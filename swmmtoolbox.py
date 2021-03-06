
# Example package with a console entry point
"""
Reads and formats data from the SWMM 5 output file.
"""
from __future__ import absolute_import
from __future__ import print_function

import sys
import struct
import datetime
import os

import mando
import pandas as pd
from six.moves import range
from six.moves import zip

from tstoolbox import tsutils

PROPCODE = {0: {1: 'Area',
               },
            1: {0: 'Type',
                2: 'Inv_elev',
                3: 'Max_depth'
               },
            2: {0: 'Type',
                4: 'Inv_offset',
                3: 'Max_depth',
                5: 'Length'
               }
           }

# Names for the 'Node type' and 'Link type' codes above
TYPECODE = {0: {1: 'Area',
               },
            1: {0: 'Junction',  # nodes
                1: 'Outfall',
                2: 'Storage',
                3: 'Divider'
               },
            2: {0: 'Conduit',   # links
                1: 'Pump',
                2: 'Orifice',
                3: 'Weir',
                4: 'Outlet'
               }
           }

VARCODE = {0: {0: 'Rainfall',
               1: 'Snow_depth',
               2: 'Evaporation_loss',
               3: 'Infiltration_loss',
               4: 'Runoff_rate',
               5: 'Groundwater_outflow',
               6: 'Groundwater_elevation',
               7: 'Soil_moisture',
               8: 'Pollutant_washoff',
              },
           1: {0: 'Depth_above_invert',
               1: 'Hydraulic_head',
               2: 'Volume_stored_ponded',
               3: 'Lateral_inflow',
               4: 'Total_inflow',
               5: 'Flow_lost_flooding',
              },
           2: {0: 'Flow_rate',
               1: 'Flow_depth',
               2: 'Flow_velocity',
               3: 'Froude_number',
               4: 'Capacity',
              },
           4: {0: 'Air_temperature',
               1: 'Rainfall',
               2: 'Snow_depth',
               3: 'Evaporation_infiltration',
               4: 'Runoff',
               5: 'Dry_weather_inflow',
               6: 'Groundwater_inflow',
               7: 'RDII_inflow',
               8: 'User_direct_inflow',
               9: 'Total_lateral_inflow',
               10: 'Flow_lost_to_flooding',
               11: 'Flow_leaving_outfalls',
               12: 'Volume_stored_water',
               13: 'Evaporation_rate',
               14: 'Potential_PET',
              }
          }

# Prior to 5.10.10
VARCODE_old = {0: {0: 'Rainfall',
                   1: 'Snow_depth',
                   2: 'Evaporation_loss',
                   3: 'Runoff_rate',
                   4: 'Groundwater_outflow',
                   5: 'Groundwater_elevation',
                  },
               1: {0: 'Depth_above_invert',
                   1: 'Hydraulic_head',
                   2: 'Volume_stored_ponded',
                   3: 'Lateral_inflow',
                   4: 'Total_inflow',
                   5: 'Flow_lost_flooding',
                  },
               2: {0: 'Flow_rate',
                   1: 'Flow_depth',
                   2: 'Flow_velocity',
                   3: 'Froude_number',
                   4: 'Capacity',
                  },
               4: {0: 'Air_temperature',
                   1: 'Rainfall',
                   2: 'Snow_depth',
                   3: 'Evaporation_infiltration',
                   4: 'Runoff',
                   5: 'Dry_weather_inflow',
                   6: 'Groundwater_inflow',
                   7: 'RDII_inflow',
                   8: 'User_direct_inflow',
                   9: 'Total_lateral_inflow',
                   10: 'Flow_lost_to_flooding',
                   11: 'Flow_leaving_outfalls',
                   12: 'Volume_stored_water',
                   13: 'Evaporation_rate',
                  }
              }

# FLOWUNITS is here, but currently not used.
FLOWUNITS = {
    0: 'CFS',
    1: 'GPM',
    2: 'MGD',
    3: 'CMS',
    4: 'LPS',
    5: 'LPD'
    }


class SwmmExtract():
    def __init__(self, filename):

        self.RECORDSIZE = 4

        self.fp = open(filename, 'rb')

        self.fp.seek(-6*self.RECORDSIZE, 2)

        self.NamesStartPos, \
            self.PropertiesStartPos, \
            self.ResultsStartPos, \
            self.nperiods, \
            errcode, \
            magic2 = struct.unpack('6i', self.fp.read(6*self.RECORDSIZE))

        self.fp.seek(0, 0)
        magic1 = struct.unpack('i', self.fp.read(self.RECORDSIZE))[0]

        if magic1 != 516114522:
            print('First magic number incorrect.')
            sys.exit(1)
        if magic2 != 516114522:
            print('Second magic number incorrect.')
            sys.exit(1)
        if errcode != 0:
            print('Error code in output file indicates a problem with the run')
            sys.exit(1)
        if self.nperiods == 0:
            print('There are zero time periods in the output file')
            sys.exit(1)

        version, \
            self.flowunits, \
            self.nsubcatch, \
            self.nnodes, \
            self.nlinks, \
            self.npolluts = struct.unpack('6i',
                                          self.fp.read(6*self.RECORDSIZE))
        if version < 5100:
            self.varcode = VARCODE_old
        else:
            self.varcode = VARCODE

        self.itemlist = ['subcatchment', 'node', 'link', 'pollutant', 'system']

        # Read in the names
        self.fp.seek(self.NamesStartPos, 0)
        self.names = {0: [], 1: [], 2: [], 3: [], 4: []}
        number_list = [self.nsubcatch,
                       self.nnodes,
                       self.nlinks,
                       self.npolluts]
        for i, j in enumerate(number_list):
            for k in range(j):
                stringsize = struct.unpack('i',
                                           self.fp.read(self.RECORDSIZE))[0]
                self.names[i].append(
                    struct.unpack('{0}s'.format(stringsize),
                                  self.fp.read(stringsize))[0])

        # Stupid Python 3
        for key in self.names:
            collect_names = []
            for name in self.names[key]:
                # Why would SMMM allow spaces in names?  Anyway...
                try:
                    rname = str(name, 'ascii', 'replace')
                except TypeError:
                    rname = name.decode('ascii', 'replace')
                try:
                    collect_names.append(rname.decode())
                except AttributeError:
                    collect_names.append(rname)
            self.names[key] = collect_names

        # Read pollutant concentration codes
        # = Number of pollutants * 4 byte integers
        self.pollutant_codes = struct.unpack(
            '{0}i'.format(self.npolluts),
            self.fp.read(self.npolluts*self.RECORDSIZE))

        self.propcode = {}
        self.prop = {0: [], 1: [], 2: []}
        nsubprop = struct.unpack('i', self.fp.read(self.RECORDSIZE))[0]
        self.propcode[0] = struct.unpack(
            '{0}i'.format(nsubprop),
            self.fp.read(nsubprop*self.RECORDSIZE))
        for i in range(self.nsubcatch):
            rprops = struct.unpack(
                '{0}f'.format(nsubprop),
                self.fp.read(nsubprop*self.RECORDSIZE))
            self.prop[0].append(list(zip(self.propcode[0], rprops)))

        nnodeprop = struct.unpack('i', self.fp.read(self.RECORDSIZE))[0]
        self.propcode[1] = struct.unpack(
            '{0}i'.format(nnodeprop),
            self.fp.read(nnodeprop*self.RECORDSIZE))
        for i in range(self.nnodes):
            rprops = struct.unpack(
                'i{0}f'.format(nnodeprop - 1),
                self.fp.read(nnodeprop*self.RECORDSIZE))
            self.prop[1].append(list(zip(self.propcode[1], rprops)))

        nlinkprop = struct.unpack('i', self.fp.read(self.RECORDSIZE))[0]
        self.propcode[2] = struct.unpack(
            '{0}i'.format(nlinkprop),
            self.fp.read(nlinkprop*self.RECORDSIZE))
        for i in range(self.nlinks):
            rprops = struct.unpack(
                'i{0}f'.format(nlinkprop - 1),
                self.fp.read(nlinkprop*self.RECORDSIZE))
            self.prop[2].append(list(zip(self.propcode[2], rprops)))

        self.vars = {}
        self.nsubcatchvars = struct.unpack(
            'i', self.fp.read(self.RECORDSIZE))[0]
        self.vars[0] = struct.unpack(
            '{0}i'.format(self.nsubcatchvars),
            self.fp.read(self.nsubcatchvars*self.RECORDSIZE))

        self.nnodevars = struct.unpack('i', self.fp.read(self.RECORDSIZE))[0]
        self.vars[1] = struct.unpack(
            '{0}i'.format(self.nnodevars),
            self.fp.read(self.nnodevars*self.RECORDSIZE))

        self.nlinkvars = struct.unpack('i', self.fp.read(self.RECORDSIZE))[0]
        self.vars[2] = struct.unpack(
            '{0}i'.format(self.nlinkvars),
            self.fp.read(self.nlinkvars*self.RECORDSIZE))

        self.vars[3] = [0]

        self.nsystemvars = struct.unpack('i', self.fp.read(self.RECORDSIZE))[0]
        self.vars[4] = struct.unpack(
            '{0}i'.format(self.nsystemvars),
            self.fp.read(self.nsystemvars*self.RECORDSIZE))

        # System vars do not have names per se, but made names = number labels
        self.names[4] = [str(i) for i in self.vars[4]]

        self.startdate = struct.unpack('d', self.fp.read(2*self.RECORDSIZE))[0]
        days = int(self.startdate)
        seconds = (self.startdate - days)*86400
        self.startdate = datetime.datetime(1899, 12, 30) + \
            datetime.timedelta(days=days, seconds=seconds)

        self.reportinterval = struct.unpack(
            'i', self.fp.read(self.RECORDSIZE))[0]
        self.reportinterval = datetime.timedelta(
            seconds=self.reportinterval)

        # Calculate the bytes for each time period when
        # reading the computed results
        self.bytesperperiod = self.RECORDSIZE*(
            2 +
            self.nsubcatch*self.nsubcatchvars +
            self.nnodes*self.nnodevars +
            self.nlinks*self.nlinkvars +
            self.nsystemvars)

    def UpdateVarCode(self, typenumber):
        start = len(self.varcode[typenumber])
        end = start + len(self.names[3])
        nlabels = list(range(start, end))
        ndict = dict(list(zip(nlabels, self.names[3])))
        self.varcode[typenumber].update(ndict)

    def TypeCheck(self, itemtype):
        if itemtype in [0, 1, 2, 3, 4]:
            return itemtype
        try:
            typenumber = self.itemlist.index(itemtype)
        except ValueError:
            print('Type argument is incorrect')
            sys.exit(1)
        return typenumber

    def NameCheck(self, itemtype, itemname):
        self.itemtype = self.TypeCheck(itemtype)
        try:
            itemindex = self.names[self.itemtype].index(itemname)
        except (ValueError, KeyError):
            print('%s was not found in %s list' % (itemname, itemtype))
            sys.exit(1)
        return (itemname, itemindex)

    def GetSwmmResults(self, itemtype, name, variableindex, period):
        if itemtype not in [0, 1, 2, 4]:
            print('Type must be one of subcatchment, node. link, or system')
            sys.exit(1)

        itemname, itemindex = self.NameCheck(itemtype, name)

        date_offset = self.ResultsStartPos + period*self.bytesperperiod

        self.fp.seek(date_offset, 0)
        date = struct.unpack('d', self.fp.read(2*self.RECORDSIZE))[0]

        offset = date_offset + 2*self.RECORDSIZE  # skip the date

        if itemtype == 0:
            offset = offset + self.RECORDSIZE*(
                itemindex*self.nsubcatchvars)
        if itemtype == 1:
            offset = offset + self.RECORDSIZE*(
                self.nsubcatch*self.nsubcatchvars +
                itemindex*self.nnodevars)
        elif itemtype == 2:
            offset = offset + self.RECORDSIZE*(
                self.nsubcatch*self.nsubcatchvars +
                self.nnodes*self.nnodevars +
                itemindex*self.nlinkvars)
        elif itemtype == 4:
            offset = offset + self.RECORDSIZE*(
                self.nsubcatch*self.nsubcatchvars +
                self.nnodes*self.nnodevars +
                self.nlinks*self.nlinkvars)

        offset = offset + self.RECORDSIZE*variableindex

        self.fp.seek(offset, 0)
        value = struct.unpack('f', self.fp.read(self.RECORDSIZE))[0]
        return (date, value)


# @mando.command()
def about():
    """Display version number and system information.
    """
    tsutils.about(__name__)


# @mando.command
def catalog(filename, itemtype=''):
    ''' List the catalog of objects in output file

    :param filename: Filename of SWMM output file.
    '''
    obj = SwmmExtract(filename)
    if itemtype:
        typenumber = obj.TypeCheck(itemtype)
        plist = [typenumber]
    else:
        plist = list(range(len(obj.itemlist)))
    print('TYPE, NAME')
    for i in plist:
        for oname in obj.names[i]:
            print('{0},{1}'.format(obj.itemlist[i], oname))


# @mando.command
def listdetail(filename, itemtype, name=''):
    ''' List nodes and metadata in output file

    :param filename: Filename of SWMM output file.
    :param itemtype: Type to print out the table of
        (subcatchment, node, or link)
    :param name: Optional specfic name to print only that entry.
    '''
    obj = SwmmExtract(filename)
    typenumber = obj.TypeCheck(itemtype)
    if name:
        objectlist = [obj.NameCheck(itemtype, name)[0]]
    else:
        objectlist = obj.names[typenumber]

    propnumbers = obj.propcode[typenumber]
    headstr = ['#Name'] + [PROPCODE[typenumber][i] for i in propnumbers]
    headfmtstr = '{0:<25},{1:<8},' + ','.join(
        ['{'+str(i)+':>10}' for i in range(2, 1+len(propnumbers))])

    print(headfmtstr.format(*tuple(headstr)))
    fmtstr = '{0:<25},{1:<8},' + ','.join(
        ['{'+str(i)+':10.2f}' for i in range(2, 1+len(propnumbers))])

    for i, oname in enumerate(objectlist):
        printvar = [oname]
        for j in obj.prop[typenumber][i]:
            if j[0] == 0:
                printvar.append(TYPECODE[typenumber][j[1]])
            else:
                printvar.append(j[1])
        print(fmtstr.format(*tuple(printvar)))


# @mando.command
def listvariables(filename):
    ''' List variables available for each type
        (subcatchment, node, link, pollutant, system)

    :param filename: Filename of SWMM output file.
    '''
    obj = SwmmExtract(filename)
    print('TYPE, DESCRIPTION, VARINDEX')
    # 'pollutant' really isn't it's own itemtype
    # but part of subcatchment, node, and link...
    for itemtype in ['subcatchment', 'node', 'link', 'system']:
        typenumber = obj.TypeCheck(itemtype)

        obj.UpdateVarCode(typenumber)

        for i in obj.vars[typenumber]:
            try:
                print('{0},{1},{2}'.format(itemtype,
                                           obj.varcode[typenumber][i].decode(),
                                           i))
            except (TypeError, AttributeError):
                print('{0},{1},{2}'.format(itemtype,
                                           str(obj.varcode[typenumber][i]),
                                           str(i)))


# @mando.command
def stdtoswmm5(start_date=None, end_date=None, input_ts='-'):
    ''' Take the toolbox standard format and return SWMM5 format.

    Toolbox standard:
    Datetime, Column_Name
    2000-01-01 00:00:00 ,  45.6
    2000-01-01 01:00:00 ,  45.2
    ...

    SWMM5 format:
    ; comment line
    01/01/2000 00:00, 45.6
    01/01/2000 01:00, 45.2
    ...

    :param input_ts: Filename with data in 'ISOdate,value' format or '-' for
        stdin.  Default is stdin.
    :param start_date: The start_date of the series in ISOdatetime format, or
        'None' for beginning.
    :param end_date: The end_date of the series in ISOdatetime format, or
        'None' for end.
    '''
    import csv
    sys.tracebacklimit = 1000
    tsd = tsutils.read_iso_ts(input_ts)[start_date:end_date]
    try:
        # Header
        print(';Datetime,', ', '.join(str(i) for i in tsd.columns))

        # Data
        cols = tsd.columns.tolist()
        tsd['date_tmp_tstoolbox'] = tsd.index.format(formatter=lambda x:
                                                     x.strftime('%m/%d/%Y'))
        tsd['time_tmp_tstoolbox'] = tsd.index.format(formatter=lambda x:
                                                     x.strftime('%H:%M:%S'))
        tsd.to_csv(sys.stdout, float_format='%g', header=False, index=False,
                   cols=['date_tmp_tstoolbox', 'time_tmp_tstoolbox'] + cols,
                   sep=' ', quoting=csv.QUOTE_NONE)
    except IOError:
        return


# @mando.command
def getdata(filename, *labels):
    ''' DEPRECATED: Use 'extract' instead.
    '''
    return extract(filename, *labels)


# @mando.command
def extract(filename, *labels):
    ''' Get the time series data for a particular object and variable

    :param filename: Filename of SWMM output file.
    :param labels: The remaining arguments uniquely identify a time-series
        in the binary file.  The format is
        'TYPE,NAME,VARINDEX'.
        For example: 'node,C64,1 node,C63,1 ...'
        TYPE and NAME can be retrieved with
            'swmmtoolbox list filename.out'
        VARINDEX can be retrieved with
            'swmmtoolbox listvariables filename.out'
    '''
    obj = SwmmExtract(filename)
    jtsd = []
    for label in labels:
        itemtype, name, variableindex = label.split(',')
        typenumber = obj.TypeCheck(itemtype)
        if itemtype != 'system':
            name = obj.NameCheck(itemtype, name)[0]

        obj.UpdateVarCode(typenumber)

        begindate = datetime.datetime(1899, 12, 30)
        dates = []
        values = []
        for time in range(obj.nperiods):
            date, value = obj.GetSwmmResults(
                typenumber, name, int(variableindex), time)
            days = int(date)
            seconds = (date - days)*86400
            date = begindate + datetime.timedelta(
                days=days, seconds=seconds)
            dates.append(date)
            values.append(value)
        jtsd.append(pd.DataFrame(
            pd.Series(values),
            columns=['{0}_{1}_{2}'.format(
                itemtype, name, obj.varcode[typenumber][int(variableindex)])]))
    result = pd.concat(jtsd, axis=1).reindex(jtsd[0].index)
    return result
    # return tsutils._printiso(result)


def main():
    if not os.path.exists('debug_swmmtoolbox'):
       sys.tracebacklimit = 0
    mando.main()


if __name__ == '__main__':
    main()

