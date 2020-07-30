"""functions used to write the swmm input file and to write modified external runoff and loading data timeseries,
    including:
    - external
    - writebio
    - writetrench
    - writerain
    - writepave
    - swmminp
    - names
    - optimset"""

import math
import matplotlib.pyplot as pl
import par
import pandas as pd
import time
from mpl_toolkits.mplot3d.axes3d import Axes3D # this is needed even though pycharm doesn't turn it orange


def external(HRUs_imp, extroot, evaproot, data, indexhrus, pollt):
    """take the SURO timeseries for each impervious HRU from the WMOST Timeseries file and format into separate text
    files required for running SWMM."""
    # note that skiprows and nrows are determined using the start and end dates in function Excel_H_dates
    # need to call this function before calling writeswmminp function to create the external SURO files for SWMM
    print("data", data.head())

    if extroot == par.extroot_flow:
        """reads the wmost hydrology file, gets the appropriate column, & creates the evaporation file for SWMM"""
        date = data.iloc[:, 0]  # grab the dates column
        evap = data[par.nameevapcol]
        df_evap = pd.concat([date, evap], axis=1)

        if max(evap) == 0:  # to avoid an error if there's no runoff
            print("zero evaporation")
            df_nz_evap = pd.concat([date.head(), df_evap.head()], axis=1)
        else:
            df_nz_evap = df_evap.iloc[evap.to_numpy().nonzero()]  # find the nonzero values of PET

        outfile_evap = evaproot + '.txt'
        df_nz_evap.to_csv(path_or_buf=outfile_evap, sep='\t', header=False, index=False)

        """creates the runoff file, input to SWMM as precip (per method)"""
        indexhrus = [x + 1 for x in indexhrus]  # shift up by 1 because of date column
        print('indexhrus', indexhrus)
        for i in range(len(indexhrus)):  # changed from HRUs_imp
            vals = data.iloc[:, indexhrus[i]]
            if max(vals) == 0:  # to avoid an error if there's no runoff
                print("zero runoff for HRU with indexhru = %d" % i)
                df_nz = pd.concat([date.head(), vals.head()], axis=1)
            else:
                df_temp = pd.concat([date, vals], axis=1)
                df_nz = df_temp.iloc[vals.to_numpy().nonzero()]

            outfile = extroot+str(i)+'.txt'
            df_nz.to_csv(path_or_buf=outfile, sep='\t', header=False, index=False)

    elif extroot == par.extroot_load:
        """creates the load file, input to SWMM at a node (per method)"""
        date = data.iloc[:, 0]  # grab the dates column
        indexhrus = [x + 1 for x in indexhrus]  # shift up by 1 because of date column
        for i in range(len(indexhrus)):  # changed from HRUs_imp
            vals = data.iloc[:, indexhrus[i]]
            if max(vals) == 0:  # to avoid an error if there's no load
                print("zero loading for HRU with indexhru = %d" % i)
                df_nz = pd.concat([date.head(), vals.head()], axis=1)
            else:
                df_temp = pd.concat([date, vals], axis=1)
                df_nz = df_temp.iloc[vals.to_numpy().nonzero()]

            outfile = extroot + pollt + '_' + str(i) + '.txt'
            df_nz.to_csv(path_or_buf=outfile, sep='\t', header=False, index=False)

    return


def writebio(m, pollt, GI_data, countvec):
    """for a GI selection of bioretention, writes the LID CONTROLS section of the swmm input file"""

    swmmfilename = 'swmm_' + str(m) + '_' + pollt + '.inp'
    file = open(swmmfilename, 'a')

    file.write('[LID_CONTROLS]\n')
    file.write(';;Name           Type/Layer Parameters\n')
    file.write(';;-------------- ---------- ----------\n')

    for i in range(par.n_imp):
        file.write('%s\tBC\n' % (GI_data[0] + '_' + str(i)))
        file.write('%s\tSURFACE\t%.1f\t0.0\t0.1\t1.0\t0\n' % ((GI_data[0] + '_' + str(i)), (GI_data['BermH'])))
        file.write('%s\tSOIL\t%.1f\t%.2f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\t%.3f\n' %
                   ((GI_data[0] + '_' + str(i)),
                    float(GI_data['SoilThick']),
                    float(GI_data['SoilPor']),
                    float(GI_data['SoilFC']),
                    float(GI_data['SoilWP']),
                    float(GI_data['SoilCond']),
                    float(GI_data['SoilCondSlope']),
                    float(GI_data['SoilSH']),
                    float(GI_data['PercentRem_' + pollt]),
                    float(GI_data[pollt + '_Decay'])))
        # (above) soil thickness, porosity, field cap, wilting point, conductivity, cond. slope, suction head

        file.write('%s\tSTORAGE\t%.1f\t%.1f\t%.4f\t0\n' % ((GI_data[0] + '_' + str(i)),
                                                           float(GI_data['StorThick']),
                                                           float(GI_data['StorVR']),
                                                           par.K_imp[i]))
        # (above) stor thickness, void ratio, seepage rate, clogging=0

        file.write('%s\tDRAIN\t1.0\t0.5\t%.1f\t6\n\n' % ((GI_data[0] + '_' + str(i)), float(GI_data['StorThick'])))
    file.close()

    return


def writetrench(m, pollt, GI_data, countvec):
    """for a GI selection of infiltration trench, writes the LID CONTROLS section of the swmm input file"""

    swmmfilename = 'swmm_' + str(m) + '_' + pollt + '.inp'
    file = open(swmmfilename, 'a')

    file.write('[LID_CONTROLS]\n')
    file.write(';;Name           Type/Layer Parameters\n')
    file.write(';;-------------- ---------- ----------\n')
    for i in range(par.n_imp):
        file.write('%s\tIT\n' % (GI_data[0] + '_' + str(i)))
        file.write('%s\tSURFACE\t%.1f\t0.0\t0.1\t1.0\t5\n' % ((GI_data[0] + '_' + str(i)), (float(GI_data['BermH']))))
        file.write('%s\tSTORAGE\t%.1f\t%.1f\t%.1f\t0\n' % ((GI_data[0] + '_' + str(i)),
                                                           float(GI_data['StorThick']),
                                                           float(GI_data['StorVR']),
                                                           par.K_imp[i]))
        file.write('%s\tDRAIN\t0\t0\t0\t6\n\n' % (GI_data[0] + '_' + str(i)))
    file.close()

    return


def writerain(m, pollt, GI_data, countvec):
    """for a GI selection of infiltration trench, writes the LID CONTROLS section of the swmm input file"""

    swmmfilename = 'swmm_' + str(m) + '_' + pollt + '.inp'
    file = open(swmmfilename, 'a')

    file.write('[LID_CONTROLS]\n')
    file.write(';;Name           Type/Layer Parameters\n')
    file.write(';;-------------- ---------- ----------\n')
    for i in range(par.n_imp):
        file.write('%s\tRG\n' % (GI_data[0] + '_' + str(i)))
        file.write('%s\tSURFACE\t%.1f\t0.0\t0.1\t1.0\t5\n' % ((GI_data[0] + '_' + str(i)), (float(GI_data['BermH']))))
        file.write('%s\tSOIL\t%.1f\t%.2f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\t%.3f\n' % ((GI_data[0] + '_' + str(i)),
                   float(GI_data['SoilThick']),
                   float(GI_data['SoilPor']),
                   float(GI_data['SoilFC']),
                   float(GI_data['SoilWP']),
                   float(GI_data['SoilCond']),
                   float(GI_data['SoilCondSlope']),
                   float(GI_data['SoilSH']),
                   float(GI_data['PercentRem_' + pollt]),
                   float(GI_data[pollt + '_Decay'])))
        file.write('%s\tSTORAGE\t%.1f\t%.2f\t%.3f\t0\n\n' % ((GI_data[0] + '_' + str(i)), float(GI_data['StorThick']),
                                                             float(GI_data['StorVR']), par.K_imp[i]))
    file.close()

    return

def writepave(m, pollt, GI_data, countvec):
    """for a GI selection of infiltration trench, writes the LID CONTROLS section of the swmm input file"""

    swmmfilename = 'swmm_' + str(m) + '_' + pollt + '.inp'
    file = open(swmmfilename, 'a')

    file.write('[LID_CONTROLS]\n')
    file.write(';;Name           Type/Layer Parameters\n')
    file.write(';;-------------- ---------- ----------\n')
    for i in range(par.n_imp):
        file.write('%s\tPP\n' % (GI_data[0] + '_' + str(i)))
        file.write('%s\tSURFACE\t%.1f\t0.0\t0.1\t1.0\t0\n' % ((GI_data[0] + '_' + str(i)), (GI_data['BermH'])))
        file.write('%s\tPAVEMENT\t%.1f\t%.1f\t%.1f\t%f\t%f\n' % ((GI_data[0] + '_' + str(i)), float(GI_data['PaveThick']),
                                                                 float(GI_data['PaveVR']), par.EIA_imp[i],
                                                                 float(GI_data['PavePerm']),
                                                                 float(GI_data['PaveCF'])))
        file.write('%s\tSOIL\t%.1f\t%.2f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\t%.3f\n' % ((GI_data[0] + '_' + str(i)),
                                                                                         float(GI_data['SoilThick']),
                                                                                         float(GI_data['SoilPor']),
                                                                                         float(GI_data['SoilFC']),
                                                                                         float(GI_data['SoilWP']),
                                                                                         float(GI_data['SoilCond']),
                                                                                         float(GI_data['SoilCondSlope']),
                                                                                         float(GI_data['SoilSH']),
                                                                                         float(GI_data['PercentRem_' + pollt]),
                                                                                         float(GI_data[pollt + '_Decay'])))
        # (above) soil thickness, porosity, field cap, wilting point, conductivity, cond. slope, suction head
        file.write('%s\tSTORAGE\t%.1f\t%.1f\t%.1f\t0\n' % ((GI_data[0] + '_' + str(i)),
                                                           float(GI_data['StorThick']),
                                                           float(GI_data['StorVR']),
                                                           par.K_imp[i]))
        # (above) stor thickness, void ratio, seepage rate, clogging=0
        file.write('%s\tDRAIN\t1.0\t0.5\t%.1f\t6\n\n' % ((GI_data[0] + '_' + str(i)), float(GI_data['StorThick'])))
    file.close()

    return



def swmminp(date_start, date_end, pollt, HRUs_imp, m, GI_wid, GI_area, GI_data, countvec):

    """This is a function to write the swmm input file for different 'managed sets' - ie GI and design depth combos.
    One SWMM input file per selected GI, ie this function will be called multiple times for each managed set, where
    b is the index of the current set."""

    """GI_par array: GI parameters are in this order: 0:selected_GI, 1:weirh, 2:orificeh, 3:soild, 4:por, 5:vegpar,
    6:soilinf, 7:udraind, 8:uvoidfr,9: useep - see calcs.py, determ_GI_params for more info"""

    # note to self: this swmm routine will be called in a loop with m (# of managed sets) as the iterating index
    # GI_par is determined in the determ_GI_params function in calcs.py -  it is a matrix with nrows = # of managed
    # sets, thus, need to do this: GI_par = GI_par[b][:], in order to input to the writeswmminp function below

    nBMP = int(len(GI_wid) / len(HRUs_imp))
    repstepstr = '0' + str(par.reportstep) + ':00:00'

    swmmfilename = 'swmm_' + str(m) + '_' + pollt +'.inp'       # numbered by managed set index
    file = open(swmmfilename, 'w')

    # HEADER
    time_of_run = time.strftime('%Y-%m-%d %H:%M')
    file.write('[TITLE]\n')
    file.write(';;SWMM input file written on %s\n\n' % time_of_run)

    # OPTIONS
    file.write('[OPTIONS]\n')
    file.write(';;Option\t\tValue\n')
    file.write('FLOW_UNITS\t\tCFS\n')
    file.write('INFILTRATION\t\tHORTON\n')
    file.write('FLOW_ROUTING\t\tKINWAVE\n')
    file.write('START_DATE\t\t%s\n' % date_start.strftime('%m/%d/%Y'))
    file.write('START_TIME\t\t%s\n' % date_start.strftime('%H:%M:%S'))
    file.write('REPORT_START_DATE\t\t%s\n' % date_start.strftime('%m/%d/%Y'))
    file.write('REPORT_START_TIME\t\t%s\n' % date_start.strftime('%H:%M:%S'))
    file.write('END_DATE\t\t%s\n' % date_end.strftime('%m/%d/%Y'))
    file.write('END_TIME\t\t%s\n' % date_end.strftime('%H:%M:%S'))
    file.write('SWEEP_START\t\t01/01\n')
    file.write('SWEEP_END\t\t12/31\n')
    file.write('DRY_DAYS\t\t0\n')
    file.write('REPORT_STEP\t\t%s\n' % repstepstr)  # might want to revisit the different timestep lengths
    file.write('WET_STEP\t\t00:01:00\n')  # add check to make sure the report step is >= the data's timestep
    file.write('DRY_STEP\t\t01:00:00\n')
    file.write('ROUTING_STEP\t\t0:01:00\n')
    file.write('ALLOW_PONDING\t\tNO\n\n')
    file.write('INERTIAL_DAMPING\t\tPARTIAL\n')
    file.write('VARIABLE_STEP\t\t0.75\n')
    file.write('LENGTHENING_STEP\t\t0\n')
    file.write('MIN_SURFAREA\t\t0\n')
    file.write('NORMAL_FLOW_LIMITED\t\tSLOPE\n')
    file.write('SKIP_STEADY_STATE\t\tNO\n')
    file.write('FORCE_MAIN_EQUATION\t\tH-W\n')
    file.write('LINK_OFFSETS\t\tDEPTH\n')
    file.write('MIN_SLOPE\t\t0\n\n')

    # EVAPORATION - already accounted for, omitting so evap is not recalculated / double-counted.

    # RAINGAGES
    file.write('[RAINGAGES]\n')
    file.write(';;               Rain      Time   Snow   Data\n')
    file.write(';;Name           Type      Intrvl Catch  Source\n')
    file.write(';;-------------- --------- ------ ------ ----------\n')
    # loop through the different HRUs and for each HRU, the different BMPs
    for i in range(len(HRUs_imp)):
        file.write('R_' + HRUs_imp[i] + '\t\tINTENSITY\t1:00\t1.0\tTIMESERIES SURO_' + HRUs_imp[i] + '\n')
        # note that a 1-hour timestep is hardcoded here
    file.write('\n')

    # SUBCATCHMENTS
    file.write('[SUBCATCHMENTS]\n')
    file.write(
        ';;Name           Rain Gage        Outlet        Area     %Imperv  Width    %Slope   CurbLen  SnowPack\n')
    file.write(
        ';;-------------- ---------------- ------------- -------- -------- -------- -------- -------- --------\n')
    for i in range(len(HRUs_imp)):
        file.write('S_' + HRUs_imp[i] + '\tR_' + HRUs_imp[i] + '\t\tL_' + HRUs_imp[i] + '\t\t%.1f\t%.1f\t%.4f\t0\t0\n'
                   %(1, 100, 208.71))  # area = 1 in acres and width = 208.71 in ft
        file.write('L_' + HRUs_imp[i] + '\tR_' + HRUs_imp[i] + '\t\tO_' + HRUs_imp[i] + '\t\t%.9f\t%.1f\t%.4f\t0\t0\n'
               % (GI_area[i * nBMP + m]/43560, 100, GI_wid[i * nBMP + m]))
    file.write('\n')

    # SUBAREAS
    file.write('[SUBAREAS]\n')
    file.write(';;Subcatchment   N-Imperv   N-Perv     S-Imperv   S-Perv     PctZero    RouteTo    PctRouted\n')
    file.write(';;-------------- ---------- ---------- ---------- ---------- ---------- ---------- ----------\n')
    for i in range(len(HRUs_imp)):
        file.write('S_' + HRUs_imp[i] + '\t0\t0\t0\t0\t100\tOUTLET\n')
        file.write('L_' + HRUs_imp[i] + '\t0\t0\t0\t0\t100\tOUTLET\n')
    file.write('\n')

    # INFILTRATION - already accounted for, omitting so infiltration is not recalculated / double-counted.

    # POLLUTANTS
    file.write('[POLLUTANTS]\n')
    file.write(';;         Mass    Rain       GW        I & I     Decay    Snow     Co-Pollut.   Co-Pollut.	 DWF\n')
    file.write(';;Name     Units   Concen.    Concen.   Concen.   Coeff.   Only     Name         Fraction    Concen.\n')
    file.write(';;------   ------  -------    ------    ------    ----	   ------   ---------	 --------    -------\n')
    file.write('%s\t%s\t%.4f\t0.0\t0.0\t0.00\t%s\t%s\t0.0\t0\n' % (pollt, 'MG/L', 0.0000, 'NO', '*'))
    file.write('\n')  # rain conc is zero bc ext load file
    file.close()

    # LID_CONTROLS
    GI_dict = {'Bioretention': writebio, 'Infiltration_Trench': writetrench, 'Rain_Garden': writerain,
               'Permeable_Pavement': writepave}
    GI_dict[GI_data[0]](m, pollt, GI_data, countvec)  # function call from dictonary. GI_data[0] corresponds to swmm names like 'Bioretention'

    # LID_USAGE
    file = open(swmmfilename, 'a')
    file.write('[LID_USAGE]\n')
    file.write(';;Subcat    LID    Number    Area    Width    InitSat    FromImp    ToPerv      RptFile\n')
    file.write(';;------    ---    ------    ----    -----    -------    -------    ------      ------\n')
    for i in range(len(HRUs_imp)):
        # note the detailed LID output file is only written for the last pollutant and last managed set
        if pollt is not par.pollts[par.poi[-1]]:
            print("should not be last pollutant", pollt)
            file.write('L_' + HRUs_imp[i] + '\t%s\t1\t%.1f\t%.2f\t0\t100\t0\n' % ((GI_data[0] + '_' + str(i)),
                                                                          GI_area[i * nBMP + m], GI_wid[i * nBMP + m]))
        else:
            print("should be last pollutant / managed set", pollt)
            file.write('L_' + HRUs_imp[i] + '\t%s\t1\t%.1f\t%.2f\t0\t100\t0\t%s\n' % ((GI_data[0] + '_' + str(i)),
                       GI_area[i * nBMP + m], GI_wid[i * nBMP + m], "detout_" + str(m) + '-' + str(i) + ".txt"))
    file.write('\n')

    # OUTFALLS
    file.write('[OUTFALLS]\n')
    file.write(';;Name           Elevation  Type       Stage Data       Gated    Route To\n')
    file.write(';;-------------- ---------- ---------- ---------------- -------- ----------------\n')
    for i in range(len(HRUs_imp)):
        file.write('O_' + HRUs_imp[i] + '\t0\tFREE\t\tNO\t\n')
    file.write('\n')

    # INFLOWS - used for external loads
    file.write('[INFLOWS]\n')
    for i in range(par.n_imp):
        file.write('O_' + HRUs_imp[i] + '\t' + pollt + '\tLOAD_' + pollt + '_' + HRUs_imp[i] + '\tMASS\t%f\n'
                   % (126.0 * 1))
    file.write('\n')

    # TIMESERIES
    file.write('[TIMESERIES]\n')
    file.write(';;Name           Date       Time       Value\n')
    file.write(';;-------------- ---------- ---------- ----------\n')
    for i in range(len(HRUs_imp)):
        file.write('SURO_' + HRUs_imp[i] + '\tFile\t%s' % par.extroot_flow + str(i) + '.txt\n')

    for i in range(len(HRUs_imp)):
        file.write('LOAD_' + pollt + '_' + HRUs_imp[i] + '\tFile\t%s' % par.extroot_load + pollt + '_' + str(i) + '.txt\n')
    file.write('\n')
    # file.write('PET\tFILE\t%s.txt\n' % evaproot)

    # REPORT
    file.write('[REPORT]\n')
    file.write(';;Reporting Options\n')
    file.write('INPUT      NO\n')
    file.write('CONTROLS   NO\n')
    file.write('SUBCATCHMENTS ALL\n')
    file.write('NODES ALL\n\n')
    file.write('LINKS ALL\n\n')

    return


def names(HRU_ID, msets, GI_ID, indexmat):
    HRU_sequence = []
    HRU_sequence_forpol = []
    HRU_sequence_wpol = []
    GI_sequence_forpol = []

    for m in range(msets):
        for n in range(len(par.poi)):
            for i in range(len(indexmat[m])):
                GI_sequence_forpol.append(GI_ID[m])
                HRU_sequence_forpol.append(HRU_ID[indexmat[m][i]] + "_m" + str(m))
                HRU_sequence_wpol.append(HRU_ID[indexmat[m][i]] + "_m" + str(m) + "_" + par.pollts[n])

    for m in range(msets):
        for i in range(len(indexmat[m])):
            HRU_sequence.append(HRU_ID[indexmat[m][i]] + "_m" + str(m))

    print("HRU sequence per m managed sets:\n", HRU_sequence)

    return HRU_sequence_forpol, HRU_sequence_wpol, GI_sequence_forpol, HRU_sequence


def optimset(result, filename, nvars, msets, indexmat):

    costvec = []; runvec = []
    Nloadvec = []; Ploadvec = []; TSSloadvec = []; Znloadvec = []
    nsol = 0

    for solution in result:
        nsol += 1
        objslist = solution.getObjectives()
        costvec[len(costvec):] = [objslist[0]]     # show cost in $ per year
        runvec[len(runvec):] = [-objslist[1]]      # show runoff reduction in inches over an acre of land
        Nloadvec[len(Nloadvec):] = [-objslist[2]]  # show loadings lbs per year
        # Ploadvec[len(Ploadvec):] = [-objslist[3]]  # show loadings lbs per year

    fileout = open(filename, 'w')
    # fileout.write('%d\n' % nsol)
    varsvec = range(nvars)
    actual_dvs = [None] * nvars
    g = 0
    for solution in result:
        # fileout.write('%.1f\t%.2f\t%.1f\t%.1f\t' % (costvec[g], runvec[g], Nloadvec[g], Ploadvec[g]))
        fileout.write('%.1f\t%.2f\t%.1f\t' % (costvec[g], runvec[g], Nloadvec[g]))
        g += 1
        dvs = solution.getVariables()
        v = 0
        for m in range(msets):
            for i in range(len(indexmat[m])):
                actual_dvs[v] = math.trunc(dvs[v]) * par.unitarea
                v += 1
        for u in range(nvars):
            fileout.write('%.1f\t' % actual_dvs[u])
        # decision space plot
        pl.plot(varsvec, actual_dvs)
        actual_dvs = [None] * nvars
        fileout.write('\n')
    fileout.close()

    # objective space plot
    fig = pl.figure(2)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(costvec, runvec, Nloadvec, marker='o')
    ax.set_xlabel('cost ($)')
    ax.set_ylabel('runoff reduction (inches)')
    ax.set_zlabel('N load reduction (lbs)')
    pl.show()

    return

