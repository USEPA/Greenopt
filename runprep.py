"""this is the main file for set-up,
- brings in user inputs from the wmost spreadsheet, including:
    - dates selected by user
    - start and end dates of available data for the watershed of interest
- reads in external data, hydrology and loading timeseries, provided through wmost
- writes swmm input file
- saves some necessary variables to pickle object for use in optimization phase"""

import calcs
import checks
import numpy as np
import os
import par
import pandas as pd
import pickle
import read
import swmmtoolbox
import write

"""initialize global variables"""
par.init()
"""read in dates selected by user"""
date_0, date_f, skiprows, nrows, countvec = read.Excel_H_dates(par.filename_wmost, par.shtname_H, par.cell_date0_avail,
                                                        par.cell_datef_avail, par.cell_date0_user, par.cell_datef_user)
print("skiprows is", skiprows)
print("countvec is", countvec)
nyears = nrows / 24 / 365.25
print("number of years in simulation is", nyears)

"""check that start & end dates are logical"""
checks.check_H_dates(date_0, date_f)

"""read in the runoff & loading timeseries data for user-selected dates"""
df_wmost_flow = read.timeseries_csv(par.catchmname, skiprows, nrows, par.fileend_H)
df_wmost_load = [None] * len(par.poi)  # is as long as the number of pollutants indicated by the user in the par file
df_conc = [None] * len(par.poi)  # is as long as the number of pollutants indicated by the user in the par file
for p in range(len(par.poi)):
    df_wmost_load[p] = read.timeseries_csv(par.catchmname, skiprows, nrows, par.fileend_L[p])

print("df_wmost_load[0].iloc[92:103]\n", df_wmost_load[0].iloc[92:103])

"""read in hydraulic response unit (HRU) user inputs"""
# HRU_ID   vec str      #  name of hydrologic response unit (HRU)
# A_base   vec double   #  baseline area for each HRU
# EIA      vec double   #  effective impervious area
# infilt   vec double   #  infiltration rate of the HRU
HRU_ID, A_base, EIA, infilt = read.Excel_LC1(par.filename_wmost, par.shtname_LC, par.col_baseA, par.row_baseA)

"""read in the stormwater (GI) inputs, the user-selected 'managed sets' ie combo of GI option and design depth"""
# GI_ID       str       names of stormwater GI
# GI_ddepth   double    design depths for GI
GI_ID, GI_ddepth = read.Excel_GI(par.filename_wmost, par.shtname_STW, par.col_GI1, par.row_GI1)
print("user-selected GI is", GI_ID)

GI_insims = []
for p in range(len(par.poi)):
    for m in range(len(GI_ID)):
        GI_insims.append(GI_ID[m])
read.GI_Config(par.filename_GI, par.tot_GIs, GI_ID)
print("df_GI\n", par.df_GI.head())

if GI_ID:
    """prep for swmm runs: determ imperv areas and calc GI size according to user-sel'd event depth for each mset"""
    # HRU_imp     hydraulic response units with nonzero effective impervious area
    # EIA_imp     the effective impervious area of those impervious subcatchments
    # A_imp       total area of the impervious subcatchments
    # indexhrus   indicies of HRUs that have nonzero impervious area
    # GI_wid      width of user-selected GI
    # GI_len      length of user-selected GI
    HRU_imp, A_imp, indexhrus = calcs.impervious(HRU_ID, A_base, EIA)
    par.n_imp = len(HRU_imp)
    print("par.n_imp", par.n_imp)

    if par.sizebyarea is True:
        GI_wid, GI_len = calcs.GI_size_byarea(par.n_imp, GI_ID)
    else:
        GI_wid, GI_len = calcs.GI_size_byvolume(par.n_imp, par.EIA_imp, GI_ID, GI_ddepth)


    GI_area = [GI_len[i] * GI_wid[i] for i in range(len(GI_len))]
    print("GI_area:", GI_area)
    # note that GI_area reflects all msets, organized as HRU0-BMP0, HRU0-BMP1, HRU0-BMP2, HRU1-BMP0, etc.
    # indexhrus is just based on whether an HRU is impervious, does not reflect the managed sets

    print("Impervious HRUs:", HRU_imp)
    print("Impervious area (in acres):", ['%.1f' % elem for elem in A_imp])

    aandc_alls, indexmat, nvars = calcs.aandc_GI(len(HRU_ID), len(GI_ID))
    print("aandc_alls\n", aandc_alls)
    print("indexmat", indexmat)

    """create names for each combination of managed set and HRU"""
    HRU_seq_wpol, HRU_seq_wpol_incl, GI_seq_wpol, HRU_seq = write.names(HRU_ID, len(GI_ID), GI_ID, indexmat)

    """write external files (SURO and PET) to use in swmm"""
    write.external(HRU_imp, par.extroot_flow, par.evaproot, df_wmost_flow, countvec, "NA")
    for p in range(len(par.poi)):
        write.external(HRU_imp, par.extroot_load, par.evaproot, df_wmost_load[p], countvec, par.pollts[p])

    """read characteristics file to get in on hydraulic conductivity for different HRU types and sort by imp HRUs"""
    hydcond = read.character_csv(par.catchmname, par.fileend_char)
    par.K_imp = [hydcond[countvec[x]] for x in range(len(indexhrus))]
    print("K_imp", par.K_imp)

    """initialize variables"""
    removal_params = []
    lengthpriorm = 0
    loadN_wGI_all = pd.DataFrame()   # holds loading timeseries (with GI) for all managed sets
    loadP_wGI_all = pd.DataFrame()   # represents the load obtained at LID subcatchment outlet O_HRUX
    df_condensedload = [None] * len(par.poi)
    flow_wGI_all_HRUL = pd.DataFrame()
    flow_wGI_all_HRUS = pd.DataFrame()
    c = 0; i = 0; a = 0
    vdf = {}
    hoursindex = {}

    for n in range(len(GI_ID)):  # msets is len(GI_ID)
        for p in range(len(par.poi)):
            print("pollutant:", par.poi[p])
            """write the swmm input file"""
            write.swmminp(date_0, date_f, par.pollts[p], HRU_imp, n, GI_wid, GI_area, par.df_GI[GI_ID[n]], countvec)

            """collect percent removal and first order decay parameters for use later"""
            removal_params.append(float(par.df_GI[GI_ID[n]]['PercentRem_' + par.pollts[p]]))
            removal_params.append(float(par.df_GI[GI_ID[n]][par.pollts[p] + '_Decay']))

            """make report file and output (binary) file names"""
            rptfile = 'swmm_' + str(n) + '_' + str(par.pollts[p]) + '.rpt'
            outfile = 'swmm_' + str(n) + '_' + str(par.pollts[p]) + '.out'

            """run swmm"""
            os.system("swmm5\t" + 'swmm_' + str(n) + '_' + str(par.pollts[p]) + '.inp\t' + rptfile + '\t' + outfile)
            # varoptions = swmmtoolbox.listvariables(outfile)

            if par.getoutlet_L is True:
                cit = 0; cend = 0
                load_wGI = pd.DataFrame()
                for i in range(len(indexmat[m])):
                    cstart = cend
                    eachlabel = 'node,' + 'O_' + HRU_ID[indexmat[m][i]] + ',6'
                    # eachlabel = 'subcatchment,' + 'S_' + HRU_ID[indexmat[m][i]] + ',8'
                    dataset = swmmtoolbox.extract(outfile, eachlabel)
                    load_wGI = pd.concat([load_wGI, dataset], axis=1)
                    dataset = []
                    cit += 1
                cend = cit
                load_wGI.columns = [HRU_seq_wpol[cstart + lengthpriorm: cend + lengthpriorm]]
                lengthpriorm += len(indexmat[m])

                if par.poi[p] is 0:
                    loadN_wGI_all = pd.concat([loadN_wGI_all, load_wGI], axis=1)
                    print("loadN_wGI_all", loadN_wGI_all.iloc[92:103])
                elif par.poi[p] is 1:
                    loadP_wGI_all = pd.concat([loadP_wGI_all, load_wGI], axis=1)
                    print("loadP_wGI_all", loadP_wGI_all.iloc[92:103])
            else: print("swmmtoolbox not used for loads")

        """read detailed LID output file for each HRU of the managed set (for the last pollutant simulated)"""
        for i in range(len(indexmat[n])):
            nHRU = len(indexmat[n])
            print("reading detailed LID report file for", HRU_seq[nHRU * n + i])
            vdf[HRU_seq[nHRU * n + i]], hoursindex[HRU_seq[nHRU * n + i]] = read.detailedLID_hourly("detout_" + str(n) + '-' + str(i) + ".txt", date_0)

        print("len(vdf)", len(vdf))
        flow_wGI_HRUL = pd.DataFrame()  # represents external runoff for 1 acre plus the impact of the GI
        flow_wGI_HRUS = pd.DataFrame()  # represents external runoff for 1 acre plus the impact of the GI
        for i in range(len(indexmat[n])):
            outfile = 'swmm_' + str(n) + '_' + str(par.pollts[0]) + '.out'
            eachlabel_HRUL = 'subcatchment,' + 'L_' + HRU_imp[i] + ',4'
            eachlabel_HRUS = 'subcatchment,' + 'S_' + HRU_imp[i] + ',4'
            dataset_HRUL = swmmtoolbox.extract(outfile, eachlabel_HRUL)
            dataset_HRUS = swmmtoolbox.extract(outfile, eachlabel_HRUS)
            flow_wGI_HRUL = pd.concat([flow_wGI_HRUL, dataset_HRUL], axis=1)
            flow_wGI_HRUS = pd.concat([flow_wGI_HRUS, dataset_HRUS], axis=1)
        flow_wGI_all_HRUL = pd.concat([flow_wGI_all_HRUL, flow_wGI_HRUL], axis=1)
        flow_wGI_all_HRUS = pd.concat([flow_wGI_all_HRUS, flow_wGI_HRUS], axis=1)

    """assign column names to indicate both HRU and managed set"""
    flow_wGI_all_HRUL.columns = HRU_seq
    flow_wGI_all_HRUS.columns = HRU_seq
    print("flow_wGI_all_HRUS\n", flow_wGI_all_HRUS.iloc[92:103])

    """sum pollutant loads for each HRU in each managed set"""
    flow_HRUL = flow_wGI_all_HRUL.sum()
    flow_HRUS = flow_wGI_all_HRUS.sum()
    print("flow_HRUS\n", flow_HRUS)
    print("flow_HRUL\n", flow_HRUL)
    flowdiff = flow_HRUS.subtract(flow_HRUL)
    print("flow reduction (cfs)\n", flowdiff)
    flowdiff = flowdiff.multiply(3600 / 43560 * 12)
    # for n in range (len(GI_ID)):
    #    nHRU = len(indexmat[n])
    #    for i in range(nHRU):
    #        flowdiff = flowdiff.multiply(3600*1/GI_area[i * len(GI_ID) + n])

    print("flow reduction (ft over the LID subcatchment)\n", flowdiff)

    df_wmost_load_relev = [None] * len(par.poi)
    df_wmost_flow_relev = pd.DataFrame()
    """reformat df_wmost_load to have the year listed first instead of day for use in condense function"""
    for p in range(len(par.poi)):
        timecol = pd.to_datetime(df_wmost_load[p].iloc[:, 0])
        timecol = pd.DatetimeIndex(timecol)
        df_wmost_load[p].index = timecol     # sets wmost timestamp col as the index col
        del df_wmost_load[p]['Timestamp']    # gets rid of now redundant wmost timestamp column
        print("sum load (just showing not using) in lbs\n", df_wmost_load[p].sum())

        # df_wmost_load includes all columns from the external load file(s) but _relev gets only the HRUs indic. by user
        for n in range(len(GI_ID)):
            for i in countvec:
                df_wmost_load_relev[p] = pd.concat([df_wmost_load_relev[p], df_wmost_load[p].iloc[:, i]], axis=1)
                print(df_wmost_load_relev[p])

        df_wmost_load_relev[p].columns = HRU_seq
        print("load in lbs", df_wmost_load_relev[p])  # .iloc[92:103]
        df_wmost_load_relev[p] *= 453592  # to convert from lbs to mg
        print("load in mg", df_wmost_load_relev[p])

    for n in range(len(GI_ID)):
        for i in countvec:
            df_wmost_flow_relev = pd.concat([df_wmost_flow_relev, df_wmost_flow.iloc[:, i+1]], axis=1)  # +1 to skip date col

    df_wmost_flow_relev.index = timecol
    df_wmost_flow_relev.columns = HRU_seq
    df_wmost_flow_relev = df_wmost_flow_relev.multiply(3600 * 28.3168)

    flow_wGI_all_HRUS.index = timecol  # put timecol as index on flow df as well (must be consistent for division below)
    flow_wGI_all_HRUS = flow_wGI_all_HRUS.multiply(3600 * 28.3168)  # convert cubic feet per second to liters per hour

    for p in range(len(par.poi)):
        df_conc[p] = df_wmost_load_relev[p].divide(df_wmost_flow_relev)
        df_conc[p] =  df_conc[p].fillna(0)
        df_conc[p] = df_conc[p].replace(np.inf, 0.0)
        # df_conc[p] = df_conc[p].apply(lambda x: [y if y < 20 else 20.0 for y in x])  # use in cases with to deal with super small denom causing unrealistically high conc.
        print("df_conc[p]", df_conc[p])

    """obtain external load for HRUs given by indexmat and then condense the load timseries based on hoursindex, which
    reflects the nonzero loads"""
    df_conccondensed = {}; v = 0;
    for n in range(len(GI_ID)):
        for p in range(len(par.poi)):
            for i in range(len(indexmat[n])):
                nHRU = len(indexmat[n])
                df_conccondensed[HRU_seq_wpol_incl[v]] = calcs.condense(df_conc[p].iloc[:, i], hoursindex[HRU_seq[nHRU * n + i]])
                print("HRU_seq_wpol_incl[v]", [HRU_seq_wpol_incl[v]])
                print("df_conccondensed", df_conccondensed[HRU_seq_wpol_incl[v]])
                v += 1

    print("percent removal and first order decay parameters:\n", removal_params)

    """simulate reductions in pollutant load in the GI due to first order decay or percent removal"""
    loadreduced = []; v = 0
    for n in range(len(GI_ID)):
        for p in range(len(par.poi)):
            nHRU = len(indexmat[n])
            for i in range(nHRU):
                print(GI_seq_wpol[v])
                pcrem = removal_params[n * 2 * len(par.poi) + p * len(par.poi)]
                kco = removal_params[n * 2 * len(par.poi) + p * len(par.poi) + 1]
                print("GI_area", GI_area[i * len(GI_ID) + n])
                print("HRU_seq[n * nHRU + i]", HRU_seq[n * nHRU + i])
                print("GI_seq_wpol[v]", GI_seq_wpol[v])
                eachloadreduced = 0
                eachloadreduced = calcs.modelLID_wq(vdf[HRU_seq[n * nHRU + i]], GI_area[i * len(GI_ID) + n],
                                                    par.df_GI[GI_seq_wpol[v]], df_conccondensed[HRU_seq_wpol_incl[v]],
                                                    pcrem, kco)
                loadreduced.append(eachloadreduced)
                print("loadreduced in loop", loadreduced)
                v += 1

    print("loadreduced (lbs)", loadreduced)

store = open("pickleobj.obj", 'wb')
pickle.dump(HRU_ID, store)
pickle.dump(GI_ID, store)
pickle.dump(nyears, store)
pickle.dump(loadreduced, store)
pickle.dump(flowdiff, store)

