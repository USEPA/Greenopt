"""functions for different calculations, including:
    -  aandc_GI
    -  bound_STW_dvs
    -  GI_size_byarea
    -  GI_size_byvolume
    -  impervious
    """

import numpy as np
import math
import par
import read
import sys


def impervious(HRUs, baseareas, EIAs):
    """get data for HRUs with nonzero impervious area"""

    HRUs_imp = []
    EIAs_imp = []
    Areas_imp = []
    indexhrus = []

    """determine which HRUs are impervious and append associated data"""
    for i in range(len(HRUs)):
        if EIAs[i] > 0:
            HRUs_imp.append(HRUs[i])
            EIAs_imp.append(EIAs[i])
            Areas_imp.append(baseareas[i])
            indexhrus.append(i)

    par.EIA_imp = EIAs_imp

    return HRUs_imp, Areas_imp, indexhrus


def GI_size_byarea(nIH, GI_ID):
    """determine the GI dimensions as the square root of the constant area (par) except in cases of fixed-width GI"""

    nGI = len(GI_ID)
    GI_wid = [[] for x in range(nIH * nGI)]
    GI_len = [[] for x in range(nIH * nGI)]

    for i in range(nIH):
        for j in range(nGI):
            widthcoeff = float(par.df_GI[GI_ID[j]]['WidthCoeff'])
            if widthcoeff > 0:
                GI_wid[i * nGI + j] = np.sqrt(par.GIconstarea)
                GI_len[i * nGI + j] = np.sqrt(par.GIconstarea)
            else:
                GI_wid[i * nGI + j] = 5  # eventually read in? requires adding a fixed-width col to GI_Config
                GI_len[i * nGI + j] = par.GIconstarea / GI_wid[i * nGI + j]
                # width and length vectors are organized as HRU0-BMP0, HRU0-BMP1, HRU0-BMP2, HRU1-BMP0, HRU1-BMP1 ...

    return GI_wid, GI_len


def GI_size_byvolume(nIH, EIAs_imp, GI_ID, GI_ddepth):
    """determine the GI dimensions based on runoff volume = user-defined design depth * acre * EIA"""

    print("nIH", nIH)
    print("EIAs_imp", EIAs_imp)
    print("GI_ddepth", GI_ddepth)

    msets = len(GI_ID)
    print("msets is", msets)

    GI_wid = [[] for x in range(nIH * msets)]
    GI_len = [[] for x in range(nIH * msets)]

    for i in range(nIH):
        for j in range(msets):

            VRu = EIAs_imp[i] * 43560 * GI_ddepth[j] / 12  # 43560 ft2 in 1 acre, 12 inches in a foot

            widthcoeff = float(par.df_GI[GI_ID[j]]['WidthCoeff'])
            lengthcoeff = float(par.df_GI[GI_ID[j]]['LengthCoeff'])

            if widthcoeff > 0:
                GI_wid[i * msets + j] = (VRu / widthcoeff)**0.5
                GI_len[i * msets + j] = (VRu / lengthcoeff)**0.5

            else:
                GI_wid[i * msets + j] = 5  # eventually read in? requires adding a fixed-width col to GI_Config
                GI_len[i * msets + j] = VRu / GI_wid[i * msets + j] / lengthcoeff
                # width and length vectors are organized as HRU0-BMP0, HRU0-BMP1, HRU0-BMP2, HRU1-BMP0, HRU1-BMP1 ...

    # HRUs_imp, EIAs_imp, and Areas_imp are arrays of length equal to # of impervious HRUs
    # GI_wid, GI_len are vectors of length = (# of imp HRUs)(# of selected BMPs)
    print("GI_wid", GI_wid)
    print("GI_len", GI_len)

    return GI_wid, GI_len


def aandc_GI(nHRU, msets):

    """gets the area and cost information for each mset, puts into the matrix aandc_alls, then collects indices of the
    HRUs on which stormwater GI can be implemented, for each mset, and stores in matrix indexmat"""

    aandc_alls = np.zeros((msets, 4, nHRU))
    for m in range(msets):
        aandc_mset = np.zeros((4, nHRU))  # reinitialize each time through mset loop
        jump = par.rowsbt * (m + 1)  # note pars.row_HRU1 is for LC not STW, hence +1 for the 1st (or rather, 0th mset)
        for f in range(len(aandc_mset)):  # f: 4 placeholders for minA, maxA, capital cost, O&M cost - in that order
            aandc_mset[f] = read.Excel_LC2(par.filename_wmost, par.shtname_LC, nHRU, par.row_HRU1 + jump, par.ocol_LUMang[f])
            for n in range(nHRU):  # this for loop/if statem are in case '-9' is entered as a string in the spreadsheet
                if aandc_mset[f][n] is not None:
                    if math.isnan(aandc_mset[f][n]) is True:
                        print("check values entered in the Land Use section of the user input file")
                    else:
                        aandc_mset[f][n] = int(float(aandc_mset[f][n]))
        aandc_alls[m] = aandc_mset
        # 3d matrix organized like this: aandc_alls[index for manag. sets][minA=0, maxA=1, cost_cp=2, cost_om=3][nHRU]

    nvars_GI = 0
    indexvec_HRUGI = []
    indexmat = []
    for m in range(msets):
        for n in range(nHRU):  # recall that the user enters -9 in spreadsheet for each HRU that isn't an option for GI
            if aandc_alls[m][2][n] > -9:  # 2 for cost_cp
                indexvec_HRUGI[len(indexvec_HRUGI):] = [n]  # vector of HRU indices, built for each managed set
                nvars_GI += 1  # counts total number of HRU on which GI can be implemented
        indexmat[len(indexmat):] = [indexvec_HRUGI]  # stores the index for each managed set into a matrix
        indexvec_HRUGI = []  # clears the current indexvec_HRUGI to start empty for the next one

    # this method is no longer used as intended because the aandc_alls are not indicating skipped HRUs
    # Made a new variable "countvec" which reads the xs off the hydro worksheet.
    # I thought the whole reason I did it this way was so that the HRUs could be different for different managed sets.
    # I'm nearly certain now that the earlier version populated the Land Use tab with all the HRUs and the user could choose.

    return aandc_alls, indexmat, nvars_GI


def bound_STW_dvs(aandc_alls, indexmat, msets):
    """loops through each managed set stormwater BMP (which have different areas assoc)"""

    napplHRU = np.zeros(msets, dtype='i4')
    for m in range(msets):
        napplHRU[m] = len(indexmat[m])
        # indexmat[m] = [x - 1 for x in indexmat[m]]  # indices are pushed by 1 since not dealing with a time column

    minA = [0] * int(max(napplHRU))
    maxA = [0] * int(max(napplHRU))
    nchoix = np.zeros((msets, int(max(napplHRU))), dtype='i4')

    maxvalGI = []
    for m in range(msets):
        n = 0
        for i in indexmat[m]:  # i is the actual numbers of the applicable HRUs which differs for each mset
            minA[n] = aandc_alls[m][0][i]  # fill in the minA from each "applicable HRU" for each manag. set
            maxA[n] = aandc_alls[m][1][i]  # n starts at zero to create vectors of lengths appropriate for the applHRUs
            nchoix[m][n] = (maxA[n] - minA[n]) / par.unitarea

            if nchoix[m][n] <= 0:
                sys.exit("Check user inputs! No difference between min, max area values in managed set #%d.\n" % (m+1))

            maxvalGI[len(maxvalGI):] = [nchoix[m][n]]
            n += 1

    return maxvalGI


def condense(df, hoursindex):
    """this function condenses the timeseries it is fed to only include the times in hoursindex, which was determined
    previously to reflect only nonzero values"""
    df_condensed = df.reindex(index=hoursindex, fill_value=0.0)

    return df_condensed


def modelLID_wq(df_vols, BMP_Ai, BMP_datacol, df_conccond, percentrem_LID, kcoeff):
    """this function models the water quality within the BMP and the determines the load reduction based on first-order
    decay and / or percent reduction as specified by the user in the pars file"""

    soildepth = BMP_datacol['SoilThick'] / 12
    BMP_Ai_ft2 = BMP_Ai

    volEnter_hr = df_vols['Vin'].divide(12).multiply(BMP_Ai_ft2).multiply(28.3168)  # final volume L/hr
    nhrs = len(volEnter_hr)
    # calculate the volume exiting the BMP)
    volExit_hr = df_vols['Vout']
    volExit_hr = volExit_hr.divide(12).multiply(BMP_Ai_ft2).multiply(28.3168)  # vol of water exiting in ft, then convert from ft3 to L
    #  reindex the df_swmm to contain only the times of the LID df (for memory savings)

    volExitbutnotInfilt_hr = df_vols['Voutnotinf']

    # calculate the change in volume in the LID, mass entering the LID:
    initvol_inBMP = soildepth * BMP_Ai_ft2 * float(BMP_datacol['SoilFC']) * 28.3168
    deltavol_inBMP = np.subtract(volEnter_hr, volExit_hr)
    massEnter = df_conccond.multiply(volEnter_hr)  # vectors should be same length now

    vol_inBMP = deltavol_inBMP.cumsum()
    vol_inBMP += initvol_inBMP

    # initialize variables and set initial conditions
    # note that if volExit at t0 is 0 then accordingly, no mass exits. If there happens to be enough flow that water
    # exits during the first timestep, then the conc of that exiting water is equal to the washoff conc, decayed
    massPrior_inBMP = 0
    mass_inBMP = np.zeros(nhrs)
    massExit = np.zeros(nhrs)
    massExit[0] = volExit_hr.iloc[0] * np.exp(-kcoeff * 1) * df_conccond.iloc[0]  # 1 for 1 hour timestep
    print("massExit[0]", massExit[0])

    # calculate mass in BMP, solve for conc in BMP, solve for new conc after 1st order decay
    for t in range(1, nhrs):
        # mass_inBMP = massPrior_inBMP + massEnter[t] - massExit[t-1]
        mass_inBMP[t] = massPrior_inBMP + massEnter[t] - massExit[t - 1]
        if vol_inBMP[t] > 0:
            conc_inBMP = mass_inBMP[t] / vol_inBMP[t]
        else:
            conc_inBMP = 0
        conc_decay = (1 - percentrem_LID) * conc_inBMP * np.exp(-kcoeff * 1)  # 1st-order decay, 1 for 1 hour timestep
        massExit[t] = conc_decay * volExitbutnotInfilt_hr[t] # volExit_hr[t]
        massPrior_inBMP = mass_inBMP[t]  # set mass prior to be the current mass for the next timestep

    loadreduced = (sum(massEnter) - sum(massExit)) / 453592  # to convert mg to lbs- mass_inBMP[nhrs-1]
    return loadreduced


def cost_GI(FPlan, msets, indexmat, aandc_alls, GI_area, nyears):

    v = 0
    c_BMP_tot = 0

    for m in range(msets):
        for n in indexmat[m]:
            # c_BMP = FPlan * (aandc_alls[m][2][n] + nyears * aandc_alls[m][3][n]) * GI_area[v]
            c_BMP = FPlan * (aandc_alls[m][2][n] + aandc_alls[m][3][n]) * GI_area[v]  # took out nyears for O&M so this isn't dependent on the # of years of simulating
            c_BMP_tot += c_BMP
            v += 1

    return c_BMP_tot


def distribload(loadreduced, msets, indexmat):
    loadreducedN   = []
    loadreducedP   = []
    loadreducedTSS = []
    loadreducedZn  = []

    if par.pollts[par.poi[0]] is 'TN':
        if len(par.poi) is 4:
            for m in range(msets):
                nHRU = len(indexmat[m])
                for i in range(nHRU):
                    loadreducedN.append(loadreduced[nHRU * m * 4 + i])
                    loadreducedP.append(loadreduced[nHRU * m * 4 + nHRU + i])
                    loadreducedTSS.append(loadreduced[nHRU * m * 4 + 2 * nHRU + i])
                    loadreducedZn.append(loadreduced[nHRU * m * 4 + 3 * nHRU + i])

        elif len(par.poi) is 3:
            if par.pollts[par.poi[1]] is 'TP':
                for m in range(msets):
                    nHRU = len(indexmat[m])
                    for i in range(nHRU):
                        if par.pollts[par.poi[2]] is 'TSS':
                            loadreducedN.append(loadreduced[nHRU * m * 3 + i])
                            loadreducedP.append(loadreduced[nHRU * m * 3 + nHRU + i])
                            loadreducedTSS.append(loadreduced[nHRU * m * 3 + 2 * nHRU + i])
                        else:
                            loadreducedN.append(loadreduced[nHRU * m * 3 + i])
                            loadreducedP.append(loadreduced[nHRU * m * 3 + nHRU + i])
                            loadreducedZn.append(loadreduced[nHRU * m * 3 + 2 * nHRU + i])
            else:
                for m in msets:
                    nHRU = len(indexmat[m])
                    for i in range(nHRU):
                        loadreducedN.append(loadreduced[nHRU * m * 3 + i])
                        loadreducedTSS.append(loadreduced[nHRU * m * 3 + nHRU + i])
                        loadreducedZn.append(loadreduced[nHRU * m * 3 + 2 * nHRU + i])

        elif len(par.poi) is 2:
            if par.pollts[par.poi[1]] is 'TP':
                for m in range(msets):
                    nHRU = len(indexmat[m])
                    for i in range(nHRU):
                        loadreducedN.append(loadreduced[nHRU * m *2 + i])
                        loadreducedP.append(loadreduced[nHRU * m * 2 + nHRU + i])
            elif par.pollts[par.poi[1]] is 'TSS':
                for m in msets:
                    nHRU = len(indexmat[m])
                    for i in range(nHRU):
                        loadreducedN.append(loadreduced[nHRU * m *2 + i])
                        loadreducedTSS.append(loadreduced[nHRU * m * 2 + nHRU + i])
            else:
                for m in range(msets):
                    nHRU = len(indexmat[m])
                    for i in range(nHRU):
                        loadreducedN.append(loadreduced[nHRU * m * 2 + i])
                        loadreducedZn.append(loadreduced[nHRU * m * 2 + nHRU + i])
        else:
            for m in range(msets):
                nHRU = len(indexmat[m])
                for i in range(nHRU):
                    loadreducedN.append(loadreduced[nHRU * m + i])

    elif par.pollts[par.poi[0]] is 'TP':
        if len(par.poi) is 3:
            for m in range(msets):
                nHRU = len(indexmat[m])
                for i in range(nHRU):
                    loadreducedP.append(loadreduced[nHRU * m * 3 + i])
                    loadreducedTSS.append(loadreduced[nHRU * m * 3 + nHRU + i])
                    loadreducedZn.append(loadreduced[nHRU * m * 3 + 2 * nHRU + i])

        elif len(par.poi) is 2:
            if par.pollts[par.poi[1]] is 'TSS':
                for m in range(msets):
                    nHRU = len(indexmat[m])
                    for i in range(nHRU):
                        loadreducedP.append(loadreduced[nHRU * m *2 + i])
                        loadreducedTSS.append(loadreduced[nHRU * m * 2 + nHRU + i])
            else:
                for m in range(msets):
                    nHRU = len(indexmat[m])
                    for i in range(nHRU):
                        loadreducedP.append(loadreduced[nHRU * m *2 + i])
                        loadreducedZn.append(loadreduced[nHRU * m * 2 + nHRU + i])
        else:
            for m in range(msets):
                nHRU = len(indexmat[m])
                for i in range(nHRU):
                    loadreducedP.append(loadreduced[nHRU * m + i])

    elif par.pollts[par.poi[0]] is 'TSS':
        if len(par.poi) is 2:
            for m in range(msets):
                nHRU = len(indexmat[m])
                for i in range(nHRU):
                    loadreducedTSS.append(loadreduced[nHRU * m * 2  + i])
                    loadreducedZn.append(loadreduced[nHRU * m * 2 + nHRU + i])
        else:
             for m in range(msets):
                nHRU = len(indexmat[m])
                for i in range(nHRU):
                    loadreducedTSS.append(loadreduced[nHRU * m + i])

    else:
        for m in range(msets):
            nHRU = len(indexmat[m])
            for i in range(nHRU):
                loadreducedZn.append(loadreduced[nHRU * m + i])

    return loadreducedN, loadreducedP, loadreducedTSS, loadreducedZn