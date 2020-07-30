"""parameter file"""

"""folder, file, and sheet names for external inputs and data"""
excelfolder    = "_excel"
catchmname   = excelfolder + "/" + "taunton"

filename_wmost = excelfolder + "/" + "WMOSTv3_00-04.xlsm"
filename_GI    = excelfolder + "/" + "BMP-Config.csv"

fileend_char   = "_Characteristics2.csv"
fileend_H = "_hydrology.csv"
fileend_L = ["_Loadings_TN.csv", "_Loadings_TP.csv"]

shtname_H      = "Hydro"
shtname_LC     = "Land Use"
shtname_STW    = "Stormwater"

"""cell locations within the wmost spreadsheet"""
# (H) hydro worksheet
cell_date0_avail = 'N58'
cell_datef_avail = 'P58'
cell_date0_user  = 'M64'
cell_datef_user  = 'M65'

# (LU) land use worksheet
ocol_LUBase  = ['A', 'B', 'C', 'D', 'E']   # ocol identifies columns occupied in the LandUse worksheet
ocol_LChange = ['G', 'H', 'I', 'J']
ocol_LUMang  = ['C', 'D', 'E', 'F']
row_HRU1     = 9                           # first row where data values are located in the LU worksheet
col_HRU1     = 'C'                         # first column " " "
rowsbt       = 53                          # number of rows between managed sets in wmost spreadsheet
row_baseA    = 9                           # for LC1 - when the number of HRUs isn't yet known by the code, counts the
col_baseA    = 'C'                         # number of baseline areas that have been populated in the wmost spreadsheet

# (STW) stormwater worksheet, plus other inputs for SWMM
col_GI1      = 'B'
row_GI1      = 23
reportstep   = 1                           # hours
extroot_flow = 'extflow_'                  # names of files created for use by swmm
extroot_load = 'extload_'
evaproot     = 'extevap'
nameevapcol  = 'PET'
pollts       = ['TN', 'TP', 'TSS', 'ZN']   # water quality parameters
GIconstarea  = 650                         # in ft2
unitarea     = 5                           # increment of HRU area that can be treated by GI

"""info for BMP_Config file"""
cols_numeric = ['WidthCoeff', 'LengthCoeff', 'BermH', 'VegVol', 'SurfRough', 'SurfSlope', 'SwaleSlope', 'SoilThick',
                'SoilPor', 'SoilFC', 'SoilWP', 'SoilCond', 'SoilCondSlope', 'SoilSH', 'StorThick', 'StorVR', 'StorSeep',
                'StorCF', 'FlowCoeff', 'FlowExp', 'FlowOH', 'FlowDelay', 'PaveThick', 'PaveVR', 'ImpSurf', 'PavePerm',
                'PaveCF', 'MatThick', 'MatVR', 'MatRough', 'PercentRem_TN', 'PercentRem_TP', 'PercentRem_TSS',
                'PercentRem_ZN', 'TN_Decay', 'TP_Decay', 'TSS_Decay', 'ZN_Decay']
tot_GIs = 18    # should be the number of rows in the BMP_Config file, minus 1

"""optimization variables"""
poi   = [0, 1]     # poi = pollt of interest indicated by user in this order: [TN, TP, TSS, ZN] - e.g. 1 for TP
nobjs = 3
ncons = 1
NFE   = 1000000
epsilons = [1000000, 1000, 1000]
reductiontarget_N = 22627


"""simulation options"""
sizebyarea = False
getoutlet_L = False
FPlan = 1


def init():
    global df_GI
    global n_imp
    global K_imp
    global EIA_imp
    df_GI = []
    n_imp = []
    K_imp = []
    EIA_imp = []

    return