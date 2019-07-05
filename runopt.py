"""this is the main file for optimization
    - loads variable values determined in runprep from the pickle object"""

from borg import *
import calcs
import numpy as np
import math
import pandas as pd
import pickle
import par
import write

use = open("pickleobj.obj", 'rb')
HRU_ID = pickle.load(use)
GI_ID = pickle.load(use)
nyears = pickle.load(use)
loadreduced_totyears = pickle.load(use)
flowdiff_totyears = pickle.load(use)

#GI_ID = ['Biofiltration w/UD', 'Sand Filter w/UD', 'Biofiltration w/UD', 'Sand Filter w/UD']
#msnames =['HRU1B_m0', 'HRU2B_m0', 'HRU3B_m0', 'HRU4B_m0', 'HRU5B_m0', 'HRU6B_m0', 'HRU7B_m0', 'HRU8B_m0', 'HRU1B_m1', 'HRU2B_m1', 'HRU3B_m1', 'HRU4B_m1', 'HRU5B_m1', 'HRU6B_m1', 'HRU7B_m1', 'HRU8B_m1', 'HRU1B_m2', 'HRU2B_m2', 'HRU3B_m2', 'HRU4B_m2', 'HRU5B_m2', 'HRU6B_m2', 'HRU7B_m2', 'HRU8B_m2', 'HRU1B_m3', 'HRU2B_m3', 'HRU3B_m3', 'HRU4B_m3', 'HRU5B_m3', 'HRU6B_m3', 'HRU7B_m3', 'HRU8B_m3']
#runred10 = [6.906139, 7.149906, 27.643884, 119.434894, 4.090303, 4.078351, 14.055731, 42.751247, 6.918703, 7.161611, 27.704689, 119.748803, 4.312163, 4.29614, 14.498194,  45.333763, 9.691565, 10.022812, 37.651649,  161.910344, 7.781426, 7.79629, 26.43928, 90.901932, 9.692647, 10.024243, 37.670989, 162.429237, 8.083145, 8.119908, 27.26221,  93.698384]
# flowdiff_totyears = pd.DataFrame(runred10, index=msnames)
#
# loadreduced_totyears = [0.0014531177688268226, 2.9381913709070697, 34.544858862676129, 197.73904417934449, 0.031631062242708473, 2.1125280119628873, 29.192159695001351, 196.19853240498887,
#                         9.3231125460020223e-05, 0.4055262973766075, 5.7702680218083788, 23.363277764201879, 0.0028549414772708294, 0.302206474245307, 4.941285162221889, 23.277988655495722,
#                         0.0014359084168929573, 2.9589257492660463, 34.789892456899928, 199.155521187571, 0.028159196603131532, 1.8280874132112894, 27.480577660843128, 193.17394439425698,
#                         9.2130763213018216e-05, 0.40860025281738827, 5.8145671939487569, 23.540932319313097, 0.00244098861287643, 0.25167866126747163, 4.584184889494459, 22.823248424702911,
#                         0.0036327944220670566, 3.5264235004757332, 41.675568874054385, 243.98480482142989, 0.12371911836087085, 3.359101678325652, 39.191517890840728, 250.20880562864915,
#                         0.00023307781365005057, 0.48039029051275167, 6.8604554697200744, 28.515326793779568, 0.010860258255956401, 0.45867978991759534, 6.447371140471458, 29.243130682933774,
#                         0.0035887605502844255, 3.5330306004422876, 41.751905261773949, 244.35279439351532, 0.1332688418541457, 3.3123077170700248, 38.960190557544117, 249.76135046883817,
#                         0.00023025262803004997, 0.48133688349786574, 6.873801355991584, 28.561414642766923, 0.011665953328123136, 0.4501173634939864, 6.3962088536079387, 29.172954946738262]
#
# print("flowdiff_totyears", flowdiff_totyears)
# print("loadreduced_totyears", loadreduced_totyears)

loadreduced = [None] * len(loadreduced_totyears)
for i in range(len(loadreduced_totyears)):
    loadreduced[i] = loadreduced_totyears[i] / nyears

flowdiff = flowdiff_totyears.divide(nyears)

print("loadreduced (div'd by nyears)", loadreduced)
print("flowdiff (div'd by nyears)", flowdiff)

nHRU = len(HRU_ID)
msets = len(GI_ID)

"""get the area and cost information for each managed set, put into the matrix aandc_alls, collect indices of the
   HRUs on which GI options can be implemented, store in the matrix indexmat"""
# aandc_alls     areas and costs for GI for each managed set, combined
# indexmat       indexvec for each managed set, combined into matrix
# nvars          number of decision variables associated with stormwater GI

aandc_alls, indexmat, nvars = calcs.aandc_GI(nHRU, msets)
print("nvars is", nvars)
print("aandc_alls", aandc_alls)

"""determine the min and max values of the decision variables"""
# maxvalGI     the maximum 'dv value' for each HRU in each managed set. actual value = dv value * unitarea
maxvalGI = calcs.bound_STW_dvs(aandc_alls, indexmat, msets)
print("maxvalGI", maxvalGI)

loadreducedN, loadreducedP, loadreducedTSS, loadreducedZn = calcs.distribload(loadreduced, msets, indexmat)
print("loadreducedN", loadreducedN)
print("loadreducedP", loadreducedP)
print("loadreducedTSS", loadreducedTSS)
print("loadreducedZn", loadreducedZn)

def greenopt(*vars):
    """runs the optimization, each candidate solution yields a set of decision variables (vars) for which objective
    function values are determined. Whether the constraints are satisfied is also evaluated."""
    """initialize variables"""
    v = 0
    runred = 0
    actarea_treated = np.zeros(nvars)
    reduction_loadN   = 0
    reduction_loadP   = 0
    reduction_loadTSS = 0
    reduction_loadZn  = 0

    """multiply decision variables by the unit increment to determine actual GI area implemented, then solve
    for objective function values"""
    for m in range(msets):
        for i in range(len(indexmat[m])):
            actarea_treated[v] = math.trunc(vars[v]) * par.unitarea   # eg for available area of 100 acres, dec vars
            runred += actarea_treated[v] * flowdiff[v] # flowdiff.iloc[v, 0]        # would be from 0-200 with 0.5-acre increments
            # runoff reduction in units of acre-ft
            if sum(loadreducedN) > 0:
                reduction_loadN   += actarea_treated[v] * loadreducedN[v]
            if sum(loadreducedP) > 0:
                reduction_loadP   += actarea_treated[v] * loadreducedP[v]
            if sum(loadreducedTSS) > 0:
                reduction_loadTSS += actarea_treated[v] * loadreducedTSS[v]
            if sum(loadreducedZn) > 0:
                reduction_loadZn  += actarea_treated[v] * loadreducedZn[v]
            v += 1

    cost_GI = calcs.cost_GI(par.FPlan, msets, indexmat, aandc_alls, actarea_treated, nyears)

    cons = [None] * par.ncons
    """evaluate constraint satisficing"""
    if reduction_loadN > par.reductiontarget_N:
        cons[0] = 0
    else:
        cons[0] = par.reductiontarget_N - reduction_loadN
        print("doesn't meet TMDL, by: ", par.reductiontarget_N - reduction_loadN)

    objs = [cost_GI, -runred, -reduction_loadN] #, -reduction_loadP]
    return objs, cons

"""use borg moea to call greenopt"""
borg = Borg(nvars, par.nobjs, par.ncons, greenopt)
borg.setEpsilons(*par.epsilons)
# borg.setBounds([0, 44], [0, 66], [0, 88], [0, 110], [0, 132], [0, 154])
borg.setBounds(*[[0, 70]] * nvars)
result = borg.solve({"maxEvaluations": par.NFE})

# print("end of optimization, results are:")
# for solution in result:
#     print(solution.getObjectives())
# for solution in result:
#     print(solution.getVariables())

write.optimset(result, "optimset.txt", nvars, msets, indexmat)


