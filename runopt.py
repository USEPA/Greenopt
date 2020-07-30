"""this is the main file for optimization
    - loads variable values determined in runprep from the pickle object"""

from borg import *
import calcs
import numpy as np
import math
import pickle
import par
import write

use = open("pickleobj_baseline.obj", 'rb')
HRU_ID = pickle.load(use)
GI_ID = pickle.load(use)
nyears = pickle.load(use)
loadreduced_totyears = pickle.load(use)
flowdiff_totyears = pickle.load(use)

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
            runred += actarea_treated[v] * flowdiff[v] # flowdiff.iloc[v, 0] # would be from 0-200 with 0.5-acre increments
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


