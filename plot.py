import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d.axes3d import Axes3D

filepath = "optimset_NFE250k_c50k-rr500-N500-P500_biobasin_infiltrench_ppave.txt"
fid = open(filepath, 'r')
nsol = int(fid.readline())
varsvec = range(24)  # nvars

with open(filepath) as f:
    data = f.readlines()[1:2 * nsol + 1]  # this file format has all objs, then all dvs

"""decision space plot"""
pl.figure(1)
for j in range(nsol, 2*nsol):
    eachd = data[j].split()
    floats = [float(y) for y in data[j].split()]
    pl.plot(varsvec, floats)

"""objective space plot"""
costvec = []; runvec = []; Nloadvec = []; Ploadvec = [];

for k in range(1, nsol):
    # eacho = data[k].split()
    eacho = [float(y) for y in data[k].split()]
    costvec[len(costvec):] = [eacho[0]]  # show cost in $ per year
    runvec[len(runvec):] = [eacho[1]]  # show runoff reduction in cubic feet per year
    Nloadvec[len(Nloadvec):] = [eacho[1]]  # show loadings lbs per year
    Ploadvec[len(Ploadvec):] = [eacho[2]]  # show loadings lbs per year

fig = pl.figure(2)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(runvec, Nloadvec, Ploadvec, marker='o')
ax.set_xlabel('cost ($)')
ax.set_ylabel('runoff reduction (ft3)')
ax.set_zlabel('N load reduction (lbs)')

pl.show()