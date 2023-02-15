#! /usr/bin/env python3

import time
from p4 import *

var.PIVEC_MIN = 1.e-6
var.RATE_MIN = 1.e-6
var.BRLEN_MIN = 1.e-5
var.GAMMA_SHAPE_MIN = 0.12


dsetNum = int(var.argvAfterDoubleDash[0])

indxs = [str(i).zfill(2) for i in range(1,100)]
indxs.append("100")
indx = indxs[dsetNum]

aDir = "1_sim_datasets"
aLen = 400  # 400, 1500, or 8000
bDir = f"{aLen}sites"
shortFNm = f"SD{aLen}_{indx}.phy"
fullFNm = f"{aDir}/{bDir}/{shortFNm}"
print(fullFNm)
assert os.path.isfile(fullFNm)
a = func.readAndPop(fullFNm)
d = Data([a])
t = func.readAndPop("simulation_tree.tree")
t.data = d
t.name = shortFNm[:-4]
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)
strt = time.time()
t.optLogLike(optBrLens=False)
t.optLogLike(optBrLens=False)
tend = time.time()
t.time = tend - strt
t.tPickle(shortFNm[:-4])
