#! /usr/bin/env python3

desc="""Do an MCMC analysis.
This script is derived from the 'Scripts for common tasks' of the P4 phyloinformatic toolkit (vers. 1.3; Foster, 2004; https://p4.nhm.ac.uk/scripts.html)

"""

import os
import sys
import time
import math
from p4 import *

THE_INPUT_MATRIX = sys.argv[1]
THE_TARGET_TREE = "simtree.tree"

read(THE_INPUT_MATRIX)
theAlignment = var.alignments[0]
d = Data()

read(THE_TARGET_TREE)

theTree = var.trees[0]
theTree.data = d

#mcmc = func.unPickleMcmc(0, d, verbose=False)

theTree.newComp(partNum=0, free=1, spec='empirical')
theTree.newRMatrix(partNum=0, free=1, spec='ones')
theTree.setNGammaCat(partNum=0, nGammaCat=4)
theTree.newGdasrv(partNum=0, free=1, val=0.5)
theTree.setPInvar(partNum=0, free=0, val=0.0)

#theTree.setModelThingsRandomly()
m = Mcmc(theTree, nChains=4, runNum=0, sampleInterval=100, checkPointInterval=200000, simulate=0, verbose=True)
m.prob.local = 0.0
m.prob.eTBR = 0.0
m.prob.root3 = 0.0
m.prob.brLen = 0.0
m.prob.allBrLens = 0.0
m.prob.rMatrixDir = 1.0
m.prob.compDir = 0.0
m.prob.allCompsDir = 1.0
m.prob.ndch2_leafCompsDir = 0.0
m.prob.ndch2_internalCompsDir = 0.0
m.prob.ndch2_internalCompsDirAlpha = 0.0
m.prob.ndch2_leafCompsDirAlpha = 0.0
m.prob.gdasrv = 1.0

#Valid chain
m.run(600000, verbose=True)
