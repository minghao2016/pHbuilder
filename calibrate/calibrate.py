#!/bin/python3

import loaddata as load

import sys
import numpy as np

dVdlList   = load.Col("lambda_1.dat", 3)
meandVdl   = np.mean(dVdlList)
sdevdVdl   = np.std(dVdlList)
lambdaInit = float(sys.argv[1])

with open('calibrate.out', 'a') as file:
    file.write('%s %s %s\n' % (lambdaInit, meandVdl, sdevdVdl))
