#!/bin/python3

import os, numpy as np
import loaddata as load

def getWallTime():
    lineNumber = load.grepLineNumber('MD.log', 'Performance:')
    wallTime   = load.Float('MD.log', lineNumber, 2)

    return wallTime

repeats = 3        # Number of runs per data point.
nsteps  = 50000    # Number of simulation steps to simulate.

for bool_pH in [0, 1]:
    for opMode in ['default', 'gpu', 'cpu']:    # default = let gromacs decide
                                                # gpu = -pme cpu, cpu = -nb cpu
        if (bool_pH == 1 and opMode == 'default'):
            continue                            # this combo is still broken so skip

        for _ in range(0, repeats):             # repeat measurement for set number of repeats

            os.system('./phbuilder.py %s %s %s' % (bool_pH, opMode, nsteps))
            os.system('./run.sh')

            wallTime = getWallTime()            # get the wallTime from MD.log

            with open('data_%s_%s.txt' % (bool_pH, opMode), 'a') as file:
                file.write('%s\n' % (wallTime)) # write walltime to output file

        wallTimes = load.Col('data_%s_%s.txt' % (bool_pH, opMode), 1)

        wallTime_mean = np.mean(wallTimes)
        wallTime_sdev = np.std(wallTimes)

        with open('data_final.txt', 'a') as file:
            file.write('%s %s %s %s %s %s\n' % (bool_pH, opMode, repeats, nsteps, wallTime_mean, wallTime_sdev))
