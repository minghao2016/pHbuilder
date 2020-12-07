#!/bin/python3

# HEADER
################################################################################

import os, re, numpy as np

# I put this function explicitly here so there are no dependencies.
def grepLineNumber(fileName, word):
    for count, line in enumerate(open(fileName)):
        if (re.search(word, line)):
            return count + 1

# I put this function explicitly here so there are no dependencies.
def Float(fileName, line, col):
    for x, y in enumerate(open(fileName)):
        if x == line - 1:
            return float(y.split()[col-1])

# I put this function explicitly here so there are no dependencies.
def Col(fileName, col, start = 0, stop = 0):
    data = []
    
    try:
        for x, y in enumerate(open(fileName)):
            if start == 0 and stop == 0 and (y.split()[0][0] not in ['@','#','&']):
                data.append(float(y.split()[col-1]))

            elif (x >= start-1) and (x <= stop-1):
                data.append(float(y.split()[col-1]))
    
    except IndexError:
        pass

    return data

# ACTUAL CODE
################################################################################

def getWallTime():
    lineNumber = grepLineNumber('MD.log', 'Performance:')
    wallTime   = Float('MD.log', lineNumber, 2)

    return wallTime

repeats = 3 # Number repeats per data point, to obtain sdev.

for NumlambdaGroups in [0]:
    
    os.chdir('lam{:02d}'.format(NumlambdaGroups))
    
    for _ in range(0, repeats):
        os.system('./run.sh')
        
        with open('performance.txt', 'a') as file:
            file.write('{}\n'.format(getWallTime()))

    performanceList = Col('performance.txt', 1)

    mean = np.mean(performanceList)
    sdev = np.std(performanceList)

    os.chdir('..')

    with open('results.txt', 'a') as file:
        file.write('{} {} {}\n'.format(NumlambdaGroups, mean, sdev))
