#!/bin/python

import loaddata as load

import sys

def mean(array):
    return sum(array) / len(array)

dVdlList   = load.Col("lambda_1.dat", 3)
meandVdl   = mean(dVdlList)
lambdaInit = float(sys.argv[1])

with open('calibrate.out', 'a') as file:
    file.write('%s %s\n' % (lambdaInit, meandVdl))
