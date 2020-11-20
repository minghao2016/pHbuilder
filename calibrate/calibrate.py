#!/bin/python

import loaddata as load

def mean(array):
    return sum(array) / len(array)

dVdlList   = load.Col("lambda_1.dat", 3)
meandVdl   = mean(dVdlList)
lambdaInit = load.Float("constant_ph_input.dat", 25, 3) # initial lambda value from .dat file

with open('calibrate.out', 'a') as file:
    file.write('%s %s\n' % (lambdaInit, meandVdl))
