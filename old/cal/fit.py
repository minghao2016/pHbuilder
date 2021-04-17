#!/bin/python3

import loaddata as load

import sys
import numpy as np
import matplotlib.pyplot as plt

initLambdaList = load.Col("calibrate.out", 1)
dvdlMeanList   = load.Col("calibrate.out", 2)
dvdlSdevList   = load.Col("calibrate.out", 3)

# print(initLambdaList) # debug
# print(dvdlMeanList)   # debug

if (len(sys.argv) == 1):
    raise Exception("Please specify fitting order (e.g. ./fit.py 3)")

order = int(sys.argv[1])

coeffs = np.polyfit(initLambdaList, dvdlMeanList, order)[::-1]
print(coeffs)

# MEASURED POINTS
plt.scatter(initLambdaList, dvdlMeanList, label="average dV/dl value")
plt.errorbar(initLambdaList, dvdlMeanList, xerr=0, yerr=dvdlSdevList, fmt='o', capsize=3, color='#1f77b4')

# OUR FIT
fit = []
for i in initLambdaList:
    value = 0
    for j in range(0, order + 1):
        value += coeffs[j] * i**j
    fit.append(value)

plt.plot(initLambdaList, fit, label="fit")

# ORIGINAL FIT
coeffsorig = [24.685, -577.05, 137.39, -172.69] # for GLU
# coeffsorig = [37.822, -566.01, 117.97, -158.79] # for ASP
# coeffsorig = [670.1, -674.4, 83.19, -150.21]    # for BUF

orig = []
for i in initLambdaList:
    value = 0
    for j in range(0, len(coeffsorig)):
        value += coeffsorig[j] * i**j
    orig.append(value)

plt.plot(initLambdaList, orig, label="Noora original fit")

# coeffsorig2 = [19.543, -596.473, 164.418, -178.376]
# orig2 = []
# for i in initLambdaList:
#     value = 0
#     for j in range(0, len(coeffsorig2)):
#         value += coeffsorig2[j] * i**j
#     orig2.append(value)
# plt.plot(initLambdaList, orig2, label="fit boxsize -d 1")

plt.ylabel(r"dV/d$\lambda$")
plt.xlabel(r"$\lambda$-coordinate")
plt.legend()
plt.grid()
plt.show()
