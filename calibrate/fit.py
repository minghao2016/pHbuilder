#!/bin/python3

import loaddata as load

import sys
import numpy as np
import matplotlib.pyplot as plt

initLambdaList = load.Col("calibrate.out", 1)
dvdlMeanList   = load.Col("calibrate.out", 2)

# print(initLambdaList) # debug
# print(dvdlMeanList)   # debug

if (len(sys.argv) == 1):
    raise Exception("Please specify fitting order (e.g. ./fit.py 3)")

order = int(sys.argv[1])

coeffs = np.polyfit(initLambdaList, dvdlMeanList, order)[::-1]
print(coeffs)

plt.scatter(initLambdaList, dvdlMeanList, label="average dV/dl value")

fit = []
for i in initLambdaList:
    value = 0
    for j in range(0, order + 1):
        value += coeffs[j] * i**j
    fit.append(value)

    # fit.append(coeffs[0] + coeffs[1]*i + coeffs[2]*i**2 + coeffs[3]*i**3 + coeffs[4]*i**4 + coeffs[5]*i**5 + coeffs[6]*i**6)
    # fit.append(coeffs[0] + coeffs[1]*i + coeffs[2]*i**2 + coeffs[3]*i**3)

plt.plot(initLambdaList, fit, label="fit")

plt.ylabel(r"dV/d$\lambda$")
plt.xlabel(r"$\lambda$-coordinate")
plt.legend()
plt.grid()
plt.show()
