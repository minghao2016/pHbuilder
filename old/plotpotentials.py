import numpy as np
import matplotlib.pyplot as plt

import loaddata as load

def plotpotentials(pH, pKa):
    R   = 8.3145 * 10**-3 # "kJ * mol⁻1 * K^-1"
    T   = 300

    lambda_i = load.Col("lambda_dwp.dat", 1, 942, 2062)
    V_bias   = load.Col("lambda_dwp.dat", 0, 942, 2062)

    V_pH = []
    for i in lambda_i:
        V_pH.append(R * T * np.log(10) * (pKa - pH) * i)

    V_comb = []
    for i in range(0, len(lambda_i)):
        V_comb.append(V_bias[i] + V_pH[i])

    plt.plot(lambda_i, V_bias, color="b", linestyle='--', label="$V_{bias}}$")
    plt.plot(lambda_i, V_pH, color="b", linestyle = ':', label="$V_{pH}$")
    plt.plot(lambda_i, V_comb, color="b", label="$V_{combined}$")

    plt.xlabel(r"$\lambda$-coordinate")
    plt.ylabel(r"$V$ (kJ/mol)")
    plt.grid()
    plt.legend()
    plt.show()
