import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

import lib
import misc
import loaddata as load

def arrhenius(pH):
    # CLEAN
    os.system("rm -f arrhenius.out")
    header = False

    # HARDCODED PARAMETERS
    pKa = {"GLU":4.25, "ASP":3.65}
    R = 8.3145 * 10**-3             # Gas constant ("kJ * mol⁻1 * K^-1")
    T = 300                         # Temperature
    nst_lambda = 1                  # Lambda output
    dt = 0.002                      # Simulation time step
    
    sim = lib.sim()                     # Load the protein.
    sim.loadpdb(misc.inferFullName())

    resNameList = []; resIdList = []    # Obtain the resnam and resdid
    for residue in sim.d_residues:      # Of the acidic residues.
        if (residue.d_resname == "GLU"):
            resNameList.append("GLU")
            resIdList.append(residue.d_resid)
        if (residue.d_resname == "ASP"):
            resNameList.append("ASP")
            resIdList.append(residue.d_resid)
                                        # Load bias potential.
    lambda_i = load.Col("lambda_dwp.dat", 1, 995, 2009)
    V_bias   = load.Col("lambda_dwp.dat", 0, 995, 2009)

    # LOOP OVER THE DIFFERENT ACIDIC RESIDUES
    for idx in range(1, len(resNameList) + 1):
        # I - Find DeltaE_01 and DeltaE_10 and associated coefficients

        V_pH = []
        for i in lambda_i: # See lambda papers for this equation.
            V_pH.append(R * T * np.log(10) * (pKa[resNameList[idx - 1]] - pH) * i)
        
        V_comb = [] # Add the potentials.
        for i in range(0, len(lambda_i)):
            V_comb.append(V_bias[i] + V_pH[i])
    
        # plt.plot(lambda_i, V_comb, color="b", label="$V_{combined}$")
        # plt.show()

        min1 = min(V_comb[:(len(V_comb)//2)])
        min2 = min(V_comb[(len(V_comb)//2):])
        maxx = max(V_comb)
        # print(min1, min2, maxx)

        DeltaE_01 = maxx - min1
        DeltaE_10 = maxx - min2
        # print(DeltaE_01, DeltaE_10)

        # coeff_01 = np.exp(- DeltaE_01 / (R * T))
        # coeff_10 = np.exp(- DeltaE_10 / (R * T))
        # print(coeff_01, coeff_10)

        # II - Find the attempt-rate A

        # Create a lambda-velocity.xvg file
        if (not os.path.isfile("lambda_{0}_velocity.xvg".format(idx))):
            t = load.Col("lambda_{0}.dat".format(idx), 1)
            v = load.Col("lambda_{0}.dat".format(idx), 5)
        
            with open("lambda_{0}_velocity.xvg".format(idx), "w+") as file:
                for i in range(0, len(t)):
                    file.write("{0}\t{1}\n".format(t[i], v[i]))
        
        # Run gmx analyze to create the velocity-autocorrelation for the lambda.
        if (not os.path.isfile("lambda_{0}_autocorr.xvg".format(idx))):
            os.system("gmx analyze -f lambda_{0}_velocity.xvg -ac lambda_{0}_autocorr.xvg".format(idx))

        Ct     = load.Col("lambda_{0}_autocorr.xvg".format(idx), 2)
        period = find_peaks(Ct)[0][1] * 10**-3 # Get location of second peak.
        # print(period)

        # A = 1 / (period * 10**-3)
        # print(A)

        # III - Find the N_01 and N_10 (number of up and down transitions)
        x = load.Col("lambda_{0}.dat".format(idx), 2)
        
        N_01 = 0; N_10 = 0; low  = False; high = False
        for coord in x:
            if (coord < 0 and (low == False and high == False)):
                low = True # the initial transition from 0.5 to 0 does not count

            if (coord > 1 and (low == False and high == False)):
                high = True # the initial transition from 0.5 to 1 does not count

            if coord < 0 and (low == False and high == True):
                N_10 += 1
                low = True; high = False
            
            if coord > 1 and (low == True and high == False):
                N_01 += 1
                low = False; high = True
        # print(N_01, N_10)

        # IV - Find the amount of time spent in the 0 and 1 states.
        time0 = 0; time1 = 0
        for coord in x:
            if (coord < 0.04 and coord > -0.04):
                time0 += 1
    
            if (coord > 0.96 and coord < 1.04):
                time1 += 1
        
        # Convert time0 and time1 from Nsteps to ns.
        time0 = (nst_lambda * time0) / (1000 / dt)
        time1 = (nst_lambda * time1) / (1000 / dt)
        # print(time0, time1)

        # Write data to file.
        with open('arrhenius.out', 'a') as file:
            if (not header):
                file.write("name, DeltaE01(kJ/mol), DeltaE10(kJ/mol), T(ps), N01, N10, time0(ns), time1(ns)\n")
                header = True
            file.write("%s-%s, %.3e, %.3e, %.3e, %.3e, %.3e, %.3e, %.3e\n" % (resNameList[idx-1], resIdList[idx-1], DeltaE_01, DeltaE_10, period, N_01, N_10, time0, time1))
