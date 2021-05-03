import os
import matplotlib.pyplot as plt
import numpy as np
import universe, load

def titrate(lambdaFileName, cutoff=0.80):
    lambda_proto   = 0
    lambda_deproto = 0
    
    for x in load.Col(lambdaFileName, 2):
        if (x > cutoff):
            lambda_deproto += 1
        if (x < 1 - cutoff):
            lambda_proto   += 1

    fraction = float(lambda_deproto) / (lambda_proto + lambda_deproto)

    return fraction

def plotlambda(plotBUF=False):
    resnameList = []    # Get the names and such of all the ASPs and GLUs.
    residList   = []
    for residue in universe.get('d_residues'):
        if residue.d_resname in ["ASP", "GLU"]:
            resnameList.append(residue.d_resname)
            residList.append(residue.d_resid)

    plt.figure()
    for idx in range(1, len(resnameList) + 1):
        t = load.Col("lambda_{0}.dat".format(idx), 1)
        x = load.Col("lambda_{0}.dat".format(idx), 2)
        
        plt.plot(t, x, label="%s-%s" % (resnameList[idx-1], residList[idx - 1]), linewidth=0.5)

    if (plotBUF):
        t = load.Col("lambda_{0}.dat".format(len(resnameList) + 1), 1)
        x = load.Col("lambda_{0}.dat".format(len(resnameList) + 1), 2)

        plt.plot(t, x, label="Buffer", linewidth=0.5)

    plt.xlabel("Time (ps)")
    plt.ylabel(r"$\lambda$-coordinate")

    plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
    
    plt.legend()
    plt.grid()

    if (not fileName == ""):
        plt.savefig("{0}.pdf".format(fileName))
        os.system("pdfcrop {0}.pdf {0}.pdf".format(fileName))
    else:
        plt.show()

def glicphstates():
    # EXPERIMENTAL DATA ON PROTONATION STATES AT VARIOUS PH ####################
    biophys = { # also prevost2012
        'ASP-13'  : 1,
        'ASP-31'  : 1,
        'ASP-32'  : 1,
        'ASP-49'  : 1,
        'ASP-55'  : 1,
        'ASP-86'  : 0,
        'ASP-88'  : 0,
        'ASP-91'  : 1,
        'ASP-97'  : 1,
        'ASP-115' : 1,
        'ASP-122' : 1,
        'ASP-136' : 1,
        'ASP-145' : 1,
        'ASP-153' : 1,
        'ASP-154' : 1,
        'ASP-161' : 1,
        'ASP-178' : 1,
        'ASP-185' : 1,
        'GLU-14'  : 1,
        'GLU-26'  : 0,
        'GLU-35'  : 0,
        'GLU-67'  : 0,
        'GLU-69'  : 1,
        'GLU-75'  : 0,
        'GLU-82'  : 0,
        'GLU-104' : 1,
        'GLU-147' : 1,
        'GLU-163' : 1,
        'GLU-177' : 0,
        'GLU-181' : 1,
        'GLU-222' : 1,
        'GLU-243' : 0,
        'GLU-272' : 1,
        'GLU-282' : 1
    }

    nury2010 = { # this is also cheng2010, calimet2013
        'ASP-13'  : 1,
        'ASP-31'  : 1,
        'ASP-32'  : 1,
        'ASP-49'  : 1,
        'ASP-55'  : 1,
        'ASP-86'  : 0,
        'ASP-88'  : 0,
        'ASP-91'  : 1,
        'ASP-97'  : 1,
        'ASP-115' : 1,
        'ASP-122' : 1,
        'ASP-136' : 1,
        'ASP-145' : 1,
        'ASP-153' : 1,
        'ASP-154' : 1,
        'ASP-161' : 1,
        'ASP-178' : 1,
        'ASP-185' : 1,
        'GLU-14'  : 1,
        'GLU-26'  : 0,
        'GLU-35'  : 0,
        'GLU-67'  : 0,
        'GLU-69'  : 0,
        'GLU-75'  : 0,
        'GLU-82'  : 0,
        'GLU-104' : 1,
        'GLU-147' : 1,
        'GLU-163' : 1,
        'GLU-177' : 0,
        'GLU-181' : 1,
        'GLU-222' : 1,
        'GLU-243' : 0,
        'GLU-272' : 1,
        'GLU-282' : 1
    }

    fritsch2011 = {
        'ASP-13'  : 0,
        'ASP-31'  : 0,
        'ASP-32'  : 1,
        'ASP-49'  : 1,
        'ASP-55'  : 0,
        'ASP-86'  : 0,
        'ASP-88'  : 0,
        'ASP-91'  : 0,
        'ASP-97'  : 0,
        'ASP-115' : 1,
        'ASP-122' : 1,
        'ASP-136' : 1,
        'ASP-145' : 0,
        'ASP-153' : 0,
        'ASP-154' : 0,
        'ASP-161' : 0,
        'ASP-178' : 0,
        'ASP-185' : 0,
        'GLU-14'  : 0,
        'GLU-26'  : 0,
        'GLU-35'  : 0,
        'GLU-67'  : 0,
        'GLU-69'  : 0,
        'GLU-75'  : 0,
        'GLU-82'  : 0,
        'GLU-104' : 1,
        'GLU-147' : 0,
        'GLU-163' : 0,
        'GLU-177' : 0,
        'GLU-181' : 0,
        'GLU-222' : 1,
        'GLU-243' : 0,
        'GLU-272' : 0,
        'GLU-282' : 0
    }

    lev2017 = {
        'ASP-13'  : 1,
        'ASP-31'  : 1,
        'ASP-32'  : 1,
        'ASP-49'  : 1,
        'ASP-55'  : 1,
        'ASP-86'  : 1,
        'ASP-88'  : 1,
        'ASP-91'  : 1,
        'ASP-97'  : 1,
        'ASP-115' : 1,
        'ASP-122' : 1,
        'ASP-136' : 1,
        'ASP-145' : 1,
        'ASP-153' : 1,
        'ASP-154' : 1,
        'ASP-161' : 1,
        'ASP-178' : 1,
        'ASP-185' : 1,
        'GLU-14'  : 1,
        'GLU-26'  : 0,
        'GLU-35'  : 0,
        'GLU-67'  : 0,
        'GLU-69'  : 0,
        'GLU-75'  : 0,
        'GLU-82'  : 0,
        'GLU-104' : 1,
        'GLU-147' : 1,
        'GLU-163' : 1,
        'GLU-177' : 0,
        'GLU-181' : 1,
        'GLU-222' : 1,
        'GLU-243' : 0,
        'GLU-272' : 1,
        'GLU-282' : 1
    }

    nemecz2017 = { # also Hu2018
        'ASP-13'  : 1,
        'ASP-31'  : 1,
        'ASP-32'  : 1,
        'ASP-49'  : 1,
        'ASP-55'  : 1,
        'ASP-86'  : 0,
        'ASP-88'  : 0,
        'ASP-91'  : 1,
        'ASP-97'  : 1,
        'ASP-115' : 1,
        'ASP-122' : 1,
        'ASP-136' : 1,
        'ASP-145' : 1,
        'ASP-153' : 1,
        'ASP-154' : 1,
        'ASP-161' : 1,
        'ASP-178' : 1,
        'ASP-185' : 1,
        'GLU-14'  : 1,
        'GLU-26'  : 0,
        'GLU-35'  : 0,
        'GLU-67'  : 1,
        'GLU-69'  : 1,
        'GLU-75'  : 1,
        'GLU-82'  : 1,
        'GLU-104' : 1,
        'GLU-147' : 1,
        'GLU-163' : 1,
        'GLU-177' : 1,
        'GLU-181' : 1,
        'GLU-222' : 0,
        'GLU-243' : 0,
        'GLU-272' : 1,
        'GLU-282' : 1
    }

    ullman = { # unpublished
        'ASP-13'  : 1,
        'ASP-31'  : 1,
        'ASP-32'  : 1,
        'ASP-49'  : 1,
        'ASP-55'  : 1,
        'ASP-86'  : 1,
        'ASP-88'  : 1,
        'ASP-91'  : 1,
        'ASP-97'  : 1,
        'ASP-115' : 1,
        'ASP-122' : 1,
        'ASP-136' : 1,
        'ASP-145' : 1,
        'ASP-153' : 1,
        'ASP-154' : 1,
        'ASP-161' : 1,
        'ASP-178' : 1,
        'ASP-185' : 1,
        'GLU-14'  : 1,
        'GLU-26'  : 0,
        'GLU-35'  : 0,
        'GLU-67'  : 0,
        'GLU-69'  : 0,
        'GLU-75'  : 0,
        'GLU-82'  : 1,
        'GLU-104' : 1,
        'GLU-147' : 0,
        'GLU-163' : 0,
        'GLU-177' : 0,
        'GLU-181' : 1,
        'GLU-222' : 1,
        'GLU-243' : 0,
        'GLU-272' : 0,
        'GLU-282' : 1
    }

    # DIRECTORY STRUCTURE ######################################################
    dirname = "lambdaplots" 
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    else:
        os.system("rm -f {0}/*.png {0}/*.pdf".format(dirname))

    # GET THE RESIDUE NUMBER, NAME, AND CHAIN OF ALL PROTO RESIDUES ############

    resnameList = []
    residList   = []
    chainList   = []
    for residue in universe.get('d_residues'):
        if residue.d_resname in ["ASP", "GLU"]:
            resnameList.append(residue.d_resname)
            residList.append(residue.d_resid)
            chainList.append(residue.d_chain)

    # CREATE LAMBDA PLOT FOR EVERY INDIVIDUAL PROTONATABLE RESIDUE #############

    # Loop through all the lambdas:
    # for idx in range(1, len(resnameList) + 1):
        
    #     plt.figure(figsize=(8, 6))

    #     # Update user
    #     print("plotting {}/{}".format(idx, len(resnameList)), end='\r')
        
    #     # Load columns from .dat files
    #     t = load.Col("lambda_{0}.dat".format(idx), 1)
    #     x = load.Col("lambda_{0}.dat".format(idx), 2)
        
    #     # Analyze a running sim not all columns will be equal long so trim:
    #     if len(t) > len(x):
    #         t.pop()
    #     elif len(t) < len(x):
    #         x.pop()

    #     plt.plot(t, x, linewidth=0.5)

    #     # Title
    #     plt.title("{0}-{1} in chain {2} in {3}.pdb\npH={4}, nstlambda={5}, deprotonation={6:.2f}\n\
    #         Experimentally determined state for {0}-{1} at this pH = {7}".format(
    #         resnameList[idx-1], 
    #         residList[idx-1],
    #         chainList[idx-1],
    #         universe.get('d_pdbName'),
    #         universe.get('ph_pH'),
    #         universe.get('ph_nstout'),
    #         titrate("lambda_{}.dat".format(idx)), 
    #         expVals40["{0}-{1}".format(resnameList[idx-1], residList[idx-1])]
    #         ))

    #     # Axes and stuff
    #     plt.ylim(-0.1, 1.1)
    #     plt.xlabel("Time (ps)")
    #     plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
    #     plt.ylabel(r"$\lambda$-coordinate")
    #     plt.grid()

    #     # Save.
    #     fileName = "{}/{}_{}-{:03d}".format(dirname, chainList[idx-1], resnameList[idx-1], residList[idx-1])
    #     # plt.savefig("{}.pdf".format(fileName)); os.system("pdfcrop {0}.pdf {0}.pdf >> /dev/null 2>&1".format(fileName))
    #     plt.savefig("{}.png".format(fileName))

    #     # clf = clear the entire current figure. close = closes a window.
    #     plt.clf(); plt.close()

    # CREATE HISTOGRAM PLOTS FOR COMBINED PROTO STATE OF ALL FIVE CHAINS #######
    number_of_chains   = len(set(chainList))
    residues_per_chain = int(len(resnameList) / number_of_chains)
    
    for ii in range(1, residues_per_chain + 1):
        data = []        
        for jj in range(0, number_of_chains):
            print(ii + residues_per_chain * jj, end=' ')
            data += (load.Col('lambda_{}.dat'.format(ii + residues_per_chain * jj), 2, 49713, 124320))
        print()

        # PLOTTING STUFF #######################################################

        plt.figure(figsize=(8, 6))
        plt.hist(data, density=True, bins=200)
        
        # Title
        plt.title("{0}-{1} (all chains) in {2}.pdb\npH={3}, nstlambda={4}, deprotonation={5:.2f}".format(
            resnameList[ii-1],
            residList[ii-1],
            universe.get('d_pdbName'),
            universe.get('ph_pH'),
            universe.get('ph_nstout'),
            titrate("lambda_{}.dat".format(ii))
            # expVals40["{0}-{1}".format(resnameList[ii-1], residList[ii-1])]
            ))

        # Axes and stuff
        plt.axis([-0.1, 1.1, -0.1, 12])
        plt.xlabel(r"$\lambda$-coordinate")
        plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
        plt.grid()

        # Add green vertical line indicating experimental value
        plt.vlines(x=biophys["{0}-{1}".format(resnameList[ii-1], residList[ii-1])], ymin=0, ymax=12, color='r', linewidth=4.0, label="biophysics.se/Prevost2012 = {}".format(biophys["{0}-{1}".format(resnameList[ii-1], residList[ii-1])]))
        plt.vlines(x=nury2010["{0}-{1}".format(resnameList[ii-1], residList[ii-1])], ymin=0, ymax=10, color='g', linewidth=4.0, label="Nury2010/Cheng2010/Calimet2013 = {}".format(nury2010["{0}-{1}".format(resnameList[ii-1], residList[ii-1])]))
        plt.vlines(x=fritsch2011["{0}-{1}".format(resnameList[ii-1], residList[ii-1])], ymin=0, ymax=8, color='b', linewidth=4.0, label="Fritsch2011 = {}".format(fritsch2011["{0}-{1}".format(resnameList[ii-1], residList[ii-1])]))
        plt.vlines(x=lev2017["{0}-{1}".format(resnameList[ii-1], residList[ii-1])], ymin=0, ymax=6, color='c', linewidth=4.0, label="Lev2017 = {}".format(lev2017["{0}-{1}".format(resnameList[ii-1], residList[ii-1])]))
        plt.vlines(x=nemecz2017["{0}-{1}".format(resnameList[ii-1], residList[ii-1])], ymin=0, ymax=4, color = 'm', linewidth=4.0, label="Nemecz2017/Hu2018 = {}".format(nemecz2017["{0}-{1}".format(resnameList[ii-1], residList[ii-1])]))
        plt.vlines(x=ullman["{0}-{1}".format(resnameList[ii-1], residList[ii-1])], ymin=0, ymax=2, color='y', linewidth=4.0, label="Ullman (unpublished) = {}".format(ullman["{0}-{1}".format(resnameList[ii-1], residList[ii-1])]))

        plt.legend()
        # Save and clear
        fileName = "{}/hist_{}-{:03d}".format(dirname, resnameList[ii-1], residList[ii-1])
        # plt.savefig("{}.pdf".format(fileName)); os.system("pdfcrop {0}.pdf {0}.pdf >> /dev/null 2>&1".format(fileName))
        plt.savefig('{}.png'.format(fileName))
        plt.clf(); plt.close()

def plotpotentials(pKa):
    R   = 8.3145 * 10**-3 # "kJ * mol⁻1 * K^-1"
    T   = 300

    lambda_i = load.Col("lambda_dwp.dat", 1, 942, 2062)
    V_bias   = load.Col("lambda_dwp.dat", 0, 942, 2062)
    
    pH = universe.get('ph_pH')

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

def plotforces(pKa):
    R   = 8.3145 * 10**-3 # "kJ * mol⁻1 * K^-1"
    T   = 300

    lambda_i = load.Col("lambda_dwp.dat", 1, 942, 2062)
    V_bias   = load.Col("lambda_dwp.dat", 0, 942, 2062)
    
    pH = universe.get('ph_pH')

    V_pH = []
    for i in lambda_i:
        V_pH.append(R * T * np.log(10) * (pKa - pH) * i)

    V_bias = np.gradient(V_bias)    # Take derivatives.
    V_pH   = np.gradient(V_pH)
    V_comb = [V_bias[i] + V_pH[i] for i in range(0, len(lambda_i))]

    plt.plot(lambda_i, V_bias, color="b", linestyle='--', label="$F_{bias}}$")
    plt.plot(lambda_i, V_pH, color="b", linestyle = ':', label="$F_{pH}$")
    plt.plot(lambda_i, V_comb, color="b", label="$F_{combined}$")

    plt.ylim(-0.1, 0.1)
    plt.xlabel(r"$\lambda$-coordinate")
    plt.ylabel("Force")
    plt.grid()
    plt.legend()
    plt.show()

def plothistogram(fname, bins=200):
    from scipy.stats import gaussian_kde
    
    data = load.Col(fname, 2)
    
    plt.hist(data, density=True, bins=bins)

    density = gaussian_kde(data)
    xs = np.linspace(-0.1, 1.1, bins)
    density.covariance_factor = lambda : .25
    density._compute_covariance()
    plt.plot(xs, density(xs), label="10 ns test")    

    plt.xlabel(r"$\lambda$-coordinate")
    # plt.axis([-0.1, 1.1, 0, 2.5])
    plt.xlim(-0.1, 1.1)
    plt.show()

def fitCalibration(order=5, compare=[]):
    # Get relevant stuff from universe.
    # Note: these data-members are only created when calibrate.py is ran.
    dVdlInitList = universe.get('ph_dvdl_initList')
    dVdlMeanList = universe.get('ph_dvdl_meanList')
    dVdlStdList  = universe.get('ph_dvdl_stdList')

    # Compute dV/dl coefficients.
    coeffs = np.polyfit(dVdlInitList, dVdlMeanList, order)[::-1]

    # Update user with the coefficients.
    print(coeffs)

    # Plot the computed values.
    plt.scatter(dVdlInitList, dVdlMeanList, label="mean dV/dl")
    plt.errorbar(dVdlInitList, dVdlMeanList, xerr=0, yerr=dVdlStdList, fmt='o', capsize=3, color='#1f77b4')

    # Our fit
    fit = []
    for i in dVdlInitList:
        value = 0
        for j in range(0, order + 1):
            value += coeffs[j] * i**j
        fit.append(value)
    plt.plot(dVdlInitList, fit, label="fit")

    # Comparison
    if len(compare) != 0:
        fit = []
        for i in dVdlInitList:
            value = 0
            for j in range(0, len(compare)):
                value += compare[j] * i**j
            fit.append(value)
        plt.plot(dVdlInitList, fit, label="compare")

    plt.title("Calibration for {}.pdb".format(universe.get('d_pdbName')))
    plt.ylabel(r"dV/d$\lambda$")
    plt.xlabel(r"$\lambda$-coordinate")
    plt.legend()
    plt.grid()
    plt.show()
