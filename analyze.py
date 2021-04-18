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

def plotlambda(fileName="", plotBUF=False):
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

def glicphstates(plotBUF=False):
    expVals70 = {
        'ASP-13'  : 1.0,
        'ASP-31'  : 1.0,
        'ASP-32'  : 1.0,
        'ASP-49'  : 1.0,
        'ASP-55'  : 1.0,
        'ASP-86'  : 1.0,
        'ASP-88'  : 1.0,
        'ASP-91'  : 1.0,
        'ASP-97'  : 1.0,
        'ASP-115' : 1.0,
        'ASP-122' : 1.0,
        'ASP-136' : 1.0,
        'ASP-145' : 1.0,
        'ASP-153' : 1.0,
        'ASP-154' : 1.0,
        'ASP-161' : 1.0,
        'ASP-178' : 1.0,
        'ASP-185' : 1.0,
        'GLU-14'  : 1.0,
        'GLU-26'  : 1.0,
        'GLU-35'  : 1.0,
        'GLU-67'  : 1.0,
        'GLU-69'  : 1.0,
        'GLU-75'  : 1.0,
        'GLU-82'  : 1.0,
        'GLU-104' : 1.0,
        'GLU-147' : 1.0,
        'GLU-163' : 1.0,
        'GLU-177' : 1.0,
        'GLU-181' : 1.0,
        'GLU-222' : 1.0,
        'GLU-243' : 1.0,
        'GLU-272' : 1.0,
        'GLU-282' : 1.0
    }

    # Got this from Paul's "gethistogrambins_ASP_dist.py" scripts
    expVals40 = {
        'ASP-13'  : 0.0,
        'ASP-31'  : 0.0,
        'ASP-32'  : 0.0,
        'ASP-49'  : 0.0,
        'ASP-55'  : 1.0,
        'ASP-86'  : 1.0,
        'ASP-88'  : 0.0,
        'ASP-91'  : 0.0,
        'ASP-97'  : 0.0,
        'ASP-115' : 0.0,
        'ASP-122' : 0.0,
        'ASP-136' : 0.0,
        'ASP-145' : 0.0,
        'ASP-153' : 0.0,
        'ASP-154' : 0.0,
        'ASP-161' : 0.0,
        'ASP-178' : 0.0,
        'ASP-185' : 0.0,
        'GLU-14'  : 1.0,
        'GLU-26'  : 1.0,
        'GLU-35'  : 1.0,
        'GLU-67'  : 1.0,
        'GLU-69'  : 0.0,
        'GLU-75'  : 1.0,
        'GLU-82'  : 1.0,
        'GLU-104' : 0.0,
        'GLU-147' : 0.0,
        'GLU-163' : 0.0,
        'GLU-177' : 1.0,
        'GLU-181' : 0.0,
        'GLU-222' : 0.0,
        'GLU-243' : 1.0,
        'GLU-272' : 0.0,
        'GLU-282' : 0.0
    }

    resnameList = []    # Get the names and such of all the ASPs and GLUs.
    residList   = []
    chainList   = []
    for residue in universe.get('d_residues'):
        if residue.d_resname in ["ASP", "GLU"]:
            resnameList.append(residue.d_resname)
            residList.append(residue.d_resid)
            chainList.append(residue.d_chain)

    # debug
    # for idx in range(0, len(resnameList)):
    #     print("{}-{} in chain {}".format(resnameList[idx], residList[idx], chainList[idx]))

    # Make directory structure
    dirname = "lambdaplots" 
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    else:
        os.system("rm -f {0}/*.png {0}/*.pdf".format(dirname))

    # Loop through all the lambdas:
    for idx in range(1, len(resnameList) + 1):
        
        plt.figure(figsize=(8, 6))

        # Update user
        print("plotting {}/{}".format(idx, len(resnameList)), end='\r')
        
        # Load columns from .dat files
        t = load.Col("lambda_{0}.dat".format(idx), 1)
        x = load.Col("lambda_{0}.dat".format(idx), 2)
        
        # Analyze a running sim not all columns will be equal long so trim:
        if len(t) > len(x):
            t.pop()
        elif len(t) < len(x):
            x.pop()

        plt.plot(t, x, linewidth=0.5)

        # Title
        plt.title("{0}-{1} in chain {2} in {3}.pdb\npH={4}, nstlambda={5}, deprotonation={6:.2f}\n\
            Experimentally determined state for {0}-{1} at this pH = {7}".format(
            resnameList[idx-1], 
            residList[idx-1],
            chainList[idx-1],
            universe.get('d_pdbName'),
            universe.get('ph_pH'),
            universe.get('ph_nstout'),
            titrate("lambda_{}.dat".format(idx)), 
            expVals40["{0}-{1}".format(resnameList[idx-1], residList[idx-1])]
            ))

        # Axes and stuff
        plt.ylim(-0.1, 1.1)
        plt.xlabel("Time (ps)")
        plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
        plt.ylabel(r"$\lambda$-coordinate")
        plt.grid()

        # Save.
        fileName = "{}/{}_{}-{:03d}".format(dirname, chainList[idx-1], resnameList[idx-1], residList[idx-1])
        # plt.savefig("{}.pdf".format(fileName)); os.system("pdfcrop {0}.pdf {0}.pdf >> /dev/null 2>&1".format(fileName))
        plt.savefig("{}.png".format(fileName))

        # clf = clear the entire current figure. close = closes a window.
        plt.clf(); plt.close()

    # To implement:
        # Throw away first 5-10 ns (because this is calibration)
        # Combine latter parts of all the residues to make average over 5 chains
        # Make a spreadsheet with correctly and incorrectly predicted and see if we can find a coincidence.

    # Other stuff we need to do to improve this:
        # Recalibrate and use higher-order polynomials.
        # Less buffers, and put buffers farher away.
        # Increase barrier energy from 5.0 to 7.5 kJ/mol.

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
    dVdlInitList = universe.get('d_dVdlInitList')
    dVdlMeanList = universe.get('d_dVdlMeanList')
    dVdlStdList  = universe.get('d_dVdlStdList')

    # Compute dV/dl coefficients.
    coeffs = np.polyfit(dVdlInitList, dVdlMeanList, order)[::-1]

    # User update.
    print(coeffs)

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
