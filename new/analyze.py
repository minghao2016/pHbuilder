import os
import matplotlib.pyplot as plt

import loaddata as load
import universe

def plotlambdacoordinates(fileName="", plotBUF=False):
    resnameList = []; residList = []
    for residue in universe.get('d_residues'):
        if (residue.d_resname == "GLU"):
            resnameList.append("GLU")
            residList.append(residue.d_resid)
        if (residue.d_resname == "ASP"):
            resnameList.append("ASP")
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

    # Get the names and such of all the ASPs and GLUs.
    resnameList = []
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

    # Loop through all the lambdas:
    for idx in range(1, len(resnameList) + 1):
        
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

        # Get/compute some additional information.
        pH = universe.get('ph_ph')
        nstlout = universe.get('ph_nstxout')

        # Title, axes, etc.
        plt.title("{}-{} chain {}\npH = {}, nstlambda = {}".format(resnameList[idx-1], residList[idx-1], chainList[idx-1], pH, nstlout))
        plt.xlabel("Time (ps)")
        plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
        plt.ylabel(r"$\lambda$-coordinate")
        plt.grid()

        # Save.
        fileName = "{}/{}_{}-{:03d}".format(dirname, chainList[idx-1], resnameList[idx-1], residList[idx-1])
        # plt.savefig("{}.pdf".format(fileName)); os.system("pdfcrop {0}.pdf {0}.pdf >> /dev/null 2>&1".format(fileName))
        plt.savefig("{}.png".format(fileName))

        # Clear the current figure
        plt.clf()
