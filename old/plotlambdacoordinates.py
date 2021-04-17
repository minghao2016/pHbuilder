import os
import matplotlib.pyplot as plt

import lib
import misc
import loaddata as load

def plotlambdacoordinates(fileName="", plotBUF=False):
    sim = lib.sim()
    sim.loadpdb(misc.inferFullName())

    resNameList = []; resIdList = []
    for residue in sim.d_residues:
        if (residue.d_resname == "GLU"):
            resNameList.append("GLU")
            resIdList.append(residue.d_resid)
        if (residue.d_resname == "ASP"):
            resNameList.append("ASP")
            resIdList.append(residue.d_resid)

    plt.figure()
    for idx in range(1, len(resNameList) + 1):
        t = load.Col("lambda_{0}.dat".format(idx), 1)
        x = load.Col("lambda_{0}.dat".format(idx), 2)
        
        plt.plot(t, x, label="%s-%s" % (resNameList[idx-1], resIdList[idx - 1]), linewidth=0.5)

    if (plotBUF):
        t = load.Col("lambda_{0}.dat".format(len(resNameList) + 1), 1)
        x = load.Col("lambda_{0}.dat".format(len(resNameList) + 1), 2)

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
