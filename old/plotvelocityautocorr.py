import os
import matplotlib.pyplot as plt
from multiprocess import Pool; THREAD_MAX = 16

import lib
import misc
import loaddata as load

def plotvelocityautocorr(fileName = ""):
    sim = lib.sim()
    sim.loadpdb(misc.inferFullName())
    
    # I - CREATE INDEX FILE
    count = 1; array = []
    file = open("autocorr.ndx", "w+")
    for residue in sim.d_residues:
        for atom in residue.d_atoms:

            if residue.d_resname == "GLU" and atom == ' OE2':
                file.write("[ GLU_{0}_OE2 ]\n".format(residue.d_resid))
                file.write("{0}\n\n".format(count))
                array.append("GLU_{0}_OE2".format(residue.d_resid))

            if residue.d_resname == "GLU" and atom == ' HE2':
                file.write("[ GLU_{0}_HE2 ]\n".format(residue.d_resid))
                file.write("{0}\n\n".format(count))
                array.append("GLU_{0}_HE2".format(residue.d_resid))

            if residue.d_resname == "ASP" and atom == ' OD2':
                file.write("[ GLU_{0}_OD2 ]\n".format(residue.d_resid))
                file.write("{0}\n\n".format(count))
                array.append("GLU_{0}_OD2".format(residue.d_resid))

            if residue.d_resname == "ASP" and atom == ' HD2':
                file.write("[ GLU_{0}_HD2 ]\n".format(residue.d_resid))
                file.write("{0}\n\n".format(count))
                array.append("GLU_{0}_HD2".format(residue.d_resid))

            count += 1

    file.close()

    # Run gmx velacc to create the velocity-autocorrelations for the atoms.
    def func(array):
        if (not os.path.isfile("{0}.xvg".format(array))):
            os.system("gmx velacc -f MD.trr -n autocorr.ndx -o {0}.xvg << EOF\n{0}\nEOF".format(array))
    
    pool = Pool(THREAD_MAX)
    pool.map(func, array, 1)
    
    countAcid = sim.protein_countRes("ASP") + sim.protein_countRes("GLU")
    
    count = 0
    for idx in range(1, countAcid + 1):

        # Create a lambda-velocity.xvg file
        t = load.Col("lambda_{0}.dat".format(idx), 1)
        v = load.Col("lambda_{0}.dat".format(idx), 5)

        with open("lambda_{0}_velocity.xvg".format(idx), "w+") as file:
            for i in range(0, len(t)):
                file.write("{0}\t{1}\n".format(t[i], v[i]))

        # Run gmx analyze to create the velocity-autocorrelation for the lambda.
        if (not os.path.isfile("lambda_{0}_autocorr.xvg".format(idx))):
            os.system("gmx analyze -f lambda_{0}_velocity.xvg -ac lambda_{0}_autocorr.xvg".format(idx))

        # Plot
        a = load.Col("lambda_{0}_autocorr.xvg".format(idx), 1)
        b = load.Col("lambda_{0}_autocorr.xvg".format(idx), 2)

        x1 = load.Col("{0}.xvg".format(array[count]), 1)
        y1 = load.Col("{0}.xvg".format(array[count]), 2)
        count += 1

        x2 = load.Col("{0}.xvg".format(array[count]), 1)
        y2 = load.Col("{0}.xvg".format(array[count]), 2)
        count += 1

        plt.figure()
        plt.plot(a,  b,  label="lambda")
        plt.plot(x1, y1, label="{0}".format(array[count-2]))
        plt.plot(x2, y2, label="{0}".format(array[count-1]))

        plt.xlabel("Time (ps)")
        plt.ylabel("C(t)")

        plt.axis([0, 1, -1, 1])
        plt.legend()
        plt.grid()
        
        if (not fileName == ""):
            plt.savefig("{0}_{1}.pdf".format(fileName, idx))
            os.system("pdfcrop {0}_{1}.pdf {0}_{1}.pdf".format(fileName, idx))
        else:
            plt.show()
