# DEFAULT
import os, scipy, numpy, fnmatch
import matplotlib.pyplot as plt
from multiprocess import Pool; THREAD_MAX = 16

# USER
from lib import sim
import loaddata as load

def inferFullName():
    filtered = fnmatch.filter(os.listdir('.'), "*_MD.pdb")
    return filtered[0]

def inferName():
    filtered = fnmatch.filter(os.listdir('.'), "*_MD.pdb")
    return filtered[0][0:len(filtered[0])-7]

################################################################################

def coordinate(fileName=""):
    # Load *_MD.pdb
    xx = sim()
    xx.loadpdb(inferFullName())

    countAcid = xx.protein_countRes("ASP") + xx.protein_countRes("GLU")
    countAcid = 1
    
    plt.figure()
    for idx in range(1, countAcid + 1):
        t = load.Col("lambda_{0}.dat".format(idx), 1)
        x = scipy.fft(load.Col("lambda_{0}.dat".format(idx), 2))
        
        plt.plot(t, x, label=idx, linewidth=0.5)

    plt.xlabel("Time (ps)")
    plt.ylabel(r"$\lambda$-coordinate")

    plt.legend()
    plt.grid()

    if (not fileName == ""):
        plt.savefig("{0}.pdf".format(fileName))
        os.system("pdfcrop {0}.pdf {0}.pdf".format(fileName))
    else:
        plt.show()

def checkArrhenius(barrier, lambdaNum):
    coefficient = numpy.exp( -(barrier*10**3)/(8.3145 * 300) )

    print("Barrier energy = %.3f kJ/mol" % barrier)
    print("Gast constant  = 8.3145 kJ/(K*mol)")
    print("Temperature    = 300 K")
    print("coefficient    = %.4f" % (coefficient))

    x = load.Col("lambda_{0}.dat".format(lambdaNum), 2)
    
    N_10 = 0; N_01 = 0
    low  = False; high = False

    # Get the number of up and down transitions
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

    # Get the amount of time in the 1 and 0 states
    time0 = 0; time1 = 0

    for coord in x:
        if (coord < 0.1):
            time0 += 1
    
        if (coord > 0.9):
            time1 += 1

    nst_lambda = 1      # hardcoded
    dt = 0.002          # hardcoded

    time0 = (nst_lambda * time0) / (1000 / dt) # in ns ipv Nsteps
    time1 = (nst_lambda * time1) / (1000 / dt) # in ns ipv Nsteps

    # # Get the amount of time a transition takes
    # timer = 0
    # for coord in x:
    #     if (coord < 0.0):
    #         timer = 0
        
    #     if (coord > 0.0 and coord < 1.0):
    #         timer += 1
        
    #     if (coord > 1.0 and timer != 0):
    #         print(timer / (1000 / dt))
    #         timer = 0

    k_01  = N_01 / time0 # rate constant
    k_10  = N_10 / time1 # rate constant

    # A_01  = 1 / time0
    # A_10  = 1 / time1

    print("N_01 (up)      = {0}".format(N_01))
    print("N_10 (down)    = {0}".format(N_10))
    print("trans_total    = {0}".format(N_01 + N_10))
    print("\n")
    print("time in 0      = {0} ns".format(time0))
    print("time in 1      = {0} ns".format(time1))

    print("k_01           = {0} ns^-1".format(k_01))
    print("k_10           = {0} ns^-1".format(k_10))

    # print("A_01           = {0} ns^-1".format(A_01))
    # print("A_10           = {0} ns^-1".format(A_10))

def velocity(fileName="", axis=[0, 1, -1, 1]):
    # Load *_MD.pdb
    x = sim()
    x.loadpdb(inferFullName())
    
    # Create index file.
    count = 1; array = []
    file = open("autocorr.ndx", "w+")
    for residue in x.d_residues:
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

    countAcid = x.protein_countRes("ASP") + x.protein_countRes("GLU")
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

        plt.xlabel("Time (ns)")
        plt.ylabel("C(t)")

        plt.axis(axis)
        plt.legend()
        plt.grid()
        
        if (not fileName == ""):
            plt.savefig("{0}_{1}.pdf".format(fileName, idx))
            os.system("pdfcrop {0}_{1}.pdf {0}_{1}.pdf".format(fileName, idx))
        else:
            plt.show()

def velocity_FT(fileName="", axis=[0, 1, -1, 1]):
    # Load *_MD.pdb
    x = sim()
    x.loadpdb(inferFullName())
    
    # Create index file.
    count = 1; array = []
    file = open("autocorr.ndx", "w+")
    for residue in x.d_residues:
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

    countAcid = x.protein_countRes("ASP") + x.protein_countRes("GLU")
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
        b = scipy.fft(load.Col("lambda_{0}_autocorr.xvg".format(idx), 2))

        x1 = load.Col("{0}.xvg".format(array[count]), 1)
        y1 = scipy.fft(load.Col("{0}.xvg".format(array[count]), 2))
        count += 1

        x2 = load.Col("{0}.xvg".format(array[count]), 1)
        y2 = scipy.fft(load.Col("{0}.xvg".format(array[count]), 2))
        count += 1

        plt.figure()
        plt.plot(a,  b,  label="lambda")
        plt.plot(x1, y1, label="{0}".format(array[count-2]))
        plt.plot(x2, y2, label="{0}".format(array[count-1]))

        plt.xlabel(r"$\omega$")
        plt.ylabel("a.u.")

        plt.xlim(left=0); plt.xlim(right=1)
        plt.legend()
        plt.grid()
        
        if (not fileName == ""):
            plt.savefig("{0}_{1}.pdf".format(fileName, idx))
            os.system("pdfcrop {0}_{1}.pdf {0}_{1}.pdf".format(fileName, idx))
        else:
            plt.show()
