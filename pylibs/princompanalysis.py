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

def runcovar():
    # Now start the principal components analysis of this trajectory.
    # This is done by building a so-called covariance matrix of the atomic
    # fluctuations. Diagonalisation of this matrix yields a set of eigenvectors
    # and eigenvalues, that describe collective modes of fluctuations of the
    # protein. The eigenvectors corresponding to the largest eigenvalues are
    # called "principal components", as they represent the largest-amplitude
    # collective motions.

    if (not os.path.isfile("eigenvec.trr")):
        os.system("gmx covar -s MD.tpr -f MD.trr << EOF\n1\n1\nEOF")

    # outputs
    #   average.pdb,        average structure (one frame)
    #   eigenval.xvg,       eigenvalues
    #   eigenvec.trr,       eigenvectors
    #   covar.log           just a log file

    # To see what type of motion the indivudual eigenvectors correspond to, we
    # filter the original trajectory and project out the part along a selected
    # eigenvector:

def runanaeig():
    os.system("gmx anaeig -v eigenvec.trr -f MD.trr -s MD.tpr -filt filter1.trr -first 1 -last 1")
