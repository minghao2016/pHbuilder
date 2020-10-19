#!/usr/bin/python3

# Default
import os, scipy, numpy
from multiprocess import Pool
from matplotlib import pyplot as plt

# Self
from lib import sim
import loaddata as load

THREAD_MAX = 4

def autocorr(x):
    result = numpy.correlate(x, x, mode='full')
    return result[int(result.size/2.0):]

def normalize(x):
    largest = max(x)
    return [i / largest for i in x]

def func(indexName):
    # os.system("gmx velacc -f MD.trr -n index.ndx -o {0}.xvg << EOF\n{0}\nEOF".format(indexName))
    pass

indexName = ["GLU_17_OH", "ASP_42_OH", "ASP_59_OH"]
pool = Pool(THREAD_MAX)
pool.map(func, indexName, 1)

for idx in range(0, 3):
    # LOAD LAMBDA VELOCITY, PERFORM AUTOCORR AND FT
    aa = load.Col("lambda_%s.dat" % (idx + 1), 1)
    bb = autocorr(load.Col("lambda_%s.dat" % (idx + 1), 5))
    bb = normalize(scipy.fft(bb))

    # LOAD OH VELOCITY, PERFORM AUTOCORR, AND FT
    x = load.Col("%s.xvg" % indexName[idx], 1)
    y = load.Col("%s.xvg" % indexName[idx], 2)
    y = normalize(scipy.fft(y))

    plt.figure()
    plt.plot(x, y, label=indexName[idx])
    plt.plot(aa, bb, label="%s-lambda" % indexName[idx])

    plt.axis([0, 0.5, -0.5, 1])

    plt.title("Normalized FT of velocity-autocorrelation of {0}".format(indexName[idx]))
    plt.xlabel("Time (ps)")
    plt.ylabel("a.u.")
    plt.legend()

    plt.savefig("%s.pdf" % indexName[idx])
