#!/bin/python3

# Default
import os, scipy, numpy
from multiprocess import Pool
from matplotlib import pyplot as plt

# Self
from lib import sim
import loaddata as load

THREAD_MAX = 16

# Load the .pdb file
# sim = sim()
# sim.loadpdb()

# Write new index file for analysis

# Run gmx velacc for every acidic residue

# Plot the FT of the 

def autocorr(x):
    result = numpy.correlate(x, x, mode='full')
    return result[int(result.size/2.0):]

def normalize(x):
    largest = max(x)
    return [i / largest for i in x]

def func(indexName):
    os.system("gmx velacc -f MD.trr -n index.ndx -o {0}.xvg << EOF\n{0}\nEOF".format(indexName))

indexName = ["GLU_17_OH", "ASP_42_OH", "ASP_59_OH"]
pool = Pool(THREAD_MAX)
pool.map(func, indexName, 1)

aa = load.Col("lambda_1.dat", 1)
bb = autocorr(load.Col("lambda_1.dat", 5))

bb = normalize(scipy.fft(bb))

x = load.Col("GLU_17_OH.xvg", 1)
y = load.Col("GLU_17_OH.xvg", 2)
y = normalize(scipy.fft(y))

plt.figure()
plt.plot(x, y, label="GLU-17 O-H")
plt.plot(aa, bb, label="GLU-17 lambda")

plt.axis([0, 0.5, -0.5, 1])

plt.title("Normalized FT of velocity-autocorrelation of {0}".format("GLU_17_OH"))
plt.xlabel("Time (ps)")
plt.ylabel("a.u.")
plt.legend()

plt.show()
