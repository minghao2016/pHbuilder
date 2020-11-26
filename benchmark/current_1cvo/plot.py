#!/bin/python3

import numpy as np
import matplotlib.pyplot as plt

import loaddata as load

meanList = load.Col('data_final.txt', 5)
sdevList = load.Col('data_final.txt', 6)

nruns = load.Int('data_final.txt', 1, 3)
nstep = load.Int('data_final.txt', 1, 4)

ind = [1, 2, 4, 3, 5] # change order of plot without changing order of data

fig, ax = plt.subplots(figsize=(9, 5))
rects   = ax.bar(ind, meanList, width=0.5, yerr=sdevList, capsize=4)

# add some text for labels, title and axes ticks
ax.set_xticks(ind)
ax.set_ylabel('Performance (ns/day)')
ax.set_title('1cvo.pdb (30\'000 atoms), rebase gromacs version, %s runs of %s steps' % (nruns, nstep))
ax.set_xticklabels(('(pH off)', '(pH off, -pme cpu)', '(pH off, -nb cpu)', '(pH on, -pme cpu)', '(pH on, -nb cpu)'))

def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 0.5*height, '%d' % int(height), ha='center', va='bottom')

autolabel(rects)

plt.savefig('barplot.pdf')
plt.show()
