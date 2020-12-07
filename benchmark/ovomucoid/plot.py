#!/bin/python3

# HEADER
################################################################################

import matplotlib.pyplot as plt

# I put this function explicitly here so there are no dependencies.
def Col(fileName, col, start = 0, stop = 0):
    data = []
    
    try:
        for x, y in enumerate(open(fileName)):
            if start == 0 and stop == 0 and (y.split()[0][0] not in ['@','#','&']):
                data.append(float(y.split()[col-1]))

            elif (x >= start-1) and (x <= stop-1):
                data.append(float(y.split()[col-1]))
    
    except IndexError:
        pass

    return data

# ACTUAL CODE
################################################################################

lambdaGroupNum = Col('results.txt', 1)
performaceMean = Col('results.txt', 2)

performaceMean = [val / performaceMean[0] for val in performaceMean]

plt.scatter(lambdaGroupNum, performaceMean)
plt.plot(lambdaGroupNum, performaceMean)

plt.axis([0, 16, 0, 1])
plt.xlabel('number of titratable groups')
plt.ylabel('relative running speed')

plt.savefig('benchmark.pdf')
plt.show()
