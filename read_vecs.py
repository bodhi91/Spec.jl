import os
import numpy as np
import matplotlib.pyplot as plt

filename = "eigenvecs.dat"

evecs = []

with open(filename, 'r') as f:
    for line in f.readlines():
        evecs.append(line.split(' '))


evec1 = [evec[0] for evec in evecs]
evec2 = [evec[1] for evec in evecs]


plt.scatter(evec1, evec2, c = "blue")
plt.show()
