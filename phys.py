import random
import numpy as np
import copy
import random
import matplotlib.pyplot as plt

def fun(pos):
    Nparticles = len(pos)
    sigma = .1 * 1e-9
    ideal_density = (1.0 * 14.0) / (0.0821 * 300)
    # computes dist_matrix[i,j] = |r_i - r_j|
    dist_matrix = np.sqrt(np.sum((pos[None, ...] - pos[None, ...].transpose((1,0,2)))**2, axis=2))

    # gives just the pairs below the main diagonal:
    #e.g. [(i,j) for i in  range(N) for j in range(0,i)]
    drs = dist_matrix[np.tril_indices_from(dist_matrix, -1)]

    bins = np.linspace(0, 3*sigma, 20)
    vols = np.diff(4/3*np.pi*bins**3)
    counts, _ = np.histogram(drs, bins=bins)
    gs = counts/(ideal_density*vols)/Nparticles
    rs = 0.5*(bins[:-1]+bins[1:])
    print gs
    print rs

def binsort(positions, dr, maxbin):
    p_arr = np.array(positions)
    c_array = copy.copy(p_arr)
    hist = dict()
    for p in p_arr:
        while len(c_array) > 1:
            c = c_array[1]
            if np.array_equal(p,c):
                continue
            px,py,pz = p
            cx,cy,cz = c
            d = np.sqrt((cx-px)**2+(cy-py)**2+(cz-pz)**2)
            b = int(d/dr) + 1
            if b <= maxbin:
                try:
                    hist[b] += 2
                except KeyError:
                    hist[b] = 2
            c_array = c_array[1:]
    return hist

def binnorm(hist, rho, dr, nparticles):
    const = 4.0 * np.pi * rho / 3.0
    gr = dict()
    for b in hist:
        rlower = (b-1)*dr
        rupper = rlower+dr
        nideal = const * (rupper**3 - rlower**3)
        gr[b+.5*dr] = hist[b] / nideal
#        gr[b] = hist[b] / sigma / nparticles / nideal 
    return gr

def pdf(positions, dr, edge):
    maxbin = int(edge / dr)
    rho = (1.0 * 14.0) / (8.314 * 300)
    hist = dict()
    hist = binsort(positions, dr, maxbin)
    hist = binnorm(hist, rho, dr, len(positions)) 
    x = []
    y = []
    for k in sorted(hist.iterkeys()):
        x.append(k)
        y.append(hist[k])
    plt.scatter(x,y)
    plt.show()



