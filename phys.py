import random
import numpy as np
import copy
import random
import matplotlib.pyplot as plt

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

def binnorm(hist, rho, dr):
    const = 4.0 * np.pi * rho / 3.0
    gr = dict()
    for b in hist:
        rlower = (b-1)*dr
        rupper = rlower+dr
        nideal = const * (rupper**3 - rlower**3)
        gr[b+.5*dr] = hist[b] / nideal
    return gr

def main():
    dr = 5
    positions = []
    box_size = 100
    maxbin = (box_size/2) / dr
    random.seed(42)
    for i in range(3000):
        r1 = random.uniform(0, box_size)
        r2 = random.uniform(0, box_size)
        r3 = random.uniform(0, box_size)
        positions.append([r1, r2, r3])
    hist = binsort(positions, dr, maxbin)
    rho = 1.0
    hist = binnorm(hist, rho, dr)
    
    x = []
    y = []
    for k in sorted(hist.iterkeys()):
        x.append(k)
        y.append(hist[k])
    plt.plot(x,y)
    plt.show()


main()
