import random
import numpy as np
import copy
import random
import matplotlib.pyplot as plt

def binsort(positions, dr):
    p_arr = np.array(positions)
    c_array = copy.copy(p_arr)
    hist = dict()
    distances = []
    for p in p_arr:
        while len(c_array) > 1:
            c = c_array[1]
            if np.array_equal(p,c):
                continue
            px,py,pz = p
            cx,cy,cz = c
            d = np.sqrt((cx-px)**2+(cy-py)**2+(cz-pz)**2)
            distances.append(d) 
            b = int(d/dr)+1
            try:
                hist[b] += 2
            except KeyError:
                hist[b] = 2
            c_array = c_array[1:]
    print "min dist is : " + str(min(distances))
    return hist

def binnorm(hist, rho, dr, nparticles):
    const = (4.0 * np.pi * rho) / 3.0
    gr = dict()
    print "normalizing keys : " + str(hist.keys())
    for b in hist:
        rlower = (b-1)*dr
        rupper = rlower+dr
        nideal = const * (rupper**3 - rlower**3)
        gr[b] = hist[b] / float(nparticles) / nideal 
    if 0 not in gr:
        gr[0] = 0.0
    return gr


def pdf(positions, dr, edge, rho, rmin_rr, rmin_rp, rmin_pp):
    hist = dict()
    hist = binsort(positions, dr)
    gr = binnorm(hist, rho, dr, len(positions)) 
    x = []
    y = []
    for k in sorted(gr.iterkeys()):
        x.append(k)
        y.append(gr[k])
 
    x = np.array(x)
    y = np.array(y)
    x = x * dr * 1e9 # nanometers
    plt.scatter(x,y, s=10, color='red')
    plt.xlabel('radius (nm)')
    plt.ylabel('G(r)')
    z = np.polyfit(x,y,12)
    p = np.poly1d(z)
    plt.plot(x, p(x),"r--", color='gray')

    # Sigmas
    line1 = plt.axvline(rmin_rr*1e9, color='red', label='RR $r_m$')
    """ 
    line4 = plt.axvline(rmin_rr*2*1e9, color='green')
    line5 = plt.axvline(rmin_rr*3*1e9, color='green')
    line6 = plt.axvline(rmin_rr*4*1e9, color='green')
    line7 = plt.axvline(rmin_rr*5*1e9, color='green')
    """ 
    #line2 = plt.axvline(rmin_rp, color='green')
    line3 = plt.axvline(rmin_pp*1e9, color='green', label='PP $r_m$')
    plt.legend(handles=[line1,line3])
    plt.show()



