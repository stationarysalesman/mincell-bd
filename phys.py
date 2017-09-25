import random
import numpy as np
import copy
import random
import matplotlib.pyplot as plt

def binsort(positions, dr):
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
            b = int(d/dr)+1
            try:
                hist[b] += 2
            except KeyError:
                hist[b] = 2
            c_array = c_array[1:]
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
    return gr


def pdf(positions, dr, edge, rho, sigma_rr, sigma_rp, sigma_pp):
    hist = dict()
    hist = binsort(positions, dr)
    hist = binnorm(hist, rho, dr, len(positions)) 
    x = []
    y = []
    for k in sorted(hist.iterkeys()):
        x.append(k)
        y.append(hist[k])
  
    x = np.array(x)
    y = np.array(y)
    x = x * dr
    print dr
    plt.scatter(x,y, s=10, color='red')
    plt.xlabel('radius')
    plt.ylabel('G(r)')
    plt.xlim((0, max(x)))
    
    z = np.polyfit(x,y,12)
    p = np.poly1d(z)
    plt.plot(x, p(x),"r--", color='gray')

    # Sigmas
    line1 = plt.axvline(sigma_rr, color='green', label='Sigma')
    line4 = plt.axvline(sigma_rr*2, color='green')
    line5 = plt.axvline(sigma_rr*3, color='green')
    line6 = plt.axvline(sigma_rr*4, color='green')
    line7 = plt.axvline(sigma_rr*5, color='green')
    #line2 = plt.axvline(sigma_rp, color='green')
    #line3 = plt.axvline(sigma_pp, color='green')
    plt.legend(handles=[line1])
    plt.show()



