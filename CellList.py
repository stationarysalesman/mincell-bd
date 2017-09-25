# Simulating ribosomes, probably
# by TC

import hoomd
import hoomd.md
import random
import math
import sys
import itertools
import numpy as np
# Utilities for setting up the simulation

class Particle:
    """Store properties needed to place particles"""

    def __init__(self, species, vdwr, pos = np.zeros(3)):
        self.species = species 
        self.vdwr = vdwr
        self.pos = pos

    def pos(self):
        return self.pos

    def vdwr(self):
        return self.vdwr


class Cell:

    def __init__(self):
        self.particles = list()

    def particles(self):
        return self.particles


class CellList:
 
   
    def __init__(self, edge):
        self.edge = edge
        self.cell_list = dict() 
        self.euc = lambda a1, a2: np.sqrt(np.sum(np.square(np.diff((a1,a2), axis=0))))
      
    def insert(self, particle):
        """Insert a particle into its cell."""
        x,y,z = particle.pos
        xi, yi, zi = int(x/edge), int(y/edge), int(z/edge)
        try:
            self.cell_list[(xi,yi,zi)].append(particle)
        except KeyError:
            self.cell_list[(xi,yi,zi)] = [particle]


    def translate(self, origin):
        """Translate all particles to a new origin."""
        for particle_lst in self.cell_list.values():
            for particle in particle_lst:
                particle.pos += origin

        return


    def export(self):
        """Export particle configuration to a file."""
        with open('particle_config.txt', 'w') as of:
            for particle_lst in self.cell_list.values():
                for particle in particle_lst:
                    of.write(particle.species + ',' + str(particle.vdwr) + ',' + str(particle.pos[0]) + ',' + str(particle.pos[1]) + ',' + str(particle.pos[2]) + '\n')


    def getNeighbors(self, particle):
        """Compute a given particle's neighbor list."""
        x,y,z = particle.pos
        i,j,k = int(x/edge), int(y/edge), int(z/edge)
        nl = list()
        for di,dj,dk in itertools.product([-1,0,1],repeat=3):
            for particle in self.cell_list.get((i+di,j+dj,k+dk), []):
               nl.append(particle)
        return nl

     

box_size = 400e-9 
split = 20
edge = float(box_size / split)
num_particles = 6000 
r_vdwr = 10e-9
p_vdwr = 1e-9

euc = lambda a1, a2: np.sqrt(np.sum(np.square(np.diff((a1,a2), axis=0))))


clist = CellList(edge)

r = 200e-9

def rand_pos():
    r1 = random.uniform(0, box_size)
    r2 = random.uniform(0, box_size)
    r3 = random.uniform(0, box_size)
    return np.array((r1, r2, r3))


def valid(particle):
    """Determine if a given position is within some distance of any 
    other particle."""
    # Find cell of particle
    p_pos = particle.pos
    p_vdwr = particle.vdwr
    nl = clist.getNeighbors(particle)
    for neighbor in nl:
        neighbor_pos = neighbor.pos
        neighbor_vdwr = neighbor.vdwr
        dist = euc(p_pos, neighbor_pos)
        thresh = (p_vdwr + neighbor_vdwr) / 2
        if dist <= thresh:
            return False

    """
    # Too close to wall?
    wall_points = list()
    for i in range(3):
        new_p = list(p_pos)
        new_p[i] = 0.0 
        wall_points.append(new_p)
        new_p = list(p_pos)
        new_p[i] = box_size
        wall_points.append(new_p)

    for p in wall_points:
        dist = euc(p_pos, p)
        if dist <= 5*p_vdwr:
            return False
    """

    # Within the sphere, and far enough away from potential?
    r = 200e-9
    a = np.array(p_pos)
    a -= (r, r, r)
    r_dist = np.sqrt(np.sum(a*a))
    if r_dist <= r * 0.9:
        return True 

    else:
        return False 


# Place the ribosomes
nr = raw_input("Number of ribosomes: ")
for i in range(int(nr)):
    print "Sampling position for ribosome " + str(i) + "."
    p = rand_pos()
    particle = Particle('R', r_vdwr, p) 
    while not valid(particle):
        #print "Resampling..."
        p = rand_pos()
        particle.pos = p 
    clist.insert(particle)


# Place the proteins
_np = raw_input("Number of proteins: ")
for i in range(int(_np)):
    print "Sampling position for particle " + str(i) + "." 
    z = rand_pos()
    particle = Particle('P', p_vdwr, z) 
    while not valid(particle):
        #print "Resampling..."
        z = rand_pos()
        particle.pos = z 
    clist.insert(particle)


origin = [-box_size/2, -box_size/2, -box_size/2]
clist.translate(origin)
clist.export()

# Acceptance test (slow: O(n^2))
#print "Test result: " + str(clist.test())
