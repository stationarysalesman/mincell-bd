# Simulating ribosomes, probably
# by TC

import hoomd
import hoomd.md
import random
import math
import sys

# Utilities for setting up the simulation

class Particle:
    """Store properties needed to place particles"""

    def __init__(self, vdwr, pos = [0.0, 0.0, 0.0]):
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
 
    euc = lambda v1, v2: math.sqrt((v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2)
    in_range = lambda v: v >= 0 and v <= n_split
    valid_idx_lst = lambda l: in_range(l[0]) and in_range(l[1]) and in_range(l[2])

   
    def __init__(self, size=100, n_split=10):
        self.n_split = n_split
        self.cells = dict() 
        self.in_range = lambda v: v >= 0 and v <= n_split
        self.valid_idx_lst = lambda l: self.in_range(l[0]) and self.in_range(l[1]) and self.in_range(l[2])
         
        for i in range(n_split):
            self.cells[i] = dict()
            for j in range(n_split):
                self.cells[i][j] = dict()
                for k in range(n_split):
                    self.cells[i][j][k] = Cell()


    def insert(self, particle):
        """Insert a particle into its cell."""
        p = particle.pos
        i = int(p[0] / (box_size / self.n_split))
        j = int(p[1] / (box_size / self.n_split))
        k = int(p[2] / (box_size / self.n_split))
        the_cell = self.cells[i][j][k]
        the_cell.particles.append(particle)


    def dump(self):
        """Dump all particles into a list. Useful, probably"""
        master_list = list()
        for cx in self.cells:
            for cy in self.cells[cx]:
                for cz in self.cells[cx][cy]:
                    for particle in self.cells[cx][cy][cz].particles:
                        master_list.append(particle)
        return master_list 


    def get_neighbors(self, particle):
        """Get the cell that holds a particle, and its neighbors."""
        p = particle.pos
        i = int(p[0] / (box_size / self.n_split))
        j = int(p[1] / (box_size / self.n_split))
        k = int(p[2] / (box_size / self.n_split))
        the_cell = self.cells[i][j][k]
        
        # Won't you be my neighbor?
        neighbors = list()
        if i != 0: 
            neighbors.append(self.cells[i-1][j][k]) 
            if j != 0:
                neighbors.append(self.cells[i-1][j-1][k])
                neighbors.append(self.cells[i][j-1][k])
                if k != 0:
                    neighbors.append(self.cells[i-1][j-1][k-1])
                    neighbors.append(self.cells[i-1][j][k-1])
                    neighbors.append(self.cells[i][j][k-1])
                    neighbors.append(self.cells[i][j-1][k-1])
                if k != self.n_split-1:
                    neighbors.append(self.cells[i-1][j-1][k+1])
                    neighbors.append(self.cells[i-1][j][k+1])
                    neighbors.append(self.cells[i][j][k+1])
                    neighbors.append(self.cells[i][j-1][k+1])
            if j != self.n_split-1:
                neighbors.append(self.cells[i-1][j+1][k])
                neighbors.append(self.cells[i][j+1][k])
                if k != 0:
                    neighbors.append(self.cells[i-1][j+1][k-1])
                    neighbors.append(self.cells[i][j+1][k-1])
                if k != self.n_split-1:
                    neighbors.append(self.cells[i-1][j+1][k+1])
                    neighbors.append(self.cells[i][j+1][k+1])

        if i != self.n_split-1:
            neighbors.append(self.cells[i+1][j][k])
            if j != self.n_split-1:
                neighbors.append(self.cells[i+1][j+1][k])
                if k != self.n_split-1:
                    neighbors.append(self.cells[i+1][j+1][k+1])
                    neighbors.append(self.cells[i+1][j][k+1])
                if k != 0:
                    neighbors.append(self.cells[i+1][j+1][k-1])
                    neighbors.append(self.cells[i+1][j][k-1])
            if j != 0:
                neighbors.append(self.cells[i+1][j-1][k])
                if k != self.n_split-1:
                    neighbors.append(self.cells[i+1][j-1][k+1])
                if k != 0:
                    neighbors.append(self.cells[i+1][j-1][k-1])

        neighbor_list = list()
        for particle in the_cell.particles:
            neighbor_list.append(particle)
        for c in neighbors:
            for particle in c.particles:
                neighbor_list.append(particle)
        return neighbor_list
   
    def computeNeighborCells(self, particle):
        """Determine which cell holds a particle, and the cells 
        holding its neighbors."""
        p = particle.pos
        i = int(p[0] / (box_size / self.n_split))
        j = int(p[1] / (box_size / self.n_split))
        k = int(p[2] / (box_size / self.n_split))
        the_cell = self.cells[i][j][k]

        # Won't you be my neighbor?
        idx_list = list() 
        idx_list.append((i, j, k))
        """ 
        if i != 0: 
            idx_list.append((i-1,j,k)) 
            if j != 0:
                idx_list.append((i-1,j-1,k))
                idx_list.append((i,j-1,k))
                if k != 0:
                    idx_list.append((i-1,j-1,k-1))
                    idx_list.append((i-1,j,k-1))
                    idx_list.append((i,j,k-1))
                    idx_list.append((i,j-1,k-1))
                if k != self.n_split-1:
                    idx_list.append((i-1,j-1,k+1))
                    idx_list.append((i-1,j,k+1))
                    idx_list.append((i,j,k+1))
                    idx_list.append((i,j-1,k+1))
            if j != self.n_split-1:
                idx_list.append((i-1,j+1,k))
                idx_list.append((i,j+1,k))
                if k != 0:
                    idx_list.append((i-1,j+1,k-1))
                    idx_list.append((i,j+1,k-1))
                if k != self.n_split-1:
                    idx_list.append((i-1,j+1,k+1))
                    idx_list.append((i,j+1,k+1))

        if i != self.n_split-1:
            idx_list.append((i+1,j,k))
            if j != self.n_split-1:
                idx_list.append((i+1,j+1,k))
                if k != self.n_split-1:
                    idx_list.append((i+1,j+1,k+1))
                    idx_list.append((i+1,j,k+1))
                if k != 0:
                    idx_list.append((i+1,j+1,k-1))
                    idx_list.append((i+1,j,k-1))
            if j != 0:
                idx_list.append((i+1,j-1,k))
                if k != self.n_split-1:
                    idx_list.append((i+1,j-1,k+1))
                if k != 0:
                    idx_list.append((i+1,j-1,k-1))
        """
        for x in range(i-1, i+2):
            for y in range(j-1, j+2):
                for z in range(k-1, k+2):
                    idx_list.append((x,y,z))
        s = set(idx_list)
        idx_list = list(s)
        return filter(self.valid_idx_lst, idx_list)
    




    def test(self):
        print "Performing acceptance test..."
        master_list = self.dump()
        N = len(master_list)
        i = 0
        while not len(master_list) == 0:
            obj = master_list.pop()
            p_pos = obj.pos
            p_vdwr = obj.vdwr
            for neighbor in master_list:
                neighbor_pos = neighbor.pos
                neighbor_vdwr = neighbor.vdwr
                dist = euc(p_pos, neighbor_pos)
                thresh = (p_vdwr + neighbor_vdwr) / 2
                if dist < thresh:
                    print "Error (particle " + str(i) + "): acceptable distance: " + str(thresh) + ", actual distance: " + str(dist)
                    return False 
        return True 

placed_particles = list() # keep track of particles already placed in simulation
neighbor_list = list() # n^2 is bad yo
box_size = 100 
split = 10 


clist = CellList(box_size, split)
rr_sigma = 10e-9 

euc = lambda v1, v2: math.sqrt((v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2)

def rand_pos():
    r1 = random.uniform(0, box_size)
    r2 = random.uniform(0, box_size)
    r3 = random.uniform(0, box_size)
    return [r1, r2, r3]

def valid(particle):
    """Determine if a given position is within some distance of any 
    other particle."""
    # Find cell of particle
    p_pos = particle.pos
    p_vdwr = particle.vdwr
    nl = clist.get_neighbors(particle)
    for neighbor in nl:
        neighbor_pos = neighbor.pos
        neighbor_vdwr = neighbor.vdwr
        dist = euc(p_pos, neighbor_pos)
        thresh = (p_vdwr + neighbor_vdwr) / 2
        if dist < thresh:
            return False
    return True


num_particles = 6000 
r_vdwr = 10e-9
p_vdwr = 1e-9

p = Particle(1, [0,0,0])
lst = clist.computeNeighborCells(p)
print lst
print "len is: " + str(len(lst))
"""
# Place the ribosomes
for i in range(700):
    print "Sampling position for ribosome " + str(i) + "."
    p = rand_pos()
    particle = Particle(r_vdwr, p) 
    while not valid(particle):
        print "Resampling..."
        p = rand_pos()
        particle.pos = p 
    clist.insert(particle)
"""
"""
# Place the proteins
for i in range(num_particles):
    print "Sampling position for particle " + str(i) + "." 
    p = rand_pos()
    particle = Particle(p_vdwr, p) 
    while not valid(particle):
        print "Resampling..."
        p = rand_pos()
        particle.pos = p 
    clist.insert(particle)
"""
print "Test result: " + str(clist.test())
