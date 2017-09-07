# Simulating ribosomes, probably
# by TC

import hoomd
import hoomd.md
import random
import math

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
    
    def __init__(self, size=100, n_split=10):
        self.n_split = n_split
        self.split = size / n_split
        self.cells = dict() 
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

    def get_neighbors(self, particle):
        """Get the cell that holds a particle, and its neighbors."""
        p = particle.pos
        i = int(p[0] / (box_size / self.n_split))
        j = int(p[1] / (box_size / self.n_split))
        k = int(p[2] / (box_size / self.n_split))
        the_cell = self.cells[i][j][k]
         
        # won't you be my neighbor
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

placed_particles = list() # keep track of particles already placed in simulation
neighbor_list = list() # n^2 is bad yo
box_size = 322e-09 
split = 60 


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
        if ((p_vdwr == neighbor_vdwr) and (dist < neighbor_vdwr)):
            return False
        else:
            r_avg = (p_vdwr + neighbor_vdwr)/2
            if dist < r_avg:
                return False
    return True


num_particles = 60000 
r_vdwr = 10e-9
p_vdwr = 1e-9

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

for i in range(700):
    print "Sampling position for ribosome " + str(i) + "."
    p = rand_pos()
    particle = Particle(r_vdwr, p) 
    while not valid(particle):
        print "Resampling..."
        p = rand_pos()
        particle.pos = p 
    clist.insert(particle)



