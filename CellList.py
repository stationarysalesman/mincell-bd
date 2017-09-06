# Simulating ribosomes, probably
# by TC

import hoomd
import hoomd.md
import random
import math
box_size = 100
split = 10

# Utilities for setting up the simulation

class Cell:

    def __init__(self):
        self.particles = list()

    def particles(self):
        return self.particles

class CellList:
    
    def __init__(self, size=100, n_split=10):
        self.n_split = n_split 
        self.cells = dict() 
        for i in range(n_split):
            self.cells[i] = dict()
            for j in range(n_split):
                self.cells[i][j] = dict()
                for k in range(n_split):
                    self.cells[i][j][k] = Cell()

    def insert(self, p):
        """Insert a particle into its cell."""
        i = int(p[0] / self.n_split)
        j = int(p[1] / self.n_split)
        k = int(p[2] / self.n_split)
        the_cell = self.cells[i][j][k]
        the_cell.particles.append(p)

    def get_neighbors(self, p):
        """Get the cell that contains point p."""
        i = int(p[0] / self.n_split)
        j = int(p[1] / self.n_split)
        k = int(p[2] / self.n_split)
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
                elif k != self.n_split-1:
                    neighbors.append(self.cells[i-1][j-1][k+1])
                    neighbors.append(self.cells[i-1][j][k+1])
                    neighbors.append(self.cells[i][j][k+1])
                    neighbors.append(self.cells[i][j-1][k+1])
            elif j != self.n_split-1:
                neighbors.append(self.cells[i-1][j+1][k])
                neighbors.append(self.cells[i][j+1][k])
                if k != 0:
                    neighbors.append(self.cells[i-1][j+1][k-1])
                    neighbors.append(self.cells[i][j+1][k-1])
                elif k != self.n_split-1:
                    neighbors.append(self.cells[i-1][j+1][k+1])
                    neighbors.append(self.cells[i][j+1][k+1])

        elif i != self.n_split-1:
            neighbors.append(self.cells[i+1][j][k])
            if j != self.n_split-1:
                neighbors.append(self.cells[i+1][j+1][k])
                if k != self.n_split-1:
                    neighbors.append(self.cells[i+1][j+1][k+1])
                    neighbors.append(self.cells[i+1][j][k+1])
                elif k != 0:
                    neighbors.append(self.cells[i+1][j+1][k-1])
                    neighbors.append(self.cells[i+1][j][k-1])
            elif j != 0:
                neighbors.append(self.cells[i+1][j-1][k])
                if k != self.n_split-1:
                    neighbors.append(self.cells[i+1][j-1][k+1])
                elif k != 0:
                    neighbors.append(self.cells[i+1][j-1][k-1])

        neighbor_list = list()
        neighbor_list.append(the_cell.particles)
        for c in neighbors:
            neighbor_list.append(c.particles)

        return neighbor_list

placed_particles = list() # keep track of particles already placed in simulation
neighbor_list = list() # n^2 is bad yo

euc = lambda v1, v2: math.sqrt((v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2)

valid = lambda v: v > rr_sigma

def rand_pos():
    r1 = random.uniform(0, box_size)
    r2 = random.uniform(0, box_size)
    r3 = random.uniform(0, box_size)
    return [r1, r2, r3]

def all_valid(v):

    # Find cell of particle

    for placed in placed_particles:
        if euc(placed, v) < rr_sigma:
            return False
    return True


clist = CellList()
num_particles = 100

for i in range(num_particles):
    print "Sampling position for particle " + str(i) + "."
    p = rand_pos()
    neighbor_list = clist.get_neighbors(p)
    clist.insert(p)
    print neighbor_list 


