# Simulating ribosomes, probably
# by TC

import hoomd
import hoomd.md
import random
import math

# All length units specified in nanometers

# Simulation params
T = 300 # Temperature in Kelvin
k_b = 1.380648e-23 # Boltzmann constant, in J/K

kT = k_b * T

seed = 42
dt = 1e-9

# For a spherical cell with diameter 400nm, use a lattice 
# with side length 322.39839nm (~322) to approximate the same cell volume.
cell_size = 32.2e-9
box_size = 322e-9 # length of any one side of simulation box

# Ribosomes
m_rib = 2700000
diff_rib =  5e-14 # ribosome diffusion constant
diam_rib = 20e-9
gamma_r = kT / diff_rib # Friction coefficient 

# Protein
m_prot = 346000
diff_prot = 10e-12 # protein diffusion constant
diam_prot = 2e-9
gamma_p = kT / diff_prot

# Lennard-Jones potential parameters
rr_sigma = diam_rib / 2
rr_epsilon = 1.69e-20
rp_sigma = (diam_rib + diam_prot) / 2 
rp_epsilon = 1.69e-20
pp_sigma = diam_prot / 2
pp_epsilon = 1.69e-20

# Initialize the context
hoomd.context.initialize("")


"""
# Create a lattice of 20nm cuboids, 8000 cells total
sim_cell = hoomd.lattice.unitcell(N=2,
                                    a1=[cell_size,0,0],
                                    a2=[0,cell_size,0],
                                    a3=[0,0,cell_size],
                                    dimensions=3,
                                    position=[[0, 0, 0], [cell_size/2, cell_size/2, cell_size/2]],
                                    type_name=['R','P'],
                                    mass=[m_rib, m_prot],
                                    diameter=[diam_rib, diam_prot])

hoomd.init.create_lattice(unitcell=sim_cell, n=10)

"""

# Utilities for setting up the simulation

class Cell:

    def __init__(self):
        self.particles = list()

    def particles(self):
        return self.particles

class CellList:
    n_split = 10
    split = box_size / n_split 
    
    def __init__(self):
        self.cells = dict() 
        for i in range(n_split):
            for j in range(n_split):
                for k in range(n_split):
                    cells[i][j][k] = Cell()

    def get_neighbors(p):
        """Get the cell that contains point p."""
        i = p[0] / split
        j = p[1] / split
        k = p[2] / split
        Cell the_cell = self.cells[i][j][k]

        # won't you be my neighbor
        neighbors = list()


        if i != 0: 
            neighbors.append(cells[i-1][j][k]) 
            if j != 0:
                neighbors.append(cells[i-1][j-1][k])
                neighbors.append(cells[i][j-1][k]
                if k != 0:
                    neighbors.append(cells[i-1][j-1][k-1])
                    neighbors.append(cells[i-1][j][k-1])
                    neighbors.append(cells[i][j][k-1])
                    neighbors.append(cells[i][j-1][k-1])
                elif k != n_split-1:
                    neighbors.append(cells[i-1][j-1][k+1])
                    neighbors.append(cells[i-1][j][k+1])
                    neighbors.append(cells[i][j][k+1])
                    neighbors.append(cells[i][j-1][k+1])
            elif j != n_split-1:
                neighbors.append(cells[i-1][j+1][k])
                neighbors.append(cells[i][j+1][k])
                if k != 0:
                    neighbors.append(cells[i-1][j+1][k-1])
                    neighbors.append(cells[i][j+1][k-1])
                elif k != n_split-1:
                    neighbors.append(cells[i-1][j+1][k+1])
                    neighbors.append(cells[i][j+1][k+1])

        elif i != n_split-1:
            neighbors.append(cells[i+1][j][k])
            if j != n_split-1:
                neighbors.append(cells[i+1][j+1][k])
                if k != n_split-1:
                    neighbors.append(cells[i+1][j+1][k+1])
                    neighbors.append(cells[i+1][j][k+1])
                elif k != 0:
                    neighbors.append(cells[i+1][j+1][k-1])
                    neighbors.append(cells[i+1][j][k-1])
            elif j != 0:
                neighbors.append(cells[i+1][j-1][k])
                if k != n_split-1:
                    neighbors.append(cells[i+1][j-1][k+1])
                elif k != 0:
                    neighbors.append(cells[i+1][j-1][k-1])

        return the_cell
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



def compute_neighbors(p):


for i in range(num_particles):
    print "Sampling position for particle " + str(i) + "."
    p = rand_pos()
    while not all_valid(p):
        p = rand_pos()
    placed_particles.append(p)
    compute_neighbors(p)


# Create a snapshot
random.seed(seed)
num_particles = 60700 


snapshot = hoomd.data.make_snapshot(N=num_particles, 
		box=hoomd.data.boxdim(L=box_size, dimensions=3), 
		particle_types=['R','P'])

for i in range(num_particles): 
    snapshot.particles.position[i] = placed_particles[i] 
    if i < 700:
        snapshot.particles.typeid[i] = 0 
    else:  
        snapshot.particles.typeid[i] = 1 
    snapshot.particles.velocity[i] = [0.0, 0.0, 0.0]	


hoomd.init.read_snapshot(snapshot)

# Make a list of cells
nl = hoomd.md.nlist.cell(r_buff=1e-9)

# Define potential energy
rr_cutoff = diam_rib
rp_cutoff = diam_rib
pp_cutoff = diam_prot
lj = hoomd.md.pair.lj(r_cut=rr_cutoff, nlist=nl)
lj.pair_coeff.set('R','R', epsilon=rr_epsilon, sigma=rr_sigma)
lj.pair_coeff.set('R','P', epsilon=rp_epsilon, sigma=rp_sigma, r_cut=rp_cutoff)
lj.pair_coeff.set('P','P', epsilon=pp_epsilon, sigma=pp_sigma, r_cut=pp_cutoff)

# Set up the simulation
hoomd.md.integrate.mode_standard(dt=dt)

# Brownian Dynamics 
all_parts = hoomd.group.all()

# Not overdamped
#hoomd.md.integrate.langevin(group=all_parts, kT=kT, seed=seed) 

# Overdamped
bd = hoomd.md.integrate.brownian(group=all_parts, kT=kT, seed=seed) 
bd.set_gamma('R',gamma_r)
bd.set_gamma('P',gamma_p)


# Log stuff
hoomd.analyze.log(filename="ribosome-protein-output.log",
                    quantities=['potential_energy', 'temperature'],
                    period=100,
                    overwrite=True)
hoomd.dump.gsd("trajectory.gsd", period=100, group=all_parts, overwrite=True)

# Run the simulation
hoomd.run(1e3)



