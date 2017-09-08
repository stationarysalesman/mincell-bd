# Simulating ribosomes, probably
# by TC

import hoomd
import hoomd.md
import random
import math
import string
import numpy

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

# Read in a particle configuration
particles = list()
with open('particle_config.txt', 'r') as infile:
    for line in infile:
        items = string.split(line, ',') 
        species, vdwr, x, y, z = items 
        pos = [float(x), float(y), float(z)]
        particles.append([species, vdwr, pos])

# Data dumping
class avg_vel:
    def __init__(self, system):
        self.system = system
    def __call__(self, timestep):
        snap = self.system.take_snapshot();
        avg_vel = numpy.mean(snap.particles.velocity, axis=0);
        print(timestep, ':', avg_vel);

# Create a snapshot
num_particles = len(particles)

snapshot = hoomd.data.make_snapshot(N=num_particles, 
		box=hoomd.data.boxdim(L=box_size, dimensions=3), 
		particle_types=['R','P'])
for i in range(num_particles):
    if (particles[i][0] == 'R'):
        snapshot.particles.typeid[i] = 0
    else:
        snapshot.particles.typeid[i] = 1
    snapshot.particles.position[i] = particles[i][2] 
system = hoomd.init.read_snapshot(snapshot)

# Make a list of cells
nl = hoomd.md.nlist.cell(r_buff=1e-9)

# Lennard-Jones potential parameters
rr_sigma = diam_rib/2 
rr_epsilon = 1.69e-20
rp_sigma = (diam_rib/2 + diam_prot/2)/2
rp_epsilon = 1.69e-20
pp_sigma = diam_prot/2
pp_epsilon = 1.69e-20
rr_cutoff = diam_rib*5
rp_cutoff = diam_rib*5
pp_cutoff = diam_prot*5
lj = hoomd.md.pair.lj(r_cut=rr_cutoff, nlist=nl)
lj.pair_coeff.set('R','R', epsilon=rr_epsilon, sigma=rr_sigma)
lj.pair_coeff.set('R','P', epsilon=rp_epsilon, sigma=rp_sigma, r_cut=rp_cutoff)
lj.pair_coeff.set('P','P', epsilon=pp_epsilon, sigma=pp_sigma, r_cut=pp_cutoff)

# Wall potentials
"""
wall_bottom = hoomd.md.wall.plane(origin=(0, 0, -box_size/2), normal=(0, 0, 1.0), inside=True)
wall_top = hoomd.md.wall.plane(origin=(0, 0, box_size/2), normal=(0, 0, -1.0), inside=True)
wall_negX = hoomd.md.wall.plane(origin=(-box_size/2, 0, 0), normal=(1.0, 0, 0), inside=True)
wall_posX = hoomd.md.wall.plane(origin=(box_size/2, 0, 0), normal=(-1.0, 0, 0), inside=True)
wall_negY = hoomd.md.wall.plane(origin=(0, -box_size/2, 0), normal=(0, 1.0, 0), inside=True)
wall_posY = hoomd.md.wall.plane(origin=(0, box_size/2, 0), normal=(0, -1.0, 0), inside=True)
wall_lst = [wall_bottom, wall_top, wall_negX, wall_negY, wall_posX, wall_posY]
all_walls = hoomd.md.wall.group(wall_lst)
wlj = hoomd.md.wall.lj(all_walls, r_cut=1e-9)
wlj.force_coeff.set('R', sigma=rr_sigma, epsilon=rr_epsilon)
wlj.force_coeff.set('P', sigma=rr_sigma, epsilon=rr_epsilon)
wlj.force_coeff.set(['R','P'], sigma=rp_sigma, epsilon=rp_epsilon)
"""

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
                    quantities=['potential_energy', 'temperature', 'kinetic_energy', 'momentum'],
                    period=1,
                    overwrite=True)
#hoomd.analyze.callback(callback=avg_vel(system), period=1)

hoomd.dump.gsd("trajectory.gsd", period=100, group=all_parts, overwrite=True)

# Run the simulation
hoomd.run(1e3)



