# Simulating ribosomes, probably
# by TC

import hoomd
import hoomd.md
import random

# All length units specified in nanometers

# Simulation params
T = 300 # Temperature in Kelvin
k_b = 1.380648e-23 # Boltzmann constant, in J/K

kT = k_b * T

seed = 42
dt = 1

# For a spherical cell with diameter 400nm, use a lattice 
# with side length 322.39839 (~322) to approximate the same cell volume.
cell_size = 32.2e-9
box_size = 322e-9 # length of any one side of simulation box

# Ribosomes
m_rib = 2700000
diff_rib =  5e-14 # ribosome diffusion constant in nm2/s
diam_rib = 20e-9
gamma_r = kT / diff_rib # Friction coefficient 

# Protein
m_prot = 346000
diff_prot = 0.01e-9 # protein diffusion constant in nm2/s
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

# Create a snapshot
random.seed(seed)

num_particles = 1 
snapshot = hoomd.data.make_snapshot(N=num_particles, 
			box=hoomd.data.boxdim(L=box_size, dimensions=3), 
			particle_types=['R','P'])
for i in range(num_particles):
	r1 = random.uniform(-box_size/2, box_size/2)
	r2 = random.uniform(-box_size/2, box_size/2)
	r3 = random.uniform(-box_size/2, box_size/2)
	snapshot.particles.position[i] = [r1, r2, r3]
	snapshot.particles.velocity[i] = [0.0, 0.0, 0.0]	
	print("particle", i, ":",snapshot.particles.position[i])

hoomd.init.read_snapshot(snapshot)

# Make a list of cells
nl = hoomd.md.nlist.cell(r_buff=1e-9)

# Define potential energy
rr_cutoff = diam_rib
rp_cutoff = diam_rib
pp_cutoff = diam_prot
lj = hoomd.md.pair.lj(r_cut=rr_cutoff, nlist=nl)
lj.pair_coeff.set('R','R', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('R','P', epsilon=1.0, sigma=1.0, r_cut=rp_cutoff)
lj.pair_coeff.set('P','P', epsilon=1.0, sigma=1.0, r_cut=pp_cutoff)

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
hoomd.run(1e5)

