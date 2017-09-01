# Simulating ribosomes, probably
# by TC

import hoomd
import hoomd.md

# Simulation params
_kT = 0.2 
_seed = 42
_dt = 1e-5
_gamma = 1 

# Lattice parameters
cell_size = 25 # nm

# Particle parameters 
# Units:
#   mass: Da
#   diameter: nm

# Ribosomes
r_mass = 2700000
r_d = 20.8

# Initialize the context
hoomd.context.initialize("")

# Create a lattice of 20nm cuboids, 8000 cells total
sim_cell = hoomd.lattice.unitcell(N=1,
                                    a1=[cell_size,0,0],
                                    a2=[0,cell_size,0],
                                    a3=[0,0,cell_size],
                                    dimensions=3,
                                    position=[[0, 0, 0]],
                                    type_name=['R'],
                                    mass=[r_mass],
                                    diameter=[r_d])

hoomd.init.create_lattice(unitcell=sim_cell, n=20)

# Make a list of cells
nl = hoomd.md.nlist.cell()

# Define potential energy
r_lj = hoomd.md.pair.lj(r_cut=2*r_d, nlist=nl)
r_lj.pair_coeff.set('R','R', epsilon=1.0, sigma=1.0)

# Set up the simulation
hoomd.md.integrate.mode_standard(dt=_dt)


# Brownian Dynamics 
all_parts = hoomd.group.all()

# Not overdamped
#hoomd.md.integrate.langevin(group=all_parts, kT=_kT, seed=_seed) 

# Overdamped
bd = hoomd.md.integrate.brownian(group=all_parts, kT=_kT, seed=_seed) 
bd.set_gamma('R',_gamma)


# Log stuff
hoomd.analyze.log(filename="ribosome-output.log",
                    quantities=['potential_energy', 'temperature'],
                    period=100,
                    overwrite=True)
hoomd.dump.gsd("trajectory.gsd", period=100, group=all_parts, overwrite=True)

# Run the simulation
hoomd.run(1e4)

