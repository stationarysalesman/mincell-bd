# Simulating ribosomes, probably
# by TC

import hoomd
import hoomd.md

# All length units specified in nanometers

# Simulation params
_T = 300 # Temperature in Kelvin
#_k_b = 1.380648e-23 # Boltzmann constant, in J/K
_k_b = 0.01

_kT = _k_b * _T

_seed = 42
_dt = 1


# Lattice parameters
# For a spherical cell with diameter 400nm, use a lattice 
# with side length 322.39839 (~322) to approximate the same cell volume.
cell_size = 32.2

# Particle parameters 
# Units:
#   mass: Da
#   diameter: nm

# Ribosomes
m_rib = 2700000
diff_rib =  5e-5 # ribosome diffusion constant in nm2/s
diam_rib = 20
_gamma_r = _kT / diff_rib # Friction coefficient 

# Protein
m_prot = 346000
diff_prot = 0.01 # protein diffusion constant in nm2/s
diam_prot = 2
_gamma_p = _kT / diff_prot

# Initialize the context
hoomd.context.initialize("")

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

# Make a list of cells
nl = hoomd.md.nlist.cell()

# Define potential energy
rr_cutoff = 2 * diam_rib
rp_cutoff = 2 * diam_rib
pp_cutoff = 2 * diam_prot
lj = hoomd.md.pair.lj(r_cut=rr_cutoff, nlist=nl)
lj.pair_coeff.set('R','R', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('R','P', epsilon=1.0, sigma=1.0, r_cut=rp_cutoff)
lj.pair_coeff.set('P','P', epsilon=1.0, sigma=1.0, r_cut=pp_cutoff)
"""
rp_lj = hoomd.md.pair.lj(r_cut=1.5*diam_rib, nlist=nl)
rp_lj.pair_coeff.set('R','P', epsilon=1.0, sigma=1.0)
pp_lj = hoomd.md.pair.lj(r_cut=2*diam_prot, nlist=nl)
pp_lj.pair_coeff.set('P','P', epsilon=1.0, sigma=1.0)
"""

# Set up the simulation
hoomd.md.integrate.mode_standard(dt=_dt)


# Brownian Dynamics 
all_parts = hoomd.group.all()

# Not overdamped
#hoomd.md.integrate.langevin(group=all_parts, kT=_kT, seed=_seed) 

# Overdamped
bd = hoomd.md.integrate.brownian(group=all_parts, kT=_kT, seed=_seed) 
bd.set_gamma('R',_gamma_r)
bd.set_gamma('P',_gamma_p)


# Log stuff
hoomd.analyze.log(filename="ribosome-protein-output.log",
                    quantities=['potential_energy', 'temperature'],
                    period=100,
                    overwrite=True)
hoomd.dump.gsd("trajectory.gsd", period=100, group=all_parts, overwrite=True)

# Run the simulation
hoomd.run(1e5)

