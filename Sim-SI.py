# Simulating ribosomes, probably
# by TC

import hoomd
import hoomd.md
import random
import math
import string
import numpy as np
import time  # is there a real life analog for this import statement?
import datetime
import re
import sys

# All length units specified in nanometers

# Simulation params
T = 300 # Temperature in Kelvin
k_b = 1.380648e-23 # Boltzmann constant, in J/K
kT = k_b * T
seed = 42
dt = 1e-15

# For a spherical cell with diameter 400nm, use a lattice 
# with side length 322.39839nm (~322) to approximate the same cell volume.
#cell_size = 32.2e-9
box_size = 322e-9 # length of any one side of simulation box

# Ribosomes
m_rib = 2700000
diff_rib =  5e-14 # ribosome diffusion constant
diam_rib = 20e-9
gamma_r = kT / diff_rib # Friction coefficient 
#gamma_r = 8.28e-08

# Protein
m_prot = 346000
diff_prot = 10e-12 # protein diffusion constant
diam_prot = 2e-9
gamma_p = kT / diff_prot
#gamma_p = 4.14e-10

# Initialize the context
hoomd.context.initialize("")

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
        avg_vel = np.mean(snap.particles.velocity, axis=0);
        print('avg',timestep, ':', avg_vel);

class max_vel:
    def __init__(self, system):
        self.system = system
    def __call__(self, timestep):
        snap = self.system.take_snapshot();
        max_vel = np.amax(snap.particles.velocity, axis=0);
        print('max',timestep, ':', max_vel);


# Create a snapshot
num_particles = len(particles)

snapshot = hoomd.data.make_snapshot(N=num_particles, 
		box=hoomd.data.boxdim(L=box_size, dimensions=3), 
		particle_types=['R','P'])

for i in range(num_particles):
    if (particles[i][0] == 'R'):
        snapshot.particles.typeid[i] = 0
        snapshot.particles.diameter[i] = diam_rib
    else:
        snapshot.particles.typeid[i] = 1
        snapshot.particles.diameter[i] = diam_prot

    snapshot.particles.position[i] = particles[i][2] 

system = hoomd.init.read_snapshot(snapshot)
system_particles = system.particles

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
lj.pair_coeff.set('R','R', epsilon=rr_epsilon, sigma=rr_sigma, alpha=0)
lj.pair_coeff.set('R','P', epsilon=rp_epsilon, sigma=rp_sigma, r_cut=rp_cutoff, alpha=0)
lj.pair_coeff.set('P','P', epsilon=pp_epsilon, sigma=pp_sigma, r_cut=pp_cutoff, alpha=0)


# Wall potentials

wall_bottom = hoomd.md.wall.plane(origin=(0, 0, -box_size/2), normal=(0, 0, 1.0), inside=True)
wall_top = hoomd.md.wall.plane(origin=(0, 0, box_size/2), normal=(0, 0, -1.0), inside=True)
wall_negX = hoomd.md.wall.plane(origin=(-box_size/2, 0, 0), normal=(1.0, 0, 0), inside=True)
wall_posX = hoomd.md.wall.plane(origin=(box_size/2, 0, 0), normal=(-1.0, 0, 0), inside=True)
wall_negY = hoomd.md.wall.plane(origin=(0, -box_size/2, 0), normal=(0, 1.0, 0), inside=True)
wall_posY = hoomd.md.wall.plane(origin=(0, box_size/2, 0), normal=(0, -1.0, 0), inside=True)
wall_lst = [wall_bottom, wall_top, wall_negX, wall_negY, wall_posX, wall_posY]
all_walls = hoomd.md.wall.group(wall_lst)
wlj = hoomd.md.wall.lj(all_walls, r_cut=diam_rib*5)
wlj.force_coeff.set('R', sigma=rr_sigma, epsilon=rr_epsilon)
wlj.force_coeff.set('P', sigma=rr_sigma, epsilon=rr_epsilon)
wlj.force_coeff.set(['R','P'], sigma=rp_sigma, epsilon=rp_epsilon)


# Set up the simulation
hoomd.md.integrate.mode_standard(dt=dt)
all_parts = hoomd.group.all()
bd = hoomd.md.integrate.brownian(group=all_parts, kT=kT, seed=seed) 
bd.set_gamma('R',gamma_r)
bd.set_gamma('P',gamma_p)

# Log stuff
hoomd.analyze.log(filename="ribosome-protein-output.log",
                    quantities=['potential_energy', 'temperature', 'kinetic_energy', 'momentum'],
                    period=1,
                    overwrite=True)
#hoomd.analyze.callback(callback=max_vel(system), period=1)
#hoomd.analyze.callback(callback=avg_vel(system), period=1)

hoomd.dump.gsd("trajectory.gsd", period=100, group=all_parts, overwrite=True)


sd = lambda v1, v2: np.array([(v1[0]-v2[0])**2, (v1[1]-v2[1])**2, (v1[2]-v2[2])**2])
initial_pos_rib = list() 
initial_pos_prot = list()
for i in system.particles:
    if i.typeid == 0:
        initial_pos_rib.append(i.position) 
    else:
        initial_pos_prot.append(i.position) 


# Run the simulation
sim_t = 1e3

t0 = time.clock()

prev_positions = [i.position for i in system.particles]
"""
pos = np.array(prev_positions)
ddr = np.sum((pos[None,:,:] - pos[None,:,:].transpose(1,0,2))**2, axis=2)
print(np.sqrt(np.min(ddr[np.tril_indices_from(ddr,-1)])))
exit()
max_dxs = np.max(np.linalg.norm(np.abs(np.diff((new_positions, prev_positions), axis=0)), axis=2))
"""
max_dxs = 0.0
for i in range(int(sim_t)):
    if not i % 100:
        print "Step " + str(i) + "..."
    try:
        hoomd.run(tsteps=1, quiet=True)
    except RuntimeError:
        print "max displacement: " + str(max_dxs) 
        sys.exit(1)
    new_positions = [i.position for i in system.particles]
    max_dxs = np.max(np.linalg.norm(np.abs(np.diff((new_positions, prev_positions), axis=0)), axis=2))
    sys.stdout.write('\rmax_dxs: {:e}'.format(max_dxs))
#    sys.stdout.flush()
#    prev_positions = new_positions

tf = time.clock()
elapsed = tf - t0

final_pos_rib = list()
final_pos_prot = list()
for i in system.particles:
    if i.typeid == 0:
        final_pos_rib.append(i.position)
    else:
        final_pos_prot.append(i.position)

# Simulated/experimental MSD
rib_sum = np.array([0.0, 0.0, 0.0]) 
for v1, v2 in zip(initial_pos_rib, final_pos_rib):
    rib_sum += sd(v2, v1)

if len(initial_pos_rib) != 0:
    rib_msd = rib_sum / len(initial_pos_rib) 
else:
    rib_msd = [0.0, 0.0, 0.0]

rib_msd_scalar = math.sqrt(rib_msd[0]**2 + rib_msd[1]**2 + rib_msd[2]**2)

prot_sum = np.array([0.0, 0.0, 0.0]) 
for v1, v2 in zip(initial_pos_prot, final_pos_prot):
    prot_sum += sd(v2, v1)

if len(initial_pos_prot) != 0:
    prot_msd = prot_sum / len(initial_pos_prot) 
else:
    prot_msd = [0.0, 0.0, 0.0]

prot_msd_scalar = math.sqrt(prot_msd[0]**2 + prot_msd[1]**2 + prot_msd[2]**2)


# Modeled MSD (via Brownian motion equations)
# Note: should be different depending on whether we simulated
# for shorter/longer than relaxation time
rib_model_msd = 6 * diff_rib * (sim_t * dt)
prot_model_msd = 6 * diff_prot * (sim_t * dt)


# Log 
fname = 'sim-'+str(datetime.datetime.now())+'.log'
fname = re.sub(r" ", "", fname)
with open(fname, 'w') as logfile:
    logfile.write('SIMULATION PARAMETERS\n\n')
    logfile.write('Simulation volume: ' + str((box_size)**3) + '\n')
    logfile.write('Number of particles: ' + str(num_particles) + '\n')
    logfile.write('Temperature: ' + str(T) + '\n')
    logfile.write('Seed: ' + str(seed) + '\n')
    logfile.write('dt: ' + str(dt) + '\n')
    logfile.write('Timesteps taken: ' + str(sim_t) + '\n')
    logfile.write('Total simulation time: ' + str(sim_t * dt) + '\n')
    logfile.write('Total wall clock time: ' + str(tf-t0) + '\n')
    logfile.write('Ribosome vibrational coefficient: ' + str(diff_rib) + '\n')
    logfile.write('Ribosome diameter: ' + str(diam_rib) + '\n')
    logfile.write('Protein vibrational coefficient: ' + str(diff_prot) + '\n')
    logfile.write('Protein diameter: ' + str(diam_rib) + '\n')
    logfile.write('LJ epsilon for RR interaction: ' + str(rr_epsilon) + '\n')
    logfile.write('LJ sigma for RR interaction: ' + str(rr_sigma) + '\n')
    logfile.write('LJ epsilon for RP interaction: ' + str(rp_epsilon) + '\n')
    logfile.write('LJ sigma for RP interaction: ' + str(rp_sigma) + '\n')
    logfile.write('LJ epsilon for PP interaction: ' + str(pp_epsilon) + '\n')
    logfile.write('LJ sigma for PP interaction: ' + str(pp_sigma) + '\n')
    logfile.write('\n\n')
    logfile.write('RESULTS\n\n')
    logfile.write("Ribosome MSD: " + str(rib_msd_scalar) + '\n')
    logfile.write("Model prediction of ribosome MSD: " + str(rib_model_msd) + '\n') 
    logfile.write("Protein MSD: " + str(prot_msd_scalar) + '\n')
    logfile.write("Model prediction of ribosome MSD: " + str(prot_model_msd) + '\n') 

