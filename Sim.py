# Simulating ribosomes and proteins, probably
# by TC

import hoomd
import hoomd.md
import random
import math
import numpy as np 
import time  # is there a real life analog for this import statement?
import datetime
import re
import sys
import os
import tmeanalysis
import itertools
import copy

# All length units specified in nanometers

# Simulation params
T = 298 # Temperature in Kelvin
k_b = 1.380648e-23 # Boltzmann constant, in J/K
kT = k_b * T
seed = 42
dt = 1e-6 # Delta t
sim_t = 1e6  # Simulation time steps
use_planar_wall = False # planar wall potentials
use_sphere_wall = False # sphere wall potential
interacting = False # enable potentials
frame_period = 10 # how often to write trajectories
box_size = 1.82e-3 # length of any one side of simulation box


# Protein
m_prot = 346000
diff_prot = 2.18e-10 # protein diffusion constant (1 nm radius)
diam_prot = 2e-9
gamma_p = kT / (diff_prot)

def run():
    """Run the simulation."""

    # Initialize the context
    hoomd.context.initialize()

    # Read in a particle configuration
    particles = list()


    with open('particle_config.txt', 'r') as infile:
        for line in infile:
            items = re.split(',', line) 
            species, vdwr, x, y, z = items 
            pos = [float(x), float(y), float(z)]
            particles.append([species, vdwr, pos])


    # Create a snapshot
    num_particles = len(particles)
    snapshot = hoomd.data.make_snapshot(N=num_particles, 
                    box=hoomd.data.boxdim(L=box_size, dimensions=3), 
                    particle_types=['P'])
    for i in range(num_particles):
        snapshot.particles.typeid[i] =0 
        snapshot.particles.diameter[i] = diam_prot
        snapshot.particles.position[i] = particles[i][2] 
    system = hoomd.init.read_snapshot(snapshot)
    print(system.particles)
    system_particles = system.particles

    # Create the neighbor lists (acceleration structure) 
    p_nl = hoomd.md.nlist.cell(r_buff=1e-10)
   
    # Specify the Lennard-Jones potential
    pp_sigma = diam_prot/2 + (diam_prot/2)*0.2
    pp_epsilon = 1.69e-20 
    pp_cutoff = diam_prot*1.5
    alpha = 1.0
    if interacting:
        # Protein-scale neighbor list 
        p_lj = hoomd.md.pair.lj(r_cut=pp_cutoff, nlist=p_nl)
        p_lj.pair_coeff.set('P','P', epsilon=pp_epsilon, sigma=pp_sigma, 
                            r_cut=pp_cutoff, alpha=alpha)

       
         
    # Wall potentials
    """
    if use_planar_wall:
        wall_bottom = hoomd.md.wall.plane(origin=(0, 0, -box_size/2), 
                                            normal=(0, 0, 1.0), inside=True)
        wall_top = hoomd.md.wall.plane(origin=(0, 0, box_size/2), 
                                        normal=(0, 0, -1.0), inside=True)
        wall_negX = hoomd.md.wall.plane(origin=(-box_size/2, 0, 0), 
                                        normal=(1.0, 0, 0), inside=True)
        wall_posX = hoomd.md.wall.plane(origin=(box_size/2, 0, 0), 
                                        normal=(-1.0, 0, 0), inside=True)
        wall_negY = hoomd.md.wall.plane(origin=(0, -box_size/2, 0), 
                                        normal=(0, 1.0, 0), inside=True)
        wall_posY = hoomd.md.wall.plane(origin=(0, box_size/2, 0), 
                                        normal=(0, -1.0, 0), inside=True)
        wall_lst = [wall_bottom, wall_top, wall_negX, wall_negY, wall_posX, wall_posY]
        all_walls = hoomd.md.wall.group(wall_lst)
    """
 
    # Sphere potential
    """
    elif use_sphere_wall:
        r_sphere = 900e-6
        sphere_potential = hoomd.md.wall.sphere(r_sphere)
        all_walls = hoomd.md.wall.group([sphere_potential])
    """
 
    # Set up the simulation
    hoomd.md.integrate.mode_standard(dt=dt)
    all_parts = hoomd.group.all()
    proteins = hoomd.group.type(name='proteins', type='P')
    bd = hoomd.md.integrate.brownian(group=all_parts, kT=kT, seed=seed) 
    bd.set_gamma('P',gamma_p)

    # Create log directory
    h = abs(hash((datetime.datetime.now())) + hash(seed))
    directory = 'validation/' + str(h) + '/'
    os.mkdir(directory)
    ptraj_fname = directory + 'ptraj' + str(h) + '.gsd'

    if proteins:
        hoomd.dump.gsd(ptraj_fname, period=frame_period, group=proteins, 
                        overwrite=True)

    # dump both for pdf
    hoomd.dump.gsd(directory + 'aggregatetraj' + str(h) + '.gsd', period=frame_period, 
                    group=all_parts, overwrite=True)
    

    # Run the simulation
    t0 = time.clock()
    hoomd.run(sim_t)
    tf = time.clock()
    elapsed = tf - t0

    # Log 
    fname = directory + 'bdsim'+str(h)
    with open(fname, 'w') as logfile:
        logfile.write('SIMULATION PARAMETERS\n\n')
        logfile.write('Simulation id: ' + str(h) + '\n') 
        logfile.write('Simulation datetime: ' + str(datetime.datetime.now()) + '\n') 
        logfile.write('Edge length: ' + str(box_size) + '\n')
        logfile.write('Number of particles: ' + str(num_particles) + '\n')
        logfile.write('Temperature: ' + str(T) + '\n')
        logfile.write('Seed: ' + str(seed) + '\n')
        logfile.write('dt: ' + str(dt) + '\n')
        logfile.write('Timesteps taken: ' + str(sim_t) + '\n')
        logfile.write('Total simulation time: ' + str(sim_t * dt) + '\n')
        logfile.write('Total wall clock time: ' + str(tf-t0) + '\n')
        logfile.write('Protein diffusion coefficient: ' + str(diff_prot) + '\n')
        logfile.write('Protein diameter: ' + str(diam_prot) + '\n')
        if interacting:
            logfile.write('LJ epsilon for PP interaction: ' + str(pp_epsilon) + '\n')
            logfile.write('LJ sigma for PP interaction: ' + str(pp_sigma) + '\n')
            logfile.write('LJ alpha: ' + str(alpha) + '\n')
        else:
            logfile.write('Simulating noninteracting particles.\n')
        if use_sphere_wall:
            logfile.write('Using LJ wall potential (sphere). ' \
                    + 'Coefficients based on particle LJ parameters.\n')
        if use_planar_wall:
            logfile.write('Using LJ wall potential (planar). ' \
                    + 'Coefficients based on particle LJ parameters.\n')
        logfile.write('Writing trajectories every ' + str(frame_period) + ' steps.\n') 
        logfile.write('\n')
    return h

def main():
    run()


main()

