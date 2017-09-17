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
import os
import tmeanalysis

# All length units specified in nanometers

# Simulation params
T = 300 # Temperature in Kelvin
k_b = 1.380648e-23 # Boltzmann constant, in J/K
kT = k_b * T
seed = 42
dt = 1e-10
use_walls = True
frame_period = 1 # how often to write trajectories

# For a spherical cell with diameter 400nm, use a lattice 
# with side length 322.39839nm (~322) to approximate the same cell volume.
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

def pdf(position, position_list):
    euc = lambda v1,v2: math.sqrt((v1[0]-v2[0])**2+(v1[1]-v2[1])**2+(v1[2]-v2[2])**2)
    bins = dict()
    for pos in position_list:
        if pos == position:
            continue
        d = int(euc(position, pos) / (box_size/10))
        try:
            bins[d] += 1
        except KeyError:
            bins[d] = 1
    pdf_dict = dict()
    for r, c in bins.items():
        print "r is " + str(r)
        print "c is " + str(c)
        if r == 0:
            v = (4/3)*math.pi
        else:
            v = (4/3)*math.pi * ((r+1) ** 3) - (4/3)*math.pi*(r**3)
        pdf_val = c / v
        pdf_dict[r] = pdf_val
    return pdf_dict



def run():

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
    #rr_epsilon = 1.5e-10 
    rp_sigma = (diam_rib/2 + diam_prot/2)/2
    rp_epsilon = 1.69e-20 
    #rp_epsilon = 1.5e-10 
    pp_sigma = diam_prot/2
    pp_epsilon = 1.69e-20 
    #pp_epsilon = 1.5e-10 
    rr_cutoff = diam_rib*5
    rp_cutoff = diam_rib*5
    pp_cutoff = diam_prot*5
    lj = hoomd.md.pair.lj(r_cut=rr_cutoff, nlist=nl)
    alpha = 1.0
    lj.pair_coeff.set('R','R', epsilon=rr_epsilon, sigma=rr_sigma, alpha=alpha)
    lj.pair_coeff.set('R','P', epsilon=rp_epsilon, sigma=rp_sigma, r_cut=rp_cutoff, alpha=alpha)
    lj.pair_coeff.set('P','P', epsilon=pp_epsilon, sigma=pp_sigma, r_cut=pp_cutoff, alpha=alpha)


    # Wall potentials
    if use_walls:
        wall_bottom = hoomd.md.wall.plane(origin=(0, 0, -box_size/2), normal=(0, 0, 1.0), inside=True)
        wall_top = hoomd.md.wall.plane(origin=(0, 0, box_size/2), normal=(0, 0, -1.0), inside=True)
        wall_negX = hoomd.md.wall.plane(origin=(-box_size/2, 0, 0), normal=(1.0, 0, 0), inside=True)
        wall_posX = hoomd.md.wall.plane(origin=(box_size/2, 0, 0), normal=(-1.0, 0, 0), inside=True)
        wall_negY = hoomd.md.wall.plane(origin=(0, -box_size/2, 0), normal=(0, 1.0, 0), inside=True)
        wall_posY = hoomd.md.wall.plane(origin=(0, box_size/2, 0), normal=(0, -1.0, 0), inside=True)
        wall_lst = [wall_bottom, wall_top, wall_negX, wall_negY, wall_posX, wall_posY]
        all_walls = hoomd.md.wall.group(wall_lst)
        wlj = hoomd.md.wall.lj(all_walls, r_cut=diam_rib*5)
        wlj.force_coeff.set('R', sigma=rr_sigma, epsilon=rr_epsilon, alpha=alpha)
        wlj.force_coeff.set('P', sigma=rr_sigma, epsilon=rr_epsilon, alpha=alpha)
        wlj.force_coeff.set(['R','P'], sigma=rp_sigma, epsilon=rp_epsilon, alpha=alpha)


    # Set up the simulation
    hoomd.md.integrate.mode_standard(dt=dt)
    all_parts = hoomd.group.all()
    ribosomes = hoomd.group.type(name='ribosomes', type='R')
    proteins = hoomd.group.type(name='proteins', type='P')
    print len(ribosomes)
    print len(proteins)
    bd = hoomd.md.integrate.brownian(group=all_parts, kT=kT, seed=seed) 
    bd.set_gamma('R',gamma_r)
    bd.set_gamma('P',gamma_p)

    # Create log directory
    h = abs(hash((datetime.datetime.now())) + hash(seed))
    directory = 'validation/' + str(h) + '/'
    os.mkdir(directory)
    rtraj_fname = directory + 'rtraj' + str(h) + '.gsd'
    ptraj_fname = directory + 'ptraj' + str(h) + '.gsd'
    if ribosomes:
        hoomd.dump.gsd(rtraj_fname, period=frame_period, group=ribosomes, overwrite=True)
    if proteins:
        hoomd.dump.gsd(ptraj_fname, period=frame_period, group=proteins, overwrite=True)


    # Run the simulation
    sim_t = 1e3
    t0 = time.clock()
    for i in range(int(sim_t)): 
        hoomd.run(1)
        if i % 100 == 0:
            p_pos = []
            for j in range(num_particles):
                p_pos.append(system.particles[j].position)
            for p in p_pos:
                d = pdf(p, p_pos)
                print "pdf of particle at " + str(p) + ": " + str(d)

    tf = time.clock()
    elapsed = tf - t0

    # Log 

    fname = directory + 'bdsim'+str(h)
    with open(fname, 'w') as logfile:
        logfile.write('SIMULATION PARAMETERS\n\n')
        logfile.write('Simulation id: ' + str(h) + '\n') 
        logfile.write('Simulation datetime: ' + str(datetime.datetime.now()) + '\n') 
        logfile.write('Simulation volume: ' + str((box_size)**3) + '\n')
        logfile.write('Number of particles: ' + str(num_particles) + '\n')
        logfile.write('Temperature: ' + str(T) + '\n')
        logfile.write('Seed: ' + str(seed) + '\n')
        logfile.write('dt: ' + str(dt) + '\n')
        logfile.write('Timesteps taken: ' + str(sim_t) + '\n')
        logfile.write('Total simulation time: ' + str(sim_t * dt) + '\n')
        logfile.write('Total wall clock time: ' + str(tf-t0) + '\n')
        logfile.write('Ribosome diffusion coefficient: ' + str(diff_rib) + '\n')
        logfile.write('Ribosome diameter: ' + str(diam_rib) + '\n')
        logfile.write('Protein diffusion coefficient: ' + str(diff_prot) + '\n')
        logfile.write('Protein diameter: ' + str(diam_prot) + '\n')
        logfile.write('LJ epsilon for RR interaction: ' + str(rr_epsilon) + '\n')
        logfile.write('LJ sigma for RR interaction: ' + str(rr_sigma) + '\n')
        logfile.write('LJ epsilon for RP interaction: ' + str(rp_epsilon) + '\n')
        logfile.write('LJ sigma for RP interaction: ' + str(rp_sigma) + '\n')
        logfile.write('LJ epsilon for PP interaction: ' + str(pp_epsilon) + '\n')
        logfile.write('LJ sigma for PP interaction: ' + str(pp_sigma) + '\n')
        logfile.write('LJ alpha: ' + str(alpha) + '\n')
        if use_walls:
            logfile.write('Using LJ wall potential. Coefficients based on particle LJ parameters.\n')
        logfile.write('Writing trajectories every ' + str(frame_period) + ' steps.\n') 
        logfile.write('\n')
    return h

def main():

    # Run the simulation
    h = run()
    
    # Analyze the simulation
    tmeanalysis.analyze('validation/' + str(h) + '/', h, dt, diff_rib, diff_prot, box_size)    

main()
