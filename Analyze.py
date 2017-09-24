import numpy as np
# GSD file reader: git clone https://bitbucket.org/glotzer/gsd.git
from gsd.hoomd import HOOMDTrajectory
from gsd.fl import GSDFile
import sys
import re
import matplotlib.pyplot as plt
import phys

def analyze(directory, h, dt, diff_rib, diff_prot, box_size, rho):

    fname = directory+'trajlog-'+str(h) 
    rfname = directory + 'rtraj' + str(h) + '.gsd'
    pfname = directory + 'ptraj' + str(h) + '.gsd'
    logfile = None 
    f,axarr = plt.subplots(2, sharex=True)
    plt.xlabel('time (ns)')
    try:
        logfile = open(fname, 'w')
        logfile.write('Calculations based on trajectories:\n')
    except EnvironmentError:
        print "404 file not found :("
        sys.exit(1)

    ts = []
    rib_msds = []

    try:
        rib_traj = HOOMDTrajectory(GSDFile(rfname, "rb", "HOOMD", "hoomd"))

        rib_pos0 = rib_traj.read_frame(0).particles.position
        frmDataRib = np.zeros((rib_traj.file.nframes, rib_pos0.shape[0], rib_pos0.shape[1]), dtype=rib_pos0.dtype)
        ts = np.zeros(rib_traj.file.nframes)
        for i in range(rib_traj.file.nframes):
            frm = rib_traj.read_frame(i)
            frmDataRib[i,...] = frm.particles.position
            ts[i] = dt*frm.configuration.step
        rib_max_dxs = np.max(np.abs(np.diff(frmDataRib,axis=0)), axis=(0,2))
        frmDataFiltRib = frmDataRib[:,rib_max_dxs < 0.1*box_size,:]
        rib_dx = frmDataFiltRib-frmDataFiltRib[0,:,:]
        rib_msd = np.mean(np.sum(rib_dx*rib_dx, axis=2), axis=1)
        rib_m,rib_b = np.polyfit(ts, rib_msd, 1)

        ymin = np.min(rib_msd)
        ymax = np.max(rib_msd)
        axarr[0].scatter(ts/dt, rib_msd, color='r', s=1)
        axarr[0].set_ylim([ymin, ymax])
        mod_ts = ts[1:]
        mod_msd = rib_msd[1:]
        rib_diffusion = np.divide(mod_msd, (6*mod_ts))
        rib_diffusion = np.mean(rib_diffusion[1:])
        logfile.write("ribosome msd linear fit: m={}, b={}, slope_err={}\n".format(rib_m,rib_b, rib_m/(6*diff_rib)))
        logfile.write("ribosome diffusion coefficient={} (error: {})\n".format(rib_diffusion, rib_diffusion/diff_rib))
    except EnvironmentError:
        logfile.write("No ribosome trajectory.\n")

    try:
        prot_traj = HOOMDTrajectory(GSDFile(pfname, "rb", "HOOMD", "hoomd"))
        prot_pos0 = prot_traj.read_frame(0).particles.position
        prot_msds = []
        frmDataProt = np.zeros((prot_traj.file.nframes, prot_pos0.shape[0], prot_pos0.shape[1]), dtype=prot_pos0.dtype)
        prot_ts = np.zeros(prot_traj.file.nframes)
        ts = np.zeros(prot_traj.file.nframes)
        for i in range(prot_traj.file.nframes):
            frm = prot_traj.read_frame(i)
            frmDataProt[i,...] = frm.particles.position
            ts[i] = dt*frm.configuration.step
        prot_max_dxs = np.max(np.abs(np.diff(frmDataProt,axis=0)), axis=(0,2))
        frmDataFiltProt = frmDataProt[:,prot_max_dxs < 0.1*box_size,:]
        prot_dx = frmDataFiltProt-frmDataFiltProt[0,:,:]
        prot_msd = np.mean(np.sum(prot_dx*prot_dx, axis=2), axis=1)
        prot_m,prot_b = np.polyfit(ts, prot_msd, 1)
        ymin_p = np.min(prot_msd)
        ymax_p = np.max(prot_msd)
        axarr[1].scatter(ts/dt, prot_msd)
        axarr[1].set_ylim([ymin_p, ymax_p])
        mod_ts = ts[1:]
        mod_msd = prot_msd[1:]
        prot_diffusion = np.divide(mod_msd, (6*mod_ts))
        prot_diffusion = np.mean(prot_diffusion[1:]) 
        logfile.write("protein linear fit: m={}, b={}, slope_err={}\n".format(prot_m,prot_b, prot_m/(6*diff_prot)))
        logfile.write("protein diffusion coefficient={} (error: {})\n".format(prot_diffusion, prot_diffusion/diff_prot))
        
    except EnvironmentError:
        logfile.write("No protein trajectory.\n")


    logfile.close()
    plt.show()
    
    # display pdf
    atfname = directory + 'aggregatetraj' + str(h) + '.gsd'
    agg_traj = HOOMDTrajectory(GSDFile(atfname, "rb", "HOOMD", "hoomd"))
    pos0 = agg_traj.read_frame(0).particles.position
    posFinal = np.zeros(pos0.shape[0])
    finalFrame = agg_traj.read_frame(agg_traj.file.nframes-1)
    posFinal = finalFrame.particles.position
    maxb = 200e-9
    dr = 2e-9
    edge = maxb/dr
    phys.pdf(posFinal, 10e-9, 200e-9, rho)


# Analyze a trajectory file



f = sys.argv[1]
x = re.split('/', f)
print x
h=x[-2]
fname = f+'bdsim'+h
d = dict()
with open(fname, 'r') as infile:
    for line in infile:
        s = re.split(': ', line)
        if len(s) > 1:
            d[s[0]] = s[1][:len(s[1])-1]

h = int(d['Simulation id'])
directory = f 
dt = float(d['dt'])
diff_rib = float(d['Ribosome diffusion coefficient'])
diff_prot = float(d['Protein diffusion coefficient'])
box_size = float(d['Edge length'])
rho = float(d['Simulation density (particle volume/simulation volume)'])

analyze(directory, h, dt, diff_rib, diff_prot, box_size, rho)
