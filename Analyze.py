import numpy as np
# GSD file reader: git clone https://bitbucket.org/glotzer/gsd.git
from gsd.hoomd import HOOMDTrajectory
from gsd.fl import GSDFile
import sys
import re
import matplotlib.pyplot as plt
import phys

def analyze(directory, h, dt, diff_rib, diff_prot, box_size, rho, sigma_rr, sigma_rp, sigma_pp):

    fname = directory+'trajlog-'+str(h) 
    rfname = directory + 'rtraj' + str(h) + '.gsd'
    pfname = directory + 'ptraj' + str(h) + '.gsd'
    logfile = None 
    rib_traj = None
    prot_traj = None
      
    f,axarr = plt.subplots(2)
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
        xmin = 0.0
        scaled_x = ts*1e9
        xmax = np.max(scaled_x)
        axarr[0].set_ylim([ymin, ymax])
        axarr[0].set_xlim([xmin, xmax])
        axarr[0].set_ylabel('Ribosome $<r^2> (m^2)$')
        axarr[0].scatter(scaled_x, rib_msd, color='blue', s=2)
        rib_trendline, = axarr[0].plot(scaled_x, rib_m*(ts)+rib_b, color='red', label='Ribosome MSD trendline')
        plt.legend(handles=[rib_trendline])
        mod_ts = ts[1:]
        mod_msd = rib_msd[1:]
        rib_diffusion = np.divide(mod_msd, (6*mod_ts))
        rib_diffusion = np.mean(rib_diffusion[1:])
        logfile.write("ribosome msd linear fit: m={}, b={}, slope_err={}\n".format(rib_m,rib_b, rib_m/(6*diff_rib)))
        funthing = rib_msd / 6
        d_m, d_b = np.polyfit(ts, funthing, 1)
        logfile.write("ribosome Diff linear fit: m={}, b={}, slope_err={}\n".format(d_m,d_b, d_m/(5e-14)))

    
    except EnvironmentError:
        logfile.write("No ribosome trajectory.\n")
        f.delaxes(axarr[0])
    
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
        scaled_x = ts * 1e9
        axarr[1].set_ylim([ymin_p, ymax_p])
        axarr[1].set_ylabel('Protein $<r^2> (m^2)$')
        xmax = np.max(scaled_x)
        xmin = 0.0
        axarr[1].set_xlim([xmin, xmax])
        axarr[1].scatter(scaled_x, prot_msd, s=2, color='blue')
        prot_trendline = axarr[1].plot(scaled_x, prot_m*(ts)+prot_b, color='red', label='Protein MSD trendline')
        mod_ts = ts[1:]
        mod_msd = prot_msd[1:]
        prot_diffusion = np.divide(mod_msd, (6*mod_ts))
        prot_diffusion = np.mean(prot_diffusion[1:]) 
        logfile.write("protein linear fit: m={}, b={}, slope_err={}\n".format(prot_m,prot_b, prot_m/(6*diff_prot)))
        logfile.write("protein diffusion coefficient={} (error: {})\n".format(prot_diffusion, prot_diffusion/diff_prot))
        
    except EnvironmentError:
        logfile.write("No protein trajectory.\n")
        print type(axarr[1])
        f.delaxes(axarr[1])
    

    logfile.close()
    plt.xlabel('time (ns)')
    plt.show()
      
    # display pdf
        
    if (rib_traj):
        pos0 = rib_traj.read_frame(0).particles.position
        posFinal = np.zeros(pos0.shape[0])
        finalFrame = rib_traj.read_frame(rib_traj.file.nframes-1)
        posFinal = finalFrame.particles.position
        maxb = 200e-9
        dr = sigma_rr * 5 
        edge = maxb/dr
        phys.pdf(posFinal, dr, 200e-9, rho, sigma_rr, sigma_rp, sigma_pp)

    if (prot_traj):
        pos0 = prot_traj.read_frame(0).particles.position
        posFinal = np.zeros(pos0.shape[0])
        finalFrame = prot_traj.read_frame(prot_traj.file.nframes-1)
        posFinal = finalFrame.particles.position
        maxb = 200e-9
        dr = sigma_pp * 5 
        edge = maxb/dr
        phys.pdf(posFinal, dr, 200e-9, rho, sigma_rr, sigma_rp, sigma_pp)


    atfname = directory + 'aggregatetraj' + str(h) + '.gsd'
    agg_traj = HOOMDTrajectory(GSDFile(atfname, "rb", "HOOMD", "hoomd"))
    pos0 = agg_traj.read_frame(0).particles.position
    posFinal = np.zeros(pos0.shape[0])
    finalFrame = agg_traj.read_frame(agg_traj.file.nframes-1)
    posFinal = finalFrame.particles.position
    maxb = 200e-9
    dr = sigma_pp * 5 
    edge = maxb/dr
    phys.pdf(posFinal, dr, 200e-9, rho, sigma_rr, sigma_rp, sigma_pp)



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
rmin_rr = float(d['LJ sigma for RR interaction']) * pow(2, 1/6.)
rmin_rp = float(d['LJ sigma for RP interaction']) * pow(2, 1/6.)
rmin_pp = float(d['LJ sigma for PP interaction']) * pow(2, 1/6.)
analyze(directory, h, dt, diff_rib, diff_prot, box_size, rho, rmin_rr, rmin_rp, rmin_pp)
