import numpy as np
# GSD file reader: git clone https://bitbucket.org/glotzer/gsd.git
from gsd.hoomd import HOOMDTrajectory
from gsd.fl import GSDFile

dt = 1e-9
diff_rib =  5e-14 # ribosome diffusion constant
diff_prot =  10e-12 # protein diffusion constant
box_size = 322e-9 # length of any one side of simulation box


ts = []
rib_msds = []
try:
    rib_traj = HOOMDTrajectory(GSDFile("ribosome-trajectory.gsd", "rb", "HOOMD", "hoomd"))
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
    print("ribosome msd linear fit: m={}, b={}, slope_err={}".format(rib_m,rib_b, rib_m/(6*diff_rib)))
    mod_ts = ts[1:]
    mod_msd = rib_msd[1:]
    rib_diffusion = np.divide(mod_msd, (6*mod_ts))
    rib_diffusion = np.mean(rib_diffusion[1:])
    print "ribosome diffusion coefficient={} (error: {})".format(rib_diffusion, rib_diffusion/diff_rib)
except IndexError:
    print "No ribosome trajectory."

try:
    prot_traj = HOOMDTrajectory(GSDFile("protein-trajectory.gsd", "rb", "HOOMD", "hoomd"))
    prot_pos0 = prot_traj.read_frame(0).particles.position
    prot_msds = []
    frmDataProt = np.zeros((prot_traj.file.nframes, prot_pos0.shape[0], prot_pos0.shape[1]), dtype=prot_pos0.dtype)
    prot_ts = np.zeros(prot_traj.file.nframes)
    for i in range(prot_traj.file.nframes):
        frm = prot_traj.read_frame(i)
        frmDataProt[i,...] = frm.particles.position
    prot_max_dxs = np.max(np.abs(np.diff(frmDataProt,axis=0)), axis=(0,2))
    frmDataFiltProt = frmDataProt[:,prot_max_dxs < 0.1*box_size,:]
    prot_dx = frmDataFiltProt-frmDataFiltProt[0,:,:]
    prot_msd = np.mean(np.sum(prot_dx*prot_dx, axis=2), axis=1)
    prot_m,prot_b = np.polyfit(ts, prot_msd, 1)
    print("protein linear fit: m={}, b={}, slope_err={}".format(prot_m,prot_b, prot_m/(6*diff_prot)))
    mod_ts = ts[1:]
    mod_msd = prot_msd[1:]
    prot_diffusion = np.divide(mod_msd, (6*mod_ts))
    prot_diffusion = np.mean(prot_diffusion[1:]) 
    print "protein diffusion coefficient={} (error: {})".format(prot_diffusion, prot_diffusion/diff_prot)
except IndexError:
    print "No protein trajectory."


