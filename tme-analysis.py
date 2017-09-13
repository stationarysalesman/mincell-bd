import numpy as np
# GSD file reader: git clone https://bitbucket.org/glotzer/gsd.git
from gsd.hoomd import HOOMDTrajectory
from gsd.fl import GSDFile

traj = HOOMDTrajectory(GSDFile("trajectory.gsd", "rb", "HOOMD", "hoomd"))

pos0 = traj.read_frame(0).particles.position
dt = 1e-9
diff_rib =  5e-14 # ribosome diffusion constant
box_size = 322e-9 # length of any one side of simulation box

ts = []
msds = []

frmData = np.zeros((traj.file.nframes, pos0.shape[0], pos0.shape[1]), dtype=pos0.dtype)
ts = np.zeros(traj.file.nframes)


for i in range(traj.file.nframes):
    frm = traj.read_frame(i)
    frmData[i,...] = frm.particles.position
    ts[i] = dt*frm.configuration.step

max_dxs = np.max(np.abs(np.diff(frmData,axis=0)), axis=(0,2))

frmDataFilt = frmData[:,max_dxs < 0.1*box_size,:]

dx = frmDataFilt-frmDataFilt[0,:,:]
msd = np.mean(np.sum(dx*dx, axis=2), axis=1)
m,b = np.polyfit(ts, msd, 1)

print("linear fit: m={}, b={}, slope_err={}".format(m,b, m/(6*diff_rib)))
