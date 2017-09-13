# this is probably not the most efficient way to do this, but it's efficient with programmer time...

positions  = [] # [ [x0,y0,z0], [x1,y1,z1], ... ]

edge = 1.0 # Length of cell edge

cell_list = dict()

# build cell list

for x,y,z in positions:
    xi, yi, zi = int(x/edge), int(y/edge), int(z/edge)
    try:
        cell_list[(xi,yi,zi)].append((x,y,z))
    except KeyError:
        cell_list[(xi,yi,zi)] = [(x,y,z)]


# compute interactions

import itertools

for (i,j,k), cell0 in cell_list.items():
    for di,dj,dk in itertools.product([-1,0,1],repeat=3):
        for (x,y,z) in cell_list.get((i+di,j+dj, k+dk), []):
            for x0,y0,z0 in cell0:
                r = np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
                if 0 < r <= edge:
                    do_something(x,y,z,x0,y0,z0)
