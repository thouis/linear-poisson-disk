import sys
import random
import time
from math import sqrt, ceil, floor, pi, cos

# This code requires Polygon 
# http://www.cs.man.ac.uk/~toby/alan/software/
# http://pypi.python.org/pypi/Polygon
from Polygon import Polygon
from Polygon.Shapes import Circle
import sample
import time

# radius of the disks to be added.  Points (aka centers) will be spaced by 2R.
R = 0.01

start_time = time.time()

if len(sys.argv) > 1:
    try:
        R = float(sys.argv[1])
    except:
        sys.stderr.write('Usage: %s <disk radius>\n'%(sys.argv[0]))
        sys.exit(1)

# Some constants
Rsquared = R**2
fourRsquared = 4 * Rsquared
grid_size = int(ceil(sqrt(2.0) / R))
grid_spacing = 1 / float(grid_size)

neighbor_offset = int(ceil(2 * R / grid_spacing))

inf = float('inf')
# This test is only available in python 3.0 and up
def isinf(v):
    return v == inf


# The grid is stored as a dictionary from (i, j) index of the grid square.
# The meaning of values in the grid:
#   no entry - a grid square that has not been requested yet
#   (p, t, r) - point p, time of arrival t, valid region r
#      if t == infinity and r is True, then the point was accepted
#      if t == infinity and r is False, then that grid square was
#        covered by other points


grid = {}
bucket = set()

def unwrap((i, j)):
    # wrap around
    i = (i + grid_size) % grid_size
    j = (j + grid_size) % grid_size
    return (i, j)
    

def fetch_grid((i, j)):
    '''on demand generation of empty grid squares.'''
    (i, j) = unwrap((i, j))
    if (i, j) not in grid:
        xloc = i * grid_spacing
        yloc = j * grid_spacing
        r = Polygon([(xloc, yloc),
                     (xloc + grid_spacing, yloc),
                     (xloc + grid_spacing, yloc + grid_spacing),
                     (xloc, yloc + grid_spacing)])
        p = sample.sample(r)
        t = random.expovariate(r.area())
        grid[i, j] = (p, t, r)
    return grid[(i, j)]

def sq_distance_point_to_grid((px, py), (i, j)):
    xlo = i * grid_spacing
    xhi = (i + 1) * grid_spacing
    ylo = j * grid_spacing
    yhi = (j + 1) * grid_spacing
    xdiff = min(abs(px - xlo), abs(px - xhi)) if not xlo < px < xhi else 0
    ydiff = min(abs(py - ylo), abs(py - yhi)) if not ylo < py < yhi else 0
    return xdiff**2 + ydiff**2

def neighbors((i, j)):
    '''the neighbor grid squares of the point at (i, j) the grid.'''
    p = fetch_grid((i, j))[0]
    nlist = []
    for offset_i in range(-neighbor_offset, neighbor_offset + 1):
        for offset_j in range(-neighbor_offset, neighbor_offset + 1):
            if offset_i == 0 and offset_j == 0:
                continue

            # this check purposefully relies on sq_distance_point_to_grid() not calling unwrap()
            if sq_distance_point_to_grid(p, (i + offset_i, j + offset_j)) > fourRsquared:
                continue

            nlist.append((i + offset_i, j + offset_j))
    return nlist


def disk((x, y)):
    ''' the exclusion disk, twice the radius'''
    
    # because we're using an approximate disk, we use a just slightly
    # larger R so that the minimum internal radius of the polygon is R.
    Radjust = R / cos(pi / 32)
    return Circle(2 * Radjust, (x, y), 32)


def too_close((x, y), (z, w)):
    return ((x - z)**2 + (y - w)**2) <= fourRsquared

def closest_p(p, n):
    # find closest p to grid n in wrapped space
    closest_d = sq_distance_point_to_grid(p, n)
    pclose = p
    for off_x in [-1, 0, 1]:
        for off_y in [-1, 0, 1]:
            off_dist = sq_distance_point_to_grid((p[0] + off_x, p[1] + off_y), n)
            if off_dist < closest_d:
                closest_d = off_dist
                pclose = (p[0] + off_x, p[1] + off_y)
    return pclose
    

def update_valid(n, p):
    n = unwrap(n)
    # this grid square must have been created already in locally_early()
    q, t, r = grid[n]

    # no change if this square is already output or covered
    if isinf(t):
        return False

    p = closest_p(p, n)

    r = r - disk(p)
    r.simplify()
    # update
    grid[n] = (q, t, r)

    # Did p invalidate q?
    if too_close(q, p):
        A = r.area()
        if A == 0:
            # q's square has been covered
            grid[n] = (q, inf, False)
        else:
            # generate a new q
            q = sample.sample(r)
            t = t + random.expovariate(A)
            grid[n] = (q, t, r)

        return True
    return False


def check_nearby((i, j)):
    ''' find all nearby points that have become locally early.'''
    # These loops have to be over all grid squares that *might* have a
    # point that has (i,j) as a neighbor.
    for offset_i in range(-neighbor_offset, neighbor_offset + 1):
        for offset_j in range(-neighbor_offset, neighbor_offset + 1):
            p, t, r = fetch_grid((i + offset_i, j + offset_j))
            
            p = closest_p(p, (i, j))
            # this check purposefully relies on sq_distance_point_to_grid() not calling unwrap()
            if sq_distance_point_to_grid(p, (i, j)) > fourRsquared:
                continue
                
            if locally_early((i + offset_i, j + offset_j)):
                bucket.add((i + offset_i, j + offset_j))


def accept((i, j)):
    (i, j) = unwrap((i, j))
    p, t, r = grid[(i, j)]

    # update all valid regions for our neighbors
    updated_neighbors = [n for n in neighbors((i, j)) if update_valid(n, p)]

    # mark p for i,j as accepted - setting region to True indicates
    # that this point should be output.
    grid[(i, j)] = (p, inf, True)

    # These checks must be after p's time is set to inf.
    check_nearby((i, j))
    for n in updated_neighbors:
        check_nearby(n)


def locally_early((i, j)):
    (i, j) = unwrap((i, j))

    if (i, j) in bucket:
        # already marked
        return True

    p, t, r = fetch_grid((i, j))

    if isinf(t):
        # already in output or grid square is covered
        return False

    for n in neighbors((i, j)):
        np, nt, nr = fetch_grid(n)
        if nt < t:
            return False

    return True



def traverse():
    '''Traverse the grid, looking for points to add.  When one is
    found, empty the bucket before continuing.'''
    for i in range(grid_size):
        for j in range(grid_size):
            if locally_early((i, j)):
                bucket.add((i, j))
            while len(bucket) > 0:
                accept(bucket.pop())

def distchk(p, q):
    p = closest_p(p, (floor(q[0] / grid_spacing), floor(q[1] / grid_spacing)))
    d =  (p[0] - q[0])**2 + (p[1] - q[1])**2
    assert d > fourRsquared
    return "%f %f   sqrt %f %f"%(d, fourRsquared, sqrt(d), 2 * R)



# for reproducibility & debugging
# random.seed(0)

traverse()

output = set([p for p, t, r in grid.values() if isinf(t) and r == True])

print R, len(output), sample.sample_count, float(sample.sample_count) / len(output), grid_size * grid_size, len(output) / (time.time() - start_time)

# Every grid square should be completely covered
for p, t, r in grid.values():
    assert r == False or r == True

# validate
for p in output:
    for q in output:
        if p != q:
            distchk(p, q)
