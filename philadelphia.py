# Copyright 2021 D-Wave Systems Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import math
import numpy as np

# Ignore errors importing matplotlib.pyplot (may not be available in
# testing framework)
try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.cm as cm
    from matplotlib.patches import RegularPolygon
except ImportError:
    pass

P1_REUSE_DISTANCES = [2*np.sqrt(3), np.sqrt(3), 1, 1, 1, 0]
P2_REUSE_DISTANCES = [np.sqrt(7), np.sqrt(3), 1, 1, 1, 0]

# Note: the original problem description states that for P2, P4, and
# P6, the first reuse distance (d0) is sqrt(7).  However, this does
# not appear to be consistent with the given diagrams (e.g., Figure 4
# at http://fap.zib.de/problems/Philadelphia/).  The diagrams indicate
# that node pairs such as (3, 20) must be separated by at least 1
# frequency.  However, the distance between such nodes is exactly
# sqrt(7), and the description states that "no interference occurs if
# and only if the centers have mutual distance >= d".  This is also
# consistent with the mathematical description given by, e.g., Koster
# (1999, pp. 28-29).
#
# In attempt to be consistent with the given diagrams, we perturb the
# first reuse distance for P2.
P2_REUSE_DISTANCES[0] += 1e-3


# Axial coordinates of the nodes/stations/cells, see
# https://www.redblobgames.com/grids/hexagons

# We take cell 1 to have coordinates q=0, r=0
qvals = [0, 1, 2, 3, 4,
         -2, -1, 0, 1, 2, 3, 4,
         -3, -2, -1, 0, 1, 2,
         -1, 0, 1]
rvals = 5*[0] + 7*[1] + 6*[2] + 3*[3]


def _euclidean_dist(cell1, cell2):
    """Get Euclidean distance between cells.

    Args:
        cell1 (int):
            Index for first cell (starting with 1)
        cell2 (int):
            Index for second cell (starting with 1)
    """
    if cell1 == cell2:
        return 0
    dq = qvals[cell2 - 1] - qvals[cell1 - 1]
    dr = rvals[cell2 - 1] - rvals[cell1 - 1]
    return np.sqrt(dq**2 + dr**2 + dq*dr)


def _ge_close(x, y):
    """Utility routine for >= with approximate equality.

    In order to avoid floating point equality comparison issues, just
    to be safe
    """
    return x > y or np.isclose(x, y, rtol=1e-8)


def get_forbidden_set(cell1, cell2, reuse_distances):
    """Get set of forbidden frequency differences for given cell numbers.

    Args:
        cell1 (int):
            Index for first cell (starting with 1)
        cell2 (int):
            Index for second cell (starting with 1)
        reuse_distances (list):
            Problem-specific list of frequency re-use distances
    
    Returns:
        set
    """
    d = _euclidean_dist(cell1, cell2)

    for j in range(1, len(reuse_distances)+1):
        if _ge_close(d, reuse_distances[j-1]):
            continue
        if _ge_close(d, reuse_distances[j]):
            return set(range(j))
    return set()


def load_problem(problem='small'):
    """Load network structure for problem definition.

    Args:
        problem (str)
    
    Returns:
        demand (dict):
            Dictionary mapping node numbers to capacities
        reuse_distances (list):
            List of reuse distances
        nfreq (int):
            Default number of frequencies to consider
    """

    # Start with default demand, which applies to P1, P2, and other
    # "small" problem versions
    demand = (8,25,8,8,8,15,18,52,77,28,13,15,31,15,36,57,28,8,10,13,8)
    if problem == 'P3' or problem == 'P4':
        demand = (5,5,5,8,12,25,30,25,30,40,40,45,20,30,25,15,15,30,20,20,25)
    elif problem == 'P5' or problem == 'P6':
        demand = (20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20)
    elif problem == 'P7':
        demand = (16,50,16,16,16,30,36,104,154,56,26,30,62,30,72,114,56,16,20,26,16)
    elif problem == 'P8':
        demand = (8,25,8,8,8,15,18,52,77,28,13,15,31,15,36,57,28,8,10,13,8)
    elif problem == 'P9':
        demand = (32,100,32,32,32,60,72,208,308,112,52,60,124,60,144,228,112,32,40,52,32)


    # Dictionary mapping node number to demand:
    demand = {i+1: demand[i] for i in range(len(demand))}

    reuse_distances = P1_REUSE_DISTANCES
    if problem == 'P2' or problem == 'P4' or problem == 'P6':
        reuse_distances = P2_REUSE_DISTANCES

    nfreq = None  # None is flag to use default based on sum(demand)

    if problem == 'trivial':
        # Further reduce problem size by reducing the demand at each node:
        demand = {1: 2}
        nfreq = 6
        
    elif problem == 'single':
        # Further reduce problem size by reducing the demand at each node:
        demand = {1: 4}
        nfreq = 16

    elif problem == 'small':
        # Simplify problem by removing nodes:
        keep_nodes = [1,2,7,8,9,15,16]

        # Further reduce problem size by reducing the demand at each node:
        demand = {i: math.ceil(demand[i] / 5) for i in keep_nodes}

    elif problem == 'very-small':
        # Simplify problem by removing nodes:
        keep_nodes = [4,5,11,12]

        # Further reduce problem size by reducing the demand at each node:
        demand = {i: math.ceil(demand[i] / 5) for i in keep_nodes}

    # Define default nfreq values for the Philadelphia benchmark
    # problems.  Start by defining opt_range, which is the reported
    # optimal value of fmax-fmin.  In theory, the problem should have
    # a feasible solution using nfreq=opt_range+1.  By default, we
    # include additional frequencies in the search.

    # opt_ranges[i] gives the range for problem Pi
    opt_ranges = (-1, 426, 426, 257, 252, 239, 179, 855, 524, 1713)
    if problem.startswith('P'):
        opt_range = opt_ranges[int(problem[1:])]
        nfreq = opt_range + 20

    if nfreq is None:
        nfreq = sum(demand.values()) + 50

    return demand, nfreq, reuse_distances


def plot_nodes(nodes, demand=None, interference_nodes=[], demand_violations=[]):
    """Plot hexagonal grid for the specified nodes.

    Create a plot of the specified nodes, which must be a subset from
    the Philadelphia instance.  The pre-defined axial coordinates are
    used to place the nodes.

    The nodes that are identified as having either interference or
    demand violations are marked using hatches.  Currently not
    distinction is made in the display of interference vs demand
    violations, but the two are separated to facilitate modifications
    that treat them differently.

    Args:
        nodes (list):
            Subset of node indices from the Philadelphia instance
        demand (dict or None):
            Mapping of node numbers to capacities.  If provided, color
            the nodes by demand.
        interference_nodes (list):
            List of nodes that have interference
        demand_violations (list):
            List of nodes where the demand violation constraint is not
            satisfied
    """

    cmap = plt.get_cmap('YlOrRd')
    # Make the colormap a little bit lighter:
    alpha = 0.3

    f, ax = plt.subplots()
    ax.set_aspect('equal')

    if demand is not None:
        # Create normalizer for mapping demand values to colors.  Set
        # vmax above the largest demand to exclude the upper portion
        # of the colormap.
        norm = mcolors.Normalize(vmin=min(demand.values()), vmax=max(demand.values())*1.1)

    for n in nodes:
        q, r = qvals[n-1], rvals[n-1]
        # Convert axial coordinates to Cartesian:
        x = np.sqrt(3) * q + np.sqrt(3)/2 * r
        y = -1.5 * r
        if demand is not None:
            color = cmap(norm(demand[n]))
        else:
            color = None
            
        h = RegularPolygon((x, y), numVertices=6, radius=1, orientation=0, facecolor=color, alpha=alpha, edgecolor='k')
        if n in interference_nodes or n in demand_violations:
            h.set_hatch('\\')

        ax.add_patch(h)
        ax.text(x, y, '{}'.format(n), ha='center', va='center')

    ax.axis('off')
    ax.autoscale_view()

    if demand is not None:
        cbar = f.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), alpha=alpha, orientation='horizontal')
        cbar.set_label('Demand')


if __name__ == '__main__':

    demand, nfreq, reuse_distances = load_problem("P1")

    plot_nodes(range(1, len(qvals)+1), demand)

    plt.show()
