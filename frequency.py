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

import numpy as np

# Ignore errors importing matplotlib.pyplot (may not be available in
# testing framework)
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass

import dimod
from dwave.system import LeapHybridSampler

from philadelphia import load_problem, get_forbidden_set, plot_nodes


def construct_bqm(demand, nfreq, reuse_distances, LAGRANGE=1.0):
    """Construct BQM for feasibility frequency assignment problem
    
    Args:
        demand (dict):
            Dictionary mapping each node number to a demand value
        nfreq (int):
            Number of frequencies to consider
        reuse_distances (list):
            List of reuse distances
        LAGRANGE (float):
            Lagrange multiplier.  Not needed in current formulation,
            which does not include an objective component of the
            problem formulation.  Retained only as a placeholder in
            case the problem is extended to include an objective.

    Returns:
        AdjVectorBQM
    """
    # Variables:
    # x_vf, v in nodes, f in frequencies: Is f assigned to node v?

    nodes = sorted(list(demand.keys()))
    n_nodes = len(nodes)

    bqm = dimod.AdjVectorBQM(dimod.BINARY)

    # Constraints to enforce demand at each node:
    # Sum_f[ (1-2C) xvf ] + Sum_j>i[ 2 xvfi xvfj ] + C^2

    # Linear parts:
    for v in nodes:
        for f in range(nfreq):
            bqm.add_variable('x_{}_{}'.format(v, f), LAGRANGE * (1.0 - 2*demand[v]))
    # Interactions:
    for v in nodes:
        for fi in range(nfreq):
            for fj in range(fi+1,nfreq):
                bqm.add_interaction('x_{}_{}'.format(v, fi), 'x_{}_{}'.format(v, fj), LAGRANGE * 2.0)


    # Define penalties associated with the interference constraints.
    # The interference constraints are represented by the inequality
    # xvf + xwg <= 1 for all combinations that would produce
    # interference.

    # First enforce the self-conflicts between frequencies in the same node:
    T = get_forbidden_set(1, 1, reuse_distances)
    for v in nodes:
        for f in range(nfreq):
            for g in range(f+1,nfreq):
                if abs(f-g) in T:
                    bqm.add_interaction('x_{}_{}'.format(v, f), 'x_{}_{}'.format(v, g), LAGRANGE)


    # Now enforce the cross-node conflicts:
    for iv,v in enumerate(nodes):
        for w in nodes[iv+1:]:
            T = get_forbidden_set(v, w, reuse_distances)
            if not T:
                # No disallowed frequencies at this distance
                continue
            for f in range(nfreq):
                # Note f and g are frequencies on different nodes, so we do need to look at all combinations
                for g in range(nfreq):
                    if abs(f-g) in T:
                        bqm.add_interaction('x_{}_{}'.format(v, f), 'x_{}_{}'.format(w, g), LAGRANGE)

    return bqm


def check_results(demand, nfreq, reuse_distances, sample, verbose=True):
    """Check whether a given solution vector satisfies the problem constraints

    Args:
        demand (dict):
            Dictionary mapping each node number to a demand value
        nfreq (int):
            Number of frequencies to consider
        reuse_distances (list):
            List of reuse distances
        sample (Sample):
            Solution vector to consider
        verbose (bool):
            If True, print additional details to screen

    Returns:
        dict:
            Dictionary that tracks counts and node lists for constraint violations
    """
    nodes = sorted(list(demand.keys()))
    frequencies = _get_frequencies(nodes, nfreq, sample)

    violation_dict = {
        'demand-count': 0,
        'self-count': 0,
        'cross-count': 0,
        'demand-nodes': [],
        'self-nodes': set(),
        'cross-nodes': set()
        }

    # Check for demand violations
    for v in nodes:
        if len(frequencies[v]) != demand[v]:
            violation_dict['demand-count'] += 1
            violation_dict['demand-nodes'].append(v)

    # Check for self-conflicts:
    T = get_forbidden_set(1, 1, reuse_distances)
    for v in nodes:
        violations = 0
        for i_f, fi in enumerate(frequencies[v]):
            for fj in frequencies[v][i_f+1:]:
                if abs(fi-fj) in T:
                    violations += 1
                    violation_dict['self-nodes'].add(v)
        violation_dict['self-count'] += violations
        if verbose:
            print('[{}]: Node {}: {} self conflicts'.format('X' if violations==0 else ' ', v, violations))

    # Check for cross-node conflicts:
    for iv, v in enumerate(nodes):
        violations = 0
        # Alternatively, to only count "new" conflicts, iterate over
        # range(iv+1,len(nodes))
        for iw, w in enumerate(nodes):
            if iw == iv:
                continue
            T = get_forbidden_set(v, w, reuse_distances)
            if not T:
                # No disallowed frequencies at this distance
                continue
            for fv in frequencies[v]:
                for fw in frequencies[w]:
                    if abs(fv-fw) in T:
                        violations += 1
                        violation_dict['cross-nodes'].add(v)
                        violation_dict['cross-nodes'].add(w)
        violation_dict['cross-count'] += violations
        if verbose:
            print('[{}]: Node {}: {} cross conflicts'.format('X' if violations==0 else ' ', v, violations))
    violation_dict['cross-count'] //= 2 # Each cross-violation was counted twice

    return violation_dict


def _get_frequencies(nodes, nfreq, sample):
    """Retrieve the frequencies selected for each node

    Returns:
        dict: Dictionary mapping the node to a list of frequencies
    """
    frequencies = {}
    for v in nodes:
        frequencies[v] = [i for i in range(nfreq) if sample.sample['x_{}_{}'.format(v,i)]]
    return frequencies


def _print_frequency_separations(reuse_distances, solution):
    """Print detailed information about frequency separation in solution"""
    if len(solution) > 10:
        # Avoid printing detail with full problem instance, which would be very lengthy
        print('Skipping interference detail for full problem instance')
        return

    print('Interference detail:')

    # Start with self-conflicts
    T = get_forbidden_set(1, 1, reuse_distances)
    print('  Station: min actual separation vs allowed separation')
    for node, frequencies in sorted(solution.items()):
        if len(frequencies) > 1:
            sep = min(frequencies[i+1] - frequencies[i] for i in range(len(frequencies)-1))
            allowed = max(T) + 1
            print('    {:2}: {:2} {:2} {:2}'.format(node, sep, '>=' if sep >= allowed else '<', allowed))

    print('  Station pair: min actual separation vs allowed separation')
    nodes = sorted(solution.keys())
    for i in range(len(nodes)):
        node1 = nodes[i]
        f1_vals = solution[node1]
        if not f1_vals:
            continue
        for j in range(i+1, len(nodes)):
            node2 = nodes[j]
            f2_vals = solution[node2]
            if not f2_vals:
                continue
            T = get_forbidden_set(node1, node2, reuse_distances)
            sep = min(abs(f1_vals[i] - f2_vals[j]) for i in range(len(f1_vals)) for j in range(len(f2_vals)))
            allowed = max(T) + 1
            print('    {:2}, {:2}: {:2} {:2} {:2}'.format(node1, node2, sep, '>=' if sep >= allowed else '<', allowed))


if __name__ == '__main__':
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(description="Run the frequency selection example on specified problem",
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=textwrap.dedent("""
                                     The Philadelphia problem instances have the following minimum
                                     span frequency ranges:

                                     - P1: 426
                                     - P2: 426
                                     - P3: 257
                                     - P4: 252
                                     - P5: 239
                                     - P6: 179
                                     - P7: 855
                                     - P8: 524
                                     - P9: 1713

                                     In theory, each problem instance has a feasible solution when
                                     NFREQ is greater than or equal to the minimum span frequency
                                     range plus 1
                                     """))
    parser.add_argument("problem", nargs="?", default="small",
                        choices=["trivial", "single", "small", "very-small"] + ["P{}".format(i) for i in range(1,10)],
                        help="problem to run (default: %(default)s)")
    parser.add_argument('-n', '--nfreq', default=None, help="number of frequencies to consider (default: problem-dependent)", type=int)
    parser.add_argument('--show-frequencies', action="store_true", help="print out selected frequencies")
    parser.add_argument('--verbose', action='store_true', help='print details about frequency separation in solution (not allowed for full problem instances)')
    parser.add_argument("--show-plot",  action='store_true', help="display plot of cell grid")
    parser.add_argument("--save-plot",  action='store_true', help="save plot of cell grid to file")

    args = parser.parse_args()

    demand, nfreq, reuse_distances = load_problem(args.problem)
    if args.nfreq:
        # Override problem-dependent default
        nfreq = args.nfreq
    print(nfreq, 'frequencies considered')

    bqm = construct_bqm(demand, nfreq, reuse_distances)

    print('{} variables'.format(bqm.num_variables))
    print('{} interactions'.format(bqm.num_interactions))
            
    sampler = LeapHybridSampler()
    results = sampler.sample(bqm)

    print('\nSolution:')
    violations = check_results(demand, nfreq, reuse_distances, results.first, verbose=False)

    print('{} demand violations'.format(violations['demand-count']))
    print('{} within-node frequency violations'.format(violations['self-count']))
    print('{} across-node frequency violations'.format(violations['cross-count']))
    print('')

    nodes = sorted(list(demand.keys()))
    frequencies = _get_frequencies(nodes, nfreq, results.first)
    if args.show_frequencies:
        for node, f in sorted(frequencies.items()):
            print('Station {}: {}'.format(node, f))
        print('')
    station_maximums = [max(freqs) for freqs in frequencies.values() if freqs]
    if station_maximums:
        print('Max frequency:', max(station_maximums))

    if args.verbose:
        print('')
        _print_frequency_separations(reuse_distances, frequencies)

    if args.show_plot or args.save_plot:
        interference = violations['self-nodes'].union(violations['cross-nodes'])
        plot_nodes(nodes, demand, interference, demand_violations=violations['demand-nodes'])

        if args.save_plot:
            filename = 'frequency_grid.png'
            plt.savefig(filename, bbox_inches='tight')
            print('Plot saved to:', filename)

        if args.show_plot:
            plt.show()
