from philadelphia import get_forbidden_set

def check_results(demand, nfreq, reuse_distances, sample, verbose=True):
    """Check whether a given solution vector satisfies the problem constraints.

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
    frequencies = get_frequencies(nodes, nfreq, sample)

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


def get_frequencies(nodes, nfreq, sample):
    """Retrieve the frequencies selected for each node.

    Returns:
        dict: Dictionary mapping the node to a list of frequencies
    """
    frequencies = {}
    for v in nodes:
        frequencies[v] = [i for i in range(nfreq) if sample.sample['x_{}_{}'.format(v,i)]]
    return frequencies


def print_frequency_separations(reuse_distances, solution):
    """Print detailed information about frequency separation in solution."""
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
