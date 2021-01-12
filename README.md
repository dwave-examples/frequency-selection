# Frequency Selection

Frequency selection problems have to do with assigning frequencies to
transmitting and receiving stations for wireless communication.  Various problem
definitions have been discussed in the literature (see, for example, [1]), but
the problems generally involve choosing the frequencies subject to constraints
associated with interference.  Interference occurs when two frequencies are both
close to each other in the electromagnetic band and assigned to connections that
are geographically close.

In this example, we consider the "feasibility frequency assignment problem" [1].
In this problem, there is a set of stations to which frequencies must be
assigned, and a set of available frequencies, and the goal is to identify
frequency assignments that satisfy two constraints:

- Constraint 1: The number of frequencies assigned to each station must equal
  the specified requirement (i.e., "demand") for that station.
- Constraint 2: No interference occurs.  Interference can be defined in terms of
  a forbidden set of frequency differences associated with each pair of
  stations.  For example, the problem may specify that frequencies assigned to
  the same station must differ by at least 5, whereas frequencies assigned to
  adjacent stations must differ by at least 2.

This example makes use of the Philadelphia benchmark instances introduced by
Ref. [2].  These problem instances involve a hexagonal grid with 21 stations, as
shown below.  Nine different problem instances are defined, labeled "P1" through
"P9", which each have associated demand and frequency reuse distance
(interference) specifications.  The schematic below is colored according to the
demand values for problem instance P1.  Further details are given in Refs [2,3].

![Philadelphia instance](_static/Philadelphia.png)

The Philadelphia instances have been widely used to study the minimum span
frequency assignment problem.  Summaries of reported results can be found at
http://fap.zib.de/problems/Philadelphia/ [3].  The minimum span frequency
assignment problem involves finding the frequency assignments that satisfy the
constraints defined above, while minimizing the difference between the largest
and smallest frequencies used.  One can view feasibility frequency assignment as
a subproblem of minimum span frequency assignment, in which the feasible
solution is searched for using increasingly smaller sets of candidate
frequencies.


## Usage

To run the demonstration, execute:

```bash
python frequency.py --show-plot
```

Or:

```bash
python frequency.py --save-plot
```

This will load the default "small" problem, which uses only a subset of the
stations from the full Philadelphia instance and reduces the demands.  The
`--show-plot` and `--save-plot` flags are optional, and they create a plot of
the stations, colored by their relative demands; any constraint violations
(either demand or interference) are indicated by a hash fill pattern in the
cells corresponding to the station in which a constraint was violated.  The
`--show-plot` flag displays an interactive plot using matplotlib, and
`--save-plot` saves the plot to a file.

To run the analysis for the full "P1" Philadelphia problem instance, execute:

```bash
python frequency.py P1
```

For each available problem definition, a default value is defined for the number
of candidate frequencies.  This can be overridden using the `-n` option.  Use
`python frequency.py -h` for a description of all command options and available
problems.


## Code Overview

The code consists of the following steps:

1. Based on the user input, the problem specification is loaded, which defines
   the set of demand values, the reuse distances, and the total number of
   available frequencies.
2. A binary quadratic model is constructed to represent the problem definition.
   The objective function is formed by encoding all of the constraints as
   penalty functions.
3. The problem is solved using a combination of classical and quantum computing
   resources via `LeapHybridSampler`.
4. Information about the feasibility of the solution is printed to the screen,
   optionally displaying a plot of the stations to indicate any constraint
   violations.


## Code Specifics

The formulation of the feasibility frequency assignment problem as a constraint
satisfaction problem with binary variables is discussed in Ref. [1].  The binary
variables are denoted by `x_{vf}`, which are indicators for whether frequency
`f` is selected for station `v`.  Mathematically, the constraints associated
with meeting the specified demand at each station are specified as:

![Eq1](_static/Eq1.png)

where `m(v)` denotes the demand at station `v`.  The constraints associated with
interference are given by:

![Eq2](_static/Eq2.png)

Here, a constraint is added whenever the combination of assignments `x_{vf}` and
`x_{gw}` would result in interference.

The constraint satisfaction problem is formulated as a binary quadratic model by
constructing an objective function that includes a penalty whenever a constraint
is violated.  The penalties associated with the demand constraints can be formed
by squaring the difference between the selected and demanded number of
frequencies at each station:

![H1](_static/H1.png)

For the interference constraints, we simply penalize interference by adding a
positive interaction coefficient for each case in which the pair `(x_{vf},
x_{gw})` would result in interference.


## References

[1] Aardal, K.I., van Hoesel, S.P.M., Koster, A.M.C.A. et al. Models and
solution techniques for frequency assignment problems. 4OR 1, 261–317,
2003. https://doi.org/10.1007/s10288-003-0022-6.

[2] Anderson, LG. A simulation study of some dynamic channel assignment
algorithms in a high capacity mobile telecommunications system. IEEE
Transactions on Communications 21, 1294– 1301, 1973.

[3] URL: http://fap.zib.de. Maintained by A. Eisenblätter and
A. M. C. A. Koster.


## License

Released under the Apache License 2.0. See [LICENSE](LICENSE) file.