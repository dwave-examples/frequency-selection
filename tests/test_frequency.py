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

import unittest
import networkx as nx
import os
import sys
import subprocess

from neal import SimulatedAnnealingSampler

import frequency
from philadelphia import load_problem, get_forbidden_set, P1_REUSE_DISTANCES, P2_REUSE_DISTANCES

example_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class TestNetwork(unittest.TestCase):
    """Test reading in the network structure."""

    def test_small_network(self):
        demand, nfreq, reuse_distances  = load_problem('small')

        self.assertEqual(len(demand), 7)

    def test_P1_definition(self):
        demand, nfreq, reuse_distances = load_problem('P1')

        self.assertEqual(len(demand), 21)
        self.assertEqual(get_forbidden_set(1, 1, P1_REUSE_DISTANCES), {0, 1, 2, 3, 4})
        self.assertEqual(get_forbidden_set(1, 2, P1_REUSE_DISTANCES), {0, 1})
        self.assertEqual(get_forbidden_set(1, 3, P1_REUSE_DISTANCES), {0})
        self.assertEqual(get_forbidden_set(5, 18, P1_REUSE_DISTANCES), {0})
        self.assertEqual(get_forbidden_set(3, 16, P1_REUSE_DISTANCES), {0})
        self.assertEqual(get_forbidden_set(1, 5, P1_REUSE_DISTANCES), set())
        self.assertEqual(get_forbidden_set(20, 5, P1_REUSE_DISTANCES), set())
        self.assertEqual(get_forbidden_set(4, 15, P1_REUSE_DISTANCES), set())
        self.assertEqual(get_forbidden_set(3, 20, P1_REUSE_DISTANCES), {0})

        # Differs from P2:
        self.assertEqual(get_forbidden_set(4, 20, P1_REUSE_DISTANCES), {0})
        
    def test_P2_definition(self):

        self.assertEqual(get_forbidden_set(1, 1, P2_REUSE_DISTANCES), {0, 1, 2, 3, 4})
        self.assertEqual(get_forbidden_set(1, 2, P2_REUSE_DISTANCES), {0, 1})
        self.assertEqual(get_forbidden_set(1, 3, P2_REUSE_DISTANCES), {0})
        self.assertEqual(get_forbidden_set(5, 18, P2_REUSE_DISTANCES), {0})
        self.assertEqual(get_forbidden_set(3, 16, P2_REUSE_DISTANCES), {0})
        self.assertEqual(get_forbidden_set(1, 5, P1_REUSE_DISTANCES), set())
        self.assertEqual(get_forbidden_set(20, 5, P1_REUSE_DISTANCES), set())
        self.assertEqual(get_forbidden_set(4, 15, P2_REUSE_DISTANCES), set())
        # This case does not match provided diagrams if reuse
        # distances taken as given.  See notes in philadelphia.py.
        self.assertEqual(get_forbidden_set(3, 20, P2_REUSE_DISTANCES), {0})

        # Differs from P1
        self.assertEqual(get_forbidden_set(4, 20, P2_REUSE_DISTANCES), set())


class TestSmallProblem(unittest.TestCase):
    """Test solution to small problem."""

    def test_small(self):
        
        demand, nfreq, reuse_distances = load_problem('small')

        bqm = frequency.construct_bqm(demand, nfreq, reuse_distances)

        sampler = SimulatedAnnealingSampler()
        results = sampler.sample(bqm)

        violation_dict = frequency.check_results(demand, nfreq, reuse_distances, results.first, verbose=False)

        self.assertEqual(violation_dict['demand-count'], 0)
        self.assertEqual(violation_dict['self-count'], 0)
        self.assertEqual(violation_dict['cross-count'], 0)


class TestIntegration(unittest.TestCase):
    """Test execution as a script."""

    @unittest.skipIf(os.getenv('SKIP_INT_TESTS'), "Skipping integration test.")
    def test_integration(self):
        file_path = os.path.join(example_dir, "frequency.py")

        output = subprocess.check_output([sys.executable, file_path, 'trivial'])
        output = output.decode('utf-8') # Bytes to str

        self.assertIn('0 demand', output.lower())
        self.assertIn('0 within', output.lower())
        self.assertIn('0 across', output.lower())
