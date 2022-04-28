#!/usr/bin/env python

import unittest
import os
import sys
import subprocess

# TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..',
#                                       'rnapolii'))

scripts_dir = '../scripts/'
modeling_scripts = ['mhr_modeling.py','mhm_modeling.py','nude_modeling.py','mhr_xl_ctrl_modeling.py']

class Tests(unittest.TestCase):
    def test_complete(self):
        """Test modeling and analysis"""
        # Run modeling
        os.chdir(f'{scripts_dir}modeling')
        for script in modeling_scripts:
            p = subprocess.check_call(["/home/shreyas/imp-clean/build/setup_environment.sh","python", script, "test", "0"])
            self.assertTrue(os.path.exists('run_0/rmfs/0.rmf3'))

        # Run clustering
        os.chdir(os.path.join(TOPDIR, 'analysis'))
        p = subprocess.check_call(["$imp python", 'clustering.py', "--test"])
        self.assertTrue(os.path.exists('kmeans_5_1/dist_matrix.pdf'))
        self.assertTrue(os.path.exists('kmeans_5_1/cluster.0/0.rmf3'))

        # Test analysis
        p = subprocess.check_call(["python", 'precision_rmsf.py'])
        self.assertTrue(os.path.exists('kmeans_5_1/precision.0.0.out'))

        p = subprocess.check_call(["python", 'accuracy.py'])

if __name__ == '__main__':
    unittest.main()
