#!/usr/bin/env python
import unittest
import os
import shutil
import sys
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '../'))

class Tests(unittest.TestCase):
    def test_mhr(self):
        """Test MHR simulations"""
        os.chdir(os.path.join(TOPDIR, 'scripts/modeling'))
        if os.path.exists('output'):
            shutil.rmtree('output')
        if os.path.exists('run_test'):
            shutil.rmtree('run_test')
        p = subprocess.check_call(["/home/shreyas/imp-clean/build/setup_environment.sh", "python", "mhr_modeling.py", "test", "mhr-test"])

        # require that the number of frames is present
        total_num_lines_stat_files = 0
        for i in range(1):
            with open("run_mhr-test/stat."+str(i)+".out", 'r') as statf:
                total_num_lines_stat_files += len(statf.readlines())
        self.assertEqual(total_num_lines_stat_files,6)

        # require that output files were produced
        for i in range(1):
            os.unlink("run_mhr-test/rmfs/"+str(i)+".rmf3")
            os.unlink("run_mhr-test/stat."+str(i)+".out")
            os.unlink("run_mhr-test/stat_replica."+str(i)+".out")


    def test_mhm(self):
        """Test MHM simulations"""
        os.chdir(os.path.join(TOPDIR, 'scripts/modeling'))
        if os.path.exists('output'):
            shutil.rmtree('output')

        p = subprocess.check_call(["/home/shreyas/imp-clean/build/setup_environment.sh", "python", "mhm_modeling.py", "test", "mhm-test"])

        # require that the number of frames is present
        total_num_lines_stat_files = 0
        for i in range(1):
            with open("run_mhm-test/stat."+str(i)+".out", 'r') as statf:
                total_num_lines_stat_files += len(statf.readlines())
        self.assertEqual(total_num_lines_stat_files,6)

        # require that output files were produced
        for i in range(1):
            os.unlink("run_mhm-test/rmfs/"+str(i)+".rmf3")
            os.unlink("run_mhm-test/stat."+str(i)+".out")
            os.unlink("run_mhm-test/stat_replica."+str(i)+".out")


    def test_nude(self):
        """Test NuDe simulations"""
        os.chdir(os.path.join(TOPDIR, 'scripts/modeling'))
        if os.path.exists('output'):
            shutil.rmtree('output')

        p = subprocess.check_call(["/home/shreyas/imp-clean/build/setup_environment.sh", "python", "nude_modeling.py", "test", "nude-test"])

        # require that the number of frames is present
        total_num_lines_stat_files = 0
        for i in range(1):
            with open("run_nude-test/stat."+str(i)+".out", 'r') as statf:
                total_num_lines_stat_files += len(statf.readlines())
        self.assertEqual(total_num_lines_stat_files,6)

        # require that output files were produced
        for i in range(1):
            os.unlink("run_nude-test/rmfs/"+str(i)+".rmf3")
            os.unlink("run_nude-test/stat."+str(i)+".out")
            os.unlink("run_nude-test/stat_replica."+str(i)+".out")


if __name__ == '__main__':
    unittest.main()
