'''
Created on Mar 31, 2017

@author: prag717
'''
import unittest
import pysheaf as ps
import numpy as np


class TestSheaf(unittest.TestCase):
    
    def setUp(self):
        self.cellA = ps.SheafCell(0,cofaces=[
         ps.SheafCoface(index=2,orientation=-1,restriction=np.matrix(1)),
         ps.SheafCoface(index=3,orientation=1,restriction=np.matrix(1))],name='A')
        self.cellB = ps.SheafCell(0,cofaces=[
         ps.SheafCoface(index=2,orientation=1,restriction=np.matrix(1)),
         ps.SheafCoface(index=3,orientation=-1,restriction=np.matrix(1))],name='B')
        self.cellAB = ps.SheafCell(1,name='AB')
        self.cellBA = ps.SheafCell(1,name='BA')
        self.CircSheaf=ps.Sheaf([self.cellA,self.cellB,self.cellAB,self.cellBA]) # Edge

    def test_sheaf_cell_instance(self):
        self.assertTrue(isinstance(self.cellA,ps.SheafCell))
        self.assertTrue(isinstance(self.CircSheaf,ps.Sheaf))

    def tearDown(self):
        pass


    def testName(self):
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()