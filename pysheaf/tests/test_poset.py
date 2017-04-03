'''
Created on Mar 31, 2017

@author: prag717
'''
import unittest
import pysheaf as ps


class TestPoset(unittest.TestCase):


    def setUp(self):
        self.poset=ps.Poset([ps.Cell(0,True,[ps.Coface(1,1),ps.Coface(2,1)]),
                ps.Cell(1,True,[ps.Coface(3,1),ps.Coface(4,1)]),
                ps.Cell(1,True,[ps.Coface(3,1)]),
                ps.Cell(2,True,[ps.Coface(5,1)]),
                ps.Cell(2,True,[]),
                ps.Cell(3,True,[])])
        
    def test_poset_instance(self):
        self.assertTrue(isinstance(self.poset,ps.Poset)) 


    def tearDown(self):
        pass

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()