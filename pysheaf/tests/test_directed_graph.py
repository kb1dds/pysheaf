'''
Created on Mar 31, 2017

@author: prag717
'''
import unittest
import pysheaf as ps


class TestDirectedGraph(unittest.TestCase):


    def setUp(self):
        self.dg = ps.DirectedGraph([(0,1),(1,2),(1,2),(0,2)])
        
    def test_directed_graph_instance(self):
        self.assertTrue(isinstance(self.dg,ps.DirectedGraph))
        
    def test_max_flow(self):
        mf = self.dg.maxFlow(4,6)        
        self.assertTrue( mf == [[4, 0, 5, 1, 6], [4, 3, 6]]  )


    def tearDown(self):
        pass


    def testName(self):
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()