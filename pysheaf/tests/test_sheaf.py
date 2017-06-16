'''
Created on Mar 31, 2017

@author: prag717
'''
import unittest
import pysheaf as ps
import numpy as np
from numpy.linalg import matrix_rank


class TestSheaf(unittest.TestCase):

    def setUp(self):
        ## From p. 98 Example 4.9:
        ## Cells=[A,B,C,D,E,F,G] inherit id from index of this list
        ## Restriction Maps:
        self.restAD = np.matrix([1,0,0,0])
        self.restAG = np.matrix([0,0,1,-1])
        self.restBD = np.matrix([1,0,0])
        self.restBE = np.matrix([[0,0,3],[0,2,0]])
        self.restBF = np.matrix([0,1,0])
        self.restCE = np.matrix([[0,0,3,0],[0,0,0,2]])
        self.restCF = np.matrix([0,0,0,1])
        self.restCG = np.matrix([1,-1,0,0])
        ## cofaces
        self.cofaceAD = ps.SheafCoface(3,1,self.restAD)
        self.cofaceAG = ps.SheafCoface(6,1,self.restAG)
        self.cofaceBD = ps.SheafCoface(3,-1,self.restBD)
        self.cofaceBE = ps.SheafCoface(4,1,self.restBE)
        self.cofaceBF = ps.SheafCoface(5,1,self.restBF)
        self.cofaceCE = ps.SheafCoface(4,-1,self.restCE)
        self.cofaceCF = ps.SheafCoface(5,-1,self.restCF)
        self.cofaceCG = ps.SheafCoface(6,-1,self.restCG)

        self.cellA = ps.SheafCell(dimension=0,id=0,stalkDim=4,name='A',
                             cofaces=[self.cofaceAD,self.cofaceAG])
        self.cellB = ps.SheafCell(dimension=0,id=1,stalkDim=3,name='B',
                             cofaces=[self.cofaceBD,self.cofaceBE,self.cofaceBF])
        self.cellC = ps.SheafCell(dimension=0,id=2,stalkDim=4,name='C',
                             cofaces=[self.cofaceCE,self.cofaceCF,self.cofaceCG])
        self.cellD = ps.SheafCell(dimension=1,id=3,stalkDim=1,name='D')
        self.cellE = ps.SheafCell(dimension=1,id=4,stalkDim=2,name='E')
        self.cellF = ps.SheafCell(dimension=1,id=5,stalkDim=1,name='F')
        self.cellG = ps.SheafCell(dimension=1,id=6,stalkDim=1,name='G')

        self.sheaf = ps.Sheaf([self.cellA,self.cellB,self.cellC,self.cellD,self.cellE,self.cellF,self.cellG])

    def test_sheaf_cell_instance(self):
        self.assertTrue(isinstance(self.cofaceAD,ps.SheafCoface))
        self.assertTrue(isinstance(self.cellA,ps.SheafCell))
        self.assertTrue(isinstance(self.sheaf,ps.Sheaf))

    def test_coboundary(self):
        d = self.sheaf.coboundary(0)
        m,n = d.shape
        gt = np.array([[ 1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0],
                       [ 0,  0,  0,  0,  0,  0,  3,  0,  0, -3,  0],
                       [ 0,  0,  0,  0,  0,  2,  0,  0,  0,  0, -2],
                       [ 0,  0,  0,  0,  0,  1,  0,  0,  0,  0, -1],
                       [ 0,  0,  1, -1,  0,  0,  0, -1,  1,  0,  0]])
        self.assertEqual(np.sum(d==gt),m*n)

    def test_cohomology(self):
        d = self.sheaf.cohomology(0).transpose()
        hk = np.matrix('0 0 1 0 0 0 0;1 0 0 0 0 0 0;0 1 0 1 -1 0 0;0 1 0 0 0 0 0;0 0 1 0 0 0 0;0 0 0 0 0 0 1; 0 0 0 0 0 1 0;0 0 0 1 0 0 0 ;0 0 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 0 1')
        gt = np.array(hk).transpose()
        tmat = np.concatenate((d,gt))

        self.assertEqual(d.shape,gt.shape)
        self.assertEqual(matrix_rank(d),matrix_rank(gt))
        self.assertEqual(matrix_rank(gt),matrix_rank(tmat))

    def test_cobetti(self):
        d = self.sheaf.coboundary(0)
        m,n = d.shape
        self.assertEqual(np.sum(d==d.astype(float)),m*n)
        d = d.astype(float)   ##### this avoids a type error when calculating matrix_rank
        r = matrix_rank(d)
        ker = d.shape[1] - r
        self.assertEqual(self.sheaf.cobetti(0),ker)


    def tearDown(self):
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
