import unittest
import pysheaf as ps
import numpy as np


class TestConsistencyRadius(unittest.TestCase):
    def setUp(self):
        self.testSheaf1=ps.Sheaf([])

        # Small sheaf
        self.testSheaf2=ps.Sheaf([ps.SheafCell(dimension=1,stalkDim=1,cofaces=[ps.SheafCoface(index=2,orientation=np.array(1),restriction=np.array([0.5]))]),
                                  ps.SheafCell(dimension=1,stalkDim=1,cofaces=[ps.SheafCoface(index=2,orientation=np.array(-1),restriction=np.array([1]))]),
                                  ps.SheafCell(dimension=2,stalkDim=1,cofaces=[]),
                                  ps.SheafCell(dimension=0,stalkDim=1,cofaces=[ps.SheafCoface(index=0,orientation=np.array(1),restriction=np.array([2])),
                                                                               ps.SheafCoface(index=1,orientation=np.array(1),restriction=np.array([1]))])])

        self.asg2_1 = ps.Assignment([ps.AssignmentCell(support=0,value=np.array(0.)),
                                  ps.AssignmentCell(support=1,value=np.array(1.))])

        self.asg2_3 = ps.Assignment([ps.AssignmentCell(support=0,value=np.array(0.)),
                                  ps.AssignmentCell(support=1,value=np.array(1.)),
                                  ps.AssignmentCell(support=2,value=np.array(0.5)),
                                  ps.AssignmentCell(support=3,value=np.array(1./3))])

    def test_consistencyRadius1(self):
        '''Test consistency of a global Assignment'''
        pass

    def test_consistencyRadius2(self):
        '''Test consistency of a non-global Assignment'''
        pass

    def test_fuseAssignment1(self):
        '''Can we replace a single missing value?'''
        pass

    def test_localConsistencyRadius(self):
        '''Test consistency radius when restricted to a subset'''
        pass

    def test_consistencyRadius3(self):
        '''Test consistency radius of a small Assignment on a small sheaf (no pre-fusion)'''
        self.assertAlmostEqual(0,self.testSheaf2.consistencyRadius(self.asg2_1,testSupport=[0,2]),4) # Local star 0
        self.assertAlmostEqual(0,self.testSheaf2.consistencyRadius(self.asg2_1,testSupport=[1,2]),4) # Local star 1
        self.assertAlmostEqual(0,self.testSheaf2.consistencyRadius(self.asg2_1,testSupport=[2]),4) # Local star 2
        self.assertAlmostEqual(0.5,self.testSheaf2.consistencyRadius(self.asg2_1),4) # Union
        
    def test_fuseAssignment2(self):
        '''Test consistency radius after fusing on a small sheaf'''
        # Do the fusion to cell 2
        asg2_2 = self.testSheaf2.fuseAssignment(self.asg2_1,activeCells=[0,1,2],testSupport=[0,1,2])
        self.assertAlmostEqual(2./3,asg2_2.assignmentCells[0].value,2) # Check that the correct value got found
        self.assertAlmostEqual(1./3,asg2_2.assignmentCells[1].value,2) # Check that the correct value got found
        self.assertAlmostEqual(1./3,asg2_2.assignmentCells[2].value,2) # Check that the correct value got found
        
        # Compute consistency radii
        self.assertAlmostEqual(0.,self.testSheaf2.consistencyRadius(asg2_2,testSupport=[0,2]),2) # Local star 0
        self.assertAlmostEqual(0.,self.testSheaf2.consistencyRadius(asg2_2,testSupport=[1,2]),2) # Local star 1
        self.assertAlmostEqual(0,self.testSheaf2.consistencyRadius(asg2_2,testSupport=[2]),4) # Local star 2
        self.assertAlmostEqual(0.,self.testSheaf2.consistencyRadius(asg2_2),4) # Union

    def test_minimizeConsistencyRadius1(self):
        '''Test consistency radius after extending assignment on a small sheaf'''
        # Do the fusion to cell 2
        asg2_2 = self.testSheaf2.minimizeConsistencyRadius(self.asg2_1,activeCells=[2],testSupport=[0,1,2])
        self.assertAlmostEqual(0.5,asg2_2.assignmentCells[2].value,2) # Check that the correct value got found

        # Compute consistency radii
        self.assertAlmostEqual(0.5,self.testSheaf2.consistencyRadius(asg2_2,testSupport=[0,2]),2) # Local star 0
        self.assertAlmostEqual(0.5,self.testSheaf2.consistencyRadius(asg2_2,testSupport=[1,2]),2) # Local star 1
        self.assertAlmostEqual(0,self.testSheaf2.consistencyRadius(asg2_2,testSupport=[2]),4) # Local star 2
        self.assertAlmostEqual(0.5,self.testSheaf2.consistencyRadius(asg2_2),4) # Union

    def test_minimizeConsistencyRadius2(self):
        '''Test consistency radius after extending assignment to unions on a small sheaf'''
        # Do the fusion to cells 2 and 3
        asg2_3 = self.testSheaf2.minimizeConsistencyRadius(self.asg2_1,activeCells=[2,3])

        self.assertLess(1./3,asg2_3.assignmentCells[2].value)
        self.assertLess(asg2_3.assignmentCells[2].value,2./3)
        self.assertAlmostEqual(1./3,asg2_3.assignmentCells[3].value,2) # Check that the correct value got found

        # Compute consistency radii
        self.assertAlmostEqual(0,self.testSheaf2.consistencyRadius(asg2_3,testSupport=[2]),4) # Local star 2
        self.assertAlmostEqual(2./3,self.testSheaf2.consistencyRadius(asg2_3),4) # Union
    

if __name__ == '__main__':
    unittest.main(verbosity=2)
