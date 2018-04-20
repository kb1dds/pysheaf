import unittest
import pysheaf as ps
import numpy as np


class TestCell(unittest.TestCase):
    def setUp(self):
        self.testSheaf=ps.Sheaf([ps.SheafCell(dimension=0,stalkDim=1,cofaces=[ps.SheafCoface(index=1,restriction=np.array([1])), ps.SheafCoface(index=2,restriction=np.array([1]))]), ps.SheafCell(dimension=1,stalkDim=1),ps.SheafCell(dimension=1,stalkDim=1)])
        self.asg=ps.Section([ps.SectionCell(support=0,value=0),ps.SectionCell(support=1,value=0),ps.SectionCell(support=2,value=1)])

    def test_consistencyRadius(self):
        crs = self.testSheaf.consistencyRadii(self.asg)
        self.assertEqual(len(crs),2)
        self.assertTrue(np.abs(crs[0])<1e-5)
        self.assertTrue(np.abs(crs[1]-1)<1e-5)
        self.assertTrue(np.abs(self.testSheaf.consistencyRadius(self.asg)-1.0) < 1e-5)
    
    def test_covering(self):
        cover=self.testSheaf.consistentCover(self.asg,threshold=0.5)
        self.assertEqual(len(cover),2)
        self.assertEqual(cover[0],{0,1,2})
        self.assertEqual(cover[1],{2})

        cover=self.testSheaf.consistentCover(self.asg,threshold=1.5)
        self.assertEqual(len(cover),1)
        self.assertEqual(cover[0],{0,1,2})

if __name__ == '__main__':
    unittest.main(verbosity=2)
