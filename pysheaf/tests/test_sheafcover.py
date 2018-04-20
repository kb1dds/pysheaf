import unittest
import pysheaf as ps
import numpy as np


class TestSheafCover(unittest.TestCase):
    def setUp(self):
        self.testSheaf=ps.Sheaf([ps.SheafCell(dimension=0,stalkDim=1,cofaces=[ps.SheafCoface(index=1,restriction=np.array([1])),
                                                                              ps.SheafCoface(index=2,restriction=np.array([1]))]),
                                 ps.SheafCell(dimension=1,stalkDim=1),
                                 ps.SheafCell(dimension=1,stalkDim=1)])
        self.asg=ps.Section([ps.SectionCell(support=0,value=0),
                             ps.SectionCell(support=1,value=0),
                             ps.SectionCell(support=2,value=1)])

        self.testSheaf2=ps.Sheaf([ps.SheafCell(dimension=0,stalkDim=1,cofaces=[ps.SheafCoface(index=1,restriction=np.array([1])),
                                                                              ps.SheafCoface(index=2,restriction=np.array([1]))]),
                                  ps.SheafCell(dimension=1,stalkDim=1),
                                  ps.SheafCell(dimension=1,stalkDim=1),
                                  ps.SheafCell(dimension=0,stalkDim=1,cofaces=[ps.SheafCoface(index=1,restriction=np.array([1]))])])
        self.asg2=ps.Section([ps.SectionCell(support=0,value=0),
                              ps.SectionCell(support=1,value=0),
                              ps.SectionCell(support=2,value=1),
                              ps.SectionCell(support=3,value=0)])
        self.asg3=ps.Section([ps.SectionCell(support=0,value=0),
                              ps.SectionCell(support=1,value=0),
                              ps.SectionCell(support=2,value=1),
                              ps.SectionCell(support=3,value=1)])

    def test_consistencyRadius1(self):
        crs = self.testSheaf.consistencyRadii(self.asg)
        self.assertEqual(len(crs),2)
        self.assertTrue(np.abs(crs[0])<1e-5)
        self.assertTrue(np.abs(crs[1]-1)<1e-5)
        self.assertTrue(np.abs(self.testSheaf.consistencyRadius(self.asg)-1.0) < 1e-5)
        
    def test_consistencyRadius2(self):
        self.assertTrue(np.abs(self.testSheaf2.consistencyRadius(self.asg2)-1.0) < 1e-5)
        self.assertTrue(np.abs(self.testSheaf2.consistencyRadius(self.asg3)-1.0) < 1e-5)
    
    def test_covering1(self):
        cover=self.testSheaf.consistentCover(self.asg,threshold=0.5)
        self.assertEqual(len(cover),2)
        self.assertEqual(cover[0],{0,1,2})
        self.assertEqual(cover[1],{2})
        
    def test_covering2(self):
        cover=self.testSheaf.consistentCover(self.asg,threshold=1.5)
        self.assertEqual(len(cover),1)
        self.assertEqual(cover[0],{0,1,2})

    def test_covering3(self):
        cover=self.testSheaf2.consistentCover(self.asg2,threshold=0.5)
        self.assertEqual(len(cover),2)
        self.assertEqual(cover[0],{0,1,2,3})
        self.assertEqual(cover[1],{2})

    def test_covering4(self):
        cover=self.testSheaf2.consistentCover(self.asg3,threshold=0.5)
        self.assertEqual(len(cover),2)
        self.assertEqual(cover[0],{0,1,2})
        self.assertEqual(cover[1],{1,2,3})

if __name__ == '__main__':
    unittest.main(verbosity=2)
