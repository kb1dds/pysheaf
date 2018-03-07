# Unit test for the set cover measurement tools

import unittest
import pysheaf.covers as covers

class CoverTest(unittest.TestCase):
    """Test case based on the collection of all irredundant set covers of a set with four elements.  The results are given in Figure 2 on page 5."""
    def setUp(self):
        self.questions=[([[1,2,3,4]], 1,0.,0.),
                        ([[1,2,3],[1,2,4],[1,3,4],[2,3,4]], 90./99,12,8),
                        ([[1,2,4],[1,3,4],[2,3,4]], 81./99,6,5),
                        ([[1,2],[1,3,4],[2,3,4]], 72./99,4,4),
                        ([[1,2,3],[2,3,4]], 63./99,2,2),([[1,2,3],[1,4],[2,4],[3,4]],63./99,6,5),
                        ([[1,2],[2,3],[1,3,4]], 54./99,3,3),([[1,2],[1,3],[2,3],[1,4],[2,4],[3,4]],54./99,12,8),
                        ([[1,2,3],[3,4]], 45./99,1,1),([[1,2],[2,3],[2,4],[1,4],[3,4]], 45./99,8,6),
                        ([[1,2,3],[4]], 36./99,0,0),([[1,2],[2,3],[3,4],[1,4]], 36./99,4,4),([[1,2],[1,3],[2,3],[3,4]],36./99,5,4),
                        ([[1,2],[2,3],[3,4]], 27./99,2,2),([[1,2],[2,3],[2,4]], 27./99,3,2),([[1,2],[2,3],[1,3],[4]],27./99,3,3),
                        ([[1,2],[3,4]], 18./99,0,0),([[1,2],[2,3],[4]], 18./99,1,1),
                        ([[1],[2,3],[4]], 9./99,0,0),
                        ([[1],[2],[3],[4]], 0,0,0)]

    def stats(self,cover):
        return (covers.normalized_coarseness(cover),covers.pairwise_overlap(cover),covers.elementwise_overlap(cover))

    def test_FourElement(self):
        for cover,coarse,pover,eover in self.questions:
            sc,sp,se = self.stats(cover)
            #print str(cover) + ' rho=' + str(sc) +', delta=' + str(sp) + ', deltaprime=' + str(se)
            self.assertTrue(abs(sc-coarse)<1e-8)
            self.assertTrue(abs(sp-pover)<1e-8)
            self.assertTrue(abs(se-eover)<1e-8)
        return True

if __name__ == "__main__":
    # Run the unit tests upon import
    suite = unittest.TestLoader().loadTestsFromTestCase(CoverTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
