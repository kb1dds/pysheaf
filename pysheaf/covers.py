#!/bin/python
#
# Cover measures test code
# Based on
#
# Cliff Joslyn, Kathleen Nowak, Emilie Purvine, "Granularity measures for irredundant set covers", unpublished, 2017.

from numpy import ceil
from scipy.special import comb
from itertools import chain, combinations

import unittest

def squash(cover):
    """Determine the unique elements in a cover"""
    return set([e for a in cover for e in a])

def powerset(iterable):
    """Return an iterator walking over power set of a given set"""
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def normalized_coarseness(cover):
    """Compute rho, the normalized coarseness measure.  Max (single blob) = 1, min (dust) = 0.  Uses formula on page 7 in Section 3.3"""
    n=len(squash(cover))
    return (len(set([e for a in cover for e in powerset(a)]))-(n+1))/(2.0**n-(n+1))

def pairwise_overlap(cover):
    """Compute delta, the pairwise overlap measure.  Uses formula on page 8."""
    return sum([comb(lm,2) for e,lm in memberdict(cover).iteritems()])

def elementwise_overlap(cover):
    """Compute deltaprime, the elementwise overlap measure.  Uses formula on page 7."""
    return sum([lm-1 for e,lm in memberdict(cover).iteritems()])

def normalized_pairwise_overlap(cover):
    """Compute deltabar, the normalized pairwise overlap measure.  Max = 1, min (partition) = 0.  Uses formula page 10."""
    n=len(squash(cover))
    return pairwise_overlap(cover)/(1.0*n*comb(comb(n-1,ceil(n/2.0)-1),2))

def normalized_elementwise_overlap(cover):
    """Compute deltaprimebar, the normalized elementwise overlap measure.  Max = 1, min (partition) = 0.  Uses formula page 10."""
    n=len(squash(cover))
    return elementwise_overlap(cover)/(1.0*n*(comb(n-1,ceil(n/2.0)-1)-1))

def memberdict(cover):
    """Determine lambda, the elementwise membership vector, defined on page 7.  This is implemented as a dictionary indexed by element."""
    return {e:len([a for a in cover if e in a]) for e in squash(cover)}

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
        return (normalized_coarseness(cover),pairwise_overlap(cover),elementwise_overlap(cover))

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
