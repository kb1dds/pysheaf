#!/bin/python
#
# Cover measures test code
# Based on
#
# Cliff Joslyn, Kathleen Nowak, Emilie Purvine, "Granularity measures for irredundant set covers", unpublished, 2017.

from numpy import ceil
from scipy.special import comb
from itertools import chain, combinations

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

def partitions_iter(elements,currentset=[]):
    """Iterate over all partitions of a given set"""
    if not elements:
        yield currentset
    else:
        for newset in powerset(elements):
            if newset:
                for p in partitions_iter([a for a in elements if a not in list(squash(currentset))+list(newset)],currentset+[list(newset)]):
                    yield p
