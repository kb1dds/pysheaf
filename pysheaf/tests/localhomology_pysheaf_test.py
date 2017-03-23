# Unit tests for local homology computations using the pysheaf library
#
# Copyright (c) 2015, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 

import numpy as np
import pysheaf as ps

# Test 1: line segment localized to the edge
toplexes=[[1,2]]
cplx=ps.AbstractSimplicialComplex(toplexes)

H0=cplx.localHomology(0,[0])
H1=cplx.localHomology(1,[0])
H2=cplx.localHomology(2,[0])

if H0.shape[1] != 0 or H1.shape[1] != 1 or H2.shape[1] != 0:
    print "Test 1 failed" 
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 1 passed"

# Test 2: line segment localized to the star over one vertex 
# (using localHomology())
toplexes=[[1,2]]
cplx=ps.AbstractSimplicialComplex(toplexes)

H0=cplx.localHomology(0,[1])
H1=cplx.localHomology(1,[1])
H2=cplx.localHomology(2,[1])

if H0.shape[1] != 0 or H1.shape[1] != 0 or H2.shape[1] != 0:
    print "Test 2 failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 2 passed"

# Test 3: line segment localized to the star over a vertex
# (using relative homology, "by hand")
toplexes=[[1,2]]
cplx=ps.AbstractSimplicialComplex(toplexes)

H0=cplx.homology(0,[2])
H1=cplx.homology(1,[2])
H2=cplx.homology(2,[2])

if H0.shape[1] != 0 or H1.shape[1] != 0 or H2.shape[1] != 0:
    print "Test 3 failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 3 passed"

# Test 4: empty triangle localized to an edge
toplexes=[[1,2],[1,3],[2,3]]
cplx=ps.AbstractSimplicialComplex(toplexes)

H0=cplx.localHomology(0,[0])
H1=cplx.localHomology(1,[0])
H2=cplx.localHomology(2,[0])

if H0.shape[1] != 0 or H1.shape[1] != 1 or H2.shape[1] != 0:
    print "Test 4 failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 4 passed"

# Test 5: solid triangle localized to star over an edge
# (using relative homology, "by hand")
toplexes=[[1,2,3]]
cplx=ps.AbstractSimplicialComplex(toplexes)

H0=cplx.homology(0,[2,3,4,5,6])
H1=cplx.homology(1,[2,3,4,5,6])
H2=cplx.homology(2,[2,3,4,5,6])

if H0.shape[1] != 0 or H1.shape[1] != 0 or H2.shape[1] != 0:
    print "Test 5 failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 5 passed"

# Test 6: solid triangle localized to the 2-face
toplexes=[[1,2,3]]
cplx=ps.AbstractSimplicialComplex(toplexes)

H0=cplx.localHomology(0,[0])
H1=cplx.localHomology(1,[0])
H2=cplx.localHomology(2,[0])

if H0.shape[1] != 0 or H1.shape[1] != 0 or H2.shape[1] != 1:
    print "Test 6 failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 6 passed"

# Test 7: solid triangle localized to the star over an edge
# (using localHomology)
toplexes=[[1,2,3]]
cplx=ps.AbstractSimplicialComplex(toplexes)

H0=cplx.localHomology(0,[1])
H1=cplx.localHomology(1,[1])
H2=cplx.localHomology(2,[1])

if H0.shape[1] != 0 or H1.shape[1] != 0 or H2.shape[1] != 0:
    print "Test 7 failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 7 passed"

# Test 8: two solid trianlges joined along an edge localized to the star over that edge
toplexes=[[1,2,3],[2,3,4]]
cplx=ps.AbstractSimplicialComplex(toplexes)

H0=cplx.localHomology(0,[4])
H1=cplx.localHomology(1,[4])
H2=cplx.localHomology(2,[4])

if H0.shape[1] != 0 or H1.shape[1] != 0 or H2.shape[1] != 1:
    print "Test 8 failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 8 passed"

# Test 9: two solid triangles joined along an edge localized to the star over some other edge
toplexes=[[1,2,3],[2,3,4]]
cplx=ps.AbstractSimplicialComplex(toplexes)

H0=cplx.localHomology(0,[2])
H1=cplx.localHomology(1,[2])
H2=cplx.localHomology(2,[2])

if H0.shape[1] != 0 or H1.shape[1] != 0 or H2.shape[1] != 0:
    print "Test 9 failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 9 passed"
