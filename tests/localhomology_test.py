# Unit tests for local homology computations
#
# Copyright (c) 2015-2017, Michael Robinson, Chris Capraro
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 

import simplicialHomology as sh
 
# Test 1: line segment localized to the edge
toplexes=[[1,2]]
   
H0=sh.localHomology(0,toplexes,[[1,2]],True)
H1=sh.localHomology(1,toplexes,[[1,2]],True)
H2=sh.localHomology(2,toplexes,[[1,2]],True)
   
if H0 != 0 or H1 != 1 or H2 != 0:
    print "Test 1a failed" 
    print "H0="+ str(H0)
    print "H1="+ str(H1)
    print "H2="+ str(H2)
else:
    print "Test 1a passed"
   
H0=sh.localHomology(0,toplexes,[[1,2]])
H1=sh.localHomology(1,toplexes,[[1,2]])
H2=sh.localHomology(2,toplexes,[[1,2]])
   
if H0.shape[1] != 0 or H1.shape[1] != 1 or H2.shape[1] != 0:
    print "Test 1b failed" 
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 1b passed"
    
# Test 2: line segment localized to the star over one vertex 
# (using localHomology())
toplexes=[[1,2]]
   
H0=sh.localHomology(0,toplexes,[[1]],True)
H1=sh.localHomology(1,toplexes,[[1]],True)
H2=sh.localHomology(2,toplexes,[[1]],True)
   
if H0 != 0 or H1 != 0 or H2 != 0:
    print "Test 2a failed" 
    print "H0="+ str(H0)
    print "H1="+ str(H1)
    print "H2="+ str(H2)
else:
    print "Test 2a passed"
    
H0=sh.localHomology(0,toplexes,[[1]])
H1=sh.localHomology(1,toplexes,[[1]])
H2=sh.localHomology(2,toplexes,[[1]])

if H0.shape[1] != 0 or H1.shape[1] != 0 or H2.shape[1] != 0:
    print "Test 2b failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 2b passed"
    
# Test 3: line segment localized to the star over a vertex
# (using relative homology, "by hand")
toplexes=[[1,2]]
   
H0=sh.simplicialHomology(0,toplexes,[[1]],True)
H1=sh.simplicialHomology(1,toplexes,[[1]],True)
H2=sh.simplicialHomology(2,toplexes,[[1]],True)
   
if H0 != 0 or H1 != 0 or H2 != 0:
    print "Test 3a failed" 
    print "H0="+ str(H0)
    print "H1="+ str(H1)
    print "H2="+ str(H2)
else:
    print "Test 3a passed"
 
H0=sh.simplicialHomology(0,toplexes,[[1]])
H1=sh.simplicialHomology(1,toplexes,[[1]])
H2=sh.simplicialHomology(2,toplexes,[[1]])
   
if H0.shape[1] != 0 or H1.shape[1] != 0 or H2.shape[1] != 0:
    print "Test 3 failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
if H0 != 0 or H1 != 0 or H2 != 0:
    print "Test 3b failed" 
    print "H0="+ str(H0)
    print "H1="+ str(H1)
    print "H2="+ str(H2)
else:
    print "Test 3b passed"
      
# Test 4: empty triangle localized to an edge
toplexes=[[1,2],[1,3],[2,3]]
   
H0=sh.localHomology(0,toplexes,[[1,2]],True)
H1=sh.localHomology(1,toplexes,[[1,2]],True)
H2=sh.localHomology(2,toplexes,[[1,2]],True)
   
if H0 != 0 or H1 != 1 or H2 != 0:
    print "Test 4a failed" 
    print "H0="+ str(H0)
    print "H1="+ str(H1)
    print "H2="+ str(H2)
else:
    print "Test 4a passed"
   
H0=sh.localHomology(0,toplexes,[[1,2]])
H1=sh.localHomology(1,toplexes,[[1,2]])
H2=sh.localHomology(2,toplexes,[[1,2]])
   
if H0.shape[1] != 0 or H1.shape[1] != 1 or H2.shape[1] != 0:
    print "Test 4b failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 4b passed"
    
# Test 5: solid triangle localized to star over an edge
# (using relative homology, "by hand")
toplexes=[[1,2,3]]
   
H0=sh.simplicialHomology(0,toplexes,[[1],[2],[3],[1,2],[2,3]],True)
H1=sh.simplicialHomology(1,toplexes,[[1],[2],[3],[1,2],[2,3]],True)
H2=sh.simplicialHomology(2,toplexes,[[1],[2],[3],[1,2],[2,3]],True)
   
if H0 != 0 or H1 != 0 or H2 != 0:
    print "Test 5a failed" 
    print "H0="+ str(H0)
    print "H1="+ str(H1)
    print "H2="+ str(H2)
else:
    print "Test 5a passed"

H0=sh.simplicialHomology(0,toplexes,[[1],[2],[3],[1,2],[2,3]])
H1=sh.simplicialHomology(1,toplexes,[[1],[2],[3],[1,2],[2,3]])
H2=sh.simplicialHomology(2,toplexes,[[1],[2],[3],[1,2],[2,3]])
   
if H0.shape[1] != 0 or H1.shape[1] != 0 or H2.shape[1] != 0:
    print "Test 5b failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 5b passed"
     
# Test 6: solid triangle localized to the 2-face
toplexes=[[1,2,3]]
   
H0=sh.localHomology(0,toplexes,[[1,2,3]],True)
H1=sh.localHomology(1,toplexes,[[1,2,3]],True)
H2=sh.localHomology(2,toplexes,[[1,2,3]],True)

if H0 != 0 or H1 != 0 or H2 != 1:
    print "Test 6a failed" 
    print "H0="+ str(H0)
    print "H1="+ str(H1)
    print "H2="+ str(H2)
else:
    print "Test 6a passed"

H0=sh.localHomology(0,toplexes,[[1,2,3]])
H1=sh.localHomology(1,toplexes,[[1,2,3]])
H2=sh.localHomology(2,toplexes,[[1,2,3]])

if H0.shape[1] != 0 or H1.shape[1] != 0 or H2.shape[1] != 1:
    print "Test 6b failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 6b passed"
    
# Test 7: solid triangle localized to the star over an edge
# (using localHomology)
toplexes=[[1,2,3]]
   
H0=sh.localHomology(0,toplexes,[[1,2]],True)
H1=sh.localHomology(1,toplexes,[[1,2]],True)
H2=sh.localHomology(2,toplexes,[[1,2]],True)
   
if H0 != 0 or H1 != 0 or H2 != 0:
    print "Test 7a failed" 
    print "H0="+ str(H0)
    print "H1="+ str(H1)
    print "H2="+ str(H2)
else:
    print "Test 7a passed"
   
H0=sh.localHomology(0,toplexes,[[1,2]])
H1=sh.localHomology(1,toplexes,[[1,2]])
H2=sh.localHomology(2,toplexes,[[1,2]])
   
if H0.shape[1] != 0 or H1.shape[1] != 0 or H2.shape[1] != 0:
    print "Test 7b failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 7b passed"
    
# Test 8: two solid triangles joined along an edge localized to the star over that edge
toplexes=[[1,2,3],[2,3,4]]
   
H0=sh.localHomology(0,toplexes,[[2,3]],True)
H1=sh.localHomology(1,toplexes,[[2,3]],True)
H2=sh.localHomology(2,toplexes,[[2,3]],True)

if H0 != 0 or H1 != 0 or H2 != 1:
    print "Test 8a failed" 
    print "H0="+ str(H0)
    print "H1="+ str(H1)
    print "H2="+ str(H2)
else:
    print "Test 8a passed"
   
H0=sh.localHomology(0,toplexes,[[2,3]])
H1=sh.localHomology(1,toplexes,[[2,3]])
H2=sh.localHomology(2,toplexes,[[2,3]])
   
if H0.shape[1] != 0 or H1.shape[1] != 0 or H2.shape[1] != 1:
    print "Test 8 failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 8b passed"
    
# Test 9: two solid triangles joined along an edge localized to the star over some other edge
toplexes=[[1,2,3],[2,3,4]]
   
H0=sh.localHomology(0,toplexes,[[1,2]],True)
H1=sh.localHomology(1,toplexes,[[1,2]],True)
H2=sh.localHomology(2,toplexes,[[1,2]],True)
   
if H0 != 0 or H1 != 0 or H2 != 0:
    print "Test 9a failed" 
    print "H0="+ str(H0)
    print "H1="+ str(H1)
    print "H2="+ str(H2)
else:
    print "Test 9a passed"
   
H0=sh.localHomology(0,toplexes,[[1,2]])
H1=sh.localHomology(1,toplexes,[[1,2]])
H2=sh.localHomology(2,toplexes,[[1,2]])
   
if H0.shape[1] != 0 or H1.shape[1] != 0 or H2.shape[1] != 0:
    print "Test 9b failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 9b passed"
    
# Test 10: a tetrahedron
toplexes= [[1,2,3],[1,3,4],[1,4,2],[2,4,3]]
    
H0=sh.simplicialHomology(0,toplexes,None,True)
H1=sh.simplicialHomology(1,toplexes,None,True)
H2=sh.simplicialHomology(2,toplexes,None,True)

if H0 != 1 or H1 != 0 or H2 != 1:
    print "Test 10a failed" 
    print "H0="+ str(H0)
    print "H1="+ str(H1)
    print "H2="+ str(H2)
else:
    print "Test 10a passed"
   
H0=sh.simplicialHomology(0,toplexes,None)
H1=sh.simplicialHomology(1,toplexes,None)
H2=sh.simplicialHomology(2,toplexes,None)
    
if H0.shape[1] != 1 or H1.shape[1] != 0 or H2.shape[1] != 1:
    print "Test 10b failed"
    print "H0="+ H0.__repr__() 
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 10b passed"
    
# Test 11: a simple torus
toplexes= [[1,4,2],[2,4,5],[2,5,3],[3,5,6],[1,6,5], \
           [1,2,6],[2,7,6],[2,3,7],[1,7,3],[1,5,7], \
           [4,7,5],[4,6,7],[1,3,4],[3,6,4]] 
   
H0=sh.simplicialHomology(0,toplexes,None,True)
H1=sh.simplicialHomology(1,toplexes,None,True)
H2=sh.simplicialHomology(2,toplexes,None,True)

if H0 != 1 or H1 != 2 or H2 != 1:
    print "Test 11a failed" 
    print "H0="+ str(H0)
    print "H1="+ str(H1)
    print "H2="+ str(H2)
else:
    print "Test 11a passed"

H0=sh.simplicialHomology(0,toplexes,None)
H1=sh.simplicialHomology(1,toplexes,None)
H2=sh.simplicialHomology(2,toplexes,None)
   
if H0.shape[1] != 1 or H1.shape[1] != 2 or H2.shape[1] != 1:
    print "Test 11b failed"
    print "H0="+ H0.__repr__()
    print "H1="+ H1.__repr__()
    print "H2="+ H2.__repr__()
else:
    print "Test 11b passed"
