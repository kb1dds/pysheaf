# Unit test for poset chain methods
#
# Copyright (c) 2013-2014, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 

import pysheaf as ps
import numpy as np

poset=ps.Poset([ps.Cell(0,True,[ps.Coface(1,1),ps.Coface(2,1)]),
                ps.Cell(1,True,[ps.Coface(3,1),ps.Coface(4,1)]),
                ps.Cell(1,True,[ps.Coface(3,1)]),
                ps.Cell(2,True,[ps.Coface(5,1)]),
                ps.Cell(2,True,[]),
                ps.Cell(3,True,[])])

print poset.maximalChains(0)

dg=ps.DirectedGraph([(0,1),(1,2),(1,2),(0,2)])
print dg.maxFlow(4,6)
