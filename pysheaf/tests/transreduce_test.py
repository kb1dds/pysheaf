# Transitive reduction for posets
#
# Copyright (c) 2013-2014, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 

import pysheaf as ps

complex1=ps.Poset([ps.Cell(dimension=0,cofaces=[ps.Coface(index=1,orientation=1),ps.Coface(index=2,orientation=1),ps.Coface(index=3,orientation=1)]),
                   ps.Cell(dimension=3,cofaces=[ps.Coface(index=2,orientation=1),ps.Coface(index=3,orientation=1)]),
                   ps.Cell(dimension=1,cofaces=[ps.Coface(index=3,orientation=1)]),
                   ps.Cell(dimension=2,cofaces=[])])

complex1.transitiveReduce()
for c in complex1.cells:
    print "Dimension: " + str(c.dimension) + ", cofaces: " + str([cf.index for cf in c.cofaces]) 
