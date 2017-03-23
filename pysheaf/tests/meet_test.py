# Test generation of the meet matrix for a poset
#
# Copyright (c) 2013-2014, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 

import pysheaf as ps

poset=ps.Poset([ps.Cell(0,cofaces=[ps.Coface(2,1)]),ps.Cell(0,cofaces=[ps.Coface(2,1)]),ps.Cell(1,[])])

poset.meetMatrix()

poset2=ps.Poset([ps.Cell(0,cofaces=[ps.Coface(1,1),ps.Coface(2,1)]),ps.Cell(1,[]),ps.Cell(1,cofaces=[ps.Coface(3,1),ps.Coface(4,1)]),ps.Cell(2,[]),ps.Cell(2,[])])
poset2.meetMatrix()
