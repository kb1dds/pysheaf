# Try to extend sections, based on examples found in Chapter 3 of
# "Topological Signal Processing", by Michael Robinson
#
# Copyright (c) 2013-2014, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 


import numpy as np
import pysheaf as ps

sh1=ps.Sheaf([ps.SheafCell(dimension=0,cofaces=[ps.SheafCoface(index=2,orientation=-1,restriction=np.matrix(1))]),
              ps.SheafCell(dimension=0,cofaces=[ps.SheafCoface(index=2,orientation=1,restriction=np.matrix(1))]),
              ps.SheafCell(dimension=1,cofaces=[],stalkDim=1)])

sec1=ps.Section([ps.SectionCell(0,1),ps.SectionCell(1,1)])
sec2=ps.Section([ps.SectionCell(0,1)])
sec3=ps.Section([ps.SectionCell(0,1),ps.SectionCell(1,2)])

# Extending section over two vertices to a common coface
if sec1.extend(sh1,2) and sec2.extend(sh1,2) and not sec3.extend(sh1,2) and sec2.extend(sh1,1):
    print "Test 1 passed"
else:
    print "Test 1 failed"

sh2=ps.Sheaf([ps.SheafCell(dimension=0,cofaces=[ps.SheafCoface(index=1,orientation=-1,restriction=np.matrix(1))]),
              ps.SheafCell(dimension=1,cofaces=[],stalkDim=1),
              ps.SheafCell(dimension=0,cofaces=[ps.SheafCoface(index=1,orientation=1,restriction=np.matrix(-1)),
                                                ps.SheafCoface(index=3,orientation=-1,restriction=np.matrix(1))]),
              ps.SheafCell(dimension=1,cofaces=[],stalkDim=1)])
sec4=ps.Section([ps.SectionCell(0,1)])

# Extending along a line
if sec4.extend(sh2,1) and sec4.extend(sh2,2) and sec4.extend(sh2,3):
    print [s.value for s in sec4.sectionCells]
    print "Test 2 passed"
else:
    print "Test 2 failed"

# Mayer-Vietoris example from tspbook
sh3=ps.FlowSheaf(ps.DirectedGraph([(None,0),(None,0),(0,None),(0,1),(None,1),(1,None)]))
sec5=ps.Section([ps.SectionCell(0,3),
                 ps.SectionCell(1,2),
                 ps.SectionCell(2,1),
                 ps.SectionCell(4,3),
                 ps.SectionCell(5,8)])
sec6=ps.Section([ps.SectionCell(0,3),
                 ps.SectionCell(1,2),
                 ps.SectionCell(2,1),
                 ps.SectionCell(4,3),
                 ps.SectionCell(5,7)])
sec7=ps.Section([ps.SectionCell(1,2),
                 ps.SectionCell(4,3)])
sec8=ps.Section([ps.SectionCell(1,2),
                 ps.SectionCell(4,3)])

# Trival extension to vertices
if sec5.extend(sh3,6) and sec6.extend(sh3,6) and sec5.extend(sh3,7) and sec6.extend(sh3,7):
    print "Test 3 passed"
else:
    print "Test 3 failed"

# Extending to a connecting edge (using the previously extended sections)
if not sec5.extend(sh3,3) and sec6.extend(sh3,3):
    print "Test 4 passed"
else:
    print "Test 4 failed"

# Extending a very underdetermined section
if sec7.extend(sh3,6) and sec7.extend(sh3,7) and sec7.extend(sh3,3):
    print "Test 5 passed"
else:
    print "Test 5 failed"

# Extending with specific values in mind; in which extension works ...
if sec7.extend(sh3,0,3) and sec7.extend(sh3,2,1) and sec7.extend(sh3,5,7) and sec7.extend(sh3,6) and sec7.extend(sh3,7) and sec7.extend(sh3,3):
    print "Test 6 passed"
else:
    print "Test 6 failed"

# ... and in which a section does not completely extend
if sec8.extend(sh3,0,3) and sec8.extend(sh3,2,1) and sec8.extend(sh3,5,8) and sec8.extend(sh3,6) and sec8.extend(sh3,7) and not sec8.extend(sh3,3):
    print "Test 7 passed"
else:
    print "Test 7 failed"
