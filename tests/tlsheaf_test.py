# Test transmission line sheaves
#
# Copyright (c) 2013-2014, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 

import pysheaf as ps
import numpy as np

f=900e6 # Operating frequency
ft2m=12*0.0254 # Conversion feet to meters
wavenumber=2*np.pi*f/3e8

dg=ps.DirectedGraph([(None,1),(None,2),(None,3),(1,2),(2,3),(3,1)])
dg.cells[3].length=70*ft2m # Edge e3
dg.cells[4].length=(150+70)*ft2m # Edge e1
dg.cells[5].length=150*ft2m # Edge e2
tl=ps.TransLineSheaf(dg,wavenumber)

print tl.cohomology(0)

sec=np.sum(tl.cohomology(0),1)

if np.allclose(tl.cells[6].cofaces[1].restriction(sec[0:3]),tl.cells[7].cofaces[1].restriction(sec[3:6])) and np.allclose(tl.cells[6].cofaces[2].restriction(sec[0:3]),tl.cells[8].cofaces[2].restriction(sec[6:9])) and np.allclose(tl.cells[7].cofaces[2].restriction(sec[3:6]),tl.cells[8].cofaces[1].restriction(sec[6:9])):
    print "Test passed"
else:
    print "Restriction to an edge is mismatched"

print "Magnitudes (dB)"
print np.log10(np.abs(sec))*10

print "Angles (degrees)"
print np.angle(sec)*180/np.pi

# Computed edge lengths
L1=-1j/wavenumber*np.log(tl.cells[7].cofaces[2].restriction(sec[3:6])[0,0]/sec[7])/ft2m
L2=1j/wavenumber*np.log(tl.cells[6].cofaces[2].restriction(sec[0:3])[0,1]/sec[8])/ft2m
L3=-1j/wavenumber*np.log(tl.cells[6].cofaces[1].restriction(sec[0:3])[0,0]/sec[4])/ft2m

print "L1 = " + str(np.real(L1)) + " ft, which is off by " + str(np.real((L1*wavenumber*ft2m-(150+70)*wavenumber*ft2m)/np.pi)) + " pi"
print "L2 = " + str(np.real(L2)) + " ft, which is off by " + str(np.real((L2*wavenumber*ft2m-150*wavenumber*ft2m)/np.pi)) + " pi"
print "L3 = " + str(np.real(L3)) + " ft, which is off by " + str(np.real((L3*wavenumber*ft2m-70*wavenumber*ft2m)/np.pi)) + " pi"
