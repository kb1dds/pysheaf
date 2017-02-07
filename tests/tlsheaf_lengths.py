# Test of using transmission line sheaves for geometry sounding as
# found in "Topological Signal Processing", by Michael Robinson
#
# Copyright (c) 2013-2014, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 


import pysheaf as ps
import numpy as np

searchstart=100
searchstop=1100
tol=1e-5

# Simulated geometry
ft2m=12*0.0254 # Conversion feet to meters
dg=ps.DirectedGraph([(None,1),(None,2),(None,3),(1,2),(2,3),(3,1)])
dg.cells[3].length=70*ft2m # Edge e3
dg.cells[4].length=(150+70)*ft2m # Edge e1
dg.cells[5].length=150*ft2m # Edge e2

frequencies=[905e6,2451e6] # Operating frequencies
wavenumbers=[2*np.pi*f/3e8 for f in frequencies]

sheaves=[ps.TransLineSheaf(dg,k) for k in wavenumbers]

L1s=[]
L2s=[]
L3s=[]
for tl,k in zip(sheaves,wavenumbers):
    print "Frequency = " + str(k*3e8/(2*np.pi)/1e6) + " MHz"

    sec=np.sum(tl.cohomology(0),1)

    print "Measurements: "
    print "v1->external: " + str(10*np.log10(np.abs(sec[0]))) + " dB " + str(np.angle(sec[0])*180/np.pi) + " degrees"
    print "v1->e2: " + str(10*np.log10(np.abs(sec[2]))) + " dB " + str(np.angle(sec[2])*180/np.pi) + " degrees"
    print "v1->e3: " + str(10*np.log10(np.abs(sec[1]))) + " dB " + str(np.angle(sec[1])*180/np.pi) + " degrees"
    print "v2->external: " + str(10*np.log10(np.abs(sec[3]))) + " dB " + str(np.angle(sec[3])*180/np.pi) + " degrees"
    print "v2->e1: " + str(10*np.log10(np.abs(sec[5]))) + " dB " + str(np.angle(sec[5])*180/np.pi) + " degrees"
    print "v2->e3: " + str(10*np.log10(np.abs(sec[4]))) + " dB " + str(np.angle(sec[4])*180/np.pi) + " degrees"
    print "v3->external: " + str(10*np.log10(np.abs(sec[6]))) + " dB " + str(np.angle(sec[6])*180/np.pi) + " degrees"
    print "v3->e1: " + str(10*np.log10(np.abs(sec[7]))) + " dB " + str(np.angle(sec[7])*180/np.pi) + " degrees"
    print "v3->e2: " + str(10*np.log10(np.abs(sec[8]))) + " dB " + str(np.angle(sec[8])*180/np.pi) + " degrees"
     
    # Computed edge lengths
    L1=np.real(-1j/k*np.log(tl.cells[7].cofaces[2].restriction(sec[3:6])[0,0]/sec[7]))
    L2=np.real(1j/k*np.log(tl.cells[6].cofaces[2].restriction(sec[0:3])[0,1]/sec[8]))
    L3=np.real(-1j/k*np.log(tl.cells[6].cofaces[1].restriction(sec[0:3])[0,0]/sec[4]))

    print "L1 = " + str(L1/ft2m) + " ft, which is off by " + str((L1-(150+70)*ft2m)*k/np.pi) + " pi"
    print "L2 = " + str(L2/ft2m) + " ft, which is off by " + str((L2-150*ft2m)*k/np.pi) + " pi"
    print "L3 = " + str(L3/ft2m) + " ft, which is off by " + str((L3-70*ft2m)*k/np.pi) + " pi"

    L1s.append(L1)
    L2s.append(L2)
    L3s.append(L3)

# Search for compatible lengths
# for which the following equation holds for integers n1, n2, ...
#  L1 - n1 * pi / k1 = L2 - n2 * pi / k2 = ...
# We suppose that the true length is this common value

L1_merge=[(L1s[0]+n*np.pi/wavenumbers[0],m,n)
          for m in range(searchstart,searchstop) 
          for n in range(searchstart,searchstop) if np.abs(L1s[0]+n*np.pi/wavenumbers[0]-L1s[1]-m*np.pi/wavenumbers[1])<tol]

L2_merge=[(L2s[0]+n*np.pi/wavenumbers[0],m,n) 
          for m in range(searchstart,searchstop) 
          for n in range(searchstart,searchstop) if np.abs(L2s[0]+n*np.pi/wavenumbers[0]-L2s[1]-m*np.pi/wavenumbers[1])<tol]

L3_merge=[(L3s[0]+n*np.pi/wavenumbers[0],m,n) 
          for m in range(searchstart,searchstop) 
          for n in range(searchstart,searchstop) if np.abs(L3s[0]+n*np.pi/wavenumbers[0]-L3s[1]-m*np.pi/wavenumbers[1])<tol]

print "Minimal merged values:"
print "L1 = " + str(L1_merge[0][0]/ft2m) + " ft"
print "L2 = " + str(L2_merge[0][0]/ft2m) + " ft"
print "L3 = " + str(L3_merge[0][0]/ft2m) + " ft"
