# Unit test for Local Simplicial Homology: Annulus
#
# Copyright (c) 2016-2017, Michael Robinson, Chris Capraro
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 

import numpy as np
import time as tm
import matplotlib.pyplot as plt
import simplicialHomology as sh
import plotComplex as pc

# Example: An annulus
locations_annulus=[]
cplx=[]
radii=np.linspace(1,10,10)
angles=np.linspace(0,2*np.pi,10)[0:-1]
for i,rad in enumerate(np.nditer(radii)):
    for j,angle in enumerate(np.nditer(angles)):
        locations_annulus.append(np.array([rad*np.cos(angle),rad*np.sin(angle)]))
        if i > 0:
            if j == len(angles)-1:
                cplx+=[[i*len(angles)+j,i*len(angles),(i-1)*len(angles)+j],
                       [(i-1)*len(angles),(i-1)*len(angles)+j,i*len(angles)]]
            else:
                cplx+=[[i*len(angles)+j,i*len(angles)+j+1,(i-1)*len(angles)+j],
                       [(i-1)*len(angles)+j+1,(i-1)*len(angles)+j,i*len(angles)+j+1]]


startTime = tm.time()

cplx_annulus=cplx+sh.ksimplices(cplx,1)+sh.ksimplices(cplx,0)

print 'colors_annulus_1'
colors_annulus_1=[sh.localHomology(1,cplx,[spx],True) for spx in cplx_annulus]
print 'colors_annulus_2'
colors_annulus_2=[sh.localHomology(2,cplx,[spx],True) for spx in cplx_annulus]

endTime = tm.time()
print "\n**** Total time = %f s ****\n" % (endTime-startTime)

plt.figure()
plt.hold(True)
plt.subplot(121)
pc.plot_complex(locations_annulus,cplx_annulus,colors_annulus_1)
plt.subplot(122)
pc.plot_complex(locations_annulus,cplx_annulus,colors_annulus_2)
plt.savefig('annulus.png')
