# Unit test for Local Simplicial Homology: A thickened random graph
#
# Copyright (c) 2016-2017, Michael Robinson, Chris Capraro
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 

import numpy as np
import time as tm
import numpy.random as random
import matplotlib.pyplot as plt
import simplicialHomology as sh
import plotComplex as pc

# Example: A "thick" graph
startTime = tm.time()
random.seed(100)

locations=[]
while len(locations)<100:
    loc=random.rand(2)
    if np.sqrt((loc[0]-0.5)**2+(loc[1]-0.5)**2)<0.5:
        locations.append(loc)
cplx=sh.vietorisRips(locations,0.1,maxdim=2)
cplx_grph=sh.ksimplices(cplx,0)+sh.ksimplices(cplx,1)+sh.ksimplices(cplx,2)
colors_grph_1=[sh.localHomology(1,cplx,[spx],True) for spx in cplx_grph]
colors_grph_2=[sh.localHomology(2,cplx,[spx],True) for spx in cplx_grph]

endTime = tm.time()
print "\n**** Total time = %f s ****\n" % (endTime-startTime)

plt.figure()
plt.hold(True)
plt.subplot(121)
pc.plot_complex(locations,cplx_grph,colors_grph_1)
plt.subplot(122)
pc.plot_complex(locations,cplx_grph,colors_grph_2)
plt.savefig('thick_graph.png')
