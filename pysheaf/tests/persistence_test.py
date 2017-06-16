# Sample persistence sheaf calculation
#
# Copyright (c) 2013-2014, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 

import numpy as np
import pysheaf as ps

if __name__ == '__main__':
    targ1=ps.DirectedGraph([(None,1),(1,2),(2,None),(None,3),(3,4),(4,None)])
    fs1=ps.FlowSheaf(targ1)
    targ2=ps.DirectedGraph([(None,1),(None,1),(1,2),(2,None),(2,None)])
    targ3=ps.DirectedGraph([(None,1),(None,1),(1,2),(2,None)])
    
    map1=[(0,0),(1,2),(2,3),(3,1),(4,2),(5,4),(6,5),(7,6),(8,5),(9,6)]
    map2=[(0,0),(1,1),(2,2),(3,3),(4,3),(5,4)]
    
    print 'Ready to compute!'

    pf1,pf1m=fs1.pushForward(targ2,map1)
    print 'pushforward FS 1 induced map ' + str(ps.inducedMap(pf1,fs1,pf1m,0))
    fs2,fsm2=pf1.flowCollapse()
    print 'collapse 1 induced map ' + str(ps.inducedMap(pf1,fs2,fsm2,0))
    pf2,pf2m=fs2.pushForward(targ3,map2)
    print 'pushforward FS 2 induced map ' + str(ps.inducedMap(pf2,fs2,pf2m,0))
    fs3,fsm3=pf2.flowCollapse()
    print 'collapse 2 induced map ' + str(ps.inducedMap(pf2,fs3,fsm3,0))

    persh=ps.PersistenceSheaf([fs1,pf1,fs2,pf2,fs3],[(1,0,pf1m),(1,2,fsm2),(3,2,pf2m),(3,4,fsm3)],0)
    
    print 'Persistence Sheaf Betti 0=' + str(persh.cobetti(0))
