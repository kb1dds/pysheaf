# Unit test for persistence sheaf library
#
# Copyright (c) 2013-2014, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 

import numpy as np
import pysheaf as ps

if __name__ == '__main__':
    # A circle
    CircSheaf=ps.Sheaf([ps.SheafCell(dimension=0,cofaces=[
         ps.SheafCoface(index=2,orientation=-1,corestriction=np.matrix(1)),
         ps.SheafCoface(index=3,orientation=1,corestriction=np.matrix(1))]), # Vertex
         ps.SheafCell(dimension=0,cofaces=[
         ps.SheafCoface(index=2,orientation=1,corestriction=np.matrix(1)),
         ps.SheafCoface(index=3,orientation=-1,corestriction=np.matrix(1))]), # Vertex
         ps.SheafCell(dimension=1,stalkDim=1), # Edge
         ps.SheafCell(dimension=1,stalkDim=1)]) # Edge
    print 'Circle Betti numbers: ' + str((CircSheaf.betti(0),CircSheaf.betti(1),CircSheaf.betti(2)))

    # A disk
    DiskSheaf=ps.Sheaf([ps.SheafCell(dimension=0,cofaces=[
              ps.SheafCoface(1,-1,np.matrix(1)),
              ps.SheafCoface(1,1,np.matrix(1))]),
              ps.SheafCell(dimension=1,cofaces=[ps.SheafCoface(2,1,np.matrix(1))]),
              ps.SheafCell(dimension=2,stalkDim=1)])
    print 'Disk Betti numbers: ' + str((DiskSheaf.betti(0),DiskSheaf.betti(1),DiskSheaf.betti(2)))

    # A 2-sphere
    SphSheaf=ps.Sheaf([ps.SheafCell(dimension=0,cofaces=[
              ps.SheafCoface(1,-1,np.matrix(1)),
              ps.SheafCoface(1,1,np.matrix(1))]),
              ps.SheafCell(dimension=1,cofaces=[
              ps.SheafCoface(2,-1,np.matrix(1)),
              ps.SheafCoface(3,1,np.matrix(1))]),
              ps.SheafCell(dimension=2,stalkDim=1),
              ps.SheafCell(dimension=2,stalkDim=1)])
    print 'S^2 Betti numbers: ' + str((SphSheaf.betti(0),SphSheaf.betti(1),SphSheaf.betti(2)))

    MorCirDisk=[ps.SheafMorphismCell(destinations=[0,1,2],maps=[np.matrix(1),np.matrix(1),np.matrix(1)]),
                ps.SheafMorphismCell(destinations=[3],maps=[np.matrix(1)]),
                ps.SheafMorphismCell(destinations=[],maps=[])]
    print 'degree 0 induced map S^1->D^2: ' + str(ps.inducedMap(DiskSheaf,CircSheaf,MorCirDisk,0))
    print 'degree 1 induced map S^1->D^2:' + str(ps.inducedMap(DiskSheaf,CircSheaf,MorCirDisk,1))
    print 'degree 2 induced map S^1->D^2:' + str(ps.inducedMap(DiskSheaf,CircSheaf,MorCirDisk,2))

    # Sheaf over a graph with an undirected loop
    LoopSheaf=ps.Sheaf([ps.SheafCell(dimension=0,cofaces=[
        ps.SheafCoface(2,1,np.matrix([1,0])),
        ps.SheafCoface(3,1,np.matrix([0,0])),
        ps.SheafCoface(4,-1,np.matrix([1,1]))]),
        ps.SheafCell(dimension=0,cofaces=[
        ps.SheafCoface(2,-1,np.matrix([1,0])),
        ps.SheafCoface(3,-1,np.matrix([0,1])),
        ps.SheafCoface(5,1,np.matrix([1,1]))]),
        ps.SheafCell(dimension=1,stalkDim=1),
        ps.SheafCell(dimension=1,stalkDim=1),
        ps.SheafCell(dimension=1,compactClosure=False,stalkDim=1),
        ps.SheafCell(dimension=1,compactClosure=False,stalkDim=1)])
    print 'Betti numbers for graph with undirected loop: ' + str((LoopSheaf.betti(0),LoopSheaf.betti(1),LoopSheaf.betti(2)))

    # Pushforward Sheaf of collapsed loop
    ColLoopSheaf=ps.Sheaf([ps.SheafCell(dimension=0,cofaces=[
        ps.SheafCoface(1,1,np.matrix([1,1])),
        ps.SheafCoface(2,-1,np.matrix([1,1]))]),
        ps.SheafCell(dimension=1,compactClosure=False,stalkDim=1),
        ps.SheafCell(dimension=1,compactClosure=False,stalkDim=1)])
    print 'Betti numbers for graph with collapsed loop: ' + str((ColLoopSheaf.betti(0),ColLoopSheaf.betti(1),ColLoopSheaf.betti(2)))

    # Sheaf on a straight line
    LineSheaf=ps.Sheaf([ps.SheafCell(dimension=0,cofaces=[
        ps.SheafCoface(1,1,np.matrix(1)),
        ps.SheafCoface(2,-1,np.matrix(1))]),
        ps.SheafCell(dimension=1,compactClosure=False,stalkDim=1),
        ps.SheafCell(dimension=1,compactClosure=False,stalkDim=1)])
    print 'Betti numbers for straight line: ' + str((LineSheaf.betti(0),LineSheaf.betti(1),LineSheaf.betti(2)))

    # Morphism from collapsed loop to loop
    MorColLoop=[ps.SheafMorphismCell(destinations=[0,1,2,3],
                                     maps=[np.matrix([[1,0],[0,1]]),
                                           np.matrix([[1,0],[0,1]]),
                                           np.matrix([1,0]),
                                           np.matrix([0,1])]),
                ps.SheafMorphismCell(destinations=[4],maps=[np.matrix(1)]),
                ps.SheafMorphismCell(destinations=[5],maps=[np.matrix(1)])]

    print 'degree 0 induced map Collapsed->Loop: ' + str(ps.inducedMap(ColLoopSheaf,LoopSheaf,MorColLoop,0))
    print 'degree 1 induced map Collapsed->Loop: ' + str(ps.inducedMap(ColLoopSheaf,LoopSheaf,MorColLoop,1))
    print 'degree 2 induced map Collapsed->Loop: ' + str(ps.inducedMap(ColLoopSheaf,LoopSheaf,MorColLoop,2))
                
    # Morphism from collapsed loop to straight line
    MorColLine=[ps.SheafMorphismCell(destinations=[0],maps=[np.matrix([1,1])]),
                ps.SheafMorphismCell(destinations=[1],maps=[np.matrix(1)]),
                ps.SheafMorphismCell(destinations=[2],maps=[np.matrix(1)])]

    print 'degree 0 induced map Collapsed->Line: ' + str(ps.inducedMap(ColLoopSheaf,LineSheaf,MorColLine,0))
    print 'degree 1 induced map Collapsed->Line: ' + str(ps.inducedMap(ColLoopSheaf,LineSheaf,MorColLine,1))
    print 'degree 2 induced map Collapsed->Line: ' + str(ps.inducedMap(ColLoopSheaf,LineSheaf,MorColLine,2))

    PerSheaf0=ps.PersistenceSheaf([ColLoopSheaf,LineSheaf,LoopSheaf],[(0,1,MorColLine),(0,2,MorColLoop)],0)

    print '0-Persistent degree 0 Betti number ' + str(PerSheaf0.betti(0))
    print '0-Persistent degree 1 Betti number ' + str(PerSheaf0.betti(1))

    PerSheaf1=ps.PersistenceSheaf([ColLoopSheaf,LineSheaf,LoopSheaf],[(0,1,MorColLine),(0,2,MorColLoop)],1)

    print '1-Persistent degree 0 Betti number ' + str(PerSheaf1.betti(0))
    print '1-Persistent degree 1 Betti number ' + str(PerSheaf1.betti(1))

    fs=ps.FlowSheaf(ps.DirectedGraph([(None,1),(None,1),(1,2),(1,None),(None,2),(2,None)]))
    print 'Flow sheaf degree 0 Betti number ' + str(fs.betti(0))
    
    fs2=fs.star([6])
    print 'Flow sheaf degree 0 Betti number ' + str(fs2.betti(0))

    fs3=fs.star([7])
    print 'Flow sheaf degree 0 Betti number ' + str(fs3.betti(0))

    print 'Dimension of local sections over a single edge ' + str(fs.localSectional([0])[0].betti(0))
    print 'Dimension of local sections over two edges ' + str(fs.localSectional([0,1])[0].betti(0))
    print 'Dimension of local sections over three edges ' + str(fs.localSectional([0,1,2])[0].betti(0))
    print 'Dimension of local sections over three edges and a common vertex ' + str(fs.localSectional([0,1,2,6])[0].betti(0))
    
    print 'Induced map on local sections from 3 edges to 2 edges ' + str(fs.localRestriction([0,1,2],[0,1]))
    print 'Induced map on local sections from 3 edges and common vertex to 3 edges ' + str(fs.localRestriction([0,1,2,6],[0,1,2]))

    pos=ps.Poset([ps.Cell(dimension=0,cofaces=[ps.Coface(index=1,orientation=1),ps.Coface(index=2,orientation=1)]),
                  ps.Cell(dimension=1,cofaces=[ps.Coface(index=3,orientation=1)]),
                  ps.Cell(dimension=1,cofaces=[ps.Coface(index=3,orientation=1)]),
                  ps.Cell(dimension=2,cofaces=[])])
                  
    dg=pos.hasseDiagram()
