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
         ps.SheafCoface(index=2,orientation=-1,restriction=np.matrix(1)),
         ps.SheafCoface(index=3,orientation=1,restriction=np.matrix(1))]), # Vertex
         ps.SheafCell(dimension=0,cofaces=[
         ps.SheafCoface(index=2,orientation=1,restriction=np.matrix(1)),
         ps.SheafCoface(index=3,orientation=-1,restriction=np.matrix(1))]), # Vertex
         ps.SheafCell(dimension=1,stalkDim=1), # Edge
         ps.SheafCell(dimension=1,stalkDim=1)]) # Edge
    print 'Circle Sheaf Betti numbers: ' + str((CircSheaf.cobetti(0),CircSheaf.cobetti(1),CircSheaf.cobetti(2)))

    # A disk
    DiskSheaf=ps.Sheaf([ps.SheafCell(dimension=0,cofaces=[
              ps.SheafCoface(1,-1,np.matrix(1)),
              ps.SheafCoface(1,1,np.matrix(1))]),
              ps.SheafCell(dimension=1,cofaces=[ps.SheafCoface(2,1,np.matrix(1))]),
              ps.SheafCell(dimension=2,stalkDim=1)])
    print 'Disk Sheaf Betti numbers: ' + str((DiskSheaf.cobetti(0),DiskSheaf.cobetti(1),DiskSheaf.cobetti(2)))

    # A 2-sphere
    SphSheaf=ps.Sheaf([ps.SheafCell(dimension=0,cofaces=[
              ps.SheafCoface(1,-1,np.matrix(1)),
              ps.SheafCoface(1,1,np.matrix(1))]),
              ps.SheafCell(dimension=1,cofaces=[
              ps.SheafCoface(2,-1,np.matrix(1)),
              ps.SheafCoface(3,1,np.matrix(1))]),
              ps.SheafCell(dimension=2,stalkDim=1),
              ps.SheafCell(dimension=2,stalkDim=1)])
    print 'S^2 Sheaf Betti numbers: ' + str((SphSheaf.cobetti(0),SphSheaf.cobetti(1),SphSheaf.cobetti(2)))

    MorCirDisk=ps.SheafMorphism([ps.SheafMorphismCell(destinations=[0,1,2],maps=[ps.LinearMorphism(np.matrix(1)),ps.LinearMorphism(np.matrix(1)),ps.LinearMorphism(np.matrix(1))]),
                                 ps.SheafMorphismCell(destinations=[3],maps=[ps.LinearMorphism(np.matrix(1))]),
                                 ps.SheafMorphismCell(destinations=[],maps=[])])
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
    print 'Sheaf Betti numbers for graph with undirected loop: ' + str((LoopSheaf.cobetti(0),LoopSheaf.cobetti(1),LoopSheaf.cobetti(2)))

    # Pushforward Sheaf of collapsed loop
    ColLoopSheaf=ps.Sheaf([ps.SheafCell(dimension=0,cofaces=[
        ps.SheafCoface(1,1,np.matrix([1,1])),
        ps.SheafCoface(2,-1,np.matrix([1,1]))]),
        ps.SheafCell(dimension=1,compactClosure=False,stalkDim=1),
        ps.SheafCell(dimension=1,compactClosure=False,stalkDim=1)])
    print 'Sheaf Betti numbers for graph with collapsed loop: ' + str((ColLoopSheaf.cobetti(0),ColLoopSheaf.cobetti(1),ColLoopSheaf.cobetti(2)))

    # Sheaf on a straight line
    LineSheaf=ps.Sheaf([ps.SheafCell(dimension=0,cofaces=[
        ps.SheafCoface(1,1,np.matrix(1)),
        ps.SheafCoface(2,-1,np.matrix(1))]),
        ps.SheafCell(dimension=1,compactClosure=False,stalkDim=1),
        ps.SheafCell(dimension=1,compactClosure=False,stalkDim=1)])
    print 'Sheaf Betti numbers for straight line: ' + str((LineSheaf.cobetti(0),LineSheaf.cobetti(1),LineSheaf.cobetti(2)))

    # Morphism from collapsed loop to loop
    MorColLoop=ps.SheafMorphism([ps.SheafMorphismCell(destinations=[0,1,2,3],
                                                      maps=[ps.LinearMorphism(np.matrix([[1,0],[0,1]])),
                                                            ps.LinearMorphism(np.matrix([[1,0],[0,1]])),
                                                            ps.LinearMorphism(np.matrix([1,0])),
                                                            ps.LinearMorphism(np.matrix([0,1]))]),
                                 ps.SheafMorphismCell(destinations=[4],maps=[ps.LinearMorphism(np.matrix(1))]),
                                 ps.SheafMorphismCell(destinations=[5],maps=[ps.LinearMorphism(np.matrix(1))])])
                                
    print 'degree 0 induced map Collapsed->Loop: ' + str(ps.inducedMap(ColLoopSheaf,LoopSheaf,MorColLoop,0))
    print 'degree 1 induced map Collapsed->Loop: ' + str(ps.inducedMap(ColLoopSheaf,LoopSheaf,MorColLoop,1))
    print 'degree 2 induced map Collapsed->Loop: ' + str(ps.inducedMap(ColLoopSheaf,LoopSheaf,MorColLoop,2))
                
    # Morphism from collapsed loop to straight line
    MorColLine=ps.SheafMorphism([ps.SheafMorphismCell(destinations=[0],maps=[ps.LinearMorphism(np.matrix([1,1]))]),
                                 ps.SheafMorphismCell(destinations=[1],maps=[ps.LinearMorphism(np.matrix(1))]),
                                 ps.SheafMorphismCell(destinations=[2],maps=[ps.LinearMorphism(np.matrix(1))])])

    print 'degree 0 induced map Collapsed->Line: ' + str(ps.inducedMap(ColLoopSheaf,LineSheaf,MorColLine,0))
    print 'degree 1 induced map Collapsed->Line: ' + str(ps.inducedMap(ColLoopSheaf,LineSheaf,MorColLine,1))
    print 'degree 2 induced map Collapsed->Line: ' + str(ps.inducedMap(ColLoopSheaf,LineSheaf,MorColLine,2))

    PerSheaf0=ps.PersistenceSheaf([ColLoopSheaf,LineSheaf,LoopSheaf],[(0,1,MorColLine),(0,2,MorColLoop)],0)

    print '0-Persistent degree 0 Sheaf Betti number ' + str(PerSheaf0.cobetti(0))
    print '0-Persistent degree 1 Sheaf Betti number ' + str(PerSheaf0.cobetti(1))

    PerSheaf1=ps.PersistenceSheaf([ColLoopSheaf,LineSheaf,LoopSheaf],[(0,1,MorColLine),(0,2,MorColLoop)],1)

    print '1-Persistent degree 0 Sheaf Betti number ' + str(PerSheaf1.cobetti(0))
    print '1-Persistent degree 1 Sheaf Betti number ' + str(PerSheaf1.cobetti(1))

    fs=ps.FlowSheaf(ps.DirectedGraph([(None,1),(None,1),(1,2),(1,None),(None,2),(2,None)]))
    print 'Flow sheaf degree 0 Sheaf Betti number ' + str(fs.cobetti(0))
    
    fs2=fs.star([6])
    print 'Flow sheaf degree 0 Sheaf Betti number ' + str(fs2.cobetti(0))

    fs3=fs.star([7])
    print 'Flow sheaf degree 0 Sheaf Betti number ' + str(fs3.cobetti(0))

    print 'Dimension of local sections over a single edge ' + str(fs.localSectional([0])[0].cobetti(0))
    print 'Dimension of local sections over two edges ' + str(fs.localSectional([0,1])[0].cobetti(0))
    print 'Dimension of local sections over three edges ' + str(fs.localSectional([0,1,2])[0].cobetti(0))
    print 'Dimension of local sections over three edges and a common vertex ' + str(fs.localSectional([0,1,2,6])[0].cobetti(0))
    
    print 'Induced map on local sections from 3 edges to 2 edges ' + str(fs.localRestriction([0,1,2],[0,1]))
    print 'Induced map on local sections from 3 edges and common vertex to 3 edges ' + str(fs.localRestriction([0,1,2,6],[0,1,2]))

    pos=ps.Poset([ps.Cell(dimension=0,cofaces=[ps.Coface(index=1,orientation=1),ps.Coface(index=2,orientation=1)]),
                  ps.Cell(dimension=1,cofaces=[ps.Coface(index=3,orientation=1)]),
                  ps.Cell(dimension=1,cofaces=[ps.Coface(index=3,orientation=1)]),
                  ps.Cell(dimension=2,cofaces=[])])
                  
    dg=pos.hasseDiagram()
