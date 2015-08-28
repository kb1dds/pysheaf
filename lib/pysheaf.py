# Persistence-capable sheaf manipulation library
#
# Copyright (c) 2013-2015, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 

import numpy as np
import random
try:
    import networkx as nx
except:
    pass

## Data structures
class Coface:
    """A coface relation"""
    def __init__(self,index,orientation):
        self.index=index
        self.orientation=orientation
        
class Cell:
    """A cell in a cell complex"""
    def __init__(self,dimension,compactClosure=True,cofaces=[]):
        self.dimension=dimension
        self.compactClosure=compactClosure
        self.cofaces=cofaces

    def cofaceList(self):
        """Return the list of cofaces of this cell"""
        return [cf.index for cf in self.cofaces]    
    
    def isCoface(self,index,orientation=None):
        """Check if a given cell index is a coface of this cell. Optionally check that the orientation is as given."""
        if orientation==None:
            return index in [cf.index for cf in self.cofaces]
        else:
            return (index,orientation) in [(cf.index,cf.orientation) for cf in self.cofaces]
        
class CellComplex:
    def __init__(self,cells):
        """Construct a cell complex from its cells"""
        self.cells=cells
        
    def isFaceOf(self,c,cells=[]):
        """Construct a list of all cells that a given cell is a face of"""
        if cells:
            cl=cells
        else:
            cl=range(len(self.cells))
        return [i for i in cl if self.cells[i].isCoface(c)]
    
    def skeleton(self,k,compactSupport=False):
        """Return the k-skeleton of a cell complex.  Optionally ensure that the complex returned is closed."""
        return [i for i in range(len(self.cells)) 
                if ((compactSupport or self.cells[i].compactClosure) and self.cells[i].dimension==k)]
                
    def faces(self,c):
        """Compute a list of all faces of a cell"""
        return [i for i in range(len(self.cells)) if self.cells[i].isCoface(c)]

    def cofaces(self,c,cells=[]):
        """Iterate over cofaces (of all dimensions) of a given cell; optional argument specifies which cells are permissible cofaces"""
        for cf in self.cells[c].cofaces:
            if cf.index in cells or not cells:
                yield cf
                
        for cf in self.cells[c].cofaces:
            if cf.index in cells or not cells:
                self.cofaces(cf.index,cells)

    def components(self,cells=[]):
        """Compute connected components; optional argument specifies permissible cells"""
        if not cells:
            cellsleft=range(len(self.cells))
        else:
            cellsleft=cells

        cpts=[]
        cpt=[]
        while cellsleft:
            cpt=self.expandComponent(cellsleft[0],cells,[cellsleft[0]])
            cpts+=[list(set(cpt))]
            cellsleft=list(set(cellsleft).difference(cpt))
        return cpts

    def expandComponent(self,start,cells=[],current_cpt=[]):
        """Compute the connected component started from a given cell.  Optional argument specifies permissible cells"""
        if not cells:
            cellsleft=list(set(range(len(self.cells))).difference(current_cpt))
        else:
            cellsleft=list(set(cells).difference(current_cpt))
	if not cellsleft:
            return current_cpt

        neighbors=self.connectedTo(start,cellsleft)
        for c in neighbors:
            current_cpt+=self.expandComponent(c,cellsleft,list(set(current_cpt+neighbors)))
            current_cpt=list(set(current_cpt))
        return current_cpt

    def connectedTo(self,start,cells=[]):
        """Which cells is a cell connected to? Optional argument specifies permissible cells"""
        if not cells:
            return list(set(self.faces(start) + self.cells[start].cofaceList()))
        else:
	    return list(set.intersection(set(self.faces(start) + self.cells[start].cofaceList()),cells))
    
    def starCells(self,cells):
        """Cells in star over a subset of a cell complex"""
        return list(set(cells+[cf.index for c in cells for cf in self.cofaces(c)]))

    def homology(self,k,subcomplex=None):
        """Compute (relative) homology of the cell complex"""
        pass

    def localHomology(self,k,cell):
        """Compute local homology localized at a (star over a) cell"""
        pass

    def inducedMapLocalHomology(self,k,cell1,cell2):
        """Compute the induced map on local homology between two cells.  It is assumed that cell2 is a coface (perhaps not dimension 1) of cell1"""
        pass

    def attachDiagram(self):
        """Draw the attachment diagram using NetworkX"""
        G=nx.DiGraph()

        G.add_nodes_from(range(len(self.cells)))

        G.add_edges_from([(i,cf.index) 
                          for i in range(len(self.cells)) 
                          for cf in self.cells[i].cofaces])

        return G

class Poset(CellComplex):
    def hasseDiagram(self):
        """Convert a poset to a directed graph"""

        # Internal edges
        graph=[(i,cf.index) for i in range(len(self.cells)) 
               for cf in self.cells[i].cofaces]

        # Attach external inputs to minimal elements
        graph += [(None,i) for i in range(len(self.cells)) 
                  if not self.faces(i)]

        # Attach external outputs to maximal elements
        graph += [(i,None) for i in range(len(self.cells)) 
                  if not self.cells[i].cofaces]

        return DirectedGraph(graph)

    def transitiveReduce(self):    
        """Remove coface relations that are redundant
        From Harry Hsu. "An algorithm for finding a minimal equivalent graph of a digraph.", Journal of the ACM, 22(1):11-16, January 1975"""

        for c1 in self.cells: # i
            for cf in c1.cofaces:
                c2=self.cells[cf.index] # j
                for cf2 in c2.cofaces: # k
                    c1.cofaces=[cf3 for cf3 in c1.cofaces 
                                if cf3.index != cf2.index]

    def maximalChains(self,start,history=[]):
        """Compute a list of maximal chains beginning at a cell"""
        if self.cells[start].cofaces:
            chains=[self.maximalChains(cf.index,history+[start]) 
                    for cf in self.cells[start].cofaces]
            lst=[]
            for ch in chains:
                lst+=ch
            return lst
        else:
            return [history+[start]]

    def maxLenChains(self,start):
        """Compute the list of chains of maximal edge length starting at a cell"""
        chains=self.maximalChains(start)
        mx=max([len(ch) for ch in chains])
        return [ch for ch in chains if len(ch)==mx]

    def mobius(self,x,y):
        """Compute the Mobius function between two elements"""
        if x==y:
            return 1

        if y not in [cf.index for cf in self.cofaces(x)]:
            return 0

        mb=-1
        for cf in self.cofaces(x):
            z=cf.index
            if z != y and y in [cf2.index for cf2 in self.cofaces(z)]:
                mb-=self.mobius(x,z)

        return mb
        
    def meet(self,c1,c2):
        """Compute the meet of two elements, if it exists"""
        if c1==c2:
            return c1
        if self.cells[c1].cofaces:
            for cf1 in self.cofaces(c1):
                if cf1.index==c2:
                    return c2
                for cf2 in self.cofaces(c2):
                    if cf2.index==c1:
                        return c1
                    if cf1.index==cf2.index:
                        return cf1.index
        else:
            for cf2 in self.cofaces(c2):
                if cf2.index==c1:
                    return c1
                for cf1 in self.cofaces(c1):
                    if cf1.index==c2:
                        return c2
                    if cf1.index==cf2.index:
                        return cf1.index
        raise ValueError('No meet exists between elements ' + str(c1) + ' and ' + str(c2))
        
    def meetMatrix(self,func=None,default=0):
        """Form the meet matrix of a poset and a function, in which pairs with no meet have a default value."""
        mat=np.zeros((len(self.cells),len(self.cells)))
        
        for i in range(len(self.cells)):
            for j in range(len(self.cells)):
                try:
                    idx=self.meet(i,j)
                    if func==None:
                        mat[i,j]=1
                    else:
                        mat[i,j]=func(idx)
                except ValueError:
                    mat[i,j]=default
                    
        return mat

class AbstractSimplicialComplex(CellComplex):
    def __init__(self,complex):
        """An abstract simplicial complex defined as a list of lists"""
        pass

class SheafCoface(Coface):
    """A coface relation"""
    def __init__(self,index,orientation,corestriction):
        self.index=index
        self.orientation=orientation
        self.corestriction=corestriction

    def __repr__(self):
        return "(index=" + str(self.index) + ",orientation="+str(self.orientation)+",corestriction="+str(self.corestriction)+")"
        
class SheafCell(Cell):
    """A cell in a cell complex with a sheaf over it"""
    def __init__(self,dimension,cofaces=[],compactClosure=True,stalkDim=1):
        if cofaces:
            try:
                self.stalkDim=cofaces[0].corestriction.shape[1]
            except AttributeError:
                self.stalkDim=0
        else:
            self.stalkDim=stalkDim
        Cell.__init__(self,dimension,compactClosure,cofaces)

    def __repr__(self):
        string="(dimension="+str(self.dimension)+",compactClosure="+str(self.compactClosure)
        if self.cofaces:
            for cf in self.cofaces:
                string+="," + cf.__repr__()
            return string+")"
        else:
            return string+",stalkdim="+str(self.stalkDim)+")"

# Sheaf class
class Sheaf(CellComplex):
    def cofaces(self,c,cells=[],currentcf=[]):
        """Iterate over cofaces (of all dimensions) of a given cell; optional argument specifies which cells are permissible cofaces"""
        
        for cf in self.cells[c].cofaces:
            if cf.index in cells or not cells:
                if currentcf:
                    cfp=SheafCoface(cf.index,
                        currentcf.orientation*cf.orientation,
                        np.dot(currentcf.corestriction,corestriction))
                else:
                    cfp=cf
                yield cfp
                
        for cf in self.cells[c].cofaces:
            if cf.index in cells or not cells:
                if currentcf:
                    cfp=SheafCoface(cf.index,
                        currentcf.orientation*cf.orientation,
                        np.dot(currentcf.corestriction,corestriction))
                else:
                    cfp=cf
                self.cofaces(cf.index,cells,cfp)

    def star(self,cells):
        """Restrict a sheaf to the star over a subset of the base space"""
        
        # Extract a list of all relevant cells in the star
        cells=CellComplex.starCells(self,cells)
        
        return Sheaf([SheafCell(dimension=self.cells[i].dimension,
                                stalkDim=self.cells[i].stalkDim,
                                compactClosure=self.cells[i].compactClosure and (not set(self.faces(i)).difference(set(cells))),
                                cofaces=[SheafCoface(cells.index(cf.index),cf.orientation,cf.corestriction) for cf in self.cells[i].cofaces]) for i in cells]) 

    def kcells(self,k,compactSupport=False):
        """Extract the compact k-cells and associated components of coboundary matrix"""
        k=CellComplex.skeleton(self,k,compactSupport)
        ksizes=[self.cells[i].stalkDim for i in k]
        kidx=list(cumulative_sum(ksizes))
        return k,ksizes,kidx
        
    def localSectional(self,cells=[]):
        """Construct a new sheaf whose global sections are the local sections of the current sheaf over the given cells, and a morphism from this new sheaf to the current one."""
        
        if not cells:
            cells=range(len(self.cells))
            
        # Edges of new sheaf = elements of S with at least one face in the list of cells
        edges=[i for i in cells if self.isFaceOf(i,cells)]

        morphism=[]
        newcells=[]                                        
        for i in edges:
            newcells.append(SheafCell(1,
                compactClosure=self.cells[i].compactClosure and (not set(self.faces(i)).difference(set(cells))),
                stalkDim=self.cells[i].stalkDim))
            morphism.append(SheafMorphismCell([i],[np.eye(self.cells[i].stalkDim)]))

        # Vertices of new sheaf = elements of S with no faces
        vert=list(set(cells).difference(edges))  

        # Corestrictions of new sheaf = compositions of corestrictions
        for i in vert:
            # starting at this vertex, do a depth-first search for the edges in the new sheaf
            cofaces=list(self.cofaces(i,cells))
            newcofaces=[]
            for cf in cofaces:
                newcofaces.append(SheafCoface(edges.index(cf.index),
                    cf.orientation,
                    cf.corestriction))
            
            if cofaces:
                newcells.append(SheafCell(0,compactClosure=True,cofaces=newcofaces))
            else:
                newcells.append(SheafCell(0,compactClosure=True,stalkDim=self.cells[i].stalkDim))
                
            morphism.append(SheafMorphismCell([i],[np.eye(self.cells[i].stalkDim)]))

        return Sheaf(newcells),morphism
        
    def localRestriction(self,cells_1,cells_2):
        """Compute the map induced on local sections by restricting from a larger set to a smaller one"""
        
        # Obtain sheaves and morphisms of local sections for both sets
        sheaf_1,mor_1=self.localSectional(cells_1)
        sheaf_2,mor_2=self.localSectional(cells_2)

        # Compute global sections of each sheaf in terms of cohomology
        H0_1=sheaf_1.cohomology(0)
        if not np.all(H0_1.shape):
            return np.zeros(H0_1.shape)

        # Extend sections of sheaf 1 to vertices of sheaf 2, if needed
        # Observe that value at each vertex of sheaf 2 either
        #  (1) comes from value at a vertex of sheaf 1 or
        #  (2) comes from value at an edge of sheaf 1, 
        #      in which case a single corestriction map obtains it 
        #      from the value at a vertex of sheaf 1
        k_1,ksizes_1,kidx_1=sheaf_1.kcells(0)
        k_2,ksizes_2,kidx_2=sheaf_2.kcells(0)
        rows=sum(ksizes_2)
        sections=np.zeros((rows,H0_1.shape[1]))
        for ss in range(H0_1.shape[1]): # Looping over sections in sheaf 1
            for i in range(len(k_2)): # Looping over vertices in sheaf 2
                # Compute compute preimages of this sheaf 2 vertex
                ms=[k for k in range(len(mor_1)) if
                    set(mor_1[k].destinations).intersection(mor_2[k_2[i]].destinations)]
                if ms:
                    if sheaf_1.cells[ms[0]].dimension==0:
                        ii=ms[0]
                        idx=k_1.index(ii)
                        map,j1,j2,j3=np.linalg.lstsq(mor_2[i].maps[0],mor_1[ii].maps[0])
                        A=np.dot(map,H0_1[kidx_1[idx]:kidx_1[idx+1],ss])
                        sections[kidx_2[i]:kidx_2[i+1],ss]=A
                    else:
                        ii=sheaf_1.faces(ms[0])[0] # parent cells
                        idx=k_1.index(ii)
                        for cf in sheaf_1.cells[ii].cofaces:
                            if cf.index==ms[0]:
                                cr=cf.corestriction
                                break
                        A=np.dot(cr,mor_1[ii].maps[0])
                        
                        map,j1,j2,j3=np.linalg.lstsq(mor_2[i].maps[0],A)
                        sections[kidx_2[i]:kidx_2[i+1],ss]=np.dot(map,H0_1[kidx_1[idx]:kidx_1[idx+1],ss])

        # Rewrite sections over sheaf 2 in terms of 0-cohomology basis
                
        map,j1,j2,j3 = np.linalg.lstsq(sections,sheaf_2.cohomology(0))
        return map.conj().T

    # Input: k = degree of cohomology to compute
    # Output: matrix
    def coboundary(self,k,compactSupport=False):
        """Compute k-th coboundary matrix"""
        # Collect the k-cells and k+1-cells
        ks,ksizes,kidx=self.kcells(k,compactSupport)
        kp1,kp1sizes,kp1idx=self.kcells(k+1,compactSupport)
        
        # Allocate output matrix
        rows=sum(kp1sizes)
        cols=sum(ksizes)
        d=np.zeros((rows,cols),dtype=np.complex)
        if rows and cols:
            # Loop over all k-cells, writing their matrices into the output matrix
            for i in range(len(ks)):
                # Loop over cofaces with compact closure
                for cf in self.cells[ks[i]].cofaces:
                    if self.cells[cf.index].compactClosure or compactSupport:
                        ridx=kp1.index(cf.index)
                        block=np.matrix(cf.orientation*cf.corestriction)
                        d[kp1idx[ridx]:kp1idx[ridx+1],kidx[i]:kidx[i+1]]+=block
            return d
        else:
            return d

    def cohomology(self,k,compactSupport=False,tol=1e-5):
        """Compute basis for k-th cohomology of the sheaf"""
        
        # Obtain coboundary matrices for the sheaf
        dm1=self.coboundary(k-1,compactSupport)
        dm1=np.compress(np.any(abs(dm1)>tol,axis=0),dm1,axis=1)
        d=self.coboundary(k,compactSupport)
        
        # Compute kernel
        if d.size:
            ker=kernel(d,tol);
        else:
            ker=np.eye(d.shape[1])

        # Remove image
        if k > 0 and dm1.any():
            map,j1,j2,j3=np.linalg.lstsq(ker,dm1)
            Hk=np.dot(ker,cokernel(map,tol));
        else:
            Hk=ker

        return Hk
        
    def betti(self,k,compactSupport=False,tol=1e-5):
        """Compute the k-th Betti number of the sheaf"""
        return self.cohomology(k,compactSupport).shape[1]

    def pushForward(self,targetComplex,map):
        """Compute the pushforward sheaf and morphism along a map"""
        
        sheafCells=[]
        mor=[]
        # Loop over cells in the target cell complex
        for cidx in range(len(targetComplex.cells)):
            c=targetComplex.cells[cidx]
            
            # Compute which cells are in the preimage of this cell
            bigPreimage=[d for d,r in map if r==cidx]
            
            # For each cell, compute map on global sections over the star
            # along each attachment
            cfs=[]
            for cf in c.cofaces:
                smallPreimage=[d for d,r in map if r==cf.index]
                corest=self.localRestriction(self.starCells(bigPreimage),
                    self.starCells(smallPreimage))
                cfs.append(SheafCoface(cf.index,
                    cf.orientation,corest))
                    
            mor.append(SheafMorphismCell(bigPreimage,
                [self.localRestriction(self.starCells(bigPreimage),[d]) for d in bigPreimage]))
            if cfs:
                sheafCells.append(SheafCell(c.dimension,cfs,c.compactClosure))
            else:
                ls,m=self.localSectional(self.starCells(bigPreimage))
                sheafCells.append(SheafCell(c.dimension,[],c.compactClosure,stalkDim=ls.betti(0)))

        return Sheaf(sheafCells),mor
        
    def flowCollapse(self):
        """Compute the sheaf morphism to collapse a sheaf to a flow sheaf over the same space"""
         
        # Generate the flow sheaf
        fs=FlowSheaf(self)
         
        mor=[]
        for i in range(len(self.cells)):
            c=self.cells[i]
             
            # If a vertex, collapse by composing edge morphism with corestrictions
            if c.dimension==0:
                map=np.zeros((0,c.stalkDim))
                for j in range(len(c.cofaces)-1):
                    cf=c.cofaces[j]
                    map=np.vstack((map,np.sum(cf.corestriction,axis=0)))

                mor.append(SheafMorphismCell([i],[map]))
            else:
                # If an edge, collapse by summing
                mor.append(SheafMorphismCell([i],[np.ones((1,c.stalkDim))]))
                 
        return fs,mor
        
class AmbiguitySheaf(Sheaf):
    def __init__(self,shf1,mor):
        """Construct an ambiguity sheaf from two sheaves (over the same base) and a morphism between them"""
        
        cellsnew=[]
        for i in range(len(shf1.cells)):
            c=shf1.cells[i]
            
            # New cell has same dimension, compactness,
            # Stalk is the kernel of the component map there
            # Corestrictions come from basis change on each corestriction
            K=kernel(mor[i].map[0])
            stalkDim=K.shape[0]
            cfnew=[]
            for cf in shf1.cells[i].cofaces:
                S=cf.corestriction
                L=kernel(mor[cf.index].map[0])
                R=np.linalg.lstsq(L,np.dot(S,K))
                cfnew.append(SheafCoface(index=cf.index,
                    orientation=cf.orientation,
                    corestriction=R))
                    
            cellsnew.append(SheafCell(dimension=c.dimension,
                compactClosure=c.compactClosure,
                stalkDim=stalkDim,       
                cofaces=cfnew))
                
        Sheaf.__init__(self,cellsnew)

class LocalHomologySheaf(Sheaf):
    def __init__(self,cellcomplex,k):
        shcells=[]

        for i,c in enumerate(cellcomplex.cells):
            shcells.append(SheafCell(c.dimension,
                                     compactClosure=c.compactClosure,
                                     stalkDim=cellcomplex.localHomology(k,i).shape[1],
                                     cofaces=[SheafCoface(index=cf.index,
                                                          orientation=cf.orientation,
                                                          corestriction=cellcomplex.inducedMapLocalHomology(k,i,cf.index))
                                              for cf in c.cofaces]))
        Sheaf.__init__(self,shcells)
        
# Poset sheaves
class ChainSheaf(Poset,Sheaf):
    def __init__(self,poset):
        """Sheaf of chains of a poset or directed graph"""
        
        shcells=[]
        for i,c in enumerate(poset.cells):
            chains=poset.maximalChains(i)
            shcells.append(SheafCell(c.dimension,
                                     compactClosure=c.compactClosure,
                                     stalkDim=len(chains),
                                     cofaces=[SheafCoface(index=cf.index,
                                                          orientation=cf.orientation,
                                                          corestriction=subchainMatrix(chains,
                                                                                       poset.maximalChains(cf.index))) 
                                              for cf in c.cofaces]))

        Sheaf.__init__(self,shcells)

# Flow sheaves
class DirectedGraph(CellComplex):
    def __init__(self,graph,vertex_capacity=-1):
        """Create a cell complex from a directed graph description, which is a list of pairs (src,dest) or triples (src,dest,capacity) of numbers representing vertices.
        The vertex labeled None is an external connection
        Cells are labeled as follows:
         First all of the edges (in the order given), 
         then all vertices (in the order they are given; not by their numerical
         values)"""

        # Construct list of vertices
        verts=[]
        for ed in graph:
            s=ed[0]
            d=ed[1]
            if s != None:
                verts.append(s)
            if d != None:
                verts.append(d)
        verts=list(set(verts))

        # Loop over edges, creating cells for each
        compcells=[]
        for i in range(len(graph)):
            compcells.append(Cell(dimension=1,
                                  compactClosure=(graph[i][0]!=None) and (graph[i][1]!=None)))
            compcells[-1].vertex_label=None
            try: # Add capacity if specified
                compcells[-1].capacity = graph[i][2]
            except:
                pass

        # Loop over vertices, creating cells for each
        for i in verts:
            # Collect cofaces
            cfs=[j for j in range(len(graph)) if graph[j][0]==i or graph[j][1]==i]
            # Compute orientations of each attachment
            orient=[]
            cofaces=[]
            for j in range(len(cfs)):
                if graph[cfs[j]][0]==i:
                    orient.append(-1)
                else:
                    orient.append(1)
                cofaces.append(Coface(cfs[j],orient[j]))
            
            compcells.append(Cell(dimension=0,
                                  compactClosure=True,
                                  cofaces=cofaces))
            compcells[-1].vertex_label=i
            compcells[-1].capacity=vertex_capacity

        CellComplex.__init__(self,compcells)

    def findPath(self,start,end,history=[]):
        """Find a path from specified start cell to end cell
Cell attribute .capacity_left specifes whether the cell can be used"""
        if start == end:
            return history+[end]
            
        # Initialize the capacities used, if unavailable
        for c in self.cells:
            if not hasattr(c,'capacity_left'):
                c.capacity_left = 1

        if self.cells[start].dimension == 0:
            # Compute list of outgoing edges
            for cf in self.cells[start].cofaces:
                if cf.orientation == -1 and cf.index not in history and self.cells[cf.index].capacity_left:
                   ch=self.findPath(cf.index,end,history+[start])
                       
                   if ch:
                       return ch
            return None
        else:
            # Locate vertices which this edge is pointing into
            fs=[i for i in range(len(self.cells)) if self.cells[i].isCoface(start,1)]
            # Is there is a vertex with remaining capcity?
            if fs and fs[0] not in history and self.cells[fs[0]].capacity_left:
                return self.findPath(fs[0],end,history+[start])
            else: # No such vertex
                return None

    def maxFlow(self,start,end):
        """Compute the maximal flow through a graph using Ford-Fulkerson algorithm.  Cell attribute .capacity specifies the number of times the cell can be used.  The default capacities are 1 for edges, infinity for vertices, which results in finding all edge-disjoint paths."""
        # Initialize the capacities on intermediate cells
        for c in self.cells:
            try: # ... to use requested capacities
                c.capacity_left = c.capacity
            except: # if no capacity specified
                if c.dimension == 0: # Vertices get infinite capacity
                    c.capacity_left=-1
                else: # Edges get capacity 1
                    c.capacity_left = 1
                
        # Initialize start/end cell capacities to be infinite
        self.cells[start].capacity_left=-1
        self.cells[end].capacity_left=-1
        
        # Search for paths
        chains=[]
        ch=self.findPath(start,end)
        while ch:
            # Add list of chains
            chains+=[ch]

            # Delete capacities from the cells in this chain
            for i in ch:
                self.cells[i].capacity_left -= 1

            # Find the next chain
            ch=self.findPath(start,end)

        return chains

    def maximalChains(self,start,history=[]):
        """Compute a list of maximal chains beginning at a cell"""

        if self.cells[start].dimension == 0:
            # Compute list of outgoing edges
            cfs=[cf for cf in self.cells[start].cofaces 
                if cf.orientation == -1 and not cf.index in history]
            if cfs: # There are outgoing edges, loop over them
                chains=[self.maximalChains(cf.index,history+[start]) 
                        for cf in cfs]
                lst=[]
                for ch in chains:
                    lst+=ch
                return lst
            else: # No outgoing edges, so this vertex is terminal
                return [history+[start]]
        else:
            # Locate vertices which this edge is pointing into
            fs=[i for i in range(len(self.cells)) if self.cells[i].isCoface(start,1)]
            # Is there is a vertex with remaining capcity
            if fs and not fs[0] in history:
                return self.maximalChains(fs[0],history+[start])
            else: # No such vertex
                return [history+[start]]

    def coveringSpace(self,sheets,partial=False):
        """Create a directed graph that is a covering space of this one with the specified number of sheets"""
        # Construct new edge set
        edges=[c for c in self.cells if c.dimension == 1]
        newcells=edges*sheets
        edgeidx=[idx for idx in range(0,len(self.cells)) 
                 if self.cells[idx].dimension==1]

        # Decompactify edges on request
        if partial:
            for i in range(0,len(edges)):
                newcells[i].compactClosure=False

        # Construct new vertex set
        for c in self.cells:
            if c.dimension == 0:
                for i in range(0,sheets):
                    if i == sheets-1 and partial:
                        break # Skip last vertex copy if requested

                    # Remap cofaces
                    cofaces=[Coface((edgeidx.index(cf.index)+
                                     (i+(1-cf.orientation)/2)*len(edges))
                                    %(len(edges)*sheets),
                                    cf.orientation)
                             for cf in c.cofaces]
                    newcells.append(Cell(dimension=0,
                                         compactClosure=True,
                                         cofaces=cofaces))

        return CellComplex(newcells)

def erdosRenyiDirectedGraph(nvert,prob):
    """Create a random graph with nvert vertices and probability of an edge prob"""
    return DirectedGraph([(a,b) for a in range(nvert)+[None]
                          for b in range(nvert)+[None]
                          if random.random() < prob and (a!=None or b!=None)])
        
class FlowSheaf(Sheaf,DirectedGraph):
    def __init__(self,graph):
        """Create a flow sheaf from a directed graph"""
        
        sheafcells=[]
        for c in graph.cells:
            cofaces=[]
            j=0
            for cf in c.cofaces:
                # Compute corestrictions
                if j in range(len(c.cofaces)-1):
                    corest=np.matrix([m==j for m in range(len(c.cofaces)-1)],dtype=int)
                else:
                    corest=np.matrix([cf.orientation for cf in c.cofaces][0:-1])
                
                cofaces.append(SheafCoface(cf.index,cf.orientation,corest))
                j+=1
            
            sheafcells.append(SheafCell(dimension=c.dimension,
                                        compactClosure=c.compactClosure,
                                        cofaces=cofaces))

        Sheaf.__init__(self,sheafcells)

class TransLineSheaf(Sheaf,DirectedGraph):
    def __init__(self,graph,wavenumber):
        """Create a transmission line sheaf from a directed graph, in which edges have been given a .length attribute"""

        # Default edge lengths
        for c in graph.cells:
            if c.dimension == 1:
                if c.compactClosure==False:
                    c.length=1
                try:
                    if c.length<0:
                        c.length=1
                except:
                    c.length=1
        
        sheafcells=[]

        for c in graph.cells:
            cofaces=[]

            if c.dimension == 0: # Edges have interesting corestrictions
                n=len(c.cofaces)
                phaselist=[2/n for i in range(n)]
                for m in range(n):
                    if c.cofaces[m].orientation == -1:
                        corest=np.matrix([[i==m for i in range(n)],
                                          phaselist],
                                         dtype=complex)
                        corest[1,m]-=1
                        corest[1,:]*=np.exp(-1j*wavenumber*graph.cells[c.cofaces[m].index].length)
                    else:
                        corest=np.matrix([phaselist,
                                          [i==m for i in range(n)]],
                                         dtype=complex)
                        corest[0,m]-=1
                        corest[0,:]*=np.exp(1j*wavenumber*graph.cells[c.cofaces[m].index].length)
                    cofaces.append(SheafCoface(c.cofaces[m].index,
                                               c.cofaces[m].orientation,
                                               corest))
            else: # All other faces have trivial corestrictions
                n=2
                cofaces=[SheafCoface(cf.index,cf.orientation,[]) for cf in c.cofaces]
            sheafcells.append(SheafCell(dimension=c.dimension,
                                        compactClosure=c.compactClosure,
                                        cofaces=cofaces,
                                        stalkDim=n))

        Sheaf.__init__(self,sheafcells)

class ConstantSheaf(Sheaf):
    def __init__(self,cells):
        """Construct a constant sheaf over a CellComplex"""
        sheafcells=[SheafCell(dimension=c.dimension,
                              compactClosure=c.compactClosure,
                              cofaces=[SheafCoface(cf.index, 
                                                   cf.orientation,
                                                   np.matrix(1))
                                       for cf in c.cofaces],
                              stalkDim=1)
                    for c in cells]

        Sheaf.__init__(self,sheafcells)
        pass

class SheafMorphismCell:
    def __init__(self,destinations=[],maps=[]):
        """Specify destinations and maps for this cell's stalk under a morphism"""
        self.destinations=destinations
        self.maps=maps

# A local section
class SectionCell:
    def __init__(self,support,value):
        """Specify support cell indices and values in each cell stalk for a local section"""
        self.support=support
        self.value=value

class Section:
    def __init__(self,sectionCells):
        self.sectionCells=sectionCells

    def support(self):
        """List the cells in the support of this section"""
        return [sc.support for sc in self.sectionCells]

    def extend(self,sheaf,cell,value=None,tol=1e-5):
        """Extend the section to another cell; returns True if successful"""
        # If the desired cell is already in the support, do nothing
        if cell in self.support():
            return True

        # Is the desired cell a coface of a cell in the support?
        for s in self.sectionCells:
            for cf in sheaf.cells[s.support].cofaces:
                if cf.index == cell:
                    # If so, extend via corestriction
                    val=np.dot(cf.corestriction,s.value)

                    # Check for consistency
                    if value != None and np.any(np.abs(val - value)>tol):
                        return False
                    value = val
                            
        # Are there are any cofaces for the desired cell in the support?
        if value == None: # Attempt to assign a new value...
            # Stack the corestrictions and values associated to existing support
            lst=[(cf.corestriction,s.value) 
                 for cf in sheaf.cells[cell].cofaces 
                 for s in self.sectionCells
                 if cf.index == s.support]
            if lst:
                crs=np.vstack([e[0] for e in lst])
                vals=np.vstack([e[1] for e in lst])

                # If the problem of solving for the value at this cell is 
                # underdetermined, refrain from solving it
                if matrixrank(crs,tol) < crs.shape[1]:
                    return True
            
                # Attempt to solve for the value at desired cell
                val,res,j2,j3=np.linalg.lstsq(crs,vals)
                if np.any(np.abs(res)>tol):
                    return False
                value = val
        else: # ...or check consistency with an old one
            for cf in sheaf.cells[cell].cofaces:
                for s in self.sectionCells:
                    if s.support == cf.index:
                        if np.any(np.abs(np.dot(cf.corestriction,value)-s.value)>tol):
                            return False
        
        # A value was successfully assigned (if no value was assigned, 
        # do nothing, but it's still possible to extend)
        if value != None:
            self.sectionCells.append(SectionCell(cell,value))

        return True

class PersistenceSheaf(Sheaf):
    # Compute k-th persistence sheaf
    # Input: list of sheaves
    #        list of triples: source sheaf index, destination sheaf index, sheaf morphism data
    def __init__(self,sheaves,morphisms,k):
        """Compute the k-th degree persistence sheaf over a graph"""
    
        persheaf=[]

        # Loop over sheaves
        for i in range(len(sheaves)):
            # Loop over morphisms initiated from this sheaf
            cofaces=[]
            for (s,d,mor) in morphisms:
                if s==i:
                    cofaces.append(SheafCoface(d,
                                               1,
                                               inducedMap(sheaves[i],sheaves[d],mor,k)))
            if cofaces:
                persheaf.append(SheafCell(dimension=0,
                                          compactClosure=True,
                                          cofaces=cofaces))
            else: # If cell does not have cofaces, compute stalk from Betti number
                persheaf.append(SheafCell(dimension=1,
                                          compactClosure=len([d for (s,d,mor) in morphisms 
                                                              if d==i])>1,
                                          stalkDim=sheaves[i].betti(k)))
        # Initialize the sheaf
        Sheaf.__init__(self,persheaf)
            
## Functions

def cumulative_sum(values, start=0):
    yield start
    for v in values:
        start += v
        yield start

def matrixrank(A,tol=1e-5):
    u, s, vh = np.linalg.svd(A)
    return sum(s > tol)

def kernel(A, tol=1e-5):
    u, s, vh = np.linalg.svd(A)
    sing=np.zeros(vh.shape[0],dtype=np.complex)
    sing[:s.size]=s
    null_mask = (sing <= tol)
    null_space = np.compress(null_mask, vh, axis=0)
    return null_space.conj().T

def cokernel(A, tol=1e-5):
    u, s, vh = np.linalg.svd(A)
    sing=np.zeros(u.shape[1],dtype=np.complex)
    sing[:s.size]=s
    null_mask = (sing <= tol)
    return np.compress(null_mask, u, axis=1)

def isSubchain(bigger,smaller):
    """Determine if smaller is a subchain of bigger"""
    try:
        idx=bigger.index(smaller[0])
    except:
        return False
    for s in smaller:
        try:
            if bigger[idx] != s:
                return False
            else:
                idx += 1
        except:
            return False
    return True

def subchainMatrix(chainsCols,chainsRows):
    """Construct a binary matrix specifying which chains in the second set are subchains of the first"""
    mat=np.zeros((len(chainsRows),len(chainsCols)))
    for i in range(len(chainsCols)):
        for j in range(len(chainsRows)):
            if isSubchain(chainsCols[i],chainsRows[j]) or isSubchain(chainsRows[j],chainsCols[i]):
                mat[j,i]=1
    return mat

# Input: domain sheaf
#        range sheaf
#        list of sheaf morphism data, one for each cell
#        k 
# Output: matrix
def inducedMap(sheaf1,sheaf2,morphism,k,compactSupport=False,tol=1e-5):
    """Compute k-th induced map on cohomology for a sheaf morphism"""
    
    # Compute cohomology basis for each sheaf
    Hk_1=sheaf1.cohomology(k,compactSupport)
    Hk_2=sheaf2.cohomology(k,compactSupport)

    if (not Hk_1.size) or (not Hk_2.size):
        return []

    # Extract the k-skeleta of each sheaf
    k_1,ksizes_1,kidx_1=sheaf1.kcells(k,compactSupport)
    k_2,ksizes_2,kidx_2=sheaf2.kcells(k,compactSupport)

    # Construct chain map
    rows=sum(ksizes_2)
    cols=sum(ksizes_1)
    m=np.zeros((rows,cols),dtype=np.complex)

    for i in range(len(k_1)):
        for j,map in zip(morphism[k_1[i]].destinations,morphism[k_1[i]].maps):
            if sheaf2.cells[j].dimension==k:
                ridx=[q for q in range(len(k_2)) if k_2[q]==j]
                if ridx:
                    ridx=ridx[0]
                    m[kidx_2[ridx]:kidx_2[ridx+1],kidx_1[i]:kidx_1[i+1]]+=map

    # Map basis for domain sheaf's cohomology through the chain map
    im=np.dot(m,Hk_1)
    
    # Expand output in new basis
    map,j1,j2,j3 = np.linalg.lstsq(im,Hk_2)

    return map.conj().T
