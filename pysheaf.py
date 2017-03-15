# Persistence-capable sheaf manipulation library
#
# Copyright (c) 2013-2017, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 

import numpy as np
import random
import matplotlib.pyplot as plt
import networkx as nx
import scipy.optimize

## Data structures
class Coface: 
    """A coface relation"""
    def __init__(self,index=None,orientation=None):
        self.index=index
        self.orientation=orientation

    def __repr__(self):
        return "(index=" + str(self.index) + ",orientation="+str(self.orientation)+")"

class Cell:
    """A cell in a cell complex.  The name attribute can be used as a key to index into the registry of cells or cofaces"""
    def __init__(self,dimension,compactClosure=True,cofaces=None, name=None):
        self.dimension=dimension
        self.compactClosure=compactClosure
        if cofaces is not None:
            self.cofaces = cofaces
        else:
            self.cofaces = []
        self.name = name

    def __repr__(self):
        string= self.name + " (dimension="+str(self.dimension)+",compactClosure="+str(self.compactClosure) 
        if self.cofaces:
            for cf in self.cofaces:
                string+="," + cf.__repr__()
        return string+")"

    def cofaceList(self):
        """Return the list of indicies of the cells which are cofaces of this cell"""
        return [cf.index for cf in self.cofaces]    
    
    def isCoface(self,index,orientation=None):
        """Check if a given cell index is a coface of this cell. Optionally check that the orientation is as given."""
        if orientation==None:
            return index in [cf.index for cf in self.cofaces]
        else:
            return (index,orientation) in [(cf.index,cf.orientation) for cf in self.cofaces]
        
class CellComplex:
    def __init__(self,cells=None):
        """Construct a cell complex from its cells"""
        if cells is None:
            self.cells=[]
        else:
            self.cells=cells

        # Register cell names with indices into the self.cells list        
        self.cell_dict={}
        self.coface_dict={}
        for i,c in enumerate(self.cells):
            if c.name is None: # Build names if not present
                c.name = str(i)
            self.cell_dict[c.name]=i
            if c.cofaces is not None: # Register cofaces if present
                for j,cf in enumerate(c.cofaces):
                    self.coface_dict[(i,cf.index)]=j

    def add_cell(self,cell):
        """Add a cell to the cell complex"""

        if cell.name is None: # Construct a name if needed
            cell.name = str(len(self.cells))

        # Register the cell
        self.cell_dict[cell.name]=len(self.cells)

        # Register its cofaces, if any
        if cell.cofaces is not None:
            for j,cf in enumerate(cell.cofaces):
                self.coface_dict[(len(self.cells),cf.index)]=j

        # Store the cell
        self.cells.append(cell)

    def add_cells_from(self,cells):
        for c in cells:
            self.add_cell(c)

    def add_coface(self,cellpair,orientation):
        """Add a coface to the cell complex, referenced by a pair of cell names.  The cellpair argument is assumed to be a pair: (face,coface).  If the cells aren't present, this will raise KeyError."""
        # Look up which cells are involved...
        source=self.cell_dict[cellpair[0]]
        target=self.cell_dict[cellpair[1]]

        # Drop in the coface
        self.coface_dict[(source,coface.index)]=len(self.cells[source].cofaces)
        self.cells[source].cofaces.append(Coface(index=target,orientation=orientation))

    def add_cofaces_from(self,cellpairs,orientations):
        for cellpair,orientation in zip(cellpairs,orientations):
            self.add_coface(cellpair,orientation)
        
    def isFaceOf(self,c,cells=None):
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

    def closure(self,cells):
        """Compute the closure of a collection of cells"""
        return [i for i in range(len(self.cells))
                if i in cells or [cf for cf in self.cofaces(i) if cf.index in cells]]

    def cofaces(self,c,cells=[]):
        """Iterate over cofaces (of all dimensions) of a given cell; optional argument specifies which cells are permissible cofaces
        Warning: duplicates are possible!"""
        for cf in self.cells[c].cofaces:
            if cf.index in cells or not cells:
                for cff in self.cofaces(cf.index,cells):
                    yield cff
                yield cf
                
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

    def homology(self,k,subcomplex=None,compactSupport=False,tol=1e-5):
        """Compute (relative) homology of the cell complex"""
                
        # Obtain boundary matrices for the complex
        d=self.boundary(k,subcomplex,compactSupport)
        dp1=self.boundary(k+1,subcomplex,compactSupport)
        dp1=np.compress(np.any(abs(dp1)>tol,axis=0),dp1,axis=1)
        
        # Compute kernel
        if d.size:
            ker=kernel(d,tol);
        else:
            ker=np.eye(d.shape[1])

        # Remove image
        if dp1.any():
            map,j1,j2,j3=np.linalg.lstsq(ker,dp1)
            Hk=np.dot(ker,cokernel(map,tol));
        else:
            Hk=ker

        return Hk

    def localPairComplex(self,cells):
        """Construct a new cell complex that consists of a cell and its boundary.  The return is the cell complex paired with a list of boundary cells"""
        raise NotImplementedError('localPairComplex is not working!  Please do not use it yet')
        # TBD fix this; it refers to the wrong cells
        # Construct the neighborhood of the cell

        star_closure = self.closure(self.starCells(cells))  ### this is the containing complex
        star=self.starCells(cells)    ###### this is the star of the set we are interested in
        bndry = list(set(star_closure) - set(star))
        starcells = [self.cells[idx] for idx in star]
        bndrycells = [self.cells[idx] for idx in bndry]
        cplx = CellComplex(starcells + bndrycells)
        bndryind = [cplx.cells.index(bdcell) for bdcell in bndrycells]
        return (cplx,bndryind)
        
    def localHomology(self,k,cells):
        """Compute local homology localized at the star over a list of cells"""
        cplx,bndry=self.localPairComplex(cells)

        # Compute the relative homology of the proxy complex
        return cplx.homology(k,subcomplex=bndry)   #######ERROR bndry are indicies of self not thet4r cmplx
        
    def boundary(self,k,subcomplex=None,compactSupport=False):
        """Compute the boundary map for the complex"""
        # Collect the k-cells and k-1-cells
        if subcomplex:
            ks=[spx for spx in self.skeleton(k,compactSupport) if not spx in subcomplex]
            km1=[spx for spx in self.skeleton(k-1,compactSupport) if not spx in subcomplex]
        else:
            ks=self.skeleton(k,compactSupport)
            km1=self.skeleton(k-1,compactSupport)

        # Allocate output matrix
#        print "ks=", ks
#        print "km1=", km1
        rows=len(km1)
        cols=len(ks)
#         d=np.zeros((rows,cols),dtype=np.complex)
        d=np.zeros((rows,cols))
        if rows and cols:
            # Loop over all k-1-cells, writing them into the output matrix
            for i in range(len(km1)):
#                print i, self.cells[km1[i]]
                # Loop over faces with compact closure
#                print "cofaces=", [cf.index for cf in self.cells[km1[i]].cofaces]
                for cf in self.cells[km1[i]].cofaces:
                    
#                    print cf.index, self.cells[cf.index], cf.orientation,
#                    if self.cells[cf.index].compactClosure or compactSupport:
                    if self.cells[cf.index].compactClosure and cf.orientation != 0:
                        d[i,ks.index(cf.index)]=cf.orientation
#                        print "ok"
            return d
        else:
            return d

    def inducedMapLocalHomology(self,k,cell1,cell2):
        """Compute the induced map on local homology between two cells.  It is assumed that cell2 is a coface (perhaps not codimension 1) of cell1"""
        # Compute local homology basis centered on each cell
        pass

    def attachDiagram(self):
        """Draw the attachment diagram using NetworkX"""
        G=nx.DiGraph()

        G.add_nodes_from(range(len(self.cells)))

        G.add_edges_from([(i,cf.index) 
                          for i in range(len(self.cells)) 
                          for cf in self.cells[i].cofaces])
        nx.draw(G)
        plt.show()
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
    def __init__(self,toplexes,maxdim=None):
        """An abstract simplicial complex defined as a list of lists; it's only necessary to pass a generating set of toplexes
        Beware: this should not be constructed for complexes involving high-dimensional simplices!
        Simplices are sorted from greatest dimension to least"""
        
        if maxdim is None:
            maxdim=max([len(tplx)-1 for tplx in toplexes])
            
        self.cells=[]
        upsimplices=[] # Simplices of greater dimension than currently being added
        upindices=[]
        for k in range(maxdim,-1,-1): 
            simplices=ksimplices(toplexes,k) # Simplices to be added
            startindex=len(self.cells)
            for s in simplices:
                cell=Cell(dimension=k,compactClosure=True,
                          cofaces=[Coface(index=upindices[i],
                                          orientation=simplexOrientation(s,upsimplices[i])) for i in range(len(upsimplices)) if set(s).issubset(upsimplices[i])])
                self.cells.append(cell)

            upsimplices=simplices
            upindices=[startindex+i for i in range(len(simplices))]


class SetMorphism():
    """A morphism in a subcategory of Set, described by a function object"""
    def __init__(self,fcn):
        self.fcn=fcn

    def __mul__(self,other): # Composition of morphisms
        return SetMorphism(lambda x : self.fcn(other.fcn(x)))

    def __call__(self,arg): # Calling the morphism on an element of the set
        return self.fcn(arg)

class LinearMorphism(SetMorphism):
    """A morphism in a category that has a matrix representation"""
    def __init__(self,matrix):
        self.matrix=matrix
        SetMorphism.__init__(self,lambda x: np.dot(matrix,x))

    def __mul__(self,other): # Composition of morphisms
        try: # Try to multiply matrices.  This might fail if the other morphism isn't a LinearMorphism
            return LinearMorphism(np.dot(self.matrix, other.matrix))
        except AttributeError:
            return SetMorphism.__mul__(self,other)

class SheafCoface(Coface):
    """A coface relation, in which the restriction is assumed to be a SetMorphism object
    If the restriction is instead a matrix, then it gets promoted to a LinearMorphism object automatically"""
    def __init__(self,index=None,orientation=None,restriction=None):
        self.index=index
        self.orientation=orientation
        if isinstance(restriction,np.ndarray):
            self.restriction=LinearMorphism(restriction)
        else:
            self.restriction=restriction

    def __repr__(self):
        if self.isLinear():
            return "(index=" + str(self.index) + ",orientation="+str(self.orientation)+",restriction="+str(self.restriction.matrix)+")"
        else:
            return "(index=" + str(self.index) + ",orientation="+str(self.orientation)+",restriction="+str(self.restriction.fcn)+")"

    def isLinear(self):
        """Does this coface relation have a LinearMorphism for a restriction?"""
        try:
            getattr(self.restriction,'matrix')
            return True
        except AttributeError:
            return False
        
class SheafCell(Cell):
    """A cell in a cell complex with a sheaf over it
    cofaces = list of SheafCoface instances, one for each coface of this cell
    stalkDim = dimension of the stalk over this cell (defaults to figuring it out from the cofaces if the restriction is a LinearMorphism) or None if stalk is not a real vector space"""
    def __init__(self,dimension,cofaces=None,compactClosure=True,stalkDim=None,metric=None,name=None):
        if stalkDim is None and cofaces:
            if cofaces[0].isLinear():
                try:  # Try to discern the stalk dimension from the matrix representation. This will fail if the matrix isn't given
                    self.stalkDim=cofaces[0].restriction.matrix.shape[1]
                except AttributeError:
                    self.stalkDim=0
            else:  # OK, not a LinearMorphism, so the stalk isn't a vector space
                self.stalkDim=stalkDim
        else:
            self.stalkDim=stalkDim

        if metric is not None:
            self.metric=metric
        else:
            self.metric=lambda x,y: np.linalg.norm(x-y)
            
        Cell.__init__(self,dimension,compactClosure,cofaces,name)

    def isLinear(self):
        """Is this cell representative of a Sheaf of vector spaces?  (All restrictions are linear maps)"""
        if self.stalkDim is None:
            return False
        for cf in self.cofaces:
            if not cf.isLinear():
                return False
        return True

    def isNumeric(self):
        """Is this cell representative of a sheaf of sets in which stalks are all real vector spaces? (restrictions may not be linear maps, though)"""
        if self.stalkDim is None:
            return False
        else:
            return True

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

    def add_coface(self,cellpair,orientation,restriction):
        """Add a coface to the sheaf, referenced by a pair of cell cellpair.  The cellpair argument is assumed to be a pair: (face,coface).  If the cells aren't present, this will raise KeyError."""
        # Look up which cells are involved...
        source=self.cell_dict[cellpair[0]]
        target=self.cell_dict[cellpair[1]]

        # Drop in the coface
        self.coface_dict[(source,coface.index)]=len(self.cells[source].cofaces)
        self.cells[source].cofaces.append(SheafCoface(index=target,orientation=orientation,restriction=restriction))

    def add_cofaces_from(self,cellpairs,orientations,restrictions):
        for cellpair,orientation,restriction in zip(cellpairs,orientations,restriction):
            self.add_coface(cellpair,orientation,restriction)

    def isLinear(self):
        """Is this a Sheaf of vector spaces?  (All restrictions are linear maps)"""
        for c in self.cells:
            if not c.isLinear():
                return False
        return True

    def isNumeric(self):
        """Is this a Sheaf of sets in which stalks are all real vector spaces? (restrictions may not be linear maps, though)"""
        for c in self.cells:
            if not c.isNumeric():
                return False
        return True
        
    def cofaces(self,c,cells=[],currentcf=[]):
        """Iterate over cofaces (of all dimensions) of a given cell c; optional argument specifies which cells are permissible cofaces"""
        if c >= len(self.cells):
            yield []

        for cf in self.cells[c].cofaces:
            if cf.index in cells or not cells:
                if currentcf: # If we've already started the iteration, there's a previous restriction to compose with
                    cfp=SheafCoface(index=cf.index,
                                    orientation=cf.orientation*currentcf.orientation,
                                    restriction=cf.restriction*currentcf.restriction)
                else: # If we're just starting this iteration, there is no restriction before this one
                    cfp=cf
                for cff in self.cofaces(cfp.index,cells,cfp): # Iterate over all higher dimensional cells
                    yield cff
                yield cfp

    def star(self,cells):
        """Restrict a sheaf to the star over a subset of the base space"""
    
        # Extract a list of all relevant cells in the star
        cells=CellComplex.starCells(self,cells)
        
        return Sheaf([SheafCell(dimension=self.cells[i].dimension,
                                stalkDim=self.cells[i].stalkDim,
                                compactClosure=self.cells[i].compactClosure and (not set(self.faces(i)).difference(set(cells))),
                                cofaces=[SheafCoface(index=cells.index(cf.index),
                                                     orientation=cf.orientation,
                                                     restriction=cf.restriction) for cf in self.cells[i].cofaces]) for i in cells]) 

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
            morphism.append(SheafMorphismCell([i],[LinearMorphism(np.eye(self.cells[i].stalkDim))]))

        # Vertices of new sheaf = elements of S with no faces
        vert=list(set(cells).difference(edges))  

        # Restrictions of new sheaf = compositions of restrictions
        for i in vert:
            # starting at this vertex, do a depth-first search for the edges in the new sheaf
            cofaces=list(self.cofaces(i,cells))
            newcofaces=[]
            for cf in cofaces:
                newcofaces.append(SheafCoface(index=edges.index(cf.index),
                    orientation=cf.orientation,
                    restriction=cf.restriction))
            
            if cofaces:
                newcells.append(SheafCell(0,compactClosure=True,cofaces=newcofaces))
            else:
                newcells.append(SheafCell(0,compactClosure=True,stalkDim=self.cells[i].stalkDim))
                
            morphism.append(SheafMorphismCell([i],[LinearMorphism(np.eye(self.cells[i].stalkDim))]))

        return Sheaf(newcells),SheafMorphism(morphism)
        
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
        #      in which case a single restriction map obtains it 
        #      from the value at a vertex of sheaf 1
        k_1,ksizes_1,kidx_1=sheaf_1.kcells(0)
        k_2,ksizes_2,kidx_2=sheaf_2.kcells(0)
        rows=sum(ksizes_2)
        sections=np.zeros((rows,H0_1.shape[1]))
        for ss in range(H0_1.shape[1]): # Looping over sections in sheaf 1
            for i in range(len(k_2)): # Looping over vertices in sheaf 2
                # Compute compute preimages of this sheaf 2 vertex
                ms=[k for k in range(len(mor_1.morphismCells)) if
                    set(mor_1.morphismCells[k].destinations).intersection(mor_2.morphismCells[k_2[i]].destinations)]
                if ms:
                    if sheaf_1.cells[ms[0]].dimension==0:
                        ii=ms[0]
                        idx=k_1.index(ii)
                        map,j1,j2,j3=np.linalg.lstsq(mor_2.morphismCells[i].maps[0].matrix,mor_1.morphismCells[ii].maps[0].matrix)
                        A=np.dot(map,H0_1[kidx_1[idx]:kidx_1[idx+1],ss])
                        sections[kidx_2[i]:kidx_2[i+1],ss]=A
                    else:
                        ii=sheaf_1.faces(ms[0])[0] # parent cells
                        idx=k_1.index(ii)
                        for cf in sheaf_1.cells[ii].cofaces:
                            if cf.index==ms[0]:
                                cr=cf.restriction
                                break
                        A=cr(mor_1.morphismCells[ii].maps[0].matrix)
                        
                        map,j1,j2,j3=np.linalg.lstsq(mor_2.morphismCells[i].maps[0].matrix,A)
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
                        block=np.matrix(cf.orientation*cf.restriction.matrix)
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

    def maximalExtend(self,assignment,multiassign=False,tol=1e-5):
        """Take a partial assignment and extend it to a maximal assignment that's non-conflicting (if multiassign=False) or one in which multiple values can be given to a given cell (if multiassign=True)"""
        for i in range(len(assignment.sectionCells)):
            for cf in self.cofaces(assignment.sectionCells[i].support):
                if not assignment.extend(self,cf.index,tol) and multiassign:
                    assignment.sectionCells.append(SectionCell(cf.index,cf.restriction(assignment.sectionCells[i].value)))
        return assignment

    def consistencyRadius(self,assignment,tol=1e-5):
        """Compute the consistency radius of an approximate section"""
        assignment=self.maximalExtend(assignment,multiassign=True,tol=tol)
        radius=0
        for c1 in assignment.sectionCells:
            for c2 in assignment.sectionCells:
                if c1.support == c2.support:
                    rad = self.cells[c1.support].metric(c1.value,c2.value)
                    if rad > radius:
                        radius = rad
        return radius

    def assignmentMetric(self,assignment1,assignment2):
        """Compute the distance between two assignments"""
        radius=0
        for c1 in assignment1.sectionCells:
            for c2 in assignment2.sectionCells:
                if c1.support == c2.support:
                    rad = self.cells[c1.support].metric(c1.value,c2.value)
                    if rad > radius:
                        radius = rad
        return radius

    def fuseAssignment(self,assignment,tol=1e-5):
        """Compute the nearest global section to a given assignment"""
        if self.isNumeric():
            # The situation where the stalks are all numeric
            res=scipy.optimize.minimize( fun = lambda sec: self.assignmentMetric(assignment,self.deserializeAssignment(sec)),
                                         x0 = self.serializeAssignment(assignment),
                                         constraints = ({'type' : 'eq',
                                                         'fun' : lambda asg: self.consistencyRadius(self.deserializeAssignment(asg))}),
                                         tol = tol )
            globalsection = self.deserializeAssignment(res.x)
        else:
            # The fallback situation, where we need to iterate over global sections manually...
            raise NotImplementedError
        
        return globalsection

    def deserializeAssignment(self,vect):
        """Transform a vector of values for a numeric-valued sheaf into an assignment as a Section instance
        (Note: this is really a helper method and should generally not be used by external callers)"""
        if not self.isNumeric():
            raise TypeError('Cannot deserialize an assignment vector for a non-numeric sheaf')

        scs=[]
        idx=0
        for i in range(len(self.cells)):
            if self.cells[i].stalkDim > 0:
                scs.append(SectionCell(support=i,value=vect[idx:idx+self.cells[i].stalkDim]))
                idx+=self.cells[i].stalkDim
                       
        return Section(scs)

    def serializeAssignment(self,assignment):
        """Transform a partial assignment in a Section instance into a vector of values"""
        if not self.isNumeric():
            raise TypeError('Cannot serialize an assignment vector for a non-numeric sheaf')

        x0 = np.zeros((sum([c.stalkDim for c in self.cells])))

        # Figure out cell boundaries in vector
        idx=0
        idxarray=[]
        for i in range(len(self.cells)):
            if self.cells[i].stalkDim > 0:
                idxarray.append(idx)
                idx+=self.cells[i].stalkDim

        # Pack data into vector.  If there are multiple values assigned, the one appearing last is used
        for cell in assignment.sectionCells:
            if self.cells[cell.support].stalkDim > 0:
                x0[idxarray[cell.support]:idxarray[cell.support]+self.cells[cell.support].stalkDim]=cell.value

        return x0
        
    def partitionAssignment(self,assignment,tol=1e-5):
        """Take an assignment to some cells of a sheaf and return a collection of disjoint maximal sets of cells on which this assignment is a local section"""
        # Extend assignment to all cofaces
        assignment=self.maximalExtend(self,assignment,multiassign=False)
       
        # Filter the cofaces into a cell complex in which the only attachments that are included as those whose data in the sheaf are consistent
        cells=[]
        for i,c in enumerate(self.cells):
            found=False
            for s in assignment.sectionCells:
                if s.support == i:
                    val=s.value
                    found=True
                    break
            if found:
                cofaces=[]
                for cf in c.cofaces:
                    vv=cf.restriction(val)
                    for s in assignment.sectionCells:
                        if s.support == cf.index and np.all(np.abs(vv-s.value)) < tol:
                            cofaces.append(cf)
                            break
                             
                cells.append(Cell(c.dimension,c.compactClosure,cofaces))
            else:
                cells.append(Cell(c.dimension,c.compactClosure,cofaces=[]))

        cplx=CellComplex(cells)
        
        # Compute the components of the complex
        return cplx.components()
    
    def smoothness(self,assignment):
        """ returns the open world smoothness Entropy/log2(number of sections) """
        partition = self.partitionAssignment(assignment)
        arr = [len(x) for x in partition]
        n = sum(arr)*1.0
        m = len(arr)
        if (n==0 or m==1):
            return 1
        else:
            f = [x/n for x in arr]
            E = sum([p*np.log2(p) for p in f])
            return -E/np.log2(m)

    def dispersion(self,assignment):
        """ returns a measure on dispersion on the number of sections to number of 0-cells """
        partition = self.partitionAssignment(assignment)
        arr = [len(x) for x in partition]
        n = sum(arr)*1.0
        m = len(arr)
        if (n==0):
            return 0
        return np.log2(m)/np.log2(n)

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
                rest=self.localRestriction(self.starCells(bigPreimage),
                    self.starCells(smallPreimage))
                cfs.append(SheafCoface(index=cf.index,
                    orientation=cf.orientation,restriction=rest))
                    
            mor.append(SheafMorphismCell(bigPreimage,
                [LinearMorphism(self.localRestriction(self.starCells(bigPreimage),[d])) for d in bigPreimage]))
            if cfs:
                sheafCells.append(SheafCell(c.dimension,cfs,c.compactClosure))
            else:
                ls,m=self.localSectional(self.starCells(bigPreimage))
                sheafCells.append(SheafCell(c.dimension,[],c.compactClosure,stalkDim=ls.betti(0)))

        return Sheaf(sheafCells),SheafMorphism(mor)
        
    def flowCollapse(self):
        """Compute the sheaf morphism to collapse a sheaf to a flow sheaf over the same space"""
         
        # Generate the flow sheaf
        fs=FlowSheaf(self)
         
        mor=[]
        for i in range(len(self.cells)):
            c=self.cells[i]
             
            # If a vertex, collapse by composing edge morphism with restrictions
            if c.dimension==0:
                map=np.zeros((0,c.stalkDim))
                for j in range(len(c.cofaces)-1):
                    cf=c.cofaces[j]
                    map=np.vstack((map,np.sum(cf.restriction.matrix,axis=0)))

                mor.append(SheafMorphismCell([i],[LinearMorphism(map)]))
            else:
                # If an edge, collapse by summing
                mor.append(SheafMorphismCell([i],[LinearMorphism(np.ones((1,c.stalkDim)))]))
                 
        return fs,SheafMorphism(mor)
        
class AmbiguitySheaf(Sheaf):
    def __init__(self,shf1,mor):
        """Construct an ambiguity sheaf from two sheaves (over the same base) and a morphism between them"""
        
        cellsnew=[]
        for i in range(len(shf1.cells)):
            c=shf1.cells[i]
            
            # New cell has same dimension, compactness,
            # Stalk is the kernel of the component map there
            # Restrictions come from basis change on each restriction
            K=kernel(mor.morphismCells[i].maps[0].matrix)
            stalkDim=K.shape[0]
            cfnew=[]
            for cf in shf1.cells[i].cofaces:
                S=cf.restriction.matrix
                L=kernel(mor.morphismCells[cf.index].maps[0].matrix)
                R=np.linalg.lstsq(L,np.dot(S,K))
                cfnew.append(SheafCoface(index=cf.index,
                    orientation=cf.orientation,
                    restriction=R))
                    
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
                                     stalkDim=cellcomplex.localHomology(k,[i]).shape[1],
                                     cofaces=[SheafCoface(index=cf.index,
                                                          orientation=cf.orientation,
                                                          restriction=cellcomplex.inducedMapLocalHomology(k,i,cf.index))
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
                                                          restriction=subchainMatrix(chains,
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
            if s is not None:
                verts.append(s)
            if d is not None:
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
                # Compute restrictions
                if j in range(len(c.cofaces)-1):
                    rest=np.matrix([m==j for m in range(len(c.cofaces)-1)],dtype=int)
                else:
                    rest=np.matrix([cf.orientation for cf in c.cofaces][0:-1])
                
                cofaces.append(SheafCoface(index=cf.index,orientation=cf.orientation,restriction=rest))
                j+=1

            if cofaces:
                sheafcells.append(SheafCell(dimension=c.dimension,
                                            compactClosure=c.compactClosure,
                                            cofaces=cofaces))
            else:
                sheafcells.append(SheafCell(dimension=c.dimension,
                                            compactClosure=c.compactClosure,
                                            cofaces=[],
                                            stalkDim=1))

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

            if c.dimension == 0: # Edges have interesting restrictions
                n=len(c.cofaces)
                phaselist=[2/n for i in range(n)]
                for m in range(n):
                    if c.cofaces[m].orientation == -1:
                        rest=np.matrix([[i==m for i in range(n)],
                                          phaselist],
                                         dtype=complex)
                        rest[1,m]-=1
                        rest[1,:]*=np.exp(-1j*wavenumber*graph.cells[c.cofaces[m].index].length)
                    else:
                        rest=np.matrix([phaselist,
                                          [i==m for i in range(n)]],
                                         dtype=complex)
                        rest[0,m]-=1
                        rest[0,:]*=np.exp(1j*wavenumber*graph.cells[c.cofaces[m].index].length)
                    cofaces.append(SheafCoface(index=c.cofaces[m].index,
                                               orientation=c.cofaces[m].orientation,
                                               restriction=rest))
            else: # All other faces have trivial restrictions
                n=2
                cofaces=[SheafCoface(index=cf.index,orientation=cf.orientation,restriction=LinearMorphism([])) for cf in c.cofaces]
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
                              cofaces=[SheafCoface(index=cf.index, 
                                                   orientation=cf.orientation,
                                                   restriction=np.matrix(1))
                                       for cf in c.cofaces],
                              stalkDim=1)
                    for c in cells]

        Sheaf.__init__(self,sheafcells)
        return

class SheafMorphismCell:
    def __init__(self,destinations=None,maps=None):
        """Specify destinations and maps for this cell's stalk under a morphism"""
        self.destinations=destinations
        self.maps=maps

class SheafMorphism:
    def __init__(self,morphismCells):
        """Construct a sheaf morphism as a list of component maps specified as SheafMorphismCells"""
        self.morphismCells=morphismCells

    def __mul__(self,other):
        """Composition of two sheaf morphisms"""
        morCells=[SheafMorphismCell(destinations=[selfdest for otherdest in mc.destinations for selfdest in self.morphismCells[otherdest].destinations],
                                    maps=[selfmap*othermap for otherdest,othermap in zip(mc.destinations,mc.maps) for selfdest,selfmap in zip(self.morphismCells[otherdest].destinations,self.morphismCells[otherdest].maps)]) for mc in other.morphismCells]
        return SheafMorphism(morCells)

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
            return 

        # Is the desired cell a coface of a cell in the support?
        for s in self.sectionCells:
            for cf in sheaf.cells[s.support].cofaces:
                if cf.index == cell:
                    # If so, extend via restriction
                    val=cf.restriction(s.value)

                    # Check for consistency
                    if value is not None and np.any(np.abs(val - value)>tol):
                        return False
                    value = val
                            
        # Are there are any cofaces for the desired cell in the support?
        if value is None: # Attempt to assign a new value...
            
            # Stack the restrictions and values associated to existing support
            lst=[(cf.restriction.matrix,s.value) 
                 for cf in sheaf.cells[cell].cofaces 
                 for s in self.sectionCells
                 if isinstance(cf.restriction,LinearMorphism) and (cf.index == s.support)]
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
                        if np.any(np.abs(cf.restriction(value)-s.value)>tol):
                            return False
        
        # A value was successfully assigned (if no value was assigned, 
        # do nothing, but it's still possible to extend)
        if value is not None:
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
                    cofaces.append(SheafCoface(index=d,
                                               orientation=1,
                                               restriction=inducedMap(sheaves[i],sheaves[d],mor,k)))
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
        return np.matrix([])

    # Extract the k-skeleta of each sheaf
    k_1,ksizes_1,kidx_1=sheaf1.kcells(k,compactSupport)
    k_2,ksizes_2,kidx_2=sheaf2.kcells(k,compactSupport)

    # Construct chain map
    rows=sum(ksizes_2)
    cols=sum(ksizes_1)
    m=np.zeros((rows,cols),dtype=np.complex)

    for i in range(len(k_1)):
        for j,map in zip(morphism.morphismCells[k_1[i]].destinations,morphism.morphismCells[k_1[i]].maps):
            if sheaf2.cells[j].dimension==k:
                ridx=[q for q in range(len(k_2)) if k_2[q]==j]
                if ridx:
                    ridx=ridx[0]
                    m[kidx_2[ridx]:kidx_2[ridx+1],kidx_1[i]:kidx_1[i+1]]+=map.matrix

    # Map basis for domain sheaf's cohomology through the chain map
    im=np.dot(m,Hk_1)
    
    # Expand output in new basis
    map,j1,j2,j3 = np.linalg.lstsq(im,Hk_2)

    return map.conj().T

def simplexOrientation(s1,s2):
    """Assuming s1 is a face of s2, what's its orientation?"""
    if not set(s1).issubset(s2):
        return 0 # s1 is not a face of s2

    # Find what element(s) are missing in s1 and the orientation arising from those missing slots
    orientation=1
    for i,s in enumerate(s2):
        if not s in s1:
            orientation *= (-1)**i

    # Find permutation of the remaining entries
    orientation*=perm_parity([s1.index(s) for s in s2 if s in s1])

    return orientation

def perm_parity(lst):
    '''\
    Given a permutation of the digits 0..N in order as a list, 
    returns its parity (or sign): +1 for even parity; -1 for odd.
    '''
    parity = 1
    for i in range(0,len(lst)-1):
        if lst[i] != i:
            parity *= -1
            mn = min(range(i,len(lst)), key=lst.__getitem__)
            lst[i],lst[mn] = lst[mn],lst[i]
    return parity

def ksublists(lst,n,sublist=[]):
    """Iterate over all ordered n-sublists of a list lst"""
    if n==0:
        yield sublist
    else:
        for idx in range(len(lst)):
            item=lst[idx]
            for tmp in ksublists(lst[idx+1:],n-1,sublist+[item]): 
                yield tmp

def ksimplices(toplexes,k,relative=None):
    """List of k-simplices in a list of toplexes"""
    simplices=[]
    for toplex in toplexes:
        for spx in ksublists(toplex,k+1):
            if not spx in simplices and (relative is None or not spx in relative):
                simplices.append(spx)
    return simplices 
