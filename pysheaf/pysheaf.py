# Persistence-capable sheaf manipulation library
#
# Copyright (c) 2013-2018, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied

import numpy as np
import random
# import matplotlib.pyplot as plt
import networkx as nx
import scipy.optimize

import warnings
import copy
import time

#import all DEAP related functions for genetic algorithms
from functools import partial

from deap import base
from deap import creator
from deap import tools
from deap import algorithms

from collections import defaultdict



## Data structures
class Coface:
    """A coface relation"""
    def __init__(self,index=None,orientation=None):
        self.index=index
        self.orientation=orientation

    def __repr__(self):
        return "(index=" + str(self.index) + ",orientation="+str(self.orientation)+")"

class Cell:
    """A cell in a cell complex.  The id attribute can be used as a key to index into the registry of cells or cofaces"""
    def __init__(self,dimension,compactClosure=True,cofaces=None, name=None, id=None):
        self.dimension=dimension
        self.compactClosure=compactClosure
        self.id=id
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

    def isCoface(self,index,orientation=None):
        """Check if a given cell index is a *listed* coface of this cell. Optionally check that the orientation is as given."""
        if orientation is None:
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

        # Register cell IDs with indices into the self.cells list
        self.cell_dict={}
        self.coface_dict={}
        for i,c in enumerate(self.cells):
            if c.id is None: # Build IDs if not present
                c.id = str(i)
            self.cell_dict[c.id]=i
            if c.cofaces is not None: # Register cofaces if present
                for j,cf in enumerate(c.cofaces):
                    self.coface_dict[(i,cf.index)]=j

        # Attribute for counting cells later on...
        self.cell_counter=len(self.cells)

    def add_cell(self,cell):
        """Add a cell to the cell complex, using the id attribute for registering the cell.  Returns the id attribute for later use"""

        if cell.id is None: # Construct an ID if needed
            cell.id = str(self.cell_counter)
            self.cell_counter += 1

        # Register the cell
        self.cell_dict[cell.id]=len(self.cells)

        # Register its cofaces, if any
        if cell.cofaces is not None:
            for j,cf in enumerate(cell.cofaces):
                self.coface_dict[(len(self.cells),cf.index)]=j

        # Store the cell
        self.cells.append(cell)

        # Return the cell's ID for later use
        return cell.id

    def add_cells_from(self,cells):
        for c in cells:
            self.add_cell(c)

    def add_coface(self,cellpair,orientation):
        """Add a coface to the cell complex, referenced by a pair of cell IDs.  The cellpair argument is assumed to be a pair: (face,coface).  If the cells aren't present, this will raise KeyError."""
        # Look up which cells are involved...
        source=self.cell_dict[cellpair[0]]
        target=self.cell_dict[cellpair[1]]

        # Drop in the coface
        self.coface_dict[(source,target)]=len(self.cells[source].cofaces)
        self.cells[source].cofaces.append(Coface(index=target,orientation=orientation))

    def add_cofaces_from(self,cellpairs,orientations):
        for cellpair,orientation in zip(cellpairs,orientations):
            self.add_coface(cellpair,orientation)

    def isFaceOf(self,c,cells=None):
        """Construct a list of indices of all cells that a given cell is a face of = indices of all cofaces"""
        return list(set(cf.index for cf in self.cofaces(c,cells)))

    def skeleton(self,k,compactSupport=False):
        """Return the k-skeleton of a cell complex.  Optionally ensure that the complex returned is closed."""
        return [i for i in range(len(self.cells))
                if ((compactSupport or self.cells[i].compactClosure) and self.cells[i].dimension==k)]

    def faces(self,c):
        """Compute a list of indices all faces of a cell"""
        return [i for i in range(len(self.cells)) if c in set(cf.index for cf in self.cofaces(i))]

    def closure(self,cells):
        """Compute the closure of a collection of cells"""
        return [i for i in range(len(self.cells))
                if i in cells or [cf for cf in self.cofaces(i) if cf.index in cells]]

    def interior(self,cells):
        """Compute the interior of a collection of cells"""
        return [i for i in cells if set(self.starCells([i])).issubset(cells)]

    def cofaces(self,c,cells=[]):
        """Iterate over cofaces (of all dimensions) of a given cell; optional argument specifies which cells are permissible cofaces
        Warning: duplicates are possible!"""
        cells = cells or []
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
        current_cpt = current_cpt or [start]
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
        cflist=[cf.index for cf in self.cofaces(start)]
        if not cells:
            return list(set(self.faces(start) + cflist))
        else:
            return list(set.intersection(set(self.faces(start) + cflist),cells))

    def starCells(self,cells):
        """Cells in star over a subset of a cell complex"""
        return set(cells).union({cf.index for c in cells for cf in self.cofaces(c)})

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
            map,j1,j2,j3=np.linalg.lstsq(ker,dp1,rcond=None)
            Hk=np.dot(ker,cokernel(map,tol));
        else:
            Hk=ker

        return Hk

    def betti(self,k,compactSupport=False,tol=1e-5):
        """Compute the k-th Betti number of the cell complex"""
        return self.homology(k,compactSupport).shape[1]

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
    def __init__(self,dimension,cofaces=None,compactClosure=True,stalkDim=None,metric=None,name=None,id=None,bounds=None):
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

        self.bounds=bounds

        if metric is not None:
            self.metric=metric
        else:
            self.metric=lambda x,y: np.linalg.norm(x-y)

        Cell.__init__(self,dimension,compactClosure,cofaces,name,id)

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
        self.coface_dict[(source,target)]=len(self.cells[source].cofaces)
        self.cells[source].cofaces.append(SheafCoface(index=target,orientation=orientation,restriction=restriction))

    def add_cofaces_from(self,cellpairs,orientations,restrictions):
        for cellpair,orientation,restriction in zip(cellpairs,orientations,restrictions):
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

    def cobetti(self,k,compactSupport=False,tol=1e-5):
        """Compute the k-th Betti number of the sheaf"""
        return self.cohomology(k,compactSupport).shape[1]

    def maximalExtend(self,assignment,multiassign=False,tol=1e-5):
        """Take a partial assignment and extend it to a maximal assignment that's non-conflicting (if multiassign=False) or one in which multiple values can be given to a given cell (if multiassign=True)"""
        for i in range(len(assignment.sectionCells)):
            for cf in self.cofaces(assignment.sectionCells[i].support):
                if multiassign or (not assignment.extend(self,cf.index,tol=tol)):
                    assignment.sectionCells.append(SectionCell(cf.index,cf.restriction(assignment.sectionCells[i].value),source=assignment.sectionCells[i].support))
        return assignment

    def consistencyRadii(self,assignment,testSupport=None,consistencyGraph=None,tol=1e-5):
        """Compute all radii for consistency across an assignment"""
        if testSupport is None:
            cellSet=set(range(len(self.cells)))
        else:
            cellSet=set(testSupport)

        if consistencyGraph is None:
            cG=self.consistencyGraph(assignment,testSupport)
        else:
            cG=consistencyGraph

        return np.unique([rad for (i,j,rad) in cG.edges(nbunch=cellSet,data='weight')])

    def consistencyRadius(self,assignment,testSupport=None,consistencyGraph=None,ord=np.inf,tol=1e-5):
        """Compute the consistency radius of an approximate section"""

        if testSupport is None:
            cellSet=set(range(len(self.cells)))
        else:
            cellSet=set(testSupport)

        if consistencyGraph is None:
            cG=self.consistencyGraph(assignment,testSupport)
        else:
            cG=consistencyGraph

        radii=[rad for (i,j,rad) in cG.edges(nbunch=cellSet,data='weight')]
        if len(radii) == 0:
            warnings.warn("No SectionCells in the assignments match, therefore nothing was compared by consistencyRadius")
            return 0.
        return np.linalg.norm(radii,ord=ord)    
    
    def maxTestSupport(self,activeCells):
        elementsTestSupport = copy.deepcopy(activeCells)
        for i in range(len(elementsTestSupport)):
            for cf in self.cofaces(elementsTestSupport[i]):
                if not cf.index in elementsTestSupport:
                    elementsTestSupport.append(cf.index)
        return elementsTestSupport

    def consistencyGraph(self,assignment,testSupport=None):
        """Construct a NetworkX graph whose vertices are cells, in which each edge connects two cells with a common coface weighted by the distance between their respective values.  Note: the assignment must be supported on the entire space. Edges also have a type attribute, explaining the kind of relationship between cells: (1 = two values assigned to this cell, 2 = one cell is a coface of the other, 3 = cells have a common coface.)"""
        
        if testSupport is None:
            cellSet=set(range(len(self.cells)))
        else:
            cellSet=set(testSupport)

        G=nx.Graph()
        G.add_nodes_from(cellSet)
        
        for c1 in assignment.sectionCells:
            if (testSupport is None) or ((c1.support in testSupport) and (c1.source in testSupport)):
                for c2 in assignment.sectionCells:
                    if (testSupport is None) or (c2.source in testSupport):
                        if c1.support == c2.support:
                            rad=self.cells[c1.support].metric(c1.value,c2.value)
                            if ((c1.support,c2.support) not in G.edges()) or G[c1.support][c2.support]['weight'] < rad:
                                G.add_edge(c1.support,c2.support,weight=rad,type=1)
                        else:
                            for cf1 in self.cofaces(c1.support):
                                if cf1.index == c2.support:
                                    rad=self.cells[cf1.index].metric(cf1.restriction(c1.value),c2.value)
                                    if ((c1.support,c2.support) not in G.edges()) or G[c1.support][c2.support]['weight'] < rad:
                                        G.add_edge(c1.support,c2.support,weight=rad,type=2)
                                else:
                                    for cf2 in self.cofaces(c2.support):
                                        if cf1.index == cf2.index:
                                            rad=0.5*self.cells[cf1.index].metric(cf1.restriction(c1.value),cf2.restriction(c2.value)) # Note the factor of 0.5
                                            if ((c1.support,c2.support) not in G.edges()) or (G[c1.support][c2.support]['type'] == 2 and G[c1.support][c2.support]['weight'] < rad):
                                                G.add_edge(c1.support,c2.support,weight=rad,type=3)

        return G

    def minimizeConsistencyRadius(self,assignment, activeCells=None, testSupport=None, method='nelder-mead', ord = np.inf, options={}, tol=1e-5):
        """
        Minimize consistency radius of an assignment given fixed cells
        Currently, any optimization supported by scipy.optimize.minimize
            Parameters:
                assignment: the partial assignment of the sheaf to fuse
                activeCells: set of cells whose values are to be changed (if None, all cells outside the support of the assignment will be changed, but nothing in the support of the assignment will be changed)
                testSupport: the set of cells over which consistency radius is assessed
                tol: the tol of numeric values to be considered the same
        """
        if activeCells is None:
            support=[sc.support for sc in assignment.sectionCells]
            ac=[idx for idx in range(len(self.cells)) if idx not in support]
        else:
            ac=activeCells
        
        if method == 'KernelProj':
            if not self.isLinear():
                raise NotImplementedError('KernelProj only works for sheaves of vector spaces')

            if ord != 2:
                warn('Kernel projection requires order 2 in minimizeConsistencyRadius')
                
            # Compile dictionary of rows
            rowstarts=dict()
            rowidx=0
            for i in support:
                rowstarts[i]=rowidx
                rowidx+=self.cells[i].stalkDim
                
            newassignment = Section([sc for sc in assignment.sectionCells])
                
            # Optimize each active cell independently
            for i in ac:
                if self.cells[i].stalkDim > 0:
                    # Matrix of all restrictions out of this cell into the support
                    mat=np.zeros((sum([self.cells[j].stalkDim for j in support]),
                                  self.cells[i].stalkDim))

                    for cf in self.cofaces(i): # Iterate over all cofaces of this activeCell
                        try:
                            supportidx=support.index(cf.index)
                            mat[rowstarts[supportidx]:rowstarts[supportidx]+self.cells[supportidx].stalkDim,:]=cf.restriction.matrix
                        except ValueError:
                            pass
            
                # Use least squares to solve for assignment rooted at this cell given the existing assignment
                asg,bnds=self.serializeAssignment(assignment,activeCells=support) # Confusingly, activeSupport here refers *only* to the support of the assignment
                result=np.linalg.lstsq(mat,asg,rcond=None)
                
                newassignment.sectionCells.append(SectionCell(i,result[0]))

            return newassignment
        elif self.isNumeric():
            initial_guess, bounds = self.serializeAssignment(assignment,ac)
            res=scipy.optimize.minimize( fun = lambda sec: self.consistencyRadius(self.deserializeAssignment(sec,ac,assignment), testSupport=testSupport, ord=ord),
                                         x0 = initial_guess,
                                         method = method, 
                                         bounds = bounds,
                                         tol = tol,
                                         options = {'maxiter' : int(100)})
            newassignment = self.deserializeAssignment(res.x,ac,assignment)
            return newassignment
        else:
            raise NotImplementedError('Non-numeric sheaf')
        return newassignment

    def consistentPartition(self,assignment,threshold,testSupport=None,consistencyGraph=None,ord=np.inf,tol=1e-5):
        """Construct a maximal collection of subsets of cells such that each subset is consistent to within the given threshold.  Note: the assignment must be supported on the entire space."""

        if testSupport is None:
            cellSet=set(range(len(self.cells)))
        else:
            cellSet=set(testSupport)

        if consistencyGraph is None:
            cG=self.consistencyGraph(assignment,testSupport)
        else:
            cG=consistencyGraph

        # Construct inconsistency graph.  Edges indicate cells that cannot be in the same cover element
        G=nx.Graph()
        G.add_nodes_from(cellSet)
        G.add_edges_from([(i,j) for (i,j,k) in cG.edges(nbunch=cellSet,data='weight') if k > threshold and i != j])
 
        # Solve graph coloring problem: nodes (cells) get colored by cover element
        color_dict=nx.coloring.greedy_color(G)
 
        # Re-sort into groups of consistent cells
        cdd=defaultdict(list)
        for cell,cover_element in color_dict.items():
            cdd[cover_element].append(cell)

        return {frozenset(s) for s in cdd.values()}

    def consistentCollection(self,assignment,threshold,testSupport=None,consistencyGraph=None,ord=np.inf,tol=1e-5):
        """Construct a maximal collection of open sets such that each subset is consistent to within the given threshold.  Note: the assignment must be supported on the entire space."""
        # First obtain a collection of consistent open sets.  These are disjoint
        initial_collection={frozenset(self.interior(s)) for s in self.consistentPartition(assignment,threshold,testSupport,consistencyGraph,ord=ord,tol=tol)}

        additions = True
        collection = set()

        while additions:
            additions=False
            for u in initial_collection:
                added_u=False
                for v in initial_collection:
                    if u is not v:
                        u_v = u.union(v)
                        if self.consistencyRadius(assignment,testSupport=u_v,ord=ord)<threshold:
                            added_u=True
                            additions=True
                            collection.add(u_v)
                if not added_u:
                    collection.add(u)
                initial_collection=collection
                collection=set()

        return initial_collection

    def assignmentMetric(self,assignment1,assignment2, testSupport=None, ord=np.inf):
        """Compute the distance between two assignments"""
        radii=[]
        for c1 in assignment1.sectionCells:
            if (testSupport is None) or ((c1.support in testSupport) and (c1.source in testSupport)):
                for c2 in assignment2.sectionCells:
                    if c1.support == c2.support and ((testSupport is None) or (c2.source in testSupport)):
                        radii.append(self.cells[c1.support].metric(c1.value,c2.value))
    
        if len(radii) == 0:
            radii.append(np.inf)
            warnings.warn("No SectionCells in the assignments match, therefore nothing was compared by assignmentMetric")
        return np.linalg.norm(radii,ord=ord)
    
    def fuseAssignment(self,assignment, activeCells=None, testSupport=None, method='SLSQP', options={}, ord=np.inf, tol=1e-5):
        """
        Compute the nearest global section to a given assignment
        Currently there are three optimization schemes to choose from
        'KernelProj': This algorithm only works for sheaves of vector spaces and uses the kernel projector for the 0-coboundary map
        'SLSQP': This algorithm is scipy.optimize.minimize's default for bounded optimization
            Parameters:
                assignment: the partial assignment of the sheaf to fuse
                activeCells: set of cells whose values are to be changed or None if all are allowed to be modified (Note: global sections may not be returned if this is changed from the default)
                testSupport: the set of cells over which consistency radius is assessed
                tol: the tol of numeric values to be considered the same
        'GA': This genetic algorithm was implemented using DEAP for optimizations over nondifferentiable functions
            Parameters:
                assignment: the partial assignment of the sheaf to fuse
                tol: the tol of numeric values to be considered the same
                options: a dictionary to store changes to parameters, the keys must be identical to the current parameters
                    keys for GA:
                        initial_pop_size - the number of individuals in the starting population
                        mutation_rate - the proportion of the offspring (newly created individuals each round) that are from mutations rather
                                        than mating
                        num_generations - the number of iterations that the genetic algorithm runs
                        num_ele_Hallfame - the number of top individuals that should be reported in the hall of fame (hof)
        """
        if activeCells is not None and testSupport is not None and not (set(testSupport) <= set(self.maxTestSupport(activeCells))):
            raise ValueError("Given testSupport is larger than the largest comparison set given activeCells")
        
        if method == 'SLSQP':
            if self.isNumeric():
                globalsection = self.optimize_SLSQP(assignment, activeCells, testSupport, ord, tol)
            else:
                # The fallback situation, where we need to iterate over global sections manually...
                raise NotImplementedError
        elif method == 'GA':
            add_parameters = {'initial_pop_size':100, 'mutation_rate':0.3, 'num_generations':100, 'num_ele_Hallfame':1, 'initial_guess_p':None}
            overlap = [st for st in add_parameters.keys() if st in set(options.keys())]
            if len(overlap) > 0:
                for st_overlap in overlap:
                    add_parameters[st_overlap] = options[st_overlap]
            globalsection = self.optimize_GA(assignment, tol, activeCells=activeCells, testSupport=testSupport, initial_pop_size=add_parameters['initial_pop_size'], mutation_rate = add_parameters['mutation_rate'], num_generations= add_parameters['num_generations'] ,num_ele_Hallfame=add_parameters['num_ele_Hallfame'], initial_guess_p = add_parameters['initial_guess_p'])
        elif method == 'KernelProj':
            if not self.isLinear():
                raise NotImplementedError('KernelProj only works for sheaves of vector spaces')

            if ord != 2:
                warn('Kernel projection requires order 2 in fuseAssignment')
            
            # Construct the coboundary map
            d = self.coboundary(k=0)
            ks,ksizes,kidx=self.kcells(k=0)
            asg,bounds = self.serializeAssignment(assignment,activeCells=ks)
            
            # Construct the kernel projector
            AAt=np.dot(d,d.transpose())
            AAti=np.linalg.pinv(AAt)
            AtAAti=np.dot(d.transpose(),AAti)
            projector=np.eye(d.shape[1])-np.dot(AtAAti,d)
            
            # Apply the kernel projector
            gs=np.dot(projector,asg)

            # Deserialize
            globalsection=self.deserializeAssignment(gs,activeCells=ks,assignment=assignment)
        else:
            raise NotImplementedError('Invalid method')
        return globalsection
    
    
    def deserializeAssignment(self,vect,activeCells=None,assignment=None):
        """Transform a vector of values for a numeric-valued sheaf into an assignment as a Section instance(Note: this is really a helper method and should generally not be used by external callers).  Inactive cells are filled from another assignment that's optionally supplied"""
        if not self.isNumeric():
            raise TypeError('Cannot deserialize an assignment vector for a non-numeric sheaf')
    
        scs=[]
        idx=0

        for i in range(len(self.cells)):
            if (activeCells is None) or (i in activeCells): # If the cell is active, pull its value from the vector
                if self.cells[i].stalkDim > 0:
                    scs.append(SectionCell(support=i,value=vect[idx:idx+self.cells[i].stalkDim]))
                    idx+=self.cells[i].stalkDim
            elif assignment is not None: # If the cell is not active, pull its value from the given assignment
                for cell in assignment.sectionCells:
                    if i == cell.support:
                        scs.append(cell)
                    
        return Section(scs)
    
    def deserializeAssignment_ga(self, vect, id_len):
        #write a new assignment from the individual so that one can use maximal extend.
        new_assignment = []
        start_index = 0
        for i in range(len(id_len)):
            new_assignment.append(SectionCell(support=id_len[i][0], value=vect[(start_index):(start_index+id_len[i][1])]))
            start_index += id_len[i][1]
        
        new_assignment = Section(new_assignment)
        
        return new_assignment

    def serializeAssignment(self,assignment,activeCells=None):
        """Transform a partial assignment in a Section instance into a vector of values"""
        if not self.isNumeric():
            raise TypeError('Cannot serialize an assignment vector for a non-numeric sheaf')

        x0 = np.zeros((sum([c.stalkDim for i,c in enumerate(self.cells) if ((activeCells is None) or (i in activeCells))])),dtype=assignment.sectionCells[0].value.dtype)

        # If any components are bounded, collect the bounds (later)
        bounded=False
        bounds=None
        if activeCells is None:
            activeCells = range(len(self.cells))
        for i in activeCells:
            if self.cells[i].bounds is not None:
                bounded=True
                bounds=[]
                break

        # Figure out cell boundaries in vector, and collect their bounds if present
        idx=0
        idxarray=[]
        for i in activeCells:
            if self.cells[i].stalkDim > 0:
                idxarray.append(idx)
                for cell in assignment.sectionCells:         # Pack data into vector.  If there are multiple values assigned, the one appearing last is used
                    if cell.support == i: 
                        x0[idx:idx+self.cells[i].stalkDim]=cell.value
                idx+=self.cells[i].stalkDim
                if bounded:
                    if self.cells[i].bounds is None:
                        bounds+= [(None,None)]*self.cells[i].stalkDim
                    else:
                        bounds.extend(self.cells[0].bounds)
        
        if bounded:
            bounds = tuple(bounds)

        return x0,bounds
                
    
    def optimize_SLSQP(self, assignment, activeCells=None, testSupport=None, ord=np.inf, tol=1e-5):
        """
        Compute the nearest global section to a given assignment using 
        scipy.optimize.minimize. When there are constraints specified, 
        scipy.optimize.minimize defaults to using SLSQP.
        Based on:
            Kraft, D. A software package for sequential quadratic programming. 1988. 
            Tech. Rep. DFVLR-FB 88-28, DLR German Aerospace Center, 
            Institute for Flight Mechanics, Koln, Germany.
        """
        initial_guess, bounds = self.serializeAssignment(assignment,activeCells)
        res=scipy.optimize.minimize( fun = lambda sec: ((self.assignmentMetric(assignment,self.deserializeAssignment(sec,activeCells,assignment), testSupport=testSupport, ord=ord))),
                                    x0 = initial_guess,
                                    method = 'SLSQP', 
                                    bounds = bounds,
                                    constraints = ({'type' : 'eq',
                                                    'fun' : lambda asg: self.consistencyRadius(self.deserializeAssignment(asg,activeCells,assignment),testSupport=testSupport, ord=ord)}), 
                                    tol = tol, 
                                    options = {'maxiter' : int(100)})
        globalsection = self.deserializeAssignment(res.x,activeCells,assignment)
        return globalsection
    
    def ga_optimization_function(self,  space_des_2_opt, assignment,individual,activeCells=None,testSupport=None, tol=1e-5):
        """Write the function for the genetic algorithm to optimize similar to fun for scipy.optimize.minimize"""
        
        #Assign a high cost to eliminate individuals with nans from the population
        if np.any(np.isnan(individual)):
            cost = 1e100
            cost = -cost
            
        else:
            #write a new assignment from the individual so that one can use maximal extend.
            multiassign = False
            
            new_assignment = []
            start_index = 0
            for i in range(len(space_des_2_opt)):
                new_assignment.append(SectionCell(support=space_des_2_opt[i][0], value=individual[(start_index):(start_index+space_des_2_opt[i][1])]))
                start_index += space_des_2_opt[i][1]
                
                
            if activeCells is not None:
                for sec in assignment.sectionCells:
                    if not (sec.support in activeCells):
                        new_assignment.append(sec)
                        multiassign=True
                
        
            new_assignment = Section(new_assignment)
        
            #start of optimization function
            new_assignment = self.maximalExtend(new_assignment,multiassign=multiassign,tol=1e-5)
        
            cost = self.assignmentMetric(assignment, new_assignment, testSupport=testSupport)
            
            #If consistency is not guarenteed force consistency check
            #TBD: verify neccessity
            #if activeCells != None:
            #Constrain value by the consistencyRadius (Note: min is 0)
            t0 = time.clock()
            radii = self.consistencyRadius(new_assignment)
            t1 = time.clock()
            
            t_finish = t1-t0
            print t_finish
            
            if radii < tol:
                pass
            else:
                cost = cost+(radii)**2   
            
            
            #Sum to formulate cost (NOTE: GA maximizes instead of minimizes so we need a negative sign)
            cost = -cost
        return (float(cost),)
        

    
    def optimize_GA(self, assignment, tol=1e-5, activeCells=None, testSupport=None, initial_pop_size=100, mutation_rate = .3, num_generations=100 ,num_ele_Hallfame=1, initial_guess_p = None):
        """
        Compute the nearest global section to a given assignment using 
        a genetic algorithm. 
        The current implementation forces the assignment to have bounds on every element
        
        Inputs:
            assignment - the initial assignment that the selection function compares guess against for consistency
            tol - the maximum distance for two items to be considered the same
            initial_pop_size - the number of individuals in the starting population
            mutation_rate - the proportion of the offspring (newly created individuals each round) that are from mutations rather
                            than mating
            num_generations - the number of iterations that the genetic algorithm runs
            num_ele_Hallfame - the number of top individuals that should be reported in the hall of fame (hof)
            
        Outputs:
            pop - the final population of the algorithm
            stats - the statistics on the fitness of the final population
            hof - the fittest individual
            opt_sp_id_len - the space the algorithm is optimizing over
        
        """
        ##########################################################################################
        # Helper Functions for Algorithm
        ##########################################################################################
        
        #add as a potential selection function to match Matlab GA
        def selnormGeom(individuals, k, prob_sel_best= 0.08, fit_attr="fitness"):
            #NormGeomSelect is a ranking selection function based on the normalized
            #geometric distribution.  
            
            #Modified from the Matlab version into the style of DEAP

            q = prob_sel_best   # Probability of selecting the best
            n = len(individuals)  # Number of individuals in pop
            
            
            chosen = [] #editted the structure of th output to reflect the structure of pysheaf
            fit = np.zeros((n,1))  #Allocates space for the prop of select
            x = np.zeros((n,2))    #Sorted list of id and rank
            x[:, 0] = range(n,0,-1) # Get original location of the element
            to_sort = zip(individuals, range(n)) #need to keep original index associated
            s_inds = sorted(to_sort, key= lambda ind: getattr(ind[0], fit_attr).values[0]) #Sort by the fitnesses
            x[:, 1] = [b for a,b in s_inds]
            r =q/(1-((1-q)**n))  # normalize the distribution, q prime
            for ind in range(n):  #Generate the probability of selection
                ind_fit = int(x[ind,1])
                fit[ind_fit] = r*((1-q)**(x[ind, 0]-1))
            fit = np.cumsum(fit)  # Calculate the cummulative prob. function
            rnums = sorted([random.random() for nn in range(n)])  # Generate n sorted random numbers
            fitIn = 0
            new_In = 0
            unique = []
            while new_In < k:
                if rnums[new_In] < fit[fitIn]:
                    unique.append(fitIn)
                    chosen.append(individuals[fitIn]) #Select the fitIn individual
                    new_In += 1  # Looking for next new individual
                else:
                    fitIn += 1 # Looking at next potential selection
            
            return chosen 
        
        #Define functions similar to Matlab optimization functionality
        def arithXover(ind1, ind2):
            for i, (x1, x2) in enumerate(zip(ind1, ind2)):
                gamma = random.random()
                ind1[i] = (1. - gamma) * x1 + gamma * x2
                ind2[i] = gamma * x1 + (1. - gamma) * x2

            return ind1, ind2
            
        #Define a function to do the mutation of the elements from a floating point perspective
        def mutUniformFloat(individual, low, up, indpb):
            """Mutate an individual by replacing attributes, with probability *indpb*,
            by a integer uniformly drawn between *low* and *up* inclusively.
            
            :param individual: :term:`Sequence <sequence>` individual to be mutated.
            :param low: The lower bound or a :term:`python:sequence` of
                        of lower bounds of the range from wich to draw the new
                        integer.
            :param up: The upper bound or a :term:`python:sequence` of
                       of upper bounds of the range from wich to draw the new
                       integer.
            :param indpb: Independent probability for each attribute to be mutated.
            :returns: A tuple of one individual.
            """
            size = len(individual)
            if len(low) < size:
                raise IndexError("low must be at least the size of individual: %d < %d" % (len(low), size))
            elif len(up) < size:
                raise IndexError("up must be at least the size of individual: %d < %d" % (len(up), size))
            
            for i, xl, xu in zip(xrange(size), low, up):
                if random.random() < indpb:
                    individual[i] = random.uniform(xl, xu)
            
            return individual,
        
        ###################################################################################################
        #Beginning of Algoithm Run
        ###################################################################################################
        
        
        
        
        #Get the spaces to optimize over while saving their indices in the complex and length as well as the bounds
        bounds = []
        opt_sp_id_len = [] #stores the cell id and 
        
        for cell in self.cells:
            if activeCells is None:
                #Use the 0-dimension cells in sheaf
                active = not any([alt_cell.isCoface(int(cell.id)) for alt_cell in self.cells if alt_cell.id != cell.id]) #should be done differently for speed
            else:
                #Use the cells denoted as active in sheaf
                active = False
                if int(cell.id) in activeCells:
                    active = True
                    
            if active:
                opt_sp_id_len.append((int(cell.id), cell.stalkDim))
                if cell.bounds !=None:
                    #Ensure that the length of the specfied bounds is equivalent to the stalkDim
                    if len(cell.bounds) == cell.stalkDim:
                        bounds.extend(cell.bounds)
                    else:
                        raise ValueError("Not all bounds specified")
                else:
                    bnds=[]
                    #This will not work (error out if stalkDim = None) however why are you trying to optimize over an empty cell
                    for x in range(cell.stalkDim):
                        bnds.append(tuple([None, None]))
                    bounds.extend(bnds)            
        
        
        #Create initial guess if that index is specfied for that section
        initial_guess = [[0 for j in range(opt_sp_id_len[i][1])] for i in range(len(opt_sp_id_len))]
        for section in assignment.sectionCells:
            id_list = [opt_sp_id_len[i][0] for i in range(len(opt_sp_id_len))]
            # all ids should be unique, enabling the the use of index
            try:
                ind = id_list.index(section.support)
                #This needs to be double checked that it is working properly
                #section.value needs to be a a numpy array
                initial_guess[ind] = section.value
                
            except ValueError:
                pass
            
        
        initial_guess_restructure = []
        for i in range(len(initial_guess)):
            initial_guess_restructure.extend(initial_guess[i])
        initial_guess = np.array(initial_guess_restructure)
        
        
        #Make the Initial guess none if section was empty
        if np.array_equal(initial_guess, np.zeros_like(initial_guess)):
            initial_guess = None
            
        #Ensure that any assigned individual to the initial population is the correct length
        #Needs to be modified to take either an array or an set of Sections
        
        #Deserialize initial_guess_p if it is expressed as a section
        if isinstance(initial_guess_p, Section):
            setOfArrays = [[] for r in range(len(opt_sp_id_len))]
            setActiveCells = [s[0] for s in opt_sp_id_len]
            for sec in initial_guess_p.sectionCells:
                if sec.support in setActiveCells:
                    ind_opt = (np.abs(np.array(setActiveCells) - sec.support)).argmin()
                    if np.size(sec.value) == opt_sp_id_len[ind_opt][1]:
                        setOfArrays[ind_opt] = sec.value
                    else:
                        initial_guess_p = None
                        break
                else:
                    warnings.warn("initial_guess_p does not contain a sectionCell for every activeCell")
                    initial_guess_p = None
                    break
            
            #Pull all the section values into a single array
            #TBD: Try to streamline
            if initial_guess_p != None:
                initial_guess_p = np.array([])
                for a_i in setOfArrays:
                    initial_guess_p = np.hstack((initial_guess_p, a_i))
                
        else:
            #Keep original functionality
            if np.any(initial_guess_p):
                if np.size(initial_guess_p) != self.cells[0].stalkDim:
                    initial_guess_p = None
        
        
        #Start of unique to the genetic algorithm
        creator.create("FitnessMax", base.Fitness, weights=(1.0,))
        creator.create("Individual", np.ndarray, fitness=creator.FitnessMax)
        
        toolbox = base.Toolbox()

        
        
        #seq_func forces the initial population to be generated randomly within the bounds
        seq_func = []
        for bnds in bounds:
            if bnds[0] != None and bnds[1] != None:
                seq_func.append(partial(random.uniform, bnds[0], bnds[1]))
            elif bnds[0] == None and bnds[1] == None:
                seq_func.extend([lambda:(1/(1-random.random()))-1]) #maps [0,1) to [0, inf)
            elif bnds[0] == None:
                multiply_bnds1 = partial(np.multiply, bnds[1])
                seq_func.extend([lambda:(-1/(1-random.random()) + 1 + multiply_bnds1(1.0))])
                #seq_func.extend([lambda:(random.randrange(copy.deepcopy(bnds[1])))]) #need to check actual opperation of randrange without a start
            else:
                multiply_bnds0 = partial(np.multiply, bnds[0])
                seq_func.extend([lambda: (1/(1-random.random()))-1 + multiply_bnds0(1.0)])
                #seq_func.extend([lambda:(-1*random.randrange(copy.deepcopy(bnds[0])) + 2*copy.deepcopy(bnds[0]))])
        
        #specify a population within the bounds
        toolbox.register("individual", tools.initCycle, creator.Individual, seq_func, n=1)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        
        

            
        #Include the initial guess in the population, otherwise register the population
        if not np.any(initial_guess) and not np.any(initial_guess_p):
            #specify a population without including an initial guess
            pop = toolbox.population(n=initial_pop_size)
        elif not np.any(initial_guess) and np.any(initial_guess_p):
            pop = toolbox.population(n=(initial_pop_size-1))
            initial_g = creator.Individual(initial_guess_p)
            pop.insert(0, initial_g)
        elif np.any(initial_guess) and not np.any(initial_guess_p):
            pop = toolbox.population(n=(initial_pop_size-1))
            initial_g = creator.Individual(initial_guess)
            pop.insert(0, initial_g)    
        else:
            #specify a population within the bounds that includes the initial guess
            pop = toolbox.population(n=(initial_pop_size-2))
            initial_g = creator.Individual(initial_guess)
            pop.insert(0, initial_g)
            initial_g = creator.Individual(initial_guess_p)
            pop.insert(0, initial_g)   
            
        #Define a function to calculate the fitness of an individual
        cost = partial(self.ga_optimization_function, opt_sp_id_len, assignment, activeCells=activeCells, testSupport=testSupport)
        toolbox.register("evaluate", cost)
        
        #Define the upper and lower bounds for each attribute in the optimization
        lower_bounds = []
        upper_bounds = []
        for i,bnds in enumerate(bounds):
            if bnds[0] is not None:
                lower_bounds.append(float(bnds[0]))
            else:
                warnings.warn("Bounds not specified on an activeCell")
                if initial_guess_p is not None:
                    lower_bounds.append(initial_guess_p[i]-100*initial_guess_p[i])
                else:
                    raise ValueError("Initial guess can't be turned into bounds")
            if bnds[1] is not None:
                upper_bounds.append(float(bnds[1]))
            else:
                warnings.warn("Bounds not specified on an activeCell")
                if initial_guess_p is not None:
                    upper_bounds.append(initial_guess_p[i]+100*initial_guess_p[i])
                else:
                    raise ValueError("Initial guess can't be turned into bounds")
                
                
         #Old Bounds Definition       
#        lower_bounds = [float(bnds[0]) for bnds in bounds]
#        upper_bounds = [float(bnds[1]) for bnds in bounds]
        
        #Define the function to do the mating between two individuals in the previous population
        #Note: Any toolbox.register that are commented out are other possibilities for the function
        
#        
        
        #Note: eta =Crowding degree of the crossover. A high eta will produce children resembling to their parents, while a small eta will produce solutions much more different
        #indpb = the probability of each attribute to be mutated
        #toolbox.register("mate", tools.cxSimulatedBinaryBounded, low=lower_bounds, up=upper_bounds, eta=20.0)
        #toolbox.register("mate", tools.cxSimulatedBinary, eta=10)
        #toolbox.register("mate", tools.cxUniform, indpb = .3)
        toolbox.register("mate", arithXover)
        
        #Define a function to do the mutation for an individual in the next generation
        toolbox.register("mutate", mutUniformFloat, low=lower_bounds, up=upper_bounds, indpb=mutation_rate)
        #toolbox.register("mutate", tools.mutUniformInt, low=lower_bounds, up=upper_bounds, indpb=mutation_rate) #Bounds must be integers to utilize
        #toolbox.register("mutate", tools.mutPolynomialBounded, low=lower_bounds, up=upper_bounds, eta=20.0, indpb=1.0/30.0)
        
        #Define a function to do the selection of individuals 
        #toolbox.register("select", tools.selNSGA2)
        #toolbox.register("select", tools.selBest)
        #toolbox.register("select", tools.selTournament, tournsize=15)
        #toolbox.register("select", tools.selAutomaticEpsilonLexicase)
        toolbox.register("select", selnormGeom, prob_sel_best=0.08)
        
        #Set up the hallof fame in order to report the best individual
        hof = tools.HallOfFame(num_ele_Hallfame, similar=np.array_equal)
    
        
        #Get statistics for each generation in the genetic algoritm
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", np.mean)
        stats.register("std", np.std)
        stats.register("min", np.min)
        stats.register("max", np.max)
        
        #Define an algorithm to structure the iteration through generations
        algorithms.eaMuPlusLambda(pop, toolbox, mu = (3*initial_pop_size), lambda_=(3*initial_pop_size), cxpb=0.5, mutpb=0.5, ngen=num_generations, stats=stats,
                        halloffame=hof) #mu and lambda are only the same because of the normGeomSelect

        globalsection = self.deserializeAssignment_ga(hof[0], opt_sp_id_len)
        return globalsection
            


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
                sheafCells.append(SheafCell(c.dimension,[],c.compactClosure,stalkDim=ls.cobetti(0)))

        return Sheaf(sheafCells),SheafMorphism(mor)

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
    def __init__(self,support,value,source=None):
        """Specify support cell indices and values in each cell stalk for a local section"""
        self.support=support
        self.value=value
        if source is None:
            self.source=support
        else:
            self.source=source

class Section:
    def __init__(self,sectionCells):
        self.sectionCells=sectionCells

    def support(self):
        """List the cells in the support of this section"""
        return {sc.support for sc in self.sectionCells}

    def extend(self,sheaf,cell,value=None,tol=1e-5):
        """Extend the section to another cell; returns True if successful"""
        
        # If the desired cell is already in the support, do nothing
        for sc in self.sectionCells:
            if cell == sc.support:
                if (value is None) or sheaf.cells[sc.support].metric(sc.value,value) < tol:
                    return True
                else:
                    return False

        # Is the desired cell a coface of a cell in the support?
        for s in self.sectionCells:
            for cf in sheaf.cells[s.support].cofaces:
                if cf.index == cell:
                    # If so, extend via restriction
                    val=cf.restriction(s.value)

                    # Check for consistency
                    if value is not None and sheaf.cells[cf.index].metric(val, value)>tol:
                    #if value is not None and np.any(np.abs(val - value)>tol):
                        return False
                    value = val
                    break
            if value is not None:
                break

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

