# MIT License

# Copyright (c) 2018 Michael Robinson & Steven Fiacco

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
import pysheaf as ps
#import dataTools
import networkx as nx

def LocallyFuseAssignment( sheaf ):
    '''
    Perform Sheaf.FuseAssignment() locally, which should use less memory (but is probably slower)
    '''
    
    # Any cell marked as not being solved for doesn't get extended assignments
    sheaf.ClearExtendedAssignments()
    
    # Locality for this is based on connected sets of optimization cells.  
    optimization_cell_indices=sheaf.GetOptimizationCellIndexList() 
    
    # To find connected components of optimization cells, build an undirected graph whose vertices are optimization cells
    # and edges are induced by the sheaf's graph
    G=sheaf.to_undirected().subgraph(optimization_cell_indices)
    
    # For each component, 
    #  1. Turn on each component as optimization cells, 
    #  2. turn off the rest, and 
    #  3. run FuseAssignment.  
    for component in nx.connected_components(G):
        # TODO TBD: The views thing is not working correctly.  I think I need to copy!
        for node in sheaf.nodes():
            if node in component:
                sheaf.GetCell(node).mOptimizationCell=True
            else:
                sheaf.GetCell(node).mOptimizationCell=False
                sheaf.GetCell(node).ClearExtendedAssignment()

        # Focus our attention on this component, only.  (This avoids flowdown dependencies from other cells)
        subsheaf=nx.classes.graphviews.subgraph_view(sheaf,filter_node=lambda x: (x in component) or [nd2 for nd2 in component if x in sheaf.successors(nd2)])
        
        # Perform the optimization
        subsheaf.FuseAssignment()

    # Put the optimization cells back to the way they were
    for node in sheaf.nodes():
        if node in optimization_cell_indices:
            sheaf.GetCell(node).mOptimizationCell=True
        else:
            sheaf.GetCell(node).mOptimizationCell=False
            
    # Propagate all extended assignments, since that's what Sheaf.FuseAssignment() does at this point
    for c in sheaf.GetCellIndexList(): 
        sheaf.MaximallyExtendCell(c)
    
    return

def BuildConstantSheaf(G, dataDimension=1):
    """Construct a constant sheaf on a graph G with a given dataDimension"""
    shf=ps.Sheaf()
    
    # Add cells for each node in the graph
    for node in G.nodes():
        shf.AddCell(node, ps.Cell('vector',dataDimension=dataDimension))
        
    # Add cofaces for each edge in the graph
    for edge in G.edges():
        shf.AddCoface(edge[0],edge[1],ps.Coface('vector','vector',dataTools.LinearMorphism(np.eye(dataDimension))))
    
    return shf # BuildConstantSheaf

class LinearAlgebraSheafAnalysisTool:
   def __init__(self, inputSheaf):
      self.mSheaf = inputSheaf
      return # __init__
   def GetEdgeMorphism(self,cellIndexStart, cellIndexTo):
      return self.mSheaf.GetCoface(cellIndexStart,cellIndexTo).mEdgeMethod # GetEdgeMorphism

   def GetAllEdgeMorphisms(self):
      all_edges_on_sheaf = self.mSheaf.edges()
      all_morphisms_list = []
      for edge in all_edges_on_sheaf:
         all_morphisms_list.append(self.GetEdgeMorphism(edge[0],edge[1]))
      return all_morphisms_list
