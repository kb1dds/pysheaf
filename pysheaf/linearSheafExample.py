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
import pysheaf.dataTools
import pysheaf.analysisTools


def CompareArrays(array1,array2):
   return np.average(array1-array2) #CompareArrays

if __name__ == '__main__':
   print("+-+-+-+-+-+-+-+-+-+-+-+-+")
   print("+-Linear Sheaf Example!-+")
   print("+-+-+-+-+-+-+-+-+-+-+-+-+")

   graph = ps.Sheaf()

   LINEAR_TYPE = "linear_morphism_type"
   graph.AddCell(0,ps.Cell(LINEAR_TYPE,CompareArrays,dataDimension=4))
   graph.AddCell(1,ps.Cell(LINEAR_TYPE,CompareArrays,dataDimension=4))
   graph.AddCell(2,ps.Cell(LINEAR_TYPE,CompareArrays,dataDimension=2))


   graph.GetCell(0).SetDataAssignment(ps.Assignment(LINEAR_TYPE,np.array([5,6,7,8])))
   graph.GetCell(0).SetBounds([(0,None),(0,None),(0,None),(0,None)])
   graph.GetCell(0).mOptimizationCell = True

   graph.GetCell(1).SetDataAssignment(ps.Assignment(LINEAR_TYPE,np.array([15,16,17,18])))
   graph.GetCell(2).SetDataAssignment(ps.Assignment(LINEAR_TYPE,np.array([30,32])))

   graph.AddCoface(0,1,ps.Coface(LINEAR_TYPE,LINEAR_TYPE,dataTools.LinearMorphism(np.array([[1,0,0,0],
                                                                                 [0,1,0,0],
                                                                                 [0,0,1,0],
                                                                                 [0,0,0,1]]))))
   graph.AddCoface(1,2,ps.Coface(LINEAR_TYPE,LINEAR_TYPE,dataTools.LinearMorphism(np.array([[2,0,0,0],
                                                                                 [0,2,0,0]]))))
   graph.MaximallyExtendCell(0)
   print("consistency radius before gobal fuse assignment:",graph.ComputeConsistencyRadius())
   graph.ClearExtendedAssignments()

   print(graph.FuseAssignment())
   print("consistency radius after global fuse assignment:",graph.ComputeConsistencyRadius())



   analysisDemo= analysisTools.LinearAlgebraSheafAnalysisTool(graph);

   all_morphisms = analysisDemo.GetAllEdgeMorphisms()
