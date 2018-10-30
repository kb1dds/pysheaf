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

def Pow2(inputNumber):
   return inputNumber**2 # Pow2
def Pow3(inputNumber):
   return inputNumber**3 # Pow3
def CompareScalars(leftValue,rightValue):
   return abs(leftValue - rightValue) # CompareScalars
def SerializeScalars(assignmentValue):
   return np.array([assignmentValue]) # SerializeScalars



if __name__ == '__main__':
   print("+-+-+-+-+-+-+-+-+-+-+")
   print("+-Unit Test Pysheaf-+")
   print("+-+-+-+-+-+-+-+-+-+-+")
   graph = ps.Sheaf()
   graph.mPreventRedundantExtendedAssignments = False

   TEST_TYPE = "test_type"
   graph.AddCell(0,ps.Cell(TEST_TYPE,CompareScalars,serializeAssignmentMethod= SerializeScalars))
   graph.AddCell(1,ps.Cell(TEST_TYPE,CompareScalars,serializeAssignmentMethod= SerializeScalars))
   graph.AddCell(2,ps.Cell(TEST_TYPE,CompareScalars,serializeAssignmentMethod= SerializeScalars))
   graph.AddCell(3,ps.Cell(TEST_TYPE,CompareScalars,serializeAssignmentMethod= SerializeScalars))
   graph.AddCell(4,ps.Cell(TEST_TYPE,CompareScalars,serializeAssignmentMethod= SerializeScalars))
   graph.AddCell(5,ps.Cell(TEST_TYPE,CompareScalars,serializeAssignmentMethod= SerializeScalars))

   assert(graph.number_of_nodes() == 6), "Incorrect number of cells"

   """
   Test Structure

         0
         |
         1
         |
         2
        / \
       4   3
        \ /
         5
   """

   coface_square = ps.Coface(TEST_TYPE,TEST_TYPE,Pow2)
   graph.AddCoface(0,1,coface_square)
   coface_square = ps.Coface(TEST_TYPE,TEST_TYPE,Pow2)
   graph.AddCoface(1,2,coface_square)
   coface_square = ps.Coface(TEST_TYPE,TEST_TYPE,Pow2)
   graph.AddCoface(2,3,coface_square)
   coface_square = ps.Coface(TEST_TYPE,TEST_TYPE,Pow2)
   graph.AddCoface(2,4,coface_square)
   coface_square = ps.Coface(TEST_TYPE,TEST_TYPE,Pow2)
   graph.AddCoface(3,5,coface_square)
   coface_square = ps.Coface(TEST_TYPE,TEST_TYPE,Pow3)
   graph.AddCoface(4,5,coface_square)

   assert(graph.number_of_edges() == 6), "Incorrect number of cofaces"

   print("number of nodes in graph",graph.number_of_nodes())
   print("number of edges in graph",graph.number_of_edges())

   known_test_values = [2,4,16,256,65536,16777216]
   graph.GetCell(0).SetDataAssignment(ps.Assignment(TEST_TYPE,2))

   print("dataAssignment: ",graph.GetCell(0).mDataAssignment)

   graph.MaximallyExtendCell(0)


   assert(graph.GetCell(0).GetDataAssignment().mValue == known_test_values[0]), "Incorrect dataAssignment"
   assert(graph.GetCell(1).GetExtendedAssignmentValueList() == [ps.Assignment(TEST_TYPE,known_test_values[1])]), "Incorrect ExtendedAssignment"
   assert(graph.GetCell(2).GetExtendedAssignmentValueList() == [ps.Assignment(TEST_TYPE,known_test_values[2])]), "Incorrect ExtendedAssignment"
   assert(graph.GetCell(3).GetExtendedAssignmentValueList() == [ps.Assignment(TEST_TYPE,known_test_values[3])]), "Incorrect ExtendedAssignment"
   assert(graph.GetCell(4).GetExtendedAssignmentValueList() == [ps.Assignment(TEST_TYPE,known_test_values[3])]), "Incorrect ExtendedAssignment"
   assert(graph.GetCell(5).GetExtendedAssignmentValueList() == [ps.Assignment(TEST_TYPE,known_test_values[4]),ps.Assignment(TEST_TYPE,known_test_values[5])]), "Incorrect ExtendedAssignment"


   cell2_test_value = 6
   cell2_test_answer = abs(known_test_values[2] - cell2_test_value)
   cell5_test_value = 12345
   cell5_test_answer = abs(known_test_values[5] - cell5_test_value)

   graph.GetCell(2).SetDataAssignment(ps.Assignment(TEST_TYPE,cell2_test_value))
   graph.GetCell(5).SetDataAssignment(ps.Assignment(TEST_TYPE,cell5_test_value))

   print("Single cell consistency:",graph.GetCell(2).ComputeConsistency(np.inf))

   assert(graph.GetCell(2).ComputeConsistency(np.inf) == cell2_test_answer), "Incorrect consistency"
   assert(graph.GetCell(5).ComputeConsistency(np.inf) == cell5_test_answer), "Incorrect consistency"

   print("consistency radius:",graph.ComputeConsistencyRadius())

   print("nodes in graph:",graph.nodes())

   print("vales for 5: ",graph.GetCell(5).GetExtendedAssignmentValueList())

   graph.GetCell(0).mOptimizationCell = True
   graph.mPreventRedundantExtendedAssignments = True

   cell2_test_value = 81
   cell5_test_value = 43046721
   result_answer = 3

   graph.GetCell(2).SetDataAssignment(ps.Assignment(TEST_TYPE,cell2_test_value))
   graph.GetCell(5).SetDataAssignment(ps.Assignment(TEST_TYPE,cell5_test_value))

   graph.ClearExtendedAssignments()
   graph.mNumpyNormType = None
   graph.MaximallyExtendCell(0)
   print("consistency radius before fuse assignment:",graph.ComputeConsistencyRadius())



   fuse_results = graph.FuseAssignment()
   assert(abs(result_answer - fuse_results.x[0]) <.1), "Incorrect fuse assignment result"
   print("consistency radius after fuse assignment:",graph.ComputeConsistencyRadius())


   print("+-+-+-+-+-+-+-+-+-+-+")
   print("+-All tests passed!-+")
   print("+-+-+-+-+-+-+-+-+-+-+")