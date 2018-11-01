""" 
debugg Python 3.6 Sheaf manipulation library
"""
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



import networkx as nx
import numpy as np
import scipy.optimize

# Assignment Types
CONST_VALUE_TYPE_SCALAR = "Scalar" #for ints, doubles, etc. one value. 

class Cell:
   """
   The cell class is the building block of a sheaf. Also a node in the graph. 
   The cell has a data assignment which serves as a \"measurement\", test or known value.
   Extended assignments are a \"guess\" by the sheaf at what the data value should be.
   
   The cell has a compare method which is used to compare two of the same types of assignment objects. 

   """
   def __init__(self,dataTagType,compareAssignmentsMethod=None, serializeAssignmentMethod=None, deSerializeAssignmentMethod=None, dataDimension=1, optimizationCell = False, extendFromThisCell=True):
      self.mDataDimension = dataDimension
      self.mDataAssignmentPresent = False
      self.mDataAssignment = 0
      self.mExtendedAssignments = {}
      self.mDataTagType = dataTagType
      self.mOptimizationCell = optimizationCell # if the cell's data assignment can be changed by the optimizer
      self.mExtendFromThisCell = extendFromThisCell # if the cell should be maximally extended from
      self.mBounds = [(None,None)] * self.mDataDimension # python magic for a list of tuples the length of dimension
      if compareAssignmentsMethod==None:
         self.Compare = self.DefaultCompareAssignments
      else:
         self.Compare = compareAssignmentsMethod

      if serializeAssignmentMethod==None:
         self.mSerialize = self.DefaultSerialize
      else:
         self.mSerialize = serializeAssignmentMethod

      if deSerializeAssignmentMethod==None:
         self.mDeserialize = self.DefaultDeserialize
      else:
         self.mDeserialize = deSerializeAssignmentMethod
      return # __init__

   def SetDataAssignment(self, dataAssignment):
      if self.mDataTagType == dataAssignment.mValueType:
         self.mDataAssignmentPresent = True
         self.mDataAssignment = dataAssignment
      else:
         print("Cell::SetDataAssignment: Error, DataAssignment is of Incorrect type name. Expected:",self.mDataTagType,"Actual:",dataAssignment.mValueType)
      return # SetDataAssignment

   def SetCompareMethod(self, methodToSet):
      self.Compare = methodToSet
      return # SetCompareMethod

   def SetSerializeMethod(self, methodToSet):
      """ 
      Sets a Serialize method that will turn a data assignment value into a one dimensional numpy array that is 1 by mDataDimension

      :returns: a numpy array of the serialized values
      """
      self.mSerialize = methodToSet
      return # SetSerializeMethod

   def SetDeserializeMethod(self, methodToSet):
      """
      methodToSet takes the place of DefaultDeserialize.
      :param methodToSet: Takes the place of DefaultDeserialize. Must populate the sheaf with the deserialized values 
         and return the same serialized argument list 
      :returns: void
      """
      self.mDeserialize = methodToSet
      return # SetDeserializeMethod

   def SetBounds(self, boundsToSet):
      """
      Sets optimization bounds on the serialized values given to the optimizer
      :param boundsToSet: A list of tuples specifying the lower and upper bound for each independent variable [(xl0, xu0),(xl1, xu1),...]
         Null means no bound for that element
      :returns: void
      """
      if len(boundsToSet) != self.mDataDimension:
         print("Cell::SetBounds: Error, bounds list length must match data dimension")
         
      self.mBounds = boundsToSet
      return #SetBounds

   def GetDataAssignment(self):
      if self.mDataAssignmentPresent == False:
         print("Cell::GetDataAssignment: Error, DataAssignment not present for cell")
      return self.mDataAssignment # GetDataAssignment

   def AddExtendedAssignment(self, cellPathTuple,extendedAssignment):
      self.mExtendedAssignments[cellPathTuple] = extendedAssignment
      return # AddExtendedAssignment

   def GetExtendedAssignmentValueList(self,cellStartIndices=None):
      if cellStartIndices is None:
         return list(self.mExtendedAssignments.values())
      else:
         return [self.mExtendedAssignments[key] for key in self.mExtendedAssignments.keys() 
                 if key[0] in cellStartIndices] # GetExtendedAssignment

   def CheckExtendedAssignmentPresent(self):
      return bool(self.mExtendedAssignments) # CheckExtendedAssignmentPresent

   def CheckDataAssignmentPresent(self):
      return self.mDataAssignmentPresent # CheckDataAssignmentPresent

   def AbleToComputeConsistency(self):
      return self.mDataAssignmentPresent and (len(self.mExtendedAssignments) != 0) # AbleToComputeConsistency

   def ComputeConsistency(self, numpyNormType=np.inf, cellStartIndices=None): 
      """
      Computes the error between the data assignments and extended assignments on the sheaf
      :param numpyNormType: Optional, how the errors across the sheaf are combined. default is take maximum error.
      :param cellStartIndices: Optional, which cells may start the data flowed into this cell
      :returns: Error between the data assignments and extended assignments on the sheaf
      """
      if self.mDataAssignmentPresent == False:
         print("Cell::ComputeConsistency: Error, DataAssignment not present for cell")
      if len(self.mExtendedAssignments) == 0:
         print("Cell::ComputeConsistency: Error, ExtendedAssignment not present for cell")
      tmp_assignments = self.GetExtendedAssignmentValueList(cellStartIndices)
      if len(tmp_assignments) == 0:
         return 0.
      assignment_comparisions = []
      for assignment in tmp_assignments:
         assignment_comparisions.append(self.Compare(self.mDataAssignment.mValue,assignment.mValue))
      return np.linalg.norm(assignment_comparisions,ord=numpyNormType) # ComputeConsistency

   def DefaultCompareAssignments(self,leftValue,rightValue):
      return abs(leftValue - rightValue) #DefaultCompareAssignments

   def ClearExtendedAssignment(self):
      self.mExtendedAssignments.clear()
      return # ClearExtendedAssignment

   def SerializeAssignment(self):
      if self.mDataAssignmentPresent == False:
         print("Cell::SerializeAssignment: Error, DataAssignment not present for cell")
      return self.mSerialize(self.mDataAssignment.mValue)# SerializeAssignment

   def DeserializeAssignment(self,serializedValueArray):
      return self.mDeserialize(serializedValueArray) # DeserializeAssignment

   def DefaultSerialize(self,assignmentValue):
      """ 
      Must return a 1-D ndarray of float.

      :param assignmentValue: The assignment value to serialize 
      :returns: Whatever is passed in by assignmentValue
      """
      return assignmentValue # DefaultSerialize

   def DefaultDeserialize(self, serializedValueArray):
      """ 
      Uses elements in serializedValueArray to populate data assignments on the sheaf.
      Uses the number equal to mDataDimension. 

      :param serializedValueArray: the serialized values for the entire sheaf
      :returns: the array given that has the elements used removed.
      """
      if self.mDataAssignmentPresent == False:
         print("Cell::DefaultDeserialize: Error, DataAssignment not present for cell")

      if self.mDataDimension > serializedValueArray.size:
         print("Cell::DefaultDeserialize: Error, not enough values in array to deserialize")

      self.mDataAssignment.mValue = serializedValueArray[:self.mDataDimension]
      return serializedValueArray[self.mDataDimension:]# DefaultDeserialize

class Coface:
   """
   This class lives on the edges between cells. Assignments are passed through
      edge methods to be placed on different cells. 
   """
   def __init__(self, inputTagType,outputTagType,edgeMethod):
      """
      edgeMethod  is called when the sheaf extends from one cell to another.
      The edge method return value will set the value on a assignment
      """
      self.mOrientation = 0 # orientation is either -1 or 1.
      self.mInputTagType = inputTagType
      self.mOutputTagType = outputTagType
      self.mEdgeMethod = edgeMethod
      return # __init__

   def SetEdgeMethod(self,methodToSet):
      """
      Sets a edge method that is called when the sheaf extends from one cell to another.
      The edge method return value will set the value on a assignment
      """
      self.mEdgeMethod = methodToSet
      return # SetEdgeMethod

   def RunEdgeMethod(self,inputAssignment):
      return Assignment(self.mOutputTagType, self.mEdgeMethod(inputAssignment.mValue)) # RunEdgeMethod

class Assignment:
   """
   This class holds the data that is placed on the sheaf. It has a type field that can be checked against
   """
   def __init__(self, valueType,value):
      self.mValue = value
      self.mValueType = valueType 
      return # __init__
   def __str__(self):
      return str(self.mValue) # __str__
   def __repr__(self):
      return str(self.mValue) # __repr__
   def __eq__(self, assignmentToCompare):
     if isinstance(assignmentToCompare, Assignment):
         return self.mValue == assignmentToCompare.mValue
     return False # __eq__


class Sheaf(nx.DiGraph):
   """ 
   This class serves to model sheafs.
   Users will build sheafs with 
   AddCell
   AddCoface

   Then evaluate sheafs with 
   FuseAssignment 
   ComputeConsistencyRadius
   """
   def __init__(self):
      """ 
      users can set the following at any time after construction 
      mNumpyNormType - what norm type the consistency radius will use. Default is np.inf
      mPreventRedundantExtendedAssignments - Gives a speed up by preventing re-extending the same assignment through different paths. Default is on
      mMaximumOptimizationIterations - limits the number of iterations the optimizer can call OptimizationIteration
      :returns: void
      """
      nx.DiGraph.__init__(self)
      self.mNumpyNormType = np.inf
      self.mPreventRedundantExtendedAssignments = False # If True, gives a speed up by preventing re-extending the same assignment through different paths.
      self.mMaximumOptimizationIterations = 100
      self.mSheafOptimizer = self.DefaultSheafOptimizer
      return # __init__

   def AddCell(self,cellIndex,cellToAdd):
      self.add_node(cellIndex,vertex = cellToAdd)
      return # AddCell

   def AddCoface(self,cellIndexFrom, cellIndexTo, cofaceToAdd):
      if (self.GetCell(cellIndexFrom).mDataTagType == cofaceToAdd.mInputTagType) and (self.GetCell(cellIndexTo).mDataTagType == cofaceToAdd.mOutputTagType):
         self.add_edge(cellIndexFrom,cellIndexTo,edge=cofaceToAdd)
      else:
         print("Sheaf::AddCoface: Error, coface types do not match cells from:",cellIndexFrom,"to:",cellIndexTo)
      return # AddCoface

   def GetCell(self,cellIndex):
      return self.node[cellIndex]['vertex'] # GetCell

   def GetCoface(self, cellIndexStart, cellIndexTo):
      return self[cellIndexStart][cellIndexTo]['edge'] # GetCoface

   def GetCellIndexList(self):
      """
      :returns: a list of indexes of all cells currently in the sheaf
      """
      return self.nodes() # GetCellIndexList

   def ExtendCellRecursiveTraverse(self,cellIndex,cellPathTuple):
      """
      This method is used by MaximallyExtendCell to populate the sheaf in a efficient manner. 

      :param cellIndex: The cell to start the traversal from. 
      :param cellPathTuple: The indexes of the cells that have been visited during this traversal
      :returns: void
      """
      dict_keyiter = self.successors(cellIndex)

      # For all children of this node
      for successor_index in dict_keyiter:
         if self.mPreventRedundantExtendedAssignments == True:
            if self.GetCell(successor_index).CheckExtendedAssignmentPresent() == True:
               continue

         # Get the path of the new assignment for the child
         next_cell_path_tuple = cellPathTuple+(successor_index,)

         # If this the the starting cell use the value in data
         if len(cellPathTuple) == 1:
            current_cell_assignment = self.GetCell(cellIndex).mDataAssignment
         else:
            current_cell_assignment = self.GetCell(cellIndex).mExtendedAssignments[cellPathTuple]

         # Calculate new assignment for the child
         next_cell_assignment = self.GetCoface(cellIndex,successor_index).RunEdgeMethod(current_cell_assignment)

         # Add the new assignment to the child cell
         self.GetCell(successor_index).AddExtendedAssignment(next_cell_path_tuple,next_cell_assignment)

         # Repeat process on the child cell
         self.ExtendCellRecursiveTraverse(successor_index,next_cell_path_tuple)
      return # ExtenCellRecursiveTraverse

   def MaximallyExtendCell(self,startCellIndex):
      """
      This method with extend a cell's data value as far as it can over the sheaf. 
      This is done depth first. 

      :param startCellIndex: cell to start extending from
      :returns: void
      """
      try:
         nx.find_cycle(self,startCellIndex)
      except:
         self.ExtendCellRecursiveTraverse(startCellIndex,(startCellIndex,))
      else:
         print("Sheaf::MaximallyExtendCell: Error, Cycle found in graph from cell index:", startCellIndex)
      return # MaximallyExtendCell

   def CheckForExtendedAssignment(self,cellPathTuple,cellIndexToCheck):
      return cellPathTuple in self.GetCell(cellIndexToCheck).mExtendedAssignments

   def ClearExtendedAssignments(self):
      cell_index_list = self.nodes()
      for cell_index in cell_index_list:
         self.GetCell(cell_index).ClearExtendedAssignment()
      return # ClearExtendedAssignments

   def ConsistentStarCollection(self,consistencyThreshold):
      """
      :returns: The indexes of cells whose stars have consistency radius less than consistencyThreshold
      """
      cell_index_list = self.nodes()
      cell_index_below_threshold = []
      for cell_index in cell_index_list:
         local_consistency_radius=self.ComputeConsistencyRadius([cell_index]+list(self.successors(cell_index)))
         if local_consistency_radius < consistencyThreshold:
            cell_index_below_threshold.append(cell_index)
            
      # Remove redundant stars in the cover
      maximal_star_list = []
      for cell_index in cell_index_below_threshold:
          not_a_duplicate=True
          for pred in self.predecessors(cell_index):
              if pred in cell_index_below_threshold:
                  not_a_duplicate=False
                  break
          if not_a_duplicate:
              maximal_star_list.append(cell_index)

      return maximal_star_list # ConsistentStarCollection
       
   def CellIndexesLessThanConsistencyThreshold(self,consistencyThreshold):
      """
      :returns: The indexes of cells that have consistencies less than consistencyThreshold
      """
      cell_index_list = self.nodes()
      cell_index_below_threshold = []
      for cell_index in cell_index_list:
         if self.GetCell(cell_index).AbleToComputeConsistency() == True:
            if self.GetCell(cell_index).ComputeConsistency(self.mNumpyNormType) < consistencyThreshold:
               cell_index_below_threshold.append(cell_index)
      return cell_index_below_threshold # CellIndexesLessThanConsistencyThreshold

   def ComputeConsistencyRadius(self,cell_indexes=None):
      """
      This method will call compute consistency on all cells (default) or those specified by cell_indexes in the sheaf. 
      The default behavior is to then return the max error. This can be changed by setting mNumpyNormType

      :returns: consistency radius of the sheaf assignment
      """
      if cell_indexes is None:
         cell_index_list = self.nodes()
      else:
         cell_index_list = cell_indexes

      cell_consistancies_list = []
      for cell_index in cell_index_list:
         if self.GetCell(cell_index).AbleToComputeConsistency() == True:
            if cell_indexes is None:
               cell_consistancies_list.append(self.GetCell(cell_index).ComputeConsistency(self.mNumpyNormType)
            else:
               cell_consistancies_list.append(self.GetCell(cell_index).ComputeConsistency(self.mNumpyNormType, cellStartIndices=cell_indexes))
      if not cell_consistancies_list:
         print("Sheaf::ComputeConsistencyRadius: Error, unable to compute consistency for any cells")
      return np.linalg.norm(cell_consistancies_list,ord=self.mNumpyNormType) # ComputeConsistencyRadius

   def SerializeAssignments(self):
      """
      This method will convert the values of a sheaf into a one dimensional numpy array that is 1 by mDataDimension

      :returns: a numpy array of the serialized values
      """
      optimization_cells_list = self.GetOptimizationCellIndexList()
      cell_value_array = np.array([])

      for cell_index in optimization_cells_list:
         cell_value_array = np.append(cell_value_array,[self.GetCell(cell_index).SerializeAssignment()])

      return cell_value_array # SerializeAssignments

   def GetOptimizationCellIndexList(self):
      """
      Used to help the optimizer know which cells have data values it can change. 
      :returns: a list of indexes of all cells that have mOptimizationCell set to true.
      """
      cell_index_list = self.nodes()
      optimization_cells_list = []

      for cell_index in cell_index_list:
         if self.GetCell(cell_index).mOptimizationCell == True:
            optimization_cells_list.append(cell_index)

      return optimization_cells_list # GetOptimizationCellIndexList

   def GetExtendableCellIndexList(self):
      """
      Used to help MaximallyExtendCell know which cells are valid to start extending from.
      :returns: list of cell indexes that have mExtendFromThisCell set to true and have a data Assignment. 
      """
      cell_index_list = self.nodes()
      extendable_cells_list = []

      for cell_index in cell_index_list:
         if self.GetCell(cell_index).mExtendFromThisCell == True and self.GetCell(cell_index).mDataAssignmentPresent == True:
            extendable_cells_list.append(cell_index)

      return extendable_cells_list # GetExtendableCellIndexList

   def GetBoundsList(self):
      """
      :returns: The a list of all bounds for the optimizer.
      """
      optimization_cells_list = self.GetOptimizationCellIndexList()
      bounds_list = []

      for optimization_cell_index in optimization_cells_list:
         bounds_list.extend(self.GetCell(optimization_cell_index).mBounds)
      return bounds_list # GetBoundsList

   def SetSheafOptimizer(self,sheafOptimizer):
      self.mSheafOptimizer = sheafOptimizer
      return # SetSheafOptimizer


   def DeserializeAssignments(self, valueArray):
      optimization_cells_list = self.GetOptimizationCellIndexList()
      tmp_value_array = valueArray
      for optimization_cell_index in optimization_cells_list:
         tmp_value_array = self.GetCell(optimization_cell_index).DeserializeAssignment(tmp_value_array)

      return # DeserializeAssignments

   def OptimizationIteration(self,optimizationCellDataAssignments):
      """
      This method is meant to be run by the genetic algorithm or optimizer. 

      :returns: The consistency radius of the sheaf
      """
      self.ClearExtendedAssignments()
      self.DeserializeAssignments(optimizationCellDataAssignments)
      extendable_cells_list = self.GetExtendableCellIndexList()
      for extendable_cell_index in extendable_cells_list:
         self.MaximallyExtendCell(extendable_cell_index)
      return self.ComputeConsistencyRadius()# OptimizationIteration

   def FuseAssignment(self):
      """
      The optimizer is set up to run on all optimization cells. 
      :returns: The return value of the employed optimization function
      """
      optimization_cell_data_values = self.SerializeAssignments()
      bounds_list = self.GetBoundsList()
      optimizer_results = self.mSheafOptimizer(self.OptimizationIteration,optimization_cell_data_values,bounds_list,self.mMaximumOptimizationIterations)
      return optimizer_results # FuseAssignment 


   def DefaultSheafOptimizer(self, functionToIterate, serializedAssignments, boundsList,maxIterations):
      return scipy.optimize.minimize(functionToIterate,serializedAssignments,method = "SLSQP",bounds=boundsList, tol = 1e-5, options = {'maxiter':maxIterations}) # DefualtSheafOptimizer
