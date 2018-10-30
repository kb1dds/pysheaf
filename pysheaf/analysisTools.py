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
import dataTools




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