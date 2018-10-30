import numpy as np


# dataTools is meant to serve as a place to put general serialize/deserialize/compare/edgeMethods that may be
# useful for other users. 


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