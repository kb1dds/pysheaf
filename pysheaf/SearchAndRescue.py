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

# Unit test for data fusion process using sheaf from
# M. Robinson, "Sheaves are the canonical data structure for information integration," Info. Fusion (36), 2017.


import pysheaf as ps
import numpy as np
import copy
# Metrics and spherical geometry
def Distance(v1,v2):
   """Compute distance in km between two lat/lon pairs"""
   x1=v1[0]
   y1=v1[1]
   x2=v2[0]
   y2=v2[1]
   return np.sqrt(((x2-x1)*np.pi/180*6356.752*np.cos((y1+y2)*np.pi/360))**2+((y2-y1)*np.pi/180*6356.752)**2)

def DistanceTime(x,y):
   return np.linalg.norm(x-y)

def DistanceLatLongAlt(v1,v2):
   """Compute distance in km between two lat/lon/alt triples"""
   x1=v1[0]
   y1=v1[1]
   z1=v1[2]
   x2=v2[0]
   y2=v2[1]
   z2=v2[2]
   return np.sqrt(((x2-x1)*np.pi/180*6356.752*np.cos((y1+y2)*np.pi/360))**2+((y2-y1)*np.pi/180*6356.752)**2 + ((z1-z2)/1000.)**2)

def Anglemetric(theta1,theta2):
   """Distance between two angles"""
   return min([abs(theta1-theta2),abs(theta1-theta2+360),abs(theta1-theta2-360)])

def MinimumOfBearingAndTime(x,y):
   return min(Anglemetric(x[0],y[0]),abs(x[1]-y[1]))

def Bearing(tx_lon,tx_lat,rx_lon,rx_lat):
   return np.arctan2((tx_lon-rx_lon)*np.cos(rx_lat*np.pi/180),tx_lat-rx_lat)*180/np.pi

# Functions from the paper, generally including additional spherical earth corrections
def A(vec):
   """The map A from the paper, but with a spherical earth correction"""
   r1x=-73.662574
   r1y=42.7338328
   vec2=E(vec) # Predict crash location by dead reckoning
   return np.array([Bearing(vec2[0],vec2[1],r1x,r1y),vec[5]])

def B(vec):
   """The map B from the paper, but with a spherical earth correction"""
   r2x=-77.0897374
   r2y=38.9352387
   vec2=E(vec) # Predict crash location by dead reckoning
   return np.array([Bearing(vec2[0],vec2[1],r2x,r2y),vec[5]])

def C(vec):
   """The map C from the paper, but with a spherical earth correction"""
   r1x=-73.662574
   r1y=42.7338328
   return Bearing(vec[0],vec[1],r1x,r1y)

def D(vec):
   """The map D from the paper, but with a spherical earth correction"""
   r2x=-77.0897374
   r2y=38.9352387
   return Bearing(vec[0],vec[1],r2x,r2y)

def E(vec):
   """The map E from the paper, but with a spherical earth correction, and accounting for units.
   Note: this map predicts crash location lon,lat using dead reckoning"""
   lon=vec[0]
   lat=vec[1]
   vlon=vec[3]
   vlat=vec[4]
   dt=vec[5]
   return np.array([lon+dt*vlon/6356.752*180/np.pi/np.cos(lat*np.pi/180),lat+dt*vlat/6356.752*180/np.pi])

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


if __name__ == '__main__':

   print("+-+-+-+-+-+-+-+-+-+-+")
   print("+-Search and rescue-+")
   print("+-+-+-+-+-+-+-+-+-+-+")

   sheaf = ps.Sheaf()

   POSITION = "lat_lon_degrees"
   POSITION_3D_VELOCITY = "lat_lon_degrees_alt_meters_velocity_kilometers_per_hour_type"
   POSITION_3D_VELOCITY_TIME = "lat_lon_degrees_alt_meters_velocity_kilometers_per_hour_time_hours_type"
   POSITION_3D = "lat_lon_alt_degrees_meters_type"
   TIME_HOURS = "time_hours"
   BEARING_DEGREES = "bearing_degrees"
   BEARING_DEGREES_AND_TIME_HOURS = "bearing_degrees_and_time_hours"



   sheaf.AddCell("X",ps.Cell(POSITION_3D_VELOCITY_TIME,DistanceLatLongAlt,dataDimension=6))
   sheaf.AddCell("U1",ps.Cell(POSITION_3D,DistanceLatLongAlt,dataDimension=3))
   sheaf.AddCell("U2",ps.Cell(POSITION_3D_VELOCITY,DistanceLatLongAlt,dataDimension=5))

   sheaf.AddCell("U3",ps.Cell(BEARING_DEGREES_AND_TIME_HOURS,MinimumOfBearingAndTime,dataDimension=2))
   sheaf.AddCell("U4",ps.Cell(BEARING_DEGREES_AND_TIME_HOURS,MinimumOfBearingAndTime,dataDimension=2))
   sheaf.AddCell("U5",ps.Cell(POSITION,Distance,dataDimension=2))
   sheaf.AddCell("V1",ps.Cell(BEARING_DEGREES,Anglemetric))
   sheaf.AddCell("V2",ps.Cell(BEARING_DEGREES,Anglemetric))
   sheaf.AddCell("V3",ps.Cell(TIME_HOURS,DistanceTime))

   """
   Search and rescue Structure

            X 
         / / \ \
       U2 U3 U5 U4
       /   /\/ \/      
      /   / /\ /\    
     U1   V1  V3 V2

   """

   sheaf.AddCoface("X","U2",ps.Coface(POSITION_3D_VELOCITY_TIME,POSITION_3D_VELOCITY,LinearMorphism(np.array([[1,0,0,0,0,0],
                                                                                          [0,1,0,0,0,0],
                                                                                          [0,0,1,0,0,0],
                                                                                          [0,0,0,1,0,0],
                                                                                          [0,0,0,0,1,0]]))))
   sheaf.AddCoface("X","U3",ps.Coface(POSITION_3D_VELOCITY_TIME,BEARING_DEGREES_AND_TIME_HOURS,A))
   sheaf.AddCoface("X","U4",ps.Coface(POSITION_3D_VELOCITY_TIME,BEARING_DEGREES_AND_TIME_HOURS,B))
   sheaf.AddCoface("X","U5",ps.Coface(POSITION_3D_VELOCITY_TIME,POSITION,E))
   sheaf.AddCoface("U2","U1",ps.Coface(POSITION_3D_VELOCITY,POSITION_3D,LinearMorphism(np.array([[1,0,0,0,0],
                                                                                                [0,1,0,0,0],
                                                                                                [0,0,1,0,0]]))))
   sheaf.AddCoface("U3","V1",ps.Coface(BEARING_DEGREES_AND_TIME_HOURS,BEARING_DEGREES,LinearMorphism(np.array([[1,0]]))))
   sheaf.AddCoface("U3","V3",ps.Coface(BEARING_DEGREES_AND_TIME_HOURS,TIME_HOURS,LinearMorphism(np.array([[0,1]]))))
   sheaf.AddCoface("U4","V3",ps.Coface(BEARING_DEGREES_AND_TIME_HOURS,TIME_HOURS,LinearMorphism(np.array([[0,1]]))))
   sheaf.AddCoface("U4","V2",ps.Coface(BEARING_DEGREES_AND_TIME_HOURS,BEARING_DEGREES,LinearMorphism(np.array([[1,0]]))))
   sheaf.AddCoface("U5","V1",ps.Coface(POSITION,BEARING_DEGREES,C))
   sheaf.AddCoface("U5","V2",ps.Coface(POSITION,BEARING_DEGREES,D))

   sheaf.GetCell("X").SetDataAssignment(ps.Assignment(POSITION_3D_VELOCITY_TIME,np.array([-70.649,42.753,11220,495,164,0.928])))
   sheaf.GetCell("U1").SetDataAssignment(ps.Assignment(POSITION_3D,np.array([-70.662,42.829,11178])))
   sheaf.GetCell("U2").SetDataAssignment(ps.Assignment(POSITION_3D_VELOCITY,np.array([-70.587,42.741,11346,495,164])))
   sheaf.GetCell("U3").SetDataAssignment(ps.Assignment(BEARING_DEGREES_AND_TIME_HOURS,np.array([77.1,0.943])))
   sheaf.GetCell("U4").SetDataAssignment(ps.Assignment(BEARING_DEGREES_AND_TIME_HOURS,np.array([61.3,0.890])))
   sheaf.GetCell("U5").SetDataAssignment(ps.Assignment(POSITION,np.array([-64.599,44.243])))

  

   sheaf.GetCell("X").mOptimizationCell = True
   sheaf.GetCell("U1").mOptimizationCell = True
   sheaf.GetCell("U2").mOptimizationCell = True
   sheaf.GetCell("U3").mOptimizationCell = True
   sheaf.GetCell("U4").mOptimizationCell = True
   sheaf.GetCell("U5").mOptimizationCell = True

   # sheaf.mNumpyNormType = None
   sheaf.MaximallyExtendCell("X")
   print("consistency radius before gobal fuse assignment:",sheaf.ComputeConsistencyRadius())
   sheaf.ClearExtendedAssignments()

   print(sheaf.FuseAssignment())
   print("consistency radius after global fuse assignment:",sheaf.ComputeConsistencyRadius())


   second_data_assignments =[
   (ps.Assignment(POSITION_3D_VELOCITY_TIME,np.array([-70.668,42.809,11431,495,164,1.05]))),
   (ps.Assignment(POSITION_3D,np.array([-70.663,42.752,11299]))),
   (ps.Assignment(POSITION_3D_VELOCITY,np.array([-70.657,42.773,11346,495,164]))),
   (ps.Assignment(BEARING_DEGREES_AND_TIME_HOURS,np.array([77.2,0.930]))),
   (ps.Assignment(BEARING_DEGREES_AND_TIME_HOURS,np.array([63.2,0.974]))),
   (ps.Assignment(POSITION,np.array([-64.630,44.287])))]

   tmp_input_assignments = copy.deepcopy(second_data_assignments)
   

   sheaf.GetCell("X").SetDataAssignment(tmp_input_assignments[0])
   sheaf.GetCell("U1").SetDataAssignment(tmp_input_assignments[1])
   sheaf.GetCell("U2").SetDataAssignment(tmp_input_assignments[2])
   sheaf.GetCell("U3").SetDataAssignment(tmp_input_assignments[3])
   sheaf.GetCell("U4").SetDataAssignment(tmp_input_assignments[4])
   sheaf.GetCell("U5").SetDataAssignment(tmp_input_assignments[5])

   sheaf.ClearExtendedAssignments()
   sheaf.MaximallyExtendCell("X")
   print("2nd data set: consistency radius before gobal fuse assignment:",sheaf.ComputeConsistencyRadius())

   sheaf.FuseAssignment()
   print("2nd data set: consistency radius after global fuse assignment:",sheaf.ComputeConsistencyRadius())

   third_data_assignments =[
   (ps.Assignment(POSITION_3D_VELOCITY_TIME,np.array([-70.626,42.814,11239,419,311,1.02]))),
   (ps.Assignment(POSITION_3D,np.array([-70.612,42.834,11237]))),
   (ps.Assignment(POSITION_3D_VELOCITY,np.array([-70.617,42.834,11236,419,310]))),
   (ps.Assignment(BEARING_DEGREES_AND_TIME_HOURS,np.array([77.2,0.985]))),
   (ps.Assignment(BEARING_DEGREES_AND_TIME_HOURS,np.array([63.3,1.05]))),
   (ps.Assignment(POSITION,np.array([-62.742,44.550])))]

   tmp_input_assignments = copy.deepcopy(third_data_assignments)
   sheaf.GetCell("X").SetDataAssignment(tmp_input_assignments[0])
   sheaf.GetCell("U1").SetDataAssignment(tmp_input_assignments[1])
   sheaf.GetCell("U2").SetDataAssignment(tmp_input_assignments[2])
   sheaf.GetCell("U3").SetDataAssignment(tmp_input_assignments[3])
   sheaf.GetCell("U4").SetDataAssignment(tmp_input_assignments[4])
   sheaf.GetCell("U5").SetDataAssignment(tmp_input_assignments[5])

   sheaf.ClearExtendedAssignments()
   sheaf.MaximallyExtendCell("X")
   print("3rd data set: consistency radius before gobal fuse assignment:",sheaf.ComputeConsistencyRadius())

   sheaf.FuseAssignment()
   print("3rd data set: consistency radius after global fuse assignment:",sheaf.ComputeConsistencyRadius())


   sheaf.GetCell("X").mOptimizationCell = False
   sheaf.GetCell("U1").mOptimizationCell = False
   sheaf.GetCell("U2").mOptimizationCell = False
   sheaf.GetCell("U3").mOptimizationCell = False
   sheaf.GetCell("U4").mOptimizationCell = False
   sheaf.GetCell("U5").mOptimizationCell = True

   sheaf.GetCell("X").SetDataAssignment(ps.Assignment(POSITION_3D_VELOCITY_TIME,np.array([-70.649,42.753,11220,495,164,0.928])))
   sheaf.GetCell("U1").SetDataAssignment(ps.Assignment(POSITION_3D,np.array([-70.662,42.829,11178])))
   sheaf.GetCell("U2").SetDataAssignment(ps.Assignment(POSITION_3D_VELOCITY,np.array([-70.587,42.741,11346,495,164])))
   sheaf.GetCell("U3").SetDataAssignment(ps.Assignment(BEARING_DEGREES_AND_TIME_HOURS,np.array([77.1,0.943])))
   sheaf.GetCell("U4").SetDataAssignment(ps.Assignment(BEARING_DEGREES_AND_TIME_HOURS,np.array([61.3,0.890])))
   sheaf.GetCell("U5").SetDataAssignment(ps.Assignment(POSITION,np.array([-64.599,44.243])))

   sheaf.ClearExtendedAssignments()
   sheaf.MaximallyExtendCell("X")
   print("consistency radius before local fuse assignment:",sheaf.ComputeConsistencyRadius())
   sheaf.FuseAssignment()
   print("consistency radius after local fuse assignment:",sheaf.ComputeConsistencyRadius())

   tmp_input_assignments = copy.deepcopy(second_data_assignments)
   sheaf.GetCell("X").SetDataAssignment(tmp_input_assignments[0])
   sheaf.GetCell("U1").SetDataAssignment(tmp_input_assignments[1])
   sheaf.GetCell("U2").SetDataAssignment(tmp_input_assignments[2])
   sheaf.GetCell("U3").SetDataAssignment(tmp_input_assignments[3])
   sheaf.GetCell("U4").SetDataAssignment(tmp_input_assignments[4])
   sheaf.GetCell("U5").SetDataAssignment(tmp_input_assignments[5])

   sheaf.ClearExtendedAssignments()
   sheaf.MaximallyExtendCell("X")
   print("2nd data set: consistency radius before local fuse assignment:",sheaf.ComputeConsistencyRadius())
   sheaf.FuseAssignment()
   print("2nd data set: consistency radius after local fuse assignment:",sheaf.ComputeConsistencyRadius())

   tmp_input_assignments = copy.deepcopy(second_data_assignments)
   sheaf.GetCell("X").SetDataAssignment(tmp_input_assignments[0])
   sheaf.GetCell("U1").SetDataAssignment(tmp_input_assignments[1])
   sheaf.GetCell("U2").SetDataAssignment(tmp_input_assignments[2])
   sheaf.GetCell("U3").SetDataAssignment(tmp_input_assignments[3])
   sheaf.GetCell("U4").SetDataAssignment(tmp_input_assignments[4])
   sheaf.GetCell("U5").SetDataAssignment(tmp_input_assignments[5])

   sheaf.ClearExtendedAssignments()
   sheaf.MaximallyExtendCell("X")
   print("3rd data set: consistency radius before local fuse assignment:",sheaf.ComputeConsistencyRadius())

   sheaf.FuseAssignment()
   print("3rd data set: consistency radius after local fuse assignment:",sheaf.ComputeConsistencyRadius())