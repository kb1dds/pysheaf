# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 10:32:46 2018

@author: jhenrich
"""

# Unit test for consistencyRadius
#
# Copyright (c) 2017, Michael Robinson, Janelle Henrich
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied

import pysheaf as ps
import numpy as np
from math import pi


# Metrics and spherical geometry

def distance(v1,v2):
    """Compute Euclidean distance between two pairs"""
    return np.linalg.norm(v2-v1)

def distance_alt(v1,v2):
    """Compute distance in km between two lat/lon/alt triples"""
    x1=v1[0]
    y1=v1[1]
    z1=v1[2]
    x2=v2[0]
    y2=v2[1]
    z2=v2[2]
    return np.sqrt(((x2-x1)*np.pi/180*6356.752*np.cos((y1+y2)*np.pi/360))**2+((y2-y1)*np.pi/180*6356.752)**2 + ((z1-z2)/1000.)**2)


# Sheaf construction (full topology)
#
# Cell: name = stalk (variables)
#
# 0: A = R^3 x R^2 x R (x,y,z,vx,vy,t) (Field Office)
# 1: B = R^3 (x,y,z) (Flight Plan)
# 2: C = R^3 x R^2 (x,y,z,vx,vy) (ATC)
# 3: D = S^1 x R (theta1, t) (RDF 1)
def construct_sheaf():
    global s1
    s1=ps.Sheaf([ps.SheafCell(dimension=0,
                              compactClosure=True,
                              stalkDim=6,
                              metric=distance, # Ignores velocity
                              cofaces=[ps.SheafCoface(index=1, 
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0,0],
                                                                                              [0,1,0,0,0,0],
                                                                                              [0,0,1,0,0,0],
                                                                                              [0,0,0,1,0,0],
                                                                                              [0,0,0,0,1,0]]))), # A->B
                                       ps.SheafCoface(index=2,
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0,0],
                                                                                              [0,1,0,0,0,0],
                                                                                              [0,0,1,0,0,0],
                                                                                              [0,0,0,2,0,0],
                                                                                              [0,0,0,0,1,0]]))), # A->C
                                       ],
                              ), 
                 ps.SheafCell(dimension=1,
                              compactClosure=True,
                              stalkDim=5,
                              metric=distance,
                              cofaces=[ps.SheafCoface(index=3, 
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0],
                                                                                              [0,1,0,0,0],
                                                                                              [0,0,1,0,0],
                                                                                              [0,0,0,2,0]]))), # B->D
                                       ]),
                 ps.SheafCell(dimension=1,
                              compactClosure=True,
                              stalkDim=5,
                              metric=distance, # Ignores velocity
                              cofaces=[ps.SheafCoface(index=3,
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0],
                                                                                              [0,1,0,0,0],
                                                                                              [0,0,1,0,0],
                                                                                              [0,0,0,1,0]]))), # C->D
                                      ]),
                 ps.SheafCell(dimension=2,
                              compactClosure=True,
                              stalkDim=4,
                              metric=distance,
                              cofaces=[]),
                 ])

# 0: A = R^3 x R^2 x R (x,y,z,vx,vy,t) (Field Office)
# 1: B = R^3 (x,y,z) (Flight Plan)
# 2: C = R^3 x R^2 (x,y,z,vx,vy) (ATC)
# 3: D = S^1 x R (theta1, t) (RDF 1)
def construct_inconsistent_sheaf():
    global s2
    s2=ps.Sheaf([ps.SheafCell(dimension=0,
                              compactClosure=True,
                              stalkDim=6,
                              metric=distance, # Ignores velocity
                              cofaces=[ps.SheafCoface(index=1, 
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0,0],
                                                                                              [0,1,0,0,0,0],
                                                                                              [0,0,1,0,0,0],
                                                                                              [0,0,0,1,0,0],
                                                                                              [0,0,0,0,1,0]]))), # A->B
                                       ps.SheafCoface(index=2,
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0,0],
                                                                                              [0,1,0,0,0,0],
                                                                                              [0,0,1,0,0,0],
                                                                                              [0,0,0,1,0,0],
                                                                                              [0,0,0,0,1,0]]))), # A->C
                                       ],
                              ), 
                 ps.SheafCell(dimension=1,
                              compactClosure=True,
                              stalkDim=5,
                              metric=distance,
                              cofaces=[ps.SheafCoface(index=3, 
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0],
                                                                                              [0,1,0,0,0],
                                                                                              [0,0,1,0,0],
                                                                                              [0,0,0,2,0]]))), # B->D
                                       ]),
                 ps.SheafCell(dimension=1,
                              compactClosure=True,
                              stalkDim=5,
                              metric=distance, # Ignores velocity
                              cofaces=[ps.SheafCoface(index=3,
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0],
                                                                                              [0,1,0,0,0],
                                                                                              [0,0,1,0,0],
                                                                                              [0,0,0,1,0]]))), # C->D
                                      ]),
                 ps.SheafCell(dimension=2,
                              compactClosure=True,
                              stalkDim=4,
                              metric=distance,
                              cofaces=[]),
                 ])



# 0: A = R^6 
# 1: B = R^5 
# 2: C = R^5 
# 3: D = R^4 
# 4: E = R^4
def construct_3_step_sheaf():
    global s3
    s3=ps.Sheaf([ps.SheafCell(dimension=0,
                              compactClosure=True,
                              stalkDim=6,
                              metric=distance, # Ignores velocity
                              cofaces=[ps.SheafCoface(index=1, 
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0,0],
                                                                                              [0,1,0,0,0,0],
                                                                                              [0,0,1,0,0,0],
                                                                                              [0,0,0,1,0,0],
                                                                                              [0,0,0,0,1,0]]))), # A->B
                                       ps.SheafCoface(index=2,
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0,0],
                                                                                              [0,1,0,0,0,0],
                                                                                              [0,0,1,0,0,0],
                                                                                              [0,0,0,2,0,0],
                                                                                              [0,0,0,0,1,0]]))), # A->C
                                       ],
                              ), 
                 ps.SheafCell(dimension=1,
                              compactClosure=True,
                              stalkDim=5,
                              metric=distance,
                              cofaces=[ps.SheafCoface(index=3, 
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0],
                                                                                              [0,1,0,0,0],
                                                                                              [0,0,1,0,0],
                                                                                              [0,0,0,2,0]]))), # B->D
                                       ]),
                 ps.SheafCell(dimension=1,
                              compactClosure=True,
                              stalkDim=5,
                              metric=distance, # Ignores velocity
                              cofaces=[ps.SheafCoface(index=3,
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0],
                                                                                              [0,1,0,0,0],
                                                                                              [0,0,1,0,0],
                                                                                              [0,0,0,1,0]]))), # C->D
                                      ]),
                 ps.SheafCell(dimension=2,
                              compactClosure=True,
                              stalkDim=4,
                              metric=distance,
                              cofaces=[ps.SheafCoface(index=4,
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0],
                                                                                              [0,1,0,0],
                                                                                              [0,0,1,0],
                                                                                              [0,0,0,1]])))]),
                 ps.SheafCell(dimension=3,
                              compactClosure=True,
                              stalkDim=4,
                              metric=distance,
                              cofaces=[]),
                 ])


def construct_3_step_sheaf():
    global s3
    s3=ps.Sheaf([ps.SheafCell(dimension=0,
                              compactClosure=True,
                              stalkDim=6,
                              metric=distance, # Ignores velocity
                              cofaces=[ps.SheafCoface(index=1, 
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0,0],
                                                                                              [0,1,0,0,0,0],
                                                                                              [0,0,1,0,0,0],
                                                                                              [0,0,0,1,0,0],
                                                                                              [0,0,0,0,1,0]]))), # A->B
                                       ps.SheafCoface(index=2,
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0,0],
                                                                                              [0,1,0,0,0,0],
                                                                                              [0,0,1,0,0,0],
                                                                                              [0,0,0,2,0,0],
                                                                                              [0,0,0,0,1,0]]))), # A->C
                                       ],
                              ), 
                 ps.SheafCell(dimension=1,
                              compactClosure=True,
                              stalkDim=5,
                              metric=distance,
                              cofaces=[ps.SheafCoface(index=3, 
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0],
                                                                                              [0,1,0,0,0],
                                                                                              [0,0,1,0,0],
                                                                                              [0,0,0,2,0]]))), # B->D
                                       ]),
                 ps.SheafCell(dimension=1,
                              compactClosure=True,
                              stalkDim=5,
                              metric=distance, # Ignores velocity
                              cofaces=[ps.SheafCoface(index=3,
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0,0],
                                                                                              [0,1,0,0,0],
                                                                                              [0,0,1,0,0],
                                                                                              [0,0,0,1,0]]))), # C->D
                                      ]),
                 ps.SheafCell(dimension=2,
                              compactClosure=True,
                              stalkDim=4,
                              metric=distance,
                              cofaces=[ps.SheafCoface(index=4,
                                                      orientation=1,
                                                      restriction=ps.LinearMorphism(np.array([[1,0,0,0],
                                                                                              [0,1,0,0],
                                                                                              [0,0,1,0],
                                                                                              [0,0,0,1]])))]),
                 ps.SheafCell(dimension=3,
                              compactClosure=True,
                              stalkDim=4,
                              metric=distance,
                              cofaces=[]),
                 ])


# Structure the data into various partially-filled Sections,
                 
#Define Input Data
input_data=[ps.Section([ps.SectionCell(support=0,value=np.array([1,1,1,1,1,1]).T)]), # A (Verify that the sheaf's maps are commutative)
            
            ps.Section([ps.SectionCell(support=0,value=np.array([1,1,1,1,1,1]).T), # A
                        ps.SectionCell(support=1,value=np.array([1,1,1,1,3]).T), # B (expected error of 2 for s1)
                        ps.SectionCell(support=2,value=np.array([1,1,1,2,2]).T), # C (expected error of 1 for s1)
                        ps.SectionCell(support=3,value=np.array([4,1,1,2]).T), # D (expected error of 3 for s1)
                        ]), 
            ps.Section([ps.SectionCell(support=0,value=np.array([1,1,1,1,1,1]).T), # A
                        ps.SectionCell(support=1,value=np.array([1,1,1,1,1]).T), # B (expected error of 0 for s2)
                        #ps.SectionCell(support=2,value=np.array([1,1,1,2,1]).T), # C (expected error of 1 for s2)
                        #ps.SectionCell(support=3,value=np.array([4,1,1,2]).T), # D (expected error of sqrt(10) for s2)
                        ]),
             ps.Section([ps.SectionCell(support=0,value=np.array([1,1,1,1,1,1]).T), # A
                        #ps.SectionCell(support=1,value=np.array([1,1,1,1,1]).T), # B (expected error of 0 for s2)
                        #ps.SectionCell(support=2,value=np.array([1,1,1,2,1]).T), # C (expected error of 1 for s2)
                        ps.SectionCell(support=3,value=np.array([4,1,1,2]).T), # D (expected error of 2 for s2)
                        ]),
            ] 
#Construct Consistent Sheaf and check its consistencyRadius
construct_sheaf()

# Exhibit the consistency radius of the partially-filled Section with the input data
consistency_radii=[(s1.isSheaf(case)) for case in input_data]
print 'Raw consistency radii for each test case: ' + str(consistency_radii)

# Demonstrate the consistency radius improves when faulty sensor (U5) is removed
print 'Case 1 consistency radius  ' + str(s1.consistencyRadius(input_data[0],testSupport=[0,1,2,3])) #(expected 0)
print 'Case 2 consistency radius  ' + str(s1.consistencyRadius(input_data[1],testSupport=[0,1,2,3])) #(expected 3)
print 'Case 4 consistency radius  ' + str(s1.consistencyRadius(input_data[3])) #(expected 3....)

#Construct Inconsistent Sheaf and check its consistencyRadius
construct_inconsistent_sheaf()

# Exhibit the consistency radius of the partially-filled Section with the input data
consistency_radii_s2=[(s2.isSheaf(case)) for case in input_data]
print 'Raw consistency radii for each test case: ' + str(consistency_radii_s2)

# Demonstrate the consistency radius improves when faulty sensor (U5) is removed
print 'Case 1 consistency radius  ' + str(s2.consistencyRadius(input_data[0],testSupport=[0,1,2,3])) #(expected 1)
print 'Case 2 consistency radius  ' + str(s2.consistencyRadius(input_data[1],testSupport=[0,1,2,3])) #(expected 3.16....)
print 'Case 3 consistency radius  ' + str(s2.consistencyRadius(input_data[2])) #(expected 1....)


#Show the issue of the 3 steps not being computed

#Redefine Input Data
input_data=[ps.Section([ps.SectionCell(support=0,value=np.array([1,1,1,1,1,1]).T)]), # A (Verify that the sheaf's maps are commutative)
            
            ps.Section([ps.SectionCell(support=0,value=np.array([1,1,1,1,1,1]).T), # A
                        ps.SectionCell(support=1,value=np.array([1,1,1,1,3]).T), # B (expected error of 2 for s1)
                        ps.SectionCell(support=2,value=np.array([1,1,1,2,2]).T), # C (expected error of 1 for s1)
                        ps.SectionCell(support=3,value=np.array([4,1,1,2]).T), # D (expected error of 3 for s1)
                        ]), 
            ps.Section([ps.SectionCell(support=0,value=np.array([1,1,1,1,1,1]).T), # A
                        ps.SectionCell(support=1,value=np.array([1,1,1,1,1]).T), # B (expected error of 0 for s2)
                        #ps.SectionCell(support=2,value=np.array([1,1,1,2,1]).T), # C (expected error of 1 for s2)
                        #ps.SectionCell(support=3,value=np.array([4,1,1,2]).T), # D (expected error of sqrt(10) for s2)
                        ]),
             ps.Section([ps.SectionCell(support=0,value=np.array([1,1,1,1,1,1]).T), # A
                        #ps.SectionCell(support=1,value=np.array([1,1,1,1,1]).T), # B (expected error of 0 for s2)
                        #ps.SectionCell(support=2,value=np.array([1,1,1,2,1]).T), # C (expected error of 1 for s2)
                        ps.SectionCell(support=4,value=np.array([4,1,1,2]).T), # D (expected error of 2 for s2)
                        ]),
            ]
             
#Construct sheaf
construct_3_step_sheaf()

# Exhibit the consistency radius of the partially-filled Section with the input data
consistency_radii=[(s3.isSheaf(case)) for case in input_data]
print 'Raw consistency radii for each test case: ' + str(consistency_radii)

# Demonstrate the consistency radius calculation
print 'Case 1 consistency radius  ' + str(s3.consistencyRadius(input_data[0],testSupport=[0,1,2,3])) #(expected 0)
print 'Case 2 consistency radius  ' + str(s3.consistencyRadius(input_data[1],testSupport=[0,1,2,3])) #(expected 3)
print 'Case 4 consistency radius  ' + str(s3.consistencyRadius(input_data[3])) #(expected 3....)