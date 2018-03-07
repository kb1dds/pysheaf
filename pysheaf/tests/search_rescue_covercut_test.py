# Unit test for data fusion process using sheaf from
# M. Robinson, "Sheaves are the canonical data structure for information integration," Info. Fusion (36), 2017.
#
# Copyright (c) 2017, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied

import pysheaf as ps
import numpy as np

import pysheaf.covers as covers

# Metrics and spherical geometry
def distance(v1,v2):
    """Compute distance in km between two lat/lon pairs"""
    x1=v1[0]
    y1=v1[1]
    x2=v2[0]
    y2=v2[1]
    return np.sqrt(((x2-x1)*np.pi/180*6356.752*np.cos((y1+y2)*np.pi/360))**2+((y2-y1)*np.pi/180*6356.752)**2)

def distance_alt(v1,v2):
    """Compute distance in km between two lat/lon/alt triples"""
    x1=v1[0]
    y1=v1[1]
    z1=v1[2]
    x2=v2[0]
    y2=v2[1]
    z2=v2[2]
    return np.sqrt(((x2-x1)*np.pi/180*6356.752*np.cos((y1+y2)*np.pi/360))**2+((y2-y1)*np.pi/180*6356.752)**2 + ((z1-z2)/1000.)**2)

def anglemetric(theta1,theta2):
    """Distance between two angles"""
    return min([abs(theta1-theta2),abs(theta1-theta2+360),abs(theta1-theta2-360)])

def bearing(tx_lon,tx_lat,rx_lon,rx_lat):
    return np.arctan2((tx_lon-rx_lon)*np.cos(rx_lat*np.pi/180),tx_lat-rx_lat)*180/np.pi

# Functions from the paper, generally including additional spherical earth corrections
def A(vec):
    """The map A from the paper, but with a spherical earth correction"""
    r1x=-73.662574
    r1y=42.7338328
    vec2=E(vec) # Predict crash location by dead reckoning
    return np.array([bearing(vec2[0],vec2[1],r1x,r1y),vec[5]])

def B(vec):
    """The map B from the paper, but with a spherical earth correction"""
    r2x=-77.0897374
    r2y=38.9352387
    vec2=E(vec) # Predict crash location by dead reckoning
    return np.array([bearing(vec2[0],vec2[1],r2x,r2y),vec[5]])

def C(vec):
    """The map C from the paper, but with a spherical earth correction"""
    r1x=-73.662574
    r1y=42.7338328
    return bearing(vec[0],vec[1],r1x,r1y)

def D(vec):
    """The map D from the paper, but with a spherical earth correction"""
    r2x=-77.0897374
    r2y=38.9352387
    return bearing(vec[0],vec[1],r2x,r2y)

def E(vec):
    """The map E from the paper, but with a spherical earth correction, and accounting for units.
    Note: this map predicts crash location lon,lat using dead reckoning"""
    lon=vec[0]
    lat=vec[1]
    vlon=vec[3]
    vlat=vec[4]
    dt=vec[5]
    return np.array([lon+dt*vlon/6356.752*180/np.pi/np.cos(lat*np.pi/180),lat+dt*vlat/6356.752*180/np.pi])

# Sheaf construction (full topology)
#
# Cell: name = stalk (variables)
#
# 0: X = R^3 x R^2 x R (x,y,z,vx,vy,t)
# 1: U1 = R^3 (x,y,z)
# 2: U2 = R^3 x R^2 (x,y,z,vx,vy)
# 3: U3 = S^1 x R (theta1, t)
# 4: U4 = S^1 x R (theta2, t)
# 5: U5 = R^2 (x,y)
# 6: V1 = S1 (theta1)
# 7: V2 = S1 (theta2)
# 8: V3 = R (t)

s1=ps.Sheaf([ps.SheafCell(dimension=0,
                          compactClosure=True,
                          stalkDim=6,
                          metric=distance_alt, # Ignores velocity
                          cofaces=[ps.SheafCoface(index=2, 
                                                  orientation=1,
                                                  restriction=ps.LinearMorphism(np.array([[1,0,0,0,0,0],
                                                                                          [0,1,0,0,0,0],
                                                                                          [0,0,1,0,0,0],
                                                                                          [0,0,0,1,0,0],
                                                                                          [0,0,0,0,1,0]]))), # X->U2
                                   ps.SheafCoface(index=3,
                                                  orientation=1,
                                                  restriction=ps.SetMorphism(A)), # X->U3
                                   ps.SheafCoface(index=4,
                                                  orientation=1,
                                                  restriction=ps.SetMorphism(B)), # X->U4
                                   ps.SheafCoface(index=5,
                                                  orientation=1,
                                                  restriction=ps.SetMorphism(E))]), # X->U5 
             ps.SheafCell(dimension=2,
                          compactClosure=True,
                          stalkDim=3,
                          metric=distance_alt,
                          cofaces=[]), # U1
             ps.SheafCell(dimension=1,
                          compactClosure=True,
                          stalkDim=5,
                          metric=distance_alt, # Ignores velocity
                          cofaces=[ps.SheafCoface(index=1,
                                                  orientation=1,
                                                  restriction=ps.LinearMorphism(np.array([[1,0,0,0,0],
                                                                                          [0,1,0,0,0],
                                                                                          [0,0,1,0,0]])))]), # U2->U1
             ps.SheafCell(dimension=1,
                          compactClosure=True,
                          stalkDim=2,
                          metric=lambda x,y: min(anglemetric(x[0],y[0]),abs(x[1]-y[1])),
                          cofaces=[ps.SheafCoface(index=6,
                                                  orientation=1,
                                                  restriction=ps.LinearMorphism(np.array([[1,0]]))), # U3->V1
                                   ps.SheafCoface(index=8,
                                                  orientation=1,
                                                  restriction=ps.LinearMorphism(np.array([[0,1]])))]), # U3->V3
             ps.SheafCell(dimension=1,
                          compactClosure=True,
                          stalkDim=2,
                          metric=lambda x,y: min(anglemetric(x[0],y[0]),abs(x[1]-y[1])),
                          cofaces=[ps.SheafCoface(index=8,
                                                  orientation=-1,
                                                  restriction=ps.LinearMorphism(np.array([[0,1]]))), # U4->V3
                                   ps.SheafCoface(index=7,
                                                  orientation=1,
                                                  restriction=ps.LinearMorphism(np.array([[1,0]])))]), # U4->V2
             ps.SheafCell(dimension=1,
                          compactClosure=True,
                          stalkDim=2,
                          metric=distance,
                          cofaces=[ps.SheafCoface(index=6,
                                                  orientation=-1,
                                                  restriction=ps.SetMorphism(C)), # U5->V1
                                   ps.SheafCoface(index=7,
                                                  orientation=-1,
                                                  restriction=ps.SetMorphism(D))]), # U5->V2
             ps.SheafCell(dimension=2,
                          compactClosure=True,
                          stalkDim=1,
                          metric=anglemetric,
                          cofaces=[]), # V1
             ps.SheafCell(dimension=2,
                          compactClosure=True,
                          stalkDim=1,
                          metric=anglemetric,
                          cofaces=[]), # V2
             ps.SheafCell(dimension=2,
                          compactClosure=True,
                          stalkDim=1,
                          cofaces=[])]) # V3

# Structure the data from Table 1 in the paper into various partially-filled Sections, one for each Case listed
input_data=[ps.Section([ps.SectionCell(support=0,value=np.array([-70.649,42.753,11220,495,164,0.928])), # X
                        ps.SectionCell(support=1,value=np.array([-70.662,42.829,11178])), # U1
                        ps.SectionCell(support=2,value=np.array([-70.587,42.741,11346,495,164])), # U2
                        ps.SectionCell(support=3,value=np.array([77.1,0.943])), # U3
                        ps.SectionCell(support=4,value=np.array([61.3,0.890])), # U4
                        ps.SectionCell(support=5,value=np.array([-64.599,44.243]))]), # U5
            
            ps.Section([ps.SectionCell(support=0,value=np.array([-70.668,42.809,11431,495,164,1.05])), # X
                        ps.SectionCell(support=1,value=np.array([-70.663,42.752,11299])), # U1
                        ps.SectionCell(support=2,value=np.array([-70.657,42.773,11346,495,164])), # U2
                        ps.SectionCell(support=3,value=np.array([77.2,0.930])), # U3
                        ps.SectionCell(support=4,value=np.array([63.2,0.974])), # U4
                        ps.SectionCell(support=5,value=np.array([-64.630,44.287]))]), # U5
            
            ps.Section([ps.SectionCell(support=0,value=np.array([-70.626,42.814,11239,419,311,1.02])), # X
                        ps.SectionCell(support=1,value=np.array([-70.612,42.834,11237])), # U1
                        ps.SectionCell(support=2,value=np.array([-70.617,42.834,11236,419,310])), # U2
                        ps.SectionCell(support=3,value=np.array([77.2,0.985])), # U3
                        ps.SectionCell(support=4,value=np.array([63.3,1.05])), # U4
                        ps.SectionCell(support=5,value=np.array([-62.742,44.550]))]), # U5

            ps.Section([ps.SectionCell(support=1,value=np.array([-70.612,42.834,11237])), # U1
                        ps.SectionCell(support=2,value=np.array([-70.617,42.834,11236,419,310])), # U2
                        ps.SectionCell(support=3,value=np.array([77.2,0.985])), # U3
                        ps.SectionCell(support=4,value=np.array([63.3,1.05])), # U4
                        ps.SectionCell(support=5,value=np.array([-62.742,44.550]))])] 


# Exhibit the consistency radius of the partially-filled Section with the input data
consistency_radii=[s1.consistencyRadius(case) for case in input_data]
print 'Raw consistency radii for each test case: ' + str(consistency_radii)

# Demonstrate the consistency radius improves when faulty sensor (U5) is removed
print 'Case 1 consistency radius after removing faulty sensor ' + str(s1.consistencyRadius(input_data[0],testSupport=[0,1,2,3,4,6,7,8]))
print 'Case 2 consistency radius after removing faulty sensor ' + str(s1.consistencyRadius(input_data[1],testSupport=[0,1,2,3,4,6,7,8]))
print 'Case 3 consistency radius after removing faulty sensor ' + str(s1.consistencyRadius(input_data[2],testSupport=[0,1,2,3,4,6,7,8]))

# Iterate over all possible covers rooted on U2, U3, U4, U5
bestCov=None
bestFOM=np.inf
for roots in covers.partitions_iter([i for i,c in enumerate(s1.cells) if c.dimension==1]):
    cov=[s1.starCells(op) for op in roots]
    print str(cov)
    FOM = s1.coverFigureofMerit(input_data[2],cov)
    if FOM < bestFOM:
        print 'Improved FOM found: ' + str(FOM) + '; ' + str(cov)
        bestCov=cov
        bestFOM=FOM

print 'Best FOM found: ' + str(bestFOM) + '; ' + str(bestCov)
