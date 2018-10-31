# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 16:16:22 2018

@author: mrobinson
"""

import pysheaf as ps

shf=ps.Sheaf()

shf.AddCell('A',ps.Cell('Scalar'))
shf.AddCell('B',ps.Cell('Scalar'))
shf.AddCell('C',ps.Cell('Scalar'))

shf.AddCoface('A','C',ps.Coface('Scalar','Scalar',lambda x: x))
shf.AddCoface('B','C',ps.Coface('Scalar','Scalar',lambda x: x))
shf.GetCell('A').SetDataAssignment(ps.Assignment('Scalar',0))
shf.GetCell('B').SetDataAssignment(ps.Assignment('Scalar',1))
shf.GetCell('C').SetDataAssignment(ps.Assignment('Scalar',0.25))

shf.mPreventRedundantExtendedAssignments=False

shf.MaximallyExtendCell('A')
shf.MaximallyExtendCell('B')
shf.MaximallyExtendCell('C')

print('Consistency radius : {}'.format(shf.ComputeConsistencyRadius()))

for thres in [0., 0.2, 0.3, 0.8]:
    print('Consistent stars at {} : {}'.format(thres,shf.ConsistentStarCollection(thres)))
