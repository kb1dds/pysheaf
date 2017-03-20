# Unit test for incremental cell creation in PySheaf
#
# Copyright (c) 2017, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 

import pysheaf as ps
import numpy as np

# Create some cells
cA=ps.SheafCell(name='A',dimension=0,stalkDim=2)
cB=ps.SheafCell(name='B',dimension=0,stalkDim=2)
cAB=ps.SheafCell(name='AB',dimension=1,stalkDim=2)

# Initialize the sheaf structure
s=ps.Sheaf()

# Add two cells, remembering their IDs assigned by the Sheaf
cA_id=s.add_cell(cA)
cAB_id=s.add_cell(cAB)

# Add a coface
s.add_coface((cA_id,cAB_id),orientation=1,restriction=np.array([[1,0],[0,1]]))

# Add another cell and coface
cB_id=s.add_cell(cB)
s.add_coface((cB_id,cAB_id),orientation=1,restriction=np.array([[0,1],[1,0]]))
