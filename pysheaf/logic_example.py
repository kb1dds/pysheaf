""" 
# Example to demonstrate the use of `Sheaf.fromNetlist()` to build a sheaf specified by a JSON file.

"""

# MIT License

# Copyright (c) 2024 Michael Robinson

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

import json
import numpy as np
import networkx as nx
import pysheaf as ps
#import matplotlib.pyplot as plt

netlist=json.load(open('pysheaf/logic_netlist.json'))

shf = ps.Sheaf()
shf.fromNetlist(netlist)
shf.mNumpyNormType = 2

#pos=nx.layout.spring_layout(shf)
#nx.draw_networkx_labels(shf,pos)
#nx.draw_networkx_edges(shf,pos)
#plt.show()

for a in [0,1]:
    for b in [0,1]:
        for c in [0,1]:
            q = min((1-a*b)+c,1)

            shf.GetCell('A').SetDataAssignment(ps.Assignment('cell',np.array((a,))))
            shf.GetCell('B').SetDataAssignment(ps.Assignment('cell',np.array((b,))))
            shf.GetCell('C').SetDataAssignment(ps.Assignment('cell',np.array((c,))))

            shf.FuseAssignment()

            qp = shf.GetCell('Q').mDataAssignment
            
            print('A={}, B={}, C={},  Q={} =?= {} '.format(a,b,c,q,qp))
            
