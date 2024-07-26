import json
import numpy as np
import networkx as nx
import pysheaf as ps
import matplotlib.pyplot as plt

netlist=json.load(open('logic_netlist.json'))

shf = ps.Sheaf()
shf.fromNetlist(netlist)
shf.mNumpyNormType = 2

pos=nx.layout.spring_layout(shf)
nx.draw_networkx_labels(shf,pos)
nx.draw_networkx_edges(shf,pos)
plt.show()

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
            
