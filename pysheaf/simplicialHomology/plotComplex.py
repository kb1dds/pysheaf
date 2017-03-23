# Unit test for persistence sheaf library
#
# Copyright (c) 2016-2017, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 

import numpy as np
import matplotlib.pyplot as plt
import simplicialHomology as sh

def plot_complex(locations,cplx,color=None):
    """Plot an abstract simplicial complex given locations of nodes"""
    colors='mbcgyr'

    # Plot simplices of dimension 2 and higher
    for idx,f in enumerate(cplx):
        if len(f)>2:
            if color != None:
                if color[idx] >= len(colors):
                    col=colors[-1]
                else:
                    col=colors[color[idx]]
            else:
                col='g'
            for face in sh.ksublists(f,3):
                x,y=zip(*[locations[i] for i in face])
                plt.fill(x,y,col,alpha=0.5)

    # Plot simplices of dimension 1 (edges)
    for idx,f in enumerate(cplx):
        if len(f)==2:
            if color != None:
                if color[idx] >= len(colors):
                    col=colors[-1]
                else:
                    col=colors[color[idx]]
            else:
                col='k'
            x,y=zip(*[locations[i] for i in f])
            plt.plot(x,y,color=col,lw=3)

    # Plot simplices of dimension 0 (vertices)
    for idx,f in enumerate(cplx):
        if len(f)==1:
            if color != None:
                if color[idx] >= len(colors):
                    col=colors[-1]
                else:
                    col=colors[color[idx]]
            else:
                col='k'
            plt.plot(locations[f[0]][0],locations[f[0]][1],color=col,ms=4*np.pi,marker='.',markersize=30)
    return
