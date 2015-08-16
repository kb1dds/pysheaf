# Relative simplicial homology library
#
# Copyright (c) 2014-2015, Michael Robinson
# Distribution of unaltered copies permitted for noncommercial use only
# All other uses require express permission of the author
# This software comes with no warrantees express or implied 

import numpy as np
try:
    import networkx as nx
except:
    pass

def kernel(A, tol=1e-5):
    u, s, vh = np.linalg.svd(A)
    sing=np.zeros(vh.shape[0],dtype=np.complex)
    sing[:s.size]=s
    null_mask = (sing <= tol)
    null_space = np.compress(null_mask, vh, axis=0)
    return null_space.conj().T
    
def cokernel(A, tol=1e-5):
    u, s, vh = np.linalg.svd(A)
    sing=np.zeros(u.shape[1],dtype=np.complex)
    sing[:s.size]=s
    null_mask = (sing <= tol)
    return np.compress(null_mask, u, axis=1)

def toplexify(simplices):
    """Reduce a simplicial complex to merely specification of its toplices"""
    simplices=sorted(simplices,key=len,reverse=True)
    return [spx for i,spx in enumerate(simplices) if not [sp2 for j,sp2 in enumerate(simplices) if (i!=j and set(spx).issubset(sp2)) and (i>j or set(spx)==set(sp2))]]
    
def ksublists(lst,n,sublist=[]):
    """Iterate over all ordered n-sublists of a list lst"""
    if n==0:
        yield sublist
    else:
        for idx in range(len(lst)):
           item=lst[idx]
           for tmp in ksublists(lst[idx+1:],n-1,sublist+[item]):
              yield tmp

def ksimplices(toplexes,k,relative=None):
    """List of k-simplices in a list of toplexes"""
    simplices=[]
    for toplex in toplexes:
        for spx in ksublists(toplex,k+1):
            if not spx in simplices and (relative == None or not spx in relative):
                simplices.append(spx)
    return simplices 

def kflag(subcomplex,k,verts=None):
    """Determine the k-cells of a flag complex"""
    if verts==None:
        verts=list(set([v for s in subcomplex for v in s]))

    return [s for s in ksublists(verts,k+1) if all([[r for r in subcomplex if (p[0] in r) and (p[1] in r)] for p in ksublists(s,2)])]

def flag(subcomplex,maxdim=None):
    """Construct the flag complex based on a given subcomplex"""
    verts=list(set([v for s in subcomplex for v in s]))
    complex=[]
    added=True
    k=0
    while added and ((maxdim == None) or (k<=maxdim)):
        added=kflag(subcomplex,k,verts=verts)
        complex+=added
        k+=1
    return complex
      
def vietorisRips(pts,diameter,maxdim=None):
    """Construct a Vietoris-Rips complex over a point cloud"""
    subcomplex=[[x,y] for x in range(len(pts)) for y in range(x,len(pts)) if np.linalg.norm(pts[x]-pts[y])<diameter]
    if maxdim == 1:
        return subcomplex
    if maxdim == None:
        G=nx.Graph()
        G.add_edges_from(subcomplex)
        return nx.find_cliques(G)
    else:
        return flag(subcomplex,maxdim)

def star(complex,face):
    """Compute the star over a given face.  Works with toplexes (it will return toplexes only, then) or full complexes."""
    return [s for s in complex if set(face).issubset(s)] 

def boundary(toplexes,k,relative=None):
    """Compute the k-boundary matrix for a simplicial complex"""
    kchain=ksimplices(toplexes,k,relative)
    
    if k <= 0:
        return np.zeros((0,len(kchain)))
    
    km1chain=ksimplices(toplexes,k-1,relative)
    bndry=np.zeros((len(km1chain),len(kchain)))
    for i,km1spx in enumerate(km1chain):
        for j,kspx in enumerate(kchain):
           parity=1
           for boundary in ksublists(kspx,k):
               parity*=-1
               if set(km1spx)==set(boundary):
                   bndry[i,j]+=parity
    return bndry    

def simplicialHomology(X,k,Y=None,tol=1e-5):
    """Compute relative k-homology of a simplicial complex"""
    dk=boundary(X,k,Y)
    dkp1=boundary(X,k+1,Y)
    return homology(dk,dkp1)

def homology(b1,b2,tol=1e-5):
    """Compute homology from two matrices whose product b1*b2=0"""
    b2p=np.compress(np.any(abs(b2)>tol,axis=0),b2,axis=1)
        
    # Compute kernel
    if b1.size:
        ker=kernel(b1,tol);
    else:
        ker=np.eye(b1.shape[1])

    # Remove image
    if b2.any():
        map,j1,j2,j3=np.linalg.lstsq(ker,b2p)
        Hk=np.dot(ker,cokernel(map,tol));
    else:
        Hk=ker

    return Hk
    
def localHomology(toplexes,k,facet):
    """Compute local homology relative to a facet"""
    rel=[spx for spx in (ksimplices(toplexes,k)+ksimplices(toplexes,k-1)) 
         if not [f for f in facet if f in spx]]
    return simplicialHomology(toplexes,k,rel)

def cone(toplexes,subcomplex,coneVertex='*'):
    """Construct the cone over a subcomplex.  The cone vertex can be renamed if desired.  The resulting complex is homotopy equivalent to the quotient by the subcomplex."""
    return toplexes + [spx+[coneVertex] for spx in subcomplex]

def integerVertices(complex):
    """Rename vertices in a complex so that they are all integers"""
    vertslist=list(set([v for s in complex for v in s]))
    return [map(lambda x: vertslist.index(x),s) for s in complex]

def vertexHoplength(toplexes,vertices,maxHoplength=None):
    """Compute the edge hoplength distance to a set of vertices within the complex.  Returns the list of tuples (vertex,distance)"""
    vertslist=list(set([v for s in toplexes for v in s if v not in vertices]))
    outlist=[(v,0) for v in vertices] # Prebuild output of the function
    outsimplices=[w for s in toplexes if [v for v in vertices if v in s] for w in s]

    # Consume the remaining vertices
    currentlen=1
    while vertslist and (maxHoplength==None or currentlen < maxHoplength):
        # Find all vertices adjacent to the current ones
        nextverts=[v for v in vertslist if v in outsimplices]

        # Rebuild the output
        outlist=outlist+[(v,currentlen) for v in nextverts]

        # Recompute implicated simplices
        outsimplices=[w for s in toplexes if [v for v,jnk in outlist if v in s] for w in s]

        # Rebuild current vertex list
        vertslist=[v for v in vertslist if v not in nextverts]
        currentlen+=1

    return outlist
    
def vertexLabels2simplices(toplexes,vertexLabels,nonConeLabel=None,coneName=None):
    """Propagate numerical vertex labels (list of pairs: (vertex, label)) to labels on simplices.  If you want the same label applied to simplices not touching a particular 'cone' vertex named coneName, set that nonConeLabel to the desired label.  Note: will use the smallest vertex label or infinity if none given"""
    if nonConeLabel == None:
        return [(s,min([label for v,label in vertexLabels if v in s]+[float("inf")])) for s in toplexes]
    else:
        return [(s,min([label for v,label in vertexLabels if v in s]+[float("inf")]) if [v for v in s if v==coneName] else nonConeLabel) for s in toplexes]

def iteratedCone(toplexes,subcomplex):
    """Prepare to compute local homology centered at subcomplex using Perseus.  In order to that, this function constructs a labeled sequence of (coned off) simplices paired with birth times. NB: We want reduced homology, but Perseus doesn't compute that..."""

    # Hoplengths to all vertices -- will be used for anti-birth times
    vertexLabels=vertexHoplength(toplexes,list(set([v for s in subcomplex for v in s])))
    
    # Maximum hoplength
    maxHopLength=max([t for v,t in vertexLabels])
    
    # All the toplexes in the final, coned off complex
    coned=cone(toplexes,[s for s in toplexes if s not in subcomplex])

    # Propagate vertex labels
    conedLabeled=vertexLabels2simplices(coned,vertexLabels,nonConeLabel=maxHopLength-1,coneName='*')

    # Renumber vertices and return
    vertslist=list(set([v for s in coned for v in s]))
    return [(map(lambda x: vertslist.index(x),s),maxHopLength-t) for s,t in conedLabeled]

def complex2perseus(filename,toplexes,labels=False):
    """Write a complex out as a file for input into Perseus.  If there are labels containing birth times, set labels=True"""
    with open(filename,'wt') as fp:
        fp.write('1\n')
        for tp in toplexes:
            fp.write(str(len(tp[0])-1)+' ')
            for v in tp[0]:
                fp.write(str(v+1)+' ')

            if labels:
                fp.write(str(tp[1])+'\n')
            else:
                fp.write('1 \n')
