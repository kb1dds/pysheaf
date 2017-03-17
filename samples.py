'''
Created on Mar 7, 2017
@author: brendapraggastis
'''
import pysheaf as ps
import plotComplex as pltx
import simplicialHomology as sh
import itertools as it
import pytest as pt
import pdb


def mkcell(*argv):
    '''
        Used for constructing test examples.
        Create a string named cell based on the character assigned to the
        nodes listed in argv
    '''
    temp = ""
    for nd in argv:
        temp += nd
    return temp
#
###### Examples: Orientation on the cell complex is implicitly defined
######           by the ordering of the vertices
v,t,r,g = ['V','T','R','G']
exWF = [v,t,r,g,mkcell(v,t),mkcell(v,r),mkcell(t,r),mkcell(r,g),mkcell(v,t,r)]

a,c,e,k,t,u,v = ['A','C','E','K','T','U','V']
ex7x6 = [a,c,e,k,t,u,v,mkcell(a,c),mkcell(a,e),mkcell(a,k),mkcell(a,u),
           mkcell(a,v),mkcell(c,e),mkcell(c,k),mkcell(k,u),mkcell(k,v),
           mkcell(u,v),mkcell(t,u),mkcell(a,c,k),mkcell(a,k,u),mkcell(a,k,v),mkcell(a,u,v),mkcell(k,u,v),mkcell(a,k,u,v)]


# example = ex7x6

example = exWF

exidx = dict( zip(example,range(len(example))) )

if __name__ == '__main__':
###### generate the cells, their cofaces, and the induced cellComplex
    exCells = []
    for cell in example:
        temp = ps.Cell(len(cell)-1,name=cell)
        exCells.append(temp)
    for idx, cell in enumerate(exCells): ## cell represents the coface
        for l in range(1,len(cell.name)):
            temp = list(it.combinations(cell.name,l))
            for cbn in temp:
                face = ("").join(cbn)
                ort = ps.simplexOrientation(face,cell.name)
                id = exidx[face]
#                 print face, "id: ",id, "ort: ", ort
#                 print "line 46 ",cells[id].cofaceList();
                exCells[id].cofaces.append(ps.Coface(exidx[cell.name],ort))

#                 print id,cells[id].name,cells[id].cofaceList()
#             print cells[id].name, cells[id].cofaceList()
    exCmplx = ps.CellComplex(exCells)

######## print skeletons of example:
#     for dm in range(len(example[-1])):
#         print dm," skeleton: ",
#         for idx in exCmplx.skeleton(dm):
#             print cells[idx].name,
#         print
#

######## print cofaces of each cell:
#     for cell in exCmplx.cells:
#         namedList = list(map(lambda cf: (cells[cf.index].name, cf.orientation), cell.cofaces))
#         print cell.name, namedList

######## print faces of a cell:  Note that self.faces() requires a cell index
    for cell in exCells:
        cdx = exidx[cell.name]
        namedList = list( map( lambda fc: exCells[fc].name, exCmplx.faces(cdx) ) )
        print "faces of ",cell.name, namedList
