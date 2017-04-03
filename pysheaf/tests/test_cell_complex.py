'''
Created on Mar 17, 2017

@author: bpraggastis
'''
import sys

sys.path.insert(0,'/Users/prag717/documents/tdm/gitrepos/pysheaf/pysheaf')

import unittest
import pysheaf as ps
import numpy as np


class TestCell(unittest.TestCase):
    def setUp(self):
        self.cell0 = ps.Cell(0,name='A')
        self.cell1 = ps.Cell(1,compactClosure=False)
        self.cell2 = ps.Cell(2,id=4)
        return

    def test_cell_instance(self):
        '''Test for cell instance'''
        self.assertTrue(isinstance(self.cell0,ps.Cell))
        self.assertTrue(isinstance(self.cell1,ps.Cell))
        self.assertTrue(isinstance(self.cell2,ps.Cell))
        return True

    def test_cell_attribute_assignment(self):
        '''Test for cell attributes'''
        self.assertEqual(self.cell0.dimension, 0)
        self.assertEqual(self.cell1.dimension, 1)
        self.assertEqual(self.cell0.name, 'A')
        self.assertEqual(self.cell2.name, None)
        self.assertEqual(self.cell2.compactClosure, True)
        self.assertEqual(self.cell1.compactClosure, False)
        self.assertEqual(self.cell0.id, None)
        self.assertEqual(self.cell2.id, 4)

    def test_cell_repr(self):
        '''Test for cell rep'''
        self.assertEqual(repr(self.cell0),'A (dimension=0,compactClosure=True)')

    def tearDown(self):
        unittest.TestCase.tearDown(self)

class TestCoface(unittest.TestCase):
    def setUp(self):
        ##   asc = [1,2,3,4,5] ####### list of cell indices
        self.coface0 = ps.Coface()
        self.coface1 = ps.Coface(2,1)
        self.coface2 = ps.Coface(3,-1)
        self.cell1 = ps.Cell(0,cofaces=[self.coface1,self.coface2])

    def test_coface_instance(self):
        '''Test for coface instance'''
        self.assertTrue(isinstance(self.coface0,ps.Coface))
        self.assertTrue(isinstance(self.coface1,ps.Coface))
        self.assertTrue(isinstance(self.coface2,ps.Coface))

    def test_coface_attribute_assignment(self):
        '''Test for coface attributes'''
        self.assertEqual(self.coface0.index, None)
        self.assertEqual(self.coface1.index, 2)
        self.assertEqual(self.coface1.orientation, 1)

    def test_coface_repr(self):
        '''Test for coface repr'''
        self.assertEqual(repr(self.coface2),'(index=3,orientation=-1)')

    def test_cell_isCoface(self):
        '''Test for cell is coface'''
        self.assertTrue(self.cell1.isCoface(3))
        self.assertTrue(self.cell1.isCoface(2,1))

    def tearDown(self):
        unittest.TestCase.tearDown(self)


class TestCellComplex(unittest.TestCase):
    def setUp(self):
        ### complex is: {A,C,E,K,AC,AE,AK,ACK}
        self.cellA = ps.Cell(0,name='A',id=0,cofaces=[ps.Coface(4,-1),ps.Coface(5, -1),ps.Coface(6, -1)])
        self.cellC = ps.Cell(0,name='C',id=1,cofaces=[ps.Coface(4, 1),ps.Coface(8, -1),ps.Coface(9, -1)])
        self.cellE = ps.Cell(0,name='E',id=2,cofaces=[ps.Coface(5, 1),ps.Coface(8, 1)])
        self.cellK = ps.Cell(0,name='K',id=3,cofaces=[ps.Coface(6, 1),ps.Coface(9, 1)])
        self.cellAC = ps.Cell(1,name='AC',id=4,cofaces=[ps.Coface(7, 1)])
        self.cellAE = ps.Cell(1,name='AE',id=5)
        self.cellAK = ps.Cell(1,name='AK',id=6,cofaces=[ps.Coface(7, -1)])
        self.cellACK = ps.Cell(2,name='ACK',id=7)
        self.cellCE = ps.Cell(1,name='CE',id=8)
        self.cellCK = ps.Cell(1,name='CK',id=9,cofaces=[ps.Coface(7, 1)])
        self.cmplx = ps.CellComplex([self.cellA,self.cellC,self.cellE,self.cellK,self.cellAC,self.cellAE,self.cellAK,self.cellACK,self.cellCE, self.cellCK])
        return

    def test_cellComplex_instance(self):
        '''Test for cell complex instance'''
        self.assertTrue(isinstance(self.cmplx,ps.CellComplex))

    def test_cellComplex_attribute_assignment(self):
        '''Test for cell complex attributes'''
        self.assertEqual(len(self.cmplx.cells),10)

    def test_cellComplex_add_cell(self):
        '''Test for cell complex add a single cell by instance returns cell id'''
        pass

    def test_cellComplex_add_cells_from(self):
        '''Test for cell complex add cells from a list of instances'''
        pass

    def test_cellComplex_add_coface(self):
        '''Test for cell complex add a single coface by index and orientation'''
        pass

    def test_add_cofaces_from(self):
        pass

    def test_isFaceOf(self):
        self.assertEqual(set(self.cmplx.isFaceOf(0)),set([4,5,6,7]))

    def test_skeleton(self):
        self.assertEqual(set(self.cmplx.skeleton(1)),set([4,5,6,8,9]))
        self.assertEqual(set(self.cmplx.skeleton(2)),set([7]))

    def test_faces(self):
        self.assertEqual(set(self.cmplx.faces(7)),set([0, 1, 3, 4, 6, 9]))
        self.assertEqual(set(self.cmplx.faces(6)),set([0,3]))

    def test_closure(self):
        self.assertEqual(set(self.cmplx.closure([6])),set([0,3,6]))
        self.assertEqual(self.cmplx.closure([0]),[0])

    def test_cofaces(self):  
        output = []
        for c in self.cmplx.cofaces(0):
            output.append(c)
        self.assertEqual(len(output),5)

    def test_components(self):
        self.assertEqual(len(self.cmplx.components()),1)
        self.assertEqual(len(self.cmplx.components([0,6,8])),2)

    def test_connectedTo(self):
        self.assertEqual(set(self.cmplx.connectedTo(0)),set([4,5,6,7]))
        self.assertEqual(set(self.cmplx.connectedTo(6)),set([0,3,7]))
        self.assertEqual(set(self.cmplx.connectedTo(6,[0,4,5])),set([0]))
        self.assertEqual(set(self.cmplx.connectedTo(0,[0,5,8,3])),set([5]))

    def test_expandComponent(self):
        self.assertEqual(len(self.cmplx.expandComponent(0)),len(self.cmplx.cells))
        self.assertEqual(set(self.cmplx.expandComponent(0,[0])),set([0]))
        self.assertEqual(set(self.cmplx.expandComponent(0,[0,2,5,8,9])),set([0,2,5,8]))

    def test_starCells(self):
        self.assertEqual(set(self.cmplx.starCells([0])),set([0,4,5,6,7]))
        self.assertEqual(set(self.cmplx.starCells([4])),set([4,7]))
        self.assertEqual(set(self.cmplx.starCells([1,3])),set([1,3,4,6,7,8,9]))

    def test_homology(self):
        self.assertEqual(self.cmplx.homology(1).shape[1],1)
        self.assertEqual(self.cmplx.homology(2).shape[1],0)

    def test_localPairComplex(self):
        pass

    def test_localHomology(self):
        pass

    def test_boundary(self):
        mat = self.cmplx.boundary(1)    
        self.assertTrue((mat[0] == [-1,-1,-1,0,0]).all())

    def test_inducedMapLocalHomology(self):
        pass

    def test_attachDiagram(self):
        pass

    def tearDown(self):
        unittest.TestCase.tearDown(self)

if __name__ == '__main__':
    unittest.main(verbosity=2)

