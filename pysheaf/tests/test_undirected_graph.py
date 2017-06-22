# -*- coding: utf-8 -*-
"""
Created on Wed May 31 18:10:17 2017

@author: jhenrich
"""

import unittest
import pysheaf as ps
import networkx as nx



class TestUndirectedGraph(unittest.TestCase):

    def setUp(self):
        ## Create the instance of the undirected graph class from a graph 
        #  whose list of edges is 
        #  [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4), (4,5), (5,6), (5,7), (6,7)] 
        #  or 
        #  [('A','B'), ('A','C'), ('A','D'), ('B','C'), ('B','D'), ('C','D'), 
        #   ('D','E'), ('E','F'), ('E','G'), ('F','G')]
        #  to check strings
        ## Edge List and Networkx Graph with nodes as integers:
        self.int_graph_edge_list = [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4), (4,5), (5,6), (5,7), (6,7)]
        self.int_graph = nx.Graph()
        self.int_graph.add_edges_from([(1,2), (1,3), (1,4), (2,3), (2,4), (3,4), (4,5), (5,6), (5,7), (6,7)])
        
        ## Edge List and Networkx Graph with nodes as strings:
        self.str_graph_edge_list = [('A','B'), ('A','C'), ('A','D'), ('B','C'), ('B','D'), ('C','D'), ('D','E'), ('E','F'), ('E','G'), ('F','G')]
        self.str_graph = nx.Graph()
        self.str_graph.add_edges_from([('A','B'), ('A','C'), ('A','D'), ('B','C'), ('B','D'), ('C','D'), ('D','E'), ('E','F'), ('E','G'), ('F','G')])
        
        ##Create undirected graphs in pysheaf:
        self.int_undir_graph_from_edge = ps.UndirectedGraph(self.int_graph_edge_list)
        self.int_undir_graph = ps.UndirectedGraph(self.int_graph)
        
        self.str_undir_graph_from_edge = ps.UndirectedGraph(self.str_graph_edge_list)
        self.str_undir_graph = ps.UndirectedGraph(self.str_graph)

    def test_undirected_graph_instance(self):
        self.assertTrue(isinstance(self.int_undir_graph_from_edge,ps.UndirectedGraph))
        self.assertTrue(isinstance(self.int_undir_graph,ps.UndirectedGraph))
        
        self.assertTrue(isinstance(self.str_undir_graph_from_edge, ps.UndirectedGraph))
        self.assertTrue(isinstance(self.str_undir_graph, ps.UndirectedGraph))
        
    def test_no_flag_int(self):
        f_edge = self.int_undir_graph_from_edge
        f_graph = self.int_undir_graph
        
        #correct answers
        c_num_cells = 17
        c_indexes_cells = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16']
        c_names_cells = [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4], [4, 5], [5, 6], [5, 7], [6, 7], [1], [2], [3], [4], [5], [6], [7]]
        c_cofaces_cells = '[[], [], [], [], [], [], [], [], [], [], [(index=0,orientation=None), (index=1,orientation=None), (index=2,orientation=None)], [(index=0,orientation=None), (index=3,orientation=None), (index=4,orientation=None)], [(index=1,orientation=None), (index=3,orientation=None), (index=5,orientation=None)], [(index=2,orientation=None), (index=4,orientation=None), (index=5,orientation=None), (index=6,orientation=None)], [(index=6,orientation=None), (index=7,orientation=None), (index=8,orientation=None)], [(index=7,orientation=None), (index=9,orientation=None)], [(index=8,orientation=None), (index=9,orientation=None)]]'
        
        #check number of cells
        self.assertEqual(len(f_edge.cells), c_num_cells)
        self.assertEqual(len(f_graph.cells), c_num_cells)
        
        #check the indices of the cells
        self.assertEqual([f_edge.cells[i].id for i in range(len(f_edge.cells))], c_indexes_cells)
        self.assertEqual([f_graph.cells[i].id for i in range(len(f_graph.cells))], c_indexes_cells)
        
        #check the cell names
        self.assertEqual([f_edge.cells[i].name for i in range(len(f_edge.cells))], c_names_cells)
        self.assertEqual([f_graph.cells[i].name for i in range(len(f_graph.cells))], c_names_cells)
        
        #check the cell cofaces
        self.assertEqual(str([f_edge.cells[i].cofaces for i in range(len(f_edge.cells))]), c_cofaces_cells)
        self.assertEqual(str([f_graph.cells[i].cofaces for i in range(len(f_graph.cells))]), c_cofaces_cells)
        
        
    def test_no_flag_str(self):
        f_edge = self.str_undir_graph_from_edge
        f_graph = self.str_undir_graph
        
        #correct answers
        c_num_cells = 17
        c_indexes_cells = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16']
        c_names_cells = [['A', 'B'], ['A', 'C'], ['A', 'D'], ['B', 'C'], ['B', 'D'], ['C', 'D'], ['D', 'E'], ['E', 'F'], ['E', 'G'], ['F', 'G'], ['A'], ['B'], ['C'], ['D'], ['E'], ['F'], ['G']]
        c_cofaces_cells = '[[], [], [], [], [], [], [], [], [], [], [(index=0,orientation=None), (index=1,orientation=None), (index=2,orientation=None)], [(index=0,orientation=None), (index=3,orientation=None), (index=4,orientation=None)], [(index=1,orientation=None), (index=3,orientation=None), (index=5,orientation=None)], [(index=2,orientation=None), (index=4,orientation=None), (index=5,orientation=None), (index=6,orientation=None)], [(index=6,orientation=None), (index=7,orientation=None), (index=8,orientation=None)], [(index=7,orientation=None), (index=9,orientation=None)], [(index=8,orientation=None), (index=9,orientation=None)]]'
        
        #check number of cells
        self.assertEqual(len(f_edge.cells), c_num_cells)
        self.assertEqual(len(f_graph.cells), c_num_cells)
        
        #check the indices of the cells
        self.assertEqual([f_edge.cells[i].id for i in range(len(f_edge.cells))], c_indexes_cells)
        self.assertEqual([f_graph.cells[i].id for i in range(len(f_graph.cells))], c_indexes_cells)
        
        #check the cell names
        self.assertEqual([f_edge.cells[i].name for i in range(len(f_edge.cells))], c_names_cells)
        self.assertEqual([f_graph.cells[i].name for i in range(len(f_graph.cells))], c_names_cells)
        
        #check the cell cofaces
        self.assertEqual(str([f_edge.cells[i].cofaces for i in range(len(f_edge.cells))]), c_cofaces_cells)
        self.assertEqual(str([f_graph.cells[i].cofaces for i in range(len(f_graph.cells))]), c_cofaces_cells)
        
        

    def tearDown(self):
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()