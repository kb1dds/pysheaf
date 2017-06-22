# -*- coding: utf-8 -*-
"""
Created on Wed Jun 07 13:48:57 2017

@author: jhenrich
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May 31 18:10:17 2017

@author: jhenrich
"""

import unittest
import pysheaf as ps
import networkx as nx



class Test_Flag_Complex(unittest.TestCase):

    def setUp(self):
        ## Create the instance of the undirected graph class from a graph 
        #  whose list of edges is 
        #  [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4), (4,5), (5,6), (5,7), (6,7)] 
        #  or 
        #  [('A','B'), ('A','C'), ('A','D'), ('B','C'), ('B','D'), ('C','D'), 
        #   ('D','E'), ('E','F'), ('E','G'), ('F','G')]
        #  to check strings
        ## Edge List and Networkx Graph with nodes as integers:
        int_graph_edge_list = [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4), (4,5), (5,6), (5,7), (6,7)]
        int_graph = nx.Graph()
        int_graph.add_edges_from([(1,2), (1,3), (1,4), (2,3), (2,4), (3,4), (4,5), (5,6), (5,7), (6,7)])
        #int_graph.add_node(8)
        
        ## Edge List and Networkx Graph with nodes as strings:
        str_graph_edge_list = [('A','B'), ('A','C'), ('A','D'), ('B','C'), ('B','D'), ('C','D'), ('D','E'), ('E','F'), ('E','G'), ('F','G')]
        str_graph = nx.Graph()
        str_graph.add_edges_from([('A','B'), ('A','C'), ('A','D'), ('B','C'), ('B','D'), ('C','D'), ('D','E'), ('E','F'), ('E','G'), ('F','G')])
        
        ##Create graphs in pysheaf to test formation of flag complex
        int_undir_graph = ps.UndirectedGraph(int_graph)
        
        str_undir_graph = ps.UndirectedGraph(str_graph)
        
        ##Create flag complexs in pysheaf:
        self.int_flag_complex_from_edges = ps.FlagComplex(int_graph_edge_list)
        self.int_flag_complex_from_nx_graph = ps.FlagComplex(int_graph)
        self.int_flag_complex_from_ps_undirectedgraph = ps.FlagComplex(int_undir_graph)
        
        
        self.str_flag_complex_from_edges = ps.FlagComplex(str_graph_edge_list)
        self.str_flag_complex_from_nx_graph = ps.FlagComplex(str_graph)
        self.str_flag_complex_from_ps_undirectedgraph = ps.FlagComplex(str_undir_graph)

    def test_flag_instance(self):
        self.assertTrue(isinstance(self.int_flag_complex_from_edges,ps.FlagComplex))
        self.assertTrue(isinstance(self.int_flag_complex_from_nx_graph, ps.FlagComplex))
        self.assertTrue(isinstance(self.int_flag_complex_from_ps_undirectedgraph,ps.FlagComplex))
        
        self.assertTrue(isinstance(self.str_flag_complex_from_edges, ps.FlagComplex))
        self.assertTrue(isinstance(self.str_flag_complex_from_nx_graph, ps.FlagComplex))
        self.assertTrue(isinstance(self.str_flag_complex_from_ps_undirectedgraph, ps.FlagComplex))
              
    def test_flag_int(self):
        f_edge = self.int_flag_complex_from_edges
        f_nx_graph = self.int_flag_complex_from_nx_graph
        f_undir_graph = self.int_flag_complex_from_ps_undirectedgraph
        
        
        #correct answers
        c_num_cells = 23
        c_indexes_cells = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
        c_names_cells = [[1, 2, 3, 4], [1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4], [5, 6, 7], [1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4], [4, 5], [5, 6], [5, 7], [6, 7], [1], [2], [3], [4], [5], [6], [7]]
        c_cofaces_cells = '[[], [(index=0,orientation=-1)], [(index=0,orientation=1)], [(index=0,orientation=-1)], [(index=0,orientation=1)], [], [(index=1,orientation=1), (index=2,orientation=1)], [(index=1,orientation=-1), (index=3,orientation=1)], [(index=2,orientation=-1), (index=3,orientation=-1)], [(index=1,orientation=1), (index=4,orientation=1)], [(index=2,orientation=1), (index=4,orientation=-1)], [(index=3,orientation=1), (index=4,orientation=1)], [], [(index=5,orientation=1)], [(index=5,orientation=-1)], [(index=5,orientation=1)], [(index=6,orientation=-1), (index=7,orientation=-1), (index=8,orientation=-1)], [(index=6,orientation=1), (index=9,orientation=-1), (index=10,orientation=-1)], [(index=7,orientation=1), (index=9,orientation=1), (index=11,orientation=-1)], [(index=8,orientation=1), (index=10,orientation=1), (index=11,orientation=1), (index=12,orientation=-1)], [(index=12,orientation=1), (index=13,orientation=-1), (index=14,orientation=-1)], [(index=13,orientation=1), (index=15,orientation=-1)], [(index=14,orientation=1), (index=15,orientation=1)]]'
        
        #check number of cells
        self.assertEqual(len(f_edge.cells), c_num_cells)
        self.assertEqual(len(f_nx_graph.cells), c_num_cells)
        self.assertEqual(len(f_undir_graph.cells), c_num_cells)
        
        
        #check the indices of the cells
        self.assertEqual([f_edge.cells[i].id for i in range(len(f_edge.cells))], c_indexes_cells)
        self.assertEqual([f_nx_graph.cells[i].id for i in range(len(f_nx_graph.cells))], c_indexes_cells)
        self.assertEqual([f_undir_graph.cells[i].id for i in range(len(f_undir_graph.cells))], c_indexes_cells)
        
        
        #check the cell names
        
        self.assertEqual([f_edge.cells[i].name for i in range(len(f_edge.cells))], c_names_cells)
        self.assertEqual([f_nx_graph.cells[i].name for i in range(len(f_nx_graph.cells))], c_names_cells)
        self.assertEqual([f_undir_graph.cells[i].name for i in range(len(f_undir_graph.cells))], c_names_cells)
        
        
        #check the cell cofaces
        self.assertEqual(str([f_edge.cells[i].cofaces for i in range(len(f_edge.cells))]), c_cofaces_cells)
        self.assertEqual(str([f_nx_graph.cells[i].cofaces for i in range(len(f_nx_graph.cells))]), c_cofaces_cells)
        self.assertEqual(str([f_undir_graph.cells[i].cofaces for i in range(len(f_undir_graph.cells))]), c_cofaces_cells)
        
        
        
    def test_flag_str(self):
        f_edge = self.str_flag_complex_from_edges
        f_nx_graph = self.str_flag_complex_from_nx_graph
        f_undir_graph = self.str_flag_complex_from_ps_undirectedgraph
        
        #correct answers
        c_num_cells = 23
        c_indexes_cells = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
        c_names_cells = [['A', 'B', 'C', 'D'], ['A', 'B', 'C'], ['A', 'B', 'D'], ['A', 'C', 'D'], ['B', 'C', 'D'], ['E', 'F', 'G'], ['A', 'B'], ['A', 'C'], ['A', 'D'], ['B', 'C'], ['B', 'D'], ['C', 'D'], ['D', 'E'], ['E', 'F'], ['E', 'G'], ['F', 'G'], ['A'], ['B'], ['C'], ['D'], ['E'], ['F'], ['G']]
        c_cofaces_cells = '[[], [(index=0,orientation=-1)], [(index=0,orientation=1)], [(index=0,orientation=-1)], [(index=0,orientation=1)], [], [(index=1,orientation=1), (index=2,orientation=1)], [(index=1,orientation=-1), (index=3,orientation=1)], [(index=2,orientation=-1), (index=3,orientation=-1)], [(index=1,orientation=1), (index=4,orientation=1)], [(index=2,orientation=1), (index=4,orientation=-1)], [(index=3,orientation=1), (index=4,orientation=1)], [], [(index=5,orientation=1)], [(index=5,orientation=-1)], [(index=5,orientation=1)], [(index=6,orientation=-1), (index=7,orientation=-1), (index=8,orientation=-1)], [(index=6,orientation=1), (index=9,orientation=-1), (index=10,orientation=-1)], [(index=7,orientation=1), (index=9,orientation=1), (index=11,orientation=-1)], [(index=8,orientation=1), (index=10,orientation=1), (index=11,orientation=1), (index=12,orientation=-1)], [(index=12,orientation=1), (index=13,orientation=-1), (index=14,orientation=-1)], [(index=13,orientation=1), (index=15,orientation=-1)], [(index=14,orientation=1), (index=15,orientation=1)]]'
        
        #check number of cells
        self.assertEqual(len(f_edge.cells), c_num_cells)
        self.assertEqual(len(f_nx_graph.cells), c_num_cells)
        self.assertEqual(len(f_undir_graph.cells), c_num_cells)
        
        #check the indices of the cells
        self.assertEqual([f_edge.cells[i].id for i in range(len(f_edge.cells))], c_indexes_cells)
        self.assertEqual([f_nx_graph.cells[i].id for i in range(len(f_nx_graph.cells))], c_indexes_cells)
        self.assertEqual([f_undir_graph.cells[i].id for i in range(len(f_undir_graph.cells))], c_indexes_cells)
        
        #check the cell names
        self.assertEqual([f_edge.cells[i].name for i in range(len(f_edge.cells))], c_names_cells)
        self.assertEqual([f_nx_graph.cells[i].name for i in range(len(f_nx_graph.cells))], c_names_cells)
        self.assertEqual([f_undir_graph.cells[i].name for i in range(len(f_undir_graph.cells))], c_names_cells)
        
        #check the cell cofaces
        self.assertEqual(str([f_edge.cells[i].cofaces for i in range(len(f_edge.cells))]), c_cofaces_cells)
        self.assertEqual(str([f_nx_graph.cells[i].cofaces for i in range(len(f_nx_graph.cells))]), c_cofaces_cells)
        self.assertEqual(str([f_undir_graph.cells[i].cofaces for i in range(len(f_undir_graph.cells))]), c_cofaces_cells)
        
    
    def test_flag_complex_betti(self):
        betti_0 = float(1) 
        betti_1 = float(0)
        betti_2 = float(0)
        betti_3 = float(0)
        
        print self.int_flag_complex_from_edges.betti(2)
        
        #Check 0th Betti
        self.assertEqual(float(self.int_flag_complex_from_edges.betti(0)), betti_0)
        self.assertEqual(float(self.int_flag_complex_from_nx_graph.betti(0)), betti_0)
        self.assertEqual(float(self.int_flag_complex_from_ps_undirectedgraph.betti(0)), betti_0)
        self.assertEqual(float(self.str_flag_complex_from_edges.betti(0)), betti_0)
        self.assertEqual(float(self.str_flag_complex_from_nx_graph.betti(0)), betti_0)
        self.assertEqual(float(self.str_flag_complex_from_ps_undirectedgraph.betti(0)), betti_0)
        
        #Check 1st Betti
        self.assertEqual(float(self.int_flag_complex_from_edges.betti(1)), betti_1)
        self.assertEqual(float(self.int_flag_complex_from_nx_graph.betti(1)), betti_1)
        self.assertEqual(float(self.int_flag_complex_from_ps_undirectedgraph.betti(1)), betti_1)
        self.assertEqual(float(self.str_flag_complex_from_edges.betti(1)), betti_1)
        self.assertEqual(float(self.str_flag_complex_from_nx_graph.betti(1)), betti_1)
        self.assertEqual(float(self.str_flag_complex_from_ps_undirectedgraph.betti(1)), betti_1)
        
        #Check 2nd Betti
        self.assertEqual(float(self.int_flag_complex_from_edges.betti(2)), betti_2)
        self.assertEqual(float(self.int_flag_complex_from_nx_graph.betti(2)), betti_2)
        self.assertEqual(float(self.int_flag_complex_from_ps_undirectedgraph.betti(2)), betti_2)
        self.assertEqual(float(self.str_flag_complex_from_edges.betti(2)), betti_2)
        self.assertEqual(float(self.str_flag_complex_from_nx_graph.betti(2)), betti_2)
        self.assertEqual(float(self.str_flag_complex_from_ps_undirectedgraph.betti(2)), betti_2)
        
        #Check 3rd Betti
        self.assertEqual(float(self.int_flag_complex_from_edges.betti(3)), betti_3)
        self.assertEqual(float(self.int_flag_complex_from_nx_graph.betti(3)), betti_3)
        self.assertEqual(float(self.int_flag_complex_from_ps_undirectedgraph.betti(3)), betti_3)
        self.assertEqual(float(self.str_flag_complex_from_edges.betti(3)), betti_3)
        self.assertEqual(float(self.str_flag_complex_from_nx_graph.betti(3)), betti_3)
        self.assertEqual(float(self.str_flag_complex_from_ps_undirectedgraph.betti(3)), betti_3)
        

    def tearDown(self):
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()