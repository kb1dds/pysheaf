# MIT License

# Copyright (c) 2018 Janelle Henrich & Michael Robinson & Steven Fiacco

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

import numpy as np
import pysheaf as ps
import copy
import random
from functools import partial
from deap import base
from deap import creator
from deap import tools
from deap import algorithms

def Pow2(inputNumber):
   return inputNumber**2 # Pow2
def Pow3(inputNumber):
   return inputNumber**3 # Pow3
def CompareScalars(leftValue,rightValue):
   return abs(leftValue - rightValue) # CompareScalars
def SerializeScalars(assignmentValue):
   return np.array([assignmentValue]) # SerializeScalars


def arithXover(ind1, ind2):
   for i, (x1, x2) in enumerate(zip(ind1, ind2)):
       gamma = random.random()
       ind1[i] = (1. - gamma) * x1 + gamma * x2
       ind2[i] = gamma * x1 + (1. - gamma) * x2

   return ind1, ind2

def selnormGeom(individuals, k, prob_sel_best= 0.08, fit_attr="fitness"):
   #NormGeomSelect is a ranking selection function based on the normalized
   #geometric distribution.  

   #Modified from the Matlab version into the style of DEAP

   probabilty_of_selecting_most_fit_individual = prob_sel_best   # Probability of selecting the best
   population_size = len(individuals)  # Number of individuals in pop


   chosen = [] #editted the structure of th output to reflect the structure of pysheaf
   fit = np.zeros((population_size,1))  #Allocates space for the prop of select
   population_index_and_rank = np.zeros((population_size,2))    #Sorted list of id and rank
   population_index_and_rank[:, 0] = list(range(population_size,0,-1)) # Get original location of the element
   to_sort = list(zip(individuals, list(range(population_size)))) #need to keep original index associated
   s_inds = sorted(to_sort, key= lambda ind: getattr(ind[0], fit_attr).values[0]) #Sort by the fitnesses
   population_index_and_rank[:, 1] = [b for a,b in s_inds]
   normalization_factor =probabilty_of_selecting_most_fit_individual/(1-((1-probabilty_of_selecting_most_fit_individual)**population_size))  # normalize the distribution, q prime
   for population_index in range(population_size):  #Generate the probability of selection
       ind_fit = int(population_index_and_rank[population_index,1])
       fit[ind_fit] = normalization_factor*((1-probabilty_of_selecting_most_fit_individual)**(population_index_and_rank[population_index, 0]-1))
   fit = np.cumsum(fit)  # Calculate the cummulative prob. function
   rnums = sorted([random.random() for nn in range(population_size)])  # Generate n sorted random numbers
   fitIn = 0
   new_In = 0
   unique = []
   while new_In < k:
       if rnums[new_In] < fit[fitIn]:
           unique.append(fitIn)
           chosen.append(individuals[fitIn]) #Select the fitIn individual
           new_In += 1  # Looking for next new individual
       else:
           fitIn += 1 # Looking at next potential selection

   return chosen 

#Define a function to do the mutation of the elements from a floating point perspective
def mutUniformFloat(individual, lowerBound, upperBound, indpb):
   """Mutate an individual by replacing attributes, with probability *indpb*,
   by a integer uniformly drawn between *low* and *up* inclusively.
   
   :param individual: :term:`Sequence <sequence>` individual to be mutated.
   :param low: The lower bound or a :term:`python:sequence` of
               of lower bounds of the range from wich to draw the new
               integer.
   :param upperBound: The upper bound or a :term:`python:sequence` of
              of upper bounds of the range from wich to draw the new
              integer.
   :param indpb: Independent probability for each attribute to be mutated.
   :returns: A tuple of one individual.
   """
   size = len(individual)
   if len(lowerBound) < size:
       raise IndexError("lowerBound must be at least the size of individual: %d < %d" % (len(lowerBound), size))
   elif len(upperBound) < size:
       raise IndexError("upperBound must be at least the size of individual: %d < %d" % (len(upperBound), size))
   
   for i, xl, xu in zip(range(size), lowerBound, upperBound):
       if random.random() < indpb:
           individual[i] = random.uniform(xl, xu)
   
   return individual,


def DeapSheafOptimizer(functionToIterate, serializedAssignments, boundsList,maxIterations):
   initial_population_size = 100
   mutation_rate =0.3
   number_of_generations = 100

   lower_bounds = []
   upper_bounds = []
   for i,bnds in enumerate(boundsList):
       if bnds[0] is not None:
           lower_bounds.append(float(bnds[0]))
       else:
           #"Bounds not specified on an activeCell"
           if serializedAssignments is not None:
               lower_bounds.append(serializedAssignments[i]-100*serializedAssignments[i])
           else:
               raise ValueError("Initial guess can't be turned into bounds")
       if bnds[1] is not None:
           upper_bounds.append(float(bnds[1]))
       else:
           #"Bounds not specified on an activeCell"
           if serializedAssignments is not None:
               upper_bounds.append(serializedAssignments[i]+100*serializedAssignments[i])
           else:
               raise ValueError("Initial guess can't be turned into bounds")

   def GetCostFromSheaf(arg):
      if np.any(np.isnan(arg)):
         return tuple([-1e100])

      return tuple([-1.0*functionToIterate(arg)])


   creator.create("FitnessMax", base.Fitness, weights=(1.0,))
   creator.create("Individual", np.ndarray, fitness=creator.FitnessMax)
   toolbox = base.Toolbox()

   #generation_function_list forces the initial population to be generated randomly within the bounds
   generation_function_list = []
   for boundIndex in range(len(boundsList)):
       if lower_bounds[boundIndex] != None and upper_bounds[boundIndex] != None:
           generation_function_list.append(partial(random.uniform, lower_bounds[boundIndex], upper_bounds[boundIndex]))
       elif lower_bounds[boundIndex] == None and upper_bounds[boundIndex] == None:
           generation_function_list.extend([lambda:(1/(1-random.random()))-1]) #maps [0,1) to [0, inf)
       elif lower_bounds[boundIndex] == None:
           multiply_bnds1 = partial(np.multiply, upper_bounds[boundIndex])
           generation_function_list.extend([lambda:(-1/(1-random.random()) + 1 + multiply_bnds1(1.0))])
           #generation_function_list.extend([lambda:(random.randrange(copy.deepcopy(bnds[1])))]) #need to check actual opperation of randrange without a start
       else:
           multiply_bnds0 = partial(np.multiply, lower_bounds[boundIndex])
           generation_function_list.extend([lambda: (1/(1-random.random()))-1 + multiply_bnds0(1.0)])
           #generation_function_list.extend([lambda:(-1*random.randrange(copy.deepcopy(bnds[0])) + 2*copy.deepcopy(bnds[0]))])
   
   #specify a population within the bounds
   toolbox.register("individual", tools.initCycle, creator.Individual, generation_function_list, n=1)
   toolbox.register("population", tools.initRepeat, list, toolbox.individual)

   toolbox.register("evaluate", GetCostFromSheaf)
   toolbox.register("mutate", mutUniformFloat,lowerBound=lower_bounds,upperBound=upper_bounds,indpb=mutation_rate)
   toolbox.register("mate", arithXover)
   toolbox.register("select", selnormGeom, prob_sel_best=0.08)

   stats = tools.Statistics(lambda ind: ind.fitness.values)
   stats.register("avg", np.mean)
   stats.register("std", np.std)
   stats.register("min", np.min)
   stats.register("max", np.max)

   mostFitIndividual = tools.HallOfFame(1,similar=np.array_equal)

   starting_pop = toolbox.population(n=initial_population_size)
   pop, log = algorithms.eaMuPlusLambda(starting_pop, toolbox, mu = (3*initial_population_size), lambda_=(3*initial_population_size),  cxpb=0.5, mutpb=0.5, ngen=number_of_generations, stats=stats, halloffame=mostFitIndividual)

   return mostFitIndividual# DeapSheafOptimizer


if __name__ == '__main__':
   graph = ps.Sheaf()
   graph.mPreventRedundantExtendedAssignments = False

   TEST_TYPE = "test_type"
   graph.AddCell(0,ps.Cell(TEST_TYPE,CompareScalars,serializeAssignmentMethod= SerializeScalars))
   graph.AddCell(1,ps.Cell(TEST_TYPE,CompareScalars,serializeAssignmentMethod= SerializeScalars))
   graph.AddCell(2,ps.Cell(TEST_TYPE,CompareScalars,serializeAssignmentMethod= SerializeScalars))
   graph.AddCell(3,ps.Cell(TEST_TYPE,CompareScalars,serializeAssignmentMethod= SerializeScalars))
   graph.AddCell(4,ps.Cell(TEST_TYPE,CompareScalars,serializeAssignmentMethod= SerializeScalars))
   graph.AddCell(5,ps.Cell(TEST_TYPE,CompareScalars,serializeAssignmentMethod= SerializeScalars))

   assert(graph.number_of_nodes() == 6), "Incorrect number of cells"

   """
   Test Structure

         0
         |
         1
         |
         2
        / \
       4   3
        \ /
         5
   """

   coface_square = ps.Coface(TEST_TYPE,TEST_TYPE,Pow2)
   graph.AddCoface(0,1,coface_square)
   coface_square = ps.Coface(TEST_TYPE,TEST_TYPE,Pow2)
   graph.AddCoface(1,2,coface_square)
   coface_square = ps.Coface(TEST_TYPE,TEST_TYPE,Pow2)
   graph.AddCoface(2,3,coface_square)
   coface_square = ps.Coface(TEST_TYPE,TEST_TYPE,Pow2)
   graph.AddCoface(2,4,coface_square)
   coface_square = ps.Coface(TEST_TYPE,TEST_TYPE,Pow2)
   graph.AddCoface(3,5,coface_square)
   coface_square = ps.Coface(TEST_TYPE,TEST_TYPE,Pow3)
   graph.AddCoface(4,5,coface_square)

   known_test_values = [2,4,16,256,65536,16777216]
   graph.GetCell(0).SetDataAssignment(ps.Assignment(TEST_TYPE,2))


   graph.GetCell(0).mOptimizationCell = True
   graph.mPreventRedundantExtendedAssignments = True

   cell2_test_value = 81
   cell5_test_value = 43046721
   result_answer = 3

   graph.GetCell(2).SetDataAssignment(ps.Assignment(TEST_TYPE,cell2_test_value))
   graph.GetCell(5).SetDataAssignment(ps.Assignment(TEST_TYPE,cell5_test_value))

   graph.ClearExtendedAssignments()
   graph.mNumpyNormType = None
   graph.MaximallyExtendCell(0)
   print("consistency radius before fuse assignment:",graph.ComputeConsistencyRadius())


   graph.SetSheafOptimizer(DeapSheafOptimizer)
   fuse_results = graph.FuseAssignment()
   print(fuse_results)