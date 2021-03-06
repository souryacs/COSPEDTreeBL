#!/usr/bin/env python

import dendropy
from dendropy import TreeList, Tree, Taxon, TaxonSet, Node
import itertools
import numpy
import time
import os
import sys
from cStringIO import StringIO
from optparse import OptionParser
import math
import scipy
from scipy import optimize 
from scipy.optimize import minimize

# we define custom edge types
RELATION_R3 = 0	# equality relationship
RELATION_R1 = 1
RELATION_R2 = 2
RELATION_R4 = 3	# no relationship
UNDEFINED_RELATION = 4

""" this is a dictionary storing cluster of nodes 
each cluster is basically a collection of nodes having equality relationship between the nodes """
Cluster_Info_Dict = dict()

""" the dictionary defines one particular taxa 
individual taxa contains the relationship information with other taxa """
Taxa_Info_Dict = dict()

""" this dictionary defines the taxa pair relations
each entry is indexed by two nodes """
TaxaPair_Reln_Dict = dict()

''' this is the list storing the costs of relations for taxa pairs depicting multi relation instance
it stores the taxa pair and the corresponding scores of different relation types
the list is used to extract the maximum score (corresponding to a particular relation) for resolving that taxa pair '''
Cost_List_Taxa_Pair_Multi_Reln = []

''' if SINGLE_EDGE_TYPE_CONN_PRIORITY option is set, then this list stores the cases 
where a pair of taxa is connected only by a single relation type  (with respect to all the candidate source trees)
cost corresponding to that single relation instance is accounted  '''
Cost_List_Taxa_Pair_Single_Reln = [] 

""" this list contains the complete set of taxa present in the input source trees """
COMPLETE_INPUT_TAXA_LIST = []

""" this list contains the current set of active cluster indices """
CURRENT_CLUST_IDX_LIST = []

# this is the debug level
# set for printing the necessary information
DEBUG_LEVEL = 0

# this text file stores all the printing output
Output_Text_File = 'complete_output_description.txt'

# this is a global string containing the objective function of the quadratic 
# programming employed for branch length determination of the tree
objective_function_string = ''

# this is a global string containing the derivative of the objective function
# employed for the above mentioned quadratic programming
deriv_objective_function_string = ''

# this dictionary stores the edge information
# key of the dictionary: (terminal_node1, terminal_node2)
# value of the element: array index (which will be used to denote corresponding variable with respect to the branch)
# variables are stored as x[0], x[1], .....
# they are employed in a quadratic programming
EdgeInfoDict = dict()

""" this variable associates weight of matrix corresponding to individual source trees """
Matrix_Weight_Val = []

# this is the corrected edge length for the output weighted supertree
# when the QP computation returns a negative edge length
CORRECTED_POSITIVE_EDGE_LEN = 0.00001

##-----------------------------------------------------
''' this class defines a leaf node of input candidate source trees
that is, it corresponds to one taxa '''
class Single_Taxa(object):
  def __init__(self):
    """ this variable signifies the cluster index that the taxa belongs 
    if the value is -1 (< 0) then it is still not inserted in a cluster 
    otherwise it is part of a valid cluster """
    self.clust_idx_part = -1

  def _Get_Taxa_Part_Clust_Idx(self):
    return self.clust_idx_part
    
  def _Set_Clust_Idx_taxa_Part(self, inp_clust_idx):
    self.clust_idx_part = inp_clust_idx
      
  # this function is called after formation of consensus tree
  def _PrintFinalTaxaInfo(self, key, Output_Text_File):
    fp = open(Output_Text_File, 'a')    
    fp.write('\n taxa key: ' + str(key))
    fp.write('\n taxa is part of the cluster ID: ' + str(self.clust_idx_part))
    fp.close()
        
    
##-----------------------------------------------------
""" this class defines the connectivity relationship between a pair of taxa
initially the information are obtained from the input source trees
later the contents of these class instances are modified according to the generation of the consensus tree
key of this class --- taxa1, taxa2  
in the class, the edge type signifies the relation between a pair of taxa """
class Reln_TaxaPair(object):
  def __init__(self):    
    
    ''' this variable denotes the no of occurrences of a particular edge type 
    there are 4 types of edges (relationship) between a pair of taxa '''
    self.edge_weight = [0] * 4    
    
    ''' a connection priority value is defined as the 
    no of occurrences of this particular edge type between these pair of nodes 
    minus the sum of no of occurrences of other edge types between these pair of nodes '''
    self.conn_pr_val = [0] * 4    
    
    ''' this cost variable denotes the cost associated with different types of edge connection between 
    these pair of nodes considered 
    this value is updated during generation of the consensus tree '''
    self.Connect_Edge_Cost = [0] * 4
    
    ''' this variable stores the number of extra lineages for the MRCA 
    node of this couplet, accumulated over all input trees '''
    self.xl_sum_all_trees = 0
    
    # for a particular taxa pair, when the tree topology is established, this array stires the branches 
    # those are in between this taxa pair
    self.Branch_Array_Idx_List = []
    
    # this list stores the tree index where the taxa pair is supported
    self.Supported_Tree_Idx_List = []

    # this list stores the distance values between these taxa pairs
    # as supported by the source trees
    self.DistMatValues = []
    
    # add - sourya
    # for one node, sum of distance to the MRCA node, corresponding to a particular relation to the other node
    self.dist_node1_mrca_sum_per_reln = [0] * 4
    self.dist_node2_mrca_sum_per_reln = [0] * 4    
    # end add - sourya
    
    # for one node, mean of its distance to the MRCA node, corresponding to a particular relation to the other node
    self.dist_node1_mrca_mean_per_reln = [0] * 4
    self.dist_node2_mrca_mean_per_reln = [0] * 4    
    
    # global mean distance of the node to the MRCA node
    self.dist_node1_mrca_global_mean = 0
    self.dist_node2_mrca_global_mean = 0
        
    # this list maintains the consensus relations for this couplet
    # with respect to the input trees
    self.consensus_reln_list = []
        
    # this list contains the union of taxa list underlying the LCA of this couplet
    # for individual input trees
    self.LCA_Underlying_Taxa_List = []

  def _AppendUnderlyingTaxonList(self, inp_list):
    self.LCA_Underlying_Taxa_List = list(set(self.LCA_Underlying_Taxa_List) | set(inp_list))

  def _GetUnderlyingTaxonList(self):
    return self.LCA_Underlying_Taxa_List
        
  def _GetConsensusRelnList(self):
    return self.consensus_reln_list
        
  #---------------------------------------------------
  def _AddLineage(self, val):
    self.xl_sum_all_trees = self.xl_sum_all_trees + val
  
  def _GetAvgTreeXL(self):
    return (self.xl_sum_all_trees * 1.0) / len(self.Supported_Tree_Idx_List)
  
  #---------------------------------------------------
  # these functions are for computing the QP 
  def _AddDistMatValue(self, val):
    self.DistMatValues.append(val)
    
  def _GetNoSupportTrees(self):
    return len(self.Supported_Tree_Idx_List)
    
  def _GetBranchArrayIdxList(self):
    return self.Branch_Array_Idx_List
    
  def _AddBranchArrayIdx(self, idx):
    if idx not in self.Branch_Array_Idx_List:
      self.Branch_Array_Idx_List.append(idx)
      
  def _AddSupportTreeIndex(self, idx):
    self.Supported_Tree_Idx_List.append(idx)
        
  def _GetSupportTreeList(self):
    return self.Supported_Tree_Idx_List
    
  # this function returns the average value of distance between current taxa pair
  def _GetAvgDistMatVal(self):
    # comment - sourya
    #return (self.dist_node1_mrca_global_mean + self.dist_node2_mrca_global_mean)
    # add - sourya
    num = 0
    denom = 0
    for i in range(len(self.DistMatValues)):
      num = num + self.DistMatValues[i] * Matrix_Weight_Val[self.Supported_Tree_Idx_List[i]]
      denom = denom + Matrix_Weight_Val[self.Supported_Tree_Idx_List[i]]
    return (num * 1.0) / denom
    # end add - sourya
  #---------------------------------------------------
  
  def _GetEdgeWeight(self, edge_type):
    return self.edge_weight[edge_type]      
    
  def _GetEdgeCost_ConnReln(self, edge_type):
    return self.Connect_Edge_Cost[edge_type]
    
  #def _IncrEdgeCost_ConnReln(self, edge_type, incr_cost):
    #self.Connect_Edge_Cost[edge_type] = self.Connect_Edge_Cost[edge_type] + incr_cost

  # this function adds one edge count (with a given input edge type)
  def _Add_Edge_Count_And_Distance(self, edge_type, node1_dist_from_mrca_node, node2_dist_from_mrca_node, val=1):
    # increment the frequency of particular relation type between this taxa pair
    self.edge_weight[edge_type] = self.edge_weight[edge_type] + val
    # add the relative branch length timings (with respect to MRCA)
    self._AddDistMatValue(node1_dist_from_mrca_node + node2_dist_from_mrca_node)
    # add respective branch lengths leading to individual taxon
    self.dist_node1_mrca_sum_per_reln[edge_type] = self.dist_node1_mrca_sum_per_reln[edge_type] + node1_dist_from_mrca_node
    self.dist_node2_mrca_sum_per_reln[edge_type] = self.dist_node2_mrca_sum_per_reln[edge_type] + node2_dist_from_mrca_node
        
  # this function prints the relationship information
  def _PrintRelnInfo(self, key, Output_Text_File):
    fp = open(Output_Text_File, 'a')    
    fp.write('\n taxa pair key: ' + str(key))
    fp.write('\n edges [type/count/conn_pr_val]: ')
    for i in range(4):
      fp.write('\n [' + str(i) + '/' + str(self.edge_weight[i]) + '/' + str(self.conn_pr_val[i]) + ']')
    fp.close()
  
  #------------------------------------
  # this function determines the mean distance between this pair of nodes, with respect to their MRCA node
  # mean distance with respect to individual couplet relation type
  # plus the overall mean distance between this couplet is computed
  def _CalcMeanDist(self):
    # mean distance with respect to individual relation types
    # between individual taxon and the MRCA node
    for edge_type in range(4):
      if (self.edge_weight[edge_type] > 0):
	self.dist_node1_mrca_mean_per_reln[edge_type] = (self.dist_node1_mrca_sum_per_reln[edge_type] * 1.0) / self.edge_weight[edge_type]
	self.dist_node2_mrca_mean_per_reln[edge_type] = (self.dist_node2_mrca_sum_per_reln[edge_type] * 1.0) / self.edge_weight[edge_type]
      else:
	self.dist_node1_mrca_mean_per_reln[edge_type] = 0
	self.dist_node2_mrca_mean_per_reln[edge_type] = 0
	
    # global mean distance between individual taxon and the MRCA node
    self.dist_node1_mrca_global_mean = (sum(self.dist_node1_mrca_sum_per_reln) * 1.0) / sum(self.edge_weight)
    self.dist_node2_mrca_global_mean = (sum(self.dist_node2_mrca_sum_per_reln) * 1.0) / sum(self.edge_weight)
    
  # this function computes the score metric value associated with individual pair of taxa 
  def _SetCostMetric(self, employ_basic_scoring):      
    # important - sourya
    # to compute the mean distance between this couplet, with respect to input source trees
    if (employ_basic_scoring == False):
      self._CalcMeanDist()   
    # now assign the score metric for individual edge type
    for edge_type in range(4):
      # if basic scoring is employed then no branch length information is accounted for couplet scoring
      if (employ_basic_scoring == True):
	self.Connect_Edge_Cost[edge_type] = self.edge_weight[edge_type] * self.conn_pr_val[edge_type]
      else:
	if (self.edge_weight[edge_type] == 0):
	  self.Connect_Edge_Cost[edge_type] = 0
	else:
	  mult_factor1 = (1.0 * self.dist_node1_mrca_global_mean) / max(math.fabs(self.dist_node1_mrca_global_mean - self.dist_node1_mrca_mean_per_reln[edge_type]), 0.00001)
	  mult_factor2 = (1.0 * self.dist_node2_mrca_global_mean) / max(math.fabs(self.dist_node2_mrca_global_mean - self.dist_node2_mrca_mean_per_reln[edge_type]), 0.00001)
	  self.Connect_Edge_Cost[edge_type] = self.edge_weight[edge_type] * self.conn_pr_val[edge_type] * (mult_factor1 * mult_factor2)
      
  # this function returns the connection priority value for input edge type
  def _GetConnPrVal(self, edge_type):
    return self.conn_pr_val[edge_type]
      
  """ 
  this function calculates connection priority value for each of the edge types, 
  for this particular connection between a pair of nodes in the final tree 
  """
  def _SetConnPrVal(self, single_edge_prior):
    # find the consensus relations 
    maxval = max(self.edge_weight)
    for i in range(4):
      if (self.edge_weight[i] == maxval):
	self.consensus_reln_list.append(i) 
    
    # this is the sum of all the edge type instances (no of occurrences)
    listsum = sum(self.edge_weight)
    # now determine the connection priority of a particular edge type with respect to other edges     
    for edge_type in range(4):
      # here we use the difference of current edge type frequency with the frequencies of all other edge types 
      #self.conn_pr_val[edge_type] = 2 * self.edge_weight[edge_type] - listsum
      self.conn_pr_val[edge_type] = ((self.edge_weight[edge_type] * 1.0) / listsum)
    
    """ 
    this code section is used when there exists NO EDGE relationship between a pair of taxa
    and we want to detect it 
    """
    if (not single_edge_prior):
      """ if there is no vote for any particular edge type other than RELATION_R4,
      (that is, corresponding settings did not occur in any of the source tree)
      then we make only the RELATION_R4 settings as valid - 
      they will only be considered for joining this pair in the final tree """
      if (self.edge_weight[RELATION_R4] != 0)\
	and (self.edge_weight[RELATION_R2] == 0)\
	and (self.edge_weight[RELATION_R1] == 0)\
	and (self.edge_weight[RELATION_R3] == 0):
	return 1
      else:
	return 0
    else:
      outlist = [0, RELATION_R4]
      for edge_type in range(4):
	if (self.edge_weight[edge_type] == listsum) and (listsum > 0):
	  outlist = [1, edge_type]
	  break
	elif (self.edge_weight[edge_type] > 0) and (self.edge_weight[edge_type] < listsum):
	  break
      return outlist

##-----------------------------------------------------
""" this class is representative of a cluster of taxa that are related via equality relationship 
according to the rule of equivalence partition """
class Cluster_node(object):
  def __init__(self, inp_taxa=None):
    # this list contains the species list of the current cluster
    self.Species_List = [] 
    # can be 0 or 1 - denote whether the cluster node has been explored during newick string construction
    self.explored = 0    
    # this variable stores the out edge list from this cluster
    # each list element is the other cluster index
    self.out_edge_list = []
    # this variable stores the in edge list from this cluster
    # each list element is the other cluster index 
    self.in_edge_list = []
    # this variable stores the NO edge list from this cluster
    # each list element is the other cluster index 
    self.no_edge_list = []
    # during initialization, append one tuple to this cluster
    if inp_taxa is not None:
      self._Append_taxa(inp_taxa)    

  # these functions keep track whether the cluster node is used during newick string formation for supertree construction
  # each of the clusters (containing a set of taxa) should be visited exactly once for supertree generation
  def _SetExploredStatus(self):
    self.explored = 1

  def _ResetExploredStatus(self):
    self.explored = 0
    
  def _GetExploredStatus(self):
    return self.explored
          
  # returns the constituent species list of this cluster
  def _GetSpeciesList(self):
    return self.Species_List
        
  # append one species information in this cluster
  def _Append_taxa(self, inp_taxa):
    if inp_taxa not in self.Species_List:
      self.Species_List.append(inp_taxa)
  
  # it returns the final cluster node connectivity (tree shape) -- in edges
  def _Get_Indegree(self):
    return len(self.in_edge_list)

  # it returns the final cluster node connectivity (tree shape) -- out edges
  def _Get_Outdegree(self):
    return len(self.out_edge_list)
      
  # it returns the out edge -- of the cluster node to the other nodes (clique formation)
  def _GetOutEdgeList(self):
    return self.out_edge_list
    
  # it returns the in edge -- of the cluster node to the other nodes (clique formation)
  def _GetInEdgeList(self):
    return self.in_edge_list    

  # it returns the in edge -- of the cluster node to the other nodes (clique formation)
  def _GetNoEdgeList(self):
    return self.no_edge_list    
    
  # it adds one out edge information to both the original connectivity (clique) and the final connectivity (tree shape)
  def _AddOutEdge(self, dest_clust_idx):
    if dest_clust_idx not in self.out_edge_list:
      self.out_edge_list.append(dest_clust_idx)
    
  # it adds one in edge information to both the original connectivity (clique) and the final connectivity (tree shape)
  def _AddInEdge(self, src_clust_idx):
    if src_clust_idx not in self.in_edge_list:
      self.in_edge_list.append(src_clust_idx)

  # it adds one in edge information to both the original connectivity (clique) and the final connectivity (tree shape)
  def _AddNoEdge(self, src_clust_idx):
    if src_clust_idx not in self.no_edge_list:
      self.no_edge_list.append(src_clust_idx)
      
  # here the final connectivity is changed (not the original clique based connectivity) -- out edge remove
  def _RemoveOutEdge(self, dest_clust_idx):
    if dest_clust_idx in self.out_edge_list:
      self.out_edge_list.remove(dest_clust_idx)    
    
  # here the final connectivity is changed (not the original clique based connectivity) -- in edge remove
  def _RemoveInEdge(self, dest_clust_idx):
    if dest_clust_idx in self.in_edge_list:
      self.in_edge_list.remove(dest_clust_idx)    
    
  # here the final connectivity is changed (not the original clique based connectivity) -- no edge remove
  def _RemoveNoEdge(self, dest_clust_idx):
    if dest_clust_idx in self.no_edge_list:
      self.no_edge_list.remove(dest_clust_idx)    
    
  def _PrintClusterInfo(self, key, Output_Text_File):
    fp = open(Output_Text_File, 'a')    
    fp.write('\n cluster key: ' + str(key))
    fp.write('\n species list: ' + str(self.Species_List))
    #print 'its indegree: ', self.indegree
    #print 'its outdegree: ', self.outdegree
    fp.write('\n out edge list: ' + str(self.out_edge_list))
    fp.write('\n in edge list: ' + str(self.in_edge_list))
    fp.close()    
    