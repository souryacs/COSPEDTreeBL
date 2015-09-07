#!/usr/bin/env python


##---------------------------------------------
''' 
this program is used to generate a supertree (consensus) from a set of constituent trees
input trees contain the branch length information as well
output needs to conform to the input source trees
so as to get the minimal least square error with respect to the source trees
the input is multiple source trees
there may be conflicts among the input tree - we have to select the consensus

Author: Sourya Bhattacharyya
Dept of CSE, IIT Kharagpur
V1.0 - 07.04.2014 - basic code
V2.0 - 01.07.2014 - modified the clustering and conflict detection routine
V3.0 - 23.03.2015 - reduced storage complexity and cleaned the code
V4.0 - 06.09.2015 - comment, code clean, reduced runtime, improved QP solver. 
		  - removed QP solvers other than BFGS.
		  - used QP execution path.
'''

## Copyright 2014, 2015, Sourya Bhattacharyya and Jayanta Mukherjee.
## All rights reserved.
##
## See "LICENSE.txt" for terms and conditions of usage.
##
##---------------------------------------------

import Header
from Header import *
import Cost_Update
from Cost_Update import *
import Process_Queues
from Process_Queues import *
import ReachGraph_Update
from ReachGraph_Update import *
import UtilFunc
from UtilFunc import *
import Edge_Len_Adjust
from Edge_Len_Adjust import *
import RefineTree
from RefineTree import *

##-----------------------------------------------------
# this function is useful to parse various options for input data processing
def parse_options():  
  parser = OptionParser()
  
  parser.add_option("-I", "--INPFILE", \
			  type="string", \
			  action="store", \
			  dest="INP_FILENAME", \
			  default="", \
			  help="name of the file containing input phylogenetic trees")
  
  parser.add_option("-p", "--inpform", \
			  type="int", \
			  action="store", \
			  dest="inp_file_format", \
			  default=1, \
			  help="1 - format of the input file (containing the input treelist) is NEWICK (default) \
			  2 - input file format is NEXUS")

  parser.add_option("-s", "--basicscore", \
			  action="store_true", \
			  dest="basic_score", \
			  default=False, \
			  help="if TRUE, we employ couplet support scores as that of COSPEDTree. \
			  There, branch length values of input trees are not accounted for couplet scores. \
			  Default FALSE.")  
    			  
  parser.add_option("-b", "--binary", \
			  action="store_false", \
			  dest="binary_suptr", \
			  default=True, \
			  help="if TRUE, it produces a strictly binary supertree. Otherwise, the tree can be non-binary. Default TRUE.")
    			  
  parser.add_option("-d", "--dfsref", \
			  action="store_true", \
			  dest="dfs_parent_refine", \
			  default=True, \
			  help="if TRUE, Multiple parent problem (C2) is tackled by deterministic parent selection. \
			  Otherwise, arbitrary parent assignment is performed. Default is TRUE.")  
    			    
  parser.add_option("-m", "--method", \
			  type="int", \
			  action="store", \
			  dest="method_of_QP", \
			  default=3, \
			  help="1 - Employes L-BFGS-B method for least square based branch length assignment \
			  2 - Employes SLSQP method for least square based branch length assignment \
			  3 - C method on GNU library (BFGS) for least square based branch length assignment (default)")
  
  parser.add_option("-l", "--loops", \
			  type="int", \
			  action="store", \
			  dest="no_of_loops", \
			  default=15, \
			  help="No of iterations (loops) employed during the branch length optimization\
			  only applicable if the method of nonlinear programming is either 1 or 2 \
			  generally more iterations yield better optimization \
			  default number of iterations is 15")
    
  parser.add_option("-Q", "--QPExec", \
			  type="string", \
			  action="store", \
			  dest="QP_Exec_Path", \
			  default="", \
			  help="Absolute path of the executable for QP solver")

  parser.add_option("-w", "--weighttaxa", \
			  action="store_false", \
			  dest="weight_taxa_subset", \
			  default=True, \
			  help="this is a boolean flag option \
				using this option toggles the existing configuration (Default TRUE) \
				if TRUE, then this option weighs couplet statistics according \
				to the size of taxa subset underlying MRCA of that couplet")  
    
  opts, args = parser.parse_args()
  return opts, args
  
  
##-----------------------------------------------------
# main function
def main():  
  opts, args = parse_options()
  
  ROOTED_TREE = True	#opts.default_rooted
  PRESERVE_UNDERSCORE = True	#opts.preserve_underscores
  if (opts.inp_file_format == 1):
    INPUT_FILE_FORMAT = 'newick'
  else:
    INPUT_FILE_FORMAT = 'nexus'
  INPUT_FILENAME = opts.INP_FILENAME
  NO_OF_QUEUES = 1	# sourya - default set
  DFS_PARENT_REFINE = opts.dfs_parent_refine
  NO_OF_LOOPS = opts.no_of_loops
  METHOD_OF_QP = opts.method_of_QP
  EMPLOY_BASIC_SCORING_COSPEDTREE = opts.basic_score
  BINARY_SUPERTREE_OPTION = opts.binary_suptr
  WEIGHT_TAXA_SUBSET = opts.weight_taxa_subset
  
  """
  abspath function converts input possibly relative path into an absoloute path
  """
  if (opts.QP_Exec_Path == ""):
    print '******** THERE IS NO PATH FOR QP SOLVER (GNU_BFGS2) IS PROVIDED - RETURN **********'
    return
  QP_EXEC_PATH = os.path.abspath(opts.QP_Exec_Path)  
  
  if (INPUT_FILENAME == ""):
    print '******** THERE IS NO INPUT FILE SPECIFIED - RETURN **********'
    return
  else:
    print 'input filename: ', INPUT_FILENAME
    
  # according to the location of input filename
  # adjust the locations of the output files as well
  k = INPUT_FILENAME.rfind("/")
  if (k == -1):
    dir_of_inp_file = './'
  else:
    dir_of_inp_file = INPUT_FILENAME[:(k+1)]
  if (DEBUG_LEVEL > 1):
    print 'dir_of_inp_file: ', dir_of_inp_file
  
  dir_of_curr_exec = dir_of_inp_file + 'COSPEDBL'
  if (DFS_PARENT_REFINE == True):
    dir_of_curr_exec = dir_of_curr_exec + '_D'    
  if (BINARY_SUPERTREE_OPTION == True):
    dir_of_curr_exec = dir_of_curr_exec + '_B'    
  if (EMPLOY_BASIC_SCORING_COSPEDTREE == True):
    dir_of_curr_exec = dir_of_curr_exec + '_BS'
  if (WEIGHT_TAXA_SUBSET == True):
    dir_of_curr_exec = dir_of_curr_exec + '_W'
  dir_of_curr_exec = dir_of_curr_exec + '_QP_' + str(METHOD_OF_QP)
  if (METHOD_OF_QP != 3):
    dir_of_curr_exec = dir_of_curr_exec + '_Loop_' + str(NO_OF_LOOPS)

  # create the directory
  if (os.path.isdir(dir_of_curr_exec) == False):
    mkdr_cmd = 'mkdir ' + dir_of_curr_exec
    os.system(mkdr_cmd)         
    
  # append the current output directory in the text file
  Output_Text_File = dir_of_curr_exec + '/' + 'Complete_Output_Description.txt'
    
  fp = open(Output_Text_File, 'w')    
        
  # note the program beginning time 
  start_timestamp = time.time()
    
  #-------------------------------------  
  """ 
  read the input treelist file
  """
  Input_Treelist = Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME)  
  
  #-------------------------------------    
  """
  from the input trees, note the number of taxa (total)
  and insert in the complete taxa list
  """
  for tr_idx in range(len(Input_Treelist)):
    taxa_labels_curr_tree = Input_Treelist[tr_idx].infer_taxa().labels()
    if (DEBUG_LEVEL > 1):
      fp.write('\n Tree no : ' + str(tr_idx+1) +  'no of leaf nodes: ' + str(len(taxa_labels_curr_tree)))
    if (DEBUG_LEVEL > 2):
      fp.write('\n taxa set belonging to current tree: ' + str(taxa_labels_curr_tree))
    for i in range(len(taxa_labels_curr_tree)):
      if taxa_labels_curr_tree[i] not in COMPLETE_INPUT_TAXA_LIST:
	COMPLETE_INPUT_TAXA_LIST.append(taxa_labels_curr_tree[i])
  
  number_of_taxa = len(COMPLETE_INPUT_TAXA_LIST)
  
  """ 
  we also define one structure "Taxa_Info_Dict" marked by a taxa
  individual taxon information is entered in this dictionary
  """
  for label in COMPLETE_INPUT_TAXA_LIST:
    Taxa_Info_Dict.setdefault(label, Single_Taxa())
  
  #---------------------------
  """
  if variable weight of individual relation frequency is considered
  then we process MRCA nodes of individual couplets
  """
  if (WEIGHT_TAXA_SUBSET == True):
    for tr_idx in range(len(Input_Treelist)):
      FindCoupletUnderlyingTaxon(Input_Treelist[tr_idx])
  #---------------------------
  """ 
  now process individual trees to find the couplet relations of those trees
  """
  for tr_idx in range(len(Input_Treelist)):
    DeriveCoupletRelations(Input_Treelist[tr_idx], tr_idx, WEIGHT_TAXA_SUBSET)
 
  # printing couplet information
  if (DEBUG_LEVEL >= 0):
    fp.write('\n  total no of taxa: ' + str(number_of_taxa))
  if (DEBUG_LEVEL > 1):
    fp.write('\n len Taxa_Info_Dict: ' + str(len(Taxa_Info_Dict)))
    fp.write('\n len COMPLETE_INPUT_TAXA_LIST: ' + str(COMPLETE_INPUT_TAXA_LIST))
    fp.write('\n len TaxaPair_Reln_Dict : ' + str(len(TaxaPair_Reln_Dict)))
  
  # close the output text file
  fp.close()
  
  #------------------------------------------------------------
  """ 
  we also calculate the connection value between each pair of nodes in the output tree
  the value defines the majority of the edge type that is between those 2 nodes 
  """
  for l in TaxaPair_Reln_Dict:
    """
    calculate the consensus and priority measures associated with each couplet for different relations
    single_edge_exist means that the couplet is non-conflicting
    only one type of relation exists between them in the input trees
    in such a case, include only that relation in the priority queue
    """
    single_edge_exist_list = TaxaPair_Reln_Dict[l]._SetConnPrVal(True)
    single_edge_exist = single_edge_exist_list[0]
    edge_type = single_edge_exist_list[1]
    """ 
    we calculate the support score value between individual couplets and for each relations
    previously the support score value was equal to the priority metric
    now we change it to make it a product of the frequency and the priority measures 
    """
    TaxaPair_Reln_Dict[l]._SetCostMetric(EMPLOY_BASIC_SCORING_COSPEDTREE)
  
    #------------------------------------------------------------
    """ also update the cost value for individual elements in the list 
    "Cost_List_Taxa_Pair_Multi_Reln or Cost_List_Taxa_Pair_Single_Reln"
    each list element contains the following values:
    1) taxa1 and taxa2    
    3) edge type (one at a time - so there will be 4 entries for each node pair)
    4) edge cost (the cost associated with one particular edge type """
    
    #-------------------
    # add - sourya - if binary supertree option is strictly maintained 
    # then we add all the conflicting relations in the source trees
    # for taxa pairs supported by only one source tree, if -a option is false
    # then only one relation is included
    # otherwise all the relations are included
    #-------------------
    
    if (single_edge_exist == 0):
      for edge_type in range(4):
	# now we add only those relations between the taxa pair which are supported by at least one source tree
	if (TaxaPair_Reln_Dict[l]._GetEdgeWeight(edge_type) > 0):
	  sublist = [l[0], l[1], edge_type, TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(edge_type)]
	  Cost_List_Taxa_Pair_Multi_Reln.append(sublist)
    else:
      # for single edge type, assign that particular connection
      sublist = [l[0], l[1], edge_type, TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(edge_type)]
      """ this connection is only possible between the current examined taxa pairs
      if we maintain different queues for storing conflicting and non-conflicting taxa pairs
      then the information is placed in Cost_List_Taxa_Pair_Single_Reln
      otherwise it is placed in Cost_List_Taxa_Pair_Multi_Reln """
      if (NO_OF_QUEUES == 2):
	Cost_List_Taxa_Pair_Single_Reln.append(sublist)
      else:
	Cost_List_Taxa_Pair_Multi_Reln.append(sublist)
  
  # we print the original connection status for all the tree nodes
  if (DEBUG_LEVEL > 2):
    for l in TaxaPair_Reln_Dict:
      #print 'printing info for the TaxaPair_Reln_Dict key: ', l
      TaxaPair_Reln_Dict[l]._PrintRelnInfo(l, Output_Text_File)  
  
  # note the timestamp
  # this will signify the time required for tree reading and couplet feature extraction
  data_read_timestamp = time.time()  
    
  #------------------------------------------------------------    
  """ here we allocate the list of clusters
  initially all the clusters contain one taxa
  each of the cluster has the index of the corresponding taxa in the COMPLETE_INPUT_TAXA_LIST """
  for i in range(len(COMPLETE_INPUT_TAXA_LIST)):
    Create_Cluster_Taxa_Label(i, COMPLETE_INPUT_TAXA_LIST[i])
  
  #------------------------------------------------------------
  """ we initialize the Reachability_Graph_Mat
  for all the clusters of nodes possible, this indicate the edges between a cluster pair
  initially the number of clusters are thought of equal as the number of taxa
  gradually, as the taxa subsets are merged in a cluster, corresponding count is decreased 
  convention -  out edge from the first node to the second node 
  this is a numpy 2D array 
  values Mat[x][y] = 1 means x->y
  Mat[x][y] = Mat[y][x] = 2 means x and y are connected via NO EDGE """
  fp = open(Output_Text_File, 'a')
  
  Reachability_Graph_Mat = numpy.zeros((number_of_taxa, number_of_taxa), dtype=numpy.int)
  if (DEBUG_LEVEL > 0):
    fp.write('\n shape of Reachability_Graph_Mat: ' + str(numpy.shape(Reachability_Graph_Mat)))
  
  # this information is printed to know the maximum possible iterations that the while loops will undergo
  if (DEBUG_LEVEL > 1):
    fp.write('\n =========== max connection pair ============= : ' + str((len(Cost_List_Taxa_Pair_Single_Reln) + len(Cost_List_Taxa_Pair_Multi_Reln))))      
    fp.write('\n len Cost_List_Taxa_Pair_Single_Reln: ' + str(len(Cost_List_Taxa_Pair_Single_Reln)))
    fp.write('\n len Cost_List_Taxa_Pair_Multi_Reln : ' + str(len(Cost_List_Taxa_Pair_Multi_Reln)))
  
  fp.close()
  #------------------------------------------------------------
  """ now we have to sort the Cost_List_Node_Pair according to the edge cost value 
  that is the 4th field of the sublist 
  we use custom sorting operation """
  Sort_Cost_List_Initial(Cost_List_Taxa_Pair_Multi_Reln)
  
  # if there is provision to include single connectivity edges then we sort that list also
  if (NO_OF_QUEUES == 2):
    Sort_Cost_List_Initial(Cost_List_Taxa_Pair_Single_Reln)
  
  # note the timestamp  
  data_initialize_timestamp = time.time()
  
  #------------------------------------------------------------
  # print the queue storing the scores of individual relations
  if (DEBUG_LEVEL > 2):
    if (NO_OF_QUEUES == 2):
      fp = open(Output_Text_File, 'a')
      fp.write('\n printing contents for the non-conflicting queue (couplet relation score)')
      fp.close()
      PrintQueueInfo(Cost_List_Taxa_Pair_Single_Reln, Output_Text_File)
      
    fp = open(Output_Text_File, 'a')
    fp.write('\n printing contents for the conflicting queue (couplet relation score)')
    fp.close()
    PrintQueueInfo(Cost_List_Taxa_Pair_Multi_Reln, Output_Text_File)
  
  #------------------------------------------------------------
  """ if we have stored taxa pairs depicting single relation instance 
  then we first process that corresponding queue """
  if (NO_OF_QUEUES == 2):
    Reachability_Graph_Mat = Proc_Queue(Reachability_Graph_Mat, 1, Output_Text_File)
  """ then we process the queue containing the taxa pairs depicting multiple relation instance """
  Reachability_Graph_Mat = Proc_Queue(Reachability_Graph_Mat, 0, Output_Text_File)
  
  # we print the final connection status for all the tree nodes
  if (DEBUG_LEVEL > 2):
    for l in Taxa_Info_Dict:
      #print 'printing information for the Taxa ', l
      Taxa_Info_Dict[l]._PrintFinalTaxaInfo(l, Output_Text_File) 
  
  # print the cluster information 
  if (DEBUG_LEVEL > 0):
    fp = open(Output_Text_File, 'a')
    fp.write('\n **** total number of clusters: ' + str(len(CURRENT_CLUST_IDX_LIST)))
    fp.write('\n CURRENT_CLUST_IDX_LIST contents: ')
    fp.write(str(CURRENT_CLUST_IDX_LIST))
    fp.write('\n ========== cluster information after reachability graph generation =============')
    fp.close()
    for i in Cluster_Info_Dict:
      #print 'printing the information for cluster node: ', i
      Cluster_Info_Dict[i]._PrintClusterInfo(i, Output_Text_File)
  
  # note the timestamp
  couplet_procssing_timestamp = time.time()    
  #------------------------------------------------------------
  """ now perform the transitive reduction of the closure formed by connection of the cluster of nodes in the above operation
  this is required to handle the following scenario:
  suppose, there exists a case such that A->C, B->C and A->B
  then in the final graph, only A->B and B->C information needs to be preserved
  in order to form the DAG """
  CompressDirectedGraph(Reachability_Graph_Mat)
    
  #------------------------------------------------------------
  # now instead of arbitrary assignment of the parent node for individual clusters 
  # we assign parent node according to the source tree relationships
  # this will solve the multiple parent problem C2 as discussed in the manuscript 
  # this is a new addition and marked under the DFS based parent refinement option 
  if (DFS_PARENT_REFINE == True):
    SolveMultipleParentC2Problem(Output_Text_File)    
    
  # print the cluster information 
  if (DEBUG_LEVEL > 0):
    fp = open(Output_Text_File, 'a')
    fp.write('\n **** total number of clusters: ' + str(len(CURRENT_CLUST_IDX_LIST)))
    fp.write('\n CURRENT_CLUST_IDX_LIST contents: ')
    fp.write(str(CURRENT_CLUST_IDX_LIST))    
    fp.write('\n ========== cluster information after transitive reduction =============')
    fp.close()
    for i in Cluster_Info_Dict:
      #print 'printing the information for cluster node: ', i
      Cluster_Info_Dict[i]._PrintClusterInfo(i, Output_Text_File)
  #----------------------------------------------
  """
  now this section constructs the supertree from the generated DAG 
  this is performed by repeatedly extracting the nodes with minimum indegree
  basically we first form a string which represents the supertree 
  """
  no_of_components = 0	# for forest
  while (1):
    root_clust_node_idx = Extract_Node_Min_Indeg(len(CURRENT_CLUST_IDX_LIST))
    if (root_clust_node_idx == -1):
      break
    Tree_Str = PrintNewick(root_clust_node_idx)
    no_of_components = no_of_components + 1
    if (no_of_components == 1):	# first component
      Final_Supertree_Str = Tree_Str
    else:
      Final_Supertree_Str = Final_Supertree_Str + ',' + Tree_Str
  
  # with the final tree string, finally generate the tree result 
  Final_Supertree_Str = '(' + Final_Supertree_Str + ')'

  fp = open(Output_Text_File, 'a')
  fp.write('\n --- original supertree as newick string --- ' + Final_Supertree_Str) 
  
  Final_Supertree_Str = Remove_Extra_Paranthesis(Final_Supertree_Str)  
  
  fp.write('\n --- after removing extra paranthesis -- supertree as newick string --- ' + Final_Supertree_Str) 
  fp.close()
  
  # read the supertree (without branch length information)
  Final_Supertree = dendropy.Tree.get_from_string(Final_Supertree_Str, schema="newick")	#, preserve_underscores=PRESERVE_UNDERSCORE, default_as_rooted=True)
    
  # add - sourya  
  if (BINARY_SUPERTREE_OPTION == True):
    # this function removes all multifurcating clusters and produces binary tree (except problem C3)
    Refine_Supertree_Binary_Form(Final_Supertree, Output_Text_File)
    fp = open(Output_Text_File, 'a')
    fp.write('\n --- after binary refinement --- output tree without branch length (in newick format): ' + Final_Supertree.as_newick_string())    
    fp.close()
  else:
    fp = open(Output_Text_File, 'a')
    fp.write('\n --- user did not provide option for producing strict binary supertree - so output tree can be non-binary')
    fp.close()
      
  # final timestamp
  Multiple_parent_solution_timestamp = time.time()    
  #-------------------------------------------------------------
  # this function assigns the branch length information on the copied tree
  AssignBranchLen(Final_Supertree, Input_Treelist, Output_Text_File, NO_OF_LOOPS, METHOD_OF_QP, QP_EXEC_PATH)
  
  # note the timestamp for counting the time elapsed in adjusting the branch length
  branch_len_timestamp = time.time() 
  #-------------------------------------------------------------
      
  # write this tree on a separate text file
  out_treefilename = dir_of_curr_exec + '/' + 'COSPEDBL_GLS_outtree_Newick.tre'  
  Final_Supertree.write_to_path(out_treefilename, 'newick')
  
  fp = open(Output_Text_File, 'a')
  fp.write('\n\n ---output tree with branch length information (in newick format): ' + Final_Supertree.as_newick_string())
  fp.close()
  
  fp = open(Output_Text_File, 'a')
  fp.write('\n \n\n ===============>>>>>>>>>>>>>>> TIME COMPLEXITY OF THE METHOD (in seconds)')
  
  fp.write('\n \n reading the data: ' + str(data_read_timestamp - start_timestamp) + \
	'\n initialization of the structure: ' + str(data_initialize_timestamp - data_read_timestamp) + \
	'\n processing all couplets and support scores: ' + \
	      str(couplet_procssing_timestamp - data_initialize_timestamp) + \
	'\n multiple parent problem + binary tree: ' + \
	      str(Multiple_parent_solution_timestamp - couplet_procssing_timestamp) + \
	    '\n adjusting branch length of the tree: ' + \
	      str(branch_len_timestamp - Multiple_parent_solution_timestamp))
  fp.write('\n \n Total time taken (in seconds) : ' + str(branch_len_timestamp - start_timestamp))
  fp.close()
    
  #--------------------------------------------------------------  
  # delete the storage variables associated with the current execution 
  
  # clear the dictionaries
  Cluster_Info_Dict.clear()
  Taxa_Info_Dict.clear()
  TaxaPair_Reln_Dict.clear()
  EdgeInfoDict.clear()
  
  # clear the lists associated
  if (len(Cost_List_Taxa_Pair_Multi_Reln) > 0):
    Cost_List_Taxa_Pair_Multi_Reln[:] = []
  if (len(Cost_List_Taxa_Pair_Single_Reln) > 0):
    Cost_List_Taxa_Pair_Single_Reln[:] = []
  if (len(COMPLETE_INPUT_TAXA_LIST) > 0):
    COMPLETE_INPUT_TAXA_LIST[:] = []
  if (len(CURRENT_CLUST_IDX_LIST) > 0):
    CURRENT_CLUST_IDX_LIST[:] = []
  if (len(Matrix_Weight_Val) > 0):
    Matrix_Weight_Val[:] = []
  
  # free the reachability graph (numpy array)
  del Reachability_Graph_Mat
    
#-----------------------------------------------------

if __name__ == "__main__":
    main() 
  
