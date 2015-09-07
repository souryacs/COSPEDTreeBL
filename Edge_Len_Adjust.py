#!/usr/bin/env python

"""
this file contains functions to adjust the branch length of the derived supertree
the branch lengths are adjusted so as to maintain the additive property 
"""

import Header
from Header import *
import UtilFunc
from UtilFunc import *

#----------------------------------------------------
# new functions used for QP based branch length assignment of the unweighted supertree
#----------------------------------------------------
"""
this function initializes EdgeInfoDict structure 
where individual edges of the final supertree are used as keys
"""
def Initialize_Edge_Dict(Inp_Tree, Output_Text_File):
  idx = 0
  for e in Inp_Tree.postorder_edge_iter():
    """
    key of edge dictionary : individual edges, represented as a pair (tail node, head node)
    tail node signifies the parent node from which the edge starts
    head node signifies the parent node to which the edge ends
    """
    key = (e.tail_node, e.head_node)
    EdgeInfoDict.setdefault(key, idx)
    idx = idx + 1
  
  if (DEBUG_LEVEL >= 2):
    fp = open(Output_Text_File, 'a')
    for key in EdgeInfoDict:
      fp.write('\n curr EdgeInfoDict key: ' + str(key) + ' val: ' + str(EdgeInfoDict[key]))
    fp.close()

#----------------------------------------------------
"""
this function scans input trees
and assigns weights for individual trees according to their constituent taxon set 
"""
def AssignMatrixWeights(Source_Treelist):
  global Matrix_Weight_Val
  
  for tr in range(len(Source_Treelist)):
    supp_taxa_count = 0
    curr_tree_taxa_set = Source_Treelist[tr].infer_taxa()
    number_of_taxa = len(curr_tree_taxa_set)
    for i in range(number_of_taxa - 1):
      t1 = curr_tree_taxa_set[i]
      for j in range(i+1, number_of_taxa):
	t2 = curr_tree_taxa_set[j]
	key1 = (t1.label, t2.label)
	key2 = (t2.label, t1.label)
	if key1 in TaxaPair_Reln_Dict:
	  target_key = key1
	elif key2 in TaxaPair_Reln_Dict:
	  target_key = key2
	else:
	  continue
	if (TaxaPair_Reln_Dict[target_key]._GetNoSupportTrees() >= 2):
	  # we find that both the taxa within this couplet should contribute to the 
	  # supporting taxa count - sourya
	  supp_taxa_count = supp_taxa_count + 2	#1
	  break
    
    # supporting taxa count related weight
    supp_taxa_related_tree_weight = (1.0 / supp_taxa_count)
   
    # now adjust the matrix weight value of the current location
    Matrix_Weight_Val.append(supp_taxa_related_tree_weight)

#----------------------------------------------------
"""
for executing QP solver, this function creates a batch file storing the required command 
in terms of an objective function
"""
def WriteObjectiveFunctionFile(Out_Text_GLS_input_file):
  fp_txt = open(Out_Text_GLS_input_file, 'w')
  
  # first write the number of variables (unknowns) of QP
  fp_txt.write(str(len(EdgeInfoDict)))

  """
  now for each row of the text file
  first write the number of branches for each taxa pair distance
  then mention those branches
  and finally mention the distance matrix entry
  """
  for taxapair_key in TaxaPair_Reln_Dict:
    br_len_list = TaxaPair_Reln_Dict[taxapair_key]._GetBranchArrayIdxList()
    fp_txt.write('\n' + str(len(br_len_list)))
    for j in range(len(br_len_list)):
      fp_txt.write('\t' + str(br_len_list[j]))
    AvgDistMatVal = TaxaPair_Reln_Dict[taxapair_key]._GetAvgDistMatVal()
    fp_txt.write('\t' + str(AvgDistMatVal))
  
  fp_txt.close()
  
  return
  
#----------------------------------------------------
def Form_ObjectiveFunction():
  # this string contains complete expression
  objective_function_string = ''
  
  for taxapair_key in TaxaPair_Reln_Dict:
    # form the variable containing string template
    # ex: (x[0]+x[1]+x[2])
    # y[] array contains the weights on individual source trees
    br_len_list = TaxaPair_Reln_Dict[taxapair_key]._GetBranchArrayIdxList()
    if (DEBUG_LEVEL >= 2):
      print 'current taxa pair key: ', taxapair_key
    # this string contains the branch length indices (variable list)
    # for one particular distance matrix entry, participating branches are accumulated
    variable_collection = ''
    for j in range(len(br_len_list)):
      variable_collection = variable_collection + 'x[' + str(br_len_list[j]) +']'
      if (j < (len(br_len_list) - 1)):
	variable_collection = variable_collection + '+'
    # for this taxa pair, extract the average distance value computed from various source matrices
    AvgDistMatVal = TaxaPair_Reln_Dict[taxapair_key]._GetAvgDistMatVal()
    if (DEBUG_LEVEL >= 2):
	print 'key : ', taxapair_key, ' AvgDistMatVal: ', AvgDistMatVal
    
    # we check whether this taxa pair is supported by at least two source trees
    # it is proved by checking the length of the branch set associated with this taxa pair
    if (objective_function_string != ''):
      objective_function_string = objective_function_string + '+'

    # extend the objective function in a string representation
    lsq_substr = '(' + str(variable_collection) + '-' + str(AvgDistMatVal) + ')**2' 
    objective_function_string = objective_function_string + lsq_substr   
      
  objective_function_string = '(' + objective_function_string + ')'
  
  return objective_function_string

#----------------------------------------------------
""" this function optimizes the edge length values according to the quadratic programming based optimization
it has inputs:
1) x: array of float type which will store the edge lengths
2) sign = 1.0 means minimization of objective function """
def edge_len_optimize_func(x, sign=1.0):
  global objective_function_string
  return sign * eval(objective_function_string)
    
#----------------------------------------------------
# comment - sourya

""" this function processes individual taxa pairs of the source trees
to note down their distance information, as well as the branches present within them
the information is stored in the dictionary TaxaPairBranchVarDict"""
"""
def Initialize_TaxaPairBranches(Inp_Tree):
  for tk in TaxaPair_Reln_Dict:
    t1 = tk[0]
    t2 = tk[1]
    # find the branches corresponding to the node pair in the Inp_Tree
    # determine the nodes in the Inp_Tree corresponding to the taxa t1 and t2
    node_inp_tree_t1 = Inp_Tree.find_node_with_taxon_label(t1)
    node_inp_tree_t2 = Inp_Tree.find_node_with_taxon_label(t2)
    # determine the MRCA of the specified nodes 
    # with respect to the Inp_Tree
    Taxa_Label_List = [t1, t2]
    mrca_inp_tree_t1_t2 = Inp_Tree.mrca(taxon_labels=Taxa_Label_List)
    if (DEBUG_LEVEL >= 2):
      print 'source tree taxa pair: ', Taxa_Label_List
    # add branch length from the first taxa (leaf node) to the MRCA node
    curr_node = node_inp_tree_t1
    parent_node = curr_node.parent_node
    while (curr_node != mrca_inp_tree_t1_t2) and (parent_node is not None):
      # add this branch information to the TaxaPairBranchVarDict
      # key follows the rule from tail node to head node
      key = (parent_node, curr_node)
      TaxaPair_Reln_Dict[tk]._AddBranchArrayIdx(EdgeInfoDict[key])
      if (DEBUG_LEVEL >= 2):
	print 'node: ', curr_node, 'parent node: ', parent_node
	print 'MRCA - added br idx: ', EdgeInfoDict[key]
      curr_node = parent_node
      parent_node = parent_node.parent_node
    # add branch length from the second taxa (leaf node) to the MRCA node
    curr_node = node_inp_tree_t2
    parent_node = curr_node.parent_node
    while (curr_node != mrca_inp_tree_t1_t2) and (parent_node is not None):
      # add this branch information to the TaxaPairBranchVarDict
      # key follows the rule from tail node to head node
      key = (parent_node, curr_node)
      TaxaPair_Reln_Dict[tk]._AddBranchArrayIdx(EdgeInfoDict[key])
      if (DEBUG_LEVEL >= 2):	  
	print 'node: ', curr_node, 'parent node: ', parent_node
	print 'MRCA - added br idx: ', EdgeInfoDict[key]
      curr_node = parent_node
      parent_node = parent_node.parent_node
    
    if (DEBUG_LEVEL >= 2):    
      print 'curr TaxaPair_Reln_Dict key: ', tk, 'branch idx list: ', TaxaPair_Reln_Dict[tk]._GetBranchArrayIdxList()
    
  return
"""

# end comment - sourya
#----------------------------------------------------
""" 
this function takes inputs of two leaf nodes and their MRCA node
branches between the couplet, with respect to their MRCA node in the derived supertree 
are assigned to their class instance
"""
def AddBranchInfo(node1, node2, mrca_node):
  """ 
  if key exists in the taxa pair information
  then assign branch information
  """
  key1 = (node1.taxon.label, node2.taxon.label)
  key2 = (node2.taxon.label, node1.taxon.label)
  if (key1 in TaxaPair_Reln_Dict):
    tk = key1
  elif (key2 in TaxaPair_Reln_Dict):
    tk = key2
  else:
    return
  
  # add branch length from the first taxa (leaf node) to the MRCA node
  curr_node = node1
  parent_node = curr_node.parent_node
  while (curr_node != mrca_node) and (parent_node is not None):
    # add this branch information to the TaxaPair_Reln_Dict
    # key follows the rule from tail node to head node
    key = (parent_node, curr_node)
    TaxaPair_Reln_Dict[tk]._AddBranchArrayIdx(EdgeInfoDict[key])
    if (DEBUG_LEVEL > 2):
      print 'node: ', curr_node, 'parent node: ', parent_node
      print 'MRCA - added br idx: ', EdgeInfoDict[key]
    curr_node = parent_node
    parent_node = parent_node.parent_node
    
  # add branch length from the second taxa (leaf node) to the MRCA node
  curr_node = node2
  parent_node = curr_node.parent_node
  while (curr_node != mrca_node) and (parent_node is not None):
    # add this branch information to the TaxaPair_Reln_Dict
    # key follows the rule from tail node to head node
    key = (parent_node, curr_node)
    TaxaPair_Reln_Dict[tk]._AddBranchArrayIdx(EdgeInfoDict[key])
    if (DEBUG_LEVEL > 2):
      print 'node: ', curr_node, 'parent node: ', parent_node
      print 'MRCA - added br idx: ', EdgeInfoDict[key]
    curr_node = parent_node
    parent_node = parent_node.parent_node
  
  if (DEBUG_LEVEL > 2):    
    print 'curr TaxaPair_Reln_Dict key: ', tk, 'branch idx list: ', TaxaPair_Reln_Dict[tk]._GetBranchArrayIdxList()
  
  return

#----------------------------------------------------
""" 
for individual couplets, this function assigns the branches between them
with respect to the output supertree (currently unweighted) 
"""
def Initialize_TaxaPairBranches(Inp_Tree):
  # traverse the internal nodes of the tree in postorder fashion
  for curr_node in Inp_Tree.postorder_internal_node_iter():
    # list the leaf and internal children of the current node
    curr_node_child_leaf_nodes = []
    curr_node_child_internal_nodes = []
    for x in curr_node.child_nodes():
      if (x.is_leaf() == True):
	curr_node_child_leaf_nodes.append(x)
      else:
	curr_node_child_internal_nodes.append(x)
    
    # pair of leaf nodes will be first processed
    if (len(curr_node_child_leaf_nodes) > 1):
      for i in range(len(curr_node_child_leaf_nodes) - 1):
	for j in range(i+1, len(curr_node_child_leaf_nodes)):
	  AddBranchInfo(curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j], curr_node)
  
    # one leaf node (direct descendant) and another leaf node (under one internal node)
    # will be related by ancestor / descendant relations
    if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
      for p in curr_node_child_leaf_nodes:
	for q in curr_node_child_internal_nodes:
	  for r in q.leaf_nodes():
	    AddBranchInfo(p, r, curr_node)

    # finally a pair of leaf nodes which are descendant of internal nodes will be related by NO_EDGE relation
    if (len(curr_node_child_internal_nodes) > 1):
      for i in range(len(curr_node_child_internal_nodes) - 1):
	for j in range(i+1, len(curr_node_child_internal_nodes)):
	  for p in curr_node_child_internal_nodes[i].leaf_nodes():
	    for q in curr_node_child_internal_nodes[j].leaf_nodes():
	      AddBranchInfo(p, q, curr_node)

#----------------------------------------------------
# this function assigns branch length information on the Inp_Tree     
def AssignBranchLen(Inp_Tree, Source_Treelist, Output_Text_File, NO_OF_LOOPS, METHOD_OF_QP, QP_Executable):
  """
  this is the objective function represented as a string format
  that need to be passed in QP optimization function
  """
  global objective_function_string  

  """
  this is the weight matrix of input phylogenetic trees
  """
  global Matrix_Weight_Val
  
  # we note the timing for branch length assignment
  start_timestamp = time.time()

  """
  function to assign individual edges of the derived unweighted supertree
  into a dictionary, 
  so that its values and attributes can be easily accessed and modified
  """
  Initialize_Edge_Dict(Inp_Tree, Output_Text_File)
  
  """
  here we process individual couplets of the output supertree
  and assign the branch indices within them 
  """
  Initialize_TaxaPairBranches(Inp_Tree)
  
  """
  assign weights of individual phylogenetic trees
  """
  AssignMatrixWeights(Source_Treelist)
  
  fp1 = open(Output_Text_File, 'a')
  for i in range(len(Matrix_Weight_Val)):
    fp1.write('\n Input tree index: ' + str(i) + ' tree weight: ' + str(Matrix_Weight_Val[i]) + \
      '  no of support taxa: ' + str(1.0 / Matrix_Weight_Val[i]) + \
      '  no of input taxa: ' + str(len(Source_Treelist[i].infer_taxa())))
  fp1.close()

  if (METHOD_OF_QP == 1) or (METHOD_OF_QP == 2):
    # these are the variables for optimization constraints
    # len(EdgeInfoDict) variables are terfor the branch length values
    Edge_Values = numpy.zeros((len(EdgeInfoDict)), dtype=numpy.float)        
    # this function generates the objective function for solving with quadratic programming
    objective_function_string = Form_ObjectiveFunction()
    if (DEBUG_LEVEL >= 2):
      print 'objective_function_string: ', objective_function_string        
    # this sequence of tuples denote the bound of individual variables
    # currently only constraint of the variables is that they must be positive
    bound_tuple = []
    for i in range(len(EdgeInfoDict)):
      t = (0, None)
      bound_tuple.append(t)    
    if (METHOD_OF_QP == 1):
      # method L-BFGS-B 
      res = minimize(edge_len_optimize_func, x0=Edge_Values, args=(1.0,), method='L-BFGS-B', \
	bounds=bound_tuple, options={'disp': True, 'maxiter': NO_OF_LOOPS})
    else:
      # method SLSQP
      res = minimize(edge_len_optimize_func, x0=Edge_Values, args=(1.0,), method='SLSQP', \
	bounds=bound_tuple, options={'disp': True, 'ftol': 0.0001, 'maxiter': NO_OF_LOOPS})
  
    #-------------------------------------------------------
    ## this function generates the constraints used for such quadratic programming 
    #cons = [{'type':'eq','fun':lambda x: (numpy.sum(x[len(EdgeInfoDict):]) - num_of_input_trees)}]
    #for i in range(len(EdgeInfoDict)+num_of_input_trees):
      #cons.append({'type':'ineq','fun':lambda x: (numpy.array(x[i]))})  
    ##for i in range(num_of_input_trees):
      ##cons.append({'type':'ineq','fun':lambda x: (numpy.array(y[i]))})
    
    #if (DEBUG_LEVEL >= 2):
      #print 'constraints: ', cons    
	
    ## now perform QP based optimization of the edge length values, 
    ## using the constraints and objective function representation

    ## method COBYLA
    ##res = minimize(edge_len_optimize_func, x0=Edge_Values, args=(1.0,), method='COBYLA', jac=deriv_edge_len_opt_func, bounds=bound_tuple, options={'disp': True})	#, 'maxiter': 15})  
    
    ## method L-BFGS-B 
    ##res = minimize(edge_len_optimize_func, x0=Edge_Values, args=(Matrix_Wt_Multiplier, 1.0,), method='L-BFGS-B', constraints=cons, options={'disp': True, 'maxiter': 15})
    #res = minimize(edge_len_optimize_func, x0=Edge_Values, args=(1.0,), method='L-BFGS-B', bounds=bound_tuple, options={'disp': True, 'maxiter': 15})
    
    ## method SLSQP
    #res = minimize(edge_len_optimize_func, x0=Edge_Values, args=(1.0,), method='SLSQP', constraints=cons, options={'disp': True, 'ftol': 0.0001, 'maxiter': 15})	#, 'eps' : 1e0})
    #res = minimize(edge_len_optimize_func, x0=Edge_Values, args=(1.0,), method='SLSQP', bounds=bound_tuple, options={'disp': True, 'ftol': 0.0001, 'maxiter': 15})	#, 'eps' : 1e0})
    #-------------------------------------------------------
    
    # printing the result and its attributes
    fp = open(Output_Text_File, 'a')
    fp.write('\n success of solution: ' + str(res.success))
    fp.write('\n status of solution: ' + str(res.status))
    fp.write('\n corresponding message: ' + str(res.message))
    #fp.write('\n no of iterations used: ' + str(res.nit))
    fp.write('\n solution - array: ' + str(res.x)) 
    fp.close()
  
    # assign the edge length values to the tree edges
    for e in Inp_Tree.postorder_edge_iter():
      # key of the dictionary 
      # tail node signifies the parent node from which the edge starts
      # head node signifies the parent node to which the edge ends
      key = (e.tail_node, e.head_node)
      e.length = res.x[EdgeInfoDict[key]]

  else:
    print '*** now starting GLS based QP optimization of the branch length values ***'
    k = Output_Text_File.rfind("/")
    Out_Text_GLS_input_file = Output_Text_File[:(k+1)] + 'GLS_input.txt'
    Out_Text_GLS_output_file = Output_Text_File[:(k+1)] + 'GLS_output.txt'
    WriteObjectiveFunctionFile(Out_Text_GLS_input_file)
    
    # call the C executable to generate QP outcome
    sys_command_str = QP_Executable + str(' ') + Out_Text_GLS_input_file + ' ' + Out_Text_GLS_output_file
    os.system(sys_command_str)
    
    # now read the edge contents from the file
    edge_value_list = []
    fp = open(Out_Text_GLS_output_file, 'r')
    for line in fp:
      if (line != ''):
	edge_value_list.append(float(line))
    fp.close()
    
    # assign the edge length values to the tree edges
    for e in Inp_Tree.postorder_edge_iter():
      # key of the dictionary 
      # tail node signifies the parent node from which the edge starts
      # head node signifies the parent node to which the edge ends
      key = (e.tail_node, e.head_node)
      # if by chance, any edge value becomes negative, this correction makes it positive
      e.length = edge_value_list[EdgeInfoDict[key]]
      #e.length = math.fabs(edge_value_list[EdgeInfoDict[key]])
      
      
  ##------------------------------------------
  print '*** end of branch length assignment ***'
  
  end_timestamp = time.time()
  fp = open(Output_Text_File, 'a')
  fp.write('\n Branch length assignment -- total time required : '+ str(end_timestamp - start_timestamp))
  fp.close()
  
  return
  
  