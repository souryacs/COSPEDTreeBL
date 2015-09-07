# COSPEDTreeBL
COSPEDTree (Couplet based supertree) with branch length assignment, to produce weighted supertree.

Description
-----------------

COSPEDTreeBL is an extended version of COSPEDTree, which produces "weighted" supertrees from 
input phylogenetic trees. 

Overview of COSPEDTree
-----------------------------------

In the basic COSPEDTree algorithm, input phylogenetic trees may contain overlapping taxa set. 
These trees often exhibit different topologies among constituent taxa set. The objective is to produce a
 supertree covering all the input taxa so that individual taxa subsets exhibit consensus relationships 
as much as possible. This is known as satisfying "Maximum agreement property".

However, given the conflicting nature of input trees, often the consensus (most freqent) relations 
among individual taxa subsets may not be reflected in the final supertree. This is because, 
consensus relation among a taxa subset may conflict with the consensus relation of another taxa subset.

In the basic COSPEDTree algorithm, supertree computation is formed using a greedy approach. 
The method is based on partitioning the input set of taxa based on equivalence relation. 
The relation is defined between individual taxa pairs (couplet), leading to the proposed couplet 
based supertree technique.

The relationship among a couplet can be of the following three types: 
1) ancestor / descendent 2) sibling 3) inconclusive relation characteristics. The sibling 
relation is in fact, an equivalence relation. Definition of these relations can be found in our paper 
(reference provided below). According to the equivalence relation, input taxa set is partitioned. After that, 
a Directed Acyclic Graph (DAG) is formed initially. From this DAG, the output tree is generated.

Input source trees can be either in NEWICK format or in NEXUS format. However, all the source trees 
should have identical input formats. Output tree is generated in the NEWICK format.

Improvements in COSPEDTreeBL
------------------------------------------

COSPEDTreeBL (COSPEDTree with Branch Length) is a program to assign the branch length of the 
supertree produced by COSPEDTree. Branch length values are obtained from 
a set of input weighted phylogenetic trees (carrying distance information).

Input
---------

A set of input phylogenetic trees with overlapping taxa set and branch length (distance) information.

Output
----------

A weighted supertree whose topology is derived from COSPEDTree algorithm (or its slight modification accounting 
input branch length information). Branch lengths of the generated supertree are assigned 
based on the weights of source trees. The branch length assignment is carried out in such a fashion that 
the difference between the output distance matrix (generated from the derived weighted supertree) and the 
input distance matrices (generated from the input phylogenetic trees) is as small as possible.

Validation
------------

For branch length prediction accuracy, use least square error between the distance matrices and the output supertree
for topological accuracy, use RF distance between the source trees and the output supertree.

Dependencies / Installation Requirements
-----------------------------------------------------

COSPEDTreeBL is developed in Linux Systems (Ubuntu 14.04), using Python 2.7.

User needs to install following before using this package:

1) Python 2.7 (available in Ubuntu, by default)

Note: We have not tested the code on Python 3. Any user having Python 3 environment need to 
check the correct execution of our code, and optionally needs to upgrade it accordingly.

********* We plan to support Python 3 environment in some future release.

2) Dendropy 3.12.0 ( available on the link: https://pythonhosted.org/DendroPy/ )

Note: there is a new release of Dendropy 4.0 but we have used 3.12.0 for the implementation. 
We did not upgrade the code for Dendropy 4.0 support, so any user having this new version of 
Dendropy might need to check the functionalities of COSPEDTreeBL and possibly 
upgrade / replace / edit few dendropy related functions. So, we recommend users to use the 
earlier version of Dendropy, to avoid any conflict.

*********** Support for Dendropy 4 and corresponding update of code will be done in a future release.

3) A binary executable file GNU_BFGS2 is provided (in a zipped archieve) along with this release. 
User needs to Download, extract the archieve and save it in the location containing the source codes.

The package requires GNU scinetific library (GSL) for its execution. The package can be installed in any of the following ways:

A) User needs to go to the link http://www.gnu.org/software/gsl/ and follow the 
download and installation instructions.

B) Otherwise, for UBUNTU 14.04 system users, please check the following link: 
http://askubuntu.com/questions/490465/install-gnu-scientific-library-gsl-on-ubuntu-14-04-via-terminal

C) Or you can use the following command for CentOS or Fedora

yum install gsl-devel

Issues (Mentioned in the following (A) and (B) points)
----------------------------------------------------------------------

A) Related Ubuntu Version
----------------------------------

For systems having Ubuntu with lower versions (lower than 14.04), please notify in case of 
any errors due to OS environments.

Note: We do not support development version corresponding to Windows XP and MacOS, although 
that will be done in some future release.

B) Related GSL
--------------------

Sometimes, after installing the GSL, following error is encountered in the system:

"error while loading shared libraries: libgsl.so.0: cannot open shared object file: No such file or directory"

To resolve this error, check out the following link and perform the steps as mentioned.

https://www.gnu.org/software/gsl/manual/html_node/Shared-Libraries.html

Otherwise, locate the file '.bashrc' within your system, and append the following line at the end of it:

LD_LIBRARY_PATH=/usr/local/lib

Execution
---------------

Download the zipped archieve and extract the codes.

First, extract the zipped archieve named GNU_BFGS2.zip, to unpack the corresponding executable. 

There is also a python file named COSPEDTreeBL.py, which is the main file of this package.

In terminal, go to the directory containing the source codes, and type the following commands:

chmod +x COSPEDTreeBL.py (To change its permission to make it an executable file)

./COSPEDTreeBL.py [options]

Details of the options are mentioned below:

-h, --help show this help message and exit

-I INP_FILENAME, --INPFILE=INP_FILENAME

                   Name of the input file containing input phylogenetic trees

-p INP_FILE_FORMAT, --inpform=INP_FILE_FORMAT

                    1 - input file (containing the input treelist) format is NEWICK (default)                         
                    2 - input file format is NEXUS

-s, --basicscore   (boolean option)

                    If TRUE, support scores of individual couplets follow the rule employed in basic COSPEDTree approach.
                    Otherwise, branch length information of the input trees are also accounted for setting the support score 
                    values of individual couplets. Default is FALSE.
			      			  
-b, --binary  (boolean option)

                    If TRUE, it produces a strictly binary supertree. Otherwise, the output supertree can be non-binary. Default FALSE.
    			  
-d, --dfsref   (boolean option)
			  
                    If TRUE, Multiple parent problem (C2) is tackled by deterministic selection of parent node. 
                    Otherwise, arbitrary selection of parent node is carried out. Default is TRUE.
    			    
-m method_val, --method method_val
			  
                  Value of this option can be 1 to 3.
                  1 - Employes L-BFGS-B method for nonlinear optimization.
                  
                  2 - Employes SLSQP method for nonlinear optimization.

                  3 - C method on GNU library (BFGS) (default).
  
-l loop_count, --loops loop_count
			  
                  No of iterations (loops) employed during the branch length optimization based on the nonlinear programming.
                  only applicable if the method of nonlinear programming is either 1 or 2.
                  generally more iterations yield better optimization \
                  default number of iterations is 15.

-Q QP_Exec_Path, --QPExec QP_Exec_Path
			  
                This is the absolute / relative path of the executable for QP solver (here GNU_BFGS2) which user needs to provide

-w, --weighttaxa  (boolean option)
			                
                Using this option toggles the existing configuration (Default TRUE) \
                if TRUE, then this option assigns different weights to individual couplet relations and corresponding frequencies. 
                The weights are computed according to the size of taxa subset underlying MRCA of that couplet (for corresponding input tree).

--------------------------------------------------------------------------------------------------
Example of a command 
(followed for the results published in the manuscript)
--------------------------------------------------------------------------------------------------

./COSPEDTreeBL.py -I source_tree_input.txt -p1-T supertree_topology_file.txt -t1 > out.txt

command descriptions:

1) -I specifies the input filename:
  
2)  source_tree_input.txt contains the input collection of trees

3) -p option is for specifying the input tree format input file contains the trees in NEWICK format, 
as specified by the option (-p1) (1 stands for newick)

4) -T option is used to specify the custom supertree topology, built from the trees in the file source_tree_input.txt

5) supertree_topology_file.txt contains the supertree topology in either newick 
(preferable) or nexus format

6) -t option is analogous to -p option, to specify the format of supertree topology file 
(1 = newick, 2 = nexus)

All the other options are put in their respective default settings.

The output texts are printed at console. User can redirect the output results to any standard text file by 
using standard redirection operation (>). For example, in the above command, all the detailed results 
(textual descriptions) are redirected to file out.txt.

In addition, one output folder 'CUSTOM_SUPERTREE_QP' is created in the current directory 
containing the supertree topology file. The folder contains the custom input unweighted supertree, 
its weighted version (computed using this package), and text files as outputs of nonlinear programming 
functions employed in the current package.

Utilities
-----------

COSPEDTreeBL requires O(MN^2+N^3) time and O(N^2) space complexity, 
for N input taxa and M input trees.

For any queries, please contact
---------------------------------------

Sourya Bhattacharyya

Department of Computer Science and Engineering

Indian Institute of Technology Kharagpur

email: sourya.bhatta@gmail.com





