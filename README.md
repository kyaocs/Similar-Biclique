# Identifying Similar-Bicliques in Bipartite Graphs

This project aims to enumerate all maximal similar-bicliques in a bipartite graph.

msbe is the executable, and is compiled on Ubuntu 18.04.5, with -O3 optimization.

Folder "datasets" contains an exmpale bipartite graph bi_github.bin, which is downloaded from [KONECT](http://konect.cc/networks/dbpedia-writer/). 


## Running Format

./msbe [1] input graph  [2] build index flag (0/1)  [3] index name (LG/SS)  [4] alpha (0.01～100)  [5] gamma (0.01～1)  [6] load index flag (0/1)  [7] index name (LG/SS)  [8] alpha (0.01～100)  [9] gamma (0.01～1) [10] epsilon (similarity constraint)  [11] tau (size constraint)

**Running example for building indexLG (with \alpha=1)**

./msbe ./datasets/bi_github.bin 1 LG 1 0

**Running example for building indexSS (with \alpha=1,\gamma=0.3)**

./msbe ./datasets/bi_github.bin 1 SS 1 0.3

**Running example for indexLG based enumeration (with \epsilon=0.5, \tau=3)**

./msbe ./datasets/bi_github.bin 0 LG 1 0 1 LG 1 0 0.5 3

**Running example for indexSS based enumeration (with \epsilon=0.5, \tau=3)**

./msbe ./datasets/bi_github.bin 0 SS 1 0.3 1 SS 1 0.3 0.5 3 

**Running example for enumeration without any index (with \epsilon=0.5, \tau=3)**

./msbe ./datasets/bi_github.bin 0 0 0 0 0 0 0 0 0.5 3 


### Note

In the following, I will explain how to set the arguments to execute the code properly. 

In summary, argument[1] is the input graph, arguments[2]-[5] control the index building, arguments[6]-[9] specify which index will be loaded, arguments[10][11] are the two important parameters \epsilon and \tau.

When argument[2]=1, MSBE will build the index according to arguments[3][4][5], and arguments[6]-[11] will be ignored. The index will be stored in the folder "datatsets" and named by the arguments, i.e., "input graph" + "alpha x 100" + "gamma x 100" + "LG/SS.bin". Note that for LG index, argument[5] (i.e., \gamma) is not needed and should be set as 0.

After the index is constructed, it is ready to make index based enumeration. Specifically, set argument[2]=0 and argument[6]=1. Set arguments[7][8][9] to load the corresponding index. (Note that, the name of the loaded index will be dirived in a same way as above, i.e., "input graph" + "alpha x 100" + "gamma x 100" + "LG/SS.bin". Thus, make sure arguments [8][9] are correctly set when loading the index. For LG index, argument[9] (i.e., \gamma) is not needed and should be set as 0.)

To run the algorithm without any index, just set arguments[2]-[9] as 0.

## Graph Format

The input graph should be in "binary" format by default. In folder "datatsets", there is an example bipartitie graph bi_github.bin. Here, edgelist2binary is the executable to transform a "txt" graph into our binary form. 

Our algorithms also support "txt" graph, the txt version should be in the following format:

number of L side vertices \t number of R side vertices \t number of edges \n

v0 \t v1

v0 \t v2

...
