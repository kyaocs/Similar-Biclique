# Enumerating Maximal Similar-Bicliques in Large Bipartite Graphs

This project aims to enumerate all maximal similar-bicliques in a bipartite graph.

MSBE is the executable, and is compiled on Ubuntu 18.04.5, with -O3 optimization.

Folder "datasets" contains an exmpale bipartite graph bi_writer.bin, which is downloaded from [KONECT](http://konect.cc/networks/dbpedia-writer/). 


## Running Format

./MSBE [1]input graph [2]build index or not (0/1) [3]alpha (0.01->100) [4]gamma (0.01 -> 1) [5]index name (LG/SS) [6]load index or not (0/1) [7] index name (LG/SS) [8]vertex reduction method (1/2) [9]epsilon (similarity threshold) [10]tau (size constraint)

**Running example for building indexLG**

./MSBE ./datasets/bi_writer.bin 1 1 0 LG

**Running example for indexLG based enumeration**

./MSBE ./datasets/bi_writer.bin 0 1 0 LG 1 LG 2 0.5 3 

**Running example for building indexSS**

./MSBE ./datasets/bi_writer.bin 1 1 0.3 SS

**Running example for indexSS based enumeration**

./MSBE ./datasets/bi_writer.bin 0 1 0.3 LG 1 SS 2 0.5 3 

**Running example for enumeration without any index**

./MSBE ./datasets/bi_writer.bin 0 0 0 0 0 0 1 0.5 3


### Note

The input graph should be in "binary" format by default. In folder "datatsets", there is an example bipartitie graph bi_writer.bin. Here, edgelist2binary is the executable to transform a "txt" graph into our binary form.

When argument[2] is 1, MSBE will build the index according to arguments[3][4][5], and arguments[6]-[9] will be ignored. The index will be stored in folder "datatsets" and named by the arguments, i.e., "input graph" + "alpha x 100" + "gamma x 100" + "LG/SS.bin". 

After the index is constructed, it is ready to make index based enumeration. Specifically, set argument[2]=0, argument[6]=1 and argument[7]=LG/SS to load the corresponding index. (Note that, the name of the loaded index will be dirived in a same way as above, i.e., "input graph" + "alpha x 100" + "gamma x 100" + "LG/SS.bin". Thus, make sure arguments [3] and [4] are correctly set when loading the corresponding index.)

As for atgument[8], "1" means vertex reduction without index, "2" means index based vertex reduction. Thus, if you set atgument[8]=2, make sure the index is loaded, i.e., argument[6]=1.


## Graph Format

### txt version

number of L side vertices \t number of R side vertices \t number of edges \n

v0 \t v1

v0 \t v2

...

### Note

In folder "datasets", edgelist2binary could transform txt file to binary file.