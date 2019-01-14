Program to identify isolated regions (such as tumors) from annotated ST-data. Clusters of spatially coherent
and isolated spots will be marked accordingly to cluster index. The number of clusters
within the sample does not have to be prespecified. Extraction is based on _connected graphs_ within treating the spots within one sample as nodes of a graph. 

The program supports multiple and single file input. In either case the output should be given as a directory. If no such directory exists this will be created.

If the filenames are "tagged" as to contain sample\_id and replicate name. This can be included in the output file as columns as well as in the output name. Use the flag "--tagged" either providing two regex patterns for (1) sample id and (2) replicate. If tagged flag is used but no arguments provided the default pattern shown below will be used.

**default pattern**
XY#####\_W#

where X,Y and W represents arbitrary capital letters, and # arbitrary digits

The x-coordinates and y-coordinates should be given in columns named "xcoord"
and "ycoords" respectively. The column containing the annotated feature can be named arbitrarily (default is "tumor"),
and should be passed as an argument (feature) if not "tumor". The annotation of interest (only support) for one as of now
should be passed as an argument if other than default ("tumor".)


Three parameters are used for cluster generation:
    
    max_distance    - the maximum manhattan distance for two spots to be considered neihbours 
    min_total_spots - the minimum number of spots that a cohort of spots must have to be considered a tumor
    norm            - Which Minkowski norm to be used in distance estimation

# Installation
## 1. Clone repo
```bash
cd /install_directory
git clone https://github.com/almaan/breast_dge.git
```
## 2. Install (makes a soft link)
```bash
cd /scripts/tumor_extraction
chmod +x install.sh
sudo install.sh
```
## 3. Run
```bash
ST_feature_extract --input /Annotated_File_Dir/single_file --output result_of_single_file.tsv --norm -1 --max_dist 1.5 --min_total_spots 4 --save_plot --tagged
```

_the above script will assign two spots with a distance 1.5 (measured by the infty-norm) as neighbours. At least 4 spots must be found within a cluster
in order to classify it as a tumor. Smaller clusters will be discarded. The result will be saved both as a .tsv file and an image._


## Tips
The positing of the spots is not perfect, thus some marginal for the distance should be used. As to illustrate; If the Von Neumann neighbourhood is to be used as neighbourhood, then rather than using the l1-norm with distance 1 use distance 1.3 or 1.5 allowing for the slight misplacement of spots to be accounted for.
