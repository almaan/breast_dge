Program to identify isolated regions (such as tumors) from annotated ST-data. Clusters of spatially coherent
and isolated spots will be marked accordingly to cluster index. The number of clusters
within the sample does not have to be prespecified (and cannot be either).

Each cluster consists of two type of spots, inner spots and boundary spots. A spot is considered an inner spot
if it has at least min_spot neighbours. An inner spot can "propagate" the cluster and will assign all of it's 
neighbors to the same cluster as itself. Boundary spots belong to a tumor but has less than min_spot neighbors
meaning that they cannot propagate a cluster. 

An option to mask clusters below a given size (number of spots) is given, which can be specified with 
the argument     .

Input can be either a directory of multiple files to be processed simultaneously or a single file. For
multiple-file mode enter a directory as output, if directory does not exist it will be created in current
working directory. For single file mode enter a file name as output.

Each file to be processed must have the sample name specified somewhere within the file according to the pattern
"XY#####" where X and Y are arbitrary alphabetical characters within the range [A-Z] (upper or lowercase) and "#"
represents a digit in the range[0-9]. The x-coordinates and y-coordinates should be given in columns named "xcoord"
and "ycoords" respectively. The column containing the annotated feature can be named arbitrarily (default is "tumor"),
and should be passed as an argument (feature) if not "tumor". The annotation of interest (only support) for one as of now
should be passed as an argument if other than default ("tumor".)


Three parameters are used for cluster generation:
    
    max_distance    - the maximum manhattan distance for two spots to be considered neihbours 
    min_spots       - the number of neihbours a spot must have to be considered an inner point of the cluster
    min_total_spots - the minimum number of spots that a cohort of spots must have to be considered a tumor

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
ST_feature_extract --input /Annotated_File_Dir/single_file --output result_of_single_file.tsv --min_spot 2 --max_dist 3 --min_total_spots 5 --save_plot
```



