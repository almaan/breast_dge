Program to identify isolated regions (such as tumors) from annotated ST-data. Clusters of spatially coherent
and isolated spots will be marked accordingly to cluster index. The number of clusters
within the sample does not have to be prespecified (and cannot be either).


Files should be provided with a 4 or 5 digit patient id, if replicates are present then this should be indicated by
appending an underscore followed by the replicate id (one capital letter followed by a digit).

 Examples are
 * feature-file-*12345*.tsv (only patient id, no replicate)
 * feature-file-*12345_A1*.tsv (both patient id and replicate id

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
ST_feature_extract --input /Annotated_File_Dir/single_file --output result_of_single_file.tsv --norm -1 --max_dist 1.5 --min_total_spots 4 --save_plot
```

_the above script will assign two spots with a distance 1.5 (measured by the infty-norm) as neighbours. At least 4 spots must be found within a cluster
in order to classify it as a tumor. Smaller clusters will be discarded. The result will be saved both as a .tsv file and an image._



