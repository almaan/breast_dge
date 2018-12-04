# TODO


## Inteded relatively easy analysis

* Tumor vs. Non-Tumor within Each patient (23 runs)
* Tumor vs- Non-Tumor within all patients (interesting to contrast with previous) (1 HUGE)
* Interpatient tumor comparision (those that have multiple tumors)
* Linnea does not think comparision based on subtype is good

**Suggestions to Linnea : Make ROC/AUC to find the perfect parameter for thresholding**

* Have suggested it through slack. Will speak more later.
* Could potentially write small software for her to do that
* Use Youden's J statistics to find a good threshold



## Code:

* Write script for pooling [x]
 * User can choose k-neighboors to use for each pooling event (important that always same numbers are used, since average is not taken)
 * Add option to filter away ambigious genes [x]
 * User can choose how many samples to pool (replacement is true)[x]
 * Add flags for k\_neighbours, n\_samples and lim [x]
 * User can choose limit (distance) withing which spots have to be to be counted as neighbord, if negative will be set to Inf  [x]
 * **Must comment on code, non-documented so far**

* Implement multifeature DGE [ ]
* Implement logger as to register how many spots that are dropped, and procedure [ ]
* Look at parallelization options, try to put this up to the cluster [ ]
* Add timer to logger as to get a sense of how long each analysis is [ ]
* enable gene-selection (would be interesting to look at say PAM50) [/] _implemented_ but not _evaluated_
* Could one subset the spots, not looking at all of them in order to increase run time [x]
* Add version print on-top [ ]

## Remember
* In DGE\_analysis.r the count matrix is taken as samples x genes
* In poolf.r the count matrix was originally taken as genes x samples. This is why a default transposition is used. Make sure to not forget this.
* Had problem of using generated idx-list for selected feature but over whole count-matrix. Have remedied this


## Ideas and thoughts (spontane)
* **number of spots** pooled seems to be more significant than samples --> less than or equal to 50 give poor results
* Run multiple DGE for same samples and take intersection of DGE-genes as to make sure pooling does not interfere
 * Have looked at intersection, out of 500 top genes 57 matched. Will do more extensive analysis
