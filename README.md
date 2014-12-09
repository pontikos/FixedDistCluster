FixedDistCluster
================

C code for memory efficient fixed distance agglomerative clustering in 3 dimensions (can be extended to multiple dimensions).
Uses exactly as much memory as datapoints thanks to its linked list implementation.
Uses region-growing to incorporate points which are within a certain distance threshold _d_ of each other within the cluster.
Obviously the choice of the _d_ is critical and should be larger for higher dimensions.

