ShapeDB - Indexed molecular shape search.

Building
========
The project is built using Eclipse.  Auto-generated makefiles can be found
in the Release directory.  You will likely have to edit these.

Usage
=====

*Important:*  Conformations are stored using the provided coordinates.
You must provide molecular conformations aligned to your desired frame of reference.
There are python scripts in the scripts directory to assist with aligning to
a canonical frame of reference (as done with VAMS) or to fragments (FOMS).

### Building a database
```
ShapeDB -Create -in selfaligned_actives_f.sdf -db foo.db
```

### Searching a database
Similarity search:
```
ShapeDB -NNSearch -k 3  -ligand salig.pdb -db foo.db -print -single-conformer -out hits.sdf
```
Distance constrain search:
```
ShapeDB -DCSearch -less 2.0 -more 2.0 -ligand salig.sdf -receptor sarec.pdb -db foo.db -print -single-conformer -out hits.sdf
```


```
USAGE: ShapeDB [options]

OPTIONS:
  Operation to perform:
    -Create                              - Create a molecule shape index.
    -CreateTrees                         - Create only the shapes with no index.
    -CreateFromTrees                     - Create an index from already created shapes.
    -CreateFromMira                      - Create an index from already created shapes.
    -MiraSearch                          - Nearest neighbor searching of mira objects
    -DCMiraSearch                        - Distance constraint searching of mira objects
    -NNSearch                            - Nearest neighbor search
    -DCSearch                            - Distance constraint search
    -MolGrid                             - Generate molecule grid and debug output
    -MergeMira                           - Merge mira files into MIV/MSV
    -BatchSearch                         - Read in a jobs file for batch processing
    -BatchMiraSearch                     - Read in a mira jobs file for batch processing
    -BatchDB                             - Read in a jobs file for batch processing of a list of directories
    -SearchAllPointCombos                - Search using all possible subsets of interaction points
  Metric for distance between shapes:
    -rel-volume                          - Relative volume difference
    -abs-volume                          - Absolute volume difference
    -hausdorff                           - Hausdorff distance
    -rel-triple                          - Triple including selectivity
    -abs-triple                          - Triple including selectivity (absolute)
    -include-exclude                     - For comparing with include/exclude constraints
    -relvol-exclude                      - Volume comparison of ligand, exclusion comparison of receptor
  -clear-cache                           - Clear file cache between each benchmarking run
  -clear-cache-first                     - Clear file cache before each benchmarking run
  -db=<string>                           - Database file
  -files=<string>                        - Files for MiraMerge
  Packing algorithm:
    -full-merge                          - Greedy full merge.
    -greedy-merge                        - Greedy iterative merge.
    -match-merge                         - Optimal matching merging.
    -spectral                            - Spectral packing
  -h                                     - Retain hydrogens in input molecules
  -help                                  - Display available options (--help-hidden for more)
  -in=<string>                           - Input file
  -interaction-distance=<number>         - Distance between ligand/receptor atom centers that are considered interacting
  -interaction-max-cluster-dist=<number> - Maximum span of interaction point cluster.
  -interaction-min-cluster=<number>      - Minimum size of interaction point cluster.
  -interaction-point-radius=<number>     - Amount to grow interaction points into ligand shape.
  -k=<uint>                              - k nearest neighbors to find for NNSearch
  -kcenters=<uint>                       - number of centers for ksample-split
  -knn=<uint>                            - K for knn graph creation
  -ksamplex=<uint>                       - multiplictive factor for ksampling
  -less=<number>                         - Distance to reduce query mol by for constraint search (default 1A).
  -ligand=<string>                       - Molecule to use for minimum included volume
  -max-dim=<number>                      - Maximum dimension.
  -mira-miv=<string>                     - Mira shape to use for minimum included volume
  -mira-msv=<string>                     - Mira shape to use for maximum surrounding volume
  -more=<number>                         - Distance to increase query mol by for constraint search (default 1A).
  -nn-search-all                         - Perform exhaustive interaction point search with NNSearch as well as DCSearch
  -nn-threshold=<number>                 - Similarity cutoff for NNSearch
  -out=<string>                          - Output file
  -pack=<uint>                           - Maximum quantities per a node
  -print                                 - Print text summary of output
  -probe-radius=<number>                 - Radius of water probe for exmol only
  -receptor=<string>                     - Molecule to use for excluded volume
  -resolution=<number>                   - Best resolution for shape database creation.
  -sentinals=<uint>                      - Number of sentinals for knn initialization (zero random)
  -single-conformer                      - Output the single best conformer
  Spectral packing sub-algorithm:
    -sort-dense                          - Simple sort followed by dense packing
    -sort-partition                      - Simple sort followed by largest separator partition packing
    -cluster-eigen                       - Cluster eigen values using greedy packer
    -relax                               - Cluster form relaxation values
  -sproxel-color=<string>                - Sproxel voxel descriptor
  -superdepth=<uint>                     - Depth to descend to create aggregrated super root
  -switch-to-pack=<uint>                 - Cutoff to trigger packing in nodes and leaves
  -time-trials=<uint>                    - Number of runs to get average for benchmarking
  Metric for cluster packing distance:
    -ave-dist                            - Use 'average' metric between MIV/MSV representations of clusters
    -complete-dist                       - Use complete linkage value between cluster members
    -single-dist                         - Use single linkage value between cluster members
    -total-dist                          - Use total (sum) linkage value between cluster members
  -use-interaction-points                - Analyze the ligand-receptor complex and generate a minimum shape centered around different interaction points.
  -use-unnorm                            - Use unnormalized laplacian in spectral packing
  -v                                     - Verbose output
```

Citation
======

Indexing volumetric shapes with matching and packing
DR Koes, CJ Camacho
Knowledge and information systems 43 (1), 157-180

Shape‐based virtual screening with volumetric aligned molecular shapes
DR Koes, CJ Camacho
Journal of computational chemistry 35 (25), 1824-1834

Fragment oriented molecular shapes
E Hain, CJ Camacho, DR Koes
Journal of Molecular Graphics and Modelling 66, 143-154
