# BFS with Intel Intrinsics

The purpose of this project was to implement algebraic BFS and then exploit SIMD to speed up BFS.
Currently there are two implementations, one using CRS format, and the other SlimSell. \
SlimSell uses Intel Intrinsics and requires AVX2 on x86 system.

## Graph generation

Graphs are generated to file using Grapg 500 benchmark generator. These graphs are Kronecker power-law graphs. \
To create a graph, move to bfs/graphgen/ and use:
1) make includes \
2) make locals \
3) make all \
 \
run with ./gen.out \
 \
Flags: \
  -s SCALE, graph scale. Determines number of vertices where nverts = 2^SCALE. \
  -e EDGEFACTOR. Determines average vertex degree where EDGEFACTOR = #edges / #vertices. \
  -f, file. Supply a file name to store the graph in \
  -r, ROOTS. Create a file containing ROOTS number of non-isolated roots to perform BFS with, these are stored in file : ../bin/sSCALE_efEDGEFACTOR_rROOTS.res \
 
 ## BFS

 To compile, move to src/ and use: \
 1) make \
  \
 Hopefully you will see an executable, bfs.out . You may need to modify makefile depending on system. \
  \
 run with ./bfs -f PATH/TO/GRAPH -S SCALE -E EDGEFACTOR \
  \
 Optional flags: \
  -m MODE , mode to run. Mode 1 = CRS only, mode 2 = SlimSell only, mode 3 = CRS and SlimSell \
  -r ROOT , user-specified root to search from. Default root is 0. \
  -R ROOTS -r NUMROOTS, ROOTS is file containing roots, NUMROOTS is maximum number of roots in file ROOTS to peform BFS with. \
  -s SIGMA , sorting degree for SlimSell. if mode > 1, then use sorting degree. Sorting degree = 2^SIGMA, where 0 <= SIGMA <= SCALE. \
 \
