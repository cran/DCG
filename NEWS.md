# DCG 0.9.2
* Fixing vignettes. Put external plots in vignettes folder. Use png::readPNG and grid::grid.raster to render plots in vignettes.

# DCG 0.91
* Importing data. When importing network data as a similarity matrix, 
network data (edgelist or matrix) first gets transformed into a undirected adjacency matrix. Then the adjacency matrix gets transformed into a similarity matrix.
