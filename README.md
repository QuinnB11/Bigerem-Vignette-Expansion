# BIGERGM Vignette Expansion

This script demonstrates the use of the [`BIGERGM`](https://www.corneliusfritz.com/bigergm/articles/bigergm.html) R package to fit a block-informed exponential random graph model (ERGM) to a large real-world network dataset from Amazon product co-purchasing communities.

## Overview

The pipeline performs the following steps:

### 1. Setup and Dependencies

* Loads necessary R packages (`bigergm`, `network`, `sna`, `microbenchmark`).
* Ensures Python dependencies are available for Infomap clustering via `bigergm::py_dep()`.

### 2. Load Community Memberships

* Reads the Amazon community membership file: `com-amazon.top5000.cmty.txt`.
* Parses each line as a group of node IDs (representing product co-purchase communities).
* Stores the size of each group.

### 3. Sample 500 Disjoint Groups

* Randomly samples 500 non-overlapping groups, each with 10–80 nodes.
* Ensures that no nodes are shared between sampled groups.
* Saves the resulting node and group indices.

### 4. Load and Filter Edge List

* Loads the full edge list from `amazon_graph.rds`.
* Filters edges to include only those between sampled nodes.
* Saves filtered edge data and node-group mapping.

### 5. Create Network and Assign Attributes

* Constructs an undirected `network` object from the filtered edge list.
* Aligns vertex attributes with group assignments (`true_partition`).
* Adds node-level covariates:

  * **Degree Level**: categorized as *low*, *medium*, or *high* based on degree quantiles.
  * **Block Size Level**: categorized as *small*, *medium*, or *large* based on block size quantiles.

### 6. Specify Model Formula

Defines an ERGM with the following components:

```r
network ~ edges + transitiveties +
          nodematch("degree_level") +
          nodematch("block_size_level")
```

### 7. Estimate BIGERGM

* Uses `bigergm()` to fit the model with two different community detection initializations:

  * **Walktrap** (in R)
  * **Infomap** (via Python)
* Extracts estimated coefficients for both *between-block* and *within-block* structures.

### 8. Simulate Networks and Benchmark

* Simulates networks using estimated parameters from both initializations.
* Measures time taken for simulations.
* Visualizes simulated networks side-by-side.

## Output Files

* `amazon_edges_500.RData`: Filtered edge list for 500 sampled groups.
* `amazon_nodes_500.RData`: Node IDs and corresponding group IDs.
* Estimated model results and simulation outputs printed to console and visualized.

## Requirements

* R packages: `bigergm`, `network`, `sna`, `microbenchmark`
* Python installation with Infomap (automatically managed by `bigergm::py_dep()`)


