# ----------------------
# Setup and Dependencies
# ----------------------
library(bigergm)
library(network)
library(sna)
library(microbenchmark)
bigergm::py_dep()  # Ensure Python dependencies for Infomap are available

# ----------------------
# Step 1: Load and Parse Community File
# ----------------------
x <- scan("com-amazon.top5000.cmty.txt", what = "", sep = "\n")
y <- strsplit(x, "[[:space:]]+")
group_sizes <- unlist(lapply(y, length))

# ----------------------
# Step 2: Sample 500 Groups of Size 10–80 Without Overlap
# ----------------------
set.seed(123456)
sampling_order <- 1:5000
iterator <- 0
node_indeces <- c()
group_indeces <- c()
group_iterator <- 0
n_nodes <- 0

while (group_iterator < 500) {
  iterator <- iterator + 1
  target_group <- y[[sampling_order[iterator]]]
  temp_length <- length(target_group)

  # Skip if group size is outside desired range
  if (temp_length < 10 || temp_length > 80) next

  # Ensure no overlapping nodes
  if (sum(target_group %in% node_indeces) == 0) {
    group_iterator <- group_iterator + 1
    node_indeces <- c(node_indeces, target_group)
    n_nodes <- n_nodes + temp_length
    print(c(iterator, n_nodes))
    group_indeces <- c(group_indeces, rep(group_iterator, temp_length))
  }
}

# ----------------------
# Step 3: Filter Edge List to Sampled Nodes
# ----------------------
t.mtx <- readRDS("amazon_graph.rds")

# Keep only edges where both nodes are in the sampled list
output <- t.mtx[(t.mtx[,1] %in% node_indeces + t.mtx[,2] %in% node_indeces) == 2, ]
rm(t.mtx)

# Save filtered data
save(output, file = "amazon_edges_500.RData")
nodes <- list(node_indeces = node_indeces, group_indeces = group_indeces)
save(nodes, file = "amazon_nodes_500.RData")

# ----------------------
# Step 4: Create Network and Assign Blocks
# ----------------------
load("amazon_edges_500.RData")
load("amazon_nodes_500.RData")

node_indeces <- nodes$node_indeces
group_indeces <- nodes$group_indeces
n_nodes <- length(node_indeces)

# Convert storage mode for compatibility
storage.mode(output) <- "character"
storage.mode(group_indeces) <- "character"
storage.mode(node_indeces) <- "character"

# Create network object
network <- as.network.matrix(output, directed = FALSE)

# Align node order with network vertex order
names <- network::get.vertex.attribute(network, "vertex.names")
new_order <- sapply(1:n_nodes, function(i) which(node_indeces == names[i]))

# Assign block (group) membership to each node
true_partition <- group_indeces[new_order]
# Uncomment if needed: set.vertex.attribute(network, "block", as.double(true_partition))

# ----------------------
# Step 5: Add Covariates
# ----------------------

# Add degree-level (low/medium/high)
deg <- sna::degree(network)
degree_level <- cut(deg, breaks = quantile(deg, probs = c(0, 0.33, 0.66, 1)),
                    labels = c("low", "medium", "high"), include.lowest = TRUE)
set.vertex.attribute(network, "degree_level", as.character(degree_level))

# Add block-size-level (small/medium/large)
group_sizes_table <- table(true_partition)
block_size <- group_sizes_table[as.character(true_partition)]
block_size_level <- cut(block_size, breaks = quantile(block_size, probs = c(0, 0.33, 0.66, 1)),
                        labels = c("small", "medium", "large"), include.lowest = TRUE)
set.vertex.attribute(network, "block_size_level", as.character(block_size_level))

# ----------------------
# Step 6: Define Model Formula
# ----------------------
model_formula <- network ~ edges + transitiveties +
                 nodematch("degree_level") +
                 nodematch("block_size_level")

# ----------------------
# Step 7: Estimate BIGERGM Model
# ----------------------

# Define shared model arguments
base_args <- list(
  object = model_formula,
  n_blocks = 2,  # NOTE: Only 2 blocks used for demonstration; increasing to 100 would improve performance accuracy
  n_MM_step_max = 100,
  tol_MM_step = 1e-6,
  estimate_parameters = TRUE,
  clustering_with_features = TRUE,
  check_block_membership = TRUE
)

# Run Walktrap initialization
res_walktrap <- do.call(bigergm, c(base_args, list(initialization = "walktrap")))

# Run Infomap initialization (Python)
res_infomap <- do.call(bigergm, c(base_args, list(
  initialization = "infomap",
  use_infomap_python = TRUE
)))

# Display estimated coefficients
print(res_walktrap$est_between$coefficients)
print(res_walktrap$est_within$coefficients)

# ----------------------
# Step 8: Simulate Networks
# ----------------------

# Time simulation using Infomap results
time_sim_infomap <- system.time({
  sim_infomap <- simulate_bigergm(
    formula = model_formula,
    coef_between = res_infomap$est_between$coefficients,
    coef_within = res_infomap$est_within$coefficients,
    nsim = 1,
    output = "network"
  )
})

# Time simulation using Walktrap results
time_sim_walktrap <- system.time({
  sim_walktrap <- simulate_bigergm(
    formula = model_formula,
    coef_between = res_walktrap$est_between$coefficients,
    coef_within = res_walktrap$est_within$coefficients,
    nsim = 1,
    output = "network"
  )
})

# Print timing results
print("Simulation time using Walktrap results:")
print(time_sim_walktrap)

print("Simulation time using Infomap results:")
print(time_sim_infomap)

# ----------------------
# Step 9: Plot Simulated Networks
# ----------------------
par(mfrow = c(1, 2))
plot(sim_infomap, main = "Simulated Network - Infomap")
plot(sim_walktrap, main = "Simulated Network - Walktrap")

# ----------------------
# Step 10: Goodness-of-Fit (GOF)
# ----------------------

# GOF for Walktrap
gof_walktrap <- gof(
  object = res_walktrap,
  nsim = 100,
  compute_geodesic_distance = TRUE,
  seed = 1234,
  start_from_observed = TRUE,
  type = "within",
  control_within = ergm::control.simulate.formula(MCMC.burnin = 1000, MCMC.interval = 1000)
)

# GOF for Infomap
gof_infomap <- gof(
  object = res_infomap,
  nsim = 100,
  compute_geodesic_distance = TRUE,
  seed = 1234,
  start_from_observed = TRUE,
  type = "within",
  control_within = ergm::control.simulate.formula(MCMC.burnin = 1000, MCMC.interval = 1000)
)

# Save GOF plots
png("gof_walktrap.png", width = 800, height = 600)
plot(gof_walktrap, main = "GOF - Walktrap Initialization")
dev.off()

png("gof_infomap.png", width = 800, height = 600)
plot(gof_infomap, main = "GOF - Infomap Initialization")
dev.off()
