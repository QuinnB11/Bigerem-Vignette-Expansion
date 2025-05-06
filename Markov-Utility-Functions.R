library("bigergm")
library("GGally")
library("igraph")
library("network")  # For conversion to igraph
library("ggplot2")

# Load the network
data("state_twitter")

# Set reproducibility
set.seed(123)

# --- Convert to igraph object ---
# The network needs to be converted to igraph for matrix operations
state_twitter_igraph <- intergraph::asIgraph(state_twitter)

# --- Visualize the Original Network ---
ggnet2(state_twitter, size = 1, color = "state", edge.alpha = 0.1) + 
  guides(color = "none") + 
  scale_colour_viridis_d()

# --- Simulate Network Evolution with Markov Chains ---

# Get the number of nodes
n_nodes <- vcount(state_twitter_igraph)

# Define the transition probability matrix
transition_matrix <- matrix(runif(n_nodes * n_nodes, 0, 1), n_nodes, n_nodes)
# Normalize rows
transition_matrix <- transition_matrix / rowSums(transition_matrix)

# Simulate the Markov chain over T steps
T <- 50
simulated_networks <- list()

# Create Markov chain evolution function
for (t in 1:T) {
  adj_matrix <- as_adjacency_matrix(state_twitter_igraph, sparse = FALSE)
  
  # Apply transition rules
  for (i in 1:n_nodes) {
    for (j in 1:n_nodes) {
      if (runif(1) < transition_matrix[i, j]) {
        adj_matrix[i, j] <- 1
      } else {
        adj_matrix[i, j] <- 0
      }
    }
  }
  
  # Store the evolved network
  evolved_graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  
  # Convert back to network object for ggnet2 visualization
  evolved_network <- intergraph::asNetwork(evolved_graph)
  simulated_networks[[t]] <- evolved_network
}

# --- Visualize the Network Evolution with ggnet2 ---
# Select 6 steps evenly distributed across the simulation
steps_to_plot <- seq(1, T, length.out = 6)

# Plot each step using ggnet2
par(mfrow = c(2, 3))  # 2x3 grid
for (i in steps_to_plot) {
  cat(paste("Plotting step", i, "\n"))
  g <- simulated_networks[[i]]
  
  # Use ggnet2 for visualization
  print(ggnet2(g, size = 2, color = "skyblue", edge.color = "grey", edge.alpha = 0.5) +
          ggtitle(paste("Step", i)))
}

# --- Fit HERGM Model on the Evolved Network ---
# Use the last evolved network
final_network <- simulated_networks[[T]]

# Define the HERGM formula
model_formula <- final_network ~ edges + triangle + nodematch("state")

# Estimate HERGM
res <- bigergm(
  object = model_formula,
  n_blocks = 4, 
  n_MM_step_max = 100,
  tol_MM_step = 1e-6,
  estimate_parameters = TRUE,
  clustering_with_features = TRUE,
  check_block_membership = TRUE, 
  initialization = "walktrap"
)

# --- Plot Convergence of the Lower Bound with ggplot2 ---
# Create a dataframe for lower bound values
df <- data.frame(
  Iteration = 1:length(res$MM_lower_bound),
  LowerBound = res$MM_lower_bound
)

# Plot convergence
ggplot(df, aes(x = Iteration, y = LowerBound)) +
  geom_line(color = "blue") +
  labs(title = "Convergence of the Lower Bound", x = "Iterations", y = "Lower Bound") +
  theme_minimal()

# --- Plot the Clustered Network with ggnet2 ---
# Extract the clustered network
clustered_graph <- res$final_network

# Convert the final igraph object back to a network object for ggnet2 visualization
clustered_network <- intergraph::asNetwork(clustered_graph)

# Plot clustered network using ggnet2
ggnet2(clustered_network, size = 2, color = "orange", edge.color = "grey", edge.alpha = 0.5) +
  ggtitle("Clustered Network")

# --- Inspect Parameter Estimates ---
summary(res$est_between)
