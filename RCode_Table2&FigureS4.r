library(SpiecEasi)
library(igraph)
library(reshape2)
library(doBy)
library(ggplot2)
library(picante)
library(agricolae)
library(ggpubr)

# Read datasets
asv <- read.csv('asv_table.csv', row.names = 1)
asv_relative <- asv/rep(colSums(asv),each = nrow(asv))

# Retain species with a frequency greater than 1 to reduce random errors
asv <- asv[which(apply(asv, 1, function(x) sum(x>0))>1), ]
asv_relative <- asv_relative[rownames(asv), ]


############## Construct SparCC network

# Estimate compositionality-robust correlations among ASVs to obtain a sparse correlation matrix
asv <- data.frame(t(asv))

set.seed(123)
asv_sparcc <- sparcc(asv, iter = 20, inner_iter = 10, th = 0.1)
sparcc0 <- asv_sparcc$Cor

colnames(sparcc0) <- colnames(asv)
rownames(sparcc0) <- colnames(asv)

# Bootstraps 100 times to obtain random correlation matrix
set.seed(123)
n = 100

for (i in 1:n) {
	asv.boot <- sample(asv, replace = TRUE)
	asv.sparcc_boot <- sparcc(asv.boot, iter = 20, inner_iter = 10, th = 0.1)
	sparcc_boot <- asv.sparcc_boot$Cor
	colnames(sparcc_boot) <- colnames(asv.boot)
	rownames(sparcc_boot) <- colnames(asv.boot)
	write.table(sparcc_boot, paste('sparcc_boot', i, '.txt', sep = ''), sep = '\t', col.names = NA, quote = FALSE)
}

# Based on the sparse correlation matrix and the results of 100 bootstraps, obtain the pseudo p-values
p <- sparcc0
p[p!=0] <- 0

for (i in 1:n) {
	p_boot <- read.delim(paste('sparcc_boot', i, '.txt', sep = ''), sep = '\t', row.names = 1)
	p[abs(p_boot)>=abs(sparcc0)] <- p[abs(p_boot)>=abs(sparcc0)] + 1
}

p <- p / n

# Remove autocorrelation
diag(sparcc0) <- 0

# Retain values where |sparse correlation| >= 0.75 and p < 0.01 (correlations that do not meet these criteria are assigned zero).
sparcc0[abs(sparcc0)<0.75 | p>=0.01] <- 0

# Remove rows and columns that are all zeros and output adjacent matrix
sparcc0 <- sparcc0[rowSums(sparcc0)!=0,rowSums(sparcc0)!=0]
write.table(sparcc0, 'neetwork.adj.txt', col.names = NA, sep = '\t', quote = FALSE)


############## Split subnetwork for each sample

adjacency_unweight <- read.delim('neetwork.adj.txt', row.names = 1, check.names = FALSE)
adjacency_unweight[adjacency_unweight>0] <- 1
adjacency_unweight[adjacency_unweight<0] <- -1
adjacency_unweight1 <- adjacency_unweight
adjacency_unweight1[adjacency_unweight1!=0] <- 1
g <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight1), mode = 'undirected', weighted = NULL, diag = FALSE)

asv <- data.frame(t(asv))
asv <- asv[rownames(adjacency_unweight1), ]
sub_graph <- list()
for (i in names(asv)) {
	sample_i <- asv[i]
	select_node <- rownames(sample_i)[which(sample_i != 0)]
	sub_graph[[i]] <- subgraph(g, select_node)
	sub_graph[[i]] <- delete.vertices(sub_graph[[i]], names(degree(sub_graph[[i]])[degree(sub_graph[[i]]) == 0]))
}
sub_graph


############## Calculate topological features for each subnetwork

sample_id <- c()
nodes_num <- c()
edges_num <- c()
Degree <- c()
average_path_length <- c()
graph_diameter <- c()
graph_density <- c()
clustering_coefficient <- c()
betweenness_centralization <- c()
degree_centralization <- c()
modularity <- c()
complementarity  <- c()
tree <- ape::read.tree('asv_seq.tre')

for(i in 1:length(sub_graph)){
	sample_id <- c(sample_id, names(sub_graph[i]))
	nodes_num <- c(nodes_num, length(V(sub_graph[[i]])))  #number of nodes
	edges_num <- c(edges_num, length(E(sub_graph[[i]])))  #number of edges
	Degree <- c(Degree, mean(degree(sub_graph[[i]])))  #degree
	clustering_coefficient <- c(clustering_coefficient, transitivity(sub_graph[[i]]))  #clustering coefficient
	modularity <- c(modularity, modularity(sub_graph[[i]], membership(cluster_fast_greedy(sub_graph[[i]]))))  #modularity
	
	# Use the total branch length of the phylogenetic tree as a proxy for community-level functional complementarity, which defined as the total branch length of the functional dendrogram of co-occurring species
	asv_i <- asv[names(V(sub_graph[[i]])), ]
	tree_i <- prune.sample(data.frame(t(asv_i)), tree)
	complementarity  <- c(complementarity , sum(ape::compute.brlen(tree_i)$edge.length))
}

sub_graph_stat <- data.frame(sample_id, nodes_num, edges_num, degree = Degree, clustering_coefficient, modularity, complementarity)
head(sub_graph_stat)

# Calculate the positive and negative connectivity of each species
g <- read.delim('neetwork.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)
g[g==0] <- NA
asv_relative <- asv_relative[rownames(g), ]

r_pos <- c()
r_neg <- c()
spe_i <- c()
for (i in names(g)) {
	co <- na.omit(g[[i]])
	r_pos <- c(r_pos, mean(co[co>0], na.rm = TRUE))
	r_neg <- c(r_neg, mean(co[co<0], na.rm = TRUE))
	spe_i <- c(spe_i, i)
}
r <- data.frame(spe_i, r_pos, r_neg)
r

# Calculate the positive and negative cohesion of each sample
C_pos <- c()
C_neg <- c()
id <- c()
 
for (j in names(asv_relative)) {
	C_pos <- c(C_pos, sum(asv_relative[[j]]*r$r_pos, na.rm = TRUE))
	C_neg <- c(C_neg, sum(asv_relative[[j]]*r$r_neg, na.rm = TRUE))
	id <- c(id, j)
}
C <- data.frame(id, C_pos, abs(C_neg)/C_pos)
names(C) <- c('sample_id', 'pos_cohesion', 'pos_neg_cohesion')
sub_graph_stat <- merge(sub_graph_stat, C, by = 'sample_id')


############## Statistical analysis

# Calculate the mean and standard deviation
env <- read.csv('env_table.csv', row.names = 1)
env <- merge(env, sub_graph_stat, by = 'sample_id')

Mean <- function(x) mean(x, na.rm = TRUE)
SD <- function(x) sd(x, na.rm = TRUE)
stat <- summaryBy(nodes_num+edges_num+degree+pos_cohesion+pos_neg_cohesion+clustering_coefficient+modularity+complementarity~watermass, data = env, FUN = c(Mean, SD))
stat  # Table 2

# Kruskal-Wallis test
kw <- kruskal(env$nodes_num, env$watermass)
kw$groups
kw <- kruskal(env$edges_num, env$watermass)
kw$groups
kw <- kruskal(env$degree, env$watermass)
kw$groups
kw <- kruskal(env$pos_cohesion, env$watermass)
kw$groups
kw <- kruskal(env$pos_neg_cohesion, env$watermass)
kw$groups
kw <- kruskal(env$clustering_coefficient, env$watermass)
kw$groups
kw <- kruskal(env$modularity, env$watermass)
kw$groups
kw <- kruskal(env$complementarity, env$watermass)
kw$groups

# Mann-Whitney U test (Wilcox test)
compare <- list(c('Plume', 'SCS'), c('Plume', 'Mixed'), c('Mixed', 'SCS'))

ggboxplot(env[c('watermass', 'nodes_num')], x = 'watermass', y = 'nodes_num', color = 'watermass') +
stat_compare_means(comparisons = compare)

ggboxplot(env[c('watermass', 'edges_num')], x = 'watermass', y = 'edges_num', color = 'watermass') +
stat_compare_means(comparisons = compare)

ggboxplot(env[c('watermass', 'degree')], x = 'watermass', y = 'degree', color = 'watermass') +
stat_compare_means(comparisons = compare)

ggboxplot(env[c('watermass', 'pos_cohesion')], x = 'watermass', y = 'pos_cohesion', color = 'watermass') +
stat_compare_means(comparisons = compare)

ggboxplot(env[c('watermass', 'pos_neg_cohesion')], x = 'watermass', y = 'pos_neg_cohesion', color = 'watermass') +
stat_compare_means(comparisons = compare)

ggboxplot(env[c('watermass', 'clustering_coefficient')], x = 'watermass', y = 'clustering_coefficient', color = 'watermass') +
stat_compare_means(comparisons = compare)

ggboxplot(env[c('watermass', 'modularity')], x = 'watermass', y = 'modularity', color = 'watermass') +
stat_compare_means(comparisons = compare)

ggboxplot(env[c('watermass', 'complementarity')], x = 'watermass', y = 'complementarity', color = 'watermass') +
stat_compare_means(comparisons = compare)


############## Export the subnetwork for visualization in Gephi

for (i in sample_id) {
	g <- sub_graph[[i]]
	g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
	edge <- data.frame(as_edgelist(g))
	edge_list <- data.frame(
		source = edge[[1]],
		target = edge[[2]],
		weight = 1
	)
	write.csv(edge_list, paste0(i, '.edge_list.csv'), row.names = FALSE, quote = FALSE)
}

