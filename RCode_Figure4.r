library(iCAMP)
library(spaa)
library(ape)
library(picante)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(patchwork)


# Read datasets
env <- read.csv('env_table.csv', row.names = 2)
tree <- read.tree('asv_seq.tre')

asv <- read.csv('asv_table.csv', row.names = 1)
asv_relative <- asv/rep(colSums(asv),each = nrow(asv))  #relative abundance

# Retain species with a frequency greater than 1 to reduce random errors
asv_relative <- asv_relative[which(apply(asv_relative, 1, function(x) sum(x>0))>1), ]
asv_relative <- data.frame(t(asv_relative))

# Prune the phylogenetic tree
asv_id <- intersect(tree$tip.label, colnames(asv_relative))
asv_relative <- asv_relative[asv_id]
tree <- prune.sample(asv_relative, tree)


############## Calculate niche differences and niche overlap

env <- env[rownames(asv_relative), ]
niche.dif <- dniche(env[c('Salinity')], asv_relative, method = 'niche.value', nworker = 6, out.dist = FALSE, bigmemo = FALSE)
niche_diff <- niche.dif$nd$Salinity

niche_overlap <- niche.overlap(asv_relative, method = 'levins')
niche_overlap <- as.matrix(niche_overlap)


############## Bin the phylogenetic tree and calculate the average phylogenetic distance and niche indices for species within each bin

env$logBP <- log10(env$BP)
env$logsBP <- log10(env$sBP)
env$logBA <- log10(env$BA)

asv_relative <- data.frame(t(asv_relative))
niche_diff <- as.data.frame(niche_diff)
niche_overlap <- as.data.frame(niche_overlap)

asv_id <- intersect(intersect(intersect(tree$tip.label, rownames(niche_overlap)), rownames(niche_diff)), rownames(asv_relative))
niche_overlap <- niche_overlap[asv_id,asv_id]
niche_diff <- niche_diff[asv_id,asv_id]
tree <- prune.sample(niche_overlap, tree)
tree_dis <- as.data.frame(as.matrix(cophenetic(tree)))
tree_dis <- tree_dis[asv_id,asv_id]

niche_overlap <- as.matrix(niche_overlap)
diag(niche_overlap) <- NA
niche_overlap <- as.data.frame(niche_overlap)
niche_diff <- as.matrix(niche_diff)
diag(niche_diff) <- NA
niche_diff <- as.data.frame(niche_diff)
tree_dis <- as.matrix(tree_dis)
diag(tree_dis) <- NA
tree_dis <- as.data.frame(tree_dis)

niche_overlap$id <- rownames(niche_overlap)
niche_overlap <- melt(niche_overlap, by = 'id')
names(niche_overlap) <- c('source', 'target', 'niche_overlap')
niche_diff$id <- rownames(niche_diff)
niche_diff <- melt(niche_diff, by = 'id')
names(niche_diff) <- c('source', 'target', 'niche_diff')
tree_dis$id <- rownames(tree_dis)
tree_dis <- melt(tree_dis, by = 'id')
names(tree_dis) <- c('source', 'target', 'tree_dis1')

dis <- cbind(tree_dis, niche_overlap['niche_overlap'], niche_diff['niche_diff'])
dis <- dis[order(dis$tree_dis1), ]

result <- NULL
for (i in seq(min(dis$tree_dis1, na.rm = TRUE), max(dis$tree_dis1, na.rm = TRUE), 0.1)) {
	dis_i <- subset(dis, tree_dis1 >= i & tree_dis1 < (i+0.1))
	asv_i <- unique(c(as.character(dis_i$source), as.character(dis_i$target)))
	otu_i <- asv_relative[asv_i, ]
	abundance_mean <- colMeans(otu_i, na.rm = TRUE)
	
	# Average phylogenetic distance, average niche differences, average niche overlap; average relative abundance and its correlation with BP, BA, and sBP
	result <- rbind(result, c(i, mean(dis_i$tree_dis1, na.rm = TRUE), mean(dis_i$niche_overlap, na.rm = TRUE), mean(dis_i$niche_diff, na.rm = TRUE), mean(abundance_mean, na.rm = TRUE), mean(cor(log10(abundance_mean[!is.na(env$logBA)]), env$logBA[!is.na(env$logBA)], method = 'pearson'), na.rm = TRUE), mean(cor(log10(abundance_mean), env$logBP, method = 'pearson'), na.rm = TRUE), mean(cor(log10(abundance_mean[!is.na(env$logsBP)]), env$logsBP[!is.na(env$logsBP)], method = 'pearson'), na.rm = TRUE)))
}

result <- as.data.frame(result)
names(result) <- c('tree', 'mean_tree_dis', 'mean_niche_overlap', 'mean_niche_diff', 'mean_abundance', 'cor_abundance_BA', 'cor_abundance_BP', 'cor_abundance_sBP')


# Plot Figure 4
cor.test(result$mean_tree_dis, result$mean_niche_overlap)
cor.test(result$mean_tree_dis, result$mean_niche_diff)
cor.test(result$mean_tree_dis, result$cor_abundance_BP)
cor.test(result$mean_tree_dis, result$cor_abundance_BA)
cor.test(result$mean_tree_dis, result$cor_abundance_sBP)

result1 <- melt(result[c('tree', 'mean_tree_dis', 'mean_niche_overlap', 'mean_niche_diff')], id = c('tree', 'mean_tree_dis'))
result2 <- melt(result[c('tree', 'mean_tree_dis', 'cor_abundance_BA', 'cor_abundance_BP', 'cor_abundance_sBP')], id = c('tree', 'mean_tree_dis'))

p1 <- ggplot(result1, aes(mean_tree_dis, value, color = variable)) +
geom_point() +
theme(panel.grid = element_blank(), 
	panel.background = element_rect(color = 'black', fill = 'white'), 
	axis.ticks = element_line(color = 'black', size = 0.5), 
	axis.text.y = element_text(color = 'black', size = 9), 
	axis.text.x = element_text(color = 'black', size = 9),
	legend.key = element_blank()) +
stat_smooth(method = 'lm', formula = y~poly(x,1), se = TRUE) +
stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'pearson', label.x.npc = 'left', label.y.npc = 'top', size = 2.7) +
#stat_poly_eq(aes(label = paste(..rr.label.., stat(p.value.label), sep = '~`,`~')), formula = y~poly(x, 2), parse = TRUE, label.x.npc = 'left', label.y.npc = 'top', size = 2.7) +
scale_y_continuous(limits = c(0, 0.4), expand = expansion(mult = c(0, 0))) +
labs(x = 'Phylogenetic distance', y = '')


p2 <- ggplot(result2, aes(mean_tree_dis, value, color = variable)) +
geom_point() +
theme(panel.grid = element_blank(), 
	panel.background = element_rect(color = 'black', fill = 'white'), 
	axis.ticks = element_line(color = 'black', size = 0.5), 
	axis.text.y = element_text(color = 'black', size = 9), 
	axis.text.x = element_text(color = 'black', size = 9),
	legend.key = element_blank()) +
stat_smooth(method = 'lm', formula = y~poly(x,1), se = TRUE) +
stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'pearson', label.x.npc = 'left', label.y.npc = 'top', size = 2.7) +
#stat_poly_eq(aes(label = paste(..rr.label.., stat(p.value.label), sep = '~`,`~')), formula = y~poly(x, 2), parse = TRUE, label.x.npc = 'left', label.y.npc = 'top', size = 2.7) +
scale_y_continuous(limits = c(-1, 1), expand = expansion(mult = c(0, 0))) +
labs(x = 'Phylogenetic distance', y = '')

layout <- c(
	area(1, 1, 1, 1),
	area(2, 1, 2, 1)
)

p1 + p2 + plot_layout(design = layout)


