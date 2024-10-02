library(reshape2)
library(doBy)
library(ggplot2)
library(agricolae)
library(ggpubr)
library(vegan)
library(picante)
library(ape)
library(GUniFrac)
library(plspm)


# Read datasets
env <- read.csv('env_table.csv', row.names = 2)
asv <- read.csv('asv_table.csv', row.names = 1)
tree <- read.tree('asv_seq.tre')

# Retain species with a frequency greater than 1 to reduce random errors
asv <- asv[which(apply(asv, 1, function(x) sum(x>0))>1), ]
asv <- data.frame(t(asv))

# Prune the phylogenetic tree
asv_id <- intersect(tree$tip.label, colnames(asv))
asv <- asv[asv_id]
tree <- prune.sample(asv, tree)

# Calculate species diversity
env <- env[rownames(asv), ]
env$richness <- rowSums(asv > 0)
env$shannon <- diversity(asv, index = 'shannon', base = exp(1))

# Calculate phylogenetic diversity
dis <- cophenetic(tree)
set.seed(123)
mpd <- ses.mpd(asv, dis, abundance.weighted = TRUE, null.model = 'taxa.labels', runs = 999)
env$mpd <- mpd$mpd.obs
env$SESmpd <- mpd$mpd.obs.z


############## Table 1

# Calculate the mean and standard deviation
Mean <- function(x) mean(x, na.rm = TRUE)
SD <- function(x) sd(x, na.rm = TRUE)
stat <- summaryBy(BP+BA+sBP+richness+shannon+mpd+SESmpd~watermass, data = env, FUN = c(Mean, SD))
stat  # Table 1

# Kruskal-Wallis test
kw <- kruskal(env$BP, env$watermass)
kw$groups
kw <- kruskal(env$BA, env$watermass)
kw$groups
kw <- kruskal(env$sBP, env$watermass)
kw$groups
kw <- kruskal(env$richness, env$watermass)
kw$groups
kw <- kruskal(env$shannon, env$watermass)
kw$groups
kw <- kruskal(env$mpd, env$watermass)
kw$groups
kw <- kruskal(env$SESmpd, env$watermass)
kw$groups

# Mann-Whitney U test (Wilcox test)
compare <- list(c('Plume', 'SCS'), c('Plume', 'Mixed'), c('Mixed', 'SCS'))

ggboxplot(na.omit(env[c('watermass', 'BP')]), x = 'watermass', y = 'BP', color = 'watermass') +
stat_compare_means(comparisons = compare)

ggboxplot(na.omit(env[c('watermass', 'BA')]), x = 'watermass', y = 'BA', color = 'watermass') +
stat_compare_means(comparisons = compare)

ggboxplot(na.omit(env[c('watermass', 'sBP')]), x = 'watermass', y = 'sBP', color = 'watermass') +
stat_compare_means(comparisons = compare)

ggboxplot(na.omit(env[c('watermass', 'richness')]), x = 'watermass', y = 'richness', color = 'watermass') +
stat_compare_means(comparisons = compare)

ggboxplot(na.omit(env[c('watermass', 'shannon')]), x = 'watermass', y = 'shannon', color = 'watermass') +
stat_compare_means(comparisons = compare)

ggboxplot(na.omit(env[c('watermass', 'mpd')]), x = 'watermass', y = 'mpd', color = 'watermass') +
stat_compare_means(comparisons = compare)

ggboxplot(na.omit(env[c('watermass', 'SESmpd')]), x = 'watermass', y = 'SESmpd', color = 'watermass') +
stat_compare_means(comparisons = compare)


############## Figure S3

env$logBP <- log10(env$BP)
env$logsBP <- log10(env$sBP)
env$logBA <- log10(env$BA)

env1 <- env[c('watermass', 'richness', 'shannon', 'mpd', 'SESmpd', 'logBP', 'logsBP', 'logBA')]
env1 <- melt(env1, id = c('watermass', 'richness', 'shannon', 'mpd', 'SESmpd'))
env1 <- melt(env1, id = c('watermass', 'variable', 'value'))
names(env1) <- c('watermass', 'variable1', 'value1', 'variable2', 'value2')

ggplot(env1, aes(value2, value1)) +
facet_grid(variable1~variable2, scale = 'free') +
geom_point(aes(color = watermass)) +
scale_color_manual(values = c('#607ebc', '#94c128', '#e74946'), limits = c('Plume', 'Mixed', 'SCS')) +
theme(panel.grid = element_blank(), 
	panel.background = element_rect(color = 'black', fill = 'white'), 
	axis.ticks = element_line(color = 'black', size = 0.5), 
	axis.text.y = element_text(color = 'black', size = 9), 
	axis.text.x = element_text(color = 'black', size = 9),
	legend.key = element_blank()) +
stat_smooth(method = 'lm', formula = y~poly(x,1), se = TRUE) +
stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'pearson', label.x.npc = 'left', label.y.npc = 'top', size = 2.7) +
#stat_poly_eq(aes(label = paste(..rr.label.., stat(p.value.label), sep = '~`,`~')), formula = y~poly(x, 2), parse = TRUE, label.x.npc = 'left', label.y.npc = 'top', size = 2.7) +
labs(x = '', y = '', color = '')


############## Figure 2A

# Calculate weighted Unifrac distance
unifrac <- GUniFrac(asv, tree)$unifracs
wei_unif_dis <- unifrac[, , 'd_1']

# Distance-based redundancy analysis (db-RDA)
env1 <- na.omit(env[c('watermass', 'Temperature', 'Salinity', 'Chla', 'BP', 'sBP', 'BA', 'HNA_LNA', 'richness', 'shannon', 'mpd', 'SESmpd')])

sample_id <- intersect(rownames(env1), rownames(wei_unif_dis))
env1 <- env1[sample_id, ]
wei_unif_dis <- wei_unif_dis[sample_id,sample_id]

wei_unif_dis <- as.dist(wei_unif_dis)
db_rda <- dbrda(wei_unif_dis~., env1[c('Temperature', 'Salinity', 'Chla', 'BP', 'sBP', 'BA', 'HNA_LNA', 'richness', 'shannon', 'mpd', 'SESmpd')], add = FALSE)

exp_adj <- RsquareAdj(db_rda)$adj.r.squared * db_rda$CCA$eig/sum(db_rda$CCA$eig)
RDA1_exp <- paste('db-RDA1 (', round(exp_adj[1]*100, 2), '% )')
RDA2_exp <- paste('db-RDA2 (', round(exp_adj[2]*100, 2), '% )')

db_rda.scaling2 <- summary(db_rda)
db_rda_site.scaling2 <- data.frame(db_rda.scaling2$site[ ,1:2])
db_rda_site.scaling2$plot <- rownames(db_rda_site.scaling2)
db_rda_env.scaling2 <- data.frame(db_rda.scaling2$biplot[ ,1:2])
db_rda_env.scaling2$name <- rownames(db_rda_env.scaling2)
db_rda_site.scaling2 <- cbind(db_rda_site.scaling2, env1)

ggplot(db_rda_site.scaling2, aes(dbRDA1, dbRDA2)) +
geom_point(aes(color = watermass), pch = 16) +
scale_color_manual(values = c('#607ebc', '#94c128', '#e74946'), limits = c('Plume', 'Mixed', 'SCS')) +
stat_ellipse(aes(fill = watermass), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
scale_fill_manual(values = c('#607ebc', '#94c128', '#e74946'), limits = c('Plume', 'Mixed', 'SCS')) +
theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5), 
	legend.key = element_blank(),  
    panel.background = element_rect(color = 'black', fill = 'transparent')) +
labs(x = RDA1_exp, y = RDA2_exp, color = '') +
geom_segment(data = db_rda_env.scaling2, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), color = 'blue', size = 0.3) +
geom_text(data = db_rda_env.scaling2, aes(x = dbRDA1, y = dbRDA2, label = name), color = 'blue', size = 3) +
geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
geom_hline(yintercept = 0, color = 'gray', size = 0.5)


############## Figure 2B and 2C

# Path analysis diagram from partial least squares path modeling (PLS-PM)
env2 <- data.frame(scale(na.omit(env[c('Temperature', 'Salinity', 'Chla', 'BA', 'BP', 'SESmpd')])))

dat_blocks <- list(
    Temperature = c('Temperature'), 
	Salinity = c('Salinity'), 
	Chla = c('Chla'), 
    PD = c('SESmpd'), 
    productivity = c('BA', 'BP')
)

Temperature <- c(0, 0, 0, 0, 0)
Salinity <- c(0, 0, 0, 0, 0)
Chla <- c(1, 1, 0, 0, 0)
PD <- c(1, 1, 1, 0, 0)
productivity <- c(1, 1, 1, 1, 0)
dat_path <- rbind(Temperature, Salinity, Chla, PD, productivity)
colnames(dat_path) <- rownames(dat_path)

dat_pls <- plspm(env2, dat_path, dat_blocks, modes = rep('A', 5))
summary(dat_pls)

# coefs
dat_pls$path_coefs

# Path
innerplot(dat_pls, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray', box.lwd = 0)
dat_pls$inner_model

# GOF
dat_pls$gof

# The direct and indirect effects
dat_pls$effects

# Plot effect
effect <- data.frame(dat_pls$effects, stringsAsFactors = FALSE)
effect <- effect[grepl('productivity', effect$relationships), ]
effect
effect <- reshape2::melt(effect, id = c('relationships', 'total'))
effect <- effect[order(effect$total, decreasing = TRUE), ]
effect$relationships <- factor(effect$relationships, levels = unique(effect$relationships))

ggplot(effect, aes(relationships, value)) +
#coord_flip() +
geom_col(aes(fill = variable), position = 'stack', width = 0.55) +
theme(panel.grid = element_blank(), panel.background = element_blank(),
    axis.line = element_line(colour = 'black')) +
labs(x = '', y = 'effects', fill = '') +
geom_hline(yintercept = 0)

