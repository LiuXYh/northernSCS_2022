library(doBy)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(agricolae)
library(patchwork)
library(ape)
library(picante)
library(hier.part)
library(rdacca.hp)


# Read datasets
tax <- read.csv('taxonomy.csv')
env <- read.csv('env_table.csv', row.names = 2)
tree <- read.tree('asv_seq.tre')

asv <- read.csv('asv_table.csv', row.names = 1)
asv_relative <- asv/rep(colSums(asv),each = nrow(asv))  #relative abundance

# Retain species with a frequency greater than 1 to reduce random errors
asv <- asv[which(apply(asv, 1, function(x) sum(x>0))>1), ]
asv_relative <- asv_relative[rownames(asv), ]

asv <- data.frame(t(asv))
asv_relative <- data.frame(t(asv_relative))

# Prune the phylogenetic tree
asv_id <- intersect(tree$tip.label, colnames(asv_relative))
asv_relative <- asv_relative[asv_id]
tree <- prune.sample(asv_relative, tree)


############## Figure 3A

# Plot
env$watermass <- factor(env$watermass, levels = c('Plume', 'Mixed', 'SCS'))

ggplot(env, aes(watermass, HNA_LNA, fill = watermass)) +
geom_boxplot(width = 0.7, size = 0.5, outlier.size = 1) +
scale_fill_manual(values = c('#607ebc', '#94c128', '#e74946'), limits = c('Plume', 'Mixed', 'SCS')) +
theme(panel.grid = element_blank(), 
	panel.background = element_rect(color = 'black', fill = 'white'), 
	axis.ticks = element_line(color = 'black', size = 0.5), 
	axis.text.y = element_text(color = 'black', size = 9), 
	axis.text.x = element_text(color = 'black', size = 9),
	legend.key = element_blank()) +
labs(x = '', y = 'HNA/LNA ratio')

# Kruskal-Wallis test
kw <- kruskal(env$HNA_LNA, env$watermass)
kw$groups

# Mann-Whitney U test (Wilcox test)
ggboxplot(env[c('watermass', 'HNA_LNA')], x = 'watermass', y = 'HNA_LNA', color = 'watermass') +
stat_compare_means(comparisons = list(c('Plume', 'SCS'), c('Plume', 'Mixed'), c('Mixed', 'SCS')))


############## Figure 3B

## Calculate the abundance-weighted mean strategy values for each ASV

# Abundance-weighted mean salinity (salinity niche)
env <- env[rownames(asv_relative), ]
asv_relative <- data.frame(t(asv_relative))
Salinity_niche_value <- as.data.frame(as.matrix(asv_relative/rowSums(asv_relative)) %*% as.matrix(env['Salinity']))

# Abundance-weighted mean cell-specific bacterial production
env_sBP <- na.omit(env[c('BA', 'sBP')])
asv_relative1 <- asv_relative[rownames(env_sBP)]
sBP_niche_value <- as.data.frame(as.matrix(asv_relative1/rowSums(asv_relative1)) %*% as.matrix(env_sBP['sBP']))

# Add taxonomic annotations (familiy level) for each ASV
# Note: Only the top 10 abundant families are displayed, with the remaining families grouped as 'Others'
tax <- tax[c('asv_id', 'family')]
names(tax) <- c('ASV', 'class')

asv_relative$average <- rowMeans(asv_relative)
asv_relative$ASV <- rownames(asv_relative)
asv_relative <- cbind(asv_relative, Salinity_niche_value, sBP_niche_value)
asv_relative <- merge(asv_relative, tax, by = 'ASV')
asv_relative <- asv_relative[order(asv_relative$average, decreasing = TRUE), ]

Mean <- function(x) mean(x, na.rm = TRUE)
stat <- summaryBy(average~class, asv_relative, FUN = Mean)
stat <- stat[order(stat$average.Mean, decreasing = TRUE), ]
stat <- subset(stat, ! class %in% c('', 'Unassigned'))
top <- stat$class[1:10]
color <- c('#67A6D0', '#F39727', '#0067B0', '#DB0119', '#0AA33E', '#82446D', '#17518D', '#FFCD0A', '#E870A3', '#857BB8', '#9D9A97')
for (i in 1:nrow(asv_relative)) asv_relative[i,'class'] <- ifelse(asv_relative[i,'class'] %in% top, asv_relative[i,'class'], 'Others')
asv_relative$class <- factor(asv_relative$class, levels = rev(c(top, 'Others')))

# Define indicator ASVs based on salinity niche
Plume <- rownames(subset(env, watermass == 'Plume'))
Mixed <- rownames(subset(env, watermass == 'Mixed'))
SCS <- rownames(subset(env, watermass == 'SCS'))

asv_relative$Plume <- apply(asv_relative[Plume], 1, Mean)
asv_relative$Mixed <- apply(asv_relative[Mixed], 1, Mean)
asv_relative$SCS <- apply(asv_relative[SCS], 1, Mean)
asv_relative$ave <- apply(asv_relative[c(Plume, Mixed, SCS)], 1, Mean)

asv_relative[which(asv_relative$Salinity <= 33), 'indicator'] <- 'Plume'
asv_relative[which(asv_relative$Salinity > 33 & asv_relative$Salinity <= 33.75), 'indicator'] <- 'Mixed'
asv_relative[which(asv_relative$Salinity > 33.75), 'indicator'] <- 'SCS'

# Calculate the mean strategy values at the family level based on the families to which the ASVs belong
asv_relative1 <- asv_relative[c('ASV', 'Salinity', 'sBP', 'average', 'class')]
SD <- function(x) sd(x, na.rm = TRUE)
asv_relative1_stat1 <- summaryBy(average~class, asv_relative1, FUN = c(Mean, SD))
asv_relative1_stat2 <- summaryBy(Salinity~class, asv_relative1, FUN = c(Mean, SD))
asv_relative1_stat3 <- summaryBy(sBP~class, asv_relative1, FUN = c(Mean, SD))
asv_relative1_stat <- cbind(asv_relative1_stat1, asv_relative1_stat2[-1], asv_relative1_stat3[-1])

ggplot(asv_relative1_stat, aes(Salinity.Mean, sBP.Mean, color = class)) +
geom_point(aes(size = average.Mean)) + 
#geom_errorbarh(aes(xmin = Salinity.Mean + Salinity.SD, xmax = Salinity.Mean - Salinity.SD)) +
#geom_errorbar(aes(ymin = average.Mean + average.SD, ymax = average.Mean - average.SD)) +
geom_linerange(aes(xmin = Salinity.Mean + Salinity.SD, xmax = Salinity.Mean - Salinity.SD)) +
geom_linerange(aes(ymin = sBP.Mean + sBP.SD, ymax = sBP.Mean - sBP.SD)) +
scale_color_manual(values = rev(color)) +
theme(panel.grid = element_blank(), 
	panel.background = element_rect(color = 'black', fill = 'white'), 
	axis.ticks = element_line(color = 'black', size = 0.5), 
	axis.text.y = element_text(color = 'black', size = 9), 
	axis.text.x = element_text(color = 'black', size = 9),
	legend.key = element_blank()) +
scale_y_continuous(limits = c(0, 5), expand = expansion(mult = c(0, 0))) +
labs(x = 'Salinity_niche_value', y = 'sBP_niche_value') +
geom_vline(xintercept = c(33, 33.75), linetype = 2)


############## Figure S5

asv_relative1 <- asv_relative[c(rownames(env), 'class')]
asv_relative1 <- melt(asv_relative1, id = 'class')

stat1 <- summaryBy(value~class+variable, data = asv_relative1, FUN = Mean)
env$variable <- rownames(env)
env <- merge(env, stat1, by = 'variable')
env <- subset(env, class != 'Others')

stat2 <- summaryBy(Salinity~class, data = asv_relative, FUN = Mean)
stat2 <- stat2[order(stat2$Salinity.Mean), ]
stat2 <- subset(stat2, class != 'Others')
env$class <- factor(env$class, levels = stat2$class)

env$logBP <- log10(env$BP)
env$logsBP <- log10(env$sBP)
env$logBA <- log10(env$BA)

p1 <- ggplot(env, aes(log10(value.Mean+0.00005), logBP)) +
facet_wrap(~class, scale = 'free_x', ncol = 10) +
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
scale_x_continuous(limits = c(-4.5, -1.4), breaks = c(-4, -3, -2)) +
labs(x = 'Relative abundance', y = 'logBP')

p2 <- ggplot(env, aes(log10(value.Mean+0.00005), logBA)) +
facet_wrap(~class, scale = 'free_x', ncol = 10) +
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
scale_x_continuous(limits = c(-4.5, -1.4), breaks = c(-4, -3, -2)) +
labs(x = 'Relative abundance', y = 'logBA')

p3 <- ggplot(env, aes(log10(value.Mean+0.00005), logsBP)) +
facet_wrap(~class, scale = 'free_x', ncol = 10) +
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
scale_x_continuous(limits = c(-4.5, -1.4), breaks = c(-4, -3, -2)) +
labs(x = 'Relative abundance', y = 'logsBP')

layout <- c(
	area(1, 1, 1, 10),
	area(2, 1, 2, 10),
	area(3, 1, 3, 10)
)

p <- p1+p2+p3 + plot_layout(design = layout)
p


############## Figure 3C

# Select indicator ASVs 
Plume <- asv_relative[which(asv_relative$indicator == 'Plume'),'ASV']
SCS <- asv_relative[which(asv_relative$indicator == 'SCS'),'ASV']
Mixed <- asv_relative[which(asv_relative$indicator == 'Mixed'),'ASV']

asv_Plume <- asv[Plume]
asv_SCS <- asv[SCS]
asv_Mixed <- asv[Mixed]

# Calculate species diversity and phylogenetic diversity for Plume indicators
env <- read.csv('env_table.csv', row.names = 2)
tree <- read.tree('asv_seq.tre')
asv_id <- intersect(tree$tip.label, colnames(asv_Plume))
asv_Plume <- asv_Plume[asv_id]
tree <- prune.sample(asv_Plume, tree)
dis <- cophenetic(tree)

set.seed(123)
mpd_Plume <- ses.mpd(asv_Plume, dis, abundance.weighted = TRUE, null.model = 'taxa.labels', runs = 999)
mpd_Plume$richness1 <- mpd_Plume$ntaxa
mpd_Plume$mpd1 <- mpd_Plume$mpd.obs
mpd_Plume$SESmpd1 <- mpd_Plume$mpd.obs.z
env <- env[rownames(mpd_Plume), ]
env <- cbind(env, mpd_Plume)
env$shannon1 <- diversity(asv_Plume, index = 'shannon', base = exp(1))

# Calculate species diversity and phylogenetic diversity for SCS indicators
tree <- read.tree('asv_seq.tre')
asv_id <- intersect(tree$tip.label, colnames(asv_SCS))
asv_SCS <- asv_SCS[asv_id]
tree <- prune.sample(asv_SCS, tree)
dis <- cophenetic(tree)

set.seed(123)
mpd_SCS <- ses.mpd(asv_SCS, dis, abundance.weighted = TRUE, null.model = 'taxa.labels', runs = 999)
mpd_SCS$richness2 <- mpd_SCS$ntaxa
mpd_SCS$mpd2 <- mpd_SCS$mpd.obs
mpd_SCS$SESmpd2 <- mpd_SCS$mpd.obs.z
env <- cbind(env, mpd_SCS)
env$shannon2 <- diversity(asv_SCS, index = 'shannon', base = exp(1))

# Calculate species diversity and phylogenetic diversity for Mixed indicators
tree <- read.tree('asv_seq.tre')
asv_id <- intersect(tree$tip.label, colnames(asv_Mixed))
asv_Mixed <- asv_Mixed[asv_id]
tree <- prune.sample(asv_Mixed, tree)
dis <- cophenetic(tree)

set.seed(123)
mpd_Mixed <- ses.mpd(asv_Mixed, dis, abundance.weighted = TRUE, null.model = 'taxa.labels', runs = 999)
mpd_Mixed$richness3 <- mpd_Mixed$ntaxa
mpd_Mixed$mpd3 <- mpd_Mixed$mpd.obs
mpd_Mixed$SESmpd3 <- mpd_Mixed$mpd.obs.z
env <- cbind(env, mpd_Mixed)
env$shannon3 <- diversity(asv_Mixed, index = 'shannon', base = exp(1))

# The contribution of these indicators to community diversity was then determined using hierarchical partitioning
env$sample_id <- rownames(env)
mpd <- env[c('sample_id', 'mpd', 'mpd1', 'mpd2', 'mpd3')]
names(mpd) <- c('sample_id', 'div', 'Plume', 'SCS', 'Mixed')
ses_mpd <- env[c('sample_id', 'SESmpd', 'SESmpd1', 'SESmpd2', 'SESmpd3')]
names(ses_mpd) <- c('sample_id', 'div', 'Plume', 'SCS', 'Mixed')
richness <- env[c('sample_id', 'richness', 'richness1', 'richness2', 'richness3')]
names(richness) <- c('sample_id', 'div', 'Plume', 'SCS', 'Mixed')
shannon <- env[c('sample_id', 'shannon', 'shannon1', 'shannon2', 'shannon3')]
names(shannon) <- c('sample_id', 'div', 'Plume', 'SCS', 'Mixed')

par(mfrow = c(2, 2))
richness_part <- hier.part(y = richness$div, xcan = richness[c('Plume', 'SCS', 'Mixed')], family = 'gaussian', gof = 'Rsqu', barplot = TRUE)
shannon_part <- hier.part(y = shannon$div, xcan = shannon[c('Plume', 'SCS', 'Mixed')], family = 'gaussian', gof = 'Rsqu', barplot = TRUE)
mpd_part <- hier.part(y = mpd$div, xcan = mpd[c('Plume', 'SCS', 'Mixed')], family = 'gaussian', gof = 'Rsqu', barplot = TRUE)
ses_mpd_part <- hier.part(y = ses_mpd$div, xcan = ses_mpd[c('Plume', 'SCS', 'Mixed')], family = 'gaussian', gof = 'Rsqu', barplot = TRUE)

# Pie plot
richness_part <- rdacca.hp(richness['div'], richness[c('Plume', 'SCS', 'Mixed')], type = 'adjR2', scale = FALSE)
shannon_part <- rdacca.hp(shannon['div'], shannon[c('Plume', 'SCS', 'Mixed')], type = 'adjR2', scale = FALSE)
mpd_part <- rdacca.hp(mpd['div'], mpd[c('Plume', 'SCS', 'Mixed')], type = 'adjR2', scale = FALSE)
ses_mpd_part <- rdacca.hp(ses_mpd['div'], ses_mpd[c('Plume', 'SCS', 'Mixed')], type = 'adjR2', scale = FALSE)

richness_part <- richness_part$Hier.part[,4]
shannon_part <- shannon_part$Hier.part[,4]
mpd_part <- mpd_part$Hier.part[,4]
ses_mpd_part <- ses_mpd_part$Hier.part[,4]

div <- data.frame(richness_part, shannon_part, mpd_part, ses_mpd_part)
div$watermass <- rownames(div)
div <- melt(div, id = 'watermass')
div$watermass <- factor(div$watermass, levels = c('Plume', 'Mixed', 'SCS'))

ggplot(div, aes(x = '', y = value, fill = watermass)) +
facet_wrap(~variable) +
geom_bar(stat = 'identity') +
coord_polar(theta = 'y') +
scale_fill_manual(values = c('#607ebc', '#94c128', '#e74946'), limits = c('Plume', 'Mixed', 'SCS')) +
theme(panel.grid = element_blank(), 
	panel.background = element_rect(color = 'black', fill = 'white'), 
	axis.text = element_blank(),
	axis.ticks = element_blank(),
	legend.key = element_blank()) +
scale_y_continuous(expand = c(0, 0)) +
labs(x = '', y = '')




