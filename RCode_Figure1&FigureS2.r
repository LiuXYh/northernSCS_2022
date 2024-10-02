library(reshape2)
library(doBy)
library(ggplot2)
library(maptools)
library(MBA)


############## Figure S2

#### Salinity

# Spatial interpolation
env <- read.csv('env_table.csv')
names(env)[which(names(env) == 'Salinity')] <- 'value'

surface <- subset(env, layer == 'surface')
surface <- mba.surf(na.omit(surface[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(surface$xyz.est$z) <- list(surface$xyz.est$x, surface$xyz.est$y)
surface <- reshape2::melt(surface$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
surface <- na.omit(surface)
surface$layer <- 'surface'

middle <- subset(env, layer == 'middle')
middle <- mba.surf(na.omit(middle[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(middle$xyz.est$z) <- list(middle$xyz.est$x, middle$xyz.est$y)
middle <- reshape2::melt(middle$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
middle <- na.omit(middle)
middle$layer <- 'middle'

bottom <- subset(env, layer == 'bottom')
bottom <- mba.surf(na.omit(bottom[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(bottom$xyz.est$z) <- list(bottom$xyz.est$x, bottom$xyz.est$y)
bottom <- reshape2::melt(bottom$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
bottom <- na.omit(bottom)
bottom$layer <- 'bottom'

env2 <- rbind(surface, middle, bottom)
env2$layer <- factor(env2$layer, levels = c('surface', 'middle', 'bottom'))
env$layer <- factor(env$layer, levels = c('surface', 'middle', 'bottom'))

# Plot map
ggplot() +
facet_wrap(~layer, ncol = 1, scale = 'free') +
geom_raster(data = env2, aes(x = longitude, y = latitude, fill = value), interpolate = TRUE) + 
scale_fill_gradientn(colors = c(colorRampPalette(c('#000084', '#004FFF'))(7), '#0EFFF0', '#8DFF70', '#F8FF05', '#FF6900', '#8D0000'), limits = c(21, 35)) +
geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group),fill = 'gray', col = 'gray') +
geom_point(data = env, aes(x = longitude, y = latitude, shape = watermass), size = 1) +
scale_shape_manual(values = c(4, 13, 1), limits = c('Plume', 'Mixed', 'SCS')) +
theme_bw() + 
theme(panel.grid.minor = element_blank(), legend.background = element_blank()) +
labs(x = 'Longitude (°E)', y = 'Latitude (°N)', fill = 'Salinity', shape = 'Watermass') +
coord_cartesian(xlim = c(113.5, 118.5), ylim = c(20.5, 23.5))


#### Temperature

# Spatial interpolation
env <- read.csv('env_table.csv')
names(env)[which(names(env) == 'Temperature')] <- 'value'

surface <- subset(env, layer == 'surface')
surface <- mba.surf(na.omit(surface[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(surface$xyz.est$z) <- list(surface$xyz.est$x, surface$xyz.est$y)
surface <- reshape2::melt(surface$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
surface <- na.omit(surface)
surface$layer <- 'surface'

middle <- subset(env, layer == 'middle')
middle <- mba.surf(na.omit(middle[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(middle$xyz.est$z) <- list(middle$xyz.est$x, middle$xyz.est$y)
middle <- reshape2::melt(middle$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
middle <- na.omit(middle)
middle$layer <- 'middle'

bottom <- subset(env, layer == 'bottom')
bottom <- mba.surf(na.omit(bottom[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(bottom$xyz.est$z) <- list(bottom$xyz.est$x, bottom$xyz.est$y)
bottom <- reshape2::melt(bottom$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
bottom <- na.omit(bottom)
bottom$layer <- 'bottom'

env2 <- rbind(surface, middle, bottom)
env2$layer <- factor(env2$layer, levels = c('surface', 'middle', 'bottom'))
env$layer <- factor(env$layer, levels = c('surface', 'middle', 'bottom'))

# Plot map
ggplot() +
facet_wrap(~layer, ncol = 1, scale = 'free') +
geom_raster(data = env2, aes(x = longitude, y = latitude, fill = value), interpolate = TRUE) + 
scale_fill_gradientn(colors = c('#000084', '#004FFF', '#0EFFF0', '#8DFF70', '#F8FF05', '#FF6900', '#8D0000'), limits = c(17, 32)) + 
geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group),fill = 'gray', col = 'gray') +
geom_point(data = env, aes(x = longitude, y = latitude, shape = watermass), size = 1) +
scale_shape_manual(values = c(4, 13, 1), limits = c('Plume', 'Mixed', 'SCS')) +
theme_bw() + 
theme(panel.grid.minor = element_blank(), legend.background = element_blank()) +
labs(x = 'Longitude (°E)', y = 'Latitude (°N)', fill = 'Temperature', shape = 'Watermass') +
coord_cartesian(xlim = c(113.5, 118.5), ylim = c(20.5, 23.5))


#### Chla

# Spatial interpolation
env <- read.csv('env_table.csv')
names(env)[which(names(env) == 'Chla')] <- 'value'

surface <- subset(env, layer == 'surface')
surface <- mba.surf(na.omit(surface[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(surface$xyz.est$z) <- list(surface$xyz.est$x, surface$xyz.est$y)
surface <- reshape2::melt(surface$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
surface <- na.omit(surface)
surface$layer <- 'surface'

middle <- subset(env, layer == 'middle')
middle <- mba.surf(na.omit(middle[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(middle$xyz.est$z) <- list(middle$xyz.est$x, middle$xyz.est$y)
middle <- reshape2::melt(middle$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
middle <- na.omit(middle)
middle$layer <- 'middle'

bottom <- subset(env, layer == 'bottom')
bottom <- mba.surf(na.omit(bottom[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(bottom$xyz.est$z) <- list(bottom$xyz.est$x, bottom$xyz.est$y)
bottom <- reshape2::melt(bottom$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
bottom <- na.omit(bottom)
bottom$layer <- 'bottom'

env2 <- rbind(surface, middle, bottom)
env2$layer <- factor(env2$layer, levels = c('surface', 'middle', 'bottom'))
env$layer <- factor(env$layer, levels = c('surface', 'middle', 'bottom'))

# Plot map
ggplot() +
facet_wrap(~layer, ncol = 1, scale = 'free') +
geom_raster(data = env2, aes(x = longitude, y = latitude, fill = value), interpolate = TRUE) + 
scale_fill_gradientn(colors = c('#000084', '#004FFF', '#0EFFF0', '#8DFF70', '#F8FF05', colorRampPalette(c('#FF6900', '#8D0000'))(6)), limits = c(0, 4.5)) + 
geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group),fill = 'gray', col = 'gray') +
geom_point(data = env, aes(x = longitude, y = latitude, shape = watermass), size = 1) +
scale_shape_manual(values = c(4, 13, 1), limits = c('Plume', 'Mixed', 'SCS')) +
theme_bw() + 
theme(panel.grid.minor = element_blank(), legend.background = element_blank()) +
labs(x = 'Longitude (°E)', y = 'Latitude (°N)', fill = 'Chla', shape = 'Watermass') +
coord_cartesian(xlim = c(113.5, 118.5), ylim = c(20.5, 23.5))


#### BP

# Spatial interpolation
env <- read.csv('env_table.csv')
names(env)[which(names(env) == 'BP')] <- 'value'

surface <- subset(env, layer == 'surface')
surface <- mba.surf(na.omit(surface[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(surface$xyz.est$z) <- list(surface$xyz.est$x, surface$xyz.est$y)
surface <- reshape2::melt(surface$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
surface <- na.omit(surface)
surface$layer <- 'surface'

middle <- subset(env, layer == 'middle')
middle <- mba.surf(na.omit(middle[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(middle$xyz.est$z) <- list(middle$xyz.est$x, middle$xyz.est$y)
middle <- reshape2::melt(middle$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
middle <- na.omit(middle)
middle$layer <- 'middle'

bottom <- subset(env, layer == 'bottom')
bottom <- mba.surf(na.omit(bottom[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(bottom$xyz.est$z) <- list(bottom$xyz.est$x, bottom$xyz.est$y)
bottom <- reshape2::melt(bottom$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
bottom <- na.omit(bottom)
bottom$layer <- 'bottom'

env2 <- rbind(surface, middle, bottom)
env2$layer <- factor(env2$layer, levels = c('surface', 'middle', 'bottom'))
env$layer <- factor(env$layer, levels = c('surface', 'middle', 'bottom'))

# Plot map
ggplot() +
facet_wrap(~layer, ncol = 1, scale = 'free') +
geom_raster(data = env2, aes(x = longitude, y = latitude, fill = value), interpolate = TRUE) + 
scale_fill_gradientn(colors = c('#000084', '#004FFF', '#0EFFF0', '#8DFF70', '#F8FF05', colorRampPalette(c('#FF6900', '#8D0000'))(6)), limits = c(0, 7)) + 
geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group),fill = 'gray', col = 'gray') +
geom_point(data = env, aes(x = longitude, y = latitude, shape = watermass), size = 1) +
scale_shape_manual(values = c(4, 13, 1), limits = c('Plume', 'Mixed', 'SCS')) +
theme_bw() + 
theme(panel.grid.minor = element_blank(), legend.background = element_blank()) +
labs(x = 'Longitude (°E)', y = 'Latitude (°N)', fill = 'BP', shape = 'Watermass') +
coord_cartesian(xlim = c(113.5, 118.5), ylim = c(20.5, 23.5))


#### sBP

# Spatial interpolation
env <- read.csv('env_table.csv')
names(env)[which(names(env) == 'sBP')] <- 'value'

surface <- subset(env, layer == 'surface')
surface <- mba.surf(na.omit(surface[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(surface$xyz.est$z) <- list(surface$xyz.est$x, surface$xyz.est$y)
surface <- reshape2::melt(surface$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
surface <- na.omit(surface)
surface$layer <- 'surface'

middle <- subset(env, layer == 'middle')
middle <- mba.surf(na.omit(middle[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(middle$xyz.est$z) <- list(middle$xyz.est$x, middle$xyz.est$y)
middle <- reshape2::melt(middle$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
middle <- na.omit(middle)
middle$layer <- 'middle'

bottom <- subset(env, layer == 'bottom')
bottom <- mba.surf(na.omit(bottom[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(bottom$xyz.est$z) <- list(bottom$xyz.est$x, bottom$xyz.est$y)
bottom <- reshape2::melt(bottom$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
bottom <- na.omit(bottom)
bottom$layer <- 'bottom'

env2 <- rbind(surface, middle, bottom)
env2$layer <- factor(env2$layer, levels = c('surface', 'middle', 'bottom'))
env$layer <- factor(env$layer, levels = c('surface', 'middle', 'bottom'))

# Plot map
ggplot() +
facet_wrap(~layer, ncol = 1, scale = 'free') +
geom_raster(data = env2, aes(x = longitude, y = latitude, fill = value), interpolate = TRUE) + 
scale_fill_gradientn(colors = c('#000084', '#004FFF', '#0EFFF0', '#8DFF70', '#F8FF05', colorRampPalette(c('#FF6900', '#8D0000'))(6)), limits = c(0, 11)) + 
geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group),fill = 'gray', col = 'gray') +
geom_point(data = env, aes(x = longitude, y = latitude, shape = watermass), size = 1) +
scale_shape_manual(values = c(4, 13, 1), limits = c('Plume', 'Mixed', 'SCS')) +
theme_bw() + 
theme(panel.grid.minor = element_blank(), legend.background = element_blank()) +
labs(x = 'Longitude (°E)', y = 'Latitude (°N)', fill = 'sBP', shape = 'Watermass') +
coord_cartesian(xlim = c(113.5, 118.5), ylim = c(20.5, 23.5))


#### BA

# Spatial interpolation
env <- read.csv('env_table.csv')
names(env)[which(names(env) == 'BA')] <- 'value'

surface <- subset(env, layer == 'surface')
surface <- mba.surf(na.omit(surface[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(surface$xyz.est$z) <- list(surface$xyz.est$x, surface$xyz.est$y)
surface <- reshape2::melt(surface$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
surface <- na.omit(surface)
surface$layer <- 'surface'

middle <- subset(env, layer == 'middle')
middle <- mba.surf(na.omit(middle[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(middle$xyz.est$z) <- list(middle$xyz.est$x, middle$xyz.est$y)
middle <- reshape2::melt(middle$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
middle <- na.omit(middle)
middle$layer <- 'middle'

bottom <- subset(env, layer == 'bottom')
bottom <- mba.surf(na.omit(bottom[c('longitude', 'latitude', 'value')]), no.X = 50, no.Y = 50, extend = FALSE)
dimnames(bottom$xyz.est$z) <- list(bottom$xyz.est$x, bottom$xyz.est$y)
bottom <- reshape2::melt(bottom$xyz.est$z, varnames = c('longitude', 'latitude'), value.name = 'value')
bottom <- na.omit(bottom)
bottom$layer <- 'bottom'

env2 <- rbind(surface, middle, bottom)
env2$layer <- factor(env2$layer, levels = c('surface', 'middle', 'bottom'))
env$layer <- factor(env$layer, levels = c('surface', 'middle', 'bottom'))

# Plot map
ggplot() +
facet_wrap(~layer, ncol = 1, scale = 'free') +
geom_raster(data = env2, aes(x = longitude, y = latitude, fill = value), interpolate = TRUE) + 
scale_fill_gradientn(colors = c('#000084', '#004FFF', '#0EFFF0', '#8DFF70', '#F8FF05', colorRampPalette(c('#FF6900', '#8D0000'))(6)), limits = c(0, 2200000)) + 
geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group),fill = 'gray', col = 'gray') +
geom_point(data = env, aes(x = longitude, y = latitude, shape = watermass), size = 1) +
scale_shape_manual(values = c(4, 13, 1), limits = c('Plume', 'Mixed', 'SCS')) +
theme_bw() + 
theme(panel.grid.minor = element_blank(), legend.background = element_blank()) +
labs(x = 'Longitude (°E)', y = 'Latitude (°N)', fill = 'BA', shape = 'Watermass') +
coord_cartesian(xlim = c(113.5, 118.5), ylim = c(20.5, 23.5))


############## Figure 1

# Sea surface salinity sourced from Copernicus Marine Data (https://data.marine.copernicus.eu/products)
Salinity <- read.delim('copernicus_salinity.txt')

# Average sea surface salinity from 21 July to 26 July 2022
Salinity <- subset(Salinity, date %in% c('2022/7/21', '2022/7/22', '2022/7/23', '2022/7/24', '2022/7/25', '2022/7/26'))
Salinity <- summaryBy(Salinity~x+y, data = Salinity, FUN = mean)
names(Salinity) <- c('lon', 'lat', 'var')

# Spatial interpolation
Salinity <- mba.surf(na.omit(Salinity[c('lon', 'lat', 'var')]), no.X = 300, no.Y = 300, extend = FALSE)
dimnames(Salinity$xyz.est$z) <- list(Salinity$xyz.est$x, Salinity$xyz.est$y)
Salinity <- melt(Salinity$xyz.est$z, varnames = c('lon', 'lat'), value.name = 'var')
Salinity <- na.omit(Salinity)
names(Salinity) <- c('x', 'y', 'Salinity.mean')

# Right panels: Plot map show sampling sites and water mass distribution
env <- read.csv('env_table.csv')
env$layer <- factor(env$layer, levels = c('surface', 'middle', 'bottom'))

ggplot() +
geom_point(data = env, aes(x = longitude, y = latitude, color = watermass), size = 1.5) +
scale_color_manual(values = c('#607ebc', '#94c128', '#e74946'), limits = c('Plume', 'Mixed', 'SCS')) +
facet_wrap(~layer, ncol = 1, scale = 'free') +
geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group),fill = 'gray', col = 'gray') +
theme_bw() + 
theme(panel.grid.minor = element_blank(), legend.background = element_blank()) +
labs(x = 'Longitude (°E)', y = 'Latitude (°N)', fill = '') +
coord_cartesian(xlim = c(113.5, 118.5), ylim = c(20.5, 23.5))

# Left panel: Plot map show average sea surface salinity
Salinity[which(Salinity$Salinity.mean<22),'Salinity.mean'] <- 21.5

ggplot() +
geom_raster(data = Salinity, aes(x = x, y = y, fill = Salinity.mean), interpolate = TRUE) + 
scale_fill_gradientn(colors = c(colorRampPalette(c('#000084', '#004FFF'))(15), '#004FFF', '#0EFFF0', '#8DFF70', '#F8FF05', '#FF6900', '#8D0000'), 
    breaks = c(22, 24, 26, 28, 30, 32, 34), labels = c('<22', '24', '26', '28', '30', '32', '34'), limits = c(21.5, 35)) + 
geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), fill = 'gray') +
geom_point(data = env, aes(x = longitude, y = latitude), size = 1.5) +
theme_bw() + 
theme(panel.grid.minor = element_blank(), legend.background = element_blank()) +
guides(fill = guide_colourbar(barwidth = 0.5, barheight = 15)) +
scale_y_continuous(breaks = c(18, 20, 22, 24, 26, 28), labels = c('18 °N', '20 °N', '22 °N', '24 °N', '26 °N', '28 °N')) +
scale_x_continuous(breaks = c(110, 112, 114, 116, 118, 120, 122), labels = c('110 °E', '112 °E', '114 °E', '116 °E', '118 °E', '120 °E', '122 °E')) +
labs(x = 'Longitude (°E)', y = 'Latitude (°N)', fill = 'Salinity') +
coord_cartesian(xlim = c(110, 122), ylim = c(18, 28))


