

library(dplyr)
library(tidyr)
library(geosphere)
library(raster)
library(ggplot2)
select <- dplyr::select

source("code/functions.r")
#setwd("E:/alphabet")

files <- list.files("input_data/turnover_matrices", full.names=T)

# load biotic distance matrices
m <- files %>% lapply(read.csv, stringsAsFactors=F)
cells <- m[[1]][,1]
m <- lapply(m, function(x)x[,2:ncol(x)])
names(m) <- sub("\\.csv", "", basename(files))


# unroll and combine
flatten <- function(x) x %>% 
  as.data.frame() %>%
  mutate(px1=colnames(x)) %>% 
  gather(px2, value, -px1)
f <- lapply(m, flatten)
for(i in names(f)) names(f[[i]]) <- c("px1", "px2", i)
f <- Reduce("full_join", f)
names(f) <- sub("California_", "", names(f))


# add geographic distances
g <- read.csv(files[1], stringsAsFactors=F) %>%
  dplyr::select(X) %>%
  separate(X, c("x", "y"), sep=":") %>%
  mutate_each(funs(as.integer))
coordinates(g) <- c("x", "y")
crs(g) <- crs(readRDS("E:/california_phylomodelling/climate/BCM2014_cwd1951-1980_wy_ave_HST.Rdata"))
g <- spTransform(g, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
gd <- distm(g)/1000
colnames(gd) <- colnames(m[[1]])
gd <- flatten(gd)
names(gd)[3] <- "geographic"
f <- full_join(gd, f)


###### add climatic distances

# load climate data
cr <- list.files("E:/california_phylomodelling/climate", full.names=T) %>%
  lapply(readRDS) %>%
  do.call("stack", .)
cr[[4]] <- log10(cr[[4]])

# load biodiverse raster to use as template
bd <- read.csv(files[1], stringsAsFactors=F) %>%
  select(X) %>%
  separate(X, c("Axis_0", "Axis_1"), sep=":") %>%
  mutate_each(funs(as.integer), Axis_0:Axis_1) %>%
  mutate(x=Axis_0, y=Axis_1, z=1) %>%
  select(x, y, z) %>%
  as.matrix() %>%
  rasterFromXYZ()

# aggregate climate data to biodiverse grid
cm <- as.data.frame(rasterToPoints(cr))
coordinates(cm) <- c("x", "y")
names(cm) <- c("cwd", "djf", "jja", "ppt")
cm <- rasterize(cm, bd, fun=mean)
cm <- mask(cm, bd)
cm <- as.data.frame(rasterToPoints(cm))
cm$cell <- paste(cm$x, cm$y, sep=":")

## climatic distances

cd <- cm %>%
  select(cwd:ppt) %>%
  scale() %>%
  dist() %>%
  as.matrix()
rownames(cd) <- colnames(cd) <- cm$cell
cd <- flatten(cd)

cdr <- cm %>%
  select(cwd:ppt) %>%
  lapply(dist) %>%
  lapply(as.matrix) %>%
  lapply(function(x){
    rownames(x) <- colnames(x) <- cm$cell
    return(x)
  }) %>%
  lapply(flatten)
for(i in 1:length(cdr)) names(cdr[[i]])[3] <- names(cdr)[i]
cdr <- Reduce("full_join", cdr)
cd <- full_join(cd, cdr)
names(cd)[names(cd)=="value"] <- "climate"



##### join climate data to mainframe
clean <- function(x){
  x <- sub("X\\.", "neg", x)
  x <- sub("X", "", x)
  x <- sub("\\.\\.", ":neg", x)
  x <- gsub("\\.", ":", x)
  x <- gsub("neg", "-", x)
  return(x)
}

d <- mutate_each(f, funs(clean), px1, px2) %>%
  #select(-x) %>%
  full_join(cd) %>%
  na.omit()



#####################


names(d)[grepl("Sorenson|phylo|RW", names(d))] <- paste0("p", names(d)[grepl("Sorenson|phylo|RW", names(d))])


# correlations

cm <- cor(as.matrix(select(d, p1_species_Sorenson:p9_cladogram_Sorenson_phylo)),
          as.matrix(select(d, climate:ppt)),
          method="spearman")

cf <- as.data.frame(cm) %>%
  mutate(bio = row.names(cm),
         bio = str_sub(bio, 2),
         bio = fix_names(bio)) %>%
  gather(clim, r, -bio)
cfbio <- cf %>%
  group_by(bio) %>%
  summarize(r=mean(r)) %>%
  arrange(desc(r))
cfclim <- cf %>%
  group_by(clim) %>%
  summarize(r=mean(r)) %>%
  arrange(desc(r))
cf <- mutate(cf,
             clim=factor(clim, levels=cfclim$clim),
             bio=factor(bio, levels=cfbio$bio))

# plot
p <- ggplot(cf, aes(bio, clim, color=r, 
                    label=substr(as.character(round(r, 2)), 2, 4))) +
  geom_point(size=18) +
  geom_text(size=6, color="white") +
  # scale_color_gradientn(colours=c("black", "red4", "orangered")) +
  scale_color_gradientn(colours=c("gray", "blue")) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.text=element_text(size=16),
        axis.title=element_blank(),
        plot.title=element_text(size=20, color="red4"),
        legend.position="none") +
  coord_fixed()# +
  #labs(title="Spearmans's correlation between\nclimatic & biotic turnover metrics")
ggsave("figures/climate/corrplot.png", p, width=8, height=6, units="in")



####

# The partial correlation is most often used when some third variable z is a plausible explanation of the correlation between X and Y.
# The semipartial is most often used when we want to show that some variable adds incremental variance in Y above and beyond other X variable

# partial correlations
library(ppcor)

stats <- names(d)[4:13]
pcm <- lapply(stats, function(x){
  pcor(cbind(d[,x], as.matrix(select(d, cwd:ppt))), method="spearman")$estimate[1,]
})
cm <- lapply(pcm, as.vector) %>%
  do.call("rbind", .) %>%
  as.data.frame() %>%
  select(-V1)
names(cm) <- c("cwd", "djf", "jja", "ppt")
cf <- cm %>%
  mutate(bio=stats,
         bio = str_sub(bio, 2),
         bio = fix_names(bio)) %>%
  gather(clim, r, -bio)

cfbio <- cf %>%
  group_by(bio) %>%
  summarize(r=mean(r)) %>%
  arrange(desc(r))
cfclim <- cf %>%
  group_by(clim) %>%
  summarize(r=mean(r)) %>%
  arrange(desc(r))
cf <- mutate(cf,
             clim=factor(clim, levels=cfclim$clim),
             bio=factor(bio, levels=cfbio$bio))
# plot
p <- ggplot(cf, aes(bio, clim, color=r, 
                    label=substr(as.character(round(r, 2)), 2, 4))) +
  geom_point(size=18) +
  geom_text(size=6, color="white") +
  # scale_color_gradientn(colours=c("black", "red4", "orangered")) +
  scale_color_gradientn(colours=c("gray", "blue")) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.text=element_text(size=16),
        axis.title=element_blank(),
        plot.title=element_text(size=20, color="red4"),
        legend.position="none") +
  coord_fixed()# +
#labs(title="Spearmans's correlation between\nclimatic & biotic turnover metrics")
ggsave("figures/climate/corrplot_partial.png", p, width=8, height=6, units="in")



# semi-partial correlations

pcm <- lapply(stats, function(x){
  # spcor return asymmetrical matrix; we want the COLUMN [,1] corresponding to community turnover
  spcor(cbind(d[,x], as.matrix(select(d, cwd:ppt))), method="spearman")$estimate[,1]
})
cm <- lapply(pcm, as.vector) %>%
  do.call("rbind", .) %>%
  as.data.frame() %>%
  select(-V1)
names(cm) <- c("cwd", "djf", "jja", "ppt")
cf <- cm %>%
  mutate(bio=stats,
         bio = str_sub(bio, 2),
         bio = fix_names(bio)) %>%
  gather(clim, r, -bio)

cfbio <- cf %>%
  group_by(bio) %>%
  summarize(r=mean(r)) %>%
  arrange(desc(r))
cfclim <- cf %>%
  group_by(clim) %>%
  summarize(r=mean(r)) %>%
  arrange(desc(r))
cf <- mutate(cf,
             clim=factor(clim, levels=cfclim$clim),
             bio=factor(bio, levels=cfbio$bio))

# plot
p <- ggplot(cf, aes(bio, clim, color=r, 
                    label=substr(as.character(round(r, 2)), 2, 4))) +
  geom_point(size=18) +
  geom_text(size=6, color="white") +
  # scale_color_gradientn(colours=c("black", "red4", "orangered")) +
  scale_color_gradientn(colours=c("gray", "blue")) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.text=element_text(size=16),
        axis.title=element_blank(),
        plot.title=element_text(size=20, color="red4"),
        legend.position="none") +
  coord_fixed()# +
#labs(title="Spearmans's correlation between\nclimatic & biotic turnover metrics")
ggsave("figures/climate/corrplot_semipartial.png", p, width=8, height=6, units="in")





# # todo: why is this showing CWD when prior analysis showed PPT?? apparently because of log ppt
# # todo: why are there negative correlations in the distance plot??
# 
# d$geo <- round(log10(d$geographic), 1)
# s <- split(d, d$geo)
# s <- lapply(s, function(x){
#   cm <- cor(as.matrix(select(x, p1_species_Sorenson:p9_cladogram_Sorenson_phylo)),
#             as.matrix(select(x, climate:ppt)),
#             method="spearman")
#   cf <- as.data.frame(cm) %>%
#     mutate(bio=row.names(cm)) %>%
#     gather(clim, r, -bio) %>%
#     mutate(geo=x$geo[1])
#   return(cf)
# })
# s <- do.call("rbind", s)
# 
# s$rw <- grepl("RW", s$bio)
# 
# p <- ggplot(mutate(s, bio=sub("_RW", "\nRW", bio),
#                    bio=sub("_Sorenson", "\nSorenson", bio)),
#             aes(geo, r, color=clim)) +
#   geom_hline(yintercept=0, color="gray") +
#   geom_line() +
#   facet_grid(.~bio) +
#   theme(legend.position="top") +
#   labs(x="log geographic distance",
#        y="spearman's r") +
#   scale_color_manual(values=c("blue", "green", "orange", "red", "black")) +
#   theme_minimal()
# ggsave("e:/alphabet/climate_betadiv/corr_distance_curves_2.png", p, width=12, height=4, units="in")
# 
# s$tree <- sub("_RW|_Sorenson", "", s$bio)
# s$tree <- sub("._", "", s$tree)
# s$tree <- sub("p", "", s$tree)
# s$tree <- sub("1", "", s$tree)
# 
# p <- ggplot(s, aes(geo, r, color=tree, linetype=rw)) +
#   geom_hline(yintercept=0, color="gray") +
#   geom_line() +
#   facet_grid(.~clim) +
#   theme(legend.position="right") +
#   labs(x="log geographic distance",
#        y="spearman's r") +
#   scale_color_manual(values=c("blue", "green", "orange", "red", "black")) +
#   theme_minimal()
# ggsave("e:/alphabet/climate_betadiv/corr_distance_curves.png", p, width=12, height=8, units="in")

