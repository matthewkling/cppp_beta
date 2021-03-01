
# spatial cluster analysis of cali flora/climate, on 15k biodiverse grid

library(tidyverse)
library(raster)
library(spdep)
library(FNN)
library(colormap)
library(gridExtra)
library(grid)
library(MASS)
library(data.table)
library(patchwork)
select <- dplyr::select

source("code/functions.R")



# find the best of the 48 possible RGB mappings
rgb_match <- function(bd, target_scheme){
   od <- bd %>%
      left_join(target_scheme) %>%
      na.omit() %>%
      as_tibble()
   target <- t(col2rgb(od$target))
   schemes <- expand.grid(order=1:6, inversion=1:8) %>%
      split(1:48) %>%
      lapply(function(x) colors3d(select(od, V1:V3), trans = "ecdf",
                                  order=x$order, inversion=x$inversion)) %>%
      lapply(col2rgb) %>%
      lapply(t)
   distance <- sapply(schemes, difff, target = target, fun = mean)
   best <- which.min(distance)
   od$color <- rgb(schemes[[best]], maxColorValue=255)
   od$metric <- tag
   od %>% select(metric, x, y, color)
}


poly_data <- function(kd, border_only=F){
   kd <- kd %>%
      select(x, y, cluster)
   names(kd)[3] <- "clust"
   if(border_only) kd$clust <- 1
   poly <- rasterFromXYZ(select(kd, x, y, clust))
   poly <- rasterToPolygons(poly, dissolve=T)
   poly <- broom::tidy(poly)
   return(poly)
}



files <- list.files("input_data/turnover_matrices", full.names=T)
files <- c("climate", files)

for(file in files){
   
   tag <- sub("California_15km_", "", basename(file))
   tag <- sub("\\.csv", "", tag)
   
   # climate
   if(file == "climate"){
      r <- readRDS("input_data/climate_data/derived/climate_pc.rds") %>%
         rasterToPoints() %>%
         as.data.frame() %>%
         mutate(ELEMENT = paste(x, y, sep = ":")) %>%
         select(V1 = pc1, V2 = pc2, V3 = pc3, ELEMENT) %>%
         mutate(V1 = rank(V1),
                V2 = rank(V2),
                V3 = rank(V3))
   } else{
      r <- read.csv(file)
      cells <- r[,1]
      r <- as.matrix(r[,2:ncol(r)])
      r[is.na(r)] <- 0
      # r <- sammon(r, k=3, niter=10000)
      # r <- data.frame(r$points)
      r <- cmdscale(r, k=3) %>% as.data.frame()
      r$ELEMENT <- cells
   }
   
   
   # pixels data
   bd <- read.csv("e:/California_phylomodelling/biodiverse/California_Clades_clean_All_final_epsg_3310_trimmed_analysed_output_SPATIAL_RESULTS.csv",
                  stringsAsFactors=F) %>%
      select(ELEMENT, Axis_0, Axis_1)
   bd <- left_join(bd, r) %>%
      na.omit()
   names(bd) <- c("pixel", "x", "y", "z1", "z2", "z3")
   
   # stash data for joint climate/phylo run
   # if(file==files[1]) phy <- bd
   
   # format data for skater
   coords <- select(bd, x, y) %>%
      as.matrix()
   phylo <- select(bd, z1:z3) %>%
      as.matrix()
   
   # construct spatial connectivity graph 
   gabn <- gabrielneigh(coords, nnmult=4)
   nb <- graph2nb(gabn, sym=T)
   
   # weight graph edges by climate dissimiliarity
   costs <- nbcosts(nb, phylo)
   nbw <- nb2listw(nb, costs)
   
   
   library(furrr)
   plan(multisession, workers = 7)
   
   # grow minimum spanning tree within connectivity graph
   # mst <- mstree(nbw, ini=1)
   niter <- 100
   msts <- 1:niter %>% map(function(x) mstree(nbw, ini=x))
   
   target <- readRDS("data/target_colors.rds")
   d3 <- bd %>% select(x, y, V1=z1, V2=z2, V3=z3) %>% rgb_match(target)
   d3 <- bind_cols(d3, col2rgb(d3$color) %>% t() %>% as.data.frame())
   # ggplot(d3, aes(x, y)) + geom_tile(fill = d3$color) + theme_void()
   
   # partition tree into clusters
   for(nclust in seq(10, 35, 5)[3]){
      
      kd <- cbind(coords, phylo) %>%
         cbind(prcomp(scale(phylo))$x) %>%
         as.data.frame() %>%
         as_tibble()
      
      sks <- msts %>%
         future_map(function(x){
            skater(x[,1:2], phylo, method="euclidean", ncuts=nclust-1)})
      saveRDS(sks, paste0("data/skater/", tag, "_", nclust, ".rds"))
      
      patterns <- sks %>% map(function(x) x$groups) %>% 
         do.call("rbind", .) %>% as.data.frame()
      unique_patterns <- distinct(patterns)
      freqs <- apply(unique_patterns, 1, 
                     function(x) apply(patterns, 1, 
                                       function(y) all.equal(x, y))) %>%
         apply(2, function(x) sum(x=="TRUE"))
      modal <- unique_patterns[which.max(freqs)[1],] %>%
         as.matrix() %>% as.vector()
      
      polys <- 1:length(freqs) %>%
         map_df(function(i){
            z <- unique_patterns[i,] %>% as.matrix() %>% as.vector()
            kd %>% mutate(cluster = z) %>% poly_data() %>%
               mutate(pattern = i,
                      freq = freqs[i])
         }) %>% 
         mutate(modal = pattern==pattern[freq==max(freq)][1])
      
      kd$cluster <- modal
      
      border <- poly_data(kd, border_only = T)
      
      kd <- kd %>% 
         left_join(d3)
      colors <- kd %>%
         arrange(cluster) %>%
         group_by(cluster) %>%
         summarize_each(funs(mean), red:blue) %>%
         select(-cluster) %>%
         na.omit() %>%
         rgb(maxColorValue=255)
      
      pg <- ggplot() +
         geom_tile(data=kd, aes(x, y, fill=factor(cluster))) +
         geom_path(data=polys %>% filter(!modal),
                   aes(long, lat, group=paste(group, pattern), alpha = freq),
                   color="black", size=.5, linetype = "dotted") +
         geom_path(data=polys %>% filter(modal),
                   aes(long, lat, group=paste(group, pattern)),
                   color="white", size=.5) +
         geom_path(data=border,
                   aes(long, lat, group=group),
                   color="white", size=.5) +
         scale_alpha_continuous() +
         scale_fill_manual(values=colors) +
         theme_void() +
         theme(legend.position="none") +
         scale_x_continuous(expand = c(0,0)) +
         scale_y_continuous(expand = c(0,0))
      
      # rgb plot
      prgb <- ggplot(kd, aes(x, y)) +
         geom_tile(fill = kd$color) +
         theme_void() +
         scale_x_continuous(expand = c(0,0)) +
         scale_y_continuous(expand = c(0,0))
      
      # ordination scatterplots
      pe12 <- ggplot(kd[sample(nrow(kd),nrow(kd)),], # permutate rows to prevent biased overplotting
                     aes(PC1, PC2, color=factor(cluster))) +
         geom_point(size=1) +
         scale_color_manual(values=colors) +
         theme_bw() +
         theme(legend.position="none",
               plot.background = element_blank(),
               title=element_text(size=10)) +
         coord_fixed()
      pe13 <- ggplot(kd[sample(nrow(kd),nrow(kd)),], # permutate rows to prevent biased overplotting
                     aes(PC1, PC3, color=factor(cluster))) +
         geom_point(size=1) +
         scale_color_manual(values=colors) +
         theme_bw() +
         theme(legend.position="none",
               plot.background = element_blank(),
               title=element_text(size=10)) +
         coord_fixed()
      
      # combined chart
      p <- prgb + inset_element(pe12, 0.425, 0.55, 1, 1) +
         pg + inset_element(pe13, 0.425, 0.55, 1, 1) +
         plot_annotation(title = fix_names(tag), theme = theme(plot.title = element_text(size = 20)))
      ggsave(paste0("figures/skater/fractures_", fix_names(tag), "_clusters_", nclust, ".png"), p,
             width=10, height=6, units = "in")
   }
}





# ################# combo climate + phylo #################
# 
# file <- files[1]
# tag <- sub("California_15km_", "", basename(file))
# tag <- sub("\\.csv", "", tag)
# 
# 
# phy <- as.data.frame(phy) # only using the first of 4 stats for now
# climate <- as.data.frame(climate)
# coords <- as.data.frame(coords)
# 
# d <- coords %>%
#    cbind(climate) %>%
#    left_join(phy) %>%
#    select(-pixel)
# 
# coords <- select(d, x, y) %>% as.matrix()
# d <- select(d, -x, -y) %>% as.matrix()
# d <- prcomp(scale(d))$x[,1:3]
# 
# # construct spatial connectivity graph 
# gabn <- gabrielneigh(coords, nnmult=4)
# nb <- graph2nb(gabn, sym=T)
# 
# # weight graph edges by dissimiliarity
# costs <- nbcosts(nb, d)
# nbw <- nb2listw(nb, costs)
# 
# # grow minimum spanning tree within connectivity graph
# mst <- mstree(nbw, ini=1)
# 
# # partition tree into clusters
# for(nclust in seq(10, 35, 5)){
#    sk <- skater(mst[,1:2], d, 
#                 method="euclidean", ncuts=nclust-1)
#    
#    # combine data
#    kd <- cbind(coords, d) %>%
#       #cbind(prcomp(scale(d))$x) %>%
#       as.data.frame()
#    kd$cluster <- sk$groups
#    
#    # geography of clusters
#    colors <- kd %>%
#       arrange(cluster) %>%
#       group_by(cluster) %>%
#       summarize_each(funs(mean), PC1:PC3) %>%
#       select(-cluster) %>%
#       as.matrix() %>%
#       colors3d(order=6, inversion=6)
#    
#    pg <- ggplot() +
#       geom_raster(data=kd, aes(x,y,fill=factor(cluster))) +
#       scale_fill_manual(values=colors) +
#       labs(title="GEOGRAPHY",
#            title=element_text(size=20)) +
#       theme_minimal() +
#       theme(legend.position="none",
#             axis.title=element_blank(), axis.text=element_blank()) +
#       coord_fixed()
#    
#    # clusters in phylo space
#    pe <- ggplot(kd[sample(nrow(kd),nrow(kd)),], # permutate rows to prevent biased overplotting
#                 aes(PC1, PC2, color=factor(cluster))) +
#       geom_point(size=3) +
#       scale_color_manual(values=colors) +
#       labs(title="CLIMATE & PHYLOGENY",
#            title=element_text(size=20)) +
#       theme_minimal() +
#       theme(legend.position="none") +
#       coord_fixed()
#    
#    # combined chart
#    p <- arrangeGrob(pg, pe, ncol=2, widths=c(1,1))
#    png(paste0("charts/combo/", tag, "_clusters_", nclust, ".png"), width=1000, height=500)
#    grid.draw(p)
#    dev.off()
# }
# 
# 
# # continuous plot in RGB space
# 
# p1 <- ggplot() +
#    geom_raster(data=kd, aes(x, y), 
#                fill=colors3d(kd[,c("PC1", "PC2", "PC3")], order=3, inversion=3)) +
#    theme_minimal() +
#    theme(legend.position="none", panel.grid=element_blank(),
#          axis.title=element_blank(), axis.text=element_blank()) +
#    coord_fixed()
# 
# p2 <- ggplot() +
#    geom_raster(data=kd, aes(x, y), 
#                fill=colors3d(kd[,c("PC1", "PC2", "PC3")], order=1, inversion=5)) +
#    theme_minimal() +
#    theme(legend.position="none", panel.grid=element_blank(),
#          axis.title=element_blank(), axis.text=element_blank()) +
#    coord_fixed()
# 
# p3 <- ggplot() +
#    geom_raster(data=kd, aes(x, y), 
#                fill=colors3d(kd[,c("PC1", "PC2", "PC3")], order=6, inversion=6)) +
#    theme_minimal() +
#    theme(legend.position="none", panel.grid=element_blank(),
#          axis.title=element_blank(), axis.text=element_blank()) +
#    coord_fixed()
# 
# p <- arrangeGrob(p1, p2, p3, nrow=1)
# 
# png(paste0("charts/combo/", tag, "_colors3d.png"), width=1000, height=600)
# grid.draw(p)
# dev.off()
# 
# 
# 
# 
# 
# # jepson ecoregions
# je <- readOGR("Geographic_Subdivisions_of_California_TJMII_v2_060415",
#               "Geographic_Subdivisions_of_California_TJMII_v2_060415")
# plot(je, col=distant_colors(length(unique(je$JEPCODE))))

