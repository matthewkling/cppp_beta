

library(raster)
library(colormap)
library(data.table)
library(vegan)
library(combinat)
library(patchwork)

library(tidyverse)

select <- dplyr::select



source("code/functions.R")


# phylo turnover matrices
files <- list.files("input_data/turnover_matrices", full.names=T)



### step 1: define the target color scheme

file <- files[6]
m <- loadm(file)
cells <- loadm(file, "cells")

o <- metaMDS(m, k=3, parallel=7)
#o <- metaMDS(m, k=3, trymax=2)

d <- cbind(cells, as.data.frame(o$points)) %>%
   separate(cells, c("x", "y"), sep=":") %>%
   mutate(x=as.integer(x), y=as.integer(y)) %>%
   mutate_each(list(rank), MDS1:MDS3) %>%
   as_tibble()

schemes <- expand.grid(order=1:6, inversion=1:8) %>%
   split(1:48) %>%
   lapply(function(x) colors3d(select(d, MDS1:MDS3), 
                               order=x$order, inversion=x$inversion)) %>%
   lapply(col2rgb) %>%
   lapply(t)

target <- schemes[[7]]


distance0 <- difff(target, target = target) # check that this is 0
distance <- sapply(schemes, difff, target = target)
best <- which.min(distance)

d$target <- rgb(schemes[[best]], maxColorValue=255)
ddd <- select(d, x, y, target)

saveRDS(ddd, "data/target_colors.rds")


### step 2: ordinate and color all 10 metrics

ord_rgb <- function(file, target_scheme){
   
   tag <- sub("California_15km_", "", basename(file))
   tag <- sub("\\.csv", "", tag)
   message(tag)
   
   m <- loadm(file)
   cells <- loadm(file, "cells")
   
   o <- cmdscale(m, k=3)
   od <- cbind(cells, as.data.frame(o)) %>%
      separate(cells, c("x", "y"), sep=":") %>%
      mutate(x=as.integer(x), y=as.integer(y)) %>%
      mutate_each(list(rank), V1:V3) %>%
      left_join(target_scheme) %>%
      na.omit() %>%
      as_tibble()
   
   
   # find the best of the 48 possible RGB mappings
   target <- t(col2rgb(od$target))
   schemes <- expand.grid(order=1:6, inversion=1:8) %>%
      split(1:48) %>%
      lapply(function(x) colors3d(select(od, V1:V3),
                                  order=x$order, inversion=x$inversion)) %>%
      lapply(col2rgb) %>%
      lapply(t)
   distance <- sapply(schemes, difff, target = target, fun = mean)
   best <- which.min(distance)
   od$color <- rgb(schemes[[best]], maxColorValue=255)
   
   od$metric <- tag
   od %>% select(metric, x, y, color)
}


ddd <- readRDS("data/target_colors.rds")

d <- map(files, ord_rgb, target_scheme = ddd)



p <- d[c(1, 4, 6, 8, 10,
         3, 5, 7, 9, 2)] %>%
   map(function(x){
      title <- x$metric[1] %>% 
         str_remove("_California_") %>%
         str_remove("Sorenson") %>%
         str_remove_all("[:digit:]") %>%
         str_replace("_", " ")
      if(str_detect(title, "clades") &
         str_detect(title, "phylo")){
         title <- str_replace(title, "clades",
                              "phylogram")
      }
      title <- str_remove(title, "_phylo")
      if(str_detect(title, "RW")) title <- paste0("\n", title)
      
      ggplot(x, aes(x, y)) + 
         geom_raster(fill = x$color) +
         theme_void() +
         labs(title = title) +
         scale_x_continuous(expand = c(0,0)) +
         scale_y_continuous(expand = c(0,0))
   }) %>%
   Reduce("+", .) +
   plot_layout(nrow = 2)

ggsave("figures/ordination_rgb_maps.png", 
       p, width = 12, height = 7, units = "in")
