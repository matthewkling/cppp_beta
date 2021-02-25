

library(raster)
library(colormap)
library(data.table)
library(vegan)
library(combinat)
library(patchwork)

library(tidyverse)

select <- dplyr::select


# utility functions

loadm <- function(file, get="dist"){
   r <- fread(file) %>%
      as.data.frame()
   if(get=="cells") return(r[,1])
   r <- r[,2:ncol(r)]
   r[is.na(r)] <- 0
   as.matrix(r)
}

rank <- function(x)ecdf(x)(x)

dst <- function(y){
   y <- matrix(y, nrow=2, byrow=T)/255
   y <- sqrt(sum((y[1,] - y[2,])^2))
   y
}

diff <- function(x, target, fun = mean){
   x <- cbind(x, target)
   x <- apply(x, 1, dst)
   fun(x)
}





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


distance0 <- diff(target, target = target) # check that this is 0
distance <- sapply(schemes, diff, target = target)
best <- which.min(distance)

d$target <- rgb(schemes[[best]], maxColorValue=255)
ddd <- select(d, x, y, target)




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
   distance <- sapply(schemes, diff, target = target, fun = mean)
   best <- which.min(distance)
   od$color <- rgb(schemes[[best]], maxColorValue=255)
   
   od$metric <- tag
   od %>% select(metric, x, y, color)
}

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
