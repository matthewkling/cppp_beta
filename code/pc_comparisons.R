

library(tidyverse)
library(data.table)
library(corrplot)
library(gplots)
library(geosphere)
library(caret)
library(raster)
library(patchwork)
library(gridExtra)
library(grid)
select <- dplyr::select


source("code/functions.R")


# load distance matrices
m <- list.files("input_data/turnover_matrices", full.names=T) %>%
        map(read_csv)
cells <- m[[1]][,1]
m <- map(m, function(x)x[,2:ncol(x)])
names(m) <- sub("\\.csv", "", list.files("input_data/turnover_matrices"))


# unroll and combine
flatten <- function(x) x %>% 
        as.data.frame() %>%
        mutate(px1=colnames(x)) %>% 
        gather(px2, value, -px1)
f <- lapply(m, flatten)
for(i in names(f)) names(f[[i]]) <- c("px1", "px2", i)
f <- Reduce("full_join", f)
names(f) <- sub("California_", "", names(f))







#### correlations and dendrogram

# pearson <- cor(select(f, -px1, -px2), use="pairwise.complete.obs", method="pearson")
spearman <- cor(select(f, -px1, -px2), use="pairwise.complete.obs", method="spearman")
colnames(spearman) <- rownames(spearman) <- fix_names(colnames(spearman))

pdf("figures/correlations_spearman.pdf", width=6, height=6)
heatmap.2(spearman, 
          dendrogram="both", symm=T, revC=T,
          key.title=NA, key.xlab="spearman correlation", key.ylab=NA,
          margins=c(10,10), 
          colsep=1:10, rowsep=1:10,
          col=colorRampPalette(c("gray", "blue")),
          trace="none", tracecol="white", hline=NA, vline=NA)
dev.off()




#### PCA of the 10 metrics

# univariate deskew
z <- f[,!grepl("geographic", names(f))]
z <- na.omit(z)
coords <- z[,c("px1", "px2")]
z <- z[,3:ncol(z)]
trans <- preProcess(sample_n(z, 100000), method=c("center", "scale", "YeoJohnson"))
z <- predict(trans, z)

# pca 
z <- prcomp(na.omit(z))

# pc loadings plot

imp <- summary(z)$importance
imp <- round(imp[2, 1:3]*100, 1)

r <- z$rotation %>% 
        as.data.frame() %>%
        mutate(var=row.names(z$rotation)) %>%
        select(var, PC1:PC3) %>%
        gather(pc, loading, PC1:PC3) %>%
        mutate(rw=grepl("RW", var))

cats <- c("species", "clade", "phylo", "chrono", "clado")
for(cat in cats) r$cat[grepl(cat, r$var)] <- cat
r$cat <- factor(r$cat, levels=cats)

r$pc <- factor(r$pc, labels=paste0(unique(r$pc), "\n(", imp, "% var)"))

p <- ggplot(r, aes(pc, loading, 
                   color=cat, group=paste(cat, rw),
                   shape=rw, linetype=rw)) +
        geom_hline(yintercept = 0, color = "gray") +
        geom_line() +
        geom_point() +
        theme_minimal() +
        scale_color_manual(values=c("black", "blue", "goldenrod1", "red", "limegreen")) +
        theme(axis.title.x=element_blank()) +
        labs(title="PCA across distance matrices for 10 turnover metrics",
             color="distance\nmetric",
             linetype="range\nweighted?",
             shape="range\nweighted?",
             y="PC loading")

ggsave("figures/pc_loadings.png", p, width=6, height=4, units="in")







# convert PCs back to distance matrices and ordinate into 3d space
pc2map <- function(i, invert=F){
        x <- data.frame(coords) %>%
                mutate(dist=z$x[,i]) %>% #,
                #dist=ecdf(dist)(dist)) %>%
                spread(px2, dist)
        row.names(x) <- x$px1
        x <- select(x, -px1) %>%
                as.matrix() %>%
                scales::rescale() # the PCA generated negative distances--remove those
        
        if(invert) x <- 1-x
        
        x[is.na(x)] <- 0
        
        
        dst <- hclust(as.dist(x))
        cluster <- cutree(dst, 10)
        
        x <- cmdscale(x, k=3)
        s <- as.data.frame(x) %>% as_tibble() %>%
                mutate(px=row.names(x),
                       pixel = sub("X\\.", "-", px),
                       pixel = sub("X", "", pixel),
                       pixel = sub("\\.\\.", ".-", pixel),
                       pixel = sub("\\.", "_", pixel)) %>%
                separate(pixel, c("x", "y"), sep=":") %>%
                mutate(x=as.integer(x), y=as.integer(y),
                       cluster=cluster)
        colors <- colormap::colors3d(s[,1:3], trans="ecdf")
        col <- t(col2rgb(colors))
        var <- apply(s[,1:3], 2, sd)
        var <- var/max(var)
        for(j in 1:ncol(col)) col[,j] <- col[,j] * var[j] # colorspace shrinkage to match variance
        colors <- rgb(col, maxColorValue=255)
        
        s <- mutate(s, color=colors) %>%
                cbind(col) %>%
                group_by(cluster) %>%
                mutate_each(funs(mean), red:blue) %>%
                mutate(cluster_color=rgb(red, green, blue, maxColorValue=255))
        
        neg <- "+"
        if(invert) neg <- "-"
        
        label <- paste0(neg, "PC", i, ":\n",
                        case_when(i == 1 & !invert ~ "overall",
                                 i == 1 & invert ~ "random",
                                 i == 2 & !invert ~ "shallow",
                                 i == 2 & invert ~ "deep",
                                 i == 3 & !invert ~ "widespread",
                                 i == 3 & invert ~ "endemic"), 
                       "\nstructure")
        
        pcont <- ggplot(s, aes(x, y)) + 
                geom_raster(fill=s$color) +
                ggmap::theme_nothing() +
                annotate(geom="text", x=100000, y=300000, label=label, 
                         size=20, hjust = 0, lineheight = .8)
        
        pclus <- ggplot(s, aes(x, y)) + 
                geom_raster(fill=s$cluster_color) +
                ggmap::theme_nothing() +
                annotate(geom="text", x=100000, y=300000, label=label, 
                         size=20, hjust = 0, lineheight = .8)
        
        return(list(pcont, pclus))
}

maps <- lapply(1:3, pc2map)
mapsi <- lapply(1:3, pc2map, invert=T)

p <- arrangeGrob(maps[[1]][[1]], maps[[2]][[1]], maps[[3]][[1]],
                 mapsi[[1]][[1]], mapsi[[2]][[1]], mapsi[[3]][[1]],
                 nrow=2)
png("figures/pc_maps.png", width=3000, height=2400)
grid.draw(p)
dev.off()


p <- arrangeGrob(maps[[1]][[2]], maps[[2]][[2]], maps[[3]][[2]],
                 mapsi[[1]][[2]], mapsi[[2]][[2]], mapsi[[3]][[2]],
                 nrow=2)
png("figures/pc_maps_clusters.png", width=3000, height=2400)
grid.draw(p)
dev.off()




