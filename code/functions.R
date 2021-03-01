fix_names <- function(x){
        x <- str_remove(x, "_California")
        recode(x, 
               "1_species_Sorenson" = "species",
               "10_cladogram_RW_phylo" = "cladogram RW",     
               "2_species_RW" = "species RW",
               "3_clades_Sorenson" = "clades",  
               "4_clades_RW" = "clades RW",
               "5_clades_Sorenson_phylo" = "phylogram",    
               "6_clades_RW_phylo" = "phylogram RW",
               "7_chronogram_Sorenson_phylo" = "chronogram",
               "8_chronogram_RW_phylo" = "chronogram RW",
               "9_cladogram_Sorenson_phylo" = "cladogram",)
}



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

difff <- function(x, target, fun = mean){
        x <- cbind(x, target)
        x <- apply(x, 1, dst)
        fun(x)
}