library("BioCor")
library("WGCNA")
library("org.Hs.eg.db")

load("~/Documents/geneSim.RData", verbose = TRUE)
load("../../Late_Network.RData", verbose = TRUE)
enableWGCNAThreads(4)
# ids <- mapIds(org.Hs.eg.db, keys = colnames(o), keytype = "ENTREZID",
#               column = "SYMBOL")
# colnames(o) <- ids
# rownames(o) <- ids
# Preparation to match the matrix
# m <- matrix(ncol = ncol(data.wgcna), nrow = ncol(data.wgcna),
#             dimnames = list(colnames(data.wgcna), colnames(data.wgcna)))
# bio_mat <- AintoB(o, m)

# sim_wgcna <- bicor(data.wgcna,
#                    use = "p",
#                    maxPOutliers = 0.10,
#                    nThreads = 6)

# sft <- pickSoftThreshold(data.wgcna,
#                          corFnc = bicor,
#                          corOptions = list(maxPOutliers = 0.10),
#                          networkType = "signed hybrid",
#                          blockSize = 5000, verbose = 1)
# adj <- adjacency(data.wgcna,
#                  type = "signed hybrid",
#                  power = sft$powerEstimate,
#                  corFnc = "bicor",
#                  corOptions = "nThreads = 6, use = 'p', maxPOutliers = 0.10")
# TOM <- TOMsimilarity(adj, TOMType = "signed")
# dissTOM <- 1 - TOM
# geneTree <- hclust(as.dist(dissTOM), method = "average")
#
# dynamicMods <- cutreeHybrid(dendro = geneTree, distM = dissTOM,
#                             deepSplit = 2, pamRespectsDendro = FALSE,
#                             minClusterSize = 30)
# y <- labels2colors(dynamicMods$labels)
# names(y) <- colnames(data.wgcna)
# modules <- y
# save(modules, file = "msimilarity_n.RData")
# First similarity then from there ####
# Function to simplify steps
pick_modules <- function(x){
  # It is the same as pickSoftThreshold.fromSimilarity
  sft <- pickSoftThreshold(x,
                           networkType = "signed hybrid",
                           corFnc = bicor,
                           corOptions = list(maxPOutliers = 0.10),
                           blockSize = 5000, verbose = 1,
                           RsquaredCut = 0.80)
  # It is the same as adjacency.fromSimilarity
  adjacency <- adjacency(x,
                         power = sft$powerEstimate,
                         type = "signed hybrid",
                         corFnc = "bicor",
                         corOptions = "nThreads = 6, use = 'p', maxPOutliers = 0.10")
  # Once we have the similarities we can calculate the TOM with TOM
  TOM <- TOMsimilarity(adjacency, TOMType = "signed")
  dissTOM <- 1 - TOM
  geneTree <- hclust(as.dist(dissTOM), method = "average")

  dynamicMods <- cutreeHybrid(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = 30)
  y <- labels2colors(dynamicMods$labels)
  names(y) <- colnames(data.wgcna)
  y
}
# sims <- list(sim = sim_wgcna, reactome = bio_mat)
# save(sims, file = "sims_late.RData")
load("sims_late.RData")
sim_mean <- similarities(sims, mean, na.rm = TRUE)
# load(file = "similarity_mean.RData")
modules <- pick_modules(sim_mean)
save(sim_mean, modules, file = "similarity_mean.RData")
save(modules, file = "msimilarity_mean.RData")

sim_09 <- similarities(sims, weighted.sum, w = c(0.9, 0.1))
# load(file = "similarity_09.RData")
modules <- pick_modules(sim_09)
save(sim_09, modules, file = "similarity_09.RData")
save(modules, file = "msimilarity_09.RData")

sim_08 <- similarities(sims, weighted.sum, w = c(0.8, 0.2))
# load(file = "similarity_08.RData")
modules <- pick_modules(sim_08)
save(sim_08, modules, file = "similarity_08.RData")
save(modules, file = "msimilarity_08.RData")

sim_07 <- similarities(sims, BioCor::weighted.sum, w = c(0.7, 0.3))
# load(file = "similarity_07.RData")
modules <- pick_modules(sim_07)
save(sim_07, modules, file = "similarity_07.RData")
save(modules, file = "msimilarity_07.RData")

sim_06 <- similarities(sims, BioCor::weighted.sum, w = c(0.6, 0.4))
# load(file = "similarity_06.RData")
modules <- pick_modules(sim_06)
save(sim_06, modules, file = "similarity_06.RData")
save(modules, file = "msimilarity_06.RData")

sim_05 <- similarities(sims, BioCor::weighted.sum, w = c(0.5, 0.5))
# load(file = "similarity_05.RData")
modules <- pick_modules(sim_05)
save(sim_05, modules, file = "similarity_05.RData")
save(modules, file = "msimilarity_05.RData")

# Adjacency and then add the similarity of BioCor
# sft <- pickSoftThreshold(data.wgcna)
# expr.adj <- adjacency(data.wgcna, power = sft$powerEstimate,
#                       type = "signed hybrid")
#
# pick2_modules <- function(x) {
#   TOM <- TOMsimilarity(x, TOMType = "signed")
#   dissTOM <- 1 - TOM
#   geneTree <- hclust(as.dist(dissTOM), method = "average")
#   # We can use a clustering tool to group the genes
#   dynamicMods <- cutreeHybrid(dendro = geneTree, distM = dissTOM,
#                               deepSplit = 2, pamRespectsDendro = FALSE,
#                               minClusterSize = 30)
#   y <- labels2colors(dynamicMods$labels)
#   names(y) <- colnames(data.wgcna)
#   y
# }
#
# sims <- list(exp = expr.adj, react = sims$reactome)
# sim_mean <- similarities(sims, mean, na.rm = TRUE)
# modules <- pick2_modules(sim_mean)
# save(sim_mean, modules, file = "similarity2_mean.RData")
# save(modules, file = "msimilarity2_mean.RData")
#
# # To combine the similarities either weighted.sum :
# sim_09 <- similarities(sims, weighted.sum, w = c(0.9, 0.1))
# modules <- pick2_modules(sim_09)
# save(sim_09, modules, file = "similarity2_09.RData")
# save(modules, file = "msimilarity2_09.RData")
#
# sim_08 <- similarities(sims, weighted.sum, w = c(0.8, 0.2))
# modules <- pick2_modules(sim_08)
# save(sim_08, modules, file = "similarity2_08.RData")
# save(modules, file = "msimilarity2_08.RData")
#
# sim_07 <- similarities(sims, weighted.sum, w = c(0.7, 0.3))
# modules <- pick2_modules(sim_07)
# save(sim_07, modules, file = "similarity2_07.RData")
# save(modules, file = "msimilarity2_07.RData")
#
# sim_06 <- similarities(sims, weighted.sum, w = c(0.6, 0.4))
# modules <- pick2_modules(sim_06)
# save(sim_06, modules, file = "similarity2_06.RData")
# save(modules, file = "msimilarity2_06.RData")
#
# sim_05 <- similarities(sims, weighted.sum, w = c(0.5, 0.5))
# modules <- pick2_modules(sim_05)
# save(sim_05, modules, file = "similarity2_05.RData")
# save(modules, file = "msimilarity2_05.RData")

paths <- list.files(pattern="msimilarity_")
out <- sapply(paths, function(x){
  y <- load(x)
  if ("modules" %in% y){
    table(modules)
  } else {
    table(moduleColors)
  }
  })

nam <- unique(unlist(sapply(out, names)))
makeDF <- function(List, Names) { # Set the element to the right column
  m <- t(vapply(List,
                FUN = function(X) unlist(X)[Names],
                FUN.VALUE = numeric(length(Names))))
  as.data.frame(m)
}

tables <- makeDF(out, nam)
colnames(tables) <- nam # In some cases it doesn't apply correctly the names
full <- cbind(type = rownames(tables), tables)
new.df <- melt(full)
new.df <- new.df[!is.na(new.df$value), ]

ggplot(new.df) +
  geom_violin(aes(x = type, y = value)) +
  theme_bw() +
  scale_y_log10() +
  ylab("# of genes by module") +
  xlab("Method to build the network") +
  geom_point(aes(x = type, y = value, color = variable, size = value)) +
  scale_color_manual(values = unique(as.character(new.df$variable))) +
  theme(legend.position = "none") +
  ggtitle("Distribution and size of the modules in each network")
