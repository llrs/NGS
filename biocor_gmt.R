library("BioCor")
library("GSEABase")
paths2Genes <- geneIds(getGmt("../data/c2.all.v5.2.symbols.gmt"))

genes <- unlist(paths2Genes, use.names = FALSE)
pathways <- rep(names(paths2Genes), lengths(paths2Genes))
genes2paths <- split(pathways, genes)

gmt_sim <- mgeneSim(names(genes2paths), genes2paths)
save(gmt_sim, file = "C2v5.2gmt.RData")
