######### Required Libraries
library(data.table)
library(reshape2)
library(trackViewer)
library(GenomicRanges)

######### Get VariantData
variants <- fread("Desktop/IGIB/TAM/variants_tetC.txt",
                  header = TRUE, stringsAsFactors = FALSE, sep = "\t")
variants <- melt.data.table(variants,
                            id.vars = c("GENE","VARIANT_TYPE","POS","PROTEIN_CHANGE"))
variants <- variants[value != 0]

######## GetDomainData
features <- fread("Desktop/IGIB/TAM/features.txt",
                  header = TRUE, stringsAsFactors = FALSE, sep = "\t")
features_gr <- GRanges("tet(C)",IRanges(start = features[,Start],
                                        width = features[,Width],
                                        names = features[,Domain],
                                        height = 0.05))

######### Prep for the Lollipop Plot
variants <- variants[,Score1 := max(value),by = PROTEIN_CHANGE]
variants <- unique(variants[,.(GENE,POS,VARIANT_TYPE,PROTEIN_CHANGE,Score1)])
variants_gr <- GRanges("tet(C)",IRanges(variants[,POS],
                                        width = 1,
                                        names = variants[,PROTEIN_CHANGE],
                                        color = "Green",
                                        score = variants[,Score1]))
lolliplot(variants_gr,
          features_gr,
          #ranges = tile(variants_gr,n = 1),
          lollipop_style_switch_limit = 25,
          cex = 1.2,
          legend = FALSE)
