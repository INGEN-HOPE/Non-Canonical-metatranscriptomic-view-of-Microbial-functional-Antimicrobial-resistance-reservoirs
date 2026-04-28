#### Required Packages
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(reshape2)
library(ggplot2)
library(dplyr)

#### Function for CIGAR String Mapping
parse_cigar_matches <- function(cigar, start_pos) {
  ops <- unlist(strsplit(cigar, "(?<=[A-Z=])", perl = TRUE))
  ref_pos <- start_pos
  out <- list()
  
  for (op in ops) {
    len <- as.integer(sub("[A-Z=]", "", op))
    #print(len)
    type <- gsub("[0-9]+", "", op)
    #print(type)
    
    if (type == "=") {
      out[[length(out)+1]] <- data.table(
        match_start = ref_pos,
        match_end   = ref_pos + len,
        type = type
      )
      ref_pos <- ref_pos + len
    } else if (type %in% c("X","D")) {
      ref_pos <- ref_pos + len
    #   out[[length(out)+1]] <- data.table(
    #     #match_start = ref_pos,
    #     #match_end   = ref_pos + len,
    #     type = type
    #   )
    } else if (type %in% c("I","S","H")) {
    #     out[[length(out)+1]] <- data.table(
    #     #match_start = ref_pos,
    #     #match_end   = ref_pos + len,
    #     type = type
    #   )
      # insertion/clip: consume read, not reference
      # do nothing to ref_pos
    }
  }
  
  if (length(out) == 0) return(data.table(match_start=integer(), match_end=integer()))
  rbindlist(out,fill = TRUE)
}

######### Get the list of files
samples <- fread("/lustre/samhita.p/tam_amr_dengue80/sampleNames.txt",
                header = FALSE, stringsAsFactors = FALSE)
bams <- data.table("bamFiles" = list.files("/lustre/samhita.p/tam_amr_dengue80/bam_amr_tam_80",
                    "*.bam$",full.names = TRUE))[-c(1,2,3,4,5,6,7)]
bamIndex <- data.table("bamIndexFiles" = list.files("/lustre/samhita.p/tam_amr_dengue80/bam_amr_tam_80",
                "*.bam.bai",full.names = TRUE))[-c(1,2,3,4,5,6,7)]
kraken <- data.table("krakenFiles" = list.files("/lustre/samhita.p/tam_amr_dengue80/kraken_80_pf",
                    "*.txt$",full.names = TRUE))
blast <- fread("/lustre/samhita.p/tam_amr_dengue80/blastData2.csv",
               sep = "\t",header = TRUE, stringsAsFactors = FALSE)                    

abundantGeneList <- fread("/lustre/samhita.p/tam_amr_dengue80/abundantGeneList.tsv",
                            header = TRUE,stringsAsFactors = FALSE)  
###########
lapply(samples$V1[1:80], function(x){
  print(x)
    blast_s <- blast[SampleName == x][,c("ReadName","ScientificName")]
    blast_s <- blast_s[,c("qname1","extra") := tstrsplit(ReadName," ", fixed = TRUE)]

    bamData <- as.data.table(readGAlignments(bams[grep(x, bamFiles)][[1]],
                                index = bamIndex[grep(x, bamIndexFiles)][[1]],
                                ScanBamParam(what = c("qname")),use.names = FALSE))[,.(qname,seqnames,start,end,cigar)]


    bamData <- bamData[,c("qname1","extra") := tstrsplit(qname," ", fixed = TRUE)]  
    bamData <- bamData[as.character(seqnames) %in% abundantGeneList[,ReferenceSequence]]
    bamData <- bamData[,seqnames := as.character(seqnames)][order(qname1)]
    bamData <- bamData[,rowId := 1:.N, by = qname1][,qname2 := paste0(qname1,"_",rowId)]
    bamData_rowWise <- bamData[,parse_cigar_matches(cigar,start),by = .(qname2,seqnames)]
    bamData_rowWise <- bamData_rowWise[,qname1 := gsub("_[0-9]$","",qname2)]
    mergedData <- merge(blast_s,bamData_rowWise,by.x = "qname1",by.y = "qname1",all.y = TRUE,all.x = FALSE)    
    mergedData <- mergedData[!is.na(ScientificName)]
    mergedData <- cbind("Sample" = x,mergedData)  

    mergedData <- mergedData[,width := match_end - match_start]

    mergedData <- mergedData[,`:=` (frag_start = min(match_start),
                                      frag_end   = max(match_end)),by = qname1]

    library(IRanges)
    splitData <- split(mergedData, by = "seqnames", keep.by = TRUE)

    # apply per-seqnames
    coverage_list <- lapply(splitData, function(dt) {
      # within each seqnames you may have multiple taxID/fragment ranges
      dt[, {
        rng <- IRanges(start = match_start, end = match_end)
        cov <- coverage(rng)

        pos <- frag_start[1]:frag_end[1]

        # extend coverage to full fragment if necessary
        if (frag_end[1] > length(cov)) {
          cov <- c(cov, Rle(0L, frag_end[1] - length(cov)))
        }

        depth <- as.integer(cov[pos])
        .(pos = pos, depth = depth)
      }, by = .(Sample,ScientificName, qname1, frag_start, frag_end)]
    })

    # bind back together
    coverageDT <- rbindlist(coverage_list, idcol = "seqnames")

    temp <- unique(coverageDT[depth >= 1][,.(Sample,seqnames,ScientificName,pos)])[,.(numBases = .N),by = .(ScientificName,seqnames)]
    temp <- merge(temp,unique(abundantGeneList),
                        by.x = "seqnames",by.y = "ReferenceSequence",
                        all.x = TRUE,all.y = FALSE,
                        allow.cartesian = TRUE)

    temp <- temp[,coverage := (numBases/ReferenceLength) * 100 ]

    fwrite(temp,
          paste0("/lustre/samhita.p/tam_amr_dengue80/results/blastBased_V3/",x,"_org_Coverage.tsv"),
          quote = FALSE, row.names = FALSE,sep = "\t")
      
})
