#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DEXSeq))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
library(magrittr)


####################################
#### parse command line options ####
####################################

option_list = list(
    make_option(c("-a", "--annotation"), type="character", default=NULL,
                help="File path to the annotation file, same format as the input for the count_exons.R script.",
                metavar="character"),
    make_option( c("-c", "--contrastFile"), type="character", default=NULL, 
                help="A contrast file specifying the differential contrast. It should be a tab separated file with at least three unnamed columns: (1) sample name (2) corresponding file resulting from the `count_exons.R` script, (3) conditions, first in line is taken as baseline. (4-n) Other covariates.",
                metavar="character" ),
    make_option(c("-o", "--outputPrefix"), type="character", default=NULL,
                help="A path with a prefix for the output file with counts.",
                metavar="character"),
    make_option(c("-n", "--ncpus"), type="integer", default=1,
                help="A numeric value specifying the number of CPUs to use.",
                metavar="integer"),
    make_option(c("-t", "--testEngine"), type="character", default="dexseq",
                help="A character string: 'dexseq' or 'diffsplice' that specifies the statistical engine.",
                metavar="character"))
opt_parser = OptionParser( option_list=option_list )
opt = parse_args(opt_parser)

annotationFile <- opt[["annotation"]]
contrastFile <- opt[["contrastFile"]]
outputPrefix <- opt[["outputPrefix"]]
nCores <- opt[["ncpus"]]
testEngine <- opt[["testEngine"]]

if( nCores > 1 ){
    bpparam <- bpstart(MulticoreParam(nCores))
}else{
    bpparam <- bpstart(SerialParam())
}

## check exon count input validity and arrange matrices

contrastDF <- read.delim( contrastFile, header=FALSE, sep="\t" )
colnames(contrastDF)[1:3] <- c("sample_name", "count_file", "condition")
rownames(contrastDF) <- contrastDF$sample_name
existFlag <- file.exists( contrastDF[["count_file"]] )

if( any(!existFlag) ){
    stop(sprintf("The following input file(s) could not be found:\n\t%s", paste( contrastDF[["count_file"]][!existFlag], collapse=", " ) ))
}
countsPerFile <- bplapply( contrastDF[["count_file"]], function(x){
    xx <- read.delim(x)
    stopifnot(all(grepl("binID|exoncount|skippingcount|psi", colnames(xx))))
    colnames(xx) <- c("binID", "exoncount", "skippingcount", "psi")
    xx
}, BPPARAM=bpparam )

if( length( unique(sapply(countsPerFile, nrow)) ) >  1 )
    stop("The count files don't have the same number of rows")

names(countsPerFile) <- contrastDF[["sample_name"]]

for( i in seq_len(length(countsPerFile)) ){
    if( i > 1 ){
        if( !all(countsPerFile[[i]]$binID == allBinIDs) ){
            stop("The exon coordinates of the count files are different")
        }
    }
    allBinIDs <- countsPerFile[[i]]$binID
}
exonCountMatrix <- vapply( countsPerFile, function(x){x$exoncount}, numeric(length(allBinIDs)) )
skippingCountMatrix <- vapply( countsPerFile, function(x){x$skippingcount}, numeric(length(allBinIDs)) )
psiMatrix <- vapply( countsPerFile, function(x){x$psi}, numeric(length(allBinIDs)) )


is.constant.factor <- sapply(contrastDF, function(x) length(unique(x)) == 1)
if ( length(is.constant.factor[is.constant.factor]) > 0) {
    sprintf("Filtering out columns:'%s' from a contrast definition because they are constant", paste(names(is.constant.factor)[is.constant.factor],collapse=" "))
}

contrastDF <- contrastDF[!is.constant.factor]
## validate annotation
covariates <- colnames(contrastDF)[!colnames( contrastDF ) %in% c("sample_name", "count_file", "condition")]

fullDesign <- as.formula(sprintf( "~ sample + exon + %scondition:exon", paste(sprintf( "%s:exon + ", covariates), collapse="" )))
reducedDesign <- as.formula(sprintf( "~ sample + exon%s", paste(sprintf( " + %s:exon", covariates), collapse="" )))
print(paste("Full design is", format(fullDesign)))
print(paste("Reduced design is", format(reducedDesign)))
if( !file.exists(annotationFile) )
    stop("The annotation file could not be found")
    
annotationData <- read.delim( annotationFile )

if( !all(c("binID", "gene_id", "tx_name", "annotated", "isIntron") %in% colnames( annotationData )) )
    stop("There are columns missing from the annotation file.")

annotationGRanges <- GRanges( annotationData$binID )
mcols(annotationGRanges) <- DataFrame(annotationData[,names(annotationData) %in% c("gene_id", "tx_name", "annotated", "isIntron", "gene_name")])
annotationGRanges$tx_name <- CharacterList(strsplit( annotationGRanges$tx_name, "\\+" ))

mcols(annotationGRanges)$exonicPartID <- sprintf("E%.3d", as.numeric(ave(annotationGRanges$gene_id, annotationGRanges$gene_id, FUN=seq_along) ))

txList <- as.list(annotationGRanges$tx_name)
annotationGRanges$tx_name <- NULL
contrastDF$count_file <- NULL

dxd <- DEXSeqDataSet( exonCountMatrix, sampleData=contrastDF,
              design=fullDesign,
              featureID=mcols(annotationGRanges)$exonicPartID,
              groupID=mcols(annotationGRanges)$gene_id,
              featureRanges=annotationGRanges,
              transcripts=txList,
              alternativeCountData=skippingCountMatrix )

dxd <- estimateSizeFactors( dxd )
colData( dxd )$condition <- relevel( colData( dxd )$condition, contrastDF$condition[1] )
keep1 <- rowSums(featureCounts(dxd) > 2) > ncol(dxd)/2/3
keep2 <- rowSums(counts(dxd)[,colData(dxd)$exon == "others"] > 2) > ncol(dxd)/2/3


dxd_sub <- dxd[keep2&keep1,]

if( !testEngine %in% c("diffsplice", "dexseq") )
    stop("The parameter 'testEngine' must be set to either 'diffsplice' or 'dexseq'")

if( testEngine == "dexseq" ){
    dxd_sub <- estimateSizeFactors(dxd_sub)
    dxd_sub <- estimateDispersions( dxd_sub, BPPARAM=bpparam )
    dxd_sub <- testForDEU( dxd_sub, reducedModel=reducedDesign, BPPARAM=bpparam )
    dxd_sub <- estimateExonFoldChanges( dxd_sub, maxRowsMF=0, denominator=contrastDF$condition[1], BPPARAM=bpparam )
    bpstop( bpparam )
    dxd_res <- DEXSeqResults( dxd_sub, independentFiltering=FALSE )
    lvs <- levels(colData(dxd_sub)$condition)
    colnames(dxd_res)[grepl("log2fold_", colnames(dxd_res))] <- "log2fc"
    dxd_res <- dxd_res[,c("log2fc", "pvalue", "padj")]
    dxd_res <- dxd_res %>%
        as.data.frame %>%
        tibble::rownames_to_column("exonID")
    names( annotationGRanges ) <- paste( mcols( annotationGRanges )$gene_id, mcols(annotationGRanges)$exonicPartID, sep=":" )
    coors <- annotationGRanges[dxd_res$exonID]
    idCoors <- sprintf("%s:%s-%s:%s", seqnames(coors), start(coors), end(coors), strand(coors) )
    dxd_res$coorsID <- idCoors
    dxd_res <- dxd_res %>%
        dplyr::rename( altsplice_log2fc=log2fc, altsplice_pvalue=pvalue, altsplice_padj=padj )
}else if( testEngine == "diffsplice" ){
    thisCols <- colData(dxd_sub)$exon == "this"
    thisCounts <- counts(dxd_sub)[,thisCols]
    otherCounts <- counts(dxd_sub)[,!thisCols]
    stopifnot(all(colData(dxd_sub)$sample[thisCols] == colData(dxd_sub)$sample[!thisCols]))
    rowData(dxd_sub)$exonID <- rownames(rowData(dxd_sub))
    dff <- data.frame(
        gene_id = rep( rowData( dxd_sub )$exonID, 2 ),
        isoform_id=paste( rowData(dxd_sub)$exonID, rep( c( "this", "others" ), each=nrow(thisCounts) ), sep=":") )
    rownames(dff) <- dff$isoform_id
    stopifnot(all(table(dff$gene_id) == 2))
    dge <- DGEList(counts=rbind( thisCounts, otherCounts ))
    dge <- calcNormFactors(dge)
    dge$genes$GeneID <- dff$gene_id
    dge$genes$ExonID <- dff$isoform_id
    dge$genes <- as.data.frame(dge$genes)
    mf <- as.data.frame( sampleAnnotation( dxd_sub ) )
    mf$condition <- relevel( mf$condition, contrastDF$condition[1] )
    fullDesign <- paste("~", gsub("sample \\+ exon \\+ |:exon", "", as.character( fullDesign )[2]))
    fullDesign <- as.formula(fullDesign)
    design <- model.matrix( fullDesign, mf )
    v <- voom( dge, design, plot=FALSE )
    fit <- lmFit( v, design )
    ex <- diffSplice( fit, geneid=dff$gene_id, exonid=dff$isoform_id )
    topT <- topSplice(ex, coef=ncol(design), test="t", number=Inf)
    topT <- topT[grepl("this$", topT$ExonID),]
    names( annotationGRanges ) <- paste( mcols( annotationGRanges )$gene_id, mcols(annotationGRanges)$exonicPartID, sep=":" )
    coors <- annotationGRanges[topT$GeneID]
    idCoors <- sprintf("%s:%s-%s:%s", seqnames(coors), start(coors), end(coors), strand(coors) )
    topT$coorsID <- idCoors
    dxd_res <- topT %>%
        dplyr::rename( exonID=GeneID, altsplice_log2fc=logFC,
                      altsplice_pvalue=`P.Value` ) %>%
        dplyr::mutate( altsplice_padj=p.adjust( altsplice_pvalue, method="BH" ) ) %>%
        dplyr::select( exonID, dplyr::matches("altsplice"), coorsID )
}

allDexseqGR <- rowRanges(dxd)
dexseqGR <- GRanges(dxd_res$coorsID)
dexseqOvl <- findOverlaps( allDexseqGR, dexseqGR, type="equal" )

resultsDF <- data.frame(
    coorsID=sprintf( "%s:%d-%d:%s", as.character( seqnames(allDexseqGR) ),
                    start( allDexseqGR ), end( allDexseqGR ), as.character( strand( allDexseqGR ) ) ) )

allDexseqGR.mtdt.df <- as.data.frame(mcols( allDexseqGR ))
resultsDF <- cbind(resultsDF, allDexseqGR.mtdt.df[,names(allDexseqGR.mtdt.df) %in% c("gene_id", "annotated", "isIntron", "gene_name")])

for( i in c("altsplice_log2fc", "altsplice_pvalue", "altsplice_padj") ){
    resultsDF[[i]] <- NA
    resultsDF[[i]][queryHits(dexseqOvl)] <- dxd_res[[i]][subjectHits(dexseqOvl)]
}

psiMatrix <- lapply(split(colData(dxd)$sample_name, colData(dxd)$condition),
       function(x){
           rowMeans( psiMatrix[,x] )
       })

deno <- names(psiMatrix)[names(psiMatrix) %in% contrastDF$condition[1]]
num <- names(psiMatrix)[!names(psiMatrix) %in% contrastDF$condition[1]]
resultsDF$altsplice_delta_psi <- psiMatrix[[num]]-psiMatrix[[deno]]

## remove exons/introns without data
resultsDF <- resultsDF %>%
    dplyr::filter( !(is.na(altsplice_log2fc) & is.na(altsplice_pvalue) & is.na(altsplice_padj) & is.na(altsplice_delta_psi) ) )
resultsDF <- resultsDF %>%
    dplyr::filter( !(is.na(altsplice_log2fc) & is.na(altsplice_pvalue) & is.na(altsplice_padj) & altsplice_delta_psi == 0 ) )

write.table(
    resultsDF, file=sprintf("%s_diff.txt", outputPrefix),
    quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE )

saveRDS( dxd, file=sprintf("%s.rds", outputPrefix)  )
