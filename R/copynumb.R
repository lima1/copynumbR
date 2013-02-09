copynumbR.gistic.read.lesions <- structure(function
### Read GISTIC all_lesions file
(filename, 
### The filename of the GISTIC output file
clinical,
### A data frame with clinical annotation for the phenoData slot of the output
### ExpressionSet
data.col=10, 
### Start column of the sample data, no need to change unless GISTIC output
### changes 
...
### Additional arguments passed to copynumbR.eset()
) {
    data <- read.delim(filename, stringsAsFactors=FALSE)
    data <- data[grep("^Actual Copy", data[,9]),]
    copynumbR.eset(data, clinical,data.col,...)
### ExpressionSet containing the significant GISTIC output    
},ex=function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <-
    copynumbR.gistic.read.lesions(system.file("extdata/gistic_stransky_bladder",
        "all_lesions.conf_95.txt", package="copynumbR"), clinical)

    # list all recurrent gains and losses
    featureData(eset)$Descriptor
    
    # show the copy number distributions of all recurrent alterations
    boxplot(t(exprs(eset)))
})

copynumbR.gistic.read.genes <- structure(function
### Read GISTIC all_data_by_genes file
(filename, 
### The filename of the GISTIC output file
clinical,
### A data frame with clinical annotation for the phenoData slot of the output
### ExpressionSet
data.col=4, 
### Start column of the sample data, no need to change unless GISTIC output
### changes 
...
### Additional arguments passed to copynumbR.eset()
) {
    data <- read.delim(filename, stringsAsFactors=FALSE)
    eset <- copynumbR.eset(data, clinical,data.col,...)
    featureNames(eset) <- featureData(eset)$Gene.Symbol
    eset
### ExpressionSet containing the copy numbers for all genes as estimated by
### GISTIC
},ex=function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <-
    copynumbR.gistic.read.genes(system.file("extdata/gistic_stransky_bladder",
        "all_data_by_genes.txt", package="copynumbR"), clinical)

    # show the copy number distribution of MYC copy numbers
    boxplot(exprs(eset)["MYC",])
})

copynumbR.gistic.write.arrayfile <- function
### Write GISTIC input arrayfile 
(labels, 
### The sample names as used in the segmented file
file
### The output filename of this function
) {
    write.table(data.frame(Array=labels),file=file, row.names=FALSE,quote=FALSE)
}


copynumbR.gistic.read.broad <- structure(function
### Read GISTIC broad_values_by_arm.txt output file
(filename="broad_values_by_arm.txt",
### The filename of the GISTIC output file
clinical,
### A data frame with clinical annotation for the phenoData slot of the output
### ExpressionSet
data.col=2, 
### Start column of the sample data, no need to change unless GISTIC output
### changes 
...
### Additional arguments passed to copynumbR.eset()
) {
    data = read.delim(filename, stringsAsFactors=FALSE)
    copynumbR.eset(data, clinical,data.col,...)
},ex=function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <-
    copynumbR.gistic.read.broad(system.file("extdata/gistic_stransky_bladder",
        "broad_values_by_arm.txt", package="copynumbR"), clinical)
    
    # display the distributions of copy numbers of all chromosome arms
    boxplot(t(exprs(eset)))
})

copynumbR.gistic.genes.region <- function
### Extract GISTIC genes in a specified region
(filename, 
### The GISTIC output file amp_genes.conf_95.txt
region
### The region to extract. 
) {
    region <- gsub("\\(probes.*$","", region)
    data <- read.delim(filename, stringsAsFactors=FALSE)
    idx <- match(region, data[3,])
    genes <- data[4:nrow(data),idx]
    genes <- genes[genes!=""]
    genes
### Genes in the specified region    
}

copynumbR.gistic.genes.band <- structure(function
### Extract the GISTIC target genes 
(filename, 
### The GISTIC output file amp_genes.conf_95.txt
band
### The chromosome band
) {
    band <- make.names(band)
    data <- read.delim(filename, stringsAsFactors=FALSE)
    lapply(data[3, grepl(band,colnames(data),fixed=TRUE)], function(region)
    copynumbR.gistic.genes.region(filename,region))
### Genes in the specified chromosome band    
},ex=function(){
    library(copynumbR)
    # extract the gene symbols in the GISTIC peak 8q24.21
    band <-
    copynumbR.gistic.genes.band(system.file("extdata/gistic_stransky_bladder",
        "amp_genes.conf_95.txt", package="copynumbR"), "8q24.21")
})

copynumbR.gistic.focal <- structure(function
### Create a data.frame with all GISTIC focal alterations. Useful for
### presentation in knitR/Sweave.
(eset, 
### The GISTIC lesions file read with copynumbR.gistic.read.lesions
gistic.lesions.file.amp="amp_genes.conf_95.txt",
### The GISTIC output file amp_genes.conf_95.txt
gistic.lesions.file.del="del_genes.conf_95.txt", 
### The GISTIC output file del_genes.conf_95.txt
gain=0.1,
### Minimum log2 ratio for a copy number gain
loss=-0.1
### Maximum log2 ratio for a copy number loss
) {
        eset <- .addGISTICregion(eset)
        df <- data.frame(Chr = featureData(eset)[[2]], Start =
        featureData(eset)$start, End = featureData(eset)$end,
        Type= sapply(featureData(eset)[[1]], function(x) strsplit(x, 
        " ")[[1]][1]), "q-value"=featureData(eset)[[6]], 
        "res. q-value"=featureData(eset)[[7]], stringsAsFactors=FALSE)
        df$Freq = sapply(1:nrow(eset), function(i) ifelse(df$Type[i] ==
        "Amplification", sum(exprs(eset)[i,] > gain ), sum(exprs(eset)[i,] < loss)))
        df$Freq.F = paste(df$Freq, " (",round(df$Freq/ncol(eset)*100,digits=1), "%)", sep="")
        .addGenes(eset, df, gistic.lesions.file.amp, gistic.lesions.file.del)

},ex=function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <-
    copynumbR.gistic.read.lesions(system.file("extdata/gistic_stransky_bladder",
        "all_lesions.conf_95.txt", package="copynumbR"), clinical)

    focal <-
        copynumbR.gistic.focal(eset, 
            system.file("extdata/gistic_stransky_bladder", 
                "amp_genes.conf_95.txt", package="copynumbR"),
            system.file("extdata/gistic_stransky_bladder",
                "del_genes.conf_95.txt", package="copynumbR")
        )
})


copynumbR.gistic.armplot <- structure(function
### Plot chromosome numbers as estimated by GISTIC
(file="broad_values_by_arm.txt",
### The GISTIC output file broad_values_by_arm.txt
subtypes=NULL,
### If not NULL, then stratify by subtype
chromosomes=c(1:12,16:20),
### Display only the specified chromosomes. 
ncol=NULL, 
### Number of columns in the plot
... ){
    data <- read.delim(file,stringsAsFactors=FALSE)
    cols <- match(paste(rep(chromosomes,2), c(rep("p",length(chromosomes)),
    rep("q", length(chromosomes))), sep=""),data[,1])
    data <- data[cols,]
    data.stack <- stack(data,select=-Chromosome.Arm)
    data$Chromosome <- gsub("p|q","", data[,1])
    data$Arm        <- gsub("^\\d+","",data[,1])

    data.stack <- cbind(
        Chromosome = rep(paste("Chr",data$Chromosome),nrow(data.stack)/nrow(data)),
        Arm        = rep(data$Arm,nrow(data.stack)/nrow(data)),
        data.stack)
    data.stack$Chromosome <- factor( data.stack$Chromosome,
    levels=unique(paste("Chr",data$Chromosome)) )
    data.stack <- cbind(data.stack[data.stack$Arm=="p",-(2:3)],
        p=data.stack[data.stack$Arm=="p",3],
        q=data.stack[data.stack$Arm=="q",3])
    if (!is.null(subtypes)) {
        data.stack$subtype = subtypes
        data.stack = data.stack[complete.cases(data.stack),]
        if (length(chromosomes)==1)
            gp <-
                ggplot(data.stack, aes(p,
                q))+geom_point(alpha=0.2)+facet_wrap(~subtype,ncol=ncol)
        else    
            gp <-
                ggplot(data.stack, aes(p,
                q))+geom_point(alpha=0.2)+facet_wrap(~Chromosome+subtype,ncol=ncol)
    } else { 
        gp <-
            ggplot(data.stack, aes(p,
            q))+geom_point(alpha=0.2)+facet_wrap(~Chromosome, ncol = ncol)

    }
    list(plot=gp, data=data.stack)
    ### A list containing the ggplot2 and the data as used in ggplot2
},ex=function(){
    library(copynumbR)

    res <-
        copynumbR.gistic.armplot(file=
            system.file("extdata/gistic_stransky_bladder", 
                "broad_values_by_arm.txt", package="copynumbR")
        )

    plot(res$plot)
})


copynumbR.tcga.write.gisticinput <- function
### Create a GISTIC input file from TCGA Level 3 data
(path=".", 
### The path containing the Level_3 folder
output.seg="tcga.seg",
### The GISTIC segemented data input file
output.alf = "tcga.alf", 
### The GISTIC array input file
hg="hg18", 
### The genome version
verbose=TRUE
### Print some additional progress information
) {
    files = dir(paste(path,"Level_3/", sep="/"), full.names=TRUE)
    files = files[grep(paste("\\.",hg,sep=""), files)]
    if (verbose) cat("Reading", files,sep="\n") 
    data = lapply(files, read.delim, stringsAsFactors=FALSE)
    
    for (i in 1:length(data)) {
        write.table(data[[i]], file=output.seg, append=i!=1, col.names=i==1,
        row.names=FALSE, quote=FALSE, sep="\t") 
    }
    arrays = read.delim(file=output.seg, header=TRUE)[,1]
    write.table(data.frame(Array=levels(arrays)), file=output.alf,
    quote=FALSE, row.names=FALSE)
}

copynumbR.read.segmented <- structure(function
### Read segmented data and turn it into an ExpressionSet
(filename, 
### The filename of segmented data
clinical, 
### A data frame with clinical annotation for the phenoData slot of the output
### ExpressionSet
gene=FALSE,
### Either segments or genes as features. 
mad=0.0,
### filter low-variance segments, typically useful for clustering
geneMap=NULL,
### The gene mapping if gene==TRUE. See ?getRS.
...
### Additional arguments passed to the copynumbR.eset function
) {
    data.col <- ifelse(gene,6,4)
    cn = read.delim(filename,as.is=TRUE)
    colnames(cn) = c("ID","chrom","loc.start","loc.end","num.mark","seg.mean") 
    seg = CNSeg(cn)
    if (gene && is.null(geneMap)) {
        data(geneInfo)
        warning("geneMap not defined. Using hg18 default provided by the CNTools package")
        geneMap <- geneInfo
    }
    if (gene)  {
        rsseg <- getRS(seg, by = "gene", imput = FALSE, XY = FALSE,
            what="mean",geneMap=geneMap)
    } else {
        rsseg <- getRS(seg, by = "region", imput = FALSE, XY = FALSE,
            what="mean",geneMap=geneMap)
    }
    if (mad > 0) {
    f1 = function(x) { sum(x == 0) == 0 }
    ffun <- filterfun(f1)
    filteredrs <- genefilter(rsseg, ffun)
    data <- madFilter(filteredrs, mad)@rs
    } else { data = rsseg@rs } 
    for (i in 1:3) data[,i] = as.numeric(as.character(data[,i]))

    eset <- copynumbR.eset(data, clinical, data.col, ...)
    
    # set gene symbols as featureName
    if (gene) {
        eset <- eset[!duplicated(featureData(eset)$genename),]
        featureNames(eset) <- featureData(eset)$genename
    }
    
    eset
},ex=function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <- copynumbR.read.segmented(system.file("extdata", "stransky_bladder.glad", package="copynumbR"), clinical)
    # Plot the distribution of copy numbers of segments in all samples
    plot(exprs(eset))

    eset.genes <- copynumbR.read.segmented(system.file("extdata",
     "stransky_bladder.glad", package="copynumbR"), clinical, gene=TRUE)

    boxplot(exprs(eset.genes)["MYC",])
})

.addGISTICregion <- function(eset) {
    featureData(eset)$chr = gsub(":.*$","",featureData(eset)[[3]])
    featureData(eset)$start =
    as.numeric(gsub("chr\\d+:|-.*$","",featureData(eset)[[3]]))
    featureData(eset)$end =
    as.numeric(gsub("^.*-|\\(.*$","",featureData(eset)[[3]]))
    eset
}

.fetchGISTICcopynumber <- function(i, eset.gistic, eset.segmented) {
    fd.gc = featureData(eset.gistic[i,])
    fd.gd = featureData(eset.segmented)
    idx = as.matrix(regionOverlap(data.frame(chr=as.numeric(gsub("chr","",fd.gc$chr)), start=fd.gc$start, end =
        fd.gc$end), data.frame(chr=fd.gd[[1]], start=fd.gd[[2]], end = fd.gd[[3]])))[1,] > 0
    cn = apply(exprs(eset.segmented[idx,]),2,function(x)
    ifelse(max(x)==max(abs(x)),max(x), min(x)) )
}

copynumbR.gistic.clone.eset <- structure(function
### Extract GISTIC peak copy numbers from segmented data 
(eset.gistic,
### An ExpressionSet with GISTIC peaks read with copynumbR.gistic.read.lesions
eset.segmented)
### An ExpressionSet with segmented data, typically from another cohort,
### read with copynumbR.read.segmented. This is useful to compare copy numbers at
### GISTIC peaks across cohorts.
{
    eset.gistic = .addGISTICregion(eset.gistic)
    M = t(mapply(rbind, lapply(1:nrow(eset.gistic), .fetchGISTICcopynumber,
    eset.gistic, eset.segmented)))
    rownames(M) = rownames(exprs(eset.gistic))
    colnames(M) = colnames(exprs(eset.segmented))
    new("ExpressionSet", exprs=M, phenoData = phenoData(eset.segmented),
    featureData= featureData(eset.gistic))
### An ExpressionSet containing copy numbers from eset.segmented of the GISTIC peaks    
### in eset.gistic
},"ex"=function(){
    library(copynumbR)

    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))

    # Read GISTIC peaks
    eset.gistic <-
    copynumbR.gistic.read.lesions(system.file("extdata/gistic_stransky_bladder",
        "all_lesions.conf_95.txt", package="copynumbR"), clinical)
    
    # Read segmented data, typically from another cohort
    eset.segmented <- copynumbR.read.segmented(system.file("extdata", "stransky_bladder.glad", package="copynumbR"), clinical)
    
    # We do not have example data of a second cohort, so we just extract the
    # GISTIC peaks from the segmented data
    eset.cloned <- copynumbR.gistic.clone.eset(eset.gistic, eset.segmented)

    plot(exprs(eset.gistic)[2,],
         exprs(eset.cloned)[2,], xlab="GISTIC", ylab="Segmented",
         main="MYC Locus")

})


.addGenes <- function(eset, df, gistic.lesions.file.amp, gistic.lesions.file.del) {
    if (!is.null(gistic.lesions.file.del)) {
        file = sapply(df$Type, function(x) ifelse(x=="Deletion",
            gistic.lesions.file.del, gistic.lesions.file.amp))
        df$Genes = 
        sapply(1:nrow(eset), function(i)
            paste(copynumbR.gistic.genes.region(file[i],
            featureData(eset[i,])[[3]]),
            collapse=", ")
            )
    }
    df
}


.loadEset <- function(...) {
    warning(".loadEset is deprecated. Use copynumbR.eset instead.")
    copynumbR.eset(...)
}

copynumbR.eset <- structure(function
### Generate an ExpressionSet from a data.frame containing the data and a
### data.frame containing clinical data
(data, 
### data.frame containing genomic data (expression, copy number, ...)
clinical, 
### data.frame containing matched clinical data
data.col=10, 
### Sample data starts at this column
id.col = 1, 
### The column of the matching sample ids in the clinical data.frame
post.process.fun = NULL,
### A function to manipulate the final ExpressionSet (adding a Surv obj for
### example)
...
### Additional arguments passed to post.process.fun
) {
    edata = data[,data.col:ncol(data)]
    edata = edata[, colnames(edata) %in% as.character(clinical[[id.col]]) ]
    if (ncol(edata) == 0) {
        edata = data[,data.col:ncol(data)]
        colnames(edata) = gsub("\\.","-",colnames(edata))
        if (sum(!substr(clinical[[id.col]],1,1) == "X")==0  &&
            sum(substr(colnames(edata),1,1) == "X")==0) {
            colnames(edata) = paste("X", colnames(edata), sep="")
        }    

        idx = sapply(clinical[[id.col]], grep, colnames(edata))
        clinical = clinical[sapply(idx, function(x) length(x) > 0),]
        idx = sapply(clinical[[id.col]], grep, colnames(edata))
        # multiple samples per patient? take first
        idx = sapply(idx, function(x) x[[1]])

        edata = edata[, idx]
        
        colnames(edata) = clinical[[id.col]]    
    }
    edata.nn = edata
    edata = apply(edata,2,as.numeric)
    colnames(edata) = colnames(edata.nn)
    rownames(edata) = rownames(edata.nn)

    pd=as(data.frame(clinical[match(colnames(edata),clinical[[id.col]]),]),"AnnotatedDataFrame")
    sampleNames(pd) = colnames(edata)
    if (data.col > 1) { 
        fd=as(data.frame(data[,1:(data.col-1)]),"AnnotatedDataFrame")
        eset = new("ExpressionSet", exprs = as.matrix(edata),
            phenoData=pd,featureData=fd)
    }    
    else eset = new("ExpressionSet", exprs = as.matrix(edata),
phenoData=pd)
    
    if (!is.null(post.process.fun)) eset <- post.process.fun(eset, ...)
    eset     
# An ExpressionSet object    
},ex=function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))

    data <- read.delim(system.file("extdata/gistic_stransky_bladder",
    "broad_values_by_arm.txt", package="copynumbR"),stringsAsFactors=FALSE)
    
    # add an example post.process.fun. makes it easy to change the annotation
    .curateGender <- function(eset, ...) {
        eset$GENDER.2 <- as.factor(as.character(ifelse(eset$GENDER=="M", "male",
            "female")))
        eset
    }

    eset <- copynumbR.eset(data, clinical, data.col=2, post.process.fun=.curateGender)
    
    eset$GENDER.2
})


copynumbR.boxplot <- structure(function
### A boxplot showing correlation of copy number and expression for matched
### data
(eset.cn, 
### ExpressionSet with copy number data. The featureNames must correspond to
### the eset.expr provided next. Typically created with
### copynumbR.read.segmented(...,gene=TRUE) or
### copynumbR.gistic.read.genes(...) 
eset.expr,
### ExpressionSet with expression data. See above, the featureNames must
### correspond to the ones of eset.cn. So for Affymetrix data for example,
#### probe sets need to be collapsed (for example with the WGCNA package). 
cutoffs=c(-Inf,-1.3,-0.1,0.1,0.9,Inf),
### Copy number cutoffs 
cutoff.labels=c("High Loss","Loss","Normal","Gain","Amplification"),
### The labels of these cutoffs
probesets=NULL, 
### Display only these genes. If null, show all genes. That default is 
### obviously only useful for already filtered ExpressionSets.
min.samples=3, 
### Minimum number of samples in each cutoff category
sqrt=FALSE,
### Square root transform the data for read counts?
highlight=NULL,
### highlight some samples (not yet implemented)
highlight.labels=NULL,
### The label of these highlighted samples shown in the legend
xlab="Copy Number",
### The label of the x-axis
ylab="Expression",
### The label of the y-axis
outlier.shape=NA
### Display outliers? Passed to geom_boxplot()
) {

    if (!is.null(probesets)) {
        eset.cn <- eset.cn[probesets,]
        eset.expr <- eset.expr[probesets,]
    }
    if (all(featureNames(eset.cn) == featureNames(eset.expr))!=TRUE) {
        stop("featureNames() of eset.cn and eset.expr do not match")
    }
    res <- lapply(2:length(cutoffs), function(i)
        apply(exprs(eset.cn),1, function(x) x>cutoffs[i-1]
        & x <=cutoffs[i]))
    names(res) <- cutoff.labels    
    res <- res[lapply(res, function(x) sum(as.vector(x))) > min.samples]

    d.f <- do.call(rbind,unlist(lapply(1:length(res), function(i)
        lapply(1:nrow(eset.cn), function(j) try(data.frame(Group=names(res)[i],
        Gene=featureNames(eset.expr)[j]
        ,Expr=exprs(eset.expr)[j,res[[i]][,j]]
        ,stringsAsFactors=FALSE)))),recursive=FALSE)) 

     d.f$Expr <- as.numeric(d.f$Expr)
     d.f$Group <- factor(d.f$Group, levels=cutoff.labels[cutoff.labels %in%
        d.f$Group])

    p <- ggplot(d.f,
    aes(Group,Expr))+geom_boxplot(outlier.shape = outlier.shape )
    if (sqrt)
    p <- p + scale_y_sqrt(breaks=trans_breaks("sqrt",function(x) x ^
    2)(c(1,1:8*500)))
    p <- p +facet_wrap(~Gene)+theme(axis.text.x=element_text(angle=45,
    hjust=1))+ylab(ylab)+xlab(xlab)
    if (!is.null(highlight)) {
        warning("Not yet implemented")    
    }
    p
### A ggplot2 object.    
},ex=function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))

    eset.genes <- copynumbR.read.segmented(system.file("extdata",
        "stransky_bladder.glad", package="copynumbR"), clinical, gene=TRUE,
        geneMap=geneMap_hg17)
    
    # load matched expression data
    data(PMID17099711.GPL91_eset)
    
    # get the samples with matching expression data
    isc <- intersect(sampleNames(eset.genes),
        sampleNames(PMID17099711.GPL91_eset))

    copynumbR.boxplot(eset.genes[,isc], PMID17099711.GPL91_eset[,isc],
        probeset=c("MYC", "ADCY8"))
})


.getCentromere <- function(eset, centromere.file) {

    if (!file.exists(centromere.file)) {
        centromere.file <- switch(centromere.file,
            "hg17" = system.file("extdata", "centromere_hg17.txt",
                package="copynumbR"),
            "hg18" = system.file("extdata", "centromere_hg18.txt",
                package="copynumbR"),
            "hg19" = system.file("extdata", "centromere_hg19.txt",
                package="copynumbR"),
            stop("Centromere file not found."))
    }

    centromere <- read.delim(centromere.file, as.is=TRUE, header=FALSE)
    centromere[,1] <- gsub("chr","",centromere[,1])
    centromere <- centromere[!duplicated(centromere[,3]),]
    centromere[,1] <- gsub("X","23",centromere[,1])
    centromere[,1] <- gsub("Y","24",centromere[,1])
    centromere[,1] <- as.numeric(centromere[,1])
    centromere <- centromere[order(centromere[,1]),]
    idx.min <- sapply(1:24,function(i) min(which(centromere[,1]==i)))
    idx.max <- sapply(1:24,function(i) max(which(centromere[,1]==i)))

    centromere <- cbind(centromere[idx.min,2], centromere[idx.max,3],0)
    centromere[13:15,1:2] <- 0

    chrs <- unique(featureData(eset)[[1]])
    chrl <- sapply(chrs, function(x)
    max(featureData(eset[featureData(eset)[[1]]==x,])[[3]]))
    chrstart <- sapply(chrs, function(i) min(featureData(eset[featureData(eset)[[1]] == i,])[[2]]))
    chrl <- c(0,chrl - chrstart)


    centromere <- cbind(cumsum(chrl)[1:24]+centromere[,1],
                    cumsum(chrl)[1:24]+centromere[,2], 0)
    list(centromere=centromere, chrstart=chrstart, chrs=chrs, chrl=chrl)                
}

.getSmoothedData <- function(esets, window, labels, gain, loss, sma,
    centromere.file) {
    hg18 <- .getCentromere(esets[[1]], centromere.file)

     .coords <- function(eset,label,gain,loss) {   
        dels <- apply(exprs(eset),1,function(x) sum(x< loss))
        amps <- apply(exprs(eset),1,function(x) sum(x> gain))
        dels <- dels/ncol(eset)
        amps <- amps/ncol(eset)
        pos <- sapply(featureData(eset)[[1]], function(x) cumsum(hg18$chrl)[x]) +
        featureData(eset)[[2]] - sapply(featureData(eset)[[1]], function(x)
         hg18$chrstart[x])
        df <- data.frame(Begin=c(0,pos),field=c(0,amps), type="Gain") 
        df <- rbind(df, data.frame(Begin=pos,field=dels*-1, type="Loss"))
        df$Begin <- round(df$Begin/window)
        x <- df$field
        if (sma>1) {
            df$field <- SMA(df$field,sma+1)
            df$field[1:sma] <- x[1:sma]
        }
        df$label <- label
        df.1 <- ddply(df[df$type=="Gain",], "Begin", function(x)
            data.frame(field=median(x[,2]),type=x[1,3],label=x[1,4]))

        df.2 <- ddply(df[df$type=="Loss",], "Begin", function(x)
            data.frame(field=median(x[,2]),type=x[1,3],label=x[1,4]))
        df <- rbind(df.1, df.2)
     }
    if (length(gain) == 1) gain <- rep(gain, length(esets))
    if (length(loss) == 1) loss <- rep(loss, length(esets))

    df <- do.call(rbind, lapply(1:length(esets), function(i)
    .coords(esets[[i]],labels[[i]],gain=gain[i],loss=loss[i])))
    

     df$label<-factor(df$label,levels=labels)
     list(df=df, hg=hg18)
}


copynumbR.plot <- structure(function
### Chromosome frequency plot
(esets,
### List of ExpressionSets
window=500000,
### Window size in base pairs. Copy numbers will be the median in each window
labels=paste(names(esets), " (n = ", sapply(esets, ncol),")", sep=""),
### The labels of the ExpressionSets.
gain=0.1,
### Minimum log2 ratio for a copy number gain
loss=-0.1,
### Maximum log2 ratio for a copy number loss
sma=10,
### Smooth copy numbers with a simple moving average
from.chr=1,
### Start plotting at this chromosome
to.chr=22,
### End plotting at this chromosome
ylab="Loss / Gain",
### The y-axis label
xlab="Chromosome",
### The x-axis label
centromere.file="hg18"
### File containing the centromere locations for each chromosome
### These files are already provided for hg17-hg19
###  curl -s
### "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBand.txt.gz" |
### gunzip -c | grep acen
) {
    max.chr <- max(sapply(esets, function(X) max(featureData(X)$chrom))) 
    min.chr <- min(sapply(esets, function(X) min(featureData(X)$chrom))) 
    if (max.chr < to.chr) to.chr <- max.chr
    if (min.chr > from.chr) from.chr <- min.chr

    res <- .getSmoothedData(esets,window,labels,gain,loss,sma,centromere.file)
    
    df <- res$df

    q <- ggplot(df, aes(Begin,field,ymin=0,
        ymax=field,fill=type))+geom_ribbon()
    q <- q+ scale_fill_brewer("",palette="Set1")+labs(fill="") +
        scale_x_continuous(minor_breaks=cumsum(res$hg$chrl)/window,
        breaks=round(res$hg$centromere[from.chr:to.chr,1]/window)
        ,labels=res$hg$chrs[from.chr:to.chr],limits=c(from.chr-1,cumsum(res$hg$chrl)[to.chr+1]/window))+
        ylab(ylab)+xlab(xlab)+facet_grid(label ~ .)
    q <- q +  scale_y_continuous(minor_breaks = NULL,
        labels = percent_format())+
        theme(panel.grid.minor =
        element_line(size = 0.4, colour =
            'white'),panel.grid.major=element_line(linetype="blank"))
### A ggplot2 object
},ex=function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <- copynumbR.read.segmented(system.file("extdata", "stransky_bladder.glad", package="copynumbR"), clinical)

    # find the non-muscle-invasive samples
    idx.noninvasive <- grepl("Ta|T1", eset$T)
    
    # create a list of ExpressionSets with the low- and high stage samples
    eset.stage <- list("Non-Muscle-Invasive"=eset[,idx.noninvasive],
        "Invasive"=eset[,!idx.noninvasive])
    
    # now compare the copy numbers of these two groups
    p <- copynumbR.plot(eset.stage, centromere.file="hg17", sma=0)

    plot(p) 
})

.fillStacked <- function(X) {
    Y <- do.call(rbind, lapply(2:nrow(X), function(i) if (X[i-1,2] ==
    X[i,2] && X$alpha[i] != 0)
    data.frame(Begin=X[i-1,1]:X[i,1],values=X[i,2], ind=X[i,3],
    values_orig=X[i,4],alpha=X[i,5])
    else X[i-1,]))
    Y[!duplicated(paste(Y$Begin, Y$ind)),]
}

copynumbR.heatmap <- structure(function
### Chromosome Heatmap
(eset,
### ExpressionSet, typically created with copynumbR.read.segmented
window=1000000,
### Window size in base pairs. Copy numbers will the median in each window
gain=0.1,
### Minimum log2 ratio for a copy number gain
loss=-0.1,
### Maximum log2 ratio for a copy number loss
Colv=NULL,  
### Cluster columns.
distfun = spearman.dist, 
### The distance function if Colv is NULL
hclustfun = hclust, 
### The cluster function if Colv is NULL
from.chr=1,
### Start plotting at this chromosome
to.chr=22,
### End plotting at this chromosome
start=NULL, 
end=NULL,
ensembl=NULL, 
plot=TRUE, 
### Plot the heatmap, otherwise the dendrogram and heatmap are returned as
### list of ggplot2 objects
hide.labels=FALSE,
### Hide the sample labels on the bottom
ylab=("Chromosome"),
### The y-axis label
centromere.file="hg18"
### File containing the centromere locations for each chromosome
### These files are already provided for hg17-hg19
###  curl -s
### "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBand.txt.gz" |
### gunzip -c | grep acen
) {
    sampleNames(eset) <- make.names(sampleNames(eset))
    eset <- eset[featureData(eset)[[1]]>= from.chr & featureData(eset)[[1]] <=
        to.chr] 

    max.chr <- max(featureData(eset)$chrom)
    min.chr <- min(featureData(eset)$chrom)
    if (max.chr < to.chr) to.chr <- max.chr
    if (min.chr > from.chr) from.chr <- min.chr

    if (is.null(Colv)) {
        Colv <- hclustfun(distfun(eset))
    }
    
    if ("hclust" %in% class(Colv)) eset <- eset[,Colv$order]

    if (!is.null(start)) {
        idx <- (featureData(eset)[[1]]== from.chr &
            featureData(eset)[[2]] >= start) & 
               (featureData(eset)[[1]]== to.chr &
            featureData(eset)[[2]] <= end)
        eset <- eset[ , order(apply(exprs(eset)[idx,],2,mean),decreasing=TRUE)]
    } 

    hg18 <- .getCentromere(eset, centromere.file)

    pos <- sapply(featureData(eset)[[1]], function(x) cumsum(hg18$chrl)[x]) +
    featureData(eset)[[2]] - sapply(featureData(eset)[[1]], function(x)
        hg18$chrstart[x])
    
    df <- data.frame(Begin=pos,exprs(eset)) 
    df$Begin <- round(df$Begin/window)
    df <- ddply(df, "Begin", function(x)
        data.frame(Begin=x[1,1],t(apply(x[,2:ncol(x)],2,median))))

    tmp <- df$Begin
    if (from.chr != to.chr) {
        df$Begin <- 1:nrow(df)
        s <- median(df$Begin/tmp,na.rm=TRUE)
    }    

    df.stack <- data.frame(Begin=rep(df$Begin, ncol(df)-1), stack(df, select = -Begin))
    df.stack$values_orig <- df.stack$values
    df.stack$values <- sapply(df.stack$values, function(x) { if(x>gain)
        return(3); if(x< loss) return(1); return(2)})
    df.stack$ind <- factor(df.stack$ind, levels=sampleNames(eset))    

    colx <- c( brewer.pal(7,"Set1")[2], "white",  brewer.pal(7,"Set1")[1] )
    df.stack$values <- colx[df.stack$values]
    df.stack$alpha <- ifelse(df.stack$values=="white",0,1)

    if (from.chr == to.chr) {
        tmp <- cumsum(hg18$chrl[1:from.chr])
        genes <- .getGenesRegion(from.chr, start, end, ensembl)
        genes$transcript_start <-  genes$transcript_start + tmp
        genes$transcript_end <-  genes$transcript_end + tmp
        range.g <- abs(  max(genes$transcript_end) -
            min(genes$transcript_start)) * 0.1

        limits <-
            c((range.g+end+tmp),(start+tmp-range.g))/window
        df.stack <- df.stack[ df.stack$Begin >= limits[2] & df.stack$Begin <=
                            limits[1], ]
        
        df.stack <- .fillStacked(df.stack)

        df.stack$alpha[ df.stack$alpha==1 ] <- df.stack$values_orig[df.stack$alpha==1]/max( abs(
         df.stack$values_orig ) )

    }    

    if (from.chr != to.chr) {     
        p1 <- ggplot(df.stack, aes(x=ind, y=Begin, alpha=alpha,
            fill=values))+geom_tile()+theme(axis.text.x = element_text(angle =
            90, hjust = 1))
        p1 <- p1 +     
            scale_y_continuous(minor_breaks=cumsum(hg18$chrl)/window*s,
        breaks=round(hg18$centromere[from.chr:to.chr,1]/window*s)
        ,labels=hg18$chrs, trans=reverse_trans(),
        limits=c(cumsum(hg18$chrl)[to.chr+1]/window*s,from.chr-1))
    } else {
        p1 <- ggplot(df.stack, aes(x=ind, y=Begin, alpha=alpha,
            fill=values))+geom_tile()+theme(axis.text.x = element_text(angle =
            90, hjust = 1))
        p1 <- p1 +     
            scale_y_continuous(
            minor_breaks=cumsum(hg18$chrl)/window,
            breaks=genes$transcript_start/window,
            labels=genes$hgnc_symbol, trans=reverse_trans(),
            limits=limits)
    }
    p1 <- p1 + xlab("")+
    ylab(ylab)+ scale_fill_identity(labels=1:3,
    breaks=colx)+theme_classic2()+theme(axis.text.x = element_text(angle = 90,
    hjust = 1),legend.position="none")
    if (!("hclust" %in% class(Colv))) {
        if (plot) print(p1)
        return(list(p1))
    }

    Colv$labels <- ""
    p2 <- ggdendrogram(dendro_data(Colv))

    if (plot) {
        tmp <- theme(plot.margin = unit(rep(0.1,4), "lines"))
        xtmp <- tmp
        if (hide.labels) xtmp <- tmp +
            theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

        grid.newpage()
        print(p2+tmp, vp=viewport(0.8, 0.2, x=0.42, y=0.92))
        print(p1+xtmp, vp=viewport(0.8, 0.84, x=0.4, y=0.47))
    }
    list(p1,p2)
### List of two ggplot2 objects (dendrogram and heatmap)    
},ex=function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <- copynumbR.read.segmented(system.file("extdata", "stransky_bladder.glad", package="copynumbR"), clinical)
    p <- copynumbR.heatmap(eset, centromere.file="hg17")
})


theme_classic2 <- function
### A simple ggplot2 theme, no fancy stuff
(base_size = 12, 
### Theme base font size
base_family = ""
### Theme base font family
) 
{
    theme_grey(base_size = base_size, base_family = base_family) %+replace% 
        theme(axis.text = element_text(size = rel(0.8)), 
              axis.ticks = element_line(colour = "black"), 
              legend.key = element_rect(colour = "grey80"), 
              panel.background = element_rect(fill = "white", colour = NA), 
              panel.border = element_rect(fill = NA, colour = "black", size=1.0), 
              panel.grid.major = element_line(linetype="blank", size = 0.2), 
              panel.grid.minor = element_line(linetype = "blank", size = 0.5), 
              panel.margin     = unit(0.75, "lines"),
              strip.background = element_rect(fill = NA, colour = "black",linetype="blank"), 
              strip.text.x = element_text(face="bold"))
}

.getGenesRegion <- function(chromosome, start, end, ensembl=NULL) {
    if (is.null(ensembl)) {
       ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    }
    x <-getBM(c("hgnc_symbol", "chromosome_name", "transcript_start",
        "transcript_end"), values=chromosome, filters="chromosome_name", ensembl)
    x <- x[x$transcript_start > start & x$transcript_end < end &
    x$hgnc_symbol != "",]
    x[!duplicated(x$hgnc_symbol),] 
}

