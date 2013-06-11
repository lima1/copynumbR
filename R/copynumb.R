copynumbR.gistic.read.lesions <- structure(function
### Read GISTIC all_lesions file.
(filename, 
### The filename of the GISTIC all_lesions.conf_95.txt output file.
clinical,
### A data frame with clinical annotation for the phenoData slot of the output
### ExpressionSet.
data.col=10, 
### Start column of the sample data, no need to change unless GISTIC output
### changes.
...
### Additional arguments passed to copynumbR.eset().
) {
    data <- read.delim(filename, stringsAsFactors=FALSE)
    data <- data[grep("^Actual Copy", data[,9]),]
    eset <- copynumbR.eset(data, clinical,data.col,...)
    # parse the wide peak region and make it accessible in the 
    # featureData slot
    .addGISTICregion(eset)
### ExpressionSet containing the significant GISTIC focal peaks.
},ex=function(){
    library(copynumbR)

    clinical <- read.csv(system.file("extdata", 
        "stransky_bladder_clinical.csv", package="copynumbR"))

    eset <-
    copynumbR.gistic.read.lesions(system.file("extdata/gistic_stransky_bladder",
        "all_lesions.conf_95.txt", package="copynumbR"), clinical)

    # list all recurrent gains and losses
    featureData(eset)$Descriptor
    
    # show the copy number distributions of all recurrent alterations
    boxplot(t(exprs(eset)))
})

copynumbR.gistic.read.genes <- structure(function
### Read GISTIC all_data_by_genes file.
(filename, 
### The filename of the GISTIC all_data_by_genes.txt output file.
clinical,
### A data frame with clinical annotation for the phenoData slot of the output
### ExpressionSet.
data.col=4, 
### Start column of the sample data, no need to change unless GISTIC output
### changes.
...
### Additional arguments passed to copynumbR.eset().
) {
    data <- read.delim(filename, stringsAsFactors=FALSE)
    eset <- copynumbR.eset(data, clinical,data.col,...)
    featureNames(eset) <- featureData(eset)$Gene.Symbol
    eset
### ExpressionSet containing the copy numbers for all genes as estimated by
### GISTIC.
},ex=function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", 
        "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <-
    copynumbR.gistic.read.genes(system.file("extdata/gistic_stransky_bladder",
        "all_data_by_genes.txt", package="copynumbR"), clinical)

    # show the copy number distribution of MYC copy numbers
    boxplot(exprs(eset)["MYC",])
})

copynumbR.gistic.write.arrayfile <- function
### Write GISTIC input arrayfile. 
(labels, 
### The sample names as used in the segmented file.
file
### The output filename of this function.
) {
    write.table(data.frame(Array=labels),file=file, row.names=FALSE,quote=FALSE)
}


copynumbR.gistic.read.broad <- structure(function
### Read GISTIC broad_values_by_arm.txt output file.
(filename="broad_values_by_arm.txt",
### The filename of the GISTIC broad_values_by_arm.txt output file.
clinical,
### A data frame with clinical annotation for the phenoData slot of the output
### ExpressionSet.
data.col=2, 
### Start column of the sample data, no need to change unless GISTIC output
### changes.
...
### Additional arguments passed to copynumbR.eset().
) {
    data <- read.delim(filename, stringsAsFactors=FALSE)
    copynumbR.eset(data, clinical,data.col,...)
},ex=function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", 
        "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <-
    copynumbR.gistic.read.broad(system.file("extdata/gistic_stransky_bladder",
        "broad_values_by_arm.txt", package="copynumbR"), clinical)
    
    # display the distributions of copy numbers of all chromosome arms
    boxplot(t(exprs(eset)))
})

copynumbR.gistic.genes.region <- function
### Extract GISTIC genes in a specified region.
(filename, 
### The filename of the GISTIC output file amp_genes.conf_95.txt.
region
### The region to extract. 
) {
    region <- gsub("\\(probes.*$","", region)
    data <- read.delim(filename, stringsAsFactors=FALSE)
    idx <- match(region, data[3,])
    genes <- data[4:nrow(data),idx]
    genes <- genes[genes!=""]
    genes
### Genes in the specified region.    
}

copynumbR.gistic.genes.band <- structure(function
### Extract the GISTIC target genes. 
(filename, 
### The filename of the GISTIC output file amp_genes.conf_95.txt.
band
### The chromosome band.
) {
    band <- make.names(band)
    data <- read.delim(filename, stringsAsFactors=FALSE)
    lapply(data[3, grepl(band,colnames(data),fixed=TRUE)], function(region)
    copynumbR.gistic.genes.region(filename,region))
### Genes in the specified chromosome band.
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
### The GISTIC lesions file read with copynumbR.gistic.read.lesions.
gistic.lesions.file.amp="amp_genes.conf_95.txt",
### The GISTIC output file amp_genes.conf_95.txt.
gistic.lesions.file.del="del_genes.conf_95.txt", 
### The GISTIC output file del_genes.conf_95.txt.
gain=0.1,
### Minimum log2 ratio for a copy number gain.
loss=-0.1
### Maximum log2 ratio for a copy number loss.
) {
        df <- data.frame(
            Chr = featureData(eset)$Descriptor, 
            Start = featureData(eset)$start, 
            End = featureData(eset)$end,
            Type= sapply(featureData(eset)[[1]], function(x) strsplit(x, 
                " ")[[1]][1]), 
            "q-value"=featureData(eset)[[6]], 
            "res. q-value"=featureData(eset)[[7]], stringsAsFactors=FALSE
        )

        df$n <- sapply(1:nrow(eset), function(i) ifelse(df$Type[i] ==
        "Amplification", sum(exprs(eset)[i,] > gain ), sum(exprs(eset)[i,] < loss)))
        df$Freq <- paste(df$n, " (",round(df$n/ncol(eset)*100,digits=1), "%)", sep="")
        .addGenes(eset, df, gistic.lesions.file.amp, gistic.lesions.file.del)
### A data.frame containing all GISTIC focal amplifcations. 
},ex=function(){
    library(copynumbR)

    clinical <- read.csv(system.file("extdata", 
        "stransky_bladder_clinical.csv", package="copynumbR"))

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
    data <- data[cols,-ncol(data)]
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
    gp
    ### A ggplot2 object.
},ex=function(){
    library(copynumbR)

    p <-
        copynumbR.gistic.armplot(file=
            system.file("extdata/gistic_stransky_bladder", 
                "broad_values_by_arm.txt", package="copynumbR")
        )

    plot(p)
})

copynumbR.read.segmented <- structure(function
### Read segmented data and turn it into an ExpressionSet.
(filename, 
### The filename of segmented data.
clinical, 
### A data frame with clinical annotation for the phenoData slot of the output
### ExpressionSet.
gene=FALSE,
### Either segments or genes as features. 
mad=0.0,
### Filter low-variance segments, typically useful for clustering.
### Calculates the mean absolute deviation across samples
### for each rows and drops rows that are not above the percentile defined.
geneMap=NULL,
### The gene mapping if gene==TRUE. See ?getRS. Maps for hg17, hg18 and hg19
### are provided by this package (geneMap_hgXX).
genes.as.features=TRUE,
### If gene==TRUE, use genes as feature names. Note that feature names
### have to be unique and only the first occurence of a gene in the geneMap
### is used if genes.as.features==TRUE.
...
### Additional arguments passed to the copynumbR.eset function.
) {
    data.col <- ifelse(gene,6,4)
    cn <- read.delim(filename,as.is=TRUE)
    colnames(cn) <- c("ID","chrom","loc.start","loc.end","num.mark","seg.mean") 
    seg <- CNSeg(cn)
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
        f1 <- function(x) { sum(x == 0) == 0 }
        ffun <- filterfun(f1)
        filteredrs <- genefilter(rsseg, ffun)
        data <- madFilter(filteredrs, mad)@rs
    } else { 
        data = rsseg@rs 
    }

    for (i in 1:3) data[,i] = as.numeric(as.character(data[,i]))

    eset <- copynumbR.eset(data, clinical, data.col, ...)
    
    # set gene symbols as featureName
    if (gene && genes.as.features) {
        eset <- eset[!duplicated(featureData(eset)$genename),]
        featureNames(eset) <- featureData(eset)$genename
    }
    
    eset
### An ExpressionSet.    
},ex=function(){
    library(copynumbR)

    clinical <- read.csv(system.file("extdata", 
        "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <- copynumbR.read.segmented(system.file("extdata", 
        "stransky_bladder.glad", package="copynumbR"), clinical)

    # Plot the distribution of copy numbers of segments in all samples
    plot(exprs(eset))

    eset.genes <- copynumbR.read.segmented(system.file("extdata",
     "stransky_bladder.glad", package="copynumbR"), clinical, gene=TRUE,
     geneMap=geneMap_hg17)

    boxplot(exprs(eset.genes)["MYC",])
})

.addGISTICregion <- function(eset) {
  
    peak <- featureData(eset)$Wide.Peak.Limits
    if (is.null(peak)) {
        warning("Could not find Wide.Peak.Limits slot.")
        return(eset)
    }

    featureData(eset)$chr <- gsub(":.*$","", peak)
    featureData(eset)$start <-
        as.numeric(gsub("chr\\d+:|-.*$","",peak))
    featureData(eset)$end <-
        as.numeric(gsub("^.*-|\\(.*$","",peak))
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
### Extract GISTIC peak copy numbers from segmented data.
(eset.gistic,
### An ExpressionSet with GISTIC peaks read with
### copynumbR.gistic.read.lesions().
eset.segmented
### An ExpressionSet with segmented data, typically from another cohort,
### read with copynumbR.read.segmented. This is useful to compare copy numbers at
### GISTIC peaks across cohorts.
){
    
    M <- t(mapply(rbind, lapply(1:nrow(eset.gistic), .fetchGISTICcopynumber,
        eset.gistic, eset.segmented)))

    rownames(M) <- rownames(exprs(eset.gistic))
    colnames(M) <- colnames(exprs(eset.segmented))

    new("ExpressionSet", exprs=M, phenoData = phenoData(eset.segmented),
    featureData= featureData(eset.gistic))
### An ExpressionSet containing copy numbers from eset.segmented of the GISTIC peaks    
### in eset.gistic.
},"ex"=function(){
    library(copynumbR)

    clinical <- read.csv(system.file("extdata", 
        "stransky_bladder_clinical.csv", package="copynumbR"))

    # Read GISTIC peaks
    eset.gistic <-
    copynumbR.gistic.read.lesions(system.file("extdata/gistic_stransky_bladder",
        "all_lesions.conf_95.txt", package="copynumbR"), clinical)
    
    # Read segmented data, typically from another cohort
    eset.segmented <- copynumbR.read.segmented(system.file("extdata", 
        "stransky_bladder.glad", package="copynumbR"), clinical)
    
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
### data.frame containing clinical data.
(data, 
### A data.frame containing genomic data (expression, copy number, ...).
clinical, 
### A data.frame containing matched clinical data.
data.col=2, 
### Sample data starts at this column.
id.col = 1, 
### The column of the matching sample ids in the clinical data.frame.
post.process.fun = NULL,
### A function to manipulate the final ExpressionSet (adding a Surv obj for
### example).
...
### Additional arguments passed to post.process.fun().
) {
    edata <- data[,data.col:ncol(data)]
    edata <- edata[, colnames(edata) %in% as.character(clinical[[id.col]]) ]
    
    # Try to be a little bit smart, this is especially useful for matching
    # TCGA barcodes. So if the colnames of the data do not match, try a few
    # things before giving up:
    if (ncol(edata) == 0) {
        warning(paste("colnames of data and sample ids in clinical data do",
        "not match. Trying some fuzzy matching."))
        edata <- data[,data.col:ncol(data)]

        # maybe the colnames got converted to valid R names?
        colnames(edata) <- make.names(colnames(edata))
        clinical[[id.col]] <- make.names(clinical[[id.col]])

        idx <- sapply(clinical[[id.col]], grep, colnames(edata))
        clinical <- clinical[sapply(idx, function(x) length(x) > 0),]
        idx <- sapply(clinical[[id.col]], grep, colnames(edata))
        # multiple samples per patient? take first
        idx <- sapply(idx, function(x) x[[1]])

        edata <- edata[, idx]
        
        colnames(edata) <- clinical[[id.col]]    
    }
    edata.nn <- edata
    edata <- apply(edata,2,as.numeric)
    colnames(edata) <- colnames(edata.nn)
    rownames(edata) <- rownames(edata.nn)

    pd <- as(data.frame(clinical[match(colnames(edata),clinical[[id.col]]),]),"AnnotatedDataFrame")
    sampleNames(pd) <- colnames(edata)

    if (data.col > 1) { 
        fd <- as(data.frame(data[,1:(data.col-1)]),"AnnotatedDataFrame")
        eset <- new("ExpressionSet", exprs = as.matrix(edata),
            phenoData=pd,featureData=fd)
    }    
    else eset <- new("ExpressionSet", exprs = as.matrix(edata),
phenoData=pd)
    
    if (!is.null(post.process.fun)) eset <- post.process.fun(eset, ...)
    eset     
# An ExpressionSet.
},ex=function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", 
        "stransky_bladder_clinical.csv", package="copynumbR"))

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

copynumbR.categorical <- structure(function
### Convert log2 copy number ratios in copy number categories
(x,
### Vector of log2 copy number ratios.
cutoffs=c(-Inf,-1.3,-0.1,0.1,0.9,Inf),
### Copy number cutoffs.
cutoff.labels=c("Homozyg. Deletion","Heterozyg. Deletion",
"Normal","Gain","Amplification")
### The labels of these cutoffs.
### Provide a summary
) {
    .hf <- function(xx) { 
        for (i in 2:length(cutoffs)) 
            if (xx > cutoffs[i-1] &
                xx <= cutoffs[i]) return(cutoff.labels[i-1])
    }
    sapply(x,.hf)
})

copynumbR.boxplot <- structure(function
### A boxplot showing correlation of copy number and expression for matched
### data.
(eset.cn, 
### ExpressionSet with copy number data. The featureNames and sampleNames
### must correspond to the eset.expr provided next. Typically created with
### copynumbR.read.segmented(...,gene=TRUE) or
### copynumbR.gistic.read.genes(...).
eset.expr,
### ExpressionSet with expression data. See above, the featureNames must
### correspond to the ones of eset.cn. So for Affymetrix data for example,
### probe sets need to be collapsed (for example with the WGCNA package). 
cutoffs=c(-Inf,-1.3,-0.1,0.1,0.9,Inf),
### Copy number cutoffs.
cutoff.labels=c("Homozyg. Deletion","Heterozyg. Deletion","Normal","Gain","Amplification"),
### The labels of these cutoffs.
probesets=NULL, 
### Display only these genes. If null, show all genes. That default is 
### obviously only useful for already filtered ExpressionSets.
min.samples=3, 
### Minimum number of samples in each cutoff category.
sqrt=FALSE,
### Square root transform the data for read counts?
xlab="Copy Number",
### The label of the x-axis.
ylab="Expression",
### The label of the y-axis.
outlier.shape=NA,
### Display outliers? Passed to geom_boxplot().
plot=TRUE
### Generate a ggplot2 object? If FALSE, then just return a list of means for
### each group and gene.
) {
    
    if (ncol(eset.cn) != ncol(eset.expr)) 
        stop("Number of samples do not match for eset.cn and eset.expr.")
    
    # if sampleNames are different, give only a warning, might be still
    # the correct order.
    if (all(sampleNames(eset.cn) == sampleNames(eset.expr))!=TRUE) 
        warning("sampleNames() of eset.cn and eset.expr do not match.")

    if (!is.null(probesets)) {
        probesets.avail <- probesets[probesets %in%
            intersect(featureNames(eset.expr), featureNames(eset.cn))]
        if (length(probesets) != length(probesets.avail)) 
            warning("Not all probesets available.")
        probesets <- probesets.avail    
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
    res <- lapply(res, function(x) x[,apply(x, 2, sum)>0, drop=FALSE])
    res <- res[lapply(res, function(x) sum(as.vector(x))) > min.samples]

    d.f <- do.call(rbind,unlist(lapply(1:length(res), function(i)
        lapply(colnames(res[[i]]), function(j) try(data.frame(Group=names(res)[i],
        Gene=j
        ,Expr=exprs(eset.expr)[j,res[[i]][,j]]
        ,id=sampleNames(eset.expr)[res[[i]][,j]]
        ,stringsAsFactors=FALSE)))),recursive=FALSE)) 

     d.f$Expr <- as.numeric(d.f$Expr)
     d.f <- d.f[!is.na(d.f$Expr),]
     d.f$Group <- factor(d.f$Group, levels=cutoff.labels[cutoff.labels %in%
        d.f$Group])

    if (plot) {
        p <- ggplot(d.f,
        aes(Group,Expr))+geom_boxplot(outlier.shape = outlier.shape )
        if (sqrt)
            p <- p + scale_y_sqrt(breaks=trans_breaks("sqrt",function(x) x ^
                2)(c(1,1:8*(max(d.f$Expr)/7))))
        p <- p +facet_wrap(~Gene)+theme(axis.text.x=element_text(angle=45,
        hjust=1))+ylab(ylab)+xlab(xlab)
    } else { 
        p=NULL 
    }
    list(plot=p, means= ddply(d.f, c("Group", "Gene"), .fun= function(x)
        mean(x$Expr)))

### A list of a ggplot2 object (plot) and the mean values in each group for each gene
### as data.frame (means)
},ex=function(){
    library(copynumbR)
    library(ggplot2)
    clinical <- read.csv(system.file("extdata", 
        "stransky_bladder_clinical.csv", package="copynumbR"))

    eset.genes <- copynumbR.read.segmented(system.file("extdata",
        "stransky_bladder.glad", package="copynumbR"), clinical, 
        gene=TRUE, geneMap=geneMap_hg17)
    
    # load matched expression data
    data(PMID17099711.GPL91_eset)
    
    # get the samples with matching expression data
    isc <- intersect(sampleNames(eset.genes),
        sampleNames(PMID17099711.GPL91_eset))

    p <- copynumbR.boxplot(eset.genes[,isc], PMID17099711.GPL91_eset[,isc],
        probeset=c("MYC", "ADCY8"))

    # Highlight samples
     plot(p$plot+geom_jitter(aes(shape=eset.genes[,id]$GENDER),size=4)+
        scale_shape_discrete(name="Gender"))
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
### Chromosome frequency plot.
(esets,
### List of ExpressionSets with segmented copy number data.
window=500000,
### Window size in base pairs. Copy numbers will be the median in each window.
labels=paste(names(esets), " (n = ", sapply(esets, ncol),")", sep=""),
### The labels of the ExpressionSets.
gain=0.1,
### Minimum log2 ratio for a copy number gain.
loss=-0.1,
### Maximum log2 ratio for a copy number loss.
sma=10,
### Smooth copy numbers with a simple moving average.
from.chr=1,
### Start plotting at this chromosome.
to.chr=22,
### End plotting at this chromosome.
ylab="Loss / Gain",
### The y-axis label.
xlab="Chromosome",
### The x-axis label.
font.size=12,
### The font size.
centromere.file="hg18"
### File containing the centromere locations for each chromosome.
### These files are already provided for hg17-hg19.
### Select these with centromere.file=c("hg17", "hg18","hg19").
##note<< The centromere files were generated with: \code{curl -s
## "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBand.txt.gz" |
## gunzip -c | grep acen}
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
        labels = .percent_format())+
        theme(text = element_text(size=font.size), panel.grid.minor =
        element_line(size = 0.4, colour =
            'white'),panel.grid.major=element_line(linetype="blank"))
### A ggplot2 object.
},ex=function(){
    library(copynumbR)

    clinical <- read.csv(system.file("extdata", 
        "stransky_bladder_clinical.csv", package="copynumbR"))

    eset <- copynumbR.read.segmented(system.file("extdata", 
        "stransky_bladder.glad", package="copynumbR"), clinical)

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
### Plot a chromosome heatmap.
(eset,
### An ExpressionSet, typically created with copynumbR.read.segmented().
window=1000000,
### Window size in base pairs. Copy numbers will be the median in each window.
gain=0.1,
### Minimum log2 ratio for a copy number gain.
loss=-0.1,
### Maximum log2 ratio for a copy number loss.
Colv=NULL,  
### Hierarchical clustering of class hclust of the columns. 
### If NULL, then calculated as specified by distfun and hclustfun.
### If NA, no ordering of columns is done.
distfun = spearman.dist, 
### The distance function if Colv is NULL.
hclustfun = hclust, 
### The cluster function if Colv is NULL.
from.chr=1,
### Start plotting at this chromosome.
to.chr=22,
### End plotting at this chromosome.
start=NULL, 
### Zoom in, start from the specified coordinates (bp) of from.chr
end=NULL,
### Zoom in, stop at the specified coordinates (bp) of to.chr. Only
### useful when from.chr and to.chr are equal
mart=NULL, 
### A biomaRt object, used for a zoom in to annotate genes.
plot=TRUE, 
### Plot the heatmap, otherwise the dendrogram and heatmap are returned as
### list of ggplot2 objects.
hide.labels=FALSE,
### Hide the sample labels on the bottom.
ylab=("Chromosome"),
### The y-axis label.
centromere.file="hg18"
### File containing the centromere locations for each chromosome.
### These files are already provided for hg17-hg19.
### Select these with centromere.file=c("hg17", "hg18","hg19").
##note<< The centromere files were generated with: \code{curl -s
## "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBand.txt.gz" |
## gunzip -c | grep acen}
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
        genes <- .getGenesRegion(from.chr, start, end, mart)
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
    list(heatmap=p1,dendrogram=p2)
### List of two ggplot2 objects (heatmap and dendrogram).
},ex=function(){
    library(copynumbR)

    clinical <- read.csv(system.file("extdata", 
        "stransky_bladder_clinical.csv", package="copynumbR"))

    eset <- copynumbR.read.segmented(system.file("extdata", 
        "stransky_bladder.glad", package="copynumbR"), clinical)

    p <- copynumbR.heatmap(eset, centromere.file="hg17")
})

copynumbR.read.maf <- function
### Parse a MAF file.
(filename,
### The filename of the MAF file. Can be a vector of filenames.
clinical=NULL,
### A data frame with clinical annotation for the phenoData slot of the output
### ExpressionSet. If NULL, return a data.frame with mutations
coding.fun = function(x) ifelse(x=="Silent",1,2),
### This function return an ExpressionSet with the mutation types coded
### numerically. This function can be used to code mutations. 0 means no
### mutation.
verbose=TRUE,
...
### Additional parameters passed to copynumbR.eset
)
{
    if (verbose) print("Reading MAF file...")
    if (length(filename) > 1) {
        sm <- lapply(filename, read.delim)
        cols <- Reduce(intersect, lapply(sm, colnames))
        sm <- do.call(rbind, lapply(sm, function(x) x[,cols]))
    } else {
        sm <- read.delim(filename)
    }
    sm$.coding <- coding.fun(sm$Variant_Classification)
    if (verbose) print("Generating mutation matrix...")
    mdf <- dcast(Hugo_Symbol ~ Tumor_Sample_Barcode, data = sm, value.var = ".coding", fun.aggregate=max)
    mdf[mdf < 0] <- 0
    if (!is.null(clinical)) {
        eset <- copynumbR.eset(mdf, clinical, ...)
        featureNames(eset) <- mdf[,1]
        return(eset)
    } 
    rownames(mdf) <- mdf[,1]
    return(mdf[,-1])
}

copynumbR.mutation.heatmap <- function
### Plot a heatmap visualizing gene mutations (somatic and copy number)
(eset.cn, 
### ExpressionSet with copy number data, typically gene-level data parsed with
### copynumbR.read.segmented( ... , gene=TRUE)
eset.maf, 
### ExpressionSet with mutation data, read with copynumbR.read.maf()
cutoffs=c(-Inf,-1.3,-0.1,0.1,0.9,Inf),
### Copy number cutoffs.
cutoff.labels=c("Homozyg. Deletion","Heterozyg. Deletion","Normal","Gain","Amplification"),
### The labels of these cutoffs.
mutation.labels=c("Silent","Non-Silent"),
### The labels of the mutations 
probesets
### Plot these probesets (genes) only
) {
    if (!is.null(probesets)) {
        probesets.avail <- probesets[probesets %in%
            intersect(featureNames(eset.maf), featureNames(eset.cn))]
        if (length(probesets) != length(probesets.avail)) 
            warning("Not all probesets available.")
        probesets <- probesets.avail    
        eset.cn <- eset.cn[probesets,]
        eset.maf <- eset.maf[probesets,]
    }

    if (all(featureNames(eset.cn) == featureNames(eset.maf))!=TRUE) {
        stop("featureNames() of eset.cn and eset.maf do not match")
    }
    
    res <- lapply(2:length(cutoffs), function(i)
        apply(exprs(eset.cn),1, function(x) x>cutoffs[i-1]
        & x <=cutoffs[i]))
    names(res) <- cutoff.labels
    res2 <- lapply(1:length(mutation.labels), function(i) apply(exprs(eset.maf),1,
    function(x) x==i))
    names(res2) <- mutation.labels
    res <- c(res2, res)

    res <- lapply(res, function(x) x[apply(x, 1, sum)>0,])
    res <- res[lapply(res, length)>3]

    d.f <- do.call(rbind,lapply(1:length(res),function(i)
        data.frame(melt(res[[i]]), type=names(res)[i])))
    
    colnames(d.f)[1:2] <- c("SampleID", "Gene")

    d.f$Gene <- factor(d.f$Gene, levels=probesets)    

    d.f$alpha= ifelse(d.f$value,1,0)
    d.f[d.f$type=="Normal","alpha"] <- 0
    d.f <- d.f[order(d.f$type),]
    .order <- sapply(levels(d.f[,1]), function(x)
        sum( (length(levels(d.f$type))-as.numeric(d.f[d.f[,1]==x,"type"]))*
             10^(as.numeric(d.f[d.f[,1]==x,"value"]))*
             100^(length(probesets)-as.numeric(d.f[d.f[,1]==x,"Gene"]))  ))
    d.f$.order <- .order[d.f[,1]]

    d.f <- d.f[order(d.f$.order,decreasing=TRUE),]
    d.f$SampleID <- factor(d.f$SampleID, levels=unique(d.f$SampleID))
    d.f <- d.f[order(d.f$alpha,decreasing=TRUE),]
    d.f <- d.f[!duplicated(paste(d.f[,1], d.f[,2])),]

    p <- ggplot(d.f, aes(SampleID,
    Gene,alpha=alpha,fill=type))+geom_tile()+ylab("")+theme_grey(16)+theme(axis.text.x=element_blank(),
    axis.ticks.x=element_blank())+xlab("")+scale_alpha_continuous(guide=FALSE)+scale_fill_discrete(name
    = "Mutation Type")+guides(fill = guide_legend(override.aes= list(alpha =
     ifelse(levels(d.f$type)=="Normal", 0, 1))))


### A ggplot2 object.    
}

copynumbR.mutation.table <- function
### Extract muations from copynumbR.mutation.heatmap
(p 
### The ggplot2 object returned by copynumbR.mutation.heatmap()
) {
    X <- p$data
    # code normal as no mutation, obviously
    X[X$type=="Normal","value"] <- FALSE
    X <- dcast(SampleID~Gene+value,data=X)
    # not present means no mutation
    X[is.na(X)] <- FALSE
    X <-X[,c(1, grep("TRUE",colnames(X)))]
    rownames(X) <- X[,1]
    X[,-1]
### A data.frame with mutations in columns and samples in rows
}

copynumbR.absolute.run <- function
(filename,
sigma.p = 0,
max.sigma.h = 0.02,
min.ploidy = 0.95,
max.ploidy = 10,
max.as.seg.count = 1500,
max.non.clonal = 0,
max.neg.genome = 0,
copy_num_type="total",
results.dir=tempdir(),
verbose=FALSE,
trans.fun=function(x) x,
...
) 
{
    x <- read.delim(filename)
    x <- lapply(levels(x[,1]), function(s) x[x[,1]==s,])
    
    
    .runAbsolute <- function(s, ...) {
        fn <- tempfile(pattern="copynumbR")
        sample.name=s[1,1]
        s <- s[,-1]
        colnames(s)[1:5] <- c("Chromosome", "Start", "End", "Num_Probes",
        "Segment_Mean")
        
        s[,5] <- trans.fun(s[,5])
        s <- s[s[,1] < 23,]

        if (verbose) cat("Writing ABSOLUTE input segmented file in", fn, 
        "\n and output files in", results.dir, "\n")

        write.table(s, file=fn, row.names=FALSE, quote=FALSE, sep="\t")
        RunAbsolute(fn, sigma.p=sigma.p, max.sigma.h=max.sigma.h,
        min.ploidy=min.ploidy, max.ploidy=max.ploidy,
        sample.name=sample.name, max.as.seg.count=max.as.seg.count,
        copy_num_type=copy_num_type, results.dir=results.dir,
        max.non.clonal=max.non.clonal, max.neg.genome=max.neg.genome,
        verbose=verbose, ...)
        load(paste(results.dir, "/", sample.name, ".ABSOLUTE.RData", sep=""))
        seg.dat
    }

    lapply(x, .runAbsolute, ...)

}


theme_classic2 <- structure(function
### A simple ggplot2 theme, no fancy stuff.
(base_size = 12, 
### Theme base font size.
base_family = ""
### Theme base font family.
) 
{
    theme_grey(base_size = base_size, base_family = base_family) %+replace% 
        theme(axis.text = element_text(size = rel(0.8)), 
              axis.ticks = element_line(colour = "black"), 
              legend.key = element_blank(),
              panel.background = element_rect(fill = "white", colour = NA), 
              panel.border = element_rect(fill = NA, colour = "black", size=1.0), 
              panel.grid.major = element_line(linetype="blank", size = 0.2), 
              panel.grid.minor = element_line(linetype = "blank", size = 0.5), 
              panel.margin     = unit(0.75, "lines"),
              strip.background = element_rect(fill = NA, colour = "black",linetype="blank"), 
              strip.text.x = element_text(face="bold"))
},ex=function(){
    library(copynumbR)

    clinical <- read.csv(system.file("extdata", 
        "stransky_bladder_clinical.csv", package="copynumbR"))

    eset <- copynumbR.read.segmented(system.file("extdata", 
        "stransky_bladder.glad", package="copynumbR"), clinical)

    # find the non-muscle-invasive samples
    idx.noninvasive <- grepl("Ta|T1", eset$T)
    
    # create a list of ExpressionSets with the low- and high stage samples
    eset.stage <- list("Non-Muscle-Invasive"=eset[,idx.noninvasive],
        "Invasive"=eset[,!idx.noninvasive])
    
    # now compare the copy numbers of these two groups
    p <- copynumbR.plot(eset.stage, centromere.file="hg17", sma=0)

    plot(p+theme_classic2(14)) 
})


.getGenesRegion <- function(chromosome, start, end, mart=NULL) {
    if (is.null(mart)) {
       mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    }
    x <-getBM(c("hgnc_symbol", "chromosome_name", "transcript_start",
        "transcript_end"), values=chromosome, filters="chromosome_name", mart)
    x <- x[x$transcript_start > start & x$transcript_end < end &
    x$hgnc_symbol != "",]
    x[!duplicated(x$hgnc_symbol),] 
}

# taken from the scales package
.percent_format <- function() {
    precision <- function(x) {
        rng <- range(x, na.rm = TRUE)
        
        span <- if (zero_range(rng)) rng[1] else diff(rng)
        10 ^ floor(log10(span))
    }
    function(x) {
        x <- abs(round_any(x, precision(x) / 100))
        str_c(comma(x * 100), "%")    
    }
}

setMethod("coerce", 
### Coerce ExpressionSet from the copynumbR package into a CNA object as
### defined in the DNAcopy package
signature(from="ExpressionSet", ##<< From ExpressionSet
          to="CNA" ##<< To CNA
), function(from, to) { 
    CNA(exprs(from), featureData(from)$chr,
    featureData(from)$start, sampleid=sampleNames(from))
})

copynumbR.cor.genes.test <- structure(function
### Find correlated probesets.
(probeset, 
### Probeset to test
eset.expr,
### ExpressionSet containing the expression data of probeset. 
eset.expr2=eset.expr, 
### ExpressionSet which we test for correlated expression, i.e, eset.expr[probeset,] will be tested for correlated
### expression with eset.expr2. Default is test for correlation in the same
### data as probeset. 
method="BH", 
### p.adjust method for adjusting for multiple testing.
cutoff=0.05, 
### P-value cutoff.
n=20, 
### Maximum number of reported correlated probesets.
direction="positive",
### Direction. Either "positive" or "negative".
annotation=NULL
### Translate probeset ids to symbols with the getSYMBOL(probeset, annotation) function if not NULL.
) {
    res <- apply(exprs(eset.expr2),1, function(x)
        cor.test(exprs(eset.expr)[probeset,],x))
    pval <- p.adjust(sapply(res, function(x) x$p.value), method=method)
    rho <- sapply(res, function(x) x$estimate)
    if (direction=="positive") {
        n <- min(n, sum(rho > 0 & pval < cutoff))
        probesets <- featureNames(eset.expr2)[order(rho, decreasing=TRUE)]
    } else {
        n <- min(n, sum(rho < 0 & pval < cutoff))
        probesets <- featureNames(eset.expr2)[order(rho, decreasing=FALSE)]
    }

    if (!is.null(annotation)) { 
        probesets <- na.omit(getSYMBOL(probesets, annotation))
    }    
    return(head(probesets,n))
},ex=function(){
    library(copynumbR)
    data(PMID17099711.GPL91_eset)

    res <- copynumbR.cor.genes.test("MYC", PMID17099711.GPL91_eset,
    method="none")

    cor.test(exprs(PMID17099711.GPL91_eset)["MYC",],
        exprs(PMID17099711.GPL91_eset)["ELL2P1///ELL2",])

})

