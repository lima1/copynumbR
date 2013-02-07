copynumbR.gistic.read.lesions <- function
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
}

attr(copynumbR.gistic.read.lesions,"ex") <- function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <-
    copynumbR.gistic.read.lesions(system.file("extdata/gistic_stransky_bladder",
        "all_lesions.conf_95.txt", package="copynumbR"), clinical)

    # list all recurrent gains and losses
    featureData(eset)$Descriptor
    
    # show the copy number distributions of all recurrent alterations
    boxplot(t(exprs(eset)))
}

copynumbR.gistic.read.genes <- function
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
}

attr(copynumbR.gistic.read.genes,"ex") <- function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <-
    copynumbR.gistic.read.genes(system.file("extdata/gistic_stransky_bladder",
        "all_data_by_genes.txt", package="copynumbR"), clinical)

    # show the copy number distribution of MYC copy numbers
    boxplot(exprs(eset)["MYC",])
}

copynumbR.gistic.write.arrayfile <- function
### Write GISTIC input arrayfile 
(labels, 
### The sample names as used in the segmented file
file
### The output filename of this function
) {
    write.table(data.frame(Array=labels),file=file, row.names=FALSE,quote=FALSE)
}


copynumbR.gistic.read.broad <- function
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
}

attr(copynumbR.gistic.read.broad,"ex") <- function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <-
    copynumbR.gistic.read.broad(system.file("extdata/gistic_stransky_bladder",
        "broad_values_by_arm.txt", package="copynumbR"), clinical)
    
    # display the distributions of copy numbers of all chromosome arms
    boxplot(t(exprs(eset)))
}

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

copynumbR.gistic.genes.band <- function
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
}

attr(copynumbR.gistic.genes.band,"ex") <- function(){
    library(copynumbR)
    # extract the gene symbols in the GISTIC peak 8q24.21
    band <-
    copynumbR.gistic.genes.band(system.file("extdata/gistic_stransky_bladder",
        "amp_genes.conf_95.txt", package="copynumbR"), "8q24.21")
}

copynumbR.gistic.focal <- function
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

}

attr(copynumbR.gistic.focal,"ex") <- function(){
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
}


copynumbR.gistic.armplot <- function
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
    require(ggplot2)
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
        gp =
            ggplot(data.stack, aes(p,
            q))+geom_point(alpha=0.2)+facet_wrap(~Chromosome, ncol = ncol)

    }
    list(plot=gp, data=data.stack)
    ### A list containing the ggplot2 and the data as used in ggplot2
}

attr(copynumbR.gistic.armplot,"ex") <- function(){
    library(copynumbR)

    res <-
        copynumbR.gistic.armplot(file=
            system.file("extdata/gistic_stransky_bladder", 
                "broad_values_by_arm.txt", package="copynumbR")
        )

    plot(res[[1]])
}


copynumbR.tcga.input <- function(path=".", output="tcga.txt", verbose=TRUE,
grepregex="\\.data.txt", row.names=1, var.column=1) {
    files = dir(paste(path, sep="/"), full.names=TRUE)
    files = files[grep(grepregex, files)]
    if (verbose) cat("Reading", files,sep="\n") 
    data <- lapply(files, read.delim, stringsAsFactors=FALSE,
        row.names=row.names, as.is=TRUE)
    cdata <- do.call(cbind, lapply(data, function(X) X[,var.column]))
    colnames(cdata) <-  gsub("^.*TCGA-","TCGA-", files)
    rownames(cdata) <-  rownames(data[[1]])
    write.table(cdata, file=output, quote=FALSE, sep="\t") 
    cdata
}

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

copynumbR.read.segmented <- function
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
    require(CNTools)
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
}

attr(copynumbR.read.segmented,"ex") <- function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <- copynumbR.read.segmented(system.file("extdata", "stransky_bladder.glad", package="copynumbR"), clinical)
    # Plot the distribution of copy numbers of segments in all samples
    plot(exprs(eset))

    eset.genes <- copynumbR.read.segmented(system.file("extdata",
     "stransky_bladder.glad", package="copynumbR"), clinical, gene=TRUE)

    boxplot(exprs(eset.genes)["MYC",])
}

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

copynumbR.clone.gistic.eset <- function
### Extract matching copy numbers 
(eset.gistic,
### An ExpressionSet with GISTIC peaks read with copynumbR.gistic.read.lesions
eset.segmented)
### An ExpressionSet with segmented data, typically from another cohort,
### read with copynumbR.segmented. So this is useful to compare copy numbers at
### GISTIC peaks across cohorts
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
}

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

copynumbR.eset <- function
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
os.col = NULL,
os.event.col=NULL,dss.col= NULL, dss.event.col=NULL) {
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

    if (!is.null(os.col)) { 
        eset$os.surv = Surv(as.numeric(eset[[os.col]]),
            eset[[os.event.col]])
        eset = .attachWaldPV(eset)    
    }        
    if (!is.null(dss.col)) {
        eset$dss.surv = Surv(eset[[dss.col]], eset[[dss.event.col]])
        eset = .attachWaldPV(eset, col="dss.surv",field="dss")    
    }    
    eset     
# An ExpressionSet object    
}    

.attachWaldPV <- function(eset, col="os.surv", field="os") {
        ft = lapply(1:nrow(eset), function(i) coxph(formula = eset[[col]] ~
            exprs(eset)[i,], data = eset))
        pv = sapply(ft, function(x) 1 - pchisq(as.vector(x$wald.test),1))
        hr = sapply(ft, function(x) x$coef)
        featureData(eset)[[field]] = data.frame(pv = pv, hr = hr)
        eset
}

copynumbR.boxplot <- function
### A boxplot showing correlation of copy number and Expression for matched
### data
(eset.cn, 
### ExpressionSet with copy number data
eset.expr,
### ExpressionSet with expression data
cutoffs=c(-Inf,-1.3,-0.1,0.1,0.9,Inf),
### Copy number cutoffs 
cutoff.labels=c("High Loss","Loss","Normal","Gain","Amplification"),
### The labels of these cutoffs
probesets=NULL, 
### Display only these genes. If null, show all genes.
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
ylab="Expression"
### The label of the y-axis
) {

    # if user utilized our copynumbR, symbols might be in the first
    # featureData slot
    if (length( intersect(featureData(eset.cn)[[1]], featureNames(eset.expr)))  > 
        length( intersect(featureNames(eset.cn), featureNames(eset.expr)) )) {
            eset.cn <- eset.cn[match(featureNames(eset.expr),
            featureData(eset.cn)[[1]]),]
            featureNames(eset.cn) <- featureNames(eset.expr)
    }    

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
    aes(Group,Expr))+geom_boxplot(outlier.shape = NA )
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
}


copynumbR.boxplot.single <- function(cn, expr,
cutoffs=c(-Inf,-1.3,-0.2,0.2,0.9,Inf),
cutoff.labels=c("High Loss","Loss",
"Normal","Gain","Amplification"),min.samples=3, gene.name,
highlight=NULL,highlight.labels=NULL) {
 res <- lapply(2:length(cutoffs), function(i)
    cn >cutoffs[i-1]
    & cn <=cutoffs[i])
names(res) <- cutoff.labels    
res <- res[lapply(res, function(x) sum(as.vector(x))) > min.samples]

   d.f <- do.call(rbind, lapply(1:length(res), function(i)
    data.frame(Group=names(res)[i],
    Gene=gene.name,Expr=expr[res[[i]]],Highlight=highlight[res[[i]]] )))
 p <- ggplot(d.f,
 aes(Group,Expr))+geom_boxplot(outlier.shape = NA )+scale_y_sqrt(breaks=trans_breaks("sqrt",
 function(x) x ^
 2)(c(1,1:8*400)))+theme(axis.text.x=element_text(angle=45,
 hjust=1))+ylab("Nanostring Read Count")+
 xlab("Copy Number")+scale_colour_discrete(name =
 "FGFR3")+scale_shape_discrete(name="FGFR3",
 solid=TRUE)+scale_size_discrete(name="FGFR3",range=c(2,5))
 plot(p+geom_point(aes(shape=Highlight,size=Highlight))+theme_grey(base_size=18))
}


.getCentromer <- function(eset, centromer.file) {
    centromer <- read.delim(centromer.file, as.is=TRUE, header=FALSE)
    centromer[,1] <- gsub("chr","",centromer[,1])
    centromer <- centromer[!duplicated(centromer[,3]),]
    centromer[,1] <- gsub("X","23",centromer[,1])
    centromer[,1] <- gsub("Y","24",centromer[,1])
    centromer[,1] <- as.numeric(centromer[,1])
    centromer <- centromer[order(centromer[,1]),]
        idx.min <- sapply(1:22,function(i) min(which(centromer[,1]==i)))
        idx.max <- sapply(1:22,function(i) max(which(centromer[,1]==i)))

    centromer <- cbind(centromer[idx.min,2], centromer[idx.max,3],0)
    centromer[13:15,1:2] <- 0

    chrs <- unique(featureData(eset)[[1]])
    chrl <- sapply(chrs, function(x)
    max(featureData(eset[featureData(eset)[[1]]==x,])[[3]]))
    chrstart <- sapply(1:22, function(i) min(featureData(eset[featureData(eset)[[1]] == i,])[[2]]))
    chrl <- c(0,chrl - chrstart)


    centromer <- cbind(cumsum(chrl)[1:22]+centromer[,1],
                    cumsum(chrl)[1:22]+centromer[,2], 0)
    list(centromer=centromer, chrstart=chrstart, chrs=chrs, chrl=chrl)                
}

.getSmoothedData <- function(esets, window, labels, gain, loss, sma,
    centromer.file) {
    require(TTR)
    require(plyr)
    hg18 <- .getCentromer(esets[[1]], centromer.file)

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
        df$field <- SMA(df$field,sma+1)
        df$field[1:sma] <- x[1:sma]
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


copynumbR.plot <- function
### Chromosome frequency plot
(esets,
### List of ExpressionSets
window=500000,
### Window size in base pairs. Copy numbers will be averaged in each window
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
centromer.file=system.file("extdata", "hg18centromer.txt",
package="copynumbR")
### File containing the centromer locations for each chromosome
) {
    require(ggplot2)
    
    res <- .getSmoothedData(esets,window,labels,gain,loss,sma,centromer.file)
    
    df <- res$df

    q <- ggplot(df, aes(Begin,field,ymin=0,
        ymax=field,fill=type))+geom_ribbon()
    q <- q+ scale_fill_brewer("",palette="Set1")+labs(fill="") +
        scale_x_continuous(minor_breaks=cumsum(res$hg$chrl)/window,
        breaks=round(res$hg$centromer[,1]/window)
        ,labels=res$hg$chrs[1:22],limits=c(from.chr-1,cumsum(res$hg$chrl)[to.chr+1]/window))+
        ylab("Loss (%) / Gain (%)")+xlab("Chromosome")+facet_grid(label ~ .)
    q <- q +  scale_y_continuous(minor_breaks = NULL,
        breaks=c(-1,-0.5,0,0.5,1),labels=c(100,50,0,50,100))+
        theme(panel.grid.minor =
        element_line(size = 0.4, colour =
            'white'),panel.grid.major=element_line(linetype="blank"))
### A ggplot2 object
}

attr(copynumbR.plot,"ex") <- function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <- copynumbR.read.segmented(system.file("extdata", "stransky_bladder.glad", package="copynumbR"), clinical)
    p <- copynumbR.plot(list(eset), labels="Stransky")
    plot(p) 
}

.fillStacked <- function(X) {
    Y <- do.call(rbind, lapply(2:nrow(X), function(i) if (X[i-1,2] ==
    X[i,2] && X$alpha[i] != 0)
    data.frame(Begin=X[i-1,1]:X[i,1],values=X[i,2], ind=X[i,3],
    values_orig=X[i,4],alpha=X[i,5])
    else X[i-1,]))
    Y[!duplicated(paste(Y$Begin, Y$ind)),]
}

copynumbR.heatmap <- function
### Chromosome Heatmap
(eset,
### ExpressionSet, typically created with copynumbR.read.segmented
window=1000000,
### Window size in base pairs. Copy numbers will be averaged in each window
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
centromer.file=system.file("extdata", "hg18centromer.txt",
package="copynumbR")
### File containing the centromer locations for each chromosome
) {
    require(TTR)
    require(RColorBrewer)
    require(ggdendro)
    require(scales)
    require(bioDist)
    require(ggplot2)
    require(plyr)
    require(grid)
    sampleNames(eset) <- make.names(sampleNames(eset))
    eset <- eset[featureData(eset)[[1]]>= from.chr & featureData(eset)[[1]] <=
        to.chr] 

    if (is.null(Colv)) {
        Colv <- hclustfun(distfun(eset))
    }

    eset <- eset[,Colv$order]

    if (!is.null(start)) {
        idx <- (featureData(eset)[[1]]== from.chr &
            featureData(eset)[[2]] >= start) & 
               (featureData(eset)[[1]]== to.chr &
            featureData(eset)[[2]] <= end)
        eset <- eset[ , order(apply(exprs(eset)[idx,],2,mean),decreasing=TRUE)]
    } 

    hg18 <- .getCentromer(eset, centromer.file)

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
        breaks=round(hg18$centromer[,1]/window*s)
        ,labels=hg18$chrs[1:22], trans=reverse_trans(),
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
}

attr(copynumbR.heatmap,"ex") <- function(){
    library(copynumbR)
    clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
    eset <- copynumbR.read.segmented(system.file("extdata", "stransky_bladder.glad", package="copynumbR"), clinical)
    p <- copynumbR.heatmap(eset)
}


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
        "transcript_end"), values=chromosome, filter="chromosome_name", ensembl)
    x <- x[x$transcript_start > start & x$transcript_end < end &
    x$hgnc_symbol != "",]
    x[!duplicated(x$hgnc_symbol),] 
}

