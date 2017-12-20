#== Args 
args<-commandArgs(TRUE)
#===mandatory parameters
phenotype.file <- args[1]
outcome.name <- args[2]
outcome.type <-  args[3]
covariate.string <- args[4]
snpinfo.file <- args[5]
genotype.files <- args[6]
output.file <- args[7]

#==optional parameters
kinship.matrix <- args[8]
pheno.id <- args[9]


# added these to JSON
BUFFER <- as.numeric(args[10]) 
gene.file <- args[11] 
agg.file <- args[12] 
snp.filter <- args[13] 
gene.filter <- args[14]
top.maf <- as.numeric(args[15]) 
test.stat <-  args[16] # Score, Wald, Firth
test.type  <-  args[17] # Burden, Single, SKAT
min.mac <- as.integer(args[18])
weights <- args[19]
conditional <- args[20]
user_cores <-  as.numeric(args[21])
het_vars <-  args[22]


weights = eval(parse(text=weights))
cat('Weights',weights,class(weights),'\n')

USE_AGG = F
# GLOBAL VARIABLES
collapsing.tests <- c("SKAT",  "Burden")
test.type.vals <- c("Single","SKAT", "Burden")
test.stat.vals <- c("Score", "Wald", "Firth")

if(!test.type %in% test.type.vals){
    stop("Test type must be one of ",paste(test.type.vals,sep=','))
}

aggregateListByAllele <- function(gds, variants, indexOnly=FALSE) {
    stopifnot(all(c("group_id", "chromosome", "position", "ref", "alt") %in% names(variants)))

    ## set filter to listed variants only
    filtOrig <- seqGetFilter(gds)
    gr <- GRanges(seqnames=variants$chromosome,
                  ranges=IRanges(variants$position, variants$position))
    seqSetFilter(gds, gr, verbose=FALSE)
  
    variants <- .expandAlleles(gds) %>%
        inner_join(variants, by=c("chromosome", "position", "ref", allele="alt"))

    seqSetFilter(gds, sample.sel=filtOrig$sample.sel,
                 variant.sel=filtOrig$variant.sel, verbose=FALSE)

    .groupVariants(variants, indexOnly)
}


GetFamilyDistribution <- function(response.type) {
               if (response.type == "Continuous"){
                      family = "gaussian"
               } else if (response.type == "Dichotomous"){
                      family = "binomial"
               } else {
                      msg = paste("Don't know how to deal with response type", response.type)
                      stop(msg)
               }
               return(family)
           }

GetKinshipMatrix <- function(kinship.matrix){
  cat('Loading Kinship Matrix:',kinship.matrix,'\n')
  if(grepl('Rda',kinship.matrix,ignore.case=TRUE)){
    kmatr = get(load(kinship.matrix))
  }
  else{
    kmatr = as.matrix(read.csv(kinship.matrix,as.is=T,check.names=F,row.names=1))
  }

  cat('Loaded Kinship NROW:',NROW(kmatr),' NCOL:',NCOL(kmatr),'\n')
  kmatr
}

getMAC <- function(gds.file){
  geno.dat <- seqGetData(gds.file,"genotype")
  apply(geno.dat,3,function(x){min(sum(x == 0,na.rm=T),sum(x == 1,na.rm=T))})
}

cat('output.file',output.file,'\n')
cat('kinship.matrix',kinship.matrix,'\n')
cat('buffer',BUFFER,'\n')
cat('genefile',gene.file,'\n')
cat('varaggfile',agg.file,'\n')
cat('snp.filter',snp.filter,'\n')
cat('top.maf',top.maf,'\n')
cat('test.stat',test.stat,'\n')
cat('test.type',test.type,'\n')
cat('outcome.type',outcome.type,'\n')
cat('het_vars',het_vars,'\n')
cat('conditional',conditional,'\n')

if(conditional != 'NA'){
  cpos = strsplit(conditional,':')[[1]][2]
}else{
  cpos = FALSE
}
cat('conditional',conditional,'\t pos',cpos,'\n')

if (!(test.stat %in% test.stat.vals)){
     msg = paste("The requested test statistic:", test.stat, "is not available (Use Firth, Score, Wald!")
     stop(msg)
}


suppressMessages(library(SeqArray))
suppressMessages(library(SeqVarTools))
suppressMessages(library(GWASTools))
suppressMessages(library(gap))
suppressMessages(library(Matrix))
suppressMessages(library(plyr))
suppressMessages(library(gdsfmt))
suppressMessages(library(bdsmatrix))
suppressMessages(library(parallel))
suppressMessages(library(GENESIS))
suppressMessages(library(data.table))
suppressMessages(library(doMC))

num_cores <- detectCores(logical=TRUE)


registerDoMC(cores=min(c(user_cores,num_cores)))
cat('Running Analysis with ',min(c(as.numeric(user_cores),num_cores)),' cores of ',num_cores,'\n')

## Setup
source("/home/dnanexus/pipelineFunctions.R")
covariates <- split.by.comma(covariate.string)  


## snp info
cat('Reading snpinfo....')
snpinfo <- fread(snpinfo.file,data.table=F)
snpinfo = as.data.frame(snpinfo)[,c('CHR','POS')]
if(!all(c("CHR","POS") %in% names(snpinfo))){
    
     msg = paste("SNPinfo file must contain column names 'CHR' and 'POS'")
     stop(msg)
}
cat('done\n')

cat('Input SNPINFO N=',nrow(snpinfo),'\n')
snpinfo = eval(parse(text= paste0('subset(snpinfo,',snp.filter,')')))
cat('Output SNPINFO N=',nrow(snpinfo),'\n')


## phenotype 
phenotype.data <- read.csv(phenotype.file, header=TRUE, as.is=TRUE)
if(NCOL(phenotype.data) < 2){
    
     msg = paste("Is the phenotype file a CSV?  Too few columns from read.csv()")
     warning(msg)
}


cat('Input pheno N=',nrow(phenotype.data),'\n')
if(het_vars != 'NA'){
    cat('prep pheno with het vars')
    pheno <- reducePheno(phenotype.data, outcome.name, covariates=covariates,hetvars=het_vars, id=pheno.id)
}else{
    cat('prep pheno without het vars\n')
    het_vars = NA
    pheno <- reducePheno(phenotype.data, outcome.name, covariates=covariates, id=pheno.id)
}
cat('Output pheno N=',nrow(pheno),'\n')

## Report dropped individuals
dropped.ids.selector <- !(phenotype.data[[pheno.id]] %in% row.names(pheno))
dropped.ids <- phenotype.data[[pheno.id]][dropped.ids.selector] 
if (NROW(dropped.ids) != 0 ) {
  cat("Dropped because of incomplete cases:", length(dropped.ids) )
}

# For GDS files
f <- seqOpen(genotype.files)
sample.ids <- seqGetData(f, "sample.id")
all.terms <- unique(c(outcome.name, covariates, het_vars))
pheno <- pheno[row.names(pheno) %in% sample.ids,na.omit(all.terms),drop=F]
cat('Output pheno after mergeing with Genos N=',nrow(pheno),'\n')
if(nrow(pheno) == 0){
    msg = paste("Phenotype ID column doesn't match IDs in GDS")
    stop(msg)
}

full.sample.ids <- sample.ids 

#subset to phenotyped samples
seqSetFilter(f,sample.id = row.names(pheno))

# order pheno to the GDS subject order
sample.ids <- seqGetData(f, "sample.id")
pheno <- pheno[match(sample.ids,row.names(pheno)),,drop=F]

# get position list before any variant filters are applied
pos = seqGetData(f, "position")
variant.ids = seqGetData(f, "variant.id")




## Conditional analaysis
if(cpos >0){
  cat('Conditioning on ',conditional,'...\n')
  cidx = which(pos == as.numeric(cpos))
  if(any(cidx)){
    
    seqSetFilter(f,variant.sel=cidx, sample.id = row.names(pheno),verbose=FALSE)
    pheno$csnp = altDosage(f,  use.names=FALSE)
  }else{
    stop('Can not find snp ',conditional,' with position ',cpos,' to condition on in data file')
  }
  
  dropConditionalCases = NROW(pheno)-NROW(pheno[complete.cases(pheno),])
  if(dropConditionalCases > 0){
    cat('Warning: Dropping ',dropConditionalCases,' samples due to missing conditional genotype calls\n')
  }

  pheno = pheno[complete.cases(pheno),]
  
  covariates[length(covariates) + 1] <- 'csnp'
}



if(agg.file != 'NO_VAR_AGG_FILE'){
    USE_AGG = T
    system.time({ agg = fread(agg.file, stringsAsFactors=F, sep=',', header=T,data.table=F) })
    if('chromosome' %in% names(agg)){
        agg$chr = sub('chr','',agg$chromosome)
    }else if('chr' %in% names(agg)){
        agg$chr = sub('chr','',agg$chr)
    }
    agg$var_long.id = paste(agg$chr,agg$pos,agg$ref,agg$alt,sep=':')

    
    pos = seqGetData(f, "position")
    
    variant.ids = seqGetData(f, "variant.id")
    chr = seqGetData(f, "chromosome")
    allele.ref = sapply(seqGetData(f, "allele"),FUN=function(x)strsplit(x,',')[[1]][1])
    allele.alt = sapply(seqGetData(f, "allele"),FUN=function(x)strsplit(x,',')[[1]][2])
    var_long.ids = paste(chr,pos,allele.ref,allele.alt,sep=':')
    genes <- unique(agg$group_id)
    cat('N AggUnits=',length(genes),'done\n')
}else{


##  Gene list
#
# Aggregation file
#  must contain 'start','stop' and 'CHR'
#  start and stop must be numeric
#
cat('reading gene file...')
if(gene.file == "NO_GENE_REGION_FILE" & test.type != 'Single'){
  stop('For aggregate tests you must provide a aggregation file.')
  
}else if(gene.file == "NO_GENE_REGION_FILE" & test.type == 'Single'){
  
  batchsize = 200000
  if(gene.file == "NO_GENE_REGION_FILE"){
    nbatch = ceiling(max(pos)/batchsize)
    kg = data.frame('name'=paste0('batch',1:nbatch),'start'=batchsize*((1:nbatch)-1),'stop'=batchsize*1:nbatch,stringsAsFactors=F)
  }
  
    genes <- kg$name
    
}else if(variantLists == 1){
  system.time({ kg = fread(gene.file, stringsAsFactors=F, sep=',', header=T) })
  kg
    
}else{
  system.time({ kg = fread(gene.file, stringsAsFactors=F, sep=',', header=T) })

  if(!(sum(names(kg) %in% c('start','stop','chr')) == 3 & is.numeric(kg$start) & is.numeric(kg$stop))){
    stop("Aggreagation file must contain columns 'start','stop' and 'chr'.  Columns start and stop must be numeric")
  }

  if(!( is.numeric(kg$start) & is.numeric(kg$stop))){
    stop("Columns start and stop in aggregation file must be numeric")
  }

  if( ! all(grepl('chr',kg$chr))){
    stop("chr column in aggregation file must be formated as 'chr#' (e.g. chr1) ")
  }

  cat('Aggregate file filter',gene.filter,'\n')
  kg = eval(parse(text= paste0('subset(kg,',gene.filter,')')))
  genes <- kg$name
}
cat('NGENEs=',NROW(kg),'done\n')
}






## Load KINSHIP matrix
## Kinship doesn't contain all samples
kmatr = GetKinshipMatrix(kinship.matrix)
pheno = pheno[row.names(pheno) %in% row.names(kmatr),,drop=F]
kmatr = kmatr[row.names(kmatr) %in% row.names(pheno),colnames(kmatr) %in% row.names(pheno)]
cat('Output pheno in Kinship N=',nrow(pheno),'\n')
kmatr = kmatr[match(row.names(pheno),row.names(kmatr)),match(row.names(pheno),colnames(kmatr))]
if(nrow(pheno) == 0){
    msg = paste("Phenotype ID column doesn't match IDs in Kinship Matrix")
    stop(msg)
}


# Get sample ids to check order 
seqSetFilter(f,sample.id = row.names(pheno))
sample.ids <- seqGetData(f, "sample.id")

if (!(identical(sample.ids,row.names(pheno)) && identical(row.names(kmatr),row.names(pheno)))){
        stop("Something is off problem with re-ordering")
}

seqClose(f)

sample.data <- data.frame(scanID = row.names(pheno),  
                    pheno, 
                    stringsAsFactors=F)
scan.annotated.frame <- ScanAnnotationDataFrame(sample.data)
modified.pheno = pheno[full.sample.ids,,drop=FALSE]
row.names(modified.pheno) <- full.sample.ids

sample.data.for.annotated <- data.frame(sample.id = full.sample.ids,
                                        modified.pheno,
                                        stringsAsFactors=F)
rm(modified.pheno)

annotated.frame <- AnnotatedDataFrame(sample.data.for.annotated)

###################
## NULL MODEL
##################
# Should depend on response type

cat('start fit....\n')
kmatr = as.matrix(kmatr)
if (test.stat == 'Firth'){
  cat('WARNING: Firth test does NOT use kinship information - unrelated only')
  nullmod <- fitNullReg(scanData = scan.annotated.frame,
                     outcome = outcome.name,
                     covars = covariates,
                     family = GetFamilyDistribution(outcome.type))

}else if (!is.na(het_vars)){
  cat('Fitting model with heterogeneous variances')
  nullmod <- fitNullMM(scanData = scan.annotated.frame,
                     outcome = outcome.name,
                     group.var = het_vars,
                     covars = covariates,
                     family = GetFamilyDistribution(outcome.type),
                     covMatList = kmatr)
}else{
  cat('Fitting model ')
  nullmod <- fitNullMM(scanData = scan.annotated.frame,
                     outcome = outcome.name,
                     covars = covariates,
                     family = GetFamilyDistribution(outcome.type),
                     covMatList = kmatr)
}

#sum( sapply(ls(),function(x){object.size(get(x))})) 
#sort( sapply(ls(),function(x){object.size(get(x))})) 
rm(kmatr)
rm(scan.annotated.frame)
sample_ids = row.names(pheno)
rm(pheno)
cat('SNPinfo format\n')
head(snpinfo,2)
snpinfo = as.data.frame(snpinfo)[,c('CHR','POS')]
x = gc()


mcoptions <- list(preschedule=FALSE, set.seed=FALSE)
sm_obj <- 
foreach (current.gene=genes, 
         .combine=function(...){rbindlist(list(...),fill=T)},
         .multicombine = TRUE,
         .inorder=FALSE,  
         .options.multicore=mcoptions) %dopar% {

  ##############
  ## Apply variant filters to get a list of variants in each gene
  ###############
  if(USE_AGG){
      geneSNPinfo = agg[agg$group_id == gene,]
      snp_idx = variant.ids[var_long.ids %in% geneSNPinfo$var_long.id]

  }else{
      gidx <- which(kg$name == current.gene)
      if (length(gidx) != 1) {
          msg=paste("The assumptions that the list of genes is unique is violated\n",
               "The gene:", cgene, "repeats", length(gidx), "times")
          stop(msg)
      } 
      geneSNPinfo = subset(snpinfo, (POS > (kg[gidx,]$start - BUFFER) & POS <= (kg[gidx,]$stop + BUFFER)))
      snp_idx = variant.ids[which(pos %in%  geneSNPinfo$POS)]
      cat(current.gene,'\t',' NSNPS:',length(snp_idx),': range(',(kg[gidx,]$start - BUFFER),'-',(kg[gidx,]$stop + BUFFER),')\n')
  }
  

    
    ## extract genotypes
    f <- seqOpen(genotype.files)
    
    seqSetFilter(f,variant.id = snp_idx,sample.id = sample_ids,verbose=FALSE)
    if(length(seqGetData(f,'variant.id')) > 0){
    
        ## filter to maf and remove monomorphic
        maf <- seqAlleleFreq(f)
        maf <- ifelse(maf < 0.5, maf, 1-maf)
        filtered.alleles = TRUE
        if (test.type %in% collapsing.tests){
            filtered.alleles <- maf < top.maf
        }
        mac <- getMAC(f)
        gc() # whole genotype matrix hanging around in memory
        if (test.type ==  'Single') filtered.alleles <- mac  >= min.mac
      
    
        seqSetFilter(f,variant.id=snp_idx[filtered.alleles & maf > 0], verbose=FALSE)
        num.polymorphic.snps <- seqSummary(f, "genotype", check="none", verbose=FALSE)[["seldim"]][3]
        num.snps = length(maf)
    
        if(num.polymorphic.snps > 0){
      
            genotype.data <- SeqVarData(f, sampleData=annotated.frame)
            ## Collapse test
            if (test.type %in% collapsing.tests) {
                xlist = list()
                xlist[[1]] = data.frame('variant.id'=seqGetData(f,"variant.id"),
                                        'allele.index'=rep(1,length(seqGetData(f,"variant.id")))
                                        )
                system.time({
                    collapse.results <- assocTestSeq( genotype.data, 
                                         nullmod, 
                                         xlist, 
                                         weight.beta = weights,
                                         test=test.type, 
                                         burden.test=test.stat,verbose=F)
                })
                generes <- cbind(data.frame(gene=current.gene, num.variants= nrow(xlist[[1]])),cmac=sum(mac[filtered.alleles & maf > 0]), collapse.results$results)
                generes$start=kg[gidx,]$start - BUFFER
                generes$stop=kg[gidx,]$stop + BUFFER
                generes$chr=kg[gidx,]$chr
                generes$cmaf=sum(maf[filtered.alleles & maf > 0])
            
                ##} else if(test.stat == 'Firth') {  # Single variant test
                ##   system.time({generes <- regression(genotype.data,
                ##                                      outcome=outcome.name,
                ##                                      model.type =  tolower(test.stat))})
                ##   generes$gene <- current.gene
                ##   generes$pos=pos[snp_idx[filtered.alleles & maf > 0]]
            } else {  # Single variant test
 
                system.time({generes <- assocTestMM(genotype.data, 
                                              nullmod, snp.block.size = 2000,
                                              test = test.stat,verbose=FALSE)})
                generes$gene <- current.gene
                generes$pos=pos[snp_idx[filtered.alleles & maf > 0]]
                generes$snpID = paste0(generes$chr,':',generes$pos)
                generes$ref = sapply(seqGetData(f, "allele"),FUN=function(x)strsplit(x,',')[[1]][1])
                generes$alt = sapply(seqGetData(f, "allele"),FUN=function(x)strsplit(x,',')[[1]][2])      
          
            }
        }
    }
    seqClose(f)
  #}
             
  if(!exists("generes")){
    generes= data.frame('gene'=current.gene);
  }
  
  generes
}
cat("\nFinished Loop\n")
cat("Total output lines: ",nrow(sm_obj),"\n")

  
out <- gzfile(output.file, "w")
write.csv(sm_obj, out, row.names=F)
close(out)
