#== Args 
args<-commandArgs(TRUE)
#===mandatory parameters
phenotype.file <- args[1]
outcome.name <- args[2]
outcome.type <-  args[3]
covariate.string <- args[4]
genotype.files <- args[5]
output.file <- args[6]

#==optional parameters
kinship.matrix <- args[7]
pheno.id <- args[8]


# added these to JSON
test.stat <-  args[9] # Score, Wald, Firth
conditional <- args[10]
het_vars <-  args[11]



# GLOBAL VARIABLES
collapsing.tests <- c("SKAT",  "Burden")
test.stat.vals <- c("Score", "Wald", "F")

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

cat('output.file',output.file,'\n')
cat('kinship.matrix',kinship.matrix,'\n')
cat('test.stat',test.stat,'\n')
cat('outcome.type',outcome.type,'\n')
cat('het_vars',het_vars,'\n')
cat('conditional',conditional,'\n')


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
suppressMessages(library(GENESIS))
suppressMessages(library(data.table))
#suppressMessages(library(parallel))
## Setup
source("/home/dnanexus/pipelineFunctions.R")

covariates <- split.by.comma(covariate.string)  

## phenotype 
phenotype.data <- read.csv(phenotype.file, header=TRUE, as.is=TRUE)
# add sex if provided, needed for X chrom MAF calcs
if(!'sex' %in% covariates && 'sex' %in% names(phenotype.data)){
    covariates = c(covariates,'sex')
}else{
    if(! 'sex' %in% names(phenotype.data)) warning("A column labeled 'sex' coded M/F should be included when performing analyses on the sex chromosomes")
} 
if(NCOL(phenotype.data) < 2) warning("Is the phenotype file a CSV?  Too few columns from read.csv()") 



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



## Conditional analaysis
if(conditional != 'NA'){
    cond_snps = strsplit(conditional,',')[[1]] 
    cond_snps = gsub(' ','',cond_snps)
    cond_snps = gsub('"','',cond_snps)

   
    chr_array = seqGetData(f, "chromosome")
    pos_array = seqGetData(f, "position")
    allele_array = seqGetData(f, "allele")
    allele_array = gsub(',',':',allele_array)
    snp_array = paste(chr_array,pos_array,allele_array,sep=':')

    
    cat('Conditioning on ',conditional,'...\n')
    cidx = which(snp_array %in% cond_snps)
    NCOND = length(cidx)
    cat('Found ',snp_array[cidx], ' in GDS file\n')
    if(NCOND < length(cond_snps)){
        warning('NOT ALL CONDITIONAL SNPS FOUND IN GDS')
    }

    
    if(any(cidx)){
        seqSetFilter(f,variant.sel=cidx, sample.id = row.names(pheno),verbose=FALSE)
        ncol = NCOL(pheno)
        pheno = cbind(pheno,as.data.frame(altDosage(f,  use.names=FALSE)))
        condheaders  = paste0('csnp',1:NCOND)
        colnames(pheno)[(ncol+1):(ncol+NCOND)] = condheaders
    }else{
        stop('Can not find snp ',conditional,' with position ',cpos,' to condition on in data file')
    }
  
    dropConditionalCases = NROW(pheno)-NROW(pheno[complete.cases(pheno),])
    if(dropConditionalCases > 0){
        cat('Warning: Dropping ',dropConditionalCases,' samples due to missing conditional genotype calls\n')
    }

    pheno = pheno[complete.cases(pheno),]
  
    covariates[(length(covariates) + 1):(length(covariates) +NCOND)] <- condheaders
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
save(nullmod,annotated.frame,file='results')
