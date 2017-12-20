#!/bin/bash

main() {

    echo "Value of phenofile: '$phenofile'"
    echo "Value of snpinfofile: '$snpinfofile'"
    echo "Value of genotypefile: '$genotypefile'"
    echo "Value of outputfilename: '$outputfilename'"

    echo "Value of kinshipmatrix: '$kinshipmatrix'"
    echo "Value of pheno_id: '$pheno_id'"
    echo "Value of snpNames: '$snpNames'"

    echo "Value of test_stat: '$test_stat'"


    covariate_list=$(echo ${covariate_list} | sed 's/ //g') 
    echo "Value of covariate_list: '$covariate_list'"
    echo "value of conditional SNP: '$conditional'"    
    echo "Value of het_vars: '$het_vars'"


    dx download "$phenofile" -o phenofile &
    dx download "$genotypefile" -o genotypefile &
    dx download "$kinshipmatrix"   &


    # install R
    # add R to path 
    
    echo "INSTALLING GENESIS"
    make >> /dev/null & 
    export PATH=/opt/R/bin/:${PATH}
    export MKL_NUM_THREADS=1
    wait

    echo "KINSHIP"
    kinshipmatrix_filename=$( dx describe --name "$kinshipmatrix" )
    
 

    sudo chmod o+rw /tmp
    # wait if debug 
    if [ ${debug} -ne 0 ]
    then
       echo "DEBUG is on sleeping for ${debug}h"
       sleep ${debug}h
    fi
    wait
    echo "Checking phenofile" 
    if [ -e phenofile ] 
    then
       head -n1 phenofile
    else
       echo "The phenofile is not ready"
    fi
 
    echo "Rscript genesis_nullmodel.R phenofile $outcome_name $outcome_type \"$covariate_list\" genotypefile results $kinshipmatrix_filename $pheno_id  $test_stat  $conditional $het_vars"
    echo "Running code"
    Rscript genesis_nullmodel.R phenofile $outcome_name $outcome_type \"$covariate_list\" genotypefile results $kinshipmatrix_filename $pheno_id $test_stat $conditional $het_vars
    echo "Finished running code"
    results=$(dx upload results --brief)
    dx-jobutil-add-output results "$results" --class=file
    dx mv ${results} ${outputfilename}.Rda
}
