#!/bin/bash 

md5_checksum_file="$1"
archive_directory="$2"


if [ ! -e ${md5_checksum_file} ]
then
   echo "File  ${md5_checksum_file} is missing"
   exit 1
fi

if [ ! -d ${archive_directory} ]
then
   echo "Directory  ${archive_directory} is missing"
   exit 1
fi

files_to_remove=$( find "${archive_directory}" -iname "*.tar.gz" | grep -v -f <( basename -a $(awk '{print $2}' ${md5_checksum_file})))
echo "The following archives ${files_to_remove} are going to be removed; continue [y/N]:" 
read var_continue
if [ "${var_continue}" == "y" ] || [ "${var_continue}" == "Y" ]
then 
   rm -v  ${files_to_remove}
else
   echo "No archives are going to be removed"
fi
