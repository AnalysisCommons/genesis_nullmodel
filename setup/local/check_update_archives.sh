#!/bin/bash 
md5_checksum_file="$1"
archive_directory=${2%/}
dna_nexus_archives="Commons:/Tools/Archive/" 

echo_if_does_not_exist() {
    file="$1"
    if [ ! -e  ${file} ]
    then
        echo "${file}"
    fi
}
export -f echo_if_does_not_exist

update_basedir() {
    file="$1"
    dirname="$2"
    sed 's|AUX_DIRNAME|'${dirname}'|' ${file}
}



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



if md5sum -c --status ${md5_checksum_file} &> /dev/null
then
    echo "Archives don't need to be updated"
    exit 0
else
    # Touch missing files
    update_basedir ${md5_checksum_file} ${archive_directory} |\
        awk '{print $2}' |\
        xargs -I{} bash -c 'echo_if_does_not_exist "$1"' _ '{}' |\
        xargs -I{} touch '{}' 
    archive_files_to_update=$(md5sum -c <(update_basedir ${md5_checksum_file} ${archive_directory}) | grep FAILED | awk -F: '{sub(".*/","",$1); print $1}')
fi

if [ "x${archive_files_to_update}" != "x" ]
then 
    if dx --help  &> /dev/null;
    then
        if dx ls "${dna_nexus_archives}" &> /dev/null;
        then
            for archive_file in ${archive_files_to_update}
            do
                dx download -f "${dna_nexus_archives}/${archive_file}" -o "${archive_directory}"
            done
        else
          echo "No access to the DNA NEXUS ARCHIVE DIR: ${dna_nexus_archives}"
          exit 1 
        fi
    else
       echo "Install DNA NEXUS DX toolkit and rerun"
       exit 1
    fi
else
    echo "The archives are current!"
fi
