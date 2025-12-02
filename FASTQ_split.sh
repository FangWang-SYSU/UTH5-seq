

#!/bin/bash

inputdir=/mRNA
# out fastq file output
outputFastq=RS240531_0619/mRNA

BarcodeFile==mRNA_index_barcode.txt



find "$inputdir" -mindepth 1 -type d | while read -r directory; do
    echo "Processing directory: $directory"
    directory_base=$(basename "$directory")
    outputdir=${outputFastq}/${directory_base}
    mkdir -p ${outputdir}
    
    file_r1=$(find "$directory" -maxdepth 1 -type f -name "*_R1.fastq.gz")
    file_r2=$(find "$directory" -maxdepth 1 -type f -name "*_R2.fastq.gz")

    # Check if both files are found
    if [[ -f "$file_r1" && -f "$file_r2" ]]; then
        echo "Files found: $file_r1 and $file_r2"

        file_r1_name=$(basename "$file_r1")
        file_r2_name=$(basename "$file_r2")
        mv ${file_r1} ${outputdir}/${file_r1_name}
        mv ${file_r2} ${outputdir}/${file_r2_name}

        barcode_splitter --bcfile $BarcodeFile $outputdir/${file_r1_name} $outputdir/${file_r2_name}  --split_all --mismatches 0 --idxread 2 --suffix .fastq --gzipout 2

    else
        echo "Files not found in $directory"
    fi

    echo "-----------------------------------------"
done


