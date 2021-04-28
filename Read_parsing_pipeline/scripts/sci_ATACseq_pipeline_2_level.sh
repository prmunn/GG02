#!/usr/bin/bash
#SBATCH -J ATACseq
#SBATCH -o %x.out
#SBATCH -n 12
#SBATCH --mem-per-cpu=18000

# Note: parts of this script are also called from the UHT_ATACseq-pipeline 
# script - if you make changes, make sure they also work for that pipeline
source ~/.bash_profile

usage(){

    echo "s c i - A T A C - S E Q   P I P E L I N E"
    echo
    echo

    echo "Usage: bash" $0 "-f arg -g arg [-h arg] [-p arg] [-d args] [-t arg] [-q arg] [-s arg] [-a arg] [-c arg] [-u arg] [-e arg]"
    echo
    echo "Example: bash" $0 "-f GIA0130-AdrianBarnyard/ -g hg38_and_mm10 -s findMESequence"
    echo
    echo "---------------------------------------------------------------------------------------------------------------"
    echo " -f  --> Folder containing data files (required)"
    echo " -g  --> Reference Genome <mm10 or hg38 or hg38_and_mm10> (required)"
    echo "[-h] --> Display Help"
    echo "[-p] --> Project Identifier Number"
    echo "[-d] --> Comma Spearated Values for Delimiter and Field <delim,field or default> default: _,5 "
    echo "[-t] --> Trimming <nextseq or nova>;"
    echo "[-q] --> Execute atacQC.R script <yes>"
    echo "[-s] --> Execute single step <step name>"
    echo "[-a] --> Add duplicate barcode(s) to 5' end of R2 read <i5 or i7i5>, default: No barcode added"
    echo "[-c] --> Define <sample>:<cell number> comma separated key / value pairs for whitelist command. e.g. sample1:2000,sample2:1600"
    echo "[-u] --> Get parameters from <filename>. Note: Any parameters specified on the command line will override the contents of this file"
    echo "[-e] --> Send email to this address when function completes"
    echo
    echo "---------------------------------------------------------------------------------------------------------------"
    echo "Assumes that the following environment variables have been set up"
    echo "For UMI_tools:"
    echo "export PYTHONPATH=/programs/UMI-tools/lib/python3.6/site-packages:/programs/UMI-tools/lib64/python3.6/site-packages"
    echo "export PATH=/programs/UMI-tools/bin:\$PATH"
    echo "For macs2:"
    echo "export PYTHONPATH=/programs/macs2-2.2.7.1/lib64/python3.6/site-packages"
    echo "export PATH=/programs/macs2-2.2.7.1/bin:\$PATH"
    echo "For annotatePeaks.pl"
    echo "export PATH=/workdir/tools/homer/bin:\$PATH"
    echo
    echo "---------------------------------------------------------------------------------------------------------------"
}


declare -A genomeDir

genomeDir=(
["mm10"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/BWAIndex/genome.fa" \
["hg38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/Homo_sapiens.GRCh38.dna.toplevel.fa" \
["hg38_and_mm10"]="/workdir/singleCellData/10x_reference_files/refdata-cellranger-atac-GRCh38-and-mm10-1.2.0/fasta/genome.fa"
)
# Removed /bwa.index from the hg38 directory reference above
# ["hg38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/bwa.index/Homo_sapiens.GRCh38.dna.toplevel.fa" \

declare -A gtfs

gtfs=(
["mm10"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/Mus_musculus.GRCm38.96.gtf" \
["hg38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/Homo_sapiens.GRCh38.96.gtf" \
["hg38_and_mm10"]="/workdir/singleCellData/10x_reference_files/refdata-cellranger-atac-GRCh38-and-mm10-1.2.0/genes/genes.gtf"
)

declare -A gAlias # for compatibility with atacQC.R
gAlias=(
["mm10"]="mouse" \
["hg38"]="human" \
["hg38_and_mm10"]="human and mouse"
)

declare -A sample_cell_number_dict

ORIG_FASTQ="orig-fastqs"

# Step 1) Use cutadapt to find the position of the ME sequence in the R2 read 
# (no modification of reads, just generate info file with ME position)
findMESequence(){

    cd "${FOLDER}"/fastqs/
    mkdir parsed_BCs
    cd parsed_BCs

    RELATIVE_ORIG_FASTQ=../"${ORIG_FASTQ}"

    # Delete info and summary files from any prior runs
    if [ -f R2_info.txt ]; then
        rm -f R2_info.txt
    fi
    if [ -f R2_info_summary.txt ]; then
        rm -f R2_info_summary.txt
    fi

    # Set file names
    input_R2_fastq_file=("${RELATIVE_ORIG_FASTQ}"/*_R2_001.fastq.gz)

    echo
    echo "Running module findMESequence()"
    echo "input_R2_fastq_file: "$input_R2_fastq_file
    echo

    # Find ME position in R2 with cutadapt (info file only)
    cutadapt -g AGATGTGTATAAGAGACAG \
    --info-file R2_info.txt \
    -e 0.11 \
    -j 10 \
    --match-read-wildcards \
    --action none \
    -o /dev/null \
    $input_R2_fastq_file

    # Summary counts for columns 2, 3, and 4
    awk -F '\t' '{a[$2]++; if($2>=0) {b[$3]++; c[$4]++;} else next;} END {print "col2=mismatch"; for(i in a) print a[i],i; print "\ncol3=startpost"; for(j in b) print b[j],j; print "\ncol4=endpos"; for(k in c) print c[k],k;}' \
    R2_info.txt > R2_info_summary.txt
 
    # Clean up
    # none

    cd ../../..

}

# Step 2) Merge I1 and I2 into R2 and pad UMIs with custom shell script 
horizontalMerge_padUMI_addDupBC(){

    cd "${FOLDER}"/fastqs/parsed_BCs/

    RELATIVE_ORIG_FASTQ=../"${ORIG_FASTQ}"

    # Set file names
    input_R2_fastq_file=("${RELATIVE_ORIG_FASTQ}"/*_R2_001.fastq.gz)
    input_R1_fastq_file=("${RELATIVE_ORIG_FASTQ}"/*_R1_001.fastq.gz)
    input_I1_fastq_file=("${RELATIVE_ORIG_FASTQ}"/*_I1_001.fastq.gz)
    input_I2_fastq_file=("${RELATIVE_ORIG_FASTQ}"/*_I2*.f*.gz)
    file_prefix=`echo $(basename $input_R2_fastq_file .gz) | cut -d "_" -f 1`

    echo
    echo "Running module horizontalMerge_padUMI_addDupBC()"
    echo "file_prefix: "$file_prefix
    echo "input_R2_fastq_file: "$input_R2_fastq_file
    echo "input_R1_fastq_file: "$input_R1_fastq_file
    echo "input_I1_fastq_file: "$input_I1_fastq_file
    echo "input_I2_fastq_file: "$input_I2_fastq_file
    echo "ADD_DUPLICATE_BC: "$ADD_DUPLICATE_BC
    echo

    # Set output file name based on barcodes added to 5' end of R2 read
    if [[ $ADD_DUPLICATE_BC == "i5" ]]; then
        merge_output_fastq_file="${file_prefix}"_I2_I2_I1_padUMI_R2.fastq
    elif [[ $ADD_DUPLICATE_BC == "i7i5" ]]; then
        merge_output_fastq_file="${file_prefix}"_i7tag_I2_I2_I1_padUMI_R2.fastq
    else
        merge_output_fastq_file="${file_prefix}"_I2_I1_padUMI_R2.fastq
    fi

    awk -v addDuplicateBC="$ADD_DUPLICATE_BC" -f ../../../horizontalMerge-padUMI-addDupBC.awk R2_info.txt \
    <(paste <(zcat $input_I2_fastq_file) <(zcat $input_I1_fastq_file) <(zcat $input_R2_fastq_file)) \
    > $merge_output_fastq_file

    # Extra demultiplexing step in needed
    if [[ $ADD_DUPLICATE_BC == "i5" ]]; then
        echo "Demultiplexing on i5tag"
        cutadapt -e 0.15 \
        --no-indels \
        -g file:"${RELATIVE_ORIG_FASTQ}"/i5tagBC_barcodes.fa 
        -o {name}_I2_I1_padUMI_R2.fastq 
        -p {name}_R1.fastq \
        $merge_output_fastq_file \
        $input_R1_fastq_file \
        > "${file_prefix}"_i5tagBC_demux.log

    elif [[ $ADD_DUPLICATE_BC == "i7i5" ]]; then
        echo "Demultiplexing on i7tag and i5tag"
        cutadapt -e 0.15 \
        --no-indels \
        -g file:"${RELATIVE_ORIG_FASTQ}"/i7tagBC_i5tagBC_barcodes.fa \
        -o {name}_I2_I1_padUMI_R2.fastq \
        -p {name}_R1.fastq \
        $merge_output_fastq_file \
        $input_R1_fastq_file \
        > "${file_prefix}"_i7tagBC_i5tagBC_demux.log

    else
        echo "No extra demultiplexing needed"
        # No demultiplexed R1 files, so copy the original R1 fastq and gzip the padded R2 fastq
        cp $input_R1_fastq_file "${file_prefix}"_R1.fastq.gz
        gzip $merge_output_fastq_file
    fi

    # If extra demultiplexing step executed, then concatenate by sample of origin
    if [[ $ADD_DUPLICATE_BC == "i5" || $ADD_DUPLICATE_BC == "i7i5" ]]; then

        mkdir demuxed_fastqs
        mv *~*.fastq demuxed_fastqs
        cd demuxed_fastqs

        ls -1 *~*_R2.fastq | cut -d "~" -f 1 | uniq > samples.list

        readarray samples < samples.list

        for i in "${samples[@]}"
        do
            # First trim blank chrs from ends of i
            trimmed_i=`echo -e $i | awk '{$1=$1;print}'`
            # Then concatenate files with the same sample name
            cat $trimmed_i*_R2.fastq | gzip > ../"${trimmed_i}"_I2_I1_padUMI_R2.fastq.gz
            cat $trimmed_i*_R1.fastq | gzip > ../"${trimmed_i}"_R1.fastq.gz
        done

        cd ..

        # Gzip the last of the demuxed files
        if [ -f unknown_I2_I1_padUMI_R2.fastq ]; then
            gzip unknown_I2_I1_padUMI_R2.fastq
            gzip unknown_R1.fastq
        fi

    fi

    cd ../../..

    # Notify user that process complete
    if [[ ! -z "${EMAIL_TO}" ]]; then
        subject="Message from "$0
        message=${FUNCNAME[0]}" function, running for "$FOLDER" dataset, completed"
        python send_email_v3.py $EMAIL_TO "${subject}" "${message}"
    fi
}

# Step 3a) Optional: generate combBC whitelist with UMI_tools whitelist
generate_whitelist() {

    # Assumes that environment variables already set up for UMI_tools
    # export PYTHONPATH=/programs/UMI-tools/lib/python3.6/site-packages:/programs/UMI-tools/lib64/python3.6/site-packages
    # export PATH=/programs/UMI-tools/bin:$PATH

    echo
    echo "Running module generate_whitelist()"
    echo

    cd "${FOLDER}"/fastqs/parsed_BCs/

    # Loop thru gziped pad UMI fastq files
    for I2_I1_padUMI_R2_fastq_file in *_R2.fastq.gz
    do
        # Set file prefix
        file_prefix=`echo $(basename $I2_I1_padUMI_R2_fastq_file .gz) | cut -d "_" -f 1`

        echo "file_prefix: "$file_prefix
        echo "I2_I1_padUMI_R2_fastq_file: "$I2_I1_padUMI_R2_fastq_file

        # Run umi_tools whitelist to generate cell barcode whitelist (specific to dataset). 
        # This step may require testing different options to get the optimal whitelist 
        # (e.g. --error-correct-threshold, --knee-method=[density/distance], --method=[reads/umis]).

        # If cell number defined for this sample then use it, otherwise use distance knee method
        # for key in "${!sample_cell_number_dict[@]}"; do
        #    echo "$key ${sample_cell_number_dict[$key]}"
        # done
        if [[ ! -z ${sample_cell_number_dict[${file_prefix}]} ]]; then
            echo "Using density knee method with cell number = "${sample_cell_number_dict[$file_prefix]}

            umi_tools whitelist --knee-method=density \
            --set-cell-number=${sample_cell_number_dict[$file_prefix]} \
            --method=reads \
            --plot-prefix "${file_prefix}"_predictBC \
            --allow-threshold-error \
            --extract-method regex \
            --bc-pattern='(?P<cell_1>.{16})(?P<umi_1>.{10})(?P<cell_2>.{10})' \
            --error-correct-threshold=2 \
            --ed-above-threshold=correct \
            -L "${file_prefix}"_predictedBCwhitelist.log \
            -I $I2_I1_padUMI_R2_fastq_file \
            -S "${file_prefix}"_predictedBCwhitelist.txt

        else
            echo "Using distance knee method with no cell number"

            umi_tools whitelist --knee-method=distance \
            --method=reads \
            --plot-prefix "${file_prefix}"_predictBC \
            --allow-threshold-error \
            --extract-method regex \
            --bc-pattern='(?P<cell_1>.{16})(?P<umi_1>.{10})(?P<cell_2>.{10})' \
            --error-correct-threshold=2 \
            --ed-above-threshold=correct \
            -L "${file_prefix}"_predictedBCwhitelist.log \
            -I $I2_I1_padUMI_R2_fastq_file \
            -S "${file_prefix}"_predictedBCwhitelist.txt
        fi
    done

    cd ../../..

    # Notify user that process complete
    if [[ ! -z "${EMAIL_TO}" ]]; then
        subject="Message from "$0
        message=${FUNCNAME[0]}" function, running for "$FOLDER" dataset, completed"
        python send_email_v3.py $EMAIL_TO "${subject}" "${message}"
    fi
}

# Step 3b) Run umi_tools extract to move cell barcode and UMI to R1/R2 headers. 
# R2 is the ‘driver’ because it contains the BC and UMI sequences
run_extract() {

    # Assumes that environment variables already set up for UMI_tools
    # export PYTHONPATH=/programs/UMI-tools/lib/python3.6/site-packages:/programs/UMI-tools/lib64/python3.6/site-packages
    # export PATH=/programs/UMI-tools/bin:$PATH

    echo
    echo "Running module run_extract()"
    echo

    cd "${FOLDER}"/fastqs/parsed_BCs/

    # Loop thru gziped pad UMI fastq files
    for I2_I1_padUMI_R2_fastq_file in *_R2.fastq.gz
    do
        # Set file prefix
        file_prefix=`echo $(basename $I2_I1_padUMI_R2_fastq_file .gz) | cut -d "_" -f 1`

        # Set file names
        input_R1_fastq_file=("${file_prefix}"_R1.fastq.gz)
        predictedBCwhitelist_file=("${file_prefix}"_predictedBCwhitelist.txt)

        echo "file_prefix: "$file_prefix
        echo "I2_I1_padUMI_R2_fastq_file: "$I2_I1_padUMI_R2_fastq_file
        echo "input_R1_fastq_file: "$input_R1_fastq_file
        echo "predictedBCwhitelist_file: "$predictedBCwhitelist_file
        echo

        umi_tools extract --extract-method=regex \
        -p '(?P<cell_1>.{16})(?P<umi_1>.{10})(?P<cell_2>.{10})' \
        --filtered-out="${file_prefix}"_extract_filtered_out.txt \
        --filtered-out2="${file_prefix}"_extract_filtered_out2.txt \
        --error-correct-cell \
        --quality-filter-mask=20 \
        --quality-encoding=phred33 \
        --whitelist=$predictedBCwhitelist_file \
        -I $I2_I1_padUMI_R2_fastq_file \
        -S "${file_prefix}"_hBC_UMI_R2.fastq.gz \
        --read2-in=$input_R1_fastq_file \
        --read2-out="${file_prefix}"_hBC_UMI_R1.fastq.gz \
        -L "${file_prefix}"_extractBC.log

        # Old code - remove after testing new code
        # umi_tools extract --extract-method=string \
        # --filter-cell-barcode \
        # --error-correct-cell \
        # --quality-filter-mask 20 \
        # --quality-encoding phred33 \
        # -p CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
        # --whitelist $predictedBCwhitelist_file \
        # -I $I2_I1_padUMI_R2_fastq_file \
        # -S "${file_prefix}"_hBC_UMI_R2.fastq.gz \
        # --read2-in=$input_R1_fastq_file \
        # --read2-out="${file_prefix}"_hBC_UMI_R1.fastq.gz \
        # -L "${file_prefix}"_extractBC.log

    done

    mkdir extract_fastqs
    mv *_hBC* extract_fastqs
    mv *_extractBC.log extract_fastqs
    mv *_extract_filtered_out* extract_fastqs
    mv extract_fastqs ..

    cd ../../..

    # Notify user that process complete
    if [[ ! -z "${EMAIL_TO}" ]]; then
        subject="Message from "$0
        message=${FUNCNAME[0]}" function, running for "$FOLDER" dataset, completed"
        python send_email_v3.py $EMAIL_TO "${subject}" "${message}"
    fi
}

# Step 4) Use cutadapt to trim (and require) 5' anchored ME in R2; write untrimmed paired reads to separate 
# file (do not use). Again, R2 is the ‘driver’ because it should now contain the ME sequence at the 5’ end
cutadapt_trim5pME(){

    echo
    echo "Running module cutadapt_trim5pME()"
    echo

    cd "${FOLDER}"/fastqs/extract_fastqs/

    # Delete output files from any prior runs
    #rm -f *_noME_hBC_UMI_*
    #rm -f *_hasME_hBC_UMI_*

    # Loop thru gziped pad UMI fastq files
    for hBC_UMI_R2_fastq_file in *_hBC_UMI_R2.fastq.gz
    do
        # Set file prefix
        file_prefix=`echo $(basename $hBC_UMI_R2_fastq_file .gz) | cut -d "_" -f 1`

        # Set file names
        hBC_UMI_R1_fastq_file=("${file_prefix}"_hBC_UMI_R1.fastq.gz)

        echo "file_prefix: "$file_prefix
        echo "hBC_UMI_R2_fastq_file: "$hBC_UMI_R2_fastq_file
        echo "hBC_UMI_R1_fastq_file: "$hBC_UMI_R1_fastq_file

        # Require 5'ME in R2
        cutadapt -g ^AGATGTGTATAAGAGACAG \
        -e 0.11 \
        -j 10 \
        --no-indels \
        --match-read-wildcards \
        --untrimmed-output "${file_prefix}"_noME_hBC_UMI_R2.fastq.gz \
        --untrimmed-paired-output "${file_prefix}"_noME_hBC_UMI_R1.fastq.gz \
        -o "${file_prefix}"_hasME_hBC_UMI_R2.fastq.gz \
        -p "${file_prefix}"_hasME_hBC_UMI_R1.fastq.gz \
        $hBC_UMI_R2_fastq_file \
        $hBC_UMI_R1_fastq_file \
        > "${file_prefix}"_trim5pME.log

    done

    mkdir hasME_hBC_fastqs
    mv *_hasME_hBC_* hasME_hBC_fastqs
    mv *_noME_hBC_* hasME_hBC_fastqs
    mv hasME_hBC_fastqs ..

    cd ../../..

    # Notify user that process complete
    if [[ ! -z "${EMAIL_TO}" ]]; then
        subject="Message from "$0
        message=${FUNCNAME[0]}" function, running for "$FOLDER" dataset, completed"
        python send_email_v3.py $EMAIL_TO "${subject}" "${message}"
    fi
}


# Step 5) Trim any internal 3' ME-REV in R1 or R2 due to small inserts, trim for quality, filter for short reads (standard cutadapt).
# Adjust 3’-quality trimming mode (-q/--nextseq-trim) based on sequencer (4-color/2-color).
# R2 is listed first because that has been the pattern in this pipeline, but the order of R1/R2 doesn’t matter (R1/R2 adaptors are the same).
cutadapt_trim3pAd(){

    echo
    echo "Running module cutadapt_trim3pAd()"
    echo

    cd "${FOLDER}"/fastqs/hasME_hBC_fastqs/

    #//--- need to add fastQC to this function
    mkdir fastQC

    # Loop thru gziped pad UMI fastq files
    for hasME_hBC_UMI_R2_fastq_file in *hasME_hBC_UMI_R2.fastq.gz
    do
        # Set file prefix
        file_prefix=`echo $(basename $hasME_hBC_UMI_R2_fastq_file .gz) | cut -d "_" -f 1`

        # Set file names
        hasME_hBC_UMI_R1_fastq_file=("${file_prefix}"_hasME_hBC_UMI_R1.fastq.gz)

        echo "file_prefix: "$file_prefix
        echo "hasME_hBC_UMI_R2_fastq_file: "$hasME_hBC_UMI_R2_fastq_file
        echo "hasME_hBC_UMI_R1_fastq_file: "$hasME_hBC_UMI_R1_fastq_file

        # Use cutadapt to trim R1/R2 for 3’ME-rev (e.g. small inserts) and for 3’quality 
        # (standard 3’adaptor/quality trimming for PE reads). Then move to mapping, peak calling steps.
        cutadapt -a CTGTCTCTTATACACATCT \
        -A CTGTCTCTTATACACATCT \
        -q 20 \
        -m 50 \
        -j 10 \
        --match-read-wildcards \
        -o "${file_prefix}"_trimmed_val_2.fastq.gz \
        -p "${file_prefix}"_trimmed_val_1.fastq.gz \
        $hasME_hBC_UMI_R2_fastq_file \
        $hasME_hBC_UMI_R1_fastq_file \
        > "${file_prefix}"_trim3pAd.log

    done

    #for i in "${fastqs[@]}"
    #do
    #    trim_galore --quality 20 --gzip --length 20  --paired --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
    #done

    mkdir TrimQC_stats trimmed_fastqs
    #mv *_trimming_report.txt TrimQC_stats
    mv *_val* trimmed_fastqs
    mv TrimQC_stats fastQC trimmed_fastqs ..

    cd ../../..

    # Notify user that process complete
    if [[ ! -z "${EMAIL_TO}" ]]; then
        subject="Message from "$0
        message=${FUNCNAME[0]}" function, running for "$FOLDER" dataset, completed"
        python send_email_v3.py $EMAIL_TO "${subject}" "${message}"
    fi
}

#trimPE(){

#    cd fastqs
#    ls -1 *_R1.fastq* > .R1
#    ls -1 *_R2.fastq* > .R2
#    paste -d " " .R1 .R2 > Reads.list

#    readarray fastqs < Reads.list
#    mkdir fastQC

#    for i in "${fastqs[@]}"
#    do
#        $TRIM --nextseq 20 --length 20  --paired --gzip --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
#    done

#    mkdir TrimQC_stats trimmed_fastqs
#    mv *_trimming_report.txt TrimQC_stats
#    mv *_val* trimmed_fastqs
#    mv TrimQC_stats fastQC trimmed_fastqs ../

#    cd ..
#}

#trimHiSeqPE(){

#    cd fastqs
#    ls -1 *_1.fq* > .R1
#    ls -1 *_2.fq* > .R2
#    paste -d " " .R1 .R2 > Reads.list

#    readarray fastqs < Reads.list
#    mkdir fastQC

#    for i in "${fastqs[@]}"
#    do
#        trim_galore --quality 20 --gzip --length 20  --paired --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
#    done

#    mkdir TrimQC_stats trimmed_fastqs
#    mv *_trimming_report.txt TrimQC_stats
#    mv *_val* trimmed_fastqs
#    mv TrimQC_stats fastQC trimmed_fastqs ..

#    cd ..
#}

# Map reads with bwa mem
alignPE(){

    cd "${FOLDER}"/fastqs/trimmed_fastqs/

    echo
    echo "Running module alignPE()"
    echo

    # Set up pairs of fastq files
    ls -1 *_val_1.f*.gz > .trR1
    ls -1 *_val_2.f*.gz > .trR2
    paste -d " " .trR1 .trR2 > Trimmed.list

    readarray trimmedFastqs < Trimmed.list

    for i in "${trimmedFastqs[@]}"

    do
        iSUB=`echo $i | cut -d "_" -f 1`
        #iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`

        echo "file_prefix: "$iSUB
        echo "Read 1 and Read 2 fastq files: "$i

        bwa mem -v 3 -t 24 -M -R "@RG\tID:${iSUB}\tSM:${iSUB}\tPL:ILLUMINA\tLB:${iSUB}\tPU:1" ${genomeDir[${DIR}]} $i \
        2> ${iSUB}_bwa.log \
        | samtools view -@ 24 -b -h -F 0x0100 -O BAM -o ${iSUB}.bam
    done

    mkdir primary_BAMS
    mv *.bam primary_BAMS
    mv primary_BAMS ..

    cd ../../..

    # Notify user that process complete
    if [[ ! -z "${EMAIL_TO}" ]]; then
        subject="Message from "$0
        message=${FUNCNAME[0]}" function, running for "$FOLDER" dataset, completed"
        python send_email_v3.py $EMAIL_TO "${subject}" "${message}"
    fi
}

sort(){

    cd "${FOLDER}"/fastqs/primary_BAMS

    echo
    echo "Running module sort()"
    echo

    for i in *.bam
    do
        echo "Sorting file: "$i
        samtools sort $i > `echo  $i | cut -d "." -f1`_sorted.bam
    done

    for i in *_sorted.bam
    do
        samtools index $i
    done

    # alignment stats etc. on raw bams
    for i in *_sorted.bam
    do
        iSUB=`echo $i | cut -d "_" -f1`
        samtools flagstat $i > ${iSUB}_primary.flagstat
        samtools idxstats $i > ${iSUB}_primary.idxstats
    done

    cd ../../..

}

# Remove mitocondrial reads
rmMT(){

    cd "${FOLDER}"/fastqs/primary_BAMS

    echo
    echo "Running module rmMT()"
    echo

    for i in *_sorted.bam
    do

        iSUB=`echo $i | cut -d "_" -f1`

        echo "file_prefix: "$iSUB
        echo "BAM file: "$i

        # //--- make sure the grep still works after changing '.' to '_' in the input file
        samtools view -H `ls -1 *_sorted.bam | head -1` | cut -f2 | grep "SN:" |  cut -d ":" -f2 | grep -v "MT\|_\|\." | xargs samtools view -b $i > ${iSUB}_noMT.bam

    done

    cd ../../..
}

# Mark duplicates
# markDups(){

#    cd "${FOLDER}"/fastqs/primary_BAMS

#    echo
#    echo "Running module markDups()"
#    echo

#    for i in *_noMT.bam
#    do
#        iSUB=`echo $i | cut -d "_" -f1`

#        echo "file_prefix: "$iSUB
#        echo "BAM file: "$i

#        java -jar /programs/bin/picard-tools/picard.jar \
#        MarkDuplicates \
#        INPUT=$i \
#        OUTPUT=${iSUB}_dupMarked_noMT.bam \
#        ASSUME_SORTED=true \
#        REMOVE_DUPLICATES=false \
#        METRICS_FILE=${iSUB}_MarkDuplicates_metrics.txt \
#        VALIDATION_STRINGENCY=LENIENT \
#        TMP_DIR=tmp

#    done
#    cd ../../..
# }

# dedupBAM(){

#    cd "${FOLDER}"/fastqs/primary_BAMS

#    echo
#    echo "Running module dedupBAM()"
#    echo

#    # alignment stats etc. on dupMarked no MT bams
#    for i in *_dupMarked_noMT.bam
#    do
#        iSUB=`echo $i | cut -d "_" -f1`

#        echo "file_prefix: "$iSUB
#        echo "BAM file: "$i

#        samtools index $i
#        samtools flagstat $i > ${iSUB}_noMT.flagstat
#        samtools idxstats $i > ${iSUB}_noMT.idxstats
#    done

#    for i in *_dupMarked_noMT.bam
#    do
#        iSUB=`echo $i | cut -d "_" -f1`
#        samtools view -b -h -F 0X400 $i > ${iSUB}_DEDUP.bam
#    done

#    for i in *_DEDUP.bam; do samtools index $i ; samtools idxstats $i > `echo $i | cut -d "_" -f1`_DEDUP.idxstats; done
#    for i in *_DEDUP.bam; do samtools flagstat $i > `echo $i | cut -d "_" -f1`_DEDUP.flagstat; done

#    multiqc -n ${PIN}_multiqc.report .

#    mkdir dedup_BAMS
#    mv *_DEDUP* dedup_BAMS/
#    mv dedup_BAMS ..

#    cd ../../..

# }

dedupBAM_with_UMI_tools(){

    cd "${FOLDER}"/fastqs/primary_BAMS

    echo
    echo "Running module dedupBAM_with_UMI_tools()"
    echo

    for noMT_BAM_file in *_noMT.bam
    do
        # Set file prefix
        file_prefix=`echo $(basename $noMT_BAM_file .bam) | cut -d "_" -f 1`

        echo "file_prefix: "$file_prefix
        echo "noMT_BAM_file: "$noMT_BAM_file

        samtools index $noMT_BAM_file

        # Dedup crashes if you attempt to sort the output
        # Added --no-sort-output and use samtools to do the sorting
        umi_tools dedup --paired \
        --per-cell \
        --no-sort-output \
        --unmapped-reads=discard \
        --chimeric-pairs=discard \
        --unpaired-reads=discard \
        --output-stats="${file_prefix}"_dedup \
        -L "${file_prefix}"_dedup.log \
        -I $noMT_BAM_file \
        -S "${file_prefix}"_DEDUP_unsorted.bam

        samtools sort "${file_prefix}"_DEDUP_unsorted.bam -o "${file_prefix}"_DEDUP.bam

    done

    # Alignment stats etc. on DEDUP bams
    #//--- Do we still need stats on dupMarked no MT bams ???
    for i in *_DEDUP.bam; do samtools index $i ; samtools idxstats $i > `echo $i | cut -d "_" -f1`_DEDUP.idxstats; done
    for i in *_DEDUP.bam; do samtools flagstat $i > `echo $i | cut -d "_" -f1`_DEDUP.flagstat; done

    multiqc -n ${PIN}_multiqc.report .

    mkdir dedup_BAMS
    mv *_DEDUP* dedup_BAMS/
    mv *_dedup* dedup_BAMS/
    mv dedup_BAMS ..

    cd ../../..

    # Notify user that process complete
    if [[ ! -z "${EMAIL_TO}" ]]; then
        subject="Message from "$0
        message=${FUNCNAME[0]}" function, running for "$FOLDER" dataset, completed"
        python send_email_v3.py $EMAIL_TO "${subject}" "${message}"
    fi
}

tagDir(){

    cd "${FOLDER}"/fastqs/dedup_BAMS

    echo
    echo "Running module tagDir()"
    echo

    for i in *_DEDUP.bam
    do
        iSUB=`echo "$i" | cut -d'_' -f1` # subset to rename

        echo "file_prefix: "$iSUB
        echo "BAM file: "$i

        /workdir/tools/homer/bin/makeTagDirectory "$iSUB"_tag.dir "$i" -genome ${genomeDir[${DIR}]}
    done
    cd ../../..
}

callPeak(){

    # Assumes that environment variables already set up for macs2
    # export PYTHONPATH=/programs/macs2-2.2.7.1/lib64/python3.6/site-packages
    # export PATH=/programs/macs2-2.2.7.1/bin:$PATH

    echo
    echo "Running module callPeak()"
    echo

    cd "${FOLDER}"/fastqs/dedup_BAMS

    echo "calling peaks on DEDUP bams"
    mkdir peaks.OUT
    for i  in *_DEDUP.bam
    do
        iSUB=`echo $i | cut -d "_" -f1`

        echo "file_prefix: "$iSUB
        echo "BAM file: "$i

        macs2 callpeak -t $i \
        -f BAMPE \
        -n ${iSUB} \
        -g hs \
        -q 0.05 \
        --outdir peaks.OUT \
        --nomodel --shift 37 --ext 73 \
        --keep-dup all
    done

    cd ../../..

}

mergedPeaks(){

    # Assumes that environment variables already set up for macs2
    # export PYTHONPATH=/programs/macs2-2.2.7.1/lib64/python3.6/site-packages
    # export PATH=/programs/macs2-2.2.7.1/bin:$PATH

   cd "${FOLDER}"/fastqs/dedup_BAMS

    echo
    echo "Running module mergedPeaks()"
    echo

    allBams=`echo *_DEDUP.bam`

    macs2 callpeak -t ${allBams} \
    -f BAMPE \
    -n allSamplesMergedPeakset \
    -g hs \
    -q 0.05 \
    --outdir peaks.OUT \
    --nomodel --shift 37 --ext 73 \
    --keep-dup all

    cd ../../..
}

saf(){

    # awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' ${sample}_peaks.narrowPeak > ${sample}_peaks.saf
    cd "${FOLDER}"/fastqs/dedup_BAMS/peaks.OUT

    echo
    echo "Running module saf()"
    echo

    awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' allSamplesMergedPeakset.peaks.narrowPeak > allSamplesMergedPeakset.saf
    cd ../../../..
}

# Count mapped reads for genomic features
frip(){

    # featureCounts -p -a ${sample}_peaks.saf -F SAF -o readCountInPeaks.txt ${sample}.sorted.marked.filtered.shifted.bam
    cd "${FOLDER}"/fastqs/dedup_BAMS/

    echo
    echo "Running module frip()"
    echo

    for i  in *_DEDUP.bam
    do
        iSUB=`echo $i | cut -d "_" -f1`

        echo "file_prefix: "$iSUB
        echo "BAM file: "$i

        /workdir/tools/subread/bin/featureCounts -p -a peaks.OUT/allSamplesMergedPeakset.saf -F SAF -o "${iSUB}"_readCountInPeaks.txt $i
    done

    cd ../../..
}

annotatePeaks(){

    cd "${FOLDER}"/fastqs/dedup_BAMS/peaks.OUT

    echo
    echo "Running module annotatePeaks()"
    echo

    /workdir/tools/homer/bin/annotatePeaks.pl allSamplesMergedPeakset.saf ${genomeDir[${DIR}]} -gtf ${gtfs[${DIR}]} > allSamplesMergedPeakset_Annotated.saf
    cd ../../../..
}

bedGraphs(){

    cd "${FOLDER}"/fastqs/dedup_BAMS

    echo
    echo "Running module bedGraphs()"
    echo

    for i in *_tag.dir
    do
        makeUCSCfile ${i} -o auto -fsize 1e10 -res 1 -color 106,42,73 -style chipseq
    done

    mkdir tagDirs
    mv *_tag.dir tagDirs
    cd tagDirs
    mkdir bedGraphs

    for i in *_tag.dir
    do
        cd $i
        zcat *.ucsc.bedGraph.gz | awk '{if(NR>1) print "chr"$0; else print $0}' | gzip > `basename *.ucsc.bedGraph.gz .ucsc.bedGraph.gz`.ucsc.bg.gz
        mv *.ucsc.bg.gz ../bedGraphs
        cd ..
    done
    cd ..

    mkdir featureCounts
    mv *.txt featureCounts

    multiqc -n ${PIN}_FRIP_multiqc.report -b "Please note that the featureCounts M Assigned Column refers to Fragments and Not Reads" --ignore tagDirs --ignore peaks.OUT .

    cd ../../..
}

atacQC(){

    cd "${FOLDER}"/fastqs/dedup_BAMS

    echo
    echo "Running module atacQC()"
    echo

    echo "genome alias" = ${gAlias[${DIR}]}
    /programs/R-3.6.3/bin/Rscript /workdir/tools/atacQC/bin/scripts/atacQC.R ${gAlias[${DIR}]}
    # ${gAlias[${DIR}]}
    /workdir/tools/atacQC/bin/scripts/html.atacQC.sh `echo ${PIN}_atacQC`

    cd ../../..

    /home/fa286/bin/tree-1.7.0/tree > folder.structure

}

# Delete files we no longer need
clean_up(){

    echo
    echo "Running module clean_up()"
    echo

    cd "${FOLDER}"/fastqs/dedup_BAMS
    rm *_DEDUP_unsorted.bam

    # if [ -f R2_info.txt ]; then
    #     rm -f R2_info.txt
    # fi

}


while getopts "hp:t:g:q:d:f:s:a:c:u:e:" opt; do
    case ${opt} in

    h)
        echo
        echo
        echo
        usage
        echo
        echo
        exit 1

    ;;

    f )

        FOLDER=$OPTARG
        echo
        echo "Project folder = " $FOLDER
        echo
    ;;

    p )

        PIN=$OPTARG
        echo "Project Identifier = " $PIN
    ;;

    t )

        T=$OPTARG

    ;;

    g)

        DIR=$OPTARG

    ;;

    q)

        QC=$OPTARG

    ;;

    d)
        DELIM=$OPTARG

    ;;

    s)
        SINGLE_STEP_LIST=$OPTARG

    ;;

    a)
        ADD_DUPLICATE_BC=$OPTARG

    ;;

    c)
        SAMPLE_CELL_NUMBERS=$OPTARG

    ;;

    u)
        PARAMETERS_FILE=$OPTARG

    ;;

    e)
        EMAIL_TO=$OPTARG

    ;;


    \? )
        echo
        echo
        echo
        usage

    ;;

    esac

done
# shift $((OPTIND -1))

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## Process parameters file, but don't overwrite parameters passed in from the command line
if [[ ! -z "${PARAMETERS_FILE}" ]]; then

    echo
    echo "Reading parameters file: " $PARAMETERS_FILE
    echo

    while IFS="=" read -r key value; do
        case "$key" in
            '#'*) ;;
            *)
                # echo $key "=" "${!key}"
                if [[ -z "${!key}" ]]; then
                    eval "$key=\"$value\""
                    echo $key "=" $value
                fi
            ;;
        esac
    done < "$PARAMETERS_FILE"
    echo
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if PIN is provided
if [[ -z "${PIN+x}" ]]; then

    PIN="PIN_Null"
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if FOLDER is provided
if [[ -z "${FOLDER}" ]]; then
    echo 
    echo "-f <folder name> is required."
    echo 
    echo 
    usage
    echo
    exit 1
fi

# PARAMETER CHECKS

                    #-------------------------------------------------------------------------------------------------------------
                    #-------------------------------------------------------------------------------------------------------------
                    ## check if delimiter parameter exists
                    if [[ ! -z "${DELIM+x}" ]]; then
                        #statements
                        if [[ $DELIM == default ]]; then

                        DELIMITER="_"
                        FIELD="5"
                        echo "file naming will be done using the default delimiter settings"
                      else

                        DELIMITER=`echo $DELIM | cut -d , -f1`
                        FIELD=`echo $DELIM | cut -d , -f2-`
                        echo "file naming will be done using the delim = $DELIMITER and field = $FIELD settings"

                      fi

                    fi

                    #-------------------------------------------------------------------------------------------------------------
                    #-------------------------------------------------------------------------------------------------------------
                    ## check if trimming parameter exists

                    if [[ ! -z "${T+x}" ]]; then

                        if [[ $T == nextseq ]]; then
                            trimPE
                        elif [[ $T == nova ]]; then
                            trimHiSeqPE
                        else
                        echo "-t only accepts nextseq or nova as arguments"
                        exit 1

                        fi
                    fi

                    ## check ADD_DUPLICATE_BC parameter

                    if [[ -z "${ADD_DUPLICATE_BC}" ]]; then
                        ADD_DUPLICATE_BC="NoAction"
                    else

                        if [[ $ADD_DUPLICATE_BC != "i5" && $ADD_DUPLICATE_BC != "i7i5" && $ADD_DUPLICATE_BC != "NoAction" ]]; then
                            echo "-a only accepts i5, i7i5, or NoAction as arguments"
                            exit 1
                        fi
                    fi

                    ## check / process SAMPLE_CELL_NUMBERS parameter
                    # No error checking: assumes correct format and no spaces
                    if [[ ! -z "${SAMPLE_CELL_NUMBERS}" ]]; then
                        for key_value in $(echo $SAMPLE_CELL_NUMBERS | sed "s/,/ /g")
                        do
                            # call your procedure/other scripts here below
                            key=`echo $key_value | cut -d ':' -f 1`
                            val=`echo $key_value | cut -d ':' -f 2`
                            sample_cell_number_dict[$key]=$val
                        done
                    fi

                    #-------------------------------------------------------------------------------------------------------------
                    #-------------------------------------------------------------------------------------------------------------
                    ## check if genomeDir provided

                    if [[ ! -z "${DIR+x}" ]]; then
                        if [ ${genomeDir[${DIR}]+_} ]; then
                            echo Reference genome selected = $DIR
                            echo
                            # If single step not set, execute all steps
                            # else execute only the step passed as a parameter
                            if [[ -z "${SINGLE_STEP_LIST}" ]]; then
                                findMESequence
                                horizontalMerge_padUMI_addDupBC
                                generate_whitelist
                                run_extract
                                cutadapt_trim5pME
                                cutadapt_trim3pAd
                                alignPE
                                sort
                                rmMT
                                # markDups
                                # dedupBAM
                                dedupBAM_with_UMI_tools
                                callPeak
                                mergedPeaks
                                saf
                                frip
                                tagDir
                                annotatePeaks
                                bedGraphs
                                clean_up
                            else
                                # No error checking: assumes correct format and no spaces
                                for SINGLE_STEP in $(echo $SINGLE_STEP_LIST | sed "s/,/ /g")
                                do
                                    # Note: I could just call the function from the parameter, but that would 
                                    # execute any string passed in by the user which cound be quite dangerous
                                    if [[ $SINGLE_STEP == findMESequence ]]; then
                                        findMESequence
                                    elif [[ $SINGLE_STEP == horizontalMerge_padUMI_addDupBC ]]; then
                                        horizontalMerge_padUMI_addDupBC
                                    elif [[ $SINGLE_STEP == generate_whitelist ]]; then
                                        generate_whitelist
                                    elif [[ $SINGLE_STEP == run_extract ]]; then
                                        run_extract
                                    elif [[ $SINGLE_STEP == cutadapt_trim5pME ]]; then
                                        cutadapt_trim5pME
                                    elif [[ $SINGLE_STEP == cutadapt_trim3pAd ]]; then
                                        cutadapt_trim3pAd
                                    elif [[ $SINGLE_STEP == alignPE ]]; then
                                        alignPE
                                    elif [[ $SINGLE_STEP == sort ]]; then
                                        sort
                                    elif [[ $SINGLE_STEP == rmMT ]]; then
                                        rmMT
                                    elif [[ $SINGLE_STEP == dedupBAM_with_UMI_tools ]]; then
                                        dedupBAM_with_UMI_tools
                                    elif [[ $SINGLE_STEP == callPeak ]]; then
                                        callPeak
                                    elif [[ $SINGLE_STEP == mergedPeaks ]]; then
                                        mergedPeaks
                                    elif [[ $SINGLE_STEP == saf ]]; then
                                        saf
                                    elif [[ $SINGLE_STEP == frip ]]; then
                                        frip
                                    elif [[ $SINGLE_STEP == tagDir ]]; then
                                        tagDir
                                    elif [[ $SINGLE_STEP == annotatePeaks ]]; then
                                        annotatePeaks
                                    elif [[ $SINGLE_STEP == bedGraphs ]]; then
                                        bedGraphs
                                    elif [[ $SINGLE_STEP == atacQC ]]; then
                                        atacQC
                                    elif [[ $SINGLE_STEP == clean_up ]]; then
                                        clean_up
                                    elif [[ $SINGLE_STEP == no_step ]]; then
                                        # This is used when we are trimming (via the -t parameter) 
                                        # or running atacQC and don't want to perform any other steps
                                        exit 1
                                    else
                                        echo
                                        echo "Error: step '"$SINGLE_STEP"' is not defined"
                                        echo "Please check your -s parameter"
                                        echo
                                        exit 1
                                    fi
                                done
                            fi
                        else
                            echo "The reference genome provided '"$DIR"' is not available"
                            exit 1

                        fi
                    fi

                    if [[ ! -z "${QC+x}" ]]; then

                        if [[ $QC == yes ]]; then
                            atacQC
                        else
                            echo "-q option only accepts yes as an argument"
                            exit 1
                        fi
                    fi






#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------

if [[ -z $1 ]] || [[  $1 = "help"  ]] ; then
    #statements
    echo
    echo
    usage
    echo
    echo
    exit 1

else
    echo >> beta5.atac.log
    echo `date -u` >> beta5.atac.log
    echo "Project Identifier Specified = " $PIN >> beta5.atac.log
    echo "Reference Genome Specified   = " $DIR >> beta5.atac.log
    echo "Trimming                     = " $T >> beta5.atac.log
    echo >> beta5.atac.log

    echo "ENV INFO: " >> beta5.atac.log
    echo >> beta5.atac.log
    # echo "STAR version:" `~/bin/STAR-2.7.0e/bin/Linux_x86_64/STAR --version` >> beta5.atac.log
    # echo "multiqc version:" `~/miniconda2/envs/RSC/bin/multiqc --version` >> beta5.atac.log
    echo "samtools version:" `/programs/bin/samtools/samtools --version` >> beta5.atac.log
    echo "macs2 version: macs2 2.1.0.20150731 " >> beta5.atac.log
    echo "HOMER version: v4.11.1" >> beta5.atac.log
    echo -------------------------------------------------------------------------------------------------- >> beta5.atac.log

fi
