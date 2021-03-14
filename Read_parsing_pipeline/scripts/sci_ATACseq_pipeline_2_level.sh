#!/usr/bin/bash
#SBATCH -J ATACseq
#SBATCH -o %x.out
#SBATCH -n 12
#SBATCH --mem-per-cpu=18000

source ~/.bash_profile

usage(){

    echo "s c i - A T A C - S E Q   P I P E L I N E"
    echo
    echo

    echo "Usage: bash" $0 "-f arg [-h arg] [-p arg] [-d args] [-t arg] [-g arg] [-q arg] [-s arg]"
    echo
    echo "---------------------------------------------------------------------------------------------------------------"
    echo " -f  --> Folder containing data files (required)"
    echo "[-h] --> Display Help"
    echo "[-p] --> Project Identifier Number"
    echo "[-d] --> Comma Spearated Values for Delimiter and Field <delim,field or default> default: _,5 "
    echo "[-t] --> Trimming <nextseq or nova>;"
    echo "[-g] --> Reference Genome <mm10 or hg38>"
    echo "[-q] --> Execute atacQC.R script <yes>"
    echo "[-s] --> Execute single step"
    echo "---------------------------------------------------------------------------------------------------------------"
}



declare -A genomeDir

genomeDir=(
["mm10"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/BWAIndex/genome.fa" \
["hg38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/bwa.index/Homo_sapiens.GRCh38.dna.toplevel.fa"
)

declare -A gtfs

gtfs=(
["mm10"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/Mus_musculus.GRCm38.96.gtf" \
["hg38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/Homo_sapiens.GRCh38.96.gtf"
)

declare -A gAlias # for compatibility with atacQC.R
gAlias=(
["mm10"]="mouse" \
["hg38"]="human"
)

# PAD_DIR="padUMI-fastq"
# I7PCR_DEMUX_FASTQ="i7PCR_demux_fastq"
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
    echo "findMESequence()"
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
horizontalMerge_padUMI(){

    cd "${FOLDER}"/fastqs/parsed_BCs/

    RELATIVE_ORIG_FASTQ=../"${ORIG_FASTQ}"

    # Set file names
    input_R2_fastq_file=("${RELATIVE_ORIG_FASTQ}"/*_R2_001.fastq.gz)
    input_I1_fastq_file=("${RELATIVE_ORIG_FASTQ}"/*_I1_001.fastq.gz)
    input_I2_fastq_file=("${RELATIVE_ORIG_FASTQ}"/*_I2.f*.gz)
    file_prefix=`echo $(basename $input_R2_fastq_file .gz) | cut -d "_" -f 1`

    echo
    echo "horizontalMerge_padUMI()"
    echo "file_prefix: "$file_prefix
    echo "input_R2_fastq_file: "$input_R2_fastq_file
    echo "input_I1_fastq_file: "$input_I1_fastq_file
    echo "input_I2_fastq_file: "$input_I2_fastq_file
    echo
        
    ../../../horizontalMerge-padUMI.awk R2_info.txt \
    <(paste <(zcat $input_I2_fastq_file) <(zcat $input_I1_fastq_file) <(zcat $input_R2_fastq_file)) \
    > "${file_prefix}"_I2_I1_padUMI_R2.fastq

    # Clean up
    # none

    cd ../../..

}

# Step 3a) Optional: generate combBC whitelist with UMI_tools whitelist
generate_whitelist() {

    # Assumes that environment variables already set up for UMI_tools
    # export PYTHONPATH=/programs/UMI-tools/lib/python3.6/site-packages:/programs/UMI-tools/lib64/python3.6/site-packages
    # export PATH=/programs/UMI-tools/bin:$PATH

    cd "${FOLDER}"/fastqs/parsed_BCs/

    # Set file names
    I2_I1_padUMI_R2_fastq_file=(*I2_I1_padUMI_R2.fastq)
    file_prefix=`echo $(basename $I2_I1_padUMI_R2_fastq_file .gz) | cut -d "_" -f 1`

    echo
    echo "generate_whitelist()"
    echo "file_prefix: "$file_prefix
    echo "I2_I1_padUMI_R2_fastq_file: "$I2_I1_padUMI_R2_fastq_file
    echo

    # Run umi_tools whitelist to generate cell barcode whitelist (specific to dataset). 
    # This step may require testing different options to get the optimal whitelist 
    # (e.g. --error-correct-threshold, --knee-method=[density/distance], --method=[reads/umis]).
    umi_tools whitelist --knee-method=density \
    --method=reads \
    --plot-prefix "${file_prefix}"_predictBC \
    --allow-threshold-error \
    --extract-method string \
    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
    --error-correct-threshold=2 \
    --ed-above-threshold=correct \
    -L "${file_prefix}"_predictedBCwhitelist.log \
    -I $I2_I1_padUMI_R2_fastq_file \
    -S "${file_prefix}"_predictedBCwhitelist.txt

    cd ../../..

}

# Step 3b) Run umi_tools extract to move cell barcode and UMI to R1/R2 headers. 
# R2 is the ‘driver’ because it contains the BC and UMI sequences
run_extract() {

    # Assumes that environment variables already set up for UMI_tools
    # export PYTHONPATH=/programs/UMI-tools/lib/python3.6/site-packages:/programs/UMI-tools/lib64/python3.6/site-packages
    # export PATH=/programs/UMI-tools/bin:$PATH

    cd "${FOLDER}"/fastqs/parsed_BCs/

    RELATIVE_ORIG_FASTQ=../"${ORIG_FASTQ}"

    # Set file names
    input_R1_fastq_file=("${RELATIVE_ORIG_FASTQ}"/*_R1_001.fastq.gz)
    I2_I1_padUMI_R2_fastq_file=(*I2_I1_padUMI_R2.fastq)
    predictedBCwhitelist_file=(*_predictedBCwhitelist.txt)
    file_prefix=`echo $(basename $input_R1_fastq_file .gz) | cut -d "_" -f 1`

    echo
    echo "run_extract()"
    echo "file_prefix: "$file_prefix
    echo "input_R1_fastq_file: "$input_R1_fastq_file
    echo "I2_I1_padUMI_R2_fastq_file: "$I2_I1_padUMI_R2_fastq_file
    echo "predictedBCwhitelist_file: "$predictedBCwhitelist_file
    echo

    umi_tools extract --extract-method=string \
    --quality-filter-mask 20 \
    --quality-encoding phred33 \
    -p CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
    --whitelist $predictedBCwhitelist_file \
    --error-correct-cell \
    -I $I2_I1_padUMI_R2_fastq_file \
    -S "${file_prefix}"_hBC_UMI_R2.fastq.gz \
    --read2-in=$input_R1_fastq_file \
    --read2-out="${file_prefix}"_hBC_UMI_R1.fastq.gz \
    -L "${file_prefix}"_extractBC.log

    cd ../../..

}

# Step 4) Use cutadapt to trim (and require) 5' anchored ME in R2; write untrimmed paired reads to separate 
# file (do not use). Again, R2 is the ‘driver’ because it should now contain the ME sequence at the 5’ end
cutadapt_trim5pME(){

    cd "${FOLDER}"/fastqs/parsed_BCs/

    # Delete output files from any prior runs
    rm -f *_noME_hBC_UMI_*
    rm -f *_hasME_hBC_UMI_*

    # Set file names
    hBC_UMI_R2_fastq_file=(*_hBC_UMI_R2.fastq.gz)
    hBC_UMI_R1_fastq_file=(*_hBC_UMI_R1.fastq.gz)
    file_prefix=`echo $(basename $hBC_UMI_R2_fastq_file .gz) | cut -d "_" -f 1`

    echo
    echo "cutadapt_trim5pME()"
    echo "file_prefix: "$file_prefix
    echo "hBC_UMI_R2_fastq_file: "$hBC_UMI_R2_fastq_file
    echo "hBC_UMI_R1_fastq_file: "$hBC_UMI_R1_fastq_file
    echo

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
    > "${file_prefix}"_trim5pME_log.out

    cd ../../..

}


# Step 5) Trim any internal 3' ME-REV in R1 or R2 due to small inserts, trim for quality, filter for short reads (standard cutadapt).
# Adjust 3’-quality trimming mode (-q/--nextseq-trim) based on sequencer (4-color/2-color).
# R2 is listed first because that has been the pattern in this pipeline, but the order of R1/R2 doesn’t matter (R1/R2 adaptors are the same).
cutadapt_trim3pAd(){

    cd "${FOLDER}"/fastqs/parsed_BCs/

    # Set file names
    hasME_hBC_UMI_R2_fastq_file=(*hasME_hBC_UMI_R2.fastq.gz)
    hasME_hBC_UMI_R1_fastq_file=(*hasME_hBC_UMI_R1.fastq.gz)
    file_prefix=`echo $(basename $hasME_hBC_UMI_R2_fastq_file .gz) | cut -d "_" -f 1`

    echo
    echo "cutadapt_trim3pAd()"
    echo "file_prefix: "$file_prefix
    echo "hasME_hBC_UMI_R2_fastq_file: "$hasME_hBC_UMI_R2_fastq_file
    echo "hasME_hBC_UMI_R1_fastq_file: "$hasME_hBC_UMI_R1_fastq_file
    echo

    # Don't need this directory if not using trim_galore - talk to Faraz about this
    mkdir fastQC

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

    #for i in "${fastqs[@]}"
    #do
    #    trim_galore --quality 20 --gzip --length 20  --paired --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
    #done

    mkdir TrimQC_stats trimmed_fastqs
    #mv *_trimming_report.txt TrimQC_stats
    mv *_val* trimmed_fastqs
    mv TrimQC_stats fastQC trimmed_fastqs ..

    cd ../../..

}


trimPE(){

    cd fastqs
    ls -1 *_R1.fastq* > .R1
    ls -1 *_R2.fastq* > .R2
    paste -d " " .R1 .R2 > Reads.list

    readarray fastqs < Reads.list
    mkdir fastQC

    for i in "${fastqs[@]}"
    do
        $TRIM --nextseq 20 --length 20  --paired --gzip --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
    done

    mkdir TrimQC_stats trimmed_fastqs
    mv *_trimming_report.txt TrimQC_stats
    mv *_val* trimmed_fastqs
    mv TrimQC_stats fastQC trimmed_fastqs ../

    cd ..
}


trimHiSeqPE(){

    cd fastqs
    ls -1 *_1.fq* > .R1
    ls -1 *_2.fq* > .R2
    paste -d " " .R1 .R2 > Reads.list

    readarray fastqs < Reads.list
    mkdir fastQC

    for i in "${fastqs[@]}"
    do
        trim_galore --quality 20 --gzip --length 20  --paired --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
    done

    mkdir TrimQC_stats trimmed_fastqs
    mv *_trimming_report.txt TrimQC_stats
    mv *_val* trimmed_fastqs
    mv TrimQC_stats fastQC trimmed_fastqs ..

    cd ..
}


# xxx
#//--- Need the mm10 / hg38 genome
alignPE(){

    cd "${FOLDER}"/fastqs/trimmed_fastqs/

    ls -1 *_val_1.fq.gz > .trR1
    ls -1 *_val_2.fq.gz > .trR2
    paste -d " " .trR1 .trR2 > Trimmed.list

    readarray trimmedFastqs < Trimmed.list

    for i in "${trimmedFastqs[@]}"

    do
        # INDEX="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/BWAIndex/genome.fa"
        iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`

        bwa mem -t 24 -M -R "@RG\tID:${iSUB}\tSM:${iSUB}\tPL:ILLUMINA\tLB:${iSUB}\tPU:1" ${genomeDir[${DIR}]} $i \
        | samtools view -@ 24 -b -h -F 0x0100 -O BAM -o ${iSUB}.bam
    done

    mkdir primary_BAMS
    mv *.bam primary_BAMS
    mv primary-BAMS ..
    cd ../../..

}

sort(){

    cd "${FOLDER}"/fastqs/primary_BAMS

    for i in *.bam
    do
        samtools sort $i > `echo  $i | cut -d "." -f1`.sorted.bam
    done

    for i in *.sorted.bam
    do
        samtools index $i
    done

    # alignment stats etc. on raw bams
    for i in *.sorted.bam
    do
        iSUB=`echo $i | cut -d "." -f1`
        samtools flagstat $i > ${iSUB}.primary.flagstat
        samtools idxstats $i > ${iSUB}.primary.idxstats
    done

    cd ../../..

}

rmMT(){

    cd "${FOLDER}"/fastqs/primary_BAMS

    for i in *.sorted.bam
    do

        iSUB=`echo $i | cut -d "." -f1`

        samtools view -H `ls -1 *.sorted.bam | head -1` | cut -f2 | grep "SN:" |  cut -d ":" -f2 | grep -v "MT\|_\|\." | xargs samtools view -b $i > ${iSUB}.noMT.bam

    done
    cd ../../..
}

markDups(){

    cd "${FOLDER}"/fastqs/primary_BAMS

    for i in *.noMT.bam
    do
        iSUB=`echo $i | cut -d "." -f1`
        java -jar /programs/bin/picard-tools/picard.jar \
        MarkDuplicates \
        INPUT=$i \
        OUTPUT=${iSUB}.dupMarked.noMT.bam \
        ASSUME_SORTED=true \
        REMOVE_DUPLICATES=false \
        METRICS_FILE=${iSUB}.MarkDuplicates.metrics.txt \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=tmp

    done
    cd ../../..
}

dedupBAM(){

    cd "${FOLDER}"/fastqs/primary_BAMS

    # alignment stats etc. on dupMarked no MT bams
    for i in *.dupMarked.noMT.bam
    do
        iSUB=`echo $i | cut -d "." -f1`
        samtools index $i
        samtools flagstat $i > ${iSUB}.noMT.flagstat
        samtools idxstats $i > ${iSUB}.noMT.idxstats
    done

    for i in *.dupMarked.noMT.bam
    do
        iSUB=`echo $i | cut -d "." -f1`
        samtools view -b -h -F 0X400 $i > ${iSUB}.DEDUP.bam
    done

    for i in *.DEDUP.bam; do samtools index $i ; samtools idxstats $i > `echo $i | cut -d "." -f1`.DEDUP.idxstats; done
    for i in *.DEDUP.bam; do samtools flagstat $i > `echo $i | cut -d "." -f1`.DEDUP.flagstat; done

    multiqc -n ${PIN}.multiqc.report .

    mkdir dedup_BAMS
    mv *.DEDUP* dedup_BAMS/
    mv dedup-BAMS ..
    cd ../../..

}

tagDir(){
    cd "${FOLDER}"/fastqs/dedup_BAMS

    for i in *.DEDUP.bam
    do
    iSUB=`echo "$i" | cut -d'.' -f1` # subset to rename
    /home/fa286/bin/HOMER/bin/makeTagDirectory "$iSUB".tag.dir "$i"
    done
    cd ../../..
}

callPeak(){

    cd "${FOLDER}"/fastqs/dedup_BAMS

    echo "calling peaks on DEDUP bams"
    mkdir peaks.OUT
    for i  in *.DEDUP.bam
    do
        iSUB=`echo $i | cut -d "." -f1`
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

    echo "running module mergedPeaks"
    cd "${FOLDER}"/fastqs/dedup_BAMS

    allBams=`echo *.DEDUP.bam`

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

    awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' allSamplesMergedPeakset_peaks.narrowPeak > allSamplesMergedPeakset.saf
    cd ../../../..
}

#saf

frip(){
    # featureCounts -p -a ${sample}_peaks.saf -F SAF -o readCountInPeaks.txt ${sample}.sorted.marked.filtered.shifted.bam
    cd "${FOLDER}"/fastqs/dedup_BAMS/

    for i  in *.DEDUP.bam
    do
        iSUB=`echo $i | cut -d "." -f1`
        featureCounts -p -a peaks.OUT/allSamplesMergedPeakset.saf -F SAF -o "${iSUB}".readCountInPeaks.txt $i
    done

    cd ../../..
}

annotatePeaks(){

    cd "${FOLDER}"/fastqs/dedup_BAMS/peaks.OUT
    /home/fa286/bin/HOMER/bin/annotatePeaks.pl allSamplesMergedPeakset.saf ${DIR} -gtf ${gtfs[${DIR}]} > allSamplesMergedPeakset.Annotated.saf
    cd ../../../..
}

bedGraphs(){

    cd "${FOLDER}"/fastqs/dedup_BAMS

    for i in *.tag.dir
    do
        makeUCSCfile ${i} -o auto -fsize 1e10 -res 1 -color 106,42,73 -style chipseq
    done

    mkdir tagDirs
    mv *.tag.dir tagDirs
    cd tagDirs
    mkdir bedGraphs

    for i in *.tag.dir
    do
        cd $i
        zcat *.ucsc.bedGraph.gz | awk '{if(NR>1) print "chr"$0; else print $0}' | gzip > `basename *.ucsc.bedGraph.gz .ucsc.bedGraph.gz`.ucsc.bg.gz
        mv *.ucsc.bg.gz ../bedGraphs
        cd ..
    done
    cd ..

    mkdir featureCounts
    mv *.txt featureCounts

    multiqc -n ${PIN}.FRIP.multiqc.report -b "Please note that the featureCounts M Assigned Column refers to Fragments and Not Reads" --ignore tagDirs --ignore peaks.OUT .

    cd ../../..
}

atacQC(){

    cd "${FOLDER}"/fastqs/dedup_BAMS

    echo "genome alias" = ${gAlias[${DIR}]}
    /programs/R-3.6.3/bin/Rscript /home/fa286/bin/scripts/atacQC.R ${gAlias[${DIR}]}
    # ${gAlias[${DIR}]}
    ~/bin/scripts/html.atacQC.sh `echo ${PIN}_atacQC`

    cd ../../..

    /home/fa286/bin/tree-1.7.0/tree > folder.structure

}


while getopts "hp:t:g:q:d:f:s:" opt; do
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
        SINGLE_STEP=$OPTARG

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

                    #-------------------------------------------------------------------------------------------------------------
                    #-------------------------------------------------------------------------------------------------------------
                    ## check if genomeDir provided

                    if [[ ! -z "${DIR+x}" ]]; then
                        if [ ${genomeDir[${DIR}]+_} ]; then
                            echo Reference genome selected = $DIR
                            echo
                            # If single step not set, execute all steps
                            # else execute only the step passed as a parameter
                            if [[ -z "${SINGLE_STEP}" ]]; then
                                findMESequence
                                horizontalMerge_padUMI
                                generate_whitelist
                                run_extract
                                cutadapt_trim5pME
                                cutadapt_trim3pAd
                                alignPE
                                sort
                                rmMT
                                markDups
                                dedupBAM
                                callPeak
                                mergedPeaks
                                saf
                                frip
                                tagDir
                                annotatePeaks
                                bedGraphs
                            else
                                # Note: I could just call the function from the parameter, but that would 
                                # execute any string passed in by the user which cound be quite dangerous
                                if [[ $SINGLE_STEP == findMESequence ]]; then
                                    findMESequence
                                elif [[ $SINGLE_STEP == horizontalMerge_padUMI ]]; then
                                    horizontalMerge_padUMI
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
                                elif [[ $SINGLE_STEP == markDups ]]; then
                                    markDups
                                elif [[ $SINGLE_STEP == dedupBAM ]]; then
                                    dedupBAM
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
    echo "HOMER version: v4.10.4" >> beta5.atac.log
    echo -------------------------------------------------------------------------------------------------- >> beta5.atac.log

fi
