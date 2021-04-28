#!/usr/bin/bash
#SBATCH -J ATACseq
#SBATCH -o %x.out
#SBATCH -n 12
#SBATCH --mem-per-cpu=18000

source ~/.bash_profile
source ../sciATAC-2level/sci_ATACseq_pipeline_2_level.sh

usage(){

    echo "U H T - A T A C - S E Q   P I P E L I N E"
    echo
    echo

    echo "Usage: bash" $0 "-f arg -g arg [-h arg] [-p arg] [-d args] [-t arg] [-q arg] [-s arg] [-a arg] [-c arg] [-u arg]"
    echo
    echo "Example: bash" $0 "-f Novogene-G201123J-noDemux/ -g hg38 -s findMESequence"
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


PAD_DIR="padUMI-fastq"
I7PCR_DEMUX_FASTQ="i7PCR_demux_fastq"

# Step 1) Use cutadapt to find the position of the ME sequence in the R2 read 
# (no modification of reads, just generate info file with ME position)
findMESequence(){

    cd fastqs/parsed-BCs/${PAD_DIR}

    RELATIVE_ORIG_FASTQ=../../"${ORIG_FASTQ}"

    # Delete summary files from any prior runs
    # rm *_R2_info_summary.txt

    for i in *_R2_info_j10.txt
    do
        iSUB=`echo $i | cut -d "_" -f 1`
        ls -1 "${RELATIVE_ORIG_FASTQ}"/"${iSUB}"*_R2_001.fastq.gz > input_fastq
        readarray input_fastq_file_list < input_fastq
        input_fastq_file=`echo $input_fastq_file_list | cut -d " " -f 9`

        echo
        echo "findMESequence()"
        echo "i: "$i
        # echo "iSUB: "$iSUB
        # echo "input_fastq_file: "$input_fastq_file
        echo

        # Start with parsing ME seq in R2 and padding UMI using cutadapt (with no reads output file, only info file)
        cutadapt -g AGATGTGTATAAGAGACAG \
        --info-file $i \
        -e 0.11 \
        -j 10 \
        --match-read-wildcards \
        --action none \
        -o /dev/null \
        $input_fastq_file

        # Summary counts for columns 2, 3, and 4
        awk -F '\t' '{a[$2]++; if($2>=0) {b[$3]++; c[$4]++;} else next;} END {print "col2=mismatch"; for(i in a) print a[i],i; print "\ncol3=startpost"; for(j in b) print b[j],j; print "\ncol4=endpos"; for(k in c) print c[k],k;}' \
        $i > "${iSUB}"_R2_info_summary.txt
    done

    # Clean up
    rm input_fastq

    cd ../../..

}

# Step 2) Custom shell script reads in the info file from step 1 and reconfigures R2 
# to build the combinatorial BC (integrate I1 and I2 reads) and ‘pad out’ the UMI
parseVarUMI(){

    cd fastqs/parsed-BCs/${PAD_DIR}
    
    RELATIVE_ORIG_FASTQ=../../"${ORIG_FASTQ}"

    for i in *.txt.gz
    do
        iSUB=`echo $i | cut -d "_" -f 1`
        ls -1 "${RELATIVE_ORIG_FASTQ}"/"${iSUB}"*_R2_001.fastq.gz > input_fastq
        readarray input_fastq_file_list < input_fastq
        input_fastq_file=`echo $input_fastq_file_list | cut -d " " -f 9`

        echo
        echo "parseVarUMI()"
        echo "i: "$i
        echo "iSUB: "$iSUB
        echo "input_fastq_file: "$input_fastq_file
        echo
        
        #../parse_varUMI.awk <(zcat $i) <(zcat $input_fastq_file) | gzip > ../"${iSUB}"_padUMI_R2.fastq.gz

    done

    cd ../../..

}

# Step 3a) Optional: generate combBC whitelist with UMI_tools whitelist
predictI7TagBCs() {

    # Assumes that environment variables already set up for UMI_tools
    # export PYTHONPATH=/programs/UMI-tools/lib/python3.6/site-packages:/programs/UMI-tools/lib64/python3.6/site-packages
    # export PATH=/programs/UMI-tools/bin:$PATH

    cd fastqs/parsed-BCs/

    umi_tools whitelist --bc-pattern=NNNNNNNNNNCCCCCCCCCC \
    --error-correct-threshold=1 \
    --ed-above-threshold=correct \
    --extract-method string \
    -L G201123J_padUMI_predictBCwhitelist.txt \
    -I G201123J_padUMI_R2.fastq.gz \
    -S G201123J_predictedBCwhitelist.txt

    cd ../..

}

# Step 3b) Use UMI_tools extract to reconfigure the R1/R2 headers to contain UMI+combBC.
moveUMItoHeader() {

    # Assumes that environment variables already set up for UMI_tools
    # export PYTHONPATH=/programs/UMI-tools/lib/python3.6/site-packages:/programs/UMI-tools/lib64/python3.6/site-packages
    # export PATH=/programs/UMI-tools/bin:$PATH

    cd fastqs/parsed-BCs/${PAD_DIR}

    umi_tools extract --extract-method=string \
    --quality-filter-mask 20 \
    --quality-encoding phred33 \
    -p NNNNNNNNNN \
    -I G201123J_padUMI_R2.fastq.gz \
    -S G201123J_UMI_R2.fq.gz \
    --read2-in=../../orig-fastq/G201123J_CKDL200167681-1a_HF5L2CCX2_S1_L005_R1_001.fastq.gz \
    --read2-out=G201123J_UMI_R1.fq.gz

    cd ../../..

}

# Step 4) Use cutadapt to trim (and require) 5’ME in R2. Reads that don’t have the ME seq in the right place 
# (after parsing/moving UMI+combBC) are written to a separate output file (maintaining R1/R2 pairing).

# Step 5) Use cutadapt to trim R1/R2 for 3’ME-rev (e.g. small inserts) and for 3’quality 
# (standard 3’adaptor/quality trimming for PE reads). Then move to mapping, peak calling steps.










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



markDups(){
                cd primary-BAMS
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
                cd ..
}

dedupBAM(){
                cd primary-BAMS
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

                mkdir dedup-BAMS
                mv *.DEDUP* dedup-BAMS/
                mv dedup-BAMS ..
                cd ..

}



while getopts "hp:t:g:q:d:f:s:a:c:u:" opt; do
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
                                parseVarUMI
                                predictI7TagBCs
                                moveUMItoHeader
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
                            else
                                # No error checking: assumes correct format and no spaces
                                for SINGLE_STEP in $(echo $SINGLE_STEP_LIST | sed "s/,/ /g")
                                do
                                    # Note: I could just call the function from the parameter, but that would 
                                    # execute any string passed in by the user which cound be quite dangerous
                                    if [[ $SINGLE_STEP == findMESequence ]]; then
                                        findMESequence
                                    elif [[ $SINGLE_STEP == parseVarUMI ]]; then
                                        parseVarUMI
                                    elif [[ $SINGLE_STEP == predictI7TagBCs ]]; then
                                        predictI7TagBCs
                                    elif [[ $SINGLE_STEP == moveUMItoHeader ]]; then
                                        moveUMItoHeader
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
