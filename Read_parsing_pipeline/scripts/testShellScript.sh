#!/usr/bin/bash
#SBATCH -J ATACseq
#SBATCH -o %x.out
#SBATCH -n 12
#SBATCH --mem-per-cpu=18000

source ~/.bash_profile

cd GIA0130-AdrianBarnyard/fastqs/parsed_BCs/

echo
echo "Running module testShellScript.sh"
echo

ADD_DUPLICATE_BC="NoAction"

../../../testAwkScript.awk R2_info.txt \
<(paste <(zcat ../orig-fastqs/CKDL200168997-1a_L6_I2.fq.gz) <(zcat ../orig-fastqs/../orig-fastqs/GIA0130_CKDL200168997-1a_HFJ5TBBXX_S1_L006_I1_001.fastq.gz) <(zcat ../orig-fastqs/GIA0130_CKDL200168997-1a_HFJ5TBBXX_S1_L006_R2_001.fastq.gz))

#awk -v addDuplicateBC="$ADD_DUPLICATE_BC" -f ../../../testAwkScript.awk R2_info.txt \
#<(paste <(zcat ../orig-fastqs/CKDL200168997-1a_L6_I2.fq.gz) <(zcat ../orig-fastqs/../orig-fastqs/GIA0130_CKDL200168997-1a_HFJ5TBBXX_S1_L006_I1_001.fastq.gz) <(zcat ../orig-fastqs/GIA0130_CKDL200168997-1a_HFJ5TBBXX_S1_L006_R2_001.fastq.gz))

#../../../testAwkScript.awk R2_info.txt \
#<(paste <(zcat ../orig-fastqs/CKDL200168997-1a_L6_I2.fq.gz) <(zcat ../orig-fastqs/../orig-fastqs/GIA0130_CKDL200168997-1a_HFJ5TBBXX_S1_L006_I1_001.fastq.gz) <(zcat ../orig-fastqs/GIA0130_CKDL200168997-1a_HFJ5TBBXX_S1_L006_R2_001.fastq.gz)) \
#$ADD_DUPLICATE_BC

