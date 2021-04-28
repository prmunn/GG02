#!/usr/bin/bash
set -x #echo on

umi_tools whitelist --knee-method=density \
--set-cell-number 2000 \
--method=reads \
--plot-prefix BPA1_predictBC \
--allow-threshold-error \
--extract-method string \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
-L BPA1_predictedBCwhitelist.log \
-I BPA1_I2_I1_padUMI_R2.fastq.gz \
-S BPA1_predictedBCwhitelist.txt

umi_tools whitelist --knee-method=density \
--set-cell-number 2200 \
--method=reads \
--plot-prefix BPA2_predictBC \
--allow-threshold-error \
--extract-method string \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
-L BPA2_predictedBCwhitelist.log \
-I BPA2_I2_I1_padUMI_R2.fastq.gz \
-S BPA2_predictedBCwhitelist.txt

umi_tools whitelist --knee-method=density \
--set-cell-number 2600 \
--method=reads \
--plot-prefix BPA3_predictBC \
--allow-threshold-error \
--extract-method string \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
-L BPA3_predictedBCwhitelist.log \
-I BPA3_I2_I1_padUMI_R2.fastq.gz \
-S BPA3_predictedBCwhitelist.txt

umi_tools whitelist --knee-method=density \
--set-cell-number 1950 \
--method=reads \
--plot-prefix BPA4_predictBC \
--allow-threshold-error \
--extract-method string \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
-L BPA4_predictedBCwhitelist.log \
-I BPA4_I2_I1_padUMI_R2.fastq.gz \
-S BPA4_predictedBCwhitelist.txt

umi_tools whitelist --knee-method=density \
--set-cell-number 2200 \
--method=reads \
--plot-prefix Ctrl1_predictBC \
--allow-threshold-error \
--extract-method string \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
-L Ctrl1_predictedBCwhitelist.log \
-I Ctrl1_I2_I1_padUMI_R2.fastq.gz \
-S Ctrl1_predictedBCwhitelist.txt

umi_tools whitelist --knee-method=density \
--set-cell-number 2600 \
--method=reads \
--plot-prefix Ctrl2_predictBC \
--allow-threshold-error \
--extract-method string \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
-L Ctrl2_predictedBCwhitelist.log \
-I Ctrl2_I2_I1_padUMI_R2.fastq.gz \
-S Ctrl2_predictedBCwhitelist.txt

umi_tools whitelist --knee-method=density \
--set-cell-number 2100 \
--method=reads \
--plot-prefix Ctrl3_predictBC \
--allow-threshold-error \
--extract-method string \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
-L Ctrl3_predictedBCwhitelist.log \
-I Ctrl3_I2_I1_padUMI_R2.fastq.gz \
-S Ctrl3_predictedBCwhitelist.txt

umi_tools whitelist --knee-method=density \
--set-cell-number 2500 \
--method=reads \
--plot-prefix Ctrl4_predictBC \
--allow-threshold-error \
--extract-method string \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
-L Ctrl4_predictedBCwhitelist.log \
-I Ctrl4_I2_I1_padUMI_R2.fastq.gz \
-S Ctrl4_predictedBCwhitelist.txt

umi_tools whitelist --knee-method=density \
--set-cell-number 2900 \
--method=reads \
--plot-prefix Mix1_predictBC \
--allow-threshold-error \
--extract-method string \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
-L Mix1_predictedBCwhitelist.log \
-I Mix1_I2_I1_padUMI_R2.fastq.gz \
-S Mix1_predictedBCwhitelist.txt

umi_tools whitelist --knee-method=density \
--set-cell-number 1700 \
--method=reads \
--plot-prefix Mix2_predictBC \
--allow-threshold-error \
--extract-method string \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
-L Mix2_predictedBCwhitelist.log \
-I Mix2_I2_I1_padUMI_R2.fastq.gz \
-S Mix2_predictedBCwhitelist.txt

umi_tools whitelist --knee-method=density \
--set-cell-number 2200 \
--method=reads \
--plot-prefix Mix3_predictBC \
--allow-threshold-error \
--extract-method string \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
-L Mix3_predictedBCwhitelist.log \
-I Mix3_I2_I1_padUMI_R2.fastq.gz \
-S Mix3_predictedBCwhitelist.txt

umi_tools whitelist --knee-method=density \
--set-cell-number 1950 \
--method=reads \
--plot-prefix Mix4_predictBC \
--allow-threshold-error \
--extract-method string \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
-L Mix4_predictedBCwhitelist.log \
-I Mix4_I2_I1_padUMI_R2.fastq.gz \
-S Mix4_predictedBCwhitelist.txt

umi_tools whitelist --knee-method=density \
--set-cell-number 750 \
--method=reads \
--plot-prefix unknown_predictBC \
--allow-threshold-error \
--extract-method string \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNCCCCCCCCCC \
-L unknown_predictedBCwhitelist.log \
-I unknown_I2_I1_padUMI_R2.fastq.gz \
-S unknown_predictedBCwhitelist.txt

echo "Done"
