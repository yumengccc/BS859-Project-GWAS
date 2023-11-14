#!/bin/bash
module load R
module load python2
module load ldsc
module load prsice
module load gcta/1.94.1
module load metal

export T1D='/projectnb/bs859/data/T1D'
export RA_UKBB='/projectnb/bs859/data/RheumatoidArthritis/ukbb'
export RA_NARAC='/projectnb/bs859/data/RheumatoidArthritis/NARAC'
export LDSCORES='/projectnb/bs859/data/ldscore_files'
export final='/projectnb/bs859/students/yumeng88/Final'

wc /projectnb/bs859/data/T2D/T2D_Xue_et_al_2018.txt

#################################################################################### 1. LD score regression
## The UK Biobank-based LD scores using hapmap3 snps 
ls $LDSCORES/UKBB.ALL.ldscore

### RA UKBB meta-analysis (1,605 cases & 359,589 controls)
cat $RA_UKBB/README
zcat $RA_UKBB/M13_RHEUMA.gwas.imputed_v3.both_sexes.tsv.gz|head
zcat $RA_UKBB/variants.tsv.gz|head
zcat $RA_UKBB/M13_RHEUMA.gwas.imputed_v3.both_sexes.tsv.gz|wc
zcat $RA_UKBB/variants.tsv.gz|wc

### merge
zcat $RA_UKBB/M13_RHEUMA.gwas.imputed_v3.both_sexes.tsv.gz | awk -F'\t' -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' > RAuk_1.tsv
zcat $RA_UKBB/variants.tsv.gz > RAuk_2.tsv
join -t $'\t' -1 1 -2 1 RAuk_1.tsv RAuk_2.tsv > RAuk_merged.tsv
wc RAuk_1.tsv
wc RAuk_2.tsv
head RAuk_merged.tsv
awk -F'\t' -v OFS='\t' '{print $13, $14, $17, $25, $15, $26, $9, $10, $12, $28}' RAuk_merged.tsv > RAuk_merged_final.tsv
head RAuk_merged_final.tsv

### T1D meta-analysis (1,445 cases & 362,050 controls)
zcat $T1D/GCST90014023_buildGRCh38.tsv.gz|head
zcat $T1D/GCST90014023_buildGRCh38.tsv.gz|wc
zcat $T1D/GCST90014023_buildGRCh38.tsv.gz|awk 'NR==1{print "SNP A1 A2 freq beta se p N"};NR>1{print $1,$5,$6,$7,$8,$9,$2,$10}' > T1D_final.tsv
awk 'NR==1{print "SNP A1 A2 freq b se p N"};$7 == ($7+0)' T1D_final.tsv > T1D_final.txt
awk 'NR==1{print "SNP A1 A2 beta se p"};NR>1{print $1,$2,$3,$5,$6,$7}' T1D_final.tsv > T1D_final_2.tsv
awk 'NR==1 || $6 < 1e300' T1D_final_2.tsv > T1D_final_3.tsv 
awk 'NR==1{print "SNP A1 A2 beta se p"};$6 == ($6+0)' T1D_final_3.tsv > T1D_final_4.txt
wc T1D_final_3.tsv T1D_final_4.txt


##Step 1:  format RA & T1D summary statistics for ldsc
munge_sumstats.py \
    --sumstats $final/RAuk_merged_final.tsv \
    --snp rsid \
    --N-cas 1605 \
    --N-con 359589 \
    --a1 minor_allele \
    --a2 ref \
    --signed-sumstats beta,0 \
    --merge-alleles $LDSCORES/w_hm3.snplist \
--out RA

munge_sumstats.py \
    --sumstats $final/T1D_final_4.txt \
    --snp SNP \
    --N-cas 1445 \
    --N-con 362050 \
    --a1 A1 \
    --a2 A2 \
    --signed-sumstats beta,0 \
    --merge-alleles $LDSCORES/w_hm3.snplist \
--out T1D

##output:  reformated, gzipped summary statistics (in your own directory)
zcat  RA.sumstats.gz|head
zcat  RA.sumstats.gz|wc
zcat  T1D.sumstats.gz|head
zcat  T1D.sumstats.gz|wc

### using UKBB EUR ldscores
ldsc.py \
    --rg RA.sumstats.gz,T1D.sumstats.gz \
    --ref-ld $LDSCORES/UKBB.ALL.ldscore/UKBB.EUR.rsid \
    --w-ld $LDSCORES/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--out RA_T1D_rg


#################################################################################### 2. PRS analysis
Rscript $SCC_PRSICE_BIN/PRSice.R --dir . \
    --prsice $SCC_PRSICE_BIN/PRSice \
    --base $final/T1D_final.txt \
    --target $RA_NARAC/narac_clean_hg19 \
    --stat b \
    --snp SNP \
    --A1 A1 \
    --A2 A2 \
    --pvalue p \
    --binary-target T \
    --quantile 5 \
    --quant-ref 1 \
   --perm 1000 \
   --seed 542386 \
--out T1D_base_RA_target

ls T1D_base_RA_target.*
cat T1D_base_RA_target.summary
head T1D_base_RA_target.prsice
wc T1D_base_RA_target.prsice
head T1D_base_RA_target.best


#################################################################################### 3. Bi-directional Mendelian randomization
##reformat the summary statistics for GSMR:
zcat $T1D/GCST90014023_buildGRCh38.tsv.gz|awk 'NR==1{print "SNP A1 A2 freq b se p N"};NR>1{print $1,$5,$6,$7,$8,$9,$2,$10}' > T1D.ss.txt
head T1D.ss.txt
##remove lines have a duplicate first column
sort -t ' ' -k 1,1 -u T1D.ss.txt > T1D.ss.dedup.txt
head T1D.ss.dedup.txt
wc T1D.ss.txt T1D.ss.dedup.txt
mv T1D.ss.dedup.txt T1D.ss.txt


head RAuk_merged.tsv
awk 'NR==1{print "SNP A1 A2 freq b se p N"};NR>1{print $3,$4,$5,$6,$7,$8,$9,$10}' RAuk_merged_final.tsv > RA.ss.txt
head RA.ss.txt
##remove lines have a duplicate first column
sort -t ' ' -k 1,1 -u RA.ss.txt > RA.ss.dedup.txt
head RA.ss.dedup.txt
wc RA.ss.txt RA.ss.dedup.txt
mv RA.ss.dedup.txt RA.ss.txt

echo "RA RA.ss.txt" > exposure.txt
echo "T1D T1D.ss.txt" > outcome.txt

gcta64 --gsmr-file exposure.txt outcome.txt --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR --gsmr-direction 2 --out RA-T1D-gsmr --effect-plot
##What output files were created?
ls RA-T1D-gsmr.*
##check log file for errors first!
more RA-T1D-gsmr.log
##SNPs removed due to HEIDI outlier procedure (pleiotropic):
cat RA-T1D-gsmr.pleio_snps
##results:
cat RA-T1D-gsmr.gsmr


#################################################################################### 4. multi-trait meta-analysis \
#awk 'NR==1{print "SNP CHR BPOS A1 A2 freq beta se pval n"};NR>1{print $3,$1,$2,$4,$5,$6,$7,$8,$9,$10}' $final/RAuk_merged_final.tsv > $final/RA.meta.txt
#zcat $T1D/GCST90014023_buildGRCh38.tsv.gz|awk 'NR==1{print "SNP CHR BPOS A1 A2 freq beta se pval n"};NR>1{print $1,$3,$4,$5,$6,$7,$8,$9,$2,$10}' > $final/T1D.metapre.txt

#awk 'BEGIN {OFS="\t"} {if (NR == 1) print $0, "z"; else if (NR > 1) {if ($8 != 0) print $0, $7/$8; else print $0, "NA"}}' /projectnb/bs859/data/T2D/T2D_Xue_et_al_2018.txt > $final/T2D.metapre.txt
#awk 'NR==1{print "SNP CHR BPOS A1 A2 freq z pval n"};NR>1{print $3,$1,$2,$4,$5,$6,$11,$9,$10}' $final/T2D.metapre.txt > $final/T2D.meta.txt

#zcat $T1D/GCST90014023_buildGRCh38.tsv.gz|awk 'NR==1{print "SNP CHR BPOS A1 A2 freq beta se pval n"};NR>1{print $1,$3,$4,$5,$6,$7,$8,$9,$2,$10}' > $final/T1D.metapre.txt
#awk 'BEGIN {OFS="\t"} {if (NR == 1) print $0, "z"; else if (NR > 1) {if ($8 != 0) print $0, $7/$8; else print $0, "NA"}}' $final/T1D.metapre.txt > $final/T1D.metapre.txt
#awk 'NR==1{print "SNP CHR BPOS A1 A2 freq z pval n"};NR>1{print $1,$2,$3,$4,$5,$6,$11,$9,$10}' $final/T1D.metapre.txt > $final/T1D.metapre.txt
#awk 'NR==1 || $8 < 1e300' $final/T1D.metapre.txt > $final/T1D.metapre.txt
#awk 'NR==1{print "SNP CHR BPOS A1 A2 freq z pval n"};$8 == ($8+0)' $final/T1D.metapre.txt > $final/T1D.metapre.txt
#awk 'NR==1 || $1!="NA"' $final/T1D.metapre.txt > $final/T1D.metapre.txt
#wc $final/T1D.metapre.txt

#zcat $T1D/GCST90014023_buildGRCh38.tsv.gz|awk 'NR==1{print "SNP CHR BPOS A1 A2 freq beta se p N"};NR>1{print $1,$5,$6,$7,$8,$9,$2,$10}' > T1D_metapre.txt
#awk 'NR==1{print "SNP A1 A2 freq b se p N"};$7 == ($7+0)' T1D_metapre.txt > T1D_metapre.txt
#awk 'NR==1{print "SNP A1 A2 beta se p"};NR>1{print $1,$2,$3,$5,$6,$7}' T1D_metapre.txt > T1D_metapre.txt
#awk 'NR==1 || $6 < 1e300' T1D_metapre.txt > T1D_metapre.txt
#awk 'NR==1{print "SNP A1 A2 beta se p"};$6 == ($6+0)' T1D_metapre.txt > T1D_metapre.txt
#head $final/T1D_metapre.txt

#cd mtag
#python mtag.py \
#    --sumstats $final/RA.meta.txt,$final/T2D.meta.txt \
#    --snp_name SNP \
#    --z_name z \
#    --n_name N \
#    --a1_name A1 \
#    --a2_name A2 \
#    --p_name pval \
#--out $final/RA_T2D_meta 

#mtag_munge.py \
#    --sumstats $final/RA.meta.txt \
#    --snp SNP \
#    --a1 A1 \
#    --a2 A2 \
#    --p pval \
#    --N-cas 1605 \
#    --N-con 359589 \
#    --maf-min 0.01 \
#--out RA

#mtag_munge.py \
#    --sumstats $final/T1D_final_4.txt \
#    --snp SNP \
#    --a1 A1 \
#    --a2 A2 \
#    --p pval \
#    --N-cas 1445 \
#    --N-con 362050 \
#    --maf-min 0.01 \
#--out T1D

#mtag.py \
#    --sumstats RA.sumstats.gz,T1D.sumstats.gz \
#	--n_min 0.0 \
#--out ./RA_T2D

#mtag.py  \
#    --sumstats $final/RA.meta.txt,$final/T2D.meta.txt \
#    --beta_name beta\
#    --se_name se\
#	--out $final/RA_T2D
#
#mtag.py  \
#	--sumstats 1_OA2016_hm3samp_NEUR.txt,1_OA2016_hm3samp_SWB.txt \
#	--out ./tutorial_results_1.1NS \
#	--n_min 0.0 \
#      --stream_stdout &  
