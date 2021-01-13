
# Pre-processing of files before running prediction - pcas and betas are stored within data folder, output
# is where all outputs from scripts are stored

# Write out minor allele report to plink.frq 
# plink --bfile ukb_cv --freq 

# Filter out all ambiguous snps (A T and C G) from plink.frq, create ambig_snps.rsids
# plink --bfile ukb_cv --chr 1-22 --exclude ambig_snps.rsids --mac 1 --make-bed --out uk_cv_filtered

# Convert bed file to vcf file
# plink --bfile uk_cv_filtered --recode vcf-iid bgz --out uk_cv_filtered_vcf
# (old) plink --bfile ../ukb_cal_chr1_v2 --recode vcf-iid bgz --out chr1_test --keep eids_inds_plink.txt

# Pre-compute PCs from genotypes
# plink --bfile /ukbb_data/ukbb_gene/ukb_cardio_cv/vcf/uk_cv_filtered --pca --out ../pcas/uk_cv_filtered


for i in 5e-3 5e-4 5e-5 5e-6 5e-7 5e-8;
do 
	echo $i

	# clump the SNPs using CEU
	plink --bfile /ukbb_data/ukbb_gene/ukb_cardio_cv/data/uk_cv_filtered --clump /ukbb_data/ukbb_gene/ukb_cardio_cv/data/cv_stage1.assoc --clump-p1 $i --clump-r2 0.2 --clump-kb 1000 --out /ukbb_data/ukbb_gene/ukb_cardio_cv/output/cv_$i
	#plink --bfile /ukbb_data/ukbb_gene/ukb_cardio_cv/vcf/uk_cv_filtered --clump /ukbb_data/ukbb_gene/ukb_cardio_cv/cv_stage1.assoc --clump-p1 5e-3 --clump-r2 0.2 --clump-kb 1000 --out /ukbb_data/ukbb_gene/ukb_cardio_cv/output/

	# get the list of rsids from 
	tail -n +2 /ukbb_data/ukbb_gene/ukb_cardio_cv/output/cv_$i.clumped | awk '{$2=$2};1' | cut -d' ' -f3 > /ukbb_data/ukbb_gene/ukb_cardio_cv/output/clumped_$i.rsids
	#tail -n +2 /extra/adni/output/igap_$i.clumped | awk '{$2=$2};1' | cut -d' ' -f3 > /extra/adni/output/clumped_$i.rsids

	# create vcf with only clumped SNPs
    plink --bfile /ukbb_data/ukbb_gene/ukb_cardio_cv/data/uk_cv_filtered --extract /ukbb_data/ukbb_gene/ukb_cardio_cv/output/clumped_$i.rsids --recode vcf --out /ukbb_data/ukbb_gene/ukb_cardio_cv/output/cv_clumped_$i
	#plink --bfile /extra/adni/data/adni_autosomes --extract /extra/adni/output/clumped_$i.rsids --recode vcf --out /extra/adni/output/adni_clumped_$i

	# generate PRS - extracted betas from CARDIoGRAM GWAS dataset
    python /ukbb_data/ukbb_gene/ukb_cardio_cv/cv_sandbox/scripts/PRS_pipeline.py -v /ukbb_data/ukbb_gene/ukb_cardio_cv/output/cv_clumped_$i.vcf -b /ukbb_data/ukbb_gene/ukb_cardio_cv/data/beta_file.txt -a /ukbb_data/ukbb_gene/ukb_cardio_cv/data/uk_cv_filtered.eigenvec  -o /ukbb_data/ukbb_gene/ukb_cardio_cv/output/prs_$i\.txt
	#python /extra/adni/scripts/PRS_pipeline.py -v  /extra/adni/output/adni_clumped_$i.vcf -b /extra/adni/data/IGAP_stage_1.txt -a /extra/adni/data/adni.eigenvec -o /extra/adni/output/prs_$i\.txt
done

# evaluate the PRS to use
#for i in 5e-3 5e-4 5e-5 5e-6 5e-7 5e-8;
#do 
#	echo $i

	#call evaluate_predictor.R
#	if [[ $i == 5e-3 ]]; then
#		Rscript evaluate_predictor_cv.R $i /ukbb_data/ukbb_gene/ukb_cardio_cv/output/prs_$i\.txt TRUE >> /ukbb_data/ukbb_gene/ukb_cardio_cv/output/cv_predictor_eval.txt
		#Rscript evaluate_predictor.R $i /extra/adni/output/prs_$i\.txt TRUE >> /extra/adni/output/ad_predictor_eval.txt
#	else
#		Rscript evaluate_predictor_cv.R $i /ukbb_data/ukbb_gene/ukb_cardio_cv/output/prs_$i\.txt FALSE >> /ukbb_data/ukbb_gene/ukb_cardio_cv/output/cv_predictor_eval.txt
		#Rscript evaluate_predictor.R $i /extra/adni/output/prs_$i\.txt FALSE >> /extra/adni/output/ad_predictor_eval.txt
#	fi
#done

 
