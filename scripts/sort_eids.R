setwd("/Users/crystalvaldez/Dropbox (Encompass Bioscience)/cv_sandbox/data")

# Extract eids from phenotypes
cv_phenos <- read.csv("cv_final.phenos", header = TRUE, sep = ",")
eids <- data.frame(cv_phenos[,1])
write.table(eids, "total_eids.txt", row.names = FALSE, col.names = FALSE)

# ethnicity specific eids
eur_eids <- read.table(file='~/Dropbox (Encompass Bioscience)/cv_sandbox/data/eur_eids.txt')
aa_eids <- read.table(file='~/Dropbox (Encompass Bioscience)/cv_sandbox/data/aa_eids.txt', sep=',', quote=F, col.names=F, row.names=F)
ind_eids <- read.table(file='~/Dropbox (Encompass Bioscience)/cv_sandbox/data/ind_eids.txt', sep=',', quote=F, col.names=F, row.names=F)

# create assoc file from cardio gwas summary statistics (http://www.cardiogramplusc4d.org/data-downloads/)
setwd("/Users/crystalvaldez/Dropbox (Encompass Bioscience)/cv_sandbox/data/cardiogram_gwas_results")
cardio_gwas <- read.table("CARDIoGRAM_GWAS_RESULTS.txt", header = TRUE, sep = '\t')
assoc <- data.frame(cardio_gwas$SNP, cardio_gwas$pvalue)
colnames(assoc) <- c("SNP", "pvalue")
write.table(assoc, "cv_stage1.assoc", quote = F, row.names = FALSE, sep = '\t')


# take log odds from logistic regression and get the probabilties
# beta file format: Chromosome	Position	RSID	Effect_allele	Non_Effect_allele	Beta	SE	Pvalue
beta_file <- data.frame(cardio_gwas$chr_pos_.b36., cardio_gwas$SNP, cardio_gwas$reference_allele, cardio_gwas$other_allele,
                        cardio_gwas$log_odds, cardio_gwas$log_odds_se, cardio_gwas$pvalue)

chr <- data.frame(cardio_gwas$chr_pos_.b36.)
chr_split <- data.frame(do.call('rbind', strsplit(as.character(chr$cardio_gwas.chr_pos_.b36.),':',fixed = TRUE)))
chr_number <- data.frame(gsub("chr","",chr_split$X1))

beta_new <- data.frame(chr_number$gsub..chr.......chr_split.X1., chr_split$X2, cardio_gwas$SNP, cardio_gwas$reference_allele, cardio_gwas$other_allele,
                       cardio_gwas$log_odds, cardio_gwas$log_odds_se, cardio_gwas$pvalue)
colnames(beta_new) <- c("Chromosome","Position","RSID","Effect_allele","Non_Effect_allele","Beta","SE","Pvalue")
write.table(beta_new, "beta_file.txt", sep = '\t', quote = F, row.names = FALSE, col.names = TRUE)

# pick best pvalue
setwd("/Users/crystalvaldez/Dropbox (Encompass Bioscience)/cv_sandbox/data/output")
prs_3 <- read.table("prs_5e-3.txt", sep = '', header = TRUE)
prs_3 <- data.frame(prs_3, type="5e-3")
prs_4 <- read.table("prs_5e-4.txt", sep = '', header = TRUE)
prs_4 <- data.frame(prs_4, type="5e-4")
prs_5 <- read.table("prs_5e-5.txt", sep = '', header = TRUE)
prs_5 <- data.frame(prs_5, type="5e-5")
prs_6 <- read.table("prs_5e-6.txt", sep = '', header = TRUE)
prs_6 <- data.frame(prs_6, type="5e-6")
prs_7 <- read.table("prs_5e-7.txt", sep = '', header = TRUE)
prs_7 <- data.frame(prs_7, type="5e-7")
prs_8 <- read.table("prs_5e-8.txt", sep = '', header = TRUE)
prs_8 <- data.frame(prs_8, type="5e-8")

pvalue_compare <- data.frame(rbind(prs_3,prs_4,prs_5,prs_6,prs_7,prs_8))
ggplot(pvalue_compare, aes(predictions, fill=type)) + geom_density(alpha = 0.2)
