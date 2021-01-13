args = commandArgs(trailingOnly=TRUE)
pval <- args[1]
prs_filename <- args[2]
print_header <- args[3]

library(survival)
library(pROC)
library(survminer)
library(ggplot2)
library(foreach)

phenos <- read.csv('~/Dropbox (Encompass Bioscience)/cv_sandbox/data/ukb11102_stroke_mi.csv', sep=',', stringsAsFactors=FALSE)
colnames(phenos) = c("eid", "sex", "birth_year","birth_month", "date1", "date2", "date3", "sys_bp10", "sys_bp11", "sys_bp20", 
                     "sys_bp21", "sys_bp30", "sys_bp31","dia_bp10", "dia_bp11", "dia_bp20", "dia_bp21", "dia_bp30", "dia_bp31", 
                     "reason_nofollow", "date_nofollow", "age_bp_dx1", "age_bp_dx2", "age_bp_dx3", "age_angina_dx1", 
                     "age_angina_dx2", "age_angina_dx3", "year_uk1", "year_uk2", "year_uk3", "mi_dx1", "mi_dx2", "mi_dx3", 
                     "stroke_dx1", "stroke_dx2", "stroke_dx3", "heart_dx11", "heart_dx12", "heart_dx13", "heart_dx14", 
                     "heart_dx21", "heart_dx22", "heart_dx23", "heart_dx24", "heart_dx31", "heart_dx32", "heart_dx33", 
                     "heart_dx34", "height", "smoke1", "smoke2", "smoke3", "e_smoke1", "e_smoke2", "e_smoke3", "race1", 
                     "race2", "race3", "bmi1", "bmi2", "bmi3", "weight1", "weight2", "weight3", "age1", "age2", "age3", 
                     "death1", "death2", "death3", "age_death1", "age_death2", "age_death3" )

# ignore for now
#phenos$dx_bl <- factor(phenos$dx_bl)
#phenos$dx_bl <- factor(phenos$dx_bl, levels(phenos$dx_bl)[c(2,3,4,1)], ordered=TRUE)

# Select a few covariates from subset with one MI to see if they are significant
phenos_clean <- data.frame(phenos$eid, phenos$age1, phenos$sex, phenos$race1, phenos$mi_dx1, phenos$mi_dx2, phenos$bmi1, phenos$age_death1, phenos$sys_bp10, phenos$smoke1, phenos$weight1, phenos$age_angina_dx1, phenos$stroke_dx1)
colnames(phenos_clean) = c("iid", "age1", "sex", "race1", "mi_dx1", "mi_dx2", "bmi1", "age_death1", "sys_bp10", 
                           "smoke1", "weight1", "age_angina_dx1", "stroke_dx1")

#prs <- read.csv(prs_filename, sep='\t')
# path_output = "~/Dropbox (Encompass Bioscience)/ukbb_analysis/ukb_cardio_cv/output"
# prs <- read.table(file.path(path_output, "prs_5e-3.txt"), sep = '\t', header = TRUE)
prs <- read.table('~/Dropbox (Encompass Bioscience)/ukbb_analysis/ukb_cardio_cv/output/prs_5e-3.txt', sep='\t', header=TRUE)
colnames(prs) <- c('iid', 'prs')

pcs <- read.table('~/Dropbox (Encompass Bioscience)/ukbb_analysis/ukb_cardio_cv/output/uk_oneMI_50.eigenvec', sep='\t', header=FALSE)
colnames(pcs) <- c('ix', 'iid', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')

# Too many people included, need to narrow down the list
mi_dx_1 <- subset(phenos_clean, phenos_clean$mi_dx1 != "NA" & phenos_clean$mi_dx1 != "-1" & phenos_clean$mi_dx1 != "-3")
compare <- subset(phenos_clean, phenos_clean$mi_dx1 == phenos_clean$mi_dx2)
no_mi <- subset(phenos_clean, is.na(phenos_clean$mi_dx1))
no_mi_sample <- no_mi[sample(nrow(no_mi), "20000"), ]

# merge new phenos_sample
phenos_samples <- rbind(mi_dx_1, compare, no_mi_sample)

t <- merge(phenos_samples, prs, on='iid')
z <- merge(t, pcs, on='iid')

z$age1 <- as.numeric(z$age1)
z$sex <- as.numeric(z$sex)
z$race1 <- as.numeric(z$race1)
z$mi_dx1 <- as.numeric(z$mi_dx1)
z$mi_dx2 <- as.numeric(z$mi_dx2)
z$bmi1 <- as.numeric(z$bmi1)
z$age_death1 <- as.numeric(z$age_death1)
z$sys_bp10 <- as.numeric(z$sys_bp10)
z$weight1 <- as.numeric(z$weight1)
z$age_angina_dx1 <- as.numeric(z$age_angina_dx1)
z$stroke_dx1 <- as.numeric(z$stroke_dx1)
z$smoke1 <- as.numeric(z$smoke1)

# create survival model based on event age
z$event_age = pmin(z$age1, z$age_angina_dx1, z$mi_dx1, z$stroke_dx1, na.rm=T)
z$event_age[z$event_age < 0] <- NA
hist(z$event_age)
z$had_event = !(is.na(z$age_angina_dx1) & is.na(z$mi_dx1) & is.na(z$stroke_dx1))
z$survival <- with(z, Surv(event_age, had_event))
#mi.coxph <- coxph(Surv( ukdatah_valid$event_age, ukdatah_valid$had_event) ~ prs_3, data=ukdatah_valid)


# time to first cardiac event
z$case_status <- 0
z[ !is.na(z$mi_dx2), ]$case_status <- 1
z[ !is.na(z$mi_dx1), ]$case_status <- 1

# too many individuals so will need to narrow it down to ~50,000 individuals 

# for individual i in N, leave i out of dataset
# with data - i, train model
# predict i, and ask how far away you are from the truth
# repeat for all individuals in N to get cross-validated model
results <- matrix(NA, nrow=dim(z)[1], ncol=10)
for (i in 1:dim(z)[1]) {
  cat(i, '\n')
	# split data up for leave one out cross validation
	training <- z[-c(i), ]
	testing <- z[i, ]
	testing[is.na(testing)] <- 0

	all.fit <- glm(case_status ~ sex + bmi1 + age1 + smoke1 + weight1 + prs + PC1 + PC2 + PC3, data=training, family=binomial(link="logit"))
	all.risk <- as.numeric(predict(all.fit,  se.fit=TRUE, type="response", newdata=testing)$fit)

	enc_prs.fit <- glm(case_status ~ prs + PC1 + PC2 + PC3, data=training, family=binomial(link="logit"))
  enc_prs.risk <- as.numeric(predict(enc_prs.fit,  se.fit=TRUE, type="response", newdata=testing)$fit)

  prs.fit <- glm(case_status ~ prs, data=training, family=binomial(link="logit"))
  prs.risk <- as.numeric(predict(prs.fit,  se.fit=TRUE, type="response", newdata=testing)$fit)

  result.row <- c(testing$iid, all.risk, enc_prs.risk, prs.risk, testing$event_age, testing$case_status, testing$prs, testing$age1, testing$mi_dx1, testing$mi_dx2)
	results[i, ] <- result.row
}

# save all the cross validation results
# results <- write.table(results, file = file.path(path_output, "LOOCV_results.csv"), sep = ',')
results <- data.frame(results, stringsAsFactors=FALSE)
colnames(results) <- c('iid', 'all_risk', 'enc_prs_risk', 'prs_risk', 'event_age', 'case_status', 'prs', 'age1', 'mi_dx1', 'mi_dx2')

results$all_risk <- as.numeric(as.numeric(results$all_risk))
results$prs_risk <- as.numeric(as.numeric(results$prs_risk))

#dementia.only <- results[ results$status == 1, ]
#ad.only <- results[ results$case_status == 1, ]

# create column with 0 - no heart attack, 1 - 1 heart attack and 2 - second heart attack
results$mi[results$mi_dx1==0] <- 0
results$mi[results$mi_dx1!=0] <- 1
results$mi[results$mi_dx2!=0] <- 2

one_mi <- results[ results$mi >= 1, ]

# analyze the results
# AD case control, diagnosis baseline versus risk
mi_corr <- c(cor(as.numeric(results$mi), results$all_risk), cor(as.numeric(results$mi), results$prs_risk))

# looking at least  one Mi
one_mi_corr <- c(cor(as.numeric(one_mi$mi), one_mi$all_risk), cor(as.numeric(one_mi$mi), one_mi$prs_risk))

# ROC 
all_roc <- roc(results$case_status, results$all_risk)
all_auc <- as.numeric(auc(all_roc))

prs_roc <- roc(results$case_status, results$prs_risk)
prs_auc <- as.numeric(auc(prs_roc))

# # graphs for time to first event
# # select top 50 percent of using all risk score
# risk_quantiles <- quantile(results$all_risk)
# risk_cutoff <- risk_quantiles[2]
# z$group <- NA
# z[ z$iid %in% results[ results$all_risk >= risk_cutoff, ]$iid, ]$group <- "Top 50% Risk"
# z[ z$iid %in% results[ results$all_risk < risk_cutoff, ]$iid, ]$group <- "Bottom 50% Risk"

# surv.cox <- survfit(Surv(event_age, case_status) ~ group, data=z)
# ggsurvplot(surv.cox, xlim=c(40, 70), conf.int.style='ribbon', conf.int=TRUE, pval=TRUE,  
#            ggtheme = theme_bw(), pval.coord=c(55, 0.25), ylab="% with First Heart Attack", xlab="Age")

# # graphs for time to second event
# # select top 50 percent of using all risk score
# results2 <- results[ results$mi_dx2 > 0, ]
# risk_quantiles <- quantile(results2$all_risk)
# risk_cutoff <- risk_quantiles[2]
# zz <- z[ !is.na(z$mi_dx2), ]
# zz$group <- NA
# zz[ zz$iid %in% results2[ results2$all_risk >= risk_cutoff, ]$iid, ]$group <- "Top 50% Risk"
# zz[ zz$iid %in% results2[ results2$all_risk < risk_cutoff, ]$iid, ]$group <- "Bottom 50% Risk"

# surv.cox2 <- survfit(Surv(event_age, case_status) ~ group, data=zz)
# ggsurvplot(surv.cox2, xlim=c(40, 70), conf.int.style='ribbon', conf.int=TRUE, pval=TRUE,  ggtheme = theme_bw(), pval.coord=c(55, 0.25), ylab="% with Second Heart Attack")


# results$mi <- as.factor((results$mi))
# ggplot(data.frame(results$mi, results$all_risk, results$prs_risk), aes(x=results.mi, y=results.prs_risk, fill=results.mi)) + geom_boxplot() + geom_smooth(method = "lm", se=FALSE, color="blue", aes(group=1)) + labs(x="# of Heart Attacks", y="Polygenic Risk Score") + annotate("text", x=2.5, y=0.2, label="Correlation: 0.112, p=2.2e-16", hjust = 1, vjust = -18) 
# results$all_risk <-as.factor((results$all_risk))
# cor.test(as.numeric(results$mi), as.numeric(results$prs_risk))
