args = commandArgs(trailingOnly=TRUE)
pval <- args[1]
prs_filename <- args[2]
print_header <- args[3]

library(survival)
library(pROC)
library(survminer)
library(ggplot2)

mace_phenos <- read.csv('/Users/crystalvaldez/Dropbox (Encompass Bioscience)/ukbb_analysis/ukb_mace_he/output/ukb_mace_hes_phenos.csv', sep=',', stringsAsFactors=FALSE)

# Select a few covariates from subset with one MI to see if they are significant
#phenos_clean <- data.frame(phenos$eid, phenos$age1, phenos$sex, phenos$race1, phenos$mi_dx1, phenos$mi_dx2, phenos$bmi1,
#                           phenos$age_death1, phenos$sys_bp10, phenos$smoke1, phenos$weight1, phenos$age_angina_dx1, phenos$stroke_dx1)
#colnames(phenos_clean) = c("iid", "age1", "sex", "race1", "mi_dx1", "mi_dx2", "bmi1", "age_death1", "sys_bp10", 
#                           "smoke1", "weight1", "age_angina_dx1", "stroke_dx1")

#prs <- read.csv(prs_filename, sep='\t')
path_output = "/Users/crystalvaldez/Dropbox (Encompass Bioscience)/ukbb_analysis/ukb_mace_he/output"
prs <- read.table(file.path(path_output, "prs_5e-3.txt"), sep = '\t', header = TRUE)
colnames(prs) <- c('iid', 'prs')

pcs <- read.table(file.path(path_output, "uk_hes_cv_pcs.eigenvec"), sep=' ', header=F)
colnames(pcs) <- c('ix', 'iid', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10',
                   'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16', 'PC17', 'PC18', 'PC19', 'PC20')



t <- merge(mace_phenos, prs, on='iid')
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
#z$case_status <- z$mi_dx1

# looking at second heart attack
#z$case_status <- 1
#z[is.na(z$mi_dx2),]$case_status <- 0
#z$second_event_age = z$mi_dx2
#z$second_event_age[z$second_event_age < 0] <- NA
#hist(z$second_event_age)
#z$had_second_event = !(is.na(z$mi_dx2)

# too many individuals so will need to narrow it down to ~50,000 individuals 

# for individual i in N, leave i out of dataset
# with data - i, train model
# predict i, and ask how far away you are from the truth
# repeat for all individuals in N to get cross-validated model
results <- matrix(NA, nrow=dim(z)[1], ncol=9)
for (i in 1:dim(z)[1]) {
  # split data up for leave one out cross validation
  training <- z[-c(i), ]
  testing <- z[i, ]
  testing[is.na(testing)] <- 0
  
  #all.fit <- glm(case_status ~ sex + bmi1 + age + smoke_status + weight1 + prs + PC1 + PC2 + PC3, data=training, family=binomial(link="logit"))
  all.fit <- glm(event_age ~ sex + bmi1 + age1 + smoke1 + weight1 + prs + PC1 + PC2 + PC3, data=training)
  all.risk <- as.numeric(predict(all.fit,  se.fit=TRUE, type="response", newdata=testing)$fit)
  
  #apoe4.fit <- glm(case_status ~ apoe4, data=training, family=binomial(link="logit"))
  #apoe4.risk <- as.numeric(predict(apoe4.fit,  se.fit=TRUE, type="response", newdata=testing)$fit)
  
  #prs.fit <- glm(case_status ~ prs + PC1 + PC2 + PC3, data=training, family=binomial(link="logit"))
  prs.fit <- glm(event_age ~ prs + PC1 + PC2 + PC3, data=training)
  prs.risk <- as.numeric(predict(prs.fit,  se.fit=TRUE, type="response", newdata=testing)$fit)
  
  result.row <- c(testing$iid, all.risk, prs.risk, testing$event_age, testing$had_event, testing$prs, testing$age1, testing$mi_dx1, testing$mi_dx2)
  #result.row <- c(testing$iid, all.risk, prs.risk, testing$mi_dx2, testing$case_status, testing$status, testing$prs, testing$age)
  #result.row <- c(testing$iid, all.risk, apoe4.risk, prs.risk, testing$dx_bl, testing$case_status, testing$status, testing$prs, testing$apoe4, testing$time)
  # cat(result.row, '\n')
  results[i, ] <- result.row
}

results <- data.frame(results, stringsAsFactors=FALSE)
colnames(results) <- c('iid', 'all_risk', 'prs_risk', 'mi_dx2', 'case_status', 'status', 'prs', 'age')

results$all_risk <- as.numeric(as.numeric(results$all_risk))
results$prs_risk <- as.numeric(as.numeric(results$prs_risk))

#dementia.only <- results[ results$status == 1, ]
#ad.only <- results[ results$case_status == 1, ]


# analyze the results

# ROC 
all_roc <- roc(results$case_status, results$all_risk)
all_auc <- as.numeric(auc(all_roc))

prs_roc <- roc(results$case_status, results$prs_risk)
prs_auc <- as.numeric(auc(prs_roc))


# select top 50 percent of using all risk score
risk_quantiles <- quantile(results$all_risk)
risk_cutoff <- risk_quantiles[2]
z$group <- NA
z[ z$iid %in% results[ results$all_risk >= risk_cutoff, ]$iid, ]$group <- "Top 50% Risk"
z[ z$iid %in% results[ results$all_risk < risk_cutoff, ]$iid, ]$group <- "Bottom 50% Risk"

surv.cox <- survfit(Surv(status) ~ group, data=z)
ggsurvplot(surv.cox, xlim=c(60, 80), conf.int.style='ribbon', conf.int=TRUE, pval=TRUE,  ggtheme = theme_bw(), pval.coord=c(65, 0.25), ylab="% with second heart attack")
survplot(miq1.npsurv, xlab="Age", xlim=c(30,72), col=c("firebrick", "steelblue"), label.curves=F, ylim=c(0.7,1), lty=1, n.risk=F, cex.text=3, lwd=2, col.fill=alpha(col=c("firebrick", "steelblue"), 0.3))


ggplot(data.frame(results$dx_bl, results$all_risk, results$prs_risk), aes(x=results.dx_bl, y=results.prs_risk, fill=results.dx_bl)) + geom_boxplot() + geom_smooth(method = "lm", se=FALSE, color="blue", aes(group=1)) + labs(x="Disease Stage", y="Polygenic Risk Score") + annotate("text", x=2.5, y=0.2, label="Correlation: 0.207, p=3e-9") + scale_x_discrete(labels=c("1" = "Cognitively Normal", "2"="Early Cognitive Impairment", "3"="Late Cognitive Impairment", "4"="Alzheimer's Disease"))