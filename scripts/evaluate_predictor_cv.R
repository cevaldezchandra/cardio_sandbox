args = commandArgs(trailingOnly=TRUE)
pval <- args[1]
prs_filename <- args[2]
print_header <- args[3]

library(survival)
library(pROC)
library(survminer)
library(ggplot2)

phenos <- read.csv('cv_final.phenos', sep=',', stringsAsFactors=FALSE)
# ignore for now
#phenos$dx_bl <- factor(phenos$dx_bl)
#phenos$dx_bl <- factor(phenos$dx_bl, levels(phenos$dx_bl)[c(2,3,4,1)], ordered=TRUE)

prs <- read.csv(prs_filename, sep='\t')
colnames(prs) <- c('iid', 'prs')

pcs <- read.csv('uk_cv_filtered.eigenvec', sep=' ', header=F)
colnames(pcs) <- c('ix', 'iid', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16', 'PC17', 'PC18', 'PC19', 'PC20')

t <- merge(phenos, prs, on='iid')
z <- merge(t, pcs, on='iid')

z$age <- as.numeric(z$age)
z$sex <- as.numeric(z$sex)
z$race <- as.numeric(z$race)
z$mi_dx1 <- as.numeric(z$mi_dx1)
z$mi_dx2 <- as.numeric(z$mi_dx2)
z$bmi1 <- as.numeric(z$bmi1)
z$age_death1 <- as.numeric(z$age_death1)
z$sys_bp <- as.numeric(z$sys_bp)
z$smoke_status <- as.numeric(z$smoke_status)
z$survival <- with(z, Surv(mi_dx1, status==1))

# time to first cardiac event
z$case_status <- z$mi_dx1

# looking at second heart attack
z$case_status <- 1
z[is.na(z$mi_dx2),]$case_status <- 0

z$second_event_age = z$mi_dx2
z$second_event_age[z$second_event_age < 0] <- NA
#hist(z$second_event_age)
z$had_second_event = !(is.na(z$mi_dx2)

mi.coxph <- coxph(Surv( ukdatah_valid$event_age, ukdatah_valid$had_event) ~ miprs0001, data=ukdatah_valid)


# for individual i in N, leave i out of dataset
# with data - i, train model
# predict i, and ask how far away you are from the truth
# repeat for all individuals in N to get cross-validated model
results <- matrix(NA, nrow=dim(z)[1], ncol=8)
for (i in 1:dim(z)[1]) {
	# split data up for leave one out cross validation
	training <- z[-c(i), ]
	testing <- z[i, ]
	testing[is.na(testing)] <- 0

	#all.fit <- glm(case_status ~ sex + bmi1 + age + smoke_status + weight1 + prs + PC1 + PC2 + PC3, data=training, family=binomial(link="logit"))
	all.fit <- glm(had_second_event ~ sex + bmi1 + age + smoke_status + weight1 + prs + PC1 + PC2 + PC3, data=training)
	all.risk <- as.numeric(predict(all.fit,  se.fit=TRUE, type="response", newdata=testing)$fit)

	#apoe4.fit <- glm(case_status ~ apoe4, data=training, family=binomial(link="logit"))
	#apoe4.risk <- as.numeric(predict(apoe4.fit,  se.fit=TRUE, type="response", newdata=testing)$fit)

	#prs.fit <- glm(case_status ~ prs + PC1 + PC2 + PC3, data=training, family=binomial(link="logit"))
	prs.fit <- glm(had_second_event ~ prs + PC1 + PC2 + PC3, data=training)

	prs.risk <- as.numeric(predict(prs.fit,  se.fit=TRUE, type="response", newdata=testing)$fit)

    result.row <- c(testing$iid, all.risk, prs.risk, testing$mi_dx2, testing$had_second_event, testing$status, testing$prs, testing$age)
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