## Uk Biobank PRS script
library(ggplot2)
library(dplyr)
library(rms)
library(pROC)
library(reshape2)

setwd("/Users/crystalvaldez/Dropbox (Encompass Bioscience)/cv_sandbox/data")
# Read in all participant data
ukdata <- read.csv("ukb11102_stroke_mi.csv")
colnames(ukdata) = c("eid", "sex", "birth_year","birth_month", "date1", "date2", "date3", "sys_bp10", "sys_bp11", "sys_bp20", "sys_bp21", "sys_bp30", "sys_bp31","dia_bp10", "dia_bp11", "dia_bp20", "dia_bp21", "dia_bp30", "dia_bp31", "reason_nofollow", "date_nofollow", "age_bp_dx1", "age_bp_dx2", "age_bp_dx3", "age_angina_dx1", "age_angina_dx2", "age_angina_dx3", "year_uk1", "year_uk2", "year_uk3", "mi_dx1", "mi_dx2", "mi_dx3", "stroke_dx1", "stroke_dx2", "stroke_dx3", "heart_dx11", "heart_dx12", "heart_dx13", "heart_dx14", "heart_dx21", "heart_dx22", "heart_dx23", "heart_dx24", "heart_dx31", "heart_dx32", "heart_dx33", "heart_dx34", "height", "smoke1", "smoke2", "smoke3", "e_smoke1", "e_smoke2", "e_smoke3", "race1", "race2", "race3", "bmi1", "bmi2", "bmi3", "weight1", "weight2", "weight3", "age1", "age2", "age3", "death1", "death2", "death3", "age_death1", "age_death2", "age_death3" )

heart_dx_map = c(-7, -3, 1, 2, 3, 4)
names = c("None", "NoAnswer", "MI", "Angina", "Stroke")

uk_twoMI <- subset(ukdata, !is.na(ukdata$mi_dx1) & !is.na(ukdata$mi_dx2))

plot(uk_twoMI$mi_dx1, uk_twoMI$mi_dx2)
sum(uk_twoMI$mi_dx1 == uk_twoMI$mi_dx2)
sum(uk_twoMI$mi_dx1 == uk_twoMI$mi_dx3, na.rm=TRUE)
sum(is.na(uk_twoMI$mi_dx3))

### read in cardio risk scores
path1 = "/Users/crystalvaldez/Dropbox (Encompass Bioscience)/ukbb_analysis/ukb_cardio_cv/output"
#prs1 <- read.table("all_prs1.txt")
#names(prs1) <- c("eid", "prs1")
#prs2 <- read.table("all_prs2.txt")
#names(prs2) <- c("eid", "prs2")

## Read in MI risk scores
prs_3 <- read.table(file.path(path1, "prs_5e-3.txt"), sep = '\t', header = TRUE)
prs_4 <- read.table(file.path(path1, "prs_5e-4.txt"), sep = '\t', header = TRUE)
prs_5 <- read.table(file.path(path1, "prs_5e-5.txt"), sep = '\t', header = TRUE)

names(prs_3) <- c("eid", "prs_5e-3")
names(prs_4) <- c("eid", "prs_5e-4")
names(prs_5) <- c("eid", "prs_5e-5")
mi_prs <- data.frame(eid = prs_3$eid, prs_3=prs_3$`prs_5e-3`, prs_4=prs_4$`prs_5e-4`, prs_5=prs_5$`prs_5e-5`)
#mi_prs <- data.frame(eid = mi_prs01$eid, miprs01=mi_prs01$miprs01, miprs005=mi_prs005$miprs005, miprs001 = mi_prs001$miprs001, miprs0005 = mi_prs0005$miprs0005, miprs0001=mi_prs0001$miprs0001, miprs14=mi_prs14$V2, miprs15=mi_prs15$V2, miprs16=mi_prs16$V2)
# Has CVD?


# Validation split

#ukdatah = subset(ukdata, ukdata$heart_dx11 != "-3" & race1%in%c(1001,1002,1003))
ukdatah = subset(ukdata, ukdata$heart_dx11 != "-3")

ukdatah$heart_problem = ukdatah$heart_dx11 %in% c(1,2,3)
ukdatah$heart_problem[ukdatah$heart_dx11=="-3"] = NA
#ukdatah = merge(ukdatah, prs1, by="eid")
ukdatah = merge(ukdatah, mi_prs, by="eid")

uk_euro = ukdatah$eid
# Sample a subset for clumping PRS
#clump_samp = sample(uk_euro, 5000, replace=F)
# Write out a fam file
#clump_frame = data.frame(clump_samp, clump_samp)
#write.table(clump_frame, "clump_ids.txt", row.names=F, col.names=F, sep="\t", quote=F)

ukdatah$heart_problem[ukdatah$heart_dx11==4] <- FALSE

# Let's get a clean dataset for heart analysis. We will remove age less than 50, BMI greater than 40 (very high risk)
ukdatah_clean <- subset(ukdatah, age1 > 50 & bmi1 < 40 & heart_dx11!=4 )

# pull out 30% for validation cohort
v_prop = 0.3
ukdatah_clean$heart_vsplit = NA
ukdatah_clean$heart_vsplit[ukdatah_clean$heart_problem == T] = sample(c(0,1), sum(ukdatah_clean$heart_problem==T, na.rm=T), prob=c(0.7,0.3), replace=T)
ukdatah_clean$heart_vsplit[ukdatah_clean$heart_problem == F] = sample(c(0,1), sum(ukdatah_clean$heart_problem==F, na.rm=T), prob=c(0.7,0.3), replace=T)

ukdatah_valid = subset(ukdatah_clean, heart_vsplit==1)

table(ukdatah_valid$heart_problem)

# AUCs not so good here
auc(roc(ukdatah_valid$heart_problem, ukdatah_valid$prs_3, algorithm=2))
auc(roc(ukdatah_valid$heart_problem, ukdatah_valid$prs_4, algorithm=2))
auc(roc(ukdatah_valid$heart_problem, ukdatah_valid$prs_5, algorithm=2))

# Without high BP samples:
ukdatah_clean_nobp <- subset(ukdatah, age1 > 50 & bmi1 < 40 & heart_dx11!=4 )

# pull out 30% for validation cohort
v_prop = 0.3
ukdatah_clean_nobp$heart_vsplit = NA
ukdatah_clean_nobp$heart_vsplit[ukdatah_clean_nobp$heart_problem == T] = sample(c(0,1), sum(ukdatah_clean_nobp$heart_problem==T, na.rm=T), prob=c(0.7,0.3), replace=T)
ukdatah_clean_nobp$heart_vsplit[ukdatah_clean_nobp$heart_problem == F] = sample(c(0,1), sum(ukdatah_clean_nobp$heart_problem==F, na.rm=T), prob=c(0.7,0.3), replace=T)

ukdatah_valid_nobp = subset(ukdatah_clean_nobp, heart_vsplit==1)

table(ukdatah_valid_nobp$heart_problem)

# AUCs still not very good here
auc(roc(ukdatah_valid_nobp$heart_problem, ukdatah_valid_nobp$prs_3, algorithm=2))
auc(roc(ukdatah_valid_nobp$heart_problem, ukdatah_valid_nobp$prs_4, algorithm=2))
auc(roc(ukdatah_valid_nobp$heart_problem, ukdatah_valid_nobp$prs_5, algorithm=2))

# Let's try a survival model
ukdatah_clean$event_age = pmin(ukdatah_clean$age1, ukdatah_clean$age_angina_dx1, ukdatah_clean$mi_dx1, ukdatah_clean$stroke_dx1, na.rm=T)
ukdatah_clean$event_age[ukdatah_clean$event_age < 0] <- NA
hist(ukdatah_clean$event_age)
ukdatah_clean$had_event = !(is.na(ukdatah_clean$age_angina_dx1) & is.na(ukdatah_clean$mi_dx1) & is.na(ukdatah_clean$stroke_dx1))

ukdatah_valid = subset(ukdatah_clean, heart_vsplit==1)

# Let's do a Cox PH model
table(ukdatah_valid$had_event)

mi.coxph <- coxph(Surv( ukdatah_valid$event_age, ukdatah_valid$had_event) ~ prs_3, data=ukdatah_valid)
mi.coxph <- coxph(Surv( ukdatah_valid$event_age, ukdatah_valid$had_event) ~ prs_3, data=ukdatah_valid)
mi.coxph <- coxph(Surv( ukdatah_valid$event_age, ukdatah_valid$had_event) ~ prs_3, data=ukdatah_valid)

# It's super effective!
# Build KM curve with top quantile
ukdatah_valid$topq_miprs <- ukdatah_valid$prs_3 > quantile(ukdatah_valid$prs_3, 0.75)

# Add in covariates from phenos
coxph(Surv( ukdatah_valid$event_age, ukdatah_valid$had_event) ~ topq_miprs + bmi1 + sex + smoke1 + age1 + weight1, data=ukdatah_valid)

# Looks pretty decent for cardiac events
miq1.npsurv <- npsurv(Surv( ukdatah_valid$event_age, ukdatah_valid$had_event) ~ ukdatah_valid$topq_miprs)

pdf("/Users/crystalvaldez/Dropbox (Encompass Bioscience)/ukbb_analysis/ukb_cardio_cv/images/prs_eventAge.pdf", height=5, width=7)
#### Plot in current deck ######
survplot(miq1.npsurv, xlab="Age", xlim=c(30,72), col=c("firebrick", "steelblue"), label.curves=F, ylim=c(0.7,1), lty=1, n.risk=F, cex.text=3, lwd=2, col.fill=alpha(col=c("firebrick", "steelblue"), 0.3))
dev.off()

pdf("/Users/crystalvaldez/Dropbox (Encompass Bioscience)/ukbb_analysis/ukb_cardio_cv/images/survival_diffprs.pdf", height=6, width=9)
survdiffplot(miq1.npsurv, xlab="Age", xlim=c(40,72), col=c("firebrick", "steelblue"), lty=1, n.risk=F)
dev.off()



#### Stopped here ######
# Let's simulate a study. 
ukdatah_sim <- subset(ukdatah_clean,  heart_vsplit==0)

study_size = 2000 # patients
start_age = 55

ukdatah_sim$study_event = ukdatah_sim$event_age - start_age
ukdatah_study_pop = subset(ukdatah_sim, study_event>0) # event after age of 45
ukdatah_study_pop$topq_miprs <- ukdatah_study_pop$prs_3 > quantile(ukdatah_valid$prs_3, 0.75) # top quartile from the validation set

# Simulate two studies and time to event
# randomly sample rows, get time to event
ukdatah_study_topq = subset(ukdatah_study_pop, topq_miprs==TRUE)
ukdatah_study_highbp = subset(ukdatah_study_pop, heart_dx11%in%c(1,2,3,4))
ukdatah_study_highbp_topq = subset(ukdatah_study_pop, heart_dx11%in%c(1,2,3,4) & topq_miprs==TRUE)

event4_year_enrich = rep(NA, 1000)
event4_year_normal = rep(NA, 1000)
for (i in 1:1000) {
  enrich_samp = ukdatah_study_topq[sample(nrow(ukdatah_study_topq), study_size),]
  event4_year_enrich[i] = sum(enrich_samp$had_event[enrich_samp$study_event < 10])
  
  normal_samp = ukdatah_study_pop[sample(nrow(ukdatah_study_pop), study_size),]
  event4_year_normal[i] = sum(normal_samp$had_event[normal_samp$study_event < 10])
  
}

enrich_frame <- data.frame(enrich = event4_year_enrich, normal= event4_year_normal)
enrich_melt <- melt(enrich_frame)

ggplot(data=enrich_melt) + geom_boxplot(aes(x=variable, y=value)) + labs(x="Enrollment Style", y="Number of Events after 10 years") + theme_bw()
ggsave("enrich_all_pop.pdf")





























