
library(survival)

# Read in participant data - taken from Julian's initial analysis
ukdata <- read.csv("~/Dropbox (Encompass Bioscience)/Encompass Shared/UKbiobank/KM_curves/ukb11102_stroke_mi.csv", check.names = FALSE)
num_ids <- data.frame(colnames(ukdata))
colnames(ukdata) = c("eid", "sex", "birth_year","birth_month", "date1", "date2", "date3", "sys_bp10", 
                     "sys_bp11", "sys_bp20", "sys_bp21", "sys_bp30", "sys_bp31","dia_bp10", "dia_bp11", 
                     "dia_bp20", "dia_bp21", "dia_bp30", "dia_bp31", "reason_nofollow", "date_nofollow", 
                     "age_bp_dx1", "age_bp_dx2", "age_bp_dx3", "age_angina_dx1", "age_angina_dx2", "age_angina_dx3",
                     "year_uk1", "year_uk2", "year_uk3", "mi_dx1", "mi_dx2", "mi_dx3", "stroke_dx1", "stroke_dx2", 
                     "stroke_dx3", "heart_dx11", "heart_dx12", "heart_dx13", "heart_dx14", "heart_dx21", "heart_dx22", 
                     "heart_dx23", "heart_dx24", "heart_dx31", "heart_dx32", "heart_dx33", "heart_dx34", "height", 
                     "smoke1", "smoke2", "smoke3", "e_smoke1", "e_smoke2", "e_smoke3", "race1", "race2", "race3", 
                     "bmi1", "bmi2", "bmi3", "weight1", "weight2", "weight3", "age1", "age2", "age3", "death1", 
                     "death2", "death3", "age_death1", "age_death2", "age_death3" )
combine <- data.frame(num_ids, colnames(ukdata))
#write.csv(combine, file="~/Dropbox (Encompass Bioscience)/cv_sandbox/data/ukb11102_stroke_mi_numeric_headers.csv")

# Test 1 -> Keep individuals age when diagnosed with a heart attack (myocardial infarction)
#uk_oneMI <- subset(ukdata, !is.na(ukdata$mi_dx1)) # 11,533 individuals

# Parse out individuals that do not know age they were diagnosed with first MI (-1) or did not answer (-3)
#uk_oneMI <- uk_oneMI[!uk_oneMI$mi_dx1 %in% c(-1,-3),] # 11,186 individuals

# Parse out individuals that did not answer their smoking status (-3)
#uk_oneMI <- uk_oneMI[!uk_oneMI$smoke1 %in% c(-3),] # 11,109

# Select a few covariates from subset with one MI to see if they are significant
#oneMI_1 <- data.frame(uk_oneMI$eid, uk_oneMI$age1, uk_oneMI$sex, uk_oneMI$race1, uk_oneMI$mi_dx1, uk_oneMI$mi_dx2, uk_oneMI$bmi1,
#                      uk_oneMI$age_death1, uk_oneMI$sys_bp10, uk_oneMI$smoke1, uk_oneMI$weight1)
#colnames(oneMI_1) <- c("iid", "age", "sex", "race", "mi_dx1", "mi_dx2", "bmi1", "age_death1", "sys_bp", "smoke_status", "weight1")

# Create status column in dataframe for survival object: 0 = alive, 1 = dead 
#oneMI_1$status <- with(oneMI_1, ifelse(oneMI_1$age_death1 == "NA", 0, 1))
#oneMI_1$status[is.na(oneMI_1$status)] <- 0

# Create a survival object
#oneMI_1$survival <- with(oneMI_1, Surv(mi_dx1, status==1))

# Create predictor with sex, bmi, race, sys_bp, smoke_status
#res.cox <- coxph(survival ~ sex + bmi1 + age + smoke_status + weight1, data=oneMI_1)
#(res.zph <- cox.zph(res.cox))
#summary(res.cox)

#setwd("/Users/crystalvaldez/Dropbox (Encompass Bioscience)/cv_sandbox/data")
#write.table(oneMI_1, file='cv_final.phenos', sep=',', quote=F, col.names=T, row.names=F)

# Test 2 -> Keep entire cohort - just clean up for survival analysis
# Keep individuals age when diagnosed with a heart attack (myocardial infarction)
uk_oneMI <- subset(ukdata, !is.na(ukdata$mi_dx1)) # 11,533 individuals

#get a subset of 50,000 individuals w/o mi_dx1
uk_noMI <- subset(ukdata, is.na(ukdata$mi_dx1))
uk_noMI_50 <- uk_noMI[sample(nrow(uk_noMI), 50000), ]
uk_clean <- rbind(uk_oneMI, uk_noMI_50)

# parse out eids
eids_uk_clean <- data.frame(uk_clean$eid)
write.table(eids_uk_clean, "eids_onMI_50.txt", row.names = FALSE, col.names = FALSE)

# Add TRUE if have 1-heart attack, 2-Angina, 3-Stroke or FALSE if -7 NONE
uk_clean$heart_problem = uk_clean$heart_dx11 %in% c(1,2,3)  
uk_clean$heart_problem[uk_clean$heart_dx11=="-3"] = NA

# Select a few covariates from subset with one MI to see if they are significant
uk_clean_sub <- data.frame(uk_clean$eid, uk_clean$age1, uk_clean$sex, uk_clean$race1, uk_clean$mi_dx1, uk_clean$mi_dx2, uk_clean$bmi1,
                           uk_clean$age_death1, uk_clean$sys_bp10, uk_clean$smoke1, uk_clean$weight1, uk_clean$heart_dx11, uk_clean$age_angina_dx1,
                           uk_clean$stroke_dx1, uk_clean$heart_problem)
colnames(uk_clean_sub) <- c("iid", "age", "sex", "race", "mi_dx1", "mi_dx2", "bmi1", "age_death1", "sys_bp", "smoke_status", "weight1", "heart_dx11", 
                       "age_angina1", "stroke_dx1", "heart_problem")

# Let's try a survival model
uk_clean_sub$event_age = pmin(uk_clean_sub$age, uk_clean_sub$age_angina1, uk_clean_sub$mi_dx1, uk_clean_sub$stroke_dx1, na.rm = T)
uk_clean_sub$event_age[uk_clean_sub$event_age < 0] <- NA
hist(uk_clean_sub$event_age)
uk_clean_sub$had_event = !(is.na(uk_clean_sub$age_angina1) & is.na(uk_clean_sub$mi_dx1) & is.na(uk_clean_sub$stroke_dx1))

# Create a survival object
uk_clean_sub$survival <- Surv(uk_clean_sub$event_age, uk_clean_sub$had_event)

# Create predictor with sex, bmi, race, sys_bp, smoke_status
res.cox <- coxph(survival ~ sex + bmi1 + age + smoke_status + weight1, data=uk_clean_sub)
summary(res.cox)

setwd("/Users/crystalvaldez/Dropbox (Encompass Bioscience)/cv_sandbox/data")
write.table(uk_clean_sub, file='uk_clean_sub.phenos', sep=',', quote=F, col.names=T, row.names=F)


##############
# split eids into europeans, african-americans and indians
europeans <- data.frame(subset(uk_clean_sub, race%in%c("1001", "1002", "1003"))$iid)
aa <- data.frame(subset(uk_clean_sub, race%in%c("2001", "2002", "4002", "2004"))$iid)
indians <- data.frame(subset(uk_clean_sub,race%in%c("3001","3002","3003", "3004"))$iid)

write.table(europeans, file='~/Dropbox (Encompass Bioscience)/cv_sandbox/data/eur_eids.txt', sep=',', quote=F, col.names=F, row.names=F)
write.table(aa, file='~/Dropbox (Encompass Bioscience)/cv_sandbox/data/aa_eids.txt', sep=',', quote=F, col.names=F, row.names=F)
write.table(indians, file='~/Dropbox (Encompass Bioscience)/cv_sandbox/data/ind_eids.txt', sep=',', quote=F, col.names=F, row.names=F)



