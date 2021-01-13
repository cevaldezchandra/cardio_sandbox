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

# Test 1 -> Keep individuals with at least one heart attack (myocardial infarction)
uk_oneMI <- subset(ukdata, !is.na(ukdata$mi_dx1)) # 11,533 individuals

# Parse out individuals that do not know age they were diagnosed with first MI (-1) or did not answer (-3)
uk_oneMI <- uk_oneMI[!uk_oneMI$mi_dx1 %in% c(-1,-3),] # 11,186 individuals

# Parse out individuals that did not answer their smoking status (-3)
uk_oneMI <- uk_oneMI[!uk_oneMI$smoke1 %in% c(-3),] # 11,109

# Select a few covariates from subset with one MI to see if they are significant
oneMI_1 <- data.frame(uk_oneMI$eid, uk_oneMI$age1, uk_oneMI$sex, uk_oneMI$race1, uk_oneMI$mi_dx1, uk_oneMI$mi_dx2, uk_oneMI$bmi1,
                      uk_oneMI$age_death1, uk_oneMI$sys_bp10, uk_oneMI$smoke1, uk_oneMI$weight1)
colnames(oneMI_1) <- c("iid", "age", "sex", "race", "mi_dx1", "mi_dx2", "bmi1", "age_death1", "sys_bp", "smoke_status", "weight1")

# Who had second heart attack
oneMI_1$status <- 0
oneMI_1[!is.na(oneMI_1$mi_dx2),]$status <- 1

# Calculate time to second heart attack
oneMI_1$time_mi2 <- abs((oneMI_1$mi_dx2 - oneMI_1$mi_dx1) * 12)

# Create a survival object
oneMI_1$survival <- with(oneMI_1, Surv(mi_dx1, status==1))

# Create predictor with sex, bmi, race, sys_bp, smoke_status
res.cox <- coxph(survival ~ sex + bmi1 + age + smoke_status + weight1, data=oneMI_1)
summary(res.cox)

