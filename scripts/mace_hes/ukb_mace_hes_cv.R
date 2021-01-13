## UK Biobank Coronary Events
library(dplyr)
library(ggplot2)
library(ukbtools) # Care this library is GPL2
library(survival)
library(rms)
library(reshape2)

hes_coronary <- read.csv("/Users/crystalvaldez/Dropbox (Encompass Bioscience)/Encompass Shared/UKbiobank/KM_curves/ukb_hes_mace.tsv", sep="\t")

# which eids are unique
length(unique(hes_coronary$eid)) # 20,854
uni <- data.frame(unique(hes_coronary$eid))
write.table(uni, file = 'ukb_hes_mace_rsids.csv', row.names = FALSE, quote = F, col.names = FALSE)

# number of individuals that were discharged due to death
sum(hes_coronary$dismeth==4, na.rm=T) # 273

# group by individuals (eids)
hes_coronary <- group_by(hes_coronary, eid)
# Data admitted to the hospital
hes_coronary$admidate = as.Date(hes_coronary$admidate, format="%Y-%m-%d")
# use with grouped data, one row per group --> have 5 columns: num of records, unique eids, difference in admit date, first time admitted
hes_coronary_admits <- summarise(hes_coronary, n_records=n(), n_uniques = length(unique(admidate)), n_diff = as.numeric(max(admidate)-min(admidate)), first_admit = min(admidate))
#write.table(unique(hes_coronary$eid), "coronary_patients.txt", quote=F, row.names=F, col.names=F)
hes_coronary_admits$first_admit = as.Date(hes_coronary_admits$first_admit, format="%Y-%m-%d")

# Analysis
head(hes_coronary)
unique(hes_coronary$diag_icd10)

# MACE icd10 codes:
mace_codes <- c("I210", "I211", "I212", "I213", "I214", "I219", "I21X", "I220", "I221", "I228", "I229", "I200", "I633", "I635", "I630", "I631", "I636")

# merge total hes_coronary to admits by row name (eid) 
hes_coronary2 <- merge(hes_coronary, hes_coronary_admits[,c("eid", 'first_admit','n_uniques')], by="eid")
# find individuals that have mace code in icd10
hes_coronary2$is_mace = hes_coronary$diag_icd10 %in% mace_codes
# summarise individuals that have their first event/admission to the hospital
hes_coronary_first_event <- summarise(hes_coronary, n_records=n(), n_uniques=length(unique(admidate)))

# summaise individuals that have >1 hosptials visits (not working)
hes_coronary_notfirst = subset(hes_coronary2, first_admit != admidate & is_mace)

#hes_coronary_notfirst <- subset(hes_coronary2, as.Date(first_admit) != as.Date(admidate))
#hes_coronary_notfirst <- subset(hes_coronary_notfirst, is_mace == TRUE)
#(old)hes_coronary_notfirst <- subset(hes_coronary2, first_admit != admidate & is_mace)

# calculate the follow-up time
hes_coronary_notfirst$follow_up = as.numeric(hes_coronary_notfirst$admidate - hes_coronary_notfirst$first_admit)

# within the first quartile, consider same event: (look at follow-up > 58 days)
hes_coronary_diff = subset(hes_coronary_notfirst, follow_up > 58)

# Get rate of outcomes - a bit high, but reasonable (5.9%)
sum(hes_coronary_diff$follow_up < 365)

# Set up survival data
ukdata <- read.csv("ukb11102_stroke_mi.csv")
colnames(ukdata) = c("eid", "sex", "birth_year","birth_month", "date1", "date2", "date3", "sys_bp10", "sys_bp11", "sys_bp20", "sys_bp21", "sys_bp30", "sys_bp31","dia_bp10", "dia_bp11", "dia_bp20", "dia_bp21", "dia_bp30", "dia_bp31", "reason_nofollow", "date_nofollow", "age_bp_dx1", "age_bp_dx2", "age_bp_dx3", "age_angina_dx1", "age_angina_dx2", "age_angina_dx3", "year_uk1", "year_uk2", "year_uk3", "mi_dx1", "mi_dx2", "mi_dx3", "stroke_dx1", "stroke_dx2", "stroke_dx3", "heart_dx11", "heart_dx12", "heart_dx13", "heart_dx14", "heart_dx21", "heart_dx22", "heart_dx23", "heart_dx24", "heart_dx31", "heart_dx32", "heart_dx33", "heart_dx34", "height", "smoke1", "smoke2", "smoke3", "e_smoke1", "e_smoke2", "e_smoke3", "race1", "race2", "race3", "bmi1", "bmi2", "bmi3", "weight1", "weight2", "weight3", "age1", "age2", "age3", "death1", "death2", "death3", "age_death1", "age_death2", "age_death3" )

## Read in MI risk scores
path_rs<- "/Users/crystalvaldez/Dropbox (Encompass Bioscience)/ukbb_analysis/ukb_mace_he/output/"
prs_3 <- read.table(file.path(path_rs,"prs_5e-3.txt"), header = TRUE)
prs_4 <- read.table(file.path(path_rs,"prs_5e-4.txt"), header = TRUE)
prs_5 <- read.table(file.path(path_rs,"prs_5e-5.txt"), header = TRUE)
prs_6 <- read.table(file.path(path_rs,"prs_5e-6.txt"), header = TRUE)
prs_7 <- read.table(file.path(path_rs,"prs_5e-7.txt"), header = TRUE)
prs_8 <- read.table(file.path(path_rs,"prs_5e-8.txt"), header = TRUE)

names(prs_3) <- c("eid", "prs_5e-3")
names(prs_3) <- c("eid", "prs_5e-3")
names(prs_4) <- c("eid", "prs_5e-4")
names(prs_5) <- c("eid", "prs_5e-5")
names(prs_6) <- c("eid", "prs_5e-6")
names(prs_7) <- c("eid", "prs_5e-7")
names(prs_8) <- c("eid", "prs_5e-8")
mi_prs <- data.frame(eid=prs_3$eid, prs_3=prs_3$`prs_5e-3`, prs_4=prs_4$`prs_5e-4`, prs_5=prs_5$`prs_5e-5`,
                     prs_6=prs_6$`prs_5e-6`, prs_7=prs_7$`prs_5e-7`, prs_8=prs_8$`prs_5e-8`)

# match survey data to hospital episodes data by eid
ukdata1 = merge(ukdata, mi_prs, by="eid")
hes_coronary3 <- merge(hes_coronary2, ukdata1, by="eid")
hes_coronary3$date_nofollow = as.Date(hes_coronary3$date_nofollow, "%Y-%m-%d") 
hes_coronary_notfirst = subset(hes_coronary3, first_admit != admidate & is_mace)

### Build the survival object. Censor at next admission (MACE event) and at last follow up
hes_coronary3$follow_up = as.numeric(hes_coronary3$admidate - hes_coronary3$first_admit)
#hes_coronary3$follow_up = as.Date(hes_coronary3$admidate) - as.Date(hes_coronary3$first_admit)
hes_coronary3$censor_date = as.numeric(hes_coronary3$date_nofollow - hes_coronary3$first_admit)

# Maximum censoring time
hes_coronary3$max_date= max(hes_coronary3$admidate, na.rm=T)

# Follow up less than 6 weeks is not good, also remove non MACE
hes_coronary3$follow_up[hes_coronary3$follow_up <42 | !hes_coronary3$is_mace] = NA

hes_coronary3$age_at_start = as.numeric(format(hes_coronary3$first_admit, "%Y")) + 1/12*as.numeric(format(hes_coronary3$first_admit, "%m")) - hes_coronary3$birth_year - 1/12*hes_coronary3$birth_month

#prs_results = read.table("ukb_cardio_small_prs")
#names(prs_results) <- c("eid", "prs")

#prs2 = read.table("prs2_ukbcardio.txt")

#plot(prs2$V2, prs_results$prs)
#names(prs2) <- c("eid", "prs")
#prs3 = read.table("prs_impute.txt")
#names(prs3) <- c("eid", "prs") 

hes_coronary4 = subset(merge(hes_coronary3, prs_3, by="eid"))
hes_coronary4 = group_by(hes_coronary4, eid)
#hes_coronary4 = subset(merge(hes_coronary3, prs3, by="eid"), race1%in%c("1001", "1002", "1003"))
#hes_coronary4 = group_by(hes_coronary4, eid)

hes_coronary4$days_to_death = (pmin(hes_coronary4$age_death1, hes_coronary4$age_death2, hes_coronary4$age_death3, na.rm=T) - hes_coronary4$age_at_start)*365
hes_coronary4$mace_death = (hes_coronary4$death1 %in% mace_codes | hes_coronary4$death2 %in% mace_codes | hes_coronary4$death3 %in% mace_codes)
#hes_coronary4$days_to_death = (pmin(hes_coronary4$age_death1, hes_coronary4$age_death2, hes_coronary4$age_death3, na.rm=T) - hes_coronary4$age_at_start)*365
#hes_coronary4$mace_death = (hes_coronary4$death1 %in% mace_codes | hes_coronary4$death2 %in% mace_codes | hes_coronary4$death3 %in% mace_codes)

hes_coronary_surv <- summarise(hes_coronary4, start_date = min(admidate), censor_days=as.numeric(min(c(censor_date, follow_up, days_to_death, as.Date(max_date)), na.rm=T)), event = (sum(!is.na(follow_up))>0) | mace_death[1], prs=prs_3[1], age_at_start= age_at_start[1], sex=sex[1], smoke=smoke1[1], bmi=bmi1[1])
hes_coronary_surv_clean = subset(hes_coronary_surv, censor_days!='Inf' & censor_days > 14) # Require 2 weeks of data (per study protocol)
#hes_coronary_surv <- summarise(hes_coronary4, start_date = min(admidate), censor_days=as.numeric(min(c(censor_date, follow_up, days_to_death, max_date-admidate), na.rm=T)), event = (sum(!is.na(follow_up))>0) | mace_death[1], prs=prs[1], age_at_start= age_at_start[1], sex = sex[1], miprs0001=miprs0001[1], smoke=smoke1[1], bmi=bmi1[1] )
#hes_coronary_surv_clean = subset(hes_coronary_surv, censor_days!='Inf' & censor_days > 14) # Require 2 weeks of data (per study protocol)

# write out phenos to csv
path_output = "/Users/crystalvaldez/Dropbox (Encompass Bioscience)/ukbb_analysis/ukb_mace_he/output/ukb_mace_hes_phenos.csv"
write.table(hes_coronary4, path_output, sep = ',')

# Get quartile definitions
q1 = summary(hes_coronary_surv_clean$prs)[2]
q2 = summary(hes_coronary_surv_clean$prs)[5]
hes_coronary_surv_clean$high_low = hes_coronary_surv_clean$prs > q2
hes_coronary_surv_clean$high_low[q1 < hes_coronary_surv_clean$prs & hes_coronary_surv_clean$prs < q2] = NA
hes_coronary_surv_clean$top_q = hes_coronary_surv_clean$prs > q2
# add smoking and bmi
prs.cph = cph(Surv(censor_days, event) ~ prs + age_at_start + sex + smoke + bmi, data=hes_coronary_surv_clean, x=T, y=T, surv=T)
prs.coxph = coxph(Surv(hes_coronary_surv_clean$censor_days, hes_coronary_surv_clean$event)~ prs + age_at_start + sex + smoke + bmi, data=hes_coronary_surv_clean)
survplot(prs.cph, prs=c(q1, q2), age_at_start=59.88457, sex=0, smoke=0, bmi=28.4091, col=c("red", "blue"), ylim = c(0.8,1), xlim = c(0,1600))
###################### STOPPED HERE #######

#prs.cph = cph(Surv(censor_days, event)~ miprs0001 + age_at_start + sex +smoke + bmi, data=hes_coronary_surv_clean, x=T, y=T, surv=T)
#prs.coxph = coxph(Surv(hes_coronary_surv_clean$censor_days, hes_coronary_surv_clean$event)~ miprs0001 + age_at_start + sex + smoke + bmi, data=hes_coronary_surv_clean)

my_sims = rep(NA, 1000)
is_sig = rep(NA, 1000)
for (k in 1:1000) {
  sframe = hes_coronary_surv_clean[sample(nrow(hes_coronary_surv_clean), 6000),]
  thiscox = coxph(Surv(censor_days, event)~ prs + age_at_start + sex + smoke + bmi, data=sframe)
  my_sims[k] = thiscox$coefficients[1] 
  is_sig[k] = (thiscox$coefficients[1] - 2 * sqrt(thiscox$var[1,1] )) > 0
}
  
survplot(prs.cph, prs=c(q1, q2), age_at_start=59.88457, sex=0, smoke=0, bmi=28.4091, col=c("red", "blue"))

prs_km <- npsurv(Surv(censor_days, event) ~ high_low, data=hes_coronary_surv_clean)
survplot(prs_km, ylim=c(0.7,1), xlim=c(0,3200), xlab="Follow Up Time in Days", lty=1, col=c("firebrick", "steelblue"), col.fill=alpha(c("firebrick", "steelblue"), 0.3), label=F, n.risk=T)

prs_topq_km <- npsurv(Surv(censor_days, event)~ top_q, data=hes_coronary_surv_clean)
survplot(prs_topq_km, ylim=c(0.7,1), xlim=c(0,3200), xlab="Follow Up Time in Days", lty=1, col=c("firebrick", "steelblue"), col.fill=alpha(c("firebrick", "steelblue"), 0.3), label=F)

# Let's censor after 4.5 years
hes_coronary_surv_clean_c4 = hes_coronary_surv_clean

hes_coronary_surv_clean_c4$censor_days[hes_coronary_surv_clean_c4$censor_days > 1642] = 1642
hes_coronary_surv_clean_c4$event[hes_coronary_surv_clean_c4$censor_days > 1642] = FALSE

prs.cph_c4 = cph(Surv(censor_days, event)~ prs + age_at_start + sex +smoke + bmi, data=hes_coronary_surv_clean_c4, x=T, y=T, surv=T)
prs.coxph_c4 = coxph(Surv(censor_days, event)~ prs + age_at_start + sex + smoke + bmi, data=hes_coronary_surv_clean_c4)
# Model concordance, 0.581(old), now 0.603
survConcordance(Surv(censor_days, event)~ predict(prs.coxph_c4), data=na.omit(hes_coronary_surv_clean_c4[,c("censor_days", "event", "prs", "age_at_start", "sex", "smoke", "bmi")]))

bmi.coxph_c4 = coxph(Surv(censor_days, event)~ age_at_start + sex + smoke + bmi, data=hes_coronary_surv_clean_c4)
survConcordance(Surv(censor_days, event)~ predict(bmi.coxph_c4), data=na.omit(hes_coronary_surv_clean_c4[,c("censor_days", "event", "prs", "age_at_start", "sex", "smoke", "bmi")]))


prs.cph_c4_topq = cph(Surv(censor_days, event)~ top_q + age_at_start + sex +smoke + bmi, data=hes_coronary_surv_clean_c4, x=T, y=T, surv=T)

my_sims = rep(NA, 1000)
is_sig = rep(NA, 1000)
for (k in 1:1000) {
  sframe = hes_coronary_surv_clean[sample(nrow(hes_coronary_surv_clean_c4), 8000),]
  thiscox = coxph(Surv(censor_days, event) ~ prs + age_at_start + sex + smoke + bmi, data=sframe)
  my_sims[k] = thiscox$coefficients[1] 
  is_sig[k] = (thiscox$coefficients[1] - 2 * sqrt(thiscox$var[1,1] )) > 0
}


study_size = 500
enrich_events_4yr = rep(NA, 1000)
normal_events_4yr = rep(NA, 1000)

km_c4 <- npsurv(Surv(censor_days, event)~ top_q, data=hes_coronary_surv_clean_c4)
km_c4_survdiff = survdiff(Surv(censor_days, event)~ top_q, data=hes_coronary_surv_clean_c4)

km_c4_survfit = survfit(Surv(censor_days, event)~ top_q, data=hes_coronary_surv_clean_c4)

km_c4_survest = stepfun(km_c4_survfit$time, c(1, km_c4_survfit$surv))

pdf("acs_c4_topQ_km.pdf", height=5, width=6)
survplot(km_c4, ylim=c(0.8,1), xlim=c(0,1600), xlab="Follow Up Time in Days", lty=1, col=c("firebrick", "steelblue"), col.fill=alpha(c("firebrick", "steelblue"), 0.3), label=F, lwd=2)
title(main="Top quartile risk vs others")
dev.off()

# ****graph on current pitch deck slides*******
km_c4_highlow <- npsurv(Surv(censor_days, event)~ high_low, data=hes_coronary_surv_clean_c4)
pdf("acs_c4_highlow_km.pdf", height=5, width=6)
survplot(km_c4_highlow, ylim=c(0.8,1), xlim=c(0,1600), xlab="Follow Up Time in Days", lty=1, col=c("firebrick", "steelblue"), col.fill=alpha(c("firebrick", "steelblue"), 0.3), label=F, lwd=2)
title(main="Top quartile risk vs bottom quartile risk")
dev.off()

# Get survival fits for each data set
km_c4_survfit = survfit(Surv(censor_days, event)~ 1, data=hes_coronary_surv_clean_c4)

km_c4_survfit_topq <- survfit(Surv(censor_days, event)~ 1, data=subset(hes_coronary_surv_clean_c4, top_q==1))
km_c4_survfit_other<- survfit(Surv(censor_days, event)~ 1, data=subset(hes_coronary_surv_clean_c4, top_q==0))
km_c4_survfit_all<- survfit(Surv(censor_days, event)~ 1, data=hes_coronary_surv_clean_c4)

topq = subset(hes_coronary_surv_clean_c4, top_q==1)

km_c4_survest_all = stepfun(km_c4_survfit_all$time, c(1, km_c4_survfit_all$surv))
km_c4_survest_topq = stepfun(km_c4_survfit_topq$time, c(1, km_c4_survfit_topq$surv))


hes_coronary_surv_clean_c4_enrich = subset(hes_coronary_surv_clean_c4, miprs0001>q2)

# Play the bootstrap game for these rates

topq_rate = rep(NA, 1000)
all_rate = rep(NA, 1000)
top_q = subset(hes_coronary_surv_clean_c4, top_q==1)
for (i in 1:1000) {
  
  topq_r = top_q[sample(nrow(top_q), nrow(top_q), replace=T),]
  all_r = hes_coronary_surv_clean_c4[sample(nrow(hes_coronary_surv_clean_c4), nrow(hes_coronary_surv_clean_c4), replace=T),]
  
  trkm <- survfit(Surv(censor_days, event)~ 1, data=topq_r)
  arkm <- survfit(Surv(censor_days, event)~ 1, data=all_r)
  t_fxn = stepfun(trkm$time, c(1, trkm$surv))
  a_fxn = stepfun(arkm$time, c(1, arkm$surv))
  
  topq_rate[i] = t_fxn(365*4)
  all_rate[i] = a_fxn(365*4)

}

# Let's make a boxplot
exp_events = data.frame( enrich = 16000*(1-topq_rate), normal = 16000*(1-all_rate))
exp_melt = melt(exp_events)

ggplot(data=exp_melt) + geom_boxplot(aes(x=variable, y=value, fill=variable), show.legend = F) + theme_bw() + labs(x="Enrollment Scheme", y="Events") + scale_fill_manual(values=alpha(c(enrich="steelblue", normal="firebrick"), 0.5))
ggsave("ACS_trialEnrichBootstrap.pdf", height=5, width=3)

enrich_frame = data.frame(var=rep(c("enrich","normal"), times=c(1000,1000)), events=c(enrich_events_4yr, normal_events_4yr))
ggplot(data=enrich_frame) + geom_boxplot(aes(x=var, y=events)) + theme_bw() + labs(y="4 Year Events", x="Sampling Scheme")
ggsave("ACS_trialEnrichmentEventSims.pdf", height=4, width=3)



