library(readxl)

## elena e. giorgi
## fred hutchinson, seattle, wa
## egiorgi@fredhutch.org

setwd("______________")

mydata <- as.data.frame(read_excel("BRAVE_under5_vax_DE-ID_dataset_080823.xlsx", sheet="BRAVE_under5_vax_dataset", range="A1:BE92"))
dim(mydata)

### no variation in 7 variables, removing them below
apply(mydata[,8:18], MARGIN=2, table)
mydata <- mydata[,-(8:18)]

table(mydata$serum_collected)
# N  Y 
# 12 79 

table(mydata$exclude)
#  N  Y 
# 75 16 

table(mydata$exclude[which(mydata$serum_collected=="N")])
#  Y 
# 12 

mydata <- mydata[-which(mydata$serum_collected=="N"),]
dim(mydata)
# [1] 79 41

mydata[which(mydata$exclude=="Y"), c("prior_infection", "prior_pcr", "vax_intervals_ok", "np_elisa")]
#    prior_infection prior_pcr vax_intervals_ok  np_elisa
# 60              NA        NA                N Equivocal
# 64              NA        NA                N  Negative
# 78              NA        NA                N  Positive
# 90              NA        NA                Y Equivocal

table(mydata$vax_intervals_ok)
#  N  Y 
#  3 76 

mydata <- mydata[-which(mydata$exclude=="Y"), ]

table(mydata$np_elisa)
# Equivocal  Negative  Positive 
#         2        54        19 

mydata[which(mydata$np_elisa=="Equivocal"), c("prior_infection", "prior_pcr", "vax_intervals_ok")]
#   prior_infection prior_pcr vax_intervals_ok
# 34               Y         Y                Y
# 85               Y         Y                Y

aggregate(vax_type ~ age_cat+sex, data=mydata, table)
#       age_cat sex vax_type.Moderna vax_type.Pfizer
# 24-59 months   F                7              11
#  6-23 months   F                9               5
# 24-59 months   M               15               7
#  6-23 months   M                9              12
 
 aggregate(vax_type ~ race+sex, data=mydata, table)
#                       race sex vax_type
#        Hispanic or Latino   F        1
#                     White   F   16, 15
#                     Asian   M     4, 2
# Black or African-American   M        1
#        Hispanic or Latino   M     1, 2
#                     Other   M        2
#                     White   M   17, 14

 aggregate(vax_type ~ age_cat, data=mydata, table)
#       age_cat vax_type.Moderna vax_type.Pfizer
# 1 24-59 months               22              18
# 2  6-23 months               18              17

 aggregate(vax_type ~ race+age_cat, data=mydata, table)
#                       race      age_cat vax_type
#                     Asian 24-59 months        2
# Black or African-American 24-59 months        1
#        Hispanic or Latino 24-59 months        1
#                     Other 24-59 months        1
#                     White 24-59 months   19, 16
#                     Asian  6-23 months     2, 2
#        Hispanic or Latino  6-23 months     1, 2
#                     Other  6-23 months        1
#                     White  6-23 months   14, 13
 
 table(mydata$race)
#                    Asian Black or African-American        Hispanic or Latino                     Other                     White 
#                        6                         1                         4                         2                        62 
 
 table(mydata$vax_type)
# Moderna  Pfizer 
#     40      35 
 
 table(mydata$age_cat)
# 24-59 months  6-23 months 
#          40           35  

 table(mydata$sex)
#  F  M 
# 32 43 

 summary(mydata$between_dose1_2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  20.00   22.50   28.00   28.88   31.00   55.00 

### Check that all NAs below are Moderna
summary(mydata$between_dose2_3)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NAs 
#  55.00   56.50   62.00   64.51   70.50   86.00      40 
  
 table(mydata$vax_type[which(is.na(mydata$between_dose2_3)==T)])
# Moderna 
#     40 

summary(mydata$timing_dose_1)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -155.00 -118.00  -87.00  -93.42  -68.00  -54.00 
 
 summary(mydata$timing_dose_2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -127.00  -96.50  -48.00  -64.55  -34.00  -26.00 
 
 summary(mydata$timing_dose_3)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -45.00  -37.00  -34.00  -33.86  -31.00  -21.00      40 
 
 summary(mydata$timing_covid)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -568.00 -226.00 -164.50 -175.50  -83.75  -12.00      51 
 
 table(mydata$prior_infection)
#  N  Y 
# 46 29 

### all prior infections were before 1-month collection
summary(mydata$timing_covid[which(mydata$prior_infection != "N")])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -568.00 -226.00 -164.50 -175.50  -83.75  -12.00       5 
 
 table(mydata$prior_pcr)
# N  Y 
# 51 24 
 
table(mydata$np_elisa)
# Equivocal  Negative  Positive 
#        2        54        19 
 
 aggregate(prior_pcr ~ vax_type, mydata, table)
#  vax_type prior_pcr.N prior_pcr.Y
#  Moderna          30          10
#   Pfizer          21          14
 
 aggregate(np_elisa ~ vax_type, mydata, table)
#  vax_type np_elisa.Equivocal np_elisa.Negative np_elisa.Positive
#  Moderna                  0                29                11
#   Pfizer                  2                25                 8
 
 aggregate(prior_infection ~ vax_type, mydata, table)
#  vax_type prior_infection.N prior_infection.Y
#  Moderna                28                  12
#   Pfizer                18                  17

aggregate(prior_infection ~ prior_pcr + np_elisa, mydata, table)
#   prior_pcr  np_elisa prior_infection
# 1         Y Equivocal               2
# 2         N  Negative              46
# 3         Y  Negative               8
# 4         N  Positive               5
# 5         Y  Positive              14

# there are 8 cases that are PCR+ and Elisa-
# there are 5 cases that are PCR- and Elisa+
### (binomial p=0.03)
table(mydata$vax_type[which(mydata$np_elisa == "Negative" & mydata$prior_pcr == "Y")])
# Moderna  Pfizer 
#      1       7 
  
table(mydata$vax_type[which(mydata$np_elisa == "Positive" & mydata$prior_pcr == "N")])
# Moderna  Pfizer 
#      2       3 

### substitute "<20" with 10 for geometric mean computations
for(i in 41:46){ mydata[,i] <- as.numeric(gsub("<20", "10", mydata[,i])) }

write.table(mydata, file="BRAVE_under5_vax_dataset_080823.cleaned.csv", sep=",", row.names=F, quote=F)