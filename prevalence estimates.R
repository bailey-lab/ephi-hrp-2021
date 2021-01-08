# E Reichert
# December 2020


##Background##

#This code reads in the cleaned full clinical dataset (n = 12,572) and the PCR dataset for which pfhrp2/3
#could be called successfully (n = 610). The proportion of those with a discordant RDT profile among all PF positives,
#both overall and by region, is then weighted by the concordance between the discordant RDT and a pfhrp2- PCR result
#(overall based on the assumption that PCR/RDT concordance will not vary by region, and to avoid small sample size bias).

#This approach estimates the true proportion of PF infections with pfhrp2- **that are leading to a false negative RDT**.
# (this excludes infections with sufficient cross-reactive pfhrp3 to be identified by HRP2-based RDTs, and therefore
# the true proportion of pfhrp2- parasites circulating may be higher)

library(dplyr)
library(asbio)
library(readr)

#read in our data
meta <- read_csv("Meta_full_clean.csv")
hrp23_pcr <- read_csv("Full_PCR_Clinical_clean.csv") %>%
  filter(pfldh_pcr_parasitemia_mean > 100)

#get proportion of Pf-positives with discordant RDT profile

#overall
RDT_overall <- meta %>%
  filter(rdt_pf_pos == 1) %>%
  mutate(discordant_RDT = rdt_discordant == 2) %>%
  group_by(discordant_RDT) %>%
  summarise(n())

disc_overall <- RDT_overall[[2,2]]/(RDT_overall[[1,2]]+RDT_overall[[2,2]])
#13.1% of Pf infections have discordant RDT overall

#by region
RDT_reg <- meta %>%
  filter(rdt_pf_pos == 1) %>%
  mutate(discordant_RDT = rdt_discordant == 2) %>%
  group_by(discordant_RDT, Region) %>%
  summarise(n())

disc_AM <- RDT_reg[[4,3]]/(RDT_reg[[1,3]]+RDT_reg[[4,3]])
#15.2% in Amhara

disc_GA <- RDT_reg[[5,3]]/(RDT_reg[[2,3]]+RDT_reg[[5,3]])
#1.6% in Gambella

disc_TG <- RDT_reg[[6,3]]/(RDT_reg[[3,3]]+RDT_reg[[6,3]])
#20.4% in Tigray

# get proportion of discordant RDTs with pfhrp2- PCR result
# here, we have a smaller subset

PCR_overall <- hrp23_pcr %>%
  filter(rdt_discordant ==2) %>%
  group_by(hrp2_neg_final_f) %>%
  summarise(n())

pcr_rdt_concord <- PCR_overall[[1,2]]/(sum(PCR_overall$`n()`))
#73.1% concordance between discordant RDT profile and pfhrp2- PCR call

## APPROACH #1

# The asbio function ci.impt estimates the product of two proportions and associated 95% CIs
# here, we want to estimate the product of (proportion of those with discordant RDT x 
# proportion with discordant RDT that have pfhrp2- PCR result), overall and by region

#Overall across 3 regions
ci.impt(RDT_overall[[2,2]], (RDT_overall[[1,2]]+RDT_overall[[2,2]]), PCR_overall[[1,2]], sum(PCR_overall$`n()`), conf = 0.95)
# 9.6% prevalence estimate (95% CI 8.4-10.9)

#for Amhara
ci.impt(RDT_reg[[4,3]], (RDT_reg[[1,3]]+RDT_reg[[4,3]]), PCR_overall[[1,2]], sum(PCR_overall$`n()`), conf = 0.95)
# 11.1% prevalence estimate (9.5-13.0)

#for Gambella
ci.impt(RDT_reg[[5,3]], (RDT_reg[[2,3]]+RDT_reg[[5,3]]), PCR_overall[[1,2]], sum(PCR_overall$`n()`), conf = 0.95)
# 1.2% prevalence estimate (0.7-2.1)

#for Tigray
ci.impt(RDT_reg[[6,3]], (RDT_reg[[3,3]]+RDT_reg[[6,3]]), PCR_overall[[1,2]], sum(PCR_overall$`n()`), conf = 0.95)
# 14.9% prevalence estimate (12.5-17.7)

## APPROACH #2

# To validate the estimates produced above, going to calculate the product of the 2 proportions directly
# and then use bootstrapping statistics (n = 1000) to generate 95% CIs and compare

#calculate product of proportions
# note: all of these estimates match those produced from ci.impt above
prev_overall <- disc_overall*pcr_rdt_concord
prev_overall_N <- prev_overall * (RDT_overall[[1,2]]+RDT_overall[[2,2]])

prev_AM <- disc_AM*pcr_rdt_concord
prev_AM_N <- prev_AM * (RDT_reg[[1,3]]+RDT_reg[[4,3]])

prev_GA <- disc_GA*pcr_rdt_concord
prev_GA_N <- prev_GA * (RDT_reg[[2,3]]+RDT_reg[[5,3]])

prev_TG <- disc_TG*pcr_rdt_concord
prev_TG_N <- prev_TG * (RDT_reg[[3,3]]+RDT_reg[[6,3]])

#write function
#source: http://pages.stat.wisc.edu/~larget/stat302/chap3.pdf
boot.mean = function(x,B,binwidth=NULL) {
  n = length(x)
  boot.samples = matrix( sample(x,size=n*B,replace=TRUE), B, n)
  boot.statistics = apply(boot.samples,1,mean)
  se = sd(boot.statistics)
  require(ggplot2)
  if ( is.null(binwidth) )
    binwidth = diff(range(boot.statistics))/30
  p = ggplot(data.frame(x=boot.statistics),aes(x=x)) +
    geom_histogram(aes(y=..density..),binwidth=binwidth) + geom_density(color="red")
  plot(p)
  interval = mean(x) + c(-1,1)*2*se
  print( interval )
  return( list(boot.statistics = boot.statistics, interval=interval, se=se, plot=p) )
}

#boostrap CIs overall
overall = c(rep(1, prev_overall_N), rep(0, (RDT_overall[[1,2]]+RDT_overall[[2,2]])-prev_overall_N))
overall.boot = boot.mean(overall, 1000, binwidth = 0.005)

#for Amhara
amhara = c(rep(1, prev_AM_N), rep(0, ((RDT_reg[[1,3]]+RDT_reg[[4,3]])-prev_AM_N)))
amhara.boot = boot.mean(amhara, 1000, binwidth = 0.005)

#for Gambella
gambella = c(rep(1, prev_GA_N), rep(0, ((RDT_reg[[2,3]]+RDT_reg[[5,3]])-prev_GA_N)))
gambella.boot = boot.mean(gambella, 1000, binwidth = 0.005)

#for Tigray
tigray = c(rep(1, prev_TG_N), rep(0, ((RDT_reg[[3,3]]+RDT_reg[[6,3]])-prev_TG_N)))
tigray.boot = boot.mean(tigray, 1000, binwidth = 0.005)

#Compared to CI's produced by ci.impt above in Approach #1, all are incredibly similar (within 0.02%)



