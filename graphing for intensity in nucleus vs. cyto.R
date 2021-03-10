nuclei<-read.csv("Nuclei.csv")
cyto<-read.csv("Cytoplasm.csv")

install.packages("dplyr")
library("dplyr")

library(ggplot2)

#summarize data
nuclei summarized <- 
  nuclei %>%
  group_by(Metadata_condition) %>%
  summarise(
    count = n(),
    mean_IntegratedIntensity = mean(Intensity_IntegratedIntensity_AHR),
    mean_Intensity = mean(Intensity_MeanIntensity_AHR), 
    sd_IntegratedIntensity = sd(Intensity_IntegratedIntensity_AHR), 
    sd_Intensity = sd(Intensity_MeanIntensity_AHR)
  )

#create function to summarize data (from R online forum)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#summarize data 
nucleiintegrateddf <- data_summary(nuclei, varname="Intensity_IntegratedIntensity_AHR", 
                    groupnames="Metadata_condition")

nucleimeandf <- data_summary(nuclei, varname="Intensity_MeanIntensity_AHR", 
                             groupnames="Metadata_condition")

#graph integrated intensity of protein in the nucleus across all nuclei per condition 
ggplot(data=nucleiintegrateddf, aes(x=Metadata_condition, y=Intensity_IntegratedIntensity_AHR)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=Intensity_IntegratedIntensity_AHR-sd, ymax=Intensity_IntegratedIntensity_AHR+sd), width=.2,
                position=position_dodge(.9)) 

#graph mean intensity of protein in nucleus across all nuclei per condition 
ggplot(data=nucleimeandf, aes(x=Metadata_condition, y=Intensity_MeanIntensity_AHR)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=Intensity_MeanIntensity_AHR-sd, ymax=Intensity_MeanIntensity_AHR+sd), width=.2,
                position=position_dodge(.9)) +
  ggtitle("Mean AHR Intensity in Nucleus")

#repeat above for cytoplasm
#summarize
cytointegrateddf <- data_summary(cyto, varname="Intensity_IntegratedIntensity_AHR", 
                                   groupnames="Metadata_condition")

cytomeandf <- data_summary(cyto, varname="Intensity_MeanIntensity_AHR", 
                             groupnames="Metadata_condition")

#graph integrated
ggplot(data=cytointegrateddf, aes(x=Metadata_condition, y=Intensity_IntegratedIntensity_AHR)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=Intensity_IntegratedIntensity_AHR-sd, ymax=Intensity_IntegratedIntensity_AHR+sd), width=.2,
                position=position_dodge(.9)) +
  ggtitle("Integrated AHR Intensity in Cytoplasm")

#graph mean
ggplot(data=cytomeandf, aes(x=Metadata_condition, y=Intensity_MeanIntensity_AHR)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=Intensity_MeanIntensity_AHR-sd, ymax=Intensity_MeanIntensity_AHR+sd), width=.2,
                position=position_dodge(.9)) +
  ggtitle("Mean AHR Intensity in Cytoplasm")

#merging data and calculating nucleus:cyto intensity
nucleimeanrenamed<-dplyr::rename(nucleimeandf, "Mean_AHR_nuclei"="Intensity_MeanIntensity_AHR")
nucleimeanrenamed<-dplyr::rename(nucleimeandf, "sd_nuclei"="sd")
cytomeanrenamed<-dplyr::rename(cytomeandf, "Mean_AHR_cyto"="Intensity_MeanIntensity_AHR")
cytomeanrenamed<-dplyr::rename(cytomeanrenamed, "sd_cyto"="sd")

mergenucleicyto<-merge(nucleimeanrenamed,cytomeanrenamed)

#create a new column for nuclei to cyto ratio
mergenucleicyto2<-dplyr::mutate(mergenucleicyto, nuclei_cyto=Mean_AHR_nuclei/Mean_AHR_cyto)
View(mergenucleicyto2)
#graph nuclei to cyto ratios 
ggplot(data=mergenucleicyto2, aes(x=Metadata_condition, y=nuclei_cyto)) + 
  geom_bar(stat='identity', fill="blue") + 
  xlab("Condition") + 
  ylab("NucleitoCyto_Ratio") +
  ggtitle("Nuclei to Cytoplasm Ratio of AHR")
 

