#install necessary packages
install.packages("dplyr")
library("dplyr")
library(ggplot2)

#import data (change this each run)
cyto<-read.csv("MyExpt_MDST8Cytoplasm.csv")
nuclei<-read.csv("MyExpt_MDST8Nuclei.csv")

#Run the following script after upload of the above files

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

#summarize data for each condition in nuclei

nucleimeanfam83h <- data_summary(nuclei, varname="Intensity_MeanIntensity_fam83h", 
                                 groupnames="Metadata_drug")

nucleimeanCK1A <- data_summary(nuclei, varname="Intensity_MeanIntensity_CK1A", 
                               groupnames="Metadata_drug")

#graph mean intensity of CK1A in nucleus across all nuclei per condition 
ggplot(data=nucleimeanCK1A, aes(x=Metadata_drug, y=Intensity_MeanIntensity_CK1A)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=Intensity_MeanIntensity_CK1A-sd, ymax=Intensity_MeanIntensity_CK1A+sd), width=.2,
                position=position_dodge(.9)) + 
  ggtitle("CK1A nuclei mean intensity")

#graph mean intensity fam83h
ggplot(data=nucleimeanfam83h, aes(x=Metadata_drug, y=Intensity_MeanIntensity_fam83h)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=Intensity_MeanIntensity_fam83h-sd, ymax=Intensity_MeanIntensity_fam83h+sd), width=.2,
                position=position_dodge(.9)) +
  ggtitle("FAM83H nuclei mean intensity")

#repeat above for cytoplasm
#summarize
cytomeanck1a <- data_summary(cyto, varname="Intensity_MeanIntensity_CK1A", 
                             groupnames="Metadata_drug")
cytomeanfam83h <- data_summary(cyto, varname="Intensity_MeanIntensity_fam83h", 
                               groupnames="Metadata_drug")

#graph mean
ggplot(data=cytomeanck1a, aes(x=Metadata_drug, y=Intensity_MeanIntensity_CK1A)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=Intensity_MeanIntensity_CK1A-sd, ymax=Intensity_MeanIntensity_CK1A+sd), width=.2,
                position=position_dodge(.9)) +
  ggtitle("CK1A mean intensity cytoplasm")

ggplot(data=cytomeanfam83h, aes(x=Metadata_drug, y=Intensity_MeanIntensity_fam83h)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=Intensity_MeanIntensity_fam83h-sd, ymax=Intensity_MeanIntensity_fam83h+sd), width=.2,
                position=position_dodge(.9)) +
  ggtitle("FAM83H mean intensity cytoplasm")

#merging data and calculating nucleus:cyto intensity of ck1a
nucleimeanCK1Arenamed<-dplyr::rename(nucleimeanCK1A, "Mean_CK1A_Nucleus"="Intensity_MeanIntensity_CK1A")
nucleimeanCK1Arenamed<-dplyr::rename(nucleimeanCK1Arenamed, "sd_nuclei"="sd")
cytomeanck1arenamed<-dplyr::rename(cytomeanck1a, "Mean_CK1A_cyto"="Intensity_MeanIntensity_CK1A")
cytomeanck1arenamed<-dplyr::rename(cytomeanck1arenamed, "sd_cyto"="sd")

mergenucleicyto<-merge(nucleimeanCK1Arenamed,cytomeanck1arenamed)

#create a new column for nuclei to cyto ratio ck1a
mergenucleicyto2<-mutate(mergenucleicyto, nuclei_cyto=Mean_CK1A_Nucleus/Mean_CK1A_cyto)

#graph nuclei to cyto ratios ck1a
ggplot(data=mergenucleicyto2, aes(x=Metadata_drug, y=nuclei_cyto)) + 
  geom_bar(stat='identity', fill="blue") + 
  xlab("Condition") + 
  ylab("NucleitoCyto_Ratio") +
  ggtitle("Nucleus to Cytoplasm Ratio of CK1A")

#repeat above steps for fam83h

nucleimeanFAM83Hrenamed<-dplyr::rename(nucleimeanfam83h, "Mean_FAM83H_Nucleus"="Intensity_MeanIntensity_fam83h")
nucleimeanFAM83Hrenamed<-dplyr::rename(nucleimeanFAM83Hrenamed, "sd_nuclei"="sd")
cytomeanFAM83Hrenamed<-dplyr::rename(cytomeanfam83h, "Mean_FAM83H_cyto"="Intensity_MeanIntensity_fam83h")
cytomeanFAM83Hrenamed<-dplyr::rename(cytomeanFAM83Hrenamed, "sd_cyto"="sd")

mergenucleicytoFAM83H<-merge(nucleimeanFAM83Hrenamed,cytomeanFAM83Hrenamed)

#create a new column for nuclei to cyto ratio ck1a
mergenucleicytoFAM83H2<-mutate(mergenucleicytoFAM83H, nuclei_cyto=Mean_FAM83H_Nucleus/Mean_FAM83H_cyto)
View(mergenucleicytoFAM83H2)

#graph nuclei to cyto ratios ck1a
ggplot(data=mergenucleicytoFAM83H2, aes(x=Metadata_drug, y=nuclei_cyto)) + 
  geom_bar(stat='identity', fill="blue") + 
  xlab("Condition") + 
  ylab("NucleitoCyto_Ratio") +
  ggtitle("Nucleus to Cytoplasm Ratio of FAM83H")

#can be used to check summary means and sd calculated by function data_summary above 
nuclei %>%
  dplyr::group_by(Metadata_drug) %>%
  dplyr::summarise(
    count = n(),
    mean_Intensity = mean(Intensity_MeanIntensity_fam83h), 
    sd_Intensity = sd(Intensity_MeanIntensity_fam83h)
  )



