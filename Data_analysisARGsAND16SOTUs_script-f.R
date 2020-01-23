###### --------------------------------------------------------- ######
#                                                                     #
#                        ARGs and MGEs                                #
#                                                                     #
###### --------------------------------------------------------- ######

######################################################################
#         R-Script for analyzing Wafergen qPCR array data            #
#            Author: Johanna Muurinen, data-analysis by              #
#                         Antti Karkman                              #
# (https://www.researchgate.net/publication/313181737_Data_analysis) #                                     
#                    was used as a template.                         #
#                                                                    #
#    Written at Dr. Tim Johnson Lab, Dept. of Animal Sciences,       #
#               Purdue University                                    #
######################################################################

#-------------- Housekeeping ------------------# 

#R-environment is cleared

rm(list=ls())

#Working directory

setwd("~/Library/Blaah blaa/blaa~blaa~blaa/Wafergen")

#Libraries and packages
#Before libraries can be loaded, the packages needs to be installed. 
#They can be installed like this:
#install.packages("reshape2")

library(reshape2)
library(dplyr)
library(tidyr)
library(MASS)
library(vegan)
library(gplots)
library(ggplot2)
library(ggthemes)
library(colorspace)
library(colorRamps)
library(RColorBrewer)
library(gridExtra)
library(fitdistrplus)
library(logspline)
library(car)
library(lme4) 
library(multcomp)
library(scales)
#library(ggridges)
#library(viridis)
#library(corrplot)
#library(PerformanceAnalytics)
library(psych)




#-------------- Data & processing ------------------# 

#We define a function to read in each plate and calculate 
#the mean and SD from the three technical replicates. 
#If only 1/3 gives signal (Ct < 27), it is considered as false positive. 
#The function will be called Array_func, but can also be called 
#something else. Just be careful that you don’t overwrite some existing functions.

Array_func <- function(PATH) {
  #read in the chip
  chip <- read.table(PATH, header=TRUE, sep="\t")
  #change all Ct > 27 to NA
  chip$Ct[chip$Ct>27] <- NA
  #make a new data drame for the results
  outArray <- data.frame()
  #loop over the whole chip and calculate the mean and SD
  
  #There are 4 samples/chip and 384 different assays. 
  #The samples have 3 replcates so that rows 1-3 have replicates of Sample1
  #and the Ct values of assay called "16S new 2_2"
  #If your samples and replicates ar in different order, you need to either 
  #change the function accordingly or reorder your .txt file.
  #reordering can be done for example in excel (just save the file as .txt file)
  
  #so first we create an index over the .txt file in every 3 rows
  for (i in seq(1,nrow(chip), 3)) {
    #if 2/3 of Ct-values are NAs, write everything as NA
    #the Ct values are in the column 6 in the .txt -file
    if(is.na(chip[i,6]) && is.na(chip[i+1,6]) || is.na(chip[i+1,6]) && 
       is.na(chip[i+2,6]) || is.na(chip[i,6]) && is.na(chip[i+2,6])) {
      #sample names are in column 4 in the .txt -file and assay names in column 3
      outArray[(i+2)/3, 1] <- chip[i,4]
      outArray[(i+2)/3, 2] <- chip[i,3]
      outArray[(i+2)/3, 3] <- NA
      outArray[(i+2)/3, 4] <- NA
    }
    #If at least 2/3 of Ct values are <27, then the mean and sd are caluculated
    #Note that the loop uses the sample name of replicate 1, so in our data set 
    #the samples are named as Sample#-Rep1. We will deal with that later. 
    else {
      outArray[(i+2)/3, 1] <- chip[i,4]
      outArray[(i+2)/3, 2] <- chip[i,3]
      outArray[(i+2)/3, 3] <- mean(c(chip[i,6], chip[i+1,6], chip[i+2,6]), na.rm=TRUE)
      outArray[(i+2)/3, 4] <- sd(c(chip[i,6], chip[i+1,6], chip[i+2,6]), na.rm=TRUE)
    }
  }
  colnames(outArray) <- c("Sample", "Assay", "mean", "sd")
  return(outArray)
}

#Then we can use the new function to read in 
#all the text files. Do this for all the chips (.txt files) you have.
# Example: chip1 <- Arrray_func("PATH/TO/YOUR/RESULTS/chip1.txt")
#if your chip data is in your working directory like in our case, 
#you dont need to specify the path. 

chip1II<- Array_func("IIJohnson_Chip1_CorrectedNames_sorted.txt")
chip2II<- Array_func("IIJohnsonChip2_CorrectedNames_sorted.txt")
chip3II<- Array_func("IIJohnsonChip3_CorrectedNames_sorted.txt")
chip4II<- Array_func("IIJohnsonChip4_CorrectedNames_sorted.txt")
chip5IIre<- Array_func("IIJohnsonChip5rerun_sorted.txt")
chip6II<- Array_func("IIJohnsonChip6_CorrectedNames_sorted.txt")

#Then combine all chips to one data frame

all_chipsII <- rbind(chip1II, chip2II, chip3II, chip4II, chip5IIre, chip6II)

#Do we have all?
View(all_chipsII)
AssaysII<-levels(all_chipsII$Assay)
length(AssaysII)
#[1] 384, yes we do.

#-------------- Data cleaning ------------------# 

# --- Unspecific amplification ---#

#Sample 1 is the negative control, but in some assays there is 
#amplification because some primers are unspecific 
#and we cannot rule out contamination either.  
#Let's see how many "problems" we have.

Pos_negC<-subset(all_chipsII, Sample=="Sample1-Rep1" & mean!="NA")
View(Pos_negC)

#These are potential problematic assays that may need to be removed 
#from the data eventually (6 target genes and the old 16S). 
#Howver we can use the other (new) 16S in normalizing.
#"The new 16S" primers did not amplify products in Sample1. 
#Let's have a look what are the Ct values with these assays
#in all samples: 

pot.problems<-Pos_negC$Assay

#I will cerate a data frame for our problems:

Unspecific_assays<-subset(all_chipsII, Assay %in% pot.problems)
Unspecific_assays<-subset(Unspecific_assays, mean!="NA")

#how many we have?

dim(Unspecific_assays)
#[1] 79  4
#Lets organize the data frame;
Unspecific_assays<-Unspecific_assays[order(Unspecific_assays$Assay,Unspecific_assays$Sample),]
#and have a look
View(Unspecific_assays)

#There probably is another way to do this but with this fuction we can 
#put the number of problematic samples & assays in a data frame

unspec_func<-function(assays) {
  TEMP<-subset(all_chipsII,  mean!="NA" & Assay==assays)
  count<-nrow(TEMP)
  return(count)
}    
#does it work?
unspec_func ("fabK_1520")
#6

#Now for all:
tmp<-unique(Unspecific_assays$Assay)
tmp2<-lapply(tmp,unspec_func)
names(tmp2)<-tmp
tmp3<-melt(tmp2)
Unspecifics<-as.data.frame(tmp3)
print(Unspecifics)

#value                     L1
#1    24            16S old 1_1
#2     4          blaOXY-1_1118
#3    19               cmlV_911
#4    12              czcA_1536
#5     6              fabK_1520
#6    13 intI1F165_clinical_359
#7     1             tetPA_1507

#tetPA is actually positive only in Sample 1 (negative control)

#In order to solve the problem with unspecific amplification, I want to have the 
# mean Ct in other samples than negative control next to the mean in sample 1
#so I can compare the Ct values. 


Pos_negC$Sample_mean<-rep(1, nrow(Pos_negC))  
Pos_negC<- within(Pos_negC, Sample_mean[Assay=="16S old 1_1"]<-mean(subset(Unspecific_assays$mean, 
                                                                           Unspecific_assays$Sample!="Sample1-Rep1" & Unspecific_assays$Assay=="16S old 1_1")))
Pos_negC<- within(Pos_negC, Sample_mean[Assay=="blaOXY-1_1118"]<-mean(subset(Unspecific_assays$mean, 
                                                                             Unspecific_assays$Sample!="Sample1-Rep1" & Unspecific_assays$Assay=="blaOXY-1_1118")))
Pos_negC<- within(Pos_negC, Sample_mean[Assay=="cmlV_911"]<-mean(subset(Unspecific_assays$mean, 
                                                                        Unspecific_assays$Sample!="Sample1-Rep1" & Unspecific_assays$Assay=="cmlV_911")))
Pos_negC<- within(Pos_negC, Sample_mean[Assay=="czcA_1536"]<-mean(subset(Unspecific_assays$mean, 
                                                                         Unspecific_assays$Sample!="Sample1-Rep1" & Unspecific_assays$Assay=="czcA_1536")))
Pos_negC<- within(Pos_negC, Sample_mean[Assay=="fabK_1520"]<-mean(subset(Unspecific_assays$mean, 
                                                                         Unspecific_assays$Sample!="Sample1-Rep1" & Unspecific_assays$Assay=="fabK_1520")))
Pos_negC<- within(Pos_negC, Sample_mean[Assay=="intI1F165_clinical_359"]<-mean(subset(Unspecific_assays$mean, 
                                                                                      Unspecific_assays$Sample!="Sample1-Rep1" & Unspecific_assays$Assay=="intI1F165_clinical_359")))
Pos_negC<- within(Pos_negC, Sample_mean[Assay=="tetPA_1507"]<-mean(subset(Unspecific_assays$mean, 
                                                                          Unspecific_assays$Sample!="Sample1-Rep1" & Unspecific_assays$Assay=="tetPA_1507")))

View(Pos_negC)

# 16S old 1_1 is not going to cause problems, we can use the other 16S because it was
# positive in all samples except negative control.
# tetPA was actually positive only in Sample 1 (negative control), so that problem is also
# solved
# czcA_1536 had lower Ct values in negative control than samples, so they results of that assay will be
# removed. 
# blaOXY-1_1118, cmlV_911b, intI1F165_clinical_359 and fabK_1520 have lower Ct values in samples than in 
# the negative control. It is possible that there is unspecific amplification (primer-dimers)
# if template DNA is not present. But to take it into account that the unspecific amplification 
# might lower the Ct values of samples I will adjust the Ct values of these assays:

# I will subtract the Ct values of samples from the mean Ct value in the negative control. Then I'm 
# going to subtract the result from the cut-off Ct value 27 and use that result as the new Ct value.  

View(all_chipsII)

all_chipsII$help_col1<-rep(1, nrow(all_chipsII))  
all_chipsII<- within(all_chipsII, help_col1[Assay=="blaOXY-1_1118"]<-subset(Pos_negC$mean, 
 Pos_negC$Assay=="blaOXY-1_1118")-subset(all_chipsII$mean, all_chipsII$Assay=="blaOXY-1_1118"))
all_chipsII<- within(all_chipsII, help_col1[Assay=="cmlV_911"]<-subset(Pos_negC$mean, 
  Pos_negC$Assay=="cmlV_911")-subset(all_chipsII$mean, all_chipsII$Assay=="cmlV_911"))
all_chipsII<- within(all_chipsII, help_col1[Assay=="intI1F165_clinical_359"]<-subset(Pos_negC$mean, 
  Pos_negC$Assay=="intI1F165_clinical_359")-subset(all_chipsII$mean, all_chipsII$Assay=="intI1F165_clinical_359"))
all_chipsII<- within(all_chipsII, help_col1[Assay=="fabK_1520"]<-subset(Pos_negC$mean, 
  Pos_negC$Assay=="fabK_1520")-subset(all_chipsII$mean, all_chipsII$Assay=="fabK_1520"))


all_chipsII$help_col2<-rep("NA", nrow(all_chipsII))  
all_chipsII<- within(all_chipsII, help_col2[Assay=="blaOXY-1_1118"]<-27-subset(all_chipsII$help_col1, 
     all_chipsII$Assay=="blaOXY-1_1118"))
all_chipsII<- within(all_chipsII, help_col2[Assay=="cmlV_911"]<-27-subset(all_chipsII$help_col1, 
    all_chipsII$Assay=="cmlV_911"))
all_chipsII<- within(all_chipsII, help_col2[Assay=="intI1F165_clinical_359"]<-27-subset(all_chipsII$help_col1, 
      all_chipsII$Assay=="intI1F165_clinical_359"))
all_chipsII<- within(all_chipsII, help_col2[Assay=="fabK_1520"]<-27-subset(all_chipsII$help_col1, 
   all_chipsII$Assay=="fabK_1520"))
#
# I will change values that are larger than 27 or 27 to NA's 
all_chipsII$help_col2[all_chipsII$help_col2>=27]<-"NA"
# new column for new new mean values
all_chipsII$meanII<-all_chipsII$mean
# and new mean values for three problematic assays
all_chipsII<- within(all_chipsII, meanII[Assay=="blaOXY-1_1118"]<-subset(all_chipsII$help_col2, 
       all_chipsII$Assay=="blaOXY-1_1118"))
all_chipsII<- within(all_chipsII, meanII[Assay=="cmlV_911"]<-subset(all_chipsII$help_col2, 
  all_chipsII$Assay=="cmlV_911"))
all_chipsII<- within(all_chipsII, meanII[Assay=="intI1F165_clinical_359"]<-subset(all_chipsII$help_col2, 
   all_chipsII$Assay=="intI1F165_clinical_359"))
all_chipsII<- within(all_chipsII, meanII[Assay=="fabK_1520"]<-subset(all_chipsII$help_col2, 
    all_chipsII$Assay=="fabK_1520"))

#new mean values need to be changed to numeric
all_chipsII$meanII<-as.numeric(all_chipsII$meanII)
#Warning message:
# NAs introduced by coercion
#Yes, there are NA's

# Next I will remove the results of assay czcA_1536 that had smaller 
#Ct values in the negative control than in samples. I don't need to remove assays tetPA_1507
# and 16S old 1_1 because they will drop out along the way. 

toBeRemoved<-which(all_chipsII$Assay=="czcA_1536")
all_chipsII<-all_chipsII[-toBeRemoved,]
str(all_chipsII)


# --- Assay names, ARG annotation and normalization to 16S Ct's ---#

#There are ususally problems with some assay names, since characters 
# like ", ' and / do not work well when some packages or when data is imported into R,
# also wafergen system changes some characters in the assay names in the chip files.
#and if the text file is opened in "non-US Excel" along the way some characters might be different
#These can be changed for example using Excel (US version) or with some other text editing probram
#and saving the files as .txt files.
#After this the all the data that has assay names has to be imported again. 
# Checking these issues:

levels(all_chipsII$Assay)

#one example: 

which(all_chips$Assay== "A. baumannii<ca>_1535") 
# I  could not find any logical explanation why "_" was changed into 
#"<ca>_" in the chip files.
#So I think the easiest way to fix these is to write a table and open it in Excel,
#see what is the problem and change the name to something R supports. 
# (I took the Assays to a separate data frame):

#Assays<-all_chipsII$Assay
#write.table(Assays,"/Users/blaa/blaa/blaa/blaa~blaa~blaa/Wafergen/AssaysInR", sep="\t")

#and then changed the assay names to all the chip files
#After this the files need to be imported and processed again.

#When all data have been processed, read in the final ARG assay list
#make sure that the assay names in your data files and in the assay list match! 
Primerset2 <- read.table("Primerset2_0_tnpAs_fixed.txt", header=TRUE, sep="\t")

#Check that changes you made are ok:
which (Primerset2 $Assay== "A.baumannii_1535") 
subset(all_chipsII, Assay=="A.baumannii_1535")

# I still need to do some corrections I missed earlier:
levels(Primerset2$Mechanism)
#[1] "bl3_cpha"   "catB3"      "catB8"      "deactivate" "efflux"     "MGE"        "none"       "other"     
#[9] "protection" "regulator"  "sugE"       "unknown" 

levels(Primerset2$Classification)
#[1] "Aminoglycoside"  "Amphenicol"      "beta Lactam"     "Beta Lactam"     "Beta lactamase" 
#[6] "fluoroquinolone" "Fluoroquinolone" "housekeeping"    "Insertional"     "Insertional "   
#[11] "Integrase"       "MDR"             "MDR-chromo"      "MDR-mobile"      "MGE"            
#[16] "MLSB"            "other"           "Phenicol"        "Phenicol_MGE"    "plasmid"        
#[21] "Plasmid-inc"     "Plasmid-inc "    "Plasmid-rep"     "rifamycin"       "Sulfonamide"    
#[26] "taxanomic"       "tetracycline"    "Tetracycline"    "Tetracycline "   "Transposase"    
#[31] "trimethoprim"    "Vancomycin"

which(Primerset2$Mechanism=="sugE")
Primerset2$Mechanism[302]<-"efflux"

which(Primerset2$Mechanism=="bl3_cpha")

Primerset2[75,]

Primerset2$Mechanism[75]<-"deactivate"

which(Primerset2$Mechanism=="catB8")

Primerset2[120,]

Primerset2$Mechanism[120]<-"deactivate"

which(Primerset2$Mechanism=="catB3")

Primerset2[119,]

Primerset2$Mechanism[119]<-"deactivate"


which(Primerset2$Mechanism=="other")
#[1]   3 111 138 159 243 308

subset(Primerset2, Mechanism=="other")
#Assay Mechanism Classification
#3   A.baumannii_1535     other      taxanomic
#111        cadC_1575     other     MDR-mobile
#138         cro_1569     other            MGE
#159   EAE_05855_1570     other            MGE
#243   merA-marko_331     other            MDR
#308        terW_1572     other     MDR-mobile

which(Primerset2$Mechanism=="MGE") #OK, a lot...

Primerset2$Mechanism<-droplevels(Primerset2$Mechanism)

levels(Primerset2$Mechanism)

#[1] "deactivate" "efflux" "MGE"  "none"  "other" "protection" "regulator"  "unknown" 


levels(Primerset2$Classification)<-c(levels(Primerset2$Classification), 
                                     "Beta-lactam")

Primerset2$Classification[Primerset2$Classification=="beta Lactam"]<-"Beta-lactam"

Primerset2$Classification[Primerset2$Classification=="Beta Lactam"]<-"Beta-lactam"

Primerset2$Classification[Primerset2$Classification=="Beta lactamase"]<-"Beta-lactam"

Primerset2$Classification[Primerset2$Classification=="Insertional "]<-"Insertional"

Primerset2$Classification[Primerset2$Classification=="Insertional "]<-"Insertional"

Primerset2$Classification[Primerset2$Classification=="Amphenicol"]<-"Phenicol"

Primerset2$Classification[Primerset2$Classification=="Plasmid-inc "]<-"Plasmid-inc"

levels(Primerset2$Classification)<-c(levels(Primerset2$Classification), 
                                     "Tetracycline")

Primerset2$Classification[Primerset2$Classification=="Tetracycline "  ]<-"Tetracycline"  

Primerset2$Classification[Primerset2$Classification=="tetracycline"  ]<-"Tetracycline"  

levels(Primerset2$Classification)<-c(levels(Primerset2$Classification), 
                                     "Trimethoprim" )

Primerset2$Classification[Primerset2$Classification=="trimethoprim"]<-"Trimethoprim"  

Primerset2$Classification[Primerset2$Classification=="fluoroquinolone"]<-"Fluoroquinolone"  

levels(Primerset2$Classification)<-c(levels(Primerset2$Classification), 
                                     "Rifamycin" )

Primerset2$Classification[Primerset2$Classification=="rifamycin"]<-"Rifamycin"  

levels(Primerset2$Classification)<-c(levels(Primerset2$Classification), 
                                     "Taxonomic" )

Primerset2$Classification[Primerset2$Classification=="taxanomic"]<-"Taxonomic"  


Primerset2$Classification[Primerset2$Classification=="Phenicol_MGE"]<-"Phenicol"  

Primerset2$Classification[Primerset2$Classification=="MDR-chromo"]<-"MDR"  

Primerset2$Classification[Primerset2$Classification=="MDR-mobile"]<-"MDR"  

levels(Primerset2$Classification)<-c(levels(Primerset2$Classification), 
                                     "Plasmid" )

Primerset2$Classification[Primerset2$Classification=="plasmid"]<-"Plasmid"  

Primerset2$Classification[Primerset2$Classification=="Plasmid-inc"]<-"Plasmid"  

Primerset2$Classification[Primerset2$Classification=="Plasmid-rep"]<-"Plasmid"  

Primerset2$Classification<-droplevels(Primerset2$Classification)

levels(Primerset2$Classification)<-c(levels(Primerset2$Classification), 
                                     "Other" )

Primerset2$Classification[Primerset2$Classification=="other"]<-"Other"  

Primerset2$Classification<-droplevels(Primerset2$Classification)

levels(Primerset2$Classification)

#[1] "Aminoglycoside"  "Fluoroquinolone" "housekeeping"    "Insertional"     "Integrase"       "MDR"             "MGE"            
#[8] "MLSB"            "Phenicol"        "Sulfonamide"     "Tetracycline"    "Transposase"     "Vancomycin"      "Beta-lactam"    
#[15] "Trimethoprim"    "Rifamycin"       "Taxonomic"       "Plasmid"         "Other"


#Gene tetU is not a resistance gene actually.
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3370814/

#the corrected insertion sequence-tnpA names need to be corrected also 
# to the chips (because I didn't remember to do it and they are not 
#in the original publication by Stedtfeld et al).

levels(all_chipsII$Assay)

levels(all_chipsII$Assay)<-c(levels(all_chipsII$Assay),"tnpA-05/IS26_201",
    "tnpA-04/IS6100_202", "tnpA-02/IS4_203","tnpA-01/IS21_204", "tnpA-07/ISEcp1B_205",
         "tnpA-06/IS1216_206", "tnpA-03/IS6_207")

all_chipsII$Assay[all_chipsII$Assay=="tnpA_201"]<-"tnpA-05/IS26_201"  
all_chipsII$Assay[all_chipsII$Assay=="tnpA_202"]<-"tnpA-04/IS6100_202"  
all_chipsII$Assay[all_chipsII$Assay=="tnpA_203"]<-"tnpA-02/IS4_203" 
all_chipsII$Assay[all_chipsII$Assay=="tnpA_204"]<-"tnpA-01/IS21_204" 
all_chipsII$Assay[all_chipsII$Assay=="tnpA_205"]<-"tnpA-07/ISEcp1B_205"
all_chipsII$Assay[all_chipsII$Assay=="tnpA_206"]<-"tnpA-06/IS1216_206"
all_chipsII$Assay[all_chipsII$Assay=="tnpA_207"]<-"tnpA-03/IS6_207" 

#After everyhting matches I'm ready to combine the ARG annotation to the results
allResults <- merge(all_chipsII, Primerset2, by="Assay")

#Now I have all Ct values in one dataframe with the annotation of 
#the different assays. Next I will make a matrix of the results 
#that can be used to make e.g. ordination plots in vegan. 

#I will combine all mean values to a matrix. 
#Rows are samples and colums assays. 
#If some assays have the same “Assay” (= gene name) 
#they will be aggregated. If this is the case, we combine the name and assay number:
#result_matrix <- dcast(final_plates, Sample~Gene.Name+Assay, value.var="mean")
#In my case the Gene.Name column has the assay name and Assay colum has the assay number
#so I dont have this problem
#################################################################
#NOTE!! If you didn't remove any assays this is what you should get:
#result_mat<- dcast(allResults, Sample~Assay, value.var="mean")
#dim(result_mat)
#[1]  24 385
#1st colum is the sample, and we have 384 assays.
#################################################################
#If you removed one assay, this is what you should have
result_matII<- dcast(allResults, Sample~Assay, value.var="meanII")
dim(result_matII)
#[1]  24 384
#1st colum is the sample, and we have 383 assays.

#I will use the "new" mean values:
result_mat <- result_matII

#Next task is to calculate relative abundances:

plot(result_mat[,2:3], main="1st 16S vs. 2nd 16S")
plot(result_mat[,2], main="1st 16S (new)")
plot(result_mat[,3], main="2nd 16S (old)")

#Seems that the first 16S is a good choice for normalization
#plus there was no unspecific amplification like in the second 16S.
#so we’ll use "the new" to calculate the deltaCt.
#outliers do not really play a role since our quantification is
#relative

#Let's calculate the relative abundance. (R = 2^-deltaCt)
#We could also calculate the relative abundances based on the old 16S.
#Just change the column number.

str(result_mat)
#I have 382 assays, assays other than 16S starting from column 4. 
#New 16S in in column 2 and the old in column 3
results_delta <- result_mat[,4:384]- (result_mat[,2])
results_relative <- 2^-(results_delta)

#Now we have all the relative abundances calculated and we can remove 
#all the assays that are notdetected in any sample can be deleted.

results_relative <- results_relative[,colSums(is.na(results_relative)) != nrow(results_relative)]

#Then check how many assays were left with dim() and after that give all NAs an artificial value 
#(NA is not the same as below detection limit, so we need a number). 
#For example deltaCt of 20, which gives relative value of 9.536743e-07 (2^-20).
#This can be anything, but probably should be smaller than the smallest real result.
#That can be checked from the result_deltaCt. So the maximum deltaCt from all results.

max(results_delta, na.rm=TRUE)
#[1] 17.43553
results_relative[is.na(results_relative)] <- 9.536743e-07
results_relative_mat<-as.matrix(results_relative)

#Now I'm almost done, just need to add some meaningful 
#row names that were lost on the way.
row.names(results_relative_mat) <- result_mat$Sample


#-------------- Data exploration ------------------# 

# I should have all the results in one matrix and I can do 
#the first ordination plots and heatmaps. I will use packages vegan and gplots, 
#that I loaded earlier.

#A basic NMDS plot is a good way to start checking that everything 
# is ok with the data

plot(metaMDS(results_relative_mat), type="text", display="sites")
#Warning message:
# In metaMDS(results_relative_mat) :
#  stress is (nearly) zero: you may have insufficient data

#Let's remember that Sample 1 was our negative control...
results_relative_mat<-results_relative_mat[-1,]

plot(metaMDS(results_relative_mat), type="text", display="sites")

#At this point importing some metadata might be a good idea. 

W_metadata <- read.table("Wafergen_metadata.txt", header=TRUE, sep="\t")

#The samples were sent for qPCR array analysis with names Sample1 - Sample 24. 
#Sample1 was the negative control and is not included in the metadata. 
#The actual sample numbers used in other analysis of this experiment 
#were according to the pens the pigs were kept.
#In order to match these results with other results, we need to change the sample names.
#Also we need shorter sample names.

levels(W_metadata$treatment)<-c(levels(W_metadata$treatment),"AB", "NTC", "M")
W_metadata$treatment[W_metadata$treatment=="carbadox"]<-"AB"
W_metadata$treatment[W_metadata$treatment=="mushroom"]<-"M"
W_metadata$treatment[W_metadata$treatment=="negative-cntl"]<-"NTC"

W_metadata$pen_treat<-interaction(W_metadata$pen, W_metadata$treatment, sep="_")

#check that the dimensions are ok (and also that the sample info is ok if you are unsure)
dim(results_relative_mat)
dim(W_metadata)

rownames(results_relative_mat)<-W_metadata$pen_treat
plot(metaMDS(results_relative_mat), type="text", display="sites")

#And heatmap is also a good tool for understanding data. 
#Some transformation is probably needed, here we used square root (sqrt()).
#You can also do it without transormation. Depends on your data.

#with transformation
heatmap.2(as.matrix(sqrt(results_relative_mat)), col = rev(heat.colors(100)), 
          trace = "none", density.info = "none",  srtCol=45, margins=c(5,8))
#without
heatmap.2(as.matrix(results_relative_mat), col = rev(heat.colors(100)), 
          trace = "none", density.info = "none",  srtCol=45, margins=c(5,8))

#Heatmap from only the most abundant genes is usually more helpful. 
#For example the ones where the maximum abundance is over 0.01. 
result_abund001_mat <- results_relative_mat[,apply(results_relative_mat, 2, max)>0.01]
result_abund0001_mat <- results_relative_mat[,apply(results_relative_mat, 2, max)>0.001]

heatmap.2(as.matrix(sqrt(result_abund001_mat)), col = rev(heat.colors(100)), 
          trace = "none", density.info = "none",  srtCol=45, margins=c(5,8))

heatmap.2(as.matrix(sqrt(result_abund0001_mat)), col = rev(heat.colors(100)), Colv=FALSE, 
          trace = "none", density.info = "none",  srtCol=45, margins=c(5,8))

#Lets put the samples in order
results_relative_matO<-results_relative_mat[order(W_metadata$treatment),]

heatmap.2(as.matrix(sqrt(results_relative_matO)), col = rev(heat.colors(100)),
          dendrogram="none",Rowv=FALSE, Colv=FALSE, 
          trace = "none", density.info = "none",  srtCol=45, margins=c(5,8))

result_abund001_mat <- results_relative_matO[,apply(results_relative_matO, 2, max)>0.01]
result_abund0001_mat <- results_relative_matO[,apply(results_relative_matO, 2, max)>0.001]

heatmap.2(as.matrix(sqrt(result_abund0001_mat)), col = rev(heat.colors(100)),
          dendrogram="none",Rowv=FALSE, Colv=FALSE, 
          trace = "none", density.info = "none",  srtCol=45, margins=c(5,8))


# We can also make a heatmap with ggplot2. 
# The data needs to be in a long format.
result_abund001_mat.m<-melt(result_abund001_mat)

#With ggplot2 the colors work better with NA's. 
result_abund001_mat.m$valueII<-result_abund001_mat.m$value
min(result_abund001_mat.m$valueII)
result_abund001_mat.m$valueII[result_abund001_mat.m$valueII==9.536743e-07]<-NA

#lets check what genes are among the most abundant ones:

levels(result_abund001_mat.m$Var2)

#[1] "aac(6)-im_417"      "aadE_174"           "ant6-ia_431"        "ant6-ib_1541"       "aph(2)-Id_104"     
#[6] "aph3-III_1540"      "aphA3_14"           "erm(B)_804"         "erm(F)_23"          "erm(F)_817"        
#[11] "erm(O)_808"         "erm(Q)_809"         "ermX_209"           "IS613_26"           "lnuA_251"          
#[16] "lnuC_1519"          "sat4_49"            "tet(32)_54"         "tet(36)_22"         "tet44_1539"        
#[21] "tetL_195"           "tetO_192"           "tetQ_185"           "tetW_191"           "tetX_196"          
#[26] "Tp614_25"           "tnpA-06/IS1216_206"

#Then lets order the ARG assays we have according to the antibiotic they give resistance

ARGord<-c("aac(6)-im_417","aadE_174","ant6-ia_431","ant6-ib_1541",  "aph(2)-Id_104", "aph3-III_1540",
          "aphA3_14", "sat4_49","erm(B)_804","erm(F)_23","erm(F)_817","erm(O)_808", "erm(Q)_809","ermX_209", 
          "lnuA_251" , "lnuC_1519","tet(32)_54",   
          "tet(36)_22",    "tet44_1539", "tetL_195", "tetO_192", "tetQ_185","tetW_191","tetX_196", 
          "IS613_26", "tnpA-06/IS1216_206","Tp614_25" )  

result_abund001_mat.m$Genes <- factor(result_abund001_mat.m$Var2, levels = ARGord)

#Lets also order the samples
result_abund001_mat.mO<-result_abund001_mat.m
levels(result_abund001_mat.m$Var1)

result_abund001_mat.mO$Var1<-factor(result_abund001_mat.mO$Var1, 
                                    levels = c("7_NTC", "9_NTC","11_NTC","19_NTC","29_NTC","31_NTC",
                                               "10_AB","13_AB", "16_AB","21_AB","23_AB","26_AB",  
                                               "4_M","15_M", "17_M","22_M","30_M",
                                               "3_ZnCu","5_ZnCu","14_ZnCu","18_ZnCu","28_ZnCu","32_ZnCu"))
#and give colors to our treatments
HMxcol<-c(rep("#924e7d", 6), rep("#749b9d", 6), rep("#8b58c8",5), rep("#aeb646",6))
#and get rid of the assay number in the gene name...
result_abund001_mat.mO$Gene_name<-result_abund001_mat.mO$Genes
result_abund001_mat.mO<- separate(data = result_abund001_mat.mO, col = Gene_name, into = c("Gene_name", "AssayNO"), sep = "_")

#For the color breaks:  
min(subset(result_abund001_mat.m$valueII, result_abund001_mat.m$valueII!="NA"))
#[1] 1.220662e-05
max(subset(result_abund001_mat.m$valueII, result_abund001_mat.m$valueII!="NA"))
#[1] 0.3530962
#https://www.rdocumentation.org/packages/colorspace/versions/1.4-1/topics/scale_colour_continuous_sequential

HM001 <- ggplot(result_abund001_mat.mO, aes(x=Var1, y=Genes, fill=valueII))+
  geom_tile()+
  scale_fill_continuous_sequential(palette =  "Lajolla",
    #I don't want the most darkest color and I want to reverse the colors
   begin = 1, end = 0.1, na.value="white", name="Relative abundance")+
  theme(panel.background=element_rect(fill="black", colour="black")) +
  theme(axis.text.y=element_text(colour= "black",size=10))+
  theme(axis.text.x = element_text(colour= HMxcol,angle = 90, hjust=0, vjust=0, size=10))+
  theme(panel.border=element_blank())+
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(legend.position="bottom")+
  theme(legend.key.size=unit(0.2, "cm"))+
  theme(legend.key.width=unit(0.5, "cm"))
HM001

#How many positive assays / group we have?

dim(results_relative_mat[grep("NTC",rownames(results_relative_mat)),]
    [, which(colSums(results_relative_mat[grep("NTC",rownames(results_relative_mat)),]) != 6*(9.536743e-07))])
#6 110

dim(results_relative_mat[grep("AB",rownames(results_relative_mat)),]
    [, which(colSums(results_relative_mat[grep("AB",rownames(results_relative_mat)),]) != 6*(9.536743e-07))])
# 6 106

dim(results_relative_mat[grep("M",rownames(results_relative_mat)),]
    [, which(colSums(results_relative_mat[grep("M",rownames(results_relative_mat)),]) != 5*(9.536743e-07))])
# 5 100

dim(results_relative_mat[grep("ZnCu",rownames(results_relative_mat)),]
    [, which(colSums(results_relative_mat[grep("ZnCu",rownames(results_relative_mat)),]) !=6*(9.536743e-07))])
# 6 103

# What genes are detected in only one treatment group and how many genes are shared?
# venn diagrams

NTC_genes<-as.matrix(results_relative_mat[grep("NTC",rownames(results_relative_mat)),]
            [, which(colSums(results_relative_mat[grep("NTC",rownames(results_relative_mat)),]) != 6*(9.536743e-07))])
str(NTC_genes)
NTC_genesU<-colnames(NTC_genes)
str(NTC_genesU)
length(NTC_genesU)

AB_genes<-as.matrix(results_relative_mat[grep("AB",rownames(results_relative_mat)),]
                     [, which(colSums(results_relative_mat[grep("AB",rownames(results_relative_mat)),]) != 6*(9.536743e-07))])
str(AB_genes)
AB_genesU<-colnames(AB_genes)
str(AB_genesU)
length(AB_genesU)

M_genes<-as.matrix(results_relative_mat[grep("M",rownames(results_relative_mat)),]
                    [, which(colSums(results_relative_mat[grep("M",rownames(results_relative_mat)),]) != 5*(9.536743e-07))])
str(M_genes)
M_genesU<-colnames(M_genes)
str(M_genesU)
length(M_genesU)

ZnCu_genes<-as.matrix(results_relative_mat[grep("ZnCu",rownames(results_relative_mat)),]
                    [, which(colSums(results_relative_mat[grep("ZnCu",rownames(results_relative_mat)),]) != 6*(9.536743e-07))])
str(ZnCu_genes)
ZnCu_genesU<-colnames(ZnCu_genes)
str(ZnCu_genesU)
length(ZnCu_genesU)

input<-list(NTC=NTC_genesU, AB=AB_genesU, M=M_genesU, ZnCu=ZnCu_genesU)

venn(input)

universe <- unique(c(NTC_genesU, AB_genesU, M_genesU, ZnCu_genesU))
NTC.l <- universe %in% NTC_genesU
AB.l <- universe %in% AB_genesU
M.l <- universe %in% M_genesU
ZnCu.l <- universe %in% ZnCu_genesU

universe[NTC.l & !AB.l & !M.l & !ZnCu.l]
#[1] "aadA2_97"         "dfrA12_59"        "dfra21_603"       "erm(E)_806"       "int1-a-marko_336" "lncF_FIC_1564"   
#[7] "mtrD_253"         "oprD_234"         "spcN_406"         "sugE_1549"        "VanB_211"         "vanHB_215"       
#[13] "vanSB_1521"       "vanYB_317

universe[!NTC.l & AB.l & !M.l & !ZnCu.l]
#[1] "aph(3)-ia_435"  "mef(B)_820"     "orf37-IS26_365" "vanTC_315"  

universe[!NTC.l & !AB.l & M.l & !ZnCu.l]
#[1] "aac(3)-Xa_1545" "Aac6-Aph2_412"  "cmlA1_127"      "cmlA5_375"      "IS1111_376"     "QnrB4_1202" 

universe[!NTC.l & !AB.l & !M.l & ZnCu.l]
#[1] "aac(6)-Iy_1502"      "floR_913"            "IncHI2-smr0018_1566" "IncN_rep_340"   

###############################
# What types of resistance genes and MGEs we have?
#Pie_charts will tell. First some data arranging:

results_relative_mat.m<-melt(results_relative_mat)
#results_relative_mat.mII<-results_relative_mat.m[,-c(5,6)]

View(results_relative_mat.m)
colnames(results_relative_mat.m)[1]<-"pen_treat"
colnames(results_relative_mat.m)[2]<-"Assay"
results_relative_mat.mII <- results_relative_mat.m

# ARG annotation data:
Assay_annot<-Primerset2[-c(1,2),]

#And those will be merged:
results_relative_mat.mII_annot<-merge(results_relative_mat.mII, Assay_annot,by = "Assay")

# ARGs and MGEs separately:
ARGresults_relative_mat.mII_annot<-subset(results_relative_mat.mII_annot, Mechanism!="MGE" )
MGEresults_relative_mat.mII_annot<-subset(results_relative_mat.mII_annot, Mechanism=="MGE" )

# unused levels need to be removed from both data frames
ARGresults_relative_mat.mII_annot$Mechanism<-droplevels(ARGresults_relative_mat.mII_annot$Mechanism)
MGEresults_relative_mat.mII_annot$Mechanism<-droplevels(MGEresults_relative_mat.mII_annot$Mechanism)
ARGresults_relative_mat.mII_annot$Classification<-droplevels(ARGresults_relative_mat.mII_annot$Classification)
MGEresults_relative_mat.mII_annot$Classification<-droplevels(MGEresults_relative_mat.mII_annot$Classification)
#checking:
levels(ARGresults_relative_mat.mII_annot$Classification)

# I will calulate how many positive assays we have per antibiotic group:
ARGresults_relative_mat.mII_annot$n_class<-rep(1, nrow(ARGresults_relative_mat.mII_annot)) 

ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_class[Classification=="Aminoglycoside"]<-length(unique(subset
  (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Classification=="Aminoglycoside"))))

ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_class[Classification=="Fluoroquinolone"]<-length(unique(subset
 (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Classification=="Fluoroquinolone"))))

ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_class[Classification== "MDR"]<-length(unique(subset
  (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Classification== "MDR"))))

ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_class[Classification==  "MLSB"]<-length(unique(subset
     (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Classification==  "MLSB"))))

ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_class[Classification=="Phenicol"]<-length(unique(subset
   (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Classification== "Phenicol" ))))

ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_class[Classification=="Sulfonamide"]<-length(unique(subset
     (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Classification== "Sulfonamide"))))

ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_class[Classification== "Tetracycline"]<-length(unique(subset
   (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Classification==  "Tetracycline" ))))

ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_class[Classification== "Vancomycin" ]<-length(unique(subset
    (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Classification==  "Vancomycin"  ))))

ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_class[Classification== "Beta-lactam" ]<-length(unique(subset
   (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Classification==  "Beta-lactam"  ))))

ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_class[Classification== "Trimethoprim"]<-length(unique(subset
    (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Classification==  "Trimethoprim" ))))

ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_class[Classification== "Other" ]<-length(unique(subset
   (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Classification==  "Other" ))))

# Organizing levels
ARGresults_relative_mat.mII_annot$Classification<-factor(ARGresults_relative_mat.mII_annot$Classification, 
        levels=c("Aminoglycoside","Beta-lactam","Fluoroquinolone", "MDR", "MLSB", "Other", "Phenicol","Sulfonamide",
         "Tetracycline","Trimethoprim", "Vancomycin" ))

# We need the unique observations for Pie chart:
Class_pie<-ARGresults_relative_mat.mII_annot
Class_pie<-Class_pie[,-c(1:4)]
Class_pie_unig<-unique(Class_pie)

Class_pie_unigII <- Class_pie_unig %>%
  arrange(desc(Classification)) %>%
  mutate(lab.ypos = cumsum(n_class) - 0.5*n_class)
# the plot:
ggplot(Class_pie_unigII, aes(x="", y=n_class, fill=Classification))+
  geom_bar(width = 1, stat = "identity")+
  theme_bw()+
  scale_fill_manual(values = topo.colors(11))+
  coord_polar("y", start=0)+
  geom_text(aes(y = lab.ypos, label = n_class), color = "Black")

# Same for resistance mechanisms:
ARGresults_relative_mat.mII_annot$Mechanism<-droplevels(ARGresults_relative_mat.mII_annot$Mechanism)
levels(ARGresults_relative_mat.mII_annot$Mechanism)

ARGresults_relative_mat.mII_annot$n_mech<-rep(1, nrow(ARGresults_relative_mat.mII_annot)) 

ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_mech[Mechanism=="deactivate"]<-length(unique(subset
  (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Mechanism=="deactivate"))))
ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_mech[Mechanism=="efflux"]<-length(unique(subset
  (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Mechanism=="efflux"))))
ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_mech[Mechanism=="other" ]<-length(unique(subset
   (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Mechanism=="other"))))
ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_mech[Mechanism=="protection" ]<-length(unique(subset
     (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Mechanism=="protection"))))
ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_mech[Mechanism== "regulator" ]<-length(unique(subset
       (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Mechanism== "regulator" ))))
ARGresults_relative_mat.mII_annot<- within(ARGresults_relative_mat.mII_annot, n_mech[Mechanism== "unknown"]<-length(unique(subset
      (ARGresults_relative_mat.mII_annot$Assay, ARGresults_relative_mat.mII_annot$Mechanism== "unknown"  ))))
# all observations:

#Again, only the unique observations:
Mech_pie<-ARGresults_relative_mat.mII_annot
Mech_pie<-Mech_pie[,-c(1:3,5,6)]
Mech_pie_unig<-unique(Mech_pie)

Mech_pie_unigII <- Mech_pie_unig %>%
  arrange(desc(Mechanism)) %>%
  mutate(lab.ypos = cumsum(n_mech) - 0.5*n_mech)

ggplot(Mech_pie_unigII, aes(x="", y=n_mech, fill=Mechanism))+
  geom_bar(width = 1, stat = "identity")+
  theme_bw()+
  scale_fill_manual(values = terrain.colors(7))+
  coord_polar("y", start=0)+
  geom_text(aes(y = lab.ypos, label = n_mech), color = "Black")

# For MGEs:
levels(MGEresults_relative_mat.mII_annot$Classification)

MGEresults_relative_mat.mII_annot$n_class<-rep(1, nrow(MGEresults_relative_mat.mII_annot)) 

MGEresults_relative_mat.mII_annot<- within(MGEresults_relative_mat.mII_annot, n_class[Classification=="Insertional"]<-length(unique(subset
    (MGEresults_relative_mat.mII_annot$Assay, MGEresults_relative_mat.mII_annot$Classification=="Insertional"))))
MGEresults_relative_mat.mII_annot<- within(MGEresults_relative_mat.mII_annot, n_class[Classification=="Integrase" ]<-length(unique(subset
    (MGEresults_relative_mat.mII_annot$Assay, MGEresults_relative_mat.mII_annot$Classification=="Integrase" ))))
MGEresults_relative_mat.mII_annot<- within(MGEresults_relative_mat.mII_annot, n_class[Classification=="Transposase"]<-length(unique(subset
     (MGEresults_relative_mat.mII_annot$Assay, MGEresults_relative_mat.mII_annot$Classification=="Transposase"))))
MGEresults_relative_mat.mII_annot<- within(MGEresults_relative_mat.mII_annot, n_class[Classification=="Plasmid"]<-length(unique(subset
   (MGEresults_relative_mat.mII_annot$Assay, MGEresults_relative_mat.mII_annot$Classification=="Plasmid"))))

#Only unique observations:
MGEClass_pie<-MGEresults_relative_mat.mII_annot
MGEClass_pie<-MGEClass_pie[,-c(1:4)]
MGEClass_pie_unig<-unique(MGEClass_pie)

MGEClass_pie_unigII <- MGEClass_pie_unig %>%
  arrange(desc(Classification)) %>%
  mutate(lab.ypos = cumsum(n_class) - 0.5*n_class)

ggplot(MGEClass_pie_unigII, aes(x="", y=n_class, fill=Classification))+
  geom_bar(width = 1, stat = "identity")+
  theme_bw()+
  scale_fill_brewer(palette="Reds")+
  coord_polar("y", start=0)+
  geom_text(aes(y = lab.ypos, label = n_class), color = "Black")


#################################################################################

#Now I know something about the data. Before I start to do statistical analyses,
# I need to check that certain assumptions are met.
#I try to follow the protocol by Zuur et al when ever possible.
#https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/j.2041-210X.2009.00001.

# I want to start from ***INDEPENDENCE***, since it is the most important assumption:

# We had only one sampling point and only one sample per pen.
# Thus we dont really have "traditional" independence problem.
# If we would have taken samples before and after some treatment from the same 
# animals, we would have to take this into account in our statistical analyses
# (use mixed models, etc.)

# --- Outliers in Y & X --- #

# We have only one X (=treatment) and several Y's. 
# Let's see what is going on. 

result_abund001_mat.mO$pen_treat<-result_abund001_mat.mO$Var1
result_abund001_mat.mO<- separate(data = result_abund001_mat.mO, col =pen_treat, into = c("pen", "treat"), sep = "_")

ggplot(result_abund001_mat.mO, aes(sqrt(value),factor(treat))) + 
  geom_point() +
  geom_text(aes(label=result_abund001_mat.mO$Gene_name),angle=45,size=2,hjust=1,vjust=1 )

#These were the most abundant genes.

#for all the genes:
results_relative_mat.m<-melt(results_relative_mat)
results_relative_mat.m$Gene_name<-results_relative_mat.m$Var2
results_relative_mat.m<- separate(data = results_relative_mat.m, 
                                  col = Gene_name, into = c("Gene_name", "AssayNO"), sep = "_")
#Warning message:
# Expected 2 pieces. Additional pieces discarded in 184 rows
#[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...]. 

#Those gene names.... Does not matter at this point.

results_relative_mat.m$pen_treat<-results_relative_mat.m$Var1
results_relative_mat.m<- separate(data =results_relative_mat.m, col =pen_treat, into = c("pen", "treat"), sep = "_")

ggplot(results_relative_mat.m, aes(sqrt(value),factor(treat))) + 
  geom_point() +
  geom_text(aes(label=results_relative_mat.m$Gene_name),angle=45,size=2,hjust=1,vjust=1 )

#I wouldnt say we have outliers... Same genes are abundant in all the treatments and the 
#differences are small. Also since we have normalized our results against 16S,
#we have taken into account the possibility that some samples have more bacterial 
#DNA than others.

#However, things are probably very different if we look at the individual genes.
#I will come back to this issue later, as it is not a good idea to plot 
#all the genes...

# --- Homogenity --- #

#Do I have the same number of observations in treatment groups?

boxplot(result_abund001_mat.mO$value ~ factor(result_abund001_mat.mO$treat), 
        varwidth = TRUE, log ="y", xlab = "Treatment" , ylab = "Relative abundance")

boxplot(results_relative_mat.m$value ~ factor(results_relative_mat.m$treat), 
        varwidth = TRUE, log ="y", xlab = "Treatment" , ylab = "Relative abundance")

# Hmmm. NTC group has the highest median.
#The command “varwidth=TRUE” means that the width of the box is relative to the amount of observations. 
#So looks quite good! And I know that we have 6 samples / group, except Mushroom-group has 5.
#However, the value we used for replacing NA's (zeros) could cause trouble. 


# --- Normality --- #

#Often people say that the data should be normally distributed, 
#if you want to apply linear models. According to Zuur et al. 
#this is not exactly true: The normality should be met at each X value. 
#This is of course difficult to check if the data does not have multiple 
#observations for each X and we don't. 
#Moreover, the Y’s (response variables = what is usually meant by “data”) 
#contain also the effects of all the explanatory variables (other X’s).
#Therefore it is misleading to assess normality according the Y data (Zuur et al).
#Better choice is to apply a model and plot the pooled residuals. 
#Since the linear models should be LINEAR, the residuals should have 
#the same distance to fitted values (predicted means). If they do not, 
#then the model is not linear and the plot of residuals against fitted 
#values will show a curve or some other pattern. This means that 
#the RESIDUALS should be normally distributed, if linear models are used.

#We know this about normal distribution:
# "It is a continuous probability distribution in which the random variable can take 
#any value". 
#Our values are "relative abundances", so something between 0 and 1.
#So we don't exactly meet the normality asumption and we cannot use a model with that asumption.


# ---  "Fixed X" --- #

#This assumption means that you know the values of X for each sample in advance. 
#You can have e.g. different fertilization rates, treatments or toxin concentrations. 
#This assumption is violated, if you happen to study some parameter of certain animals,
#which happen to be found somewhere and your X is for example the age of animals. In other words, 
#this assumption is violated if the X is somehow random. 

#With this data (qPCR array) we are safe, but with the 16S data we might be in trouble
# since for example the library size could be a random x-factor. 

#Zero trouble Y (Fitting data to a distribution)

#Are there lots of zeros in the data? We replaced the NA's (zeros) with value 9.536743e-07.
# So we dont have zeros at all, but our smallest (artificial) value could 
#cause trouble since we have definetely most of these observations.
#Usually zeros are a problem with count data, which is also not normally 
#distributed. Choices for probability distributions with count data are poisson or negative binomial.
#So what distribution we should use?

#Similar question was asked here:
#https://stats.stackexchange.com/questions/99425/distribution-for-percentage-data
# The answer:
#You are right that the binomial distribution is for discrete proportions that arise 
#from the number of 'successes' from a finite number of Bernoulli trials, 
#and that this makes the distribution inappropriate for your data. 
#You should use the Gamma distribution divided by the sum of that Gamma plus another Gamma. 
#That is, you should use the beta distribution to model continuous proportions.

#Problem is that there are not that many R-packages suitable for beta distribution.
#So would Gamma work? This can be tested after fitting a Gamma GLM, but lets see in advance:

library(fitdistrplus)
library(logspline)

?fitdist

#All values:
fit_gammaARGs<-fitdist(results_relative_mat.m$value,"gamma")
plot(fit_gammaARGs)
summary(fit_gammaARGs)

#Fitting of the distribution ' gamma ' by maximum likelihood 
#Parameters : 
#estimate  Std. Error
#shape  0.129399 0.002456243
#rate  15.312416 0.817346134
#Loglikelihood:  23842.34   AIC:  -47680.68   BIC:  -47668.6 
#Correlation matrix:
 # shape      rate
#shape 1.0000000 0.3556131
#rate  0.3556131 1.0000000

#I think that all the genes together do not match with the gamma distribution so well? 
#It might also be because of the value I used to replace the NA's.

qqp(results_relative_mat.m$value, "gamma", 
    shape = fit_gammaARGs$estimate[[1]], rate = fit_gammaARGs$estimate[[2]])
#Yeap, very poor fit

#With only detected ARGs:
qqp(subset(results_relative_mat.m$value, results_relative_mat.m$value>9.536743e-07), "gamma", 
    shape = fit_gammaARGs$estimate[[1]], rate = fit_gammaARGs$estimate[[2]])
#Pretty good!!

#Let's try with some individual genes:
fit_gamma_tetws<-fitdist(subset(results_relative_mat.m$value, 
                                results_relative_mat.m$Gene_name=="tetW"),"gamma")

plot(fit_gamma_tetws)

summary(fit_gamma_tetws)

#Fitting of the distribution ' gamma ' by maximum likelihood 
#Parameters : 
# estimate Std. Error
#shape 13.92350   4.057587
#rate  61.49422  18.247141
#Loglikelihood:  32.37661   AIC:  -60.75323   BIC:  -58.48224 
#Correlation matrix:
# shape      rate
#shape 1.0000000 0.9821073
#rate  0.9821073 1.0000000

qqp(subset(results_relative_mat.m$value, results_relative_mat.m$Gene_name=="tetW"), "gamma", 
    shape = fit_gamma_tetws$estimate[[1]], rate = fit_gamma_tetws$estimate[[2]])

#good!

fit_gamma_aadAs<-fitdist(subset(results_relative_mat.m$value, 
                                results_relative_mat.m$Gene_name=="aadA"),"gamma")

plot(fit_gamma_aadAs)

summary(fit_gamma_aadAs)

#Fitting of the distribution ' gamma ' by maximum likelihood 
#Parameters : 
# estimate  Std. Error
#shape    0.6484776   0.1617893
#rate  2610.1823697 938.7749041
#Loglikelihood:  169.6124   AIC:  -335.2248   BIC:  -332.9538 
#Correlation matrix:
# shape      rate
#shape 1.0000000 0.6939129
#rate  0.6939129 1.0000000

qqp(subset(results_relative_mat.m$value, results_relative_mat.m$Gene_name=="aadA"), "gamma", 
    shape = fit_gamma_tetws$estimate[[1]], rate = fit_gamma_tetws$estimate[[2]])

#also OK.

# So seems that when using all the data the undetected ARGs with value 9.536743e-07
#gives me trouble with the fit but individual genes might do just fine. Which will be 
#enough for what I have in mind. 

#Although I need to check this again afterwards with more genes. 

# --- Collinearity X --- #
#If X-observations correlate. We have only one X so no problems with this.

# --- Relationships Y & X --- #
#again we have only one X.

# --- Interactions --- #
#Same. Only one X.

########################################################################
#                 End of part 1 for ARGs and MGEs                      #
########################################################################

#----------------------------------------------------------------------#

###### --------------------------------------------------------- ######
#                                                                     #
#                   16S OTUs & two different methods                  #
#                                                                     #
###### --------------------------------------------------------- ######

# In this script I'm analyzing the results from 16S sequencing from the 
# same samples. I will be using two different approaches, that is, OTU-
# table that are rarefied and subsampled and and OTU table that is 
# "relative abundance normalized".
# For both, the sequences were clustered into OTUs using Mothur version xxxx
#


#--- DATA ---#

setwd("~/Library/Blaah blaa/blaa~blaa~blaa/Wafergen/ZnCu")
metadata <- "ZnCu_metadata.txt" #metadata, because the different numbering of the samples,
# I need to have different metadata file as with ARGs and MGEs
sharedfile <- "Phylo.tx.1.pick.1.subsample.shared_correctedID.txt" #Phylo table that is subsampled (rarefied)
phylotx_table <-"Phylo.tx.1.pick.shared_correctedID.txt" #not subsampled for relative abundance normalization
taxfile <- "Phylo.cons.taxonomy" #taxonomy file

#Read in OTU tables
phylOTUtable<- read.table(phylotx_table, header = TRUE)
View(phylOTUtable) #This file is used for TSS normalization

phyl_OTUsubsample<- read.table(sharedfile, header = TRUE)
View(phyl_OTUsubsample) #This is the rarefied and subsampled data frame

#Remove sequence ID and leave just sample names in both data frames
phylOTUtable <- separate(data = phylOTUtable, col = Group, into = c("sequencing_id", "Group"), sep = "_")
phylOTUtable <- phylOTUtable[,-2] #remove seq ID column

phyl_OTUsubsample <- separate(data = phyl_OTUsubsample, col = Group, into = c("sequencing_id", "Group"), sep = "_")
phyl_OTUsubsample <- phyl_OTUsubsample[,-2] #remove seq ID column

#remove the negative control from both data frames
phyl_OTUsubsampleII <- subset(phyl_OTUsubsample, Group!="WATER")
phylOTUtableII <- subset(phylOTUtable, Group!="WATER")

#I will move the sample number to the rownames of the dataframe 
rownames(phyl_OTUsubsampleII) <-phyl_OTUsubsampleII$Group
rownames(phylOTUtableII) <- phylOTUtableII$Group

# I was just being cautios, the previous data frame can be overwritten now:
phylOTUtable <- phylOTUtableII

#Read in metadata
meta <- read.table(file = metadata, sep = '\t', header = TRUE)
# I don't want to include negative controls:
View(phyl_OTUsubsampleII)
#easy way to get rid of samples that were excluded due to low sequencing quality: 
metaII<- meta[meta$sample_id %in% rownames(phyl_OTUsubsampleII),]

# I was just being cautios, the previous data frame can be overwritten now:
phyl_OTUsubsample<-phyl_OTUsubsampleII

# There is one treatmentgroup that we did not send for ARG and MGE 
# analysis and we will not include it in this analysis.
# So for clarity I will store the data with all the samples here:
ALL_samples_subsampled<-phyl_OTUsubsampleII

# I will remove that treatment group from the metadata: 
meta_noMC <- subset(meta, treatment != "Mushroom-carbadox")

#Remove the mushroom-carbadox group from data files also:
noMC_phylOTUtable <- phylOTUtable[rownames(phylOTUtable) %in% meta_noMC$sample_id,]
noMC_phyl_OTUsubsubsample <- phyl_OTUsubsample[rownames(phyl_OTUsubsample) %in% meta_noMC$sample_id,]
#Two samples were filtered out in quality control and they need to be removed from the NOT subsampled file:
noMC_phylOTUtable <- noMC_phylOTUtable[rownames(noMC_phylOTUtable) %in% rownames(noMC_phyl_OTUsubsubsample),]
#and they need to removed also from the metadata:
meta_noMC <- meta_noMC[meta_noMC$sample_id %in% rownames(noMC_phyl_OTUsubsubsample),] 

#check that these all match:
meta_noMC$sample_id
rownames(noMC_phylOTUtable)
rownames(noMC_phyl_OTUsubsubsample)

# next lines will remove the extra info that mothur includes 
#in OTU tables and outlier points
noMC_phylOTUtable <- noMC_phylOTUtable[,-c(1:3)] 
noMC_phyl_OTUsubsTable  <- noMC_phyl_OTUsubsubsample[,-c(1:3)] 

# I will need matrices instead of data frames for NMDS's.
noMC_phyl_OTUsubsTable_mat<-as.matrix(noMC_phyl_OTUsubsTable)

noMC_phylOTUtable_mat<-as.matrix(noMC_phylOTUtable)

# The singletons and doubletons were removed in mothur. 
#However, since we dont have the mushroom-carbadox group anymore, 
# it looks like we have those present in the data.
#If you would like to remove them, this is how:
##Alllargerthan3subs<- noMC_phyl_OTUsubsTable_mat[,apply(noMC_phyl_OTUsubsTable_mat, 2, max)>2]
##Alllargerthan3<- noMC_phylOTUtable_mat[,apply(noMC_phylOTUtable_mat, 2, max)>2]

####
#next I need to normalize the not subsampled data frame by dividing the OTU counts
#with the rowsum (total 16S count, which is also the library size)
rowSums(noMC_phylOTUtable_mat)

noMC_phylOTUtable_RA_mat <- apply(noMC_phylOTUtable_mat, 1, function(i) i/sum(i))
# the table just needs to be turned 
noMC_phylOTUtable_RA_mat <- t(noMC_phylOTUtable_RA_mat)
# RA stands for relative abundance

# I need to get rid of those columns (OTUs) where the column sum is zero
# these were detected only in mushroom-carbadox -group
noMC_phylOTUtable_RA_NoZero<-noMC_phylOTUtable_RA_mat[, which(colSums(noMC_phylOTUtable_RA_mat) != 0)]
noMC_phylOTUtable_RA_NoZero_mat<-noMC_phylOTUtable_RA_NoZero
min(noMC_phylOTUtable_RA_NoZero_mat)
dim(noMC_phylOTUtable_RA_NoZero_mat)

str(noMC_phyl_OTUsubsTable_mat)
str(noMC_phylOTUtable_RA_NoZero_mat)

# Now both OTU tables have been processed.

##--- Data exploration ---##

#Do the OTU tables correlate? 
#The data tables need to be identical.

setdiff(colnames(noMC_phylOTUtable_RA_NoZero_mat),colnames(noMC_phyl_OTUsubsTable_mat))
#[1] "Otu119"

setdiff(colnames(noMC_phyl_OTUsubsTable_mat), colnames(noMC_phylOTUtable_RA_NoZero_mat))
#[1] "Otu130" "Otu132" "Otu171"

OTU_RA_for_cor_mat<-noMC_phylOTUtable_RA_NoZero_mat[,-119]
OTU_Subs_for_cor_mat<-noMC_phyl_OTUsubsTable_mat[,-c(125,127,134)]

plot(OTU_RA_for_cor_mat,OTU_Subs_for_cor_mat)
#looks pretty good!
# But I want to have a nice figure for the paper
OTU_RA_for_cor_mat.m<-melt(OTU_RA_for_cor_mat)
OTU_Subs_for_cor_mat.m<-melt(OTU_Subs_for_cor_mat)
colnames(OTU_RA_for_cor_mat.m)<-c("SampleRA","OTURA","TSS_OTUs")
colnames(OTU_Subs_for_cor_mat.m)<-c("SampleSubs","OTUSubs","Rarefied_OTUs")
# the two data frames need to be in the same file:
OtusForScat<-cbind(OTU_RA_for_cor_mat.m,OTU_Subs_for_cor_mat.m)

library(ggpubr)

ggscatter(OtusForScat, x= "TSS_OTUs", y="Rarefied_OTUs", shape = 20,
          add = "reg.line",  # Add regression line
          conf.int = TRUE,   # Add confidence interval
          add.params = list(color = "chocolate3",
                            fill = "goldenrod"))+
  stat_cor(method = "spearman", label.x = 0.15, label.y = 3000)  # Add correlation coefficient


## Next step is the barplots. 

#We need the taxonomy data;
taxonomy <- read.table(taxfile, sep = "\t", header = T) #Read the taxonomy file into R
taxonomy <- separate(data = taxonomy, col = Taxonomy, 
  into = c("kingdom", "phylum", "class", "family", "order", "genus", "species"), sep = ";")
str(taxonomy)

# TSS normalized (Relative Abundance = RA) OTU table:
View(noMC_phylOTUtable_RA_NoZero_mat)
otu.summary <- prop.table(as.matrix(noMC_phylOTUtable_RA_NoZero_mat), 1) 
otu_abund <- colSums(otu.summary)
otu.summary <- rbind(otu_abund, otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]

num_genera <- 16 # 16 is the number of genera I want to have in the plots

melt_otu <- melt(otu.summary_sorted[,c(1:num_genera)])
colnames(melt_otu) <- c("Sample", "OTU", "Abundance")
tail(melt_otu)

#Same for noMC_phyl_OTUsubsubsample: 
Subs.otu.summary <- prop.table(as.matrix(noMC_phyl_OTUsubsTable ), 1) 
Subs.otu_abund <- colSums(Subs.otu.summary)
Subs.otu.summary <- rbind(Subs.otu_abund, Subs.otu.summary)
Subs.otu.summary_sorted <- Subs.otu.summary[,order(Subs.otu.summary[1,], decreasing = TRUE)]

Subs.melt_otu <- melt(Subs.otu.summary_sorted[,c(1:num_genera)])
colnames(Subs.melt_otu) <- c("Sample", "OTU", "Abundance")
tail(Subs.melt_otu)

# Metadata processing: 
meta_noMCII<- meta_noMC[-c(2,4,6)]
levels(meta_noMCII$treatment)
#[1] "carbadox"  "mushroom" "Mushroom-carbadox" "negative-cntl""ZnCu"   
# Shorter treatment names:
levels(meta_noMCII$treatment)<-c(levels(meta_noMCII$treatment),"AB", "NTC", "M")
meta_noMCII$treatment[meta_noMCII$treatment=="carbadox"]<-"AB"
meta_noMCII$treatment[meta_noMCII$treatment=="mushroom"]<-"M"
meta_noMCII$treatment[meta_noMCII$treatment=="negative-cntl"]<-"NTC"

#Colors that can be used later on:
meta_noMCII$color<-rep(1, nrow(meta_noMCII))
meta_noMCII<- within(meta_noMCII, color[treatment=="M"]<-"#8b58c8")
meta_noMCII<- within(meta_noMCII, color[treatment=="ZnCu"]<-"#aeb646")
meta_noMCII<- within(meta_noMCII, color[treatment=="NTC"]<-"#924e7d")
meta_noMCII<- within(meta_noMCII, color[treatment=="AB"]<-"#749b9d")
# I will drop the unused levels
meta_noMCII$treatment<-droplevels(meta_noMCII$treatment)
# this column will have the sample names:
meta_noMCII$pen_treat<-interaction(meta_noMCII$pen, meta_noMCII$treatment,sep = "_")
meta_noMCII$pen_treat<-droplevels(meta_noMCII$pen_treat)

# For barplot I need to merge the metadata and OTU tables

#Relative abundance OTUs (TSS normalized) first:
meta_otu_RA <- merge(meta_noMCII, melt_otu, by.x = "sample_id", by.y = "Sample")
meta_otu_RA_tax <- merge(meta_otu_RA, taxonomy)
str(meta_otu_RA_tax)
#I want to have the samples ordered
levels(meta_otu_RA_tax$treatment)
#[1] "ZnCu" "AB"   "NTC"  "M" 
unique(meta_otu_RA_tax$pen_treat)
#[1] 2_M  3_AB   19_NTC  3_ZnCu  32_ZnCu 7_NTC  4_M 14_ZnCu 17_M    
#29_NTC  10_AB   18_ZnCu 13_AB   5_ZnCu  9_NTC  
#[16] 31_NTC  15_M    22_M    28_ZnCu 30_M    21_AB   26_AB  

meta_otu_RA_tax$pen_treat<-factor(meta_otu_RA_tax$pen_treat, 
levels = c("7_NTC", "9_NTC","19_NTC","29_NTC","31_NTC",
 "10_AB","13_AB","21_AB","23_AB","26_AB",
"2_M", "4_M","15_M", "17_M","22_M","30_M",
"3_ZnCu","5_ZnCu","14_ZnCu","18_ZnCu","28_ZnCu","32_ZnCu"))

#Colors that dont look too bad next to the ARG heatmap:
Laj_16 <- c("#885f00", "#3a2100","#4a0000", "#ffc73a","#e87100", "#f2cea1",
            "#3c0031", "#cd2540", "#ff845e", "#ffa1d9","#ffde84",
            "#d11a67","#ff727e", "#863647", "#9e5d3c", "#f79534")
#Also I want to get rid of the "(100)" in the genus
#and I want to have them ordered according to abundance
meta_otu_RA_tax$Genus<- substr(meta_otu_RA_tax$genus, 1,nchar(meta_otu_RA_tax$genus)-5)
unique(meta_otu_RA_tax$Genus)
meta_otu_RA_tax$Genus<-factor(meta_otu_RA_tax$Genus, 
levels = rev(c("Prevotella","Lactobacillus","Megasphaera","Ruminococcaceae_unclassified",
   "Veillonellaceae_unclassified","Streptococcus",
     "Lachnospiraceae_unclassified","Dialister", "Blautia",
  "Clostridiales_unclassified","Roseburia","Phascolarctobacterium","Alloprevotella",
  "Faecalibacterium","Clostridium_sensu_stricto", "Bacteria_unclassified")))

#I will again give colors to our treatments:
Barxcol<-c(rep("#924e7d", 5), rep("#749b9d", 5), rep("#8b58c8",6), rep("#aeb646",6))

#The plot:
Nosubs<-ggplot(meta_otu_RA_tax, aes(x = pen_treat, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity",color="black") +
  scale_fill_manual(values = Laj_16)+
  theme_bw ()+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(colour= Barxcol,angle = 90, hjust=0, vjust=0, size=10))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  #ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(axis.text.y=element_text(colour= "black",size=10))+
  theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  #theme(axis.text.x = element_blank()) +
  ylab(paste0("Relative Abundance (top ", num_genera, " genera)"))
Nosubs

#Same for subsampled
Subs.meta_otu <- merge(meta_noMCII, Subs.melt_otu, by.x = "sample_id", by.y = "Sample")
Subs.meta_otu_tax <- merge(Subs.meta_otu, taxonomy, by="OTU")
str(Subs.meta_otu_tax)
unique(Subs.meta_otu_tax$pen_treat)
Subs.meta_otu_tax$pen_treat<-factor(Subs.meta_otu_tax$pen_treat, 
levels = c("7_NTC",
 "9_NTC","19_NTC","29_NTC","31_NTC",
 "10_AB","13_AB","21_AB","23_AB","26_AB",  
"2_M", "4_M","15_M", "17_M","22_M","30_M",
"3_ZnCu","5_ZnCu","14_ZnCu","18_ZnCu","28_ZnCu","32_ZnCu"))

# Colors for subsampled data:
Laj_16subs <- c("#3a2100", "#ffc73a","#4a0000" ,"#601078" ,
                "#cd2540", "#f2cea1",
                "#ff845e","#e31bbe", "#e87100","#ffde84", "#ffa1d9",
                "#d11a67","#ff727e", "#863647", "#9e5d3c", "#f79534")
#genuses:
Subs.meta_otu_tax$Genus<- substr(Subs.meta_otu_tax$genus, 1,nchar(Subs.meta_otu_tax$genus)-5)
unique(Subs.meta_otu_tax$Genus)
Subs.meta_otu_tax$Genus<-factor(Subs.meta_otu_tax$Genus,
 levels = rev(c( "Prevotella","Lactobacillus","Megasphaera",
 "Ruminococcaceae_unclassified", "Veillonellaceae_unclassified","Lachnospiraceae_unclassified",
  "Streptococcus", "Phascolarctobacterium","Treponema", "Dialister","Roseburia",
  "Blautia","Clostridiales_unclassified","Faecalibacterium",               
  "Alloprevotella","Clostridium_sensu_stricto"))) 

#The plot:
Subs<-ggplot(Subs.meta_otu_tax, aes(x = pen_treat, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", color= "black") +
  scale_fill_manual(values = Laj_16subs) +
  # Remove x axis title
  theme_bw()+
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(colour= Barxcol,angle = 90, hjust=0, vjust=0, size=10))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  #ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(axis.text.y=element_text(colour= "black",size=10))+
  theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  #theme(axis.text.x = element_blank()) +
  ylab(paste0("Relative Abundance subsampled (top ", num_genera, " genera)"))
Subs

# I want to have a panel plot:
#library(gridExtra)
?grid.arrange
grid.arrange(Nosubs,Subs, HM001, nrow=2,ncol=2)

##### OTU venns ######

View(noMC_phyl_OTUsubsTable_mat)
View(meta_noMCII)

#Vens for TSS normalized OTUs first:

rownames(noMC_phylOTUtable_RA_NoZero_mat)<-meta_noMCII$pen_treat
min(noMC_phylOTUtable_RA_NoZero_mat)

NTC_OTUs<-as.matrix(noMC_phylOTUtable_RA_NoZero_mat[grep("NTC",rownames(noMC_phylOTUtable_RA_NoZero_mat)),]
   [, which(colSums(noMC_phylOTUtable_RA_NoZero_mat[grep("NTC",rownames(noMC_phylOTUtable_RA_NoZero_mat)),]) != 0)])

str(NTC_OTUs)
NTC_OTUsU<-colnames(NTC_OTUs)
length(NTC_OTUsU)
#[1] 121

AB_OTUs<-as.matrix(noMC_phylOTUtable_RA_NoZero_mat[grep("AB",rownames(noMC_phylOTUtable_RA_NoZero_mat)),]
     [, which(colSums(noMC_phylOTUtable_RA_NoZero_mat[grep("AB",rownames(noMC_phylOTUtable_RA_NoZero_mat)),]) != 0)])

str(AB_OTUs)
AB_OTUsU<-colnames(AB_OTUs)
length(AB_OTUsU)
#[1] 119

M_OTUs<-as.matrix(noMC_phylOTUtable_RA_NoZero_mat[grep("M",rownames(noMC_phylOTUtable_RA_NoZero_mat)),]
  [, which(colSums(noMC_phylOTUtable_RA_NoZero_mat[grep("M",rownames(noMC_phylOTUtable_RA_NoZero_mat)),]) != 0)])

str(M_OTUs)
M_OTUsU<-colnames(M_OTUs)
length(M_OTUsU)
#[1] 123

ZnCu_OTUs<-as.matrix(noMC_phylOTUtable_RA_NoZero_mat[grep("ZnCu",rownames(noMC_phylOTUtable_RA_NoZero_mat)),]
   [, which(colSums(noMC_phylOTUtable_RA_NoZero_mat[grep("ZnCu",rownames(noMC_phylOTUtable_RA_NoZero_mat)),]) != 0)])

str(ZnCu_OTUs)
ZnCu_OTUsU<-colnames(ZnCu_OTUs)
length(ZnCu_OTUsU)
#[1] 109

OTUrel_input<-list(NTC=NTC_OTUsU, AB=AB_OTUsU, M=M_OTUsU, ZnCu=ZnCu_OTUsU)
venn(OTUrel_input)
OTUrel_universe <- unique(c(NTC_OTUsU, AB_OTUsU, M_OTUsU,ZnCu_OTUsU))
OTUrelNTC.l <- OTUrel_universe %in% NTC_OTUsU
OTUrelAB.l <- OTUrel_universe  %in% AB_OTUsU
OTUrelM.l <- OTUrel_universe  %in%  M_OTUsU
OTUrelZnCu.l <- OTUrel_universe %in% ZnCu_OTUsU

OTUrel_universe [OTUrelNTC.l & !OTUrelAB.l & !OTUrelM.l & !OTUrelZnCu.l ]
#[1] "Otu123" "Otu127"
OTUrel_universe [!OTUrelNTC.l & OTUrelAB.l & !OTUrelM.l & !OTUrelZnCu.l ]
#[1] "Otu168"
OTUrel_universe [!OTUrelNTC.l & !OTUrelAB.l & OTUrelM.l & !OTUrelZnCu.l ]
#[1] "Otu112"
OTUrel_universe [!OTUrelNTC.l & !OTUrelAB.l & !OTUrelM.l & OTUrelZnCu.l ]
#[1] "Otu117" "Otu139"

#Same for rarefied and subsampled OTUs

rownames(noMC_phyl_OTUsubsTable_mat)<-meta_noMCII$pen_treat

min(noMC_phyl_OTUsubsTable_mat)

sNTC_OTUs<-as.matrix(noMC_phyl_OTUsubsTable_mat[grep("NTC",rownames(noMC_phyl_OTUsubsTable_mat)),]
                    [, which(colSums(noMC_phyl_OTUsubsTable_mat[grep("NTC",rownames(noMC_phyl_OTUsubsTable_mat)),]) != 0)])

str(sNTC_OTUs)
sNTC_OTUsU<-colnames(sNTC_OTUs)
length(sNTC_OTUsU)
#[1] 112

sAB_OTUs<-as.matrix(noMC_phyl_OTUsubsTable_mat[grep("AB",rownames(noMC_phyl_OTUsubsTable_mat)),]
                     [, which(colSums(noMC_phyl_OTUsubsTable_mat[grep("AB",rownames(noMC_phyl_OTUsubsTable_mat)),]) != 0)])

str(sAB_OTUs)
sAB_OTUsU<-colnames(sAB_OTUs)
length(sAB_OTUsU)
#[1] 111

sM_OTUs<-as.matrix(noMC_phyl_OTUsubsTable_mat[grep("M",rownames(noMC_phyl_OTUsubsTable_mat)),]
                    [, which(colSums(noMC_phyl_OTUsubsTable_mat[grep("M",rownames(noMC_phyl_OTUsubsTable_mat)),]) != 0)])

str(sM_OTUs)
sM_OTUsU<-colnames(sM_OTUs)
length(sM_OTUsU)
#[1] 108

sZnCu_OTUs<-as.matrix(noMC_phyl_OTUsubsTable_mat[grep("ZnCu",rownames(noMC_phyl_OTUsubsTable_mat)),]
                   [, which(colSums(noMC_phyl_OTUsubsTable_mat[grep("ZnCu",rownames(noMC_phyl_OTUsubsTable_mat)),]) != 0)])

str(sZnCu_OTUs)
sZnCu_OTUsU<-colnames(sZnCu_OTUs)
length(sZnCu_OTUsU)
#[1] 98

OTUs_input<-list(NTC=sNTC_OTUsU, AB=sAB_OTUsU, M=sM_OTUsU, ZnCu=sZnCu_OTUsU)
venn(OTUs_input)
OTUs_universe <- unique(c(sNTC_OTUsU, sAB_OTUsU,sM_OTUsU, sZnCu_OTUsU))
OTUsNTC.l <- OTUs_universe %in% sNTC_OTUsU
OTUsAB.l <- OTUs_universe  %in% sAB_OTUsU
OTUsM.l <- OTUs_universe  %in%  sM_OTUsU
OTUsZnCu.l <- OTUs_universe %in% sZnCu_OTUsU

OTUs_universe [OTUsNTC.l & !OTUsAB.l & !OTUsM.l & !OTUsZnCu.l ]
#[1] "Otu098" "Otu114" "Otu123"
OTUs_universe [!OTUsNTC.l & OTUsAB.l & !OTUsM.l & !OTUsZnCu.l ]
#[1] "Otu078" "Otu104" "Otu129" "Otu168"
OTUs_universe [!OTUsNTC.l & !OTUsAB.l & OTUsM.l & !OTUsZnCu.l ]
#[[1] "Otu055" "Otu112"
OTUs_universe [!OTUsNTC.l & !OTUsAB.l & !OTUsM.l & OTUsZnCu.l ]
#[[1] "Otu109" "Otu111" "Otu117" "Otu139"


###########################################################################
#               ARGs and MGEs part 2 + stats for OTUs                     #
###########################################################################

############################################################################
#                                                                          #
#--------------  It's time for the statistical analyses! ------------------# 
#                                                                          #
############################################################################

##### ARGs and MGEs first #######

#PERMANOVA

#We use function adonis from the "vegan" package
#See tutorial:
#http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf

#How much of the variances in ARG abundance is explained by the treatment?
adonis(results_relative_mat ~ treatment, data=W_metadata, permutations=9999)
#Call:
# adonis(formula = results_relative_mat ~ treatment, data = W_metadata,      permutations = 9999) 

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#          Df SumsOfSqs  MeanSqs  F.Model R2       Pr(>F)   
#treatment  3   0.19549  0.065162  2.3094 0.26721 0.0081 **
# Residuals 19   0.53609 0.028215         0.73279          
#Total     22   0.73158                  1.00000          
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Our treatments explain ~27% of the differences in ARG & MGE abundances 
# in different groups.

#We can also do pairwise comparisons with a fuction*:

pairwise.adonis <- function(x,factors, sim.method, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method, permutations = 9999);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}
#*https://github.com/pmartinezarbizu/pairwiseAdonis, 
#original function “published” in a ResearchGate thread in 2015 
#(https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R)

pairwise.adonis(results_relative_mat, W_metadata$treatment, sim.method = "bray", p.adjust.m = "BH")
#p.adjust.m = "BH" is False discovery rate

#       pairs     F.Model  R2        p.value p.adjusted
#1     ZnCu vs M 2.574984 0.2224611  0.0401     0.0838
#2    ZnCu vs NTC 1.487533 0.1294911  0.2340     0.2583
#3    ZnCu vs AB 1.479772 0.1289026  0.2583     0.2583
#4    M vs NTC 2.213385 0.1973877  0.0419     0.0838
#5     M vs AB 4.264392 0.3214917  0.0045     0.0270
#6   NTC vs AB 2.121798 0.1750399  0.0664     0.0996

# So significat ~ 30% difference in M vs AB 

# NMDS (plots)
#lets do the NMDS again

NMDS1 <- metaMDS(results_relative_mat)
#NMDS1 <- metaMDS(results_relative_mat, trymax = 100)

#And plot it with colors
ordiplot(NMDS1,type="n", xlim=c(-0.4,0.5), ylim=c(-0.3,0.4))
with(W_metadata, points(NMDS1$points, col=W_metadata$treatment))
with(W_metadata, ordiellipse(NMDS1, treatment, kind = "se", conf = 0.95, col="goldenrod1"))
with(W_metadata, ordispider(NMDS1, treatment, col = "black", label= TRUE))

NMDS1

# Call:
# metaMDS(comm = results_relative_mat) 

#global Multidimensional Scaling using monoMDS

#Data:     results_relative_mat 
#  Distance: bray 

# Dimensions: 2 
#  Stress:     0.1579451 
#  Stress type 1, weak ties
#  Two convergent solutions found after 20 tries
#  Scaling: centring, PC rotation, halfchange scaling 
#  Species: expanded scores based on ‘results_relative_mat’ 

#Let's change the colors. Website I Want Hue is handy.
#https://medialab.github.io/iwanthue/
#These colors should be colorblind friendly
# (the colors were actually used before already)

W_metadata$color<-rep(1, nrow(W_metadata))   
W_metadata<- within(W_metadata, color[treatment=="M"]<-"#8b58c8")
W_metadata<- within(W_metadata, color[treatment=="ZnCu"]<-"#aeb646")
W_metadata<- within(W_metadata, color[treatment=="NTC"]<-"#924e7d")
W_metadata<- within(W_metadata, color[treatment=="AB"]<-"#749b9d")
# The order of colors will not work if the unused levels are not dropped:
W_metadata$treatment <- droplevels(W_metadata$treatment)

ordiplot(NMDS1,type="n", xlim=c(-0.4,0.5),
         ylim=c(-0.3,0.4),cex.axis = 1.5, cex.lab = 1.5, 
         main = "ARGs and MGEs")
with(W_metadata, ordiellipse(NMDS1, treatment, 
 kind = "se", conf = 0.95, col=c("#aeb646","#749b9d","#924e7d","#8b58c8")))
with(W_metadata, ordispider(NMDS1, treatment,col=c("#aeb646","#749b9d","#924e7d","#8b58c8"), label= TRUE, cex=.9))
with(W_metadata, points(NMDS1$points,pch=17, cex=1.7, col=W_metadata$color))

# The ARG & MGE NMDS plot supports the result from adonis.

#######  16S OTUs #############

## First I want to check that the NMDSes work

#NMDS with rarefied & subsampled OTU-table that does not have mushroom-carbadox group:
mdsII<-metaMDS(noMC_phyl_OTUsubsTable_mat)
#mdsII<-metaMDS(noMC_phyl_OTUsubsTable_mat, trymax = 100)

mdsII
#global Multidimensional Scaling using monoMDS

#Data:     wisconsin(sqrt(noMC_phyl_OTUsubsTable_mat)) 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.1550644 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘wisconsin(sqrt(noMC_phyl_OTUsubsTable_mat))’

#everything is ok

########

# I will store the NMDS data and merge it with metadata in case I need it later. 
# Since both NMDs:es have the same samples I can use the same metadata for both plots
nmdsII <-as.data.frame(mdsII$points)
nmdsII$sample <- rownames(nmdsII)
#before I merge the NMDS data with metadata, I need to make sure the samples will be in the correct order
nmdsII$sample_id <- meta_noMCII$sample_id
nmdsII$pen_treat <- meta_noMCII$pen_treat
# Check that they match 
View(nmdsII)
View(meta_noMCII)
#Ok
meta_noMCnmdsII <- merge(meta_noMCII, nmdsII, by="sample_id")
#now there are extal columns but it doesnt matter
#reason I didn't want to change the order is this:
# I want to check how much the library size (number of sequences/sample) explains the results
#I will add it to the column in the metadata:
meta_noMCnmdsII$lib_size<-rowSums(noMC_phylOTUtable_mat)

#NMDS with "relative abundance normalized" (TSS normalized) OTU table that does not have mushroom-carbadox group
mdsIII.f<-metaMDS(noMC_phylOTUtable_RA_NoZero_mat)

mdsIII.f
#Call:
 # metaMDS(comm = noMC_phylOTUtable_RA_NoZero_mat) 

#global Multidimensional Scaling using monoMDS

#Data:     noMC_phylOTUtable_RA_NoZero_mat 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.1711675 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘noMC_phylOTUtable_RA_NoZero_mat’ 

#everything is ok

# I will again store the data
mdsIII.f_df<-as.data.frame(mdsIII.f$points)
#Also to metadata
meta_noMCnmdsII$RA_mds1<-mdsIII.f_df$MDS1
meta_noMCnmdsII$RA_mds2<-mdsIII.f_df$MDS2

#Next I will do pairwise adonis comparisons for both OTU-tables

pairwise.adonis(noMC_phylOTUtable_RA_NoZero_mat, meta_noMCnmdsII$treatment, sim.method = "bray", p.adjust.m = "BH")
#        pairs   F.Model  R2         p.value p.adjusted
#1   M vs ZnCu 2.3578974 0.19080086  0.0357     0.1071
#2    M vs NTC 0.6697361 0.06926106  0.7828     0.7828
#3     M vs AB 1.7451701 0.16241438  0.0940     0.1410
#4 ZnCu vs NTC 2.1935702 0.19596699  0.0809     0.1410
#5  ZnCu vs AB 2.8060522 0.23767913  0.0317     0.1071
#6   NTC vs AB 1.4359403 0.15217777  0.1995     0.2394

pairwise.adonis(noMC_phyl_OTUsubsTable_mat, meta_noMCnmdsII$treatment, sim.method = "bray", p.adjust.m = "BH")
#      pairs   F.Model    R2       p.value p.adjusted
#1   M vs ZnCu 2.3619814 0.19106819  0.0376    0.11280
#2    M vs NTC 0.6337811 0.06578737  0.8084    0.80840
#3     M vs AB 1.7347814 0.16160379  0.0921    0.13815
#4 ZnCu vs NTC 2.2305950 0.19861771  0.0844    0.13815
#5  ZnCu vs AB 2.8470539 0.24031746  0.0271    0.11280
#6   NTC vs AB 1.4552948 0.15391321  0.1672    0.20064

#General adonis for relative abundance normalized OTUs
#that has the also the variable library size:
adonis(noMC_phylOTUtable_RA_NoZero_mat ~ meta_noMCnmdsII$treatment + meta_noMCnmdsII$lib_size, 
       permutations = 9999)
#Call:
 # adonis(formula = noMC_phylOTUtable_RA_NoZero_mat ~ meta_noMCnmdsII$treatment +
#meta_noMCnmdsII$lib_size, permutations = 9999) 
#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#                           Df SumsOfSqs  MeanSqs F.Model R2      Pr(>F)   
#meta_noMCnmdsII$treatment  3   0.20371 0.067904  1.9739 0.23861 0.0075 **
# meta_noMCnmdsII$lib_size  1   0.06523 0.065233  1.8963 0.07641 0.0725. 
#Residuals                 17   0.58481 0.034401         0.68499          
#Total                     21   0.85375                  1.00000          
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Treatment explains roughly 24 % and number of sequences would explain 8 % but P>0.05 so the result is not reliable.

adonis(noMC_phyl_OTUsubsTable_mat~ meta_noMCnmdsII$treatment + meta_noMCnmdsII$lib_size, permutations = 9999)

#Call:
 # adonis(formula = noMC_phyl_OTUsubsTable_mat ~ meta_noMCnmdsII$treatment +      meta_noMCnmdsII$lib_size, permutations = 9999) 

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#                           Df SumsOfSqs  MeanSqs F.Model R2      Pr(>F)   
#meta_noMCnmdsII$treatment  3   0.20722 0.069074  1.9717 0.23930 0.0085 **
# meta_noMCnmdsII$lib_size  1   0.06316 0.063160  1.8029 0.07294 0.0677 . 
#Residuals                 17   0.59556 0.035033         0.68776          
#Total                     21   0.86594                  1.00000          
#---
 # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Treatment explains roughly 24 % and number of sequences would explain 8 % but P>0.05 so the result is not reliable.

#NMDS plots:

#relative abundance (TSS normalized)
ordiplot(mdsIII.f,type="n", xlim=c(-0.4,0.4), ylim=c(-0.3,0.3), main = "OTUs, relative abundance")
with(meta_noMCnmdsII, ordiellipse(mdsIII.f, treatment, kind = "se", conf = 0.95, col=c("#aeb646","#749b9d","#924e7d","#8b58c8")))
with(meta_noMCnmdsII, ordispider(mdsIII.f, treatment, col =c("#aeb646","#749b9d","#924e7d","#8b58c8"), label= TRUE, cex=.7))
with(meta_noMCnmdsII, points(mdsIII.f$points,pch=19, cex=1.7, col=meta_noMCnmdsII$color))

# rarefied and subsampled
ordiplot(mdsII,type="n", xlim=c(-0.5,0.5), ylim=c(-0.4,0.4), main = "OTUs rarefied & subsampled")
with(meta_noMCnmdsII, ordiellipse(mdsII, treatment, kind = "se", conf = 0.95, col=c("#aeb646","#749b9d","#924e7d","#8b58c8")))
with(meta_noMCnmdsII, ordispider(mdsII, treatment, col=c("#aeb646","#749b9d","#924e7d","#8b58c8"), label= TRUE, cex=.7))
with(meta_noMCnmdsII, points(mdsII$points,pch=15,cex=1.7, col=meta_noMCnmdsII$color))

# I this section I studied the influence of treatments to abundances and compositions of ARGs, MGEs and OTUs
# Next I will study the changes in gene and OTU level.

##########################################
# Differentially abundant ARGs and OTUs
#########################################

# I know that the relative abundance data
# would most likely fit to gamma distribution. 
# and since we want to compare the differentially abundant ARGs and MGEs to 
# differentially abundant OTUs, it would make sense to use glm models 
# for all of them. 

# I will start from ARGs and MGEs first

View(results_relative_mat.m)
View(W_metadata)
colnames(results_relative_mat.m)[1]<-"pen_treat"
colnames(results_relative_mat.m)[2]<-"Assay"

ARGs_meta_glm<-merge(results_relative_mat.m,W_metadata, by="pen_treat")

ARGs_meta_glm<-ARGs_meta_glm[,-c(6,7,12)]
#some checking...
unique(results_relative_mat.m$Gene_name)
unique(ARGs_meta_glm$Gene_name)
unique(ARGs_meta_glm$pen_treat)

min(ARGs_meta_glm$value)
#[1] 9.536743e-07
min(ARGs_meta_glm$value[ARGs_meta_glm$value!=min(ARGs_meta_glm$value)])
#[1] 5.64135e-06

?glm

#Lets first check that the gamma model works technically

str(ARGs_meta_glm)

fit.ARGgamma1<-glm(value ~ treatment, data=ARGs_meta_glm, family=Gamma (link="log"))
#Error in glm.fit(x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  : 
#                        NA/NaN/Inf in 'x
fit.ARGgamma1<-glm(value ~ treatment, data=ARGs_meta_glm, family=Gamma (link="log"),maxit=100000)
#Error in glm.fit(x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  : 
#                        NA/NaN/Inf in 'x'

#it doesn't :D
# for figuring out why I will try some genes one by one:

tetWglm<-subset(ARGs_meta_glm, Gene_name=="tetW")

fit.tetWgamma1<-glm(value ~ treatment, data=tetWglm, family=Gamma (link="log"),maxit=10000)
summary(fit.tetWgamma1)

#Call:
 # glm(formula = value ~ treatment, family = Gamma(link = "log"), 
  #    data = tetWglm, maxit = 10000)

#Deviance Residuals: 
 # Min        1Q    Median        3Q       Max  
#-0.29099  -0.10653  -0.01989   0.09981   0.40479  

#Coefficients:
  #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  -1.56353    0.07745 -20.188 2.69e-14 ***
#  treatmentAB  -0.19714    0.10953  -1.800  0.08778 .  
#treatmentNTC  0.08247    0.10953   0.753  0.46072    
#treatmentM    0.39580    0.11487   3.446  0.00271 ** 
 # ---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for Gamma family taken to be 0.03598921)

#Null deviance: 1.67131  on 22  degrees of freedom
#Residual deviance: 0.65732  on 19  degrees of freedom
#AIC: -76.384

#Number of Fisher Scoring iterations: 4

## #This works and even some significant differences.

#I'm guessing I get the error beacuse the variance is (almost) zero in somewhere
summary(ARGs_meta_glm)
#it might also be that mushroom group has 5 samples and others have 6...

#lets check this
table(subset(ARGs_meta_glm$Gene_name,ARGs_meta_glm$value>9.536743e-07))
#yes there are several ARGs that have detected in only one sample.

VanBWglm<-subset(ARGs_meta_glm, Gene_name=="VanB")

fit.VanBgamma1<-glm(value ~ treatment, data=VanBWglm, family=Gamma (link="log"))
summary(fit.VanBgamma1)
#hmm it works...

tet32glm<-subset(ARGs_meta_glm, Gene_name=="tet(32)")

fit.tet32gamma1<-glm(value ~ treatment, data=tet32glm, family=Gamma (link="log"))
summary(fit.tet32gamma1)

IS6100glm<-subset(ARGs_meta_glm, Gene_name=="tnpA-04/IS6100")

fit.IS6100gamma1<-glm(value ~ treatment, data=IS6100glm, family=Gamma (link="log"))
summary(fit.IS6100gamma1)

acrFglm<-subset(ARGs_meta_glm, Gene_name=="acrF")

fit.acrFgamma1<-glm(value ~ treatment, data=acrFglm, family=Gamma (link="log"))
summary(fit.acrFgamma1)

strAglm<-subset(ARGs_meta_glm, Gene_name=="strA")

fit.strAgamma1<-glm(value ~ treatment, data=strAglm, family=Gamma (link="log"))
summary(fit.strAgamma1)

vanTCglm<-subset(ARGs_meta_glm, Gene_name=="vanTC")

fit.vanTCgamma1<-glm(value ~ treatment, data=vanTCglm, family=Gamma (link="log"))
summary(fit.vanTCgamma1)

#I'm not sure what is going on but the genes individually seem to work
#which is weird. But I cannot find the problem this way and I will find out the reason for
#this error pretty soon :)

#so to the next step:

#library(multcomp)

?glht
#lets test it:
glht.acrFgamma1 <- glht(fit.acrFgamma1, mcp(treatment = "Tukey"))

summary(glht(glht.acrFgamma1))

#Simultaneous Tests for General Linear Hypotheses

#Linear Hypotheses:
# Estimate Std. Error z value Pr(>|z|)
#AB - ZnCu == 0    0.1801     0.9883   0.182    0.998
#NTC - ZnCu == 0   0.3344     0.9883   0.338    0.987
#M - ZnCu == 0     1.5278     1.0365   1.474    0.453
#NTC - AB == 0     0.1543     0.9883   0.156    0.999
#M - AB == 0       1.3477     1.0365   1.300    0.563
#M - NTC == 0      1.1935     1.0365   1.151    0.657
#(Adjusted p values reported -- single-step method)

#Warning messages:
#1: In chkdots(...) : Argument(s) ‘complete’ passed to ‘...’ are ignored
#2: In chkdots(...) : Argument(s) ‘complete’ passed to ‘...’ are ignored

#I need to store the results
acrF_gamma<-summary(glht(glht.acrFgamma1))
#And test if I can print them 
acrF_gamma$test[-(1:2)]
#I would like to have them in a data frame
acrF_gammadf<-as.data.frame(acrF_gamma$test[-(1:2)])
#yay!


#since this works, lets make a fuction 
# I'm feeling great so I will name the function accordingly :)
JohannaIsTheBest<-function(ARGx) {
  TEMP<-subset(ARGs_meta_glm,Gene_name==ARGx)
  gammafit <-glm(value~ treatment, data=TEMP, family=Gamma (link="log"),maxit = 10000)
  glht.gamma <- glht(gammafit, mcp(treatment = "Tukey"), p.adjust="BH")
  gamma_glht_df<-as.data.frame(summary(glht(glht.gamma))$test[-(1:2)])
  return(gamma_glht_df) 
}

JohannaIsTheBest("tetW")

#             coefficients     sigma      tstat      pvalues        type
#AB - ZnCu   -0.19713759 0.1095281 -1.7998815 2.732429e-01 single-step
#NTC - ZnCu   0.08246763 0.1095281  0.7529359 8.753653e-01 single-step
#M - ZnCu     0.39580496 0.1148740  3.4455565 3.043499e-03 single-step
#NTC - AB     0.27960522 0.1095281  2.5528175 5.241734e-02 single-step
#M - AB       0.59294255 0.1148740  5.1616762 1.409397e-06 single-step
#M - NTC      0.31333733 0.1148740  2.7276602 3.191257e-02 single-step

#Warning messages:
 # 1: In chkdots(...) : Argument(s) ‘complete’ passed to ‘...’ are ignored
#2: In chkdots(...) : Argument(s) ‘complete’ passed to ‘...’ are ignored

# before I continue, I want to explain how to calculate fold changes, they are not the "coefficients"
# above (which are actually estimates). 

######################################################################################

######################################
# How to calculate the fold changes? #
######################################

# The output above shows "delta estimates" under the column coefficients. 
# By these I mean for example that the estimate of treatment "M" minus 
#the estimate of treatment "AB" is 0.59294255 
#we can check this:
summary(fit.tetWgamma1)
#Coefficients:
               # Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  -1.56353    0.07745 -20.188 2.69e-14 ***
# treatmentAB  -0.19714    0.10953  -1.800  0.08778 .  
#treatmentNTC  0.08247    0.10953   0.753  0.46072    
#treatmentM    0.39580    0.11487   3.446  0.00271 ** 

0.39580 -(-0.19714)
# [1] 0.59294

# Gamma distribution is a member of the exponental family and I used the log link. 
# in other words, to backtransform the "delta estimate",
# I need to do this:
exp(0.59294255)
#[1] 1.809305
#This number is the fold change. BUT WE ARE NOT DONE!

# gene tetW was way more abundant in treatment group M than in AB
# but doesn't positive fold change more than 1 usually means increasing abundance? 

# to understand what is going on, I will do some calculations with the data: 
mean(with (ARGs_meta_glm, subset(value, Gene_name=="tetW" & treatment=="M")))
#[1] 0.3110727
mean(with (ARGs_meta_glm, subset(value, Gene_name=="tetW" & treatment=="AB")))
#[1] 0.1719294
# in our glht output we had M - AB = 0.59294255. This was in log scale.
# the corresponding in normal scale would be:
0.311072/0.1719294
#[1] 1.809305
# It's the same number that we got with exp(0.59294255)!
# so gene tetW is ~ 1.8 times more abundant in treatment group M than in AB

# these link should provide further information:
#https://stats.stackexchange.com/questions/96972/how-to-interpret-parameters-in-glm-with-family-gamma/126225#126225
#https://stats.stackexchange.com/questions/161216/backtransform-coefficients-of-a-gamma-log-glmm

#################################################################################

# back to work.

# now I want to se if I can put the "summary output" in a data frame
gamma_glhtdf<-as.data.frame(JohannaIsTheBest("tetW"))

# It works so not I will test it with all the gene names.
# First a vector that has all the gene names.
gamma_tmpARG <- unique(ARGs_meta_glm$Gene_name)
# then I will put the function to work on that
gamma_tmpARG_test <- lapply(gamma_tmpARG, JohannaIsTheBest)
#Error in glm.fit(x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  : 
#NA/NaN/Inf in 'x'
#In addition: There were 50 or more warnings (use warnings() to see the first 50)

# I guess this means that all I need to go through 
# the ARGs in batches to discover what is wrong...
gamma_fewARG_test <- lapply(gamma_tmpARG[1:12], JohannaIsTheBest)
gamma_fewARG_test <- lapply(gamma_tmpARG[1:24], JohannaIsTheBest)
gamma_fewARG_test <- lapply(gamma_tmpARG[1:36], JohannaIsTheBest)
#Error in glm.fit(x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  : 
#NA/NaN/Inf in 'x'
#so somewhere here:
gamma_fewARG_test <- lapply(gamma_tmpARG[24:36], JohannaIsTheBest)
#Error in glm.fit(x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  : 
#NA/NaN/Inf in 'x'
gamma_tmpARG[24:36]
#[1] "IS613"        "int1-a-marko" "vanYB"        "IS6100"       "lnuB"         "IS1247"       "erm(O)"       "tet(36)"      "vat(E)"       "mef(B)"       "aphA3"       
#[12] "intI1F165"    "tetM" 
# It's trial and error...
JohannaIsTheBest("tet(36)")
#Error in glm.fit(x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  : 
#NA/NaN/Inf in 'x'
# Why tet36?? 
with (ARGs_meta_glm, subset(value, Gene_name=="tet(36)"))
#[1] 9.536743e-07 9.536743e-07 9.536743e-07 9.536743e-07 9.536743e-07 9.536743e-07 9.536743e-07 9.536743e-07 9.536743e-07 9.536743e-07 9.536743e-07 1.525376e-02 1.351940e-02
#[14] 1.422495e-02 9.536743e-07 9.536743e-07 9.536743e-07 9.536743e-07 9.536743e-07 9.536743e-07 9.536743e-07 9.536743e-07 9.536743e-07
#looks like it has been detected only in 3 samples and in one treatment group.

# lets test without tet(36)
ARGs_meta_glm_NOtet36<-subset(ARGs_meta_glm, Gene_name!="tet(36)")

JohannaIsTheBest<-function(ARGx) {
  TEMP<-subset(ARGs_meta_glm_NOtet36,Gene_name==ARGx)
  gammafit <-glm(value~ treatment, data=TEMP, family=Gamma (link="log"),maxit = 10000)
  glht.gamma <- glht(gammafit, mcp(treatment = "Tukey"), p.adjust="BH")
  gamma_glht_df<-as.data.frame(summary(glht(glht.gamma))$test[-(1:2)])
  return(gamma_glht_df) 
}

# I'm checking again that everything works
JohannaIsTheBest("tetW")
gamma_glhtdf<-as.data.frame(JohannaIsTheBest("tetW"))

gamma_tmpARG <- unique(ARGs_meta_glm_NOtet36$Gene_name)
gamma_tmpARG_test <- lapply(gamma_tmpARG, JohannaIsTheBest)
#There were 50 or more warnings (use warnings() to see the first 50)
#Great!! Only warnings!

# In case you would also like to control the order of our comparisons, see:
levels(ARGs_meta_glm_NOtet36$treatment)
#[1] "ZnCu" "AB"   "NTC"  "M" 

#If you want to change the order,
#make a column called for example treatment.f and with levels "NTC", "AB", "M", ZnCu".
#then rerun the function:
#JohannaIsTheBest<-function(ARGx) {
 # TEMP<-subset(ARGs_meta_glm_NOtet36,Gene_name==ARGx)
 # gammafit <-glm(value~ treatment.f, data=TEMP, family=Gamma (link="log"),maxit = 10000)
 # glht.gamma <- glht(gammafit, mcp(treatment.f = "Tukey"), p.adjust="BH")
 # gamma_glht_df<-as.data.frame(summary(glht(glht.gamma))$test[-(1:2)])
 # return(gamma_glht_df) 
#}
#gamma_tmpARG <- unique(ARGs_meta_glm_NOtet36$Gene_name)
#gamma_tmpARG_test <- lapply(gamma_tmpARG, JohannaIsTheBest)
head(gamma_tmpARG_test)
# We have a list of data frames for each ARG
names(gamma_tmpARG_test)
# they don't have names, lets give them names
names(gamma_tmpARG_test) <- gamma_tmpARG
head(gamma_tmpARG_test)
#I will use the melt fubction to get a single data frame
gamma_tmpARG_test.m<- melt(gamma_tmpARG_test)
head(gamma_tmpARG_test.m)
dim(gamma_tmpARG_test.m)
#[1] 3192   4
levels(gamma_tmpARG_test.m$variable)
#[1] "coefficients" "sigma"        "tstat"        "pvalues" 
#Now I will collect what I need to different data frames
coefDF<-subset(gamma_tmpARG_test.m, variable=="coefficients")
stderDF<-subset(gamma_tmpARG_test.m, variable=="sigma")
tstatDF<-subset(gamma_tmpARG_test.m, variable=="tstat")
pvalDF<-subset(gamma_tmpARG_test.m, variable=="pvalues")
#and remove unneccessary columns and rename one column for checking
coefDF<-coefDF[,-c(1,2)]
colnames(coefDF)<-c("coefficients","ARG_c")
stderDF<-stderDF[,-c(1,2)]
colnames(stderDF)<-c("std_error","ARG_s")
tstatDF<-tstatDF[,-c(1,2)]
colnames(tstatDF)<-c("tstat","ARG_t")
pvalDF<-pvalDF[,-c(1,2)]
colnames(pvalDF)<-c("adj_p-value","ARG_p")
# I'll collect the results in a single data frame in the order I want to havew them
Gamma_ARGs_results<-cbind(coefDF,stderDF,tstatDF,pvalDF)
#I need also the comparisons that were lost along the way
df_ARGcomp_gamma<-as.data.frame(lapply(gamma_tmpARG, JohannaIsTheBest))
#There were 50 or more warnings (use warnings() to see the first 50)
getARGcomparisons<-rownames(df_ARGcomp_gamma)
getARGcomparisons
#[1] "AB - ZnCu"  "NTC - ZnCu" "M - ZnCu"   "NTC - AB"   "M - AB"     "M - NTC" 
#Is everyhting ok?
dim(Gamma_ARGs_results)
length(unique(ARGs_meta_glm_NOtet36$Gene_name))
798 /133
#[1] 6, good, we have 6 comparisons
#I will add a column that has the comparison
Gamma_ARGs_results$comparison<-rep(getARGcomparisons,133)
#And I want to remove columns that were used for checking the data but save the gene name
Gamma_ARGs_results$Gene<-Gamma_ARGs_results$ARG_c
Gamma_ARGs_results<-Gamma_ARGs_results[,-c(2,4,6,8)]
#Renaming columns
colnames(Gamma_ARGs_results)<-c("Delta.Estimate", "Std.error", "z-value", "p.adjusted","Comparison", "Gene")
#And finally the rownames
rownames(Gamma_ARGs_results)<-interaction(Gamma_ARGs_results$Comparison,Gamma_ARGs_results$Gene, sep="_")
head(rownames(Gamma_ARGs_results))
#[1] "[1] "AB - ZnCu_ere(A)"  "NTC - ZnCu_ere(A)" "M - ZnCu_ere(A)"   "NTC - AB_ere(A)"   "M - AB_ere(A)"     "M - NTC_ere(A)"   

#And I want to filter those that are significant. 
Gamma_ARGs_results_sigP<-subset(Gamma_ARGs_results, p.adjusted <= 0.05)
# the differentially abundant genes were
unique(Gamma_ARGs_results_sigP$Gene)
#[1] "ere(A)"        "mdth"          "dfra21"        "erm(A)"        "IS1247"        "dfrA15"        "lncF"         
#[8] "tetM"          "int1-a-marko"  "ISEfm1-Entero" "erm(B)"        "IS613"         "IS200"         "tetU"         
#[15] "mcr-1"         "blaSFO"        "vanRB"         "vanHB"         "tetW"          "vat(E)"        "tetA"         
#[22] "bacA"          "cmr"           "dfrA12"        "orf37-IS26"    "aadA2"         "cmlA1"         "ermT"         
#[29] "mphA"          "IS1111"        "vanYB"         "aac(6)-Iy"     "aac3ia"        "cmlA5"         "strA"         
#[36] "QnrB4"         "aadD"          "tet(32)"       "spcN"          "tetC"          "aac(3)-Xa"     "ant6-ia"      
#[43] "Aac6-Aph2"     "aadA17"        "aph(3)-ia"     "sugE"          "VanB"          "acc3-iva"      "aph(2)-Id"    

# I can calculate the fold changes:
Gamma_ARGs_results_sigP$Fold.Change<-exp(Gamma_ARGs_results_sigP$Delta.Estimate)

#now I want to plot the results:
plot_tableGammaARG<- ARGs_meta_glm_NOtet36[ARGs_meta_glm_NOtet36$Gene_name %in% Gamma_ARGs_results_sigP$Gene,]
unique(plot_tableGammaARG$Gene_name)
str(plot_tableGammaARG)
#The data frame needs to be ordered
plot_tableGammaARG_o<-plot_tableGammaARG[order(-plot_tableGammaARG$value),]
unique(plot_tableGammaARG_o$Gene_name)
#Gene names are turned into factors so I can control the order in the plot
plot_tableGammaARG_o$Gene<-factor(plot_tableGammaARG_o$Gene_name, levels=unique(plot_tableGammaARG_o$Gene_name))

levels(plot_tableGammaARG_o$Gene)

#Now lets plot
dif_abundARGs<-ggplot(plot_tableGammaARG_o, aes(x = Gene, y = sqrt(value), fill = treatment)) + 
  geom_bar(position="dodge",stat = "identity") +
  scale_fill_manual(values=c("#924e7d","#749b9d","#8b58c8","#aeb646")) +
  theme_bw()+
  scale_x_discrete(limits = rev(levels(plot_tableGammaARG_o$Gene)))+
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend( keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  coord_flip()
dif_abundARGs
#So many genes...
#Difficult to see anything. 

# Lets see a boxplot. 
ggplot(plot_tableGammaARG_o, aes(x = Gene, y =value)) + 
  geom_boxplot (aes(color=treatment)) +
  scale_color_manual(values=c("#924e7d","#749b9d","#8b58c8","#aeb646")) +
  theme_bw()+
  scale_x_discrete(limits = levels(plot_tableGammaARG_o$Gene))+
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0.95, vjust = 1)) +
  theme(axis.text.x=element_text(colour="black")) +
  scale_y_log10(expand = c(0.001, 0.01))+
  theme(axis.text.y=element_text(colour="black")) +
  guides(color = guide_legend(title="Treatment", keywidth = .5, keyheight = .5, ncol = 1))

#I suspected this. Some (or even most) of the "significant differences" might be due to outliers. 
#I checked the outliers in the beginning, but I did it for the whole data, not for every gene one by one.
# I'm maybe not confident in saying that all these genes were differentially abbundant (P<0.05).
#removing the outliers is not a good option since we have only 6 observations per group.

#I will zoom in.
# 25 most abundant genes:
View(plot_tableGammaARG_o)
levels(plot_tableGammaARG_o$Gene)[1:25]

plot_tableGammaARGabund<- plot_tableGammaARG_o[plot_tableGammaARG_o$Gene %in% levels(plot_tableGammaARG_o$Gene)[1:25],]
plot_tableGammaARGabund$Gene<-droplevels(plot_tableGammaARGabund$Gene)

ggplot(plot_tableGammaARGabund, aes(x = Gene, y =sqrt(value), fill = treatment)) + 
  geom_bar(position="dodge",stat = "identity") +
  scale_fill_manual(values=c("#924e7d","#749b9d","#8b58c8","#aeb646")) +
  theme_bw()+
  scale_x_discrete(limits = levels(plot_tableGammaARGabund$Gene))+
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 0.95, vjust = 1)) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  guides(fill = guide_legend(reverse = T, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=8))

#and boxplot:

ggplot(plot_tableGammaARGabund, aes(x = Gene, y =value)) + 
  geom_boxplot (aes(color=treatment)) +
  scale_color_manual(values=c("#924e7d","#749b9d","#8b58c8","#aeb646")) +
  theme_bw()+
  scale_x_discrete(limits = levels(plot_tableGammaARGabund$Gene))+
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 0.95, vjust = 1)) +
  theme(axis.text.x=element_text(colour="black")) +
  scale_y_log10(expand =c(0, 0.1))+
  theme(axis.text.y=element_text(colour="black")) +
  guides(color = guide_legend(title="Treatment", keywidth = .5, keyheight = .5, ncol = 1)) 

# I'm thinking that tetM is probably the "last" gene with reliable result,
#but let's still have one more closer look.

plot_tableGammaARGabund2<- plot_tableGammaARG_o[plot_tableGammaARG_o$Gene %in% levels(plot_tableGammaARG_o$Gene)[1:17],]
plot_tableGammaARGabund2$Gene<-droplevels(plot_tableGammaARGabund2$Gene)

#Only the "detected ARGs"
ggplot(plot_tableGammaARGabund2, aes(x = Gene, y =value)) + 
  geom_boxplot (aes(color=treatment)) +
  scale_color_manual(values=c("#924e7d","#749b9d","#8b58c8","#aeb646")) +
  theme_bw()+
  scale_x_discrete(limits = levels(plot_tableGammaARGabund2$Gene))+
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 0.95, vjust = 1)) +
  theme(axis.text.x=element_text(colour="black")) +
  scale_y_log10(expand =c(0, 0.1), limits = c(1e-06, 0.36))+
  theme(axis.text.y=element_text(colour="black")) +
  guides(color = guide_legend(title="Treatment", keywidth = .5, keyheight = .5, ncol = 1)) 
#Warning message:
#Removed 151 rows containing non-finite values (stat_boxplot).
# I got the warning because I filtered out the genes that were not actually detected. 
# In my opinion these should not be show3n in the plot. 

#With 12 genes:
plot_tableGammaARGabund3<- plot_tableGammaARG_o[plot_tableGammaARG_o$Gene %in% levels(plot_tableGammaARG_o$Gene)[1:12],]
plot_tableGammaARGabund3$Gene<-droplevels(plot_tableGammaARGabund3$Gene)

ggplot(plot_tableGammaARGabund3, aes(x = Gene, y =value)) + 
  geom_boxplot (aes(color=treatment)) +
  scale_color_manual(values=c("#924e7d","#749b9d","#8b58c8","#aeb646")) +
  theme_bw()+
  scale_x_discrete(limits = rev(levels(plot_tableGammaARGabund3$Gene)))+
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 0.95, vjust = 1)) +
  theme(axis.text.x=element_text(colour="black")) +
  scale_y_log10(expand =c(0, 0.2), limits = c(0.98e-06, 0.36))+
  theme(axis.text.y=element_text(colour="black")) +
  guides(color = guide_legend(title="Treatment", keywidth = .5, keyheight = .5, ncol = 1))
# Warning message:
#  Removed 58 rows containing non-finite values (stat_boxplot). 

#With 10 genes:
plot_tableGammaARGabund4<- plot_tableGammaARG_o[plot_tableGammaARG_o$Gene %in% levels(plot_tableGammaARG_o$Gene)[1:10],]
plot_tableGammaARGabund4$Gene<-droplevels(plot_tableGammaARGabund4$Gene)

# I want to have the treatment groups in certain order and match the color vector:
plot_tableGammaARGabund4$treatment.f <- factor(plot_tableGammaARGabund4$treatment, levels = c("NTC", "AB", "M", "ZnCu"))

Dif_abund_ARGs10<-ggplot(plot_tableGammaARGabund4, aes(x = Gene, y =value)) + 
  geom_boxplot (aes(color=treatment.f)) +
  scale_color_manual(values=c("#924e7d","#749b9d","#8b58c8","#aeb646")) +
  theme_bw()+
  scale_x_discrete(limits = rev(levels(plot_tableGammaARGabund4$Gene)))+
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 0.95, vjust = 1)) +
  theme(axis.text.x=element_text(colour="black")) +
  scale_y_log10(expand =c(0, 0.2), limits = c(9.536743e-07, 0.36))+
  theme(axis.text.y=element_text(colour="black")) +
  guides(color = guide_legend(title="Treatment", keywidth = .5, keyheight = .5, ncol = 1))

Dif_abund_ARGs10
# in case I would want to remove the undetected genes from the plot
ggplot(plot_tableGammaARGabund4, aes(x = Gene, y =value)) + 
  geom_boxplot (aes(color=treatment.f)) +
  scale_color_manual(values=c("#924e7d","#749b9d","#8b58c8","#aeb646")) +
  theme_bw()+
  scale_x_discrete(limits = rev(levels(plot_tableGammaARGabund4$Gene)))+
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 0.95, vjust = 1)) +
  theme(axis.text.x=element_text(colour="black")) +
  scale_y_log10(expand =c(0, 0.2), limits = c(9.536743e-06, 0.36))+
  theme(axis.text.y=element_text(colour="black")) +
  guides(color = guide_legend(title="Treatment", keywidth = .5, keyheight = .5, ncol = 1))
# Warning message:
#  Removed 20 rows containing non-finite values (stat_boxplot). 

#I will save all the results:  
write.table(Gamma_ARGs_results_sigP, file="Gamma_ARGs_results_sigP_forXLSx011420.txt",sep="\t")

# One more thing: I want to have two more heatmaps. 
# One with the genes that had significant differences according to the models
# and one with all the other (detected) genes.

#for the first one I will use this data
View(plot_tableGammaARG_o)
# and I want to order the genes according to antibiotic groups
levels(plot_tableGammaARG_o$Gene)

allSig_ARGord<-c("aac(3)-Xa","acc3-iva","aac3ia","aac(6)-Iy","Aac6-Aph2","aadA2","aadA17","aadD","aph(2)-Id","aph(3)-ia", 
                 "ant6-ia", "strA","spcN","blaSFO",  "cmlA1","cmlA5","cmr", "ere(A)","erm(A)", "erm(B)","ermT","mphA","vat(E)",    
                 "bacA","fabK" ,"mcr-1","mdth" ,"sugE", "QnrB4","tetA", "tetC","tet(32)", "tetM","tetU", "tetW",
                 "dfrA12","dfrA15","dfra21",  "VanB", "vanHB", "vanRB", "vanYB", 
                 "lncF","int1-a-marko","IS1111", "IS200","IS613", "IS1247","ISEfm1-Entero","orf37-IS26")  
  
plot_tableGammaARG_o$GeneHM <- factor(plot_tableGammaARG_o$Gene, levels = allSig_ARGord)

#And I want to also order the samples
levels(plot_tableGammaARG_o$pen_treat)

plot_tableGammaARG_o$pen_treatHM<-factor(plot_tableGammaARG_o$pen_treat, 
                                    levels = c("7_NTC", "9_NTC","11_NTC","19_NTC","29_NTC","31_NTC",
                                               "10_AB","13_AB", "16_AB","21_AB","23_AB","26_AB",  
                                               "4_M","15_M", "17_M","22_M","30_M",
                                               "3_ZnCu","5_ZnCu","14_ZnCu","18_ZnCu","28_ZnCu","32_ZnCu"))
# and for heatmap with ggplot I need to change the minumum value back to NA
plot_tableGammaARG_o$valueII<-plot_tableGammaARG_o$value
plot_tableGammaARG_o$valueII[plot_tableGammaARG_o$valueII==9.536743e-07]<-NA
# Now I can plot
 ggplot(plot_tableGammaARG_o, aes(x=pen_treatHM, y=GeneHM, fill=valueII))+
  geom_tile()+
  scale_fill_continuous_sequential(palette =  "Lajolla",
    #I don't want the most darkest color and I want to reverse the colors
begin = 1, end = 0.1, na.value="white", name="Relative abundance")+
  theme(panel.background=element_rect(fill="black", colour="black")) +
  theme(axis.text.y=element_text(colour= "black",size=10))+
  theme(axis.text.x = element_text(colour="black",angle=90,size=10))+
  theme(panel.border=element_blank())+
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(legend.position="bottom")+
  theme(legend.key.size=unit(0.2, "cm"))+
  theme(legend.key.width=unit(1, "cm"))
 
 
 #Now the genes that didnt have significant differences:
 plot_tableINSIGF_ARG<- ARGs_meta_glm_NOtet36[!(ARGs_meta_glm_NOtet36$Gene_name %in% Gamma_ARGs_results_sigP$Gene),]
 unique(plot_tableINSIGF_ARG$Gene_name)
# [1] "mdtE-yhiU"       "tetG"            "IS26"            "IS91"            "tetL"            "IncN"            "vanTC"           "ISEcp1"         
# [9] "lnuA"            "vanSB"           "tnpA-02/IS4"     "floR"            "IS6100"          "lnuB"            "erm(O)"          "mef(B)"         
# [17] "aphA3"           "tetR"            "lnuC"            "vanTG"           "erm(E)"          "tnpA-07/ISEcp1B" "mdtA"            "IncHI2-smr0018" 
# [25] "intI1F165"       "vanYD"           "cfxA"            "erm(Q)"          "tet44"           "tnpA-05/IS26"    "mefA"            "intl2"          
# [33] "tetO"            "vanWB"           "blaCMY"          "erm(F)"          "tnpA-06/IS1216"  "oprD"            "aph6ic"          "aac3-Via"       
# [41] "pbrT"            "tnpA-04/IS6100"  "mepA"            "bla-ACT"         "intl3"           "tetQ"            "sul2"            "cphA"           
# [49] "vanXB"           "blaOXY-1"        "IS1133"          "aac(3)-iid"      "tnpA-03/IS6"     "aphA1"           "tetX"            "aadA7"          
# [57] "pBS228-IncP-1"   "acrF"            "ermX"            "tetB"            "str"             "qepA"            "acrA"            "terW"           
# [65] "aac(6)-im"       "sat4"            "tetbP"           "mtrD"            "aadE"            "acrB"            "TN5"             "cmlV"           
# [73] "pcoA"            "tolC"            "Tp614"           "aadA"            "tetD"            "aacA-aphD"       "strB"            "trb-C"          
# [81] "ant6-ib"         "aph3-III"        "aph4ib"             
 
# I will first order them alphabetically:
 
INSIGF_ARGs_a<-sort(unique(plot_tableINSIGF_ARG$Gene_name))
INSIGF_ARGs_a
INSIGF_ARGs_o<-c("aac(3)-iid","aac(6)-im","aac3-Via","aacA-aphD","aadA","aadA7", "aadE","ant6-ib" ,  "aph3-III" , "aph4ib", "aph6ic",        
"aphA1","aphA3","sat4" , "str","strB","bla-ACT","blaCMY", "blaOXY-1", "oprD","cfxA", "cphA", "cmlV","floR", "qepA", 
"erm(E)","erm(F)","erm(O)" , "erm(Q)","ermX", "lnuA" , "lnuB","lnuC",                 
"acrA" ,"acrB", "acrF",  "mdtA", "mdtE-yhiU", "mef(B)", "mefA","mepA", "mtrD",  "pbrT", "pcoA","tolC" ,"sul2", "terW", "tet44" , "tetB" , "tetbP",
"tetD","tetG", "tetL", "tetO" , "tetQ","tetR", "tetX",  
"vanSB","vanTC" , "vanTG" ,"vanWB" , "vanXB", "vanYD", 
"IncHI2-smr0018","IncN" , "intI1F165" ,"intl2","intl3" ,"IS1133", "IS26",           
"IS6100" ,"IS91","ISEcp1","pBS228-IncP-1" ,"TN5","tnpA-03/IS6","tnpA-02/IS4","tnpA-04/IS6100","tnpA-05/IS26", "tnpA-06/IS1216" ,"tnpA-07/ISEcp1B","Tp614","trb-C")   

plot_tableINSIGF_ARG$GeneHM <- factor(plot_tableINSIGF_ARG$Gene_name, levels = INSIGF_ARGs_o)

#And I want to also order the samples
levels(plot_tableINSIGF_ARG$pen_treat)

plot_tableINSIGF_ARG$pen_treatHM<-factor(plot_tableINSIGF_ARG$pen_treat, 
                                         levels = c("7_NTC", "9_NTC","11_NTC","19_NTC","29_NTC","31_NTC",
                                                    "10_AB","13_AB", "16_AB","21_AB","23_AB","26_AB",  
                                                    "4_M","15_M", "17_M","22_M","30_M",
                                                    "3_ZnCu","5_ZnCu","14_ZnCu","18_ZnCu","28_ZnCu","32_ZnCu"))
# and for heatmap with ggplot I need to change the minumum value back to NA
plot_tableINSIGF_ARG$valueII<-plot_tableINSIGF_ARG$value
plot_tableINSIGF_ARG$valueII[plot_tableINSIGF_ARG$valueII==9.536743e-07]<-NA
# Plot:
ggplot(plot_tableINSIGF_ARG, aes(x=pen_treatHM, y=GeneHM, fill=valueII))+
  geom_tile()+
  scale_fill_continuous_sequential(palette =  "Light Grays", 
                    na.value="white", name="Relative abundance")+#,
  theme(panel.background=element_rect(fill="black", colour="black")) +
  theme(axis.text.y=element_text(colour= "black",size=10))+
  theme(axis.text.x = element_text(colour="black",angle=90,size=10))+
  theme(panel.border=element_blank())+
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(legend.position="bottom")+
  theme(legend.key.size=unit(0.2, "cm"))+
  theme(legend.key.width=unit(1, "cm"))


###########################
# Now the same for OTUs
###########################
# I will start with the TSS normalized OTUs (relative abundance):

View(noMC_phylOTUtable_RA_NoZero)
noMC_phylOTUtable_RA_NoZeroII<-as.data.frame(noMC_phylOTUtable_RA_NoZero)
noMC_phylOTUtable_RA_NoZeroII$SampleID<-rownames(noMC_phylOTUtable_RA_NoZeroII)
noMC_phylOTUtable_RA_NoZeroII.m<-melt(noMC_phylOTUtable_RA_NoZeroII)
colnames(noMC_phylOTUtable_RA_NoZeroII.m)<-c("SampleID", "OTU", "Abundance")

View(taxonomy)
View(meta_noMCII)
#Just in case I need to change something
meta_noMC_glm<-meta_noMCII

RA_OTUs_meta_glm<-merge(noMC_phylOTUtable_RA_NoZeroII.m, meta_noMC_glm, by.x="SampleID",by.y="sample_id")
RA_OTUs_meta_tax_glm<-merge(RA_OTUs_meta_glm, taxonomy, by="OTU")

#If I'm using gamma model, I cannot have zero's in my data. 
# Also if we think about it, 0 is not the same thing as not detected, if the data is not count data.
# Like with ARGs I need a value that is smaller than any "true" result.
min(RA_OTUs_meta_tax_glm$Abundance)
# [1] 0
min(RA_OTUs_meta_tax_glm$Abundance[RA_OTUs_meta_tax_glm$Abundance!=min(RA_OTUs_meta_tax_glm$Abundance)])
#[1] 1.490935e-05
# 100 times lower than the last true observation would be nice beacuse it is also 
# close to the value I used for ARGs. 
RA_OTUs_meta_tax_glm$Abundance[RA_OTUs_meta_tax_glm$Abundance==0]<-1.490935e-07

min(RA_OTUs_meta_tax_glm$Abundance)

#Lets first again check that the gamma model works technically...
fit.gamma1<-glm(Abundance ~ treatment, data=RA_OTUs_meta_tax_glm, family=Gamma (link="log"),maxit = 1000)
#This time it does!!!
fit.gamma1

#Call:  glm(formula = Abundance ~ treatment, family = Gamma(link = "log"), 
#          data = RA_OTUs_meta_tax_glm, maxit = 1000)

#Coefficients:
#(Intercept)  treatmentNTC    treatmentM   treatmentAB  
#-4.883e+00    -2.669e-06    -2.187e-06    -1.566e-06  

#Degrees of Freedom: 2903 Total (i.e. Null);  2900 Residual
#Null Deviance:	    28100 
#Residual Deviance: 28100 	AIC: -38700

summary(fit.gamma1)

#Call:
#glm(formula = Abundance ~ treatment, family = Gamma(link = "log"), 
#data = RA_OTUs_meta_tax_glm, maxit = 1000)

#Deviance Residuals: 
# Min       1Q   Median       3Q      Max  
#-4.4353  -4.4353  -2.4390  -0.6799   8.6391  

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  -4.883e+00  1.179e-01  -41.41   <2e-16 ***
# treatmentNTC -2.669e-06  1.749e-01    0.00        1    
#treatmentM   -2.187e-06  1.668e-01    0.00        1    
#treatmentAB  -1.566e-06  1.749e-01    0.00        1    
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for Gamma family taken to be 11.01104)

#Null deviance: 28096  on 2903  degrees of freedom
#Residual deviance: 28096  on 2900  degrees of freedom
#AIC: -38695

#Number of Fisher Scoring iterations: 239

#Now I can try check if the fit is ok
plot(fit.gamma1)
plot(residuals(fit.gamma1),fitted(fit.gamma1))
#Hit <Return> to see next plot:
#Hit <Return> to see next plot:
#Hit <Return> to see next plot:
#Hit <Return> to see next plot:

#Not really, but I didn't expect it to be with all OTUs 

#What happens with only one OTU?

gamma_OTU_test<-subset(RA_OTUs_meta_tax_glm, OTU=="Otu001")
gamma_OTU19_test<-subset(RA_OTUs_meta_tax_glm, OTU=="Otu019")
gamma_OTU19_test$treatment.f <- factor(gamma_OTU19_test$treatment, levels=c("AB","M","NTC","ZnCu"))

gamma_peekaboo<-glm(Abundance ~ treatment, data=gamma_OTU_test, family=Gamma (link="log"))
gamma_OTU19<-glm(Abundance ~ treatment.f, data=gamma_OTU19_test, family=Gamma (link="log"))

summary(gamma_peekaboo)

#Call:
# glm(formula = Abundance ~ treatment, family = Gamma(link = "log"), 
#    data = gamma_OTU_test)

#Deviance Residuals: 
# Min        1Q    Median        3Q       Max  
#-0.33614  -0.18388   0.04582   0.12071   0.32239  

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  -1.31225    0.08058 -16.286 3.23e-12 ***
# treatmentNTC -0.21913    0.11952  -1.834   0.0833 .  
#treatmentM   -0.16329    0.11395  -1.433   0.1690    
#treatmentAB  -0.22287    0.11952  -1.865   0.0786 .  
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for Gamma family taken to be 0.03895604)

#Null deviance: 0.93298  on 21  degrees of freedom
#Residual deviance: 0.74257  on 18  degrees of freedom
#AIC: -66.855

#Number of Fisher Scoring iterations: 4

# 1-(0.74257 /0.93298 )
#[1] 0.204088

#I will use again the packages for checking the fit 

#library(fitdistrplus)
#library(logspline)

fit_gamma_allOTUs<-fitdist(RA_OTUs_meta_tax_glm$Abundance, "gamma")

plot(fit_gamma_allOTUs)

summary(fit_gamma_allOTUs)

#Fitting of the distribution ' gamma ' by maximum likelihood 
#Parameters : 
# estimate  Std. Error
#shape  0.1576799 0.003127083
#rate  20.8246070 1.057179291
#Loglikelihood:  19553.84   AIC:  -39103.68   BIC:  -39091.73 
#Correlation matrix:
# shape      rate
#shape 1.0000000 0.3906531
#rate  0.3906531 1.0000000

qqp(RA_OTUs_meta_tax_glm$Abundance, "gamma", 
    shape = fit_gamma_allOTUs$estimate[[1]], rate = fit_gamma_allOTUs$estimate[[2]])

min(RA_OTUs_meta_tax_glm$Abundance)
#[1] 1.490935e-07

qqp(subset(RA_OTUs_meta_tax_glm$Abundance, RA_OTUs_meta_tax_glm$Abundance>1.490935e-07), "gamma", 
    shape = fit_gamma_allOTUs$estimate[[1]], rate = fit_gamma_allOTUs$estimate[[2]])
# Without the minimum value the fit is not perfect but OK!

#I will check again with one OTU:

fit_gammaOTU<-fitdist(gamma_OTU_test$Abundance,"gamma")

plot(fit_gammaOTU)

summary(fit_gammaOTU)

#Fitting of the distribution ' gamma ' by maximum likelihood 
#Parameters : 
# estimate Std. Error
#shape  23.74313    7.10910
#rate  101.51317   30.71759
#Loglikelihood:  35.90088   AIC:  -67.80175   BIC:  -65.61967 
#Correlation matrix:
# shape      rate
#shape 1.0000000 0.9894908
#rate  0.9894908 1.0000000

qqp(gamma_OTU_test$Abundance, "gamma", 
    shape = fit_gammaOTU$estimate[[1]], rate = fit_gammaOTU$estimate[[2]])
#Good!
#So seems the same as with ARGs, one OTU at a time is good, as well as the
#"abundant" OTUs.  

#I will move to the next step but after I'm done with all the OTUs 
# I will need to check the fit for those that are signifficant. 

glht.gammaOTU1 <- glht(gamma_peekaboo, mcp(treatment = "Tukey"))

summary(glht(glht.gammaOTU1))

#Simultaneous Tests for General Linear Hypotheses

#Linear Hypotheses:
#Estimate Std. Error z value Pr(>|z|)
#NTC - ZnCu == 0 -0.219134   0.119515  -1.834    0.257
#M - ZnCu == 0   -0.163287   0.113953  -1.433    0.478
#AB - ZnCu == 0  -0.222866   0.119515  -1.865    0.243
#M - NTC == 0     0.055847   0.119515   0.467    0.966
#AB - NTC == 0   -0.003732   0.124830  -0.030    1.000
#AB - M == 0     -0.059579   0.119515  -0.499    0.959
#(Adjusted p values reported -- single-step method)


glht.gammaOTU19 <- glht(gamma_OTU19, mcp(treatment.f = "Tukey"))
summary(glht(glht.gammaOTU19))


# I will use the function again. This time I discovered that
# Otu051 is for some reason causing problems so I will remove that

RA_OTUs_meta_tax_glmno51<-subset(RA_OTUs_meta_tax_glm, OTU!="Otu051")

RA_OTUs_meta_tax_glmno51$OTU<-droplevels(RA_OTUs_meta_tax_glmno51$OTU)

#And I could again order the treatments:
#RA_OTUs_meta_tax_glmno51$treatment2<-factor(RA_OTUs_meta_tax_glmno51$treatment, levels = c("NTC", "AB","M","ZnCu" ))

#The function for gamma OTUs:
gamma_JohannaIsTheBest<-function(Otux) {
  TEMP<-subset(RA_OTUs_meta_tax_glmno51,OTU==Otux)
  #TEMP<-subset(var_test,OTU==Otux)
  gammafit_x <-glm(Abundance ~ treatment, data=TEMP, family=Gamma (link="log"),maxit = 10000)
  # if the treatments were ordered, change the above line to below:
  #gammafit_x <-glm(Abundance ~ treatment2, data=TEMP, family=Gamma (link="log"),maxit = 10000)
  glht.gammax <- glht(gammafit_x, mcp(treatment= "Tukey"), p.adjust="BH")
  # same here:
  #glht.gammax <- glht(gammafit_x, mcp(treatment2 = "Tukey"), p.adjust="BH")
  gamma_glhtSUMdf<-as.data.frame(summary(glht(glht.gammax))$test[-(1:2)])
  return(gamma_glhtSUMdf) 
}
# I will run the function for all the OTUs
gamma_tmpOTU <- levels(RA_OTUs_meta_tax_glmno51$OTU)
#There were 50 or more warnings (use warnings() to see the first 50)
gamma_tmpOTU_test <- lapply(gamma_tmpOTU, gamma_JohannaIsTheBest)
head(gamma_tmpOTU_test)
names(gamma_tmpOTU_test)
names(gamma_tmpOTU_test) <- gamma_tmpOTU
head(gamma_tmpOTU_test)
gamma_tmpOTU_test.m<- melt(gamma_tmpOTU_test)
head(gamma_tmpOTU_test.m)

dim(gamma_tmpOTU_test.m)
#and store the results:
levels(gamma_tmpOTU_test.m$variable)
coefDF<-subset(gamma_tmpOTU_test.m, variable=="coefficients")
stderDF<-subset(gamma_tmpOTU_test.m, variable=="sigma")
tstatDF<-subset(gamma_tmpOTU_test.m, variable=="tstat")
pvalDF<-subset(gamma_tmpOTU_test.m, variable=="pvalues")
coefDF<-coefDF[,-c(1,2)]
colnames(coefDF)<-c("coefficients","OTU_c")
stderDF<-stderDF[,-c(1,2)]
colnames(stderDF)<-c("std_error","OTU_s")
tstatDF<-tstatDF[,-c(1,2)]
colnames(tstatDF)<-c("tstat","OTU_t")
pvalDF<-pvalDF[,-c(1,2)]
colnames(pvalDF)<-c("adj_p-value","OTU_p")

Gamma_OTUs_results<-cbind(coefDF,stderDF,tstatDF,pvalDF)

# I need the comparisons:
df_comp_gamma<-as.data.frame(lapply(gamma_tmpOTU, gamma_JohannaIsTheBest))
gamma_comparisons<-rownames(df_comp_gamma)
dim(Gamma_OTUs_results)
unique(RA_OTUs_meta_tax_glmno51$OTU)
786 /131
#[1] 6

Gamma_OTUs_results$comparison<-rep(gamma_comparisons,131)

Gamma_OTUs_results$OTU<-Gamma_OTUs_results$OTU_c
Gamma_OTUs_results<-Gamma_OTUs_results[,-c(2,4,6,8)]

colnames(Gamma_OTUs_results)<-c("Delta.Estimate", "Std.Error", "z-value", "p.adjusted","Comparison", "OTU")
rownames(Gamma_OTUs_results)<-interaction(Gamma_OTUs_results$Comparison,Gamma_OTUs_results$OTU, sep="_")

Gamma_OTUs_results_sigP<-subset(Gamma_OTUs_results, p.adjusted <= 0.05)
unique(Gamma_OTUs_results_sigP$OTU)
#[1] "Otu002" "Otu006" "Otu010" "Otu013" "Otu014" "Otu017" "Otu018" "Otu019" "Otu024" "Otu025" "Otu027"
#[12] "Otu029" "Otu033" "Otu034" "Otu036" "Otu040" "Otu044" "Otu046" "Otu050" "Otu052" "Otu059" "Otu062"
#[23] "Otu065" "Otu066" "Otu069" "Otu071" "Otu074" "Otu075" "Otu076" "Otu077" "Otu083" "Otu085" "Otu087"
#[34] "Otu089" "Otu090" "Otu091" "Otu093" "Otu094" "Otu095" "Otu097" "Otu098" "Otu099" "Otu100" "Otu102"
#[45] "Otu104" "Otu106" "Otu109" "Otu110" "Otu111" "Otu112" "Otu113" "Otu114" "Otu116" "Otu117" "Otu118"
#[56] "Otu119" "Otu122" "Otu123" "Otu124" "Otu127" "Otu128" "Otu129" "Otu131" "Otu134" "Otu139" "Otu141"
#[67] "Otu158" "Otu168"

# I will add the fold change column:

Gamma_OTUs_results_sigP$Fold.Change<-exp(Gamma_OTUs_results_sigP$Delta.Estimate)

mean(with (RA_OTUs_meta_tax_glm, subset(Abundance, OTU=="Otu002" & treatment=="M")))
#[1] 0.0489911
mean(with (RA_OTUs_meta_tax_glm, subset(Abundance, OTU=="Otu002" & treatment=="ZnCu")))
#[1] 0.07401708
0.0489911/0.07401708
#0.6618891
#Corresponding number in Gamma_OTUs_results_sigP$Fold.Change is 6.618891e-01.

# I will again save the results:
Gamma_OTUs_results_sigP_forXLSx <- merge(Gamma_OTUs_results_sigP, taxonomy, by= 'OTU')
rownames(Gamma_OTUs_results_sigP_forXLSx)<-rownames(Gamma_OTUs_results_sigP)
#I want to remove the (100) at the end of each genera
Gamma_OTUs_results_sigP_forXLSx$genus<- substr(Gamma_OTUs_results_sigP_forXLSx$genus, 1,nchar(Gamma_OTUs_results_sigP_forXLSx$genus)-5)
write.table(Gamma_OTUs_results_sigP_forXLSx, file="Gamma_OTUs_results_sigP_forXLSx011420.txt",sep="\t")


# Plots:
plot_tableGammaOTU<- RA_OTUs_meta_tax_glm[RA_OTUs_meta_tax_glm$OTU %in% Gamma_OTUs_results_sigP$OTU,]
unique(plot_tableGammaOTU$OTU)
plot_tableGammaOTU$OTU<-droplevels(plot_tableGammaOTU$OTU)

str(plot_tableGammaOTU)
levels(plot_tableGammaOTU$genus)
plot_tableGammaOTU$genus.f<-as.factor(plot_tableGammaOTU$genus)
levels(plot_tableGammaOTU$genus.f)

plot_tableGammaOTU_o<-plot_tableGammaOTU[order(-plot_tableGammaOTU$Abundance),]
unique(plot_tableGammaOTU_o$genus)

plot_tableGammaOTU_o$Genus<-factor(plot_tableGammaOTU_o$genus, levels=unique(plot_tableGammaOTU_o$genus))

levels(plot_tableGammaOTU_o$Genus)

str(plot_tableGammaOTU_o)

rownames(plot_tableGammaOTU_o)<-1:nrow(plot_tableGammaOTU_o)

levels(plot_tableGammaOTU_o$Genus)

AbundGenus<-c("Veillonellaceae_unclassified(100)", "Ruminococcaceae_unclassified(100)",            
              "Streptococcus(100)","Treponema(100)",                               
              "Phascolarctobacterium(100)","Roseburia(100)","Porphyromonadaceae_unclassified(100)",
              "Selenomonas(100)","Bifidobacterium(100)","Erysipelotrichaceae_unclassified(100)",        
              "Olsenella(100)","Oscillibacter(100)",   "Butyricicoccus(100)" ,"Campylobacter(100)",                           
              "Prevotellaceae_unclassified(100)",  "Ruminococcus(100)" ,                           
              "Gammaproteobacteria_unclassified(100)", "Anaerovibrio(100)" ,                           
              "Proteobacteria_unclassified(100)" , "Elusimicrobium(100)",                          
              "Fibrobacter(100)","Clostridia_unclassified(100)",                 
              "Clostridium_XI(100)",  "Pseudobutyrivibrio(100)",                      
              "Peptococcus(100)", "Pseudoscardovia(100)",                         
              "Desulfovibrio(100)",  "Faecalicoccus(100)",                           
              "Veillonella(100)" ,"Alphaproteobacteria_unclassified(100)")    

plot_tableGammaOTU_abund<- plot_tableGammaOTU_o[plot_tableGammaOTU_o$Genus %in% AbundGenus,]

dim(plot_tableGammaOTU_abund)
dim(plot_tableGammaOTU_o)

plot_tableGammaOTU_abund$Genus<-droplevels(plot_tableGammaOTU_abund$Genus)

dif_abundRA<-ggplot(plot_tableGammaOTU_abund, aes(x = Genus, y = Abundance, fill = treatment)) + 
  #ggplot(plot_tableGammaOTU_o, aes(x = Genus, y = Abundance, fill = treatment)) + 
  geom_bar(position="dodge",stat = "identity") +
  #facet_grid(treatment ~ ., scales="free")
  #scale_fill_manual(values = HUE_16subs) +
  # Remove x axis title
  scale_fill_manual(values=c("#aeb646","#924e7d","#8b58c8","#749b9d")) +
  theme_bw()+
  #scale_x_discrete(limits = rev(levels(plot_tableGammaOTU_o$Genus)))+
  scale_x_discrete(limits = rev(levels(plot_tableGammaOTU_abund$Genus)))+
  theme(axis.title.x = element_blank()) + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.12)) +
  #ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = T, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  coord_flip()

dif_abundRA


plot_tableGammaOTU_abund$Genus2<- substr(plot_tableGammaOTU_abund$genus, 1,nchar( plot_tableGammaOTU_abund$genus)-5)
plot_tableGammaOTU_abund$Genus2<-factor(plot_tableGammaOTU_abund$Genus2, levels=unique(plot_tableGammaOTU_abund$Genus2))

# I also forgot to order the treatments:
plot_tableGammaOTU_abund$treatment.f <- factor(plot_tableGammaOTU_abund$treatment, levels = c("NTC","AB","M","ZnCu"))

ggplot(plot_tableGammaOTU_abund, aes(x = Genus2, y =Abundance)) + 
  geom_boxplot (aes(color=treatment.f)) +
  #facet_grid(treatment ~ ., scales="free")
  #scale_fill_manual(values = HUE_16subs) +
  # Remove x axis title
  scale_color_manual(values=c("#924e7d","#749b9d","#8b58c8","#aeb646")) +
  theme_bw()+
  #scale_x_discrete(limits = rev(levels(plot_tableGammaOTU_o$Genus)))+
  scale_x_discrete(limits = levels(plot_tableGammaOTU_abund$Genus2))+
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 0.95, vjust = 1)) +
  theme(axis.text.x=element_text(colour="black")) +
  scale_y_log10(expand =c(0, 0.1), limits = c(1e-06, 0.36))+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 0.36)) +
  theme(axis.text.y=element_text(colour="black")) +
  #ylim(c(0,1)) +
  guides(color = guide_legend(title="Treatment", keywidth = .5, keyheight = .5, ncol = 1)) 
#Warning message:
 # Removed 100 rows containing non-finite values (stat_boxplot). 

#Maybe 14 first? 

#With 14 genera:
plot_tableGammaOTU_abund14<- plot_tableGammaOTU_abund[plot_tableGammaOTU_abund$Genus2 %in% levels(plot_tableGammaOTU_abund$Genus2)[1:14],]
plot_tableGammaOTU_abund14$Genus2<-droplevels(plot_tableGammaOTU_abund14$Genus2)

OTUs14RA<-ggplot(plot_tableGammaOTU_abund14, aes(x = Genus2, y =Abundance)) + 
  geom_boxplot (aes(color=treatment.f)) +
  #facet_grid(treatment ~ ., scales="free")
  #scale_fill_manual(values = HUE_16subs) +
  # Remove x axis title
  scale_color_manual(values=c("#924e7d","#749b9d","#8b58c8","#aeb646")) +
  theme_bw()+
  #scale_x_discrete(limits = rev(levels(plot_tableGammaOTU_o$Genus)))+
  #scale_x_discrete(limits = levels(plot_tableGammaOTU_abund14$Genus2))+
  scale_x_discrete(limits = rev(levels(plot_tableGammaOTU_abund14$Genus2)))+
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0.95, vjust = 1)) +
  theme(axis.text.x=element_text(colour="black")) +
  scale_y_log10(expand =c(0, 0.1), limits = c( 1.490935e-07, 0.36))+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 0.36)) +
  theme(axis.text.y=element_text(colour="black")) +
  #ylim(c(0,1)) +
  guides(color = guide_legend(title="Treatment", keywidth = .5, keyheight = .5, ncol = 1)) 
OTUs14RA

# Again I want to have two more heatmaps. 
# One with the Otus that had significant differences according to the models
# and one with all the other (detected) Otus.

#for the first one I will use this data
View(plot_tableGammaOTU_o)
dim(plot_tableGammaOTU_o)

#And I want to also order the samples
levels(plot_tableGammaOTU_o$pen_treat)

plot_tableGammaOTU_o$pen_treatHM<-factor(plot_tableGammaOTU_o$pen_treat, 
                                         levels = c("7_NTC", "9_NTC","19_NTC","29_NTC","31_NTC",
                                                    "10_AB","13_AB", "21_AB","23_AB","26_AB",  
                                                    "2_M","4_M","15_M", "17_M","22_M","30_M",
                                                    "3_ZnCu","5_ZnCu","14_ZnCu","18_ZnCu","28_ZnCu","32_ZnCu"))
# and for heatmap with ggplot I need to change the minumum value back to NA
plot_tableGammaOTU_o$valueII<-plot_tableGammaOTU_o$Abundance
min(plot_tableGammaOTU_o$Abundance)
#[1] 1.490935e-07
plot_tableGammaOTU_o$valueII[plot_tableGammaOTU_o$valueII==1.490935e-07]<-NA

plot_tableGammaOTU_o$genusII<- substr(plot_tableGammaOTU_o$genus, 1,nchar(plot_tableGammaOTU_o$genus)-5)

# Now I can plot
ggplot(plot_tableGammaOTU_o, aes(x=pen_treatHM, y=genusII, fill=valueII))+
  geom_tile()+
  scale_fill_continuous_sequential(palette =  "Lajolla",
 #I don't want the most darkest color and I want to reverse the colors
  begin = 1, end = 0.1, na.value="white", name="Relative abundance")+
  theme(panel.background=element_rect(fill="black", colour="black")) +
  theme(axis.text.y=element_text(colour= "black",size=10))+
  theme(axis.text.x = element_text(colour="black",angle=90,size=10))+
  theme(panel.border=element_blank())+
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(legend.position="bottom")+
  theme(legend.key.size=unit(0.2, "cm"))+
  theme(legend.key.width=unit(1, "cm"))


#Now the Otus that didnt have significant differences:
plot_tableINSIGF_Otus<- RA_OTUs_meta_tax_glmno51[!(RA_OTUs_meta_tax_glmno51$OTU %in% Gamma_OTUs_results_sigP$OTU),]
unique(plot_tableINSIGF_Otus$OTU)
#[1] Otu001 Otu003 Otu004 Otu005 Otu007 Otu008 Otu009 Otu011 Otu012 Otu015 Otu016 Otu020 Otu021 Otu022 Otu023 Otu026 Otu028 Otu030 Otu031 Otu032 Otu035
#[22] Otu037 Otu038 Otu039 Otu041 Otu042 Otu043 Otu045 Otu047 Otu048 Otu049 Otu053 Otu054 Otu055 Otu056 Otu057 Otu058 Otu060 Otu061 Otu063 Otu064 Otu067
#[43] Otu068 Otu070 Otu072 Otu073 Otu078 Otu079 Otu080 Otu081 Otu082 Otu084 Otu086 Otu088 Otu092 Otu096 Otu101 Otu103 Otu105 Otu107 Otu108 Otu115 Otu149
#131 Levels: Otu001 Otu002 Otu003 Otu004 Otu005 Otu006 Otu007 Otu008 Otu009 Otu010 Otu011 Otu012 Otu013 Otu014 Otu015 Otu016 Otu017 Otu018 Otu019 ... Otu168

#And I want to order the samples
levels(plot_tableINSIGF_Otus$pen_treat)

plot_tableINSIGF_Otus$pen_treatHM<-factor(plot_tableINSIGF_Otus$pen_treat, 
                                         levels = c("7_NTC", "9_NTC","19_NTC","29_NTC","31_NTC",
                                                    "10_AB","13_AB", "21_AB","23_AB","26_AB",  
                                                    "2_M","4_M","15_M", "17_M","22_M","30_M",
                                                    "3_ZnCu","5_ZnCu","14_ZnCu","18_ZnCu","28_ZnCu","32_ZnCu"))
# and for heatmap with ggplot I need to change the minumum value back to NA
plot_tableINSIGF_Otus$valueII<-plot_tableINSIGF_Otus$Abundance
min(plot_tableINSIGF_Otus$Abundance)
plot_tableINSIGF_Otus$valueII[plot_tableINSIGF_Otus$valueII==1.490935e-07]<-NA
plot_tableINSIGF_Otus$genusII<- substr(plot_tableINSIGF_Otus$genus, 1,nchar(plot_tableINSIGF_Otus$genus)-5)
# Plot:
ggplot(plot_tableINSIGF_Otus, aes(x=pen_treatHM, y=genusII, fill=valueII))+
  geom_tile()+
  scale_fill_continuous_sequential(palette =  "Light Grays", 
                                   na.value="white", name="Relative abundance")+#,
  theme(panel.background=element_rect(fill="black", colour="black")) +
  theme(axis.text.y=element_text(colour= "black",size=10))+
  theme(axis.text.x = element_text(colour="black",angle=90,size=10))+
  theme(panel.border=element_blank())+
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(legend.position="bottom")+
  theme(legend.key.size=unit(0.2, "cm"))+
  theme(legend.key.width=unit(1, "cm"))

############################

# Now the rarefied and subsampled OTUs
# Which are counts, so I cannot use gamma distribution model

#data processing:
View(noMC_phyl_OTUsubsTable_mat)
df.noMC_phyl_OTUsubsTable_matII<-as.data.frame(noMC_phyl_OTUsubsTable_mat)
df.noMC_phyl_OTUsubsTable_matII$SampleID<-rownames(noMC_phyl_OTUsubsTable_mat)
df.noMC_phyl_OTUsubsTable_matII.m<-melt(df.noMC_phyl_OTUsubsTable_matII)

colnames(df.noMC_phyl_OTUsubsTable_matII.m)<-c("pen_treat", "OTU", "Abundance")

meta_noMC_glm<-meta_noMCnmdsII[,-c(2,4,6:11)]

View(meta_noMC_glm)
Subs_OTUs_meta_glm<-merge(df.noMC_phyl_OTUsubsTable_matII.m, meta_noMC_glm, by.x="pen_treat",by.y="pen_treat.x")
Subs_OTUs_meta_tax_glm<-merge(Subs_OTUs_meta_glm, taxonomy, by="OTU")

#I want to fit an negative binomial model since the subsampled OTU data are counts
# and I have a lot of zero observations and I might have overdispersion problem.
#library(MASS)
?glm.nb 
View(Subs_OTUs_meta_tax_glm)

# I will first again try if the model works technically
nb.fit1<-glm.nb(Abundance ~ treatment, data = Subs_OTUs_meta_tax_glm, link = log)

summary( nb.fit1)

# Call:
#  glm.nb(formula = Abundance ~ treatment, data = Subs_OTUs_meta_tax_glm, 
#€        link = log, init.theta = 0.1267667015)

# Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-1.2718  -1.2718  -0.9101  -0.2556   3.0861  

#Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   4.313e+00  9.914e-02    43.5   <2e-16 ***
#  treatmentNTC -1.086e-13  1.470e-01     0.0        1    
#treatmentM   -5.685e-14  1.402e-01     0.0        1    
#treatmentAB  -7.209e-14  1.470e-01     0.0        1    
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for Negative Binomial(0.1268) family taken to be 1)

#Null deviance: 2908.7  on 2947  degrees of freedom
#Residual deviance: 2908.7  on 2944  degrees of freedom
#AIC: 21240

# Number of Fisher Scoring iterations: 1

# Theta:  0.12677 
#Std. Err.:  0.00348 

#Seems so!

# I need to check is if the deviance is below critical value
qchisq(0.95, df.residual(nb.fit1))
#[1] 3071.342
deviance(nb.fit1)
#[1] 2908.691
# It is.

# Although I will use the model for OTUs one by one, it is always useful
# to check the fit iin multiple ways.
plot(nb.fit1)
plot(residuals(nb.fit1),fitted(nb.fit1))

# Again it is possible to order the treatment levels
#Subs_OTUs_meta_tax_glm$treatment.f<-factor(Subs_OTUs_meta_tax_glm$treatment, levels = c("NTC", "AB","M","ZnCu" ))

# I will check with one OTU
OTU1_test<-subset(Subs_OTUs_meta_tax_glm, OTU=="Otu001")

nb.OTU1test<-glm.nb(Abundance ~ treatment, data = OTU1_test, link = log)

summary(nb.OTU1test)

#Call:
 # glm.nb(formula = Abundance ~ treatment, data = OTU1_test, link = log, 
      #   init.theta = 30.05371596)

#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-1.8124  -1.0566   0.1995   0.7078   1.7961  

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)   7.89922    0.07488 105.487   <2e-16 ***
#  treatmentAB  -0.22935    0.11116  -2.063   0.0391 *  
#  treatmentNTC -0.22311    0.11115  -2.007   0.0447 *  
#  treatmentM   -0.16509    0.10595  -1.558   0.1192    
#---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for Negative Binomial(30.0537) family taken to be 1)

#Null deviance: 28.032  on 21  degrees of freedom
#Residual deviance: 22.122  on 18  degrees of freedom
#AIC: 338.42

#Number of Fisher Scoring iterations: 1


#Theta:  30.05 
#Std. Err.:  9.13 
# 2 x log-likelihood:  -328.419 

#The deviance:
qchisq(0.95, df.residual(nb.OTU1test))
#[1] 28.8693
deviance(nb.OTU1test)
#[1] 22.12181

# I again want to check the fit to the distribution

#library(fitdistrplus)
#library(logspline)

#all data:
fitnbinom<-fitdist(Subs_OTUs_meta_tax_glm$Abundance,"nbinom")
plot(fitnbinom)

qqp(Subs_OTUs_meta_tax_glm$Abundance,"nbinom",
    size=fitnbinom$estimate[[1]], mu=fitnbinom$estimate[[2]])
#does not fit very well but this is all data. 

# With one OTU:

fnbinomOTU1<-fitdist(OTU1_test$Abundance,"nbinom")
plot(fnbinomOTU1)

OTU2_test<-subset(Subs_OTUs_meta_tax_glm, OTU=="Otu002")

fnbinomOTU2<-fitdist(OTU2_test$Abundance,"nbinom")
plot(fnbinomOTU2)

qqp(OTU1_test$Abundance,"nbinom",
    size=fnbinomOTU1$estimate[[1]], mu=fnbinomOTU1$estimate[[2]])
# good fit.
qqp(OTU2_test$Abundance,"nbinom",
    size=fnbinomOTU2$estimate[[1]], mu=fnbinomOTU2$estimate[[2]])
#1 outlier
OTU_test<-subset(Subs_OTUs_meta_tax_glm, OTU=="Otu003")
qqp(OTU_test$Abundance,"nbinom",
    size=fnbinomOTU2$estimate[[1]], mu=fnbinomOTU2$estimate[[2]])
##Also good here

OTU12_test<-subset(Subs_OTUs_meta_tax_glm, OTU=="Otu012")
fnbinomOTU12<-fitdist(OTU12_test$Abundance,"nbinom")
plot(fnbinomOTU12)

qqp(OTU12_test$Abundance,"nbinom",
    size=fnbinomOTU12$estimate[[1]], mu=fnbinomOTU12$estimate[[2]])
#And ok here

OTU33_test<-subset(Subs_OTUs_meta_tax_glm, OTU=="Otu033")
fnbinomOTU33<-fitdist(OTU33_test$Abundance,"nbinom")

qqp(OTU33_test$Abundance,"nbinom",
    size=fnbinomOTU33$estimate[[1]], mu=fnbinomOTU33$estimate[[2]])
# One outlier

#Now I will test the glht function

#library(multcomp)
#?glht

# all data first:
glht.nb.fit1 <- glht(nb.fit1, mcp(treatment = "Tukey"))

summary(glht( glht.nb.fit1))

#Simultaneous Tests for General Linear Hypotheses

#Linear Hypotheses:
#  Estimate Std. Error z value Pr(>|z|)
#AB - ZnCu == 0  -4.554e-15  1.470e-01       0        1
#NTC - ZnCu == 0 -1.706e-14  1.470e-01       0        1
#M - ZnCu == 0   -4.284e-14  1.402e-01       0        1
#NTC - AB == 0   -1.250e-14  1.536e-01       0        1
#M - AB == 0     -3.828e-14  1.470e-01       0        1
#M - NTC == 0    -2.578e-14  1.470e-01       0        1
#(Adjusted p values reported -- single-step method)
#Warning messages:
# 1: In chkdots(...) : Argument(s) ‘complete’ passed to ‘...’ are ignored
#2: In chkdots(...) : Argument(s) ‘complete’ passed to ‘...’ are ignored

# same with one OTU
nb.OTU1test.glht1 <- glht(nb.OTU1test, mcp(treatment = "Tukey"), p.adjust="BH")
summary(glht(nb.OTU1test.glht1))

#Simultaneous Tests for General Linear Hypotheses

#Linear Hypotheses:
 # Estimate Std. Error z value Pr(>|z|)
#AB - ZnCu == 0  -0.229347   0.111156  -2.063    0.165
#NTC - ZnCu == 0 -0.223113   0.111153  -2.007    0.185
#M - ZnCu == 0   -0.165094   0.105953  -1.558    0.403
#NTC - AB == 0    0.006234   0.116171   0.054    1.000
#M - AB == 0      0.064253   0.111205   0.578    0.939
#M - NTC == 0     0.058019   0.111203   0.522    0.954
#(Adjusted p values reported -- single-step method)

#The function:

JohannaIsTheBest_nb<-function(Otux) {
  #TEMP<-subset(EVERYTHING,OTU==Otux)
  #TEMP<-subset(targetOTUs,OTU==Otux)
  TEMP<-subset(Subs_OTUs_meta_tax_glm, OTU==Otux)
  fit_nbx <- glm.nb(Abundance ~ treatment, data = TEMP, link = log)
  glht.fit_nbx <- glht(fit_nbx, mcp(treatment = "Tukey"), p.adjust="BH")
  #print(summary(glht.mod1))
  #glhtSUM<-summary(glht(glht.mod1))
  nb.glhtSUMdf<-as.data.frame(summary(glht(glht.fit_nbx))$test[-(1:2)])
  return(nb.glhtSUMdf) 
}
#test with one OTU
nb.glhtSUMdf<-as.data.frame(JohannaIsTheBest_nb("Otu001"))
str(nb.glhtSUMdf)

#'data.frame':	6 obs. of  5 variables:
 # $ coefficients: num  -0.22935 -0.22311 -0.16509 0.00623 0.06425 ...
#$ sigma       : num  0.111 0.111 0.106 0.116 0.111 ...
#$ tstat       : num  -2.0633 -2.0073 -1.5582 0.0537 0.5778 ...
#$ pvalues     : num  0.165 0.185 0.402 1 0.939 ...
#..- attr(*, "error")= num 0.000328
#$ type        : Factor w/ 1 level "single-step": 1 1 1 1 1 1

head(unique(Subs_OTUs_meta_tax_glm$OTU))
nb_tmpOTU <- unique(Subs_OTUs_meta_tax_glm$OTU)
nb_tmpOTU_test <- lapply(nb_tmpOTU, JohannaIsTheBest_nb)
#Error in while ((it <- it + 1) < limit && abs(del) > eps) { : 
#missing value where TRUE/FALSE needed

#new kind of error but the problem is probably the same.
# I found this answer:
#http://r.789695.n4.nabble.com/glm-nb-error-td4668949.html
#seems that somewhere in pairwise comparisons there are 
#observations only in one treatment or something.
#maybe lot of zeros.

#looking for the problem:
table(Subs_OTUs_meta_tax_glm$OTU)
Subs_OTUs_meta_tax_glm_noZero<-subset(Subs_OTUs_meta_tax_glm, Abundance!=0)
dim(Subs_OTUs_meta_tax_glm_noZero)
#[1] 1668   16
dim(Subs_OTUs_meta_tax_glm)
#[1] 2948   16
# a lot of zeros.
table(Subs_OTUs_meta_tax_glm_noZero$OTU)

#"Otu132" = 0
#"Otu139" = 1
#"Otu095" = 3
#"Otu069" = 11
#"Otu080" = 7
#"Otu089" = 4
#"Otu090" = 5

JohannaIsTheBest_nb("Otu132")
#Error in while ((it <- it + 1) < limit && abs(del) > eps) { : 
# missing value where TRUE/FALSE needed
JohannaIsTheBest_nb("Otu139")
#Warning messages:
# 1: In theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace = control$trace >  :
# iteration limit reached
# 2: In sqrt(1/i) : NaNs produced
#3: In theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace = control$trace >  :
#iteration limit reached
#4: In sqrt(1/i) : NaNs produced

JohannaIsTheBest_nb("Otu095")
#Warning messages:
#1: In theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace = control$trace >  :
#iteration limit reached
#  2: In theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace = control$trace >  :
# iteration limit reached

# I wonder how many above zero observations is needed for glht not to work

JohannaIsTheBest_nb("Otu069")
#coefficients        sigma       tstat   pvalues        type
#AB - ZnCu   21.33220453 6344.9393458 0.003362082 1.0000000 single-step
#NTC - ZnCu  22.40671927 6344.9393421 0.003531432 1.0000000 single-step
#M - ZnCu    22.49980969 6344.9393376 0.003546103 1.0000000 single-step
#NTC - AB     1.07451474    0.8440909 1.272984683 0.5299104 single-step
#M - AB       1.16760516    0.8094618 1.442446244 0.4204273 single-step
#M - NTC      0.09309042    0.7798654 0.119367289 0.9992557 single-step

JohannaIsTheBest_nb("Otu080")
#coefficients     sigma       tstat   pvalues        type
#AB - ZnCu   -0.91629073 1.3517349 -0.67786276 0.9038998 single-step
#NTC - ZnCu   0.47000363 1.0378763  0.45285131 0.9685145 single-step
#M - ZnCu     0.51082562 0.9911061  0.51540961 0.9546165 single-step
#NTC - AB     1.38629436 1.3374414  1.03652715 0.7241608 single-step
#M - AB       1.42711636 1.3014814  1.09653228 0.6875573 single-step
#M - NTC      0.04082199 0.9715214  0.04201863 0.9999728 single-step
JohannaIsTheBest_nb("Otu089")

#Warning messages:
#1: In theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace = control$trace >  :
#iteration limit reached
#2: In theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace = control$trace >  :
# iteration limit reached

#So this is the warning that produces the error eventually

JohannaIsTheBest_nb("Otu090")

#coefficients        sigma        tstat   pvalues        type
#AB - NTC     -0.6931472    1.2827567 -0.540357476 0.9379592 single-step
#M - NTC       0.7339692    0.9128761  0.804018419 0.8243779 single-step
#ZnCu - NTC  -19.3862944 6344.9393396 -0.003055395 1.0000000 single-step
#M - AB        1.4271164    1.1547046  1.235914667 0.5507770 single-step
#ZnCu - AB   -18.6931472 6344.9393790 -0.002946151 1.0000000 single-step
#ZnCu - M    -20.1202635 6344.9393150 -0.003171073 1.0000000 single-step
#Warning messages:
# 1: In chkdots(...) : Argument(s) ‘complete’ passed to ‘...’ are ignored
#2: In chkdots(...) : Argument(s) ‘complete’ passed to ‘...’ are ignored

JohannaIsTheBest_nb("Otu062")
#works; 6 x > 0

JohannaIsTheBest_nb("Otu091")
#works; 7 x > 0

JohannaIsTheBest_nb("Otu056")
# no; 6 x > 0
JohannaIsTheBest_nb("Otu064")
#works; 8 x > 0

JohannaIsTheBest_nb("Otu048")
#works; 6 x > 0

table(Subs_OTUs_meta_tax_glm_noZero$OTU)

#Otu001 Otu002 Otu003 Otu004 Otu005 Otu006 Otu007 Otu008 Otu009 Otu010 Otu011 Otu012 Otu013 Otu014 Otu015 Otu016 
#22     22     22     22     22     22     22     22     22     22     22     22     22     20     22     22 
#Otu017 Otu018 Otu019 Otu020 Otu021 Otu022 Otu023 Otu024 Otu025 Otu026 Otu027 Otu028 Otu029 Otu030 Otu031 Otu032 
#22     22     22     19     22     22     22     22     22     19     22     22     19     22     22     22 
#Otu033 Otu034 Otu035 Otu036 Otu037 Otu038 Otu039 Otu040 Otu041 Otu042 Otu043 Otu044 Otu045 Otu046 Otu047 Otu048 
#22     21     22     20     21     22     16     21     21      8     22     22     21     20     22      6 
#Otu049 Otu050 Otu051 Otu052 Otu053 Otu054 Otu055 Otu056 Otu057 Otu058 Otu059 Otu060 Otu061 Otu062 Otu063 Otu064 
#5     15      4     18     19      0      1      6     22     14      4     21     17      6      7      8 
#Otu065 Otu066 Otu067 Otu068 Otu069 Otu070 Otu071 Otu072 Otu073 Otu074 Otu075 Otu076 Otu077 Otu078 Otu079 Otu080 
#11     14     20     18     11     19     11     22     10     14      6      7     20      1     17      7 
#Otu081 Otu082 Otu083 Otu084 Otu085 Otu086 Otu087 Otu088 Otu089 Otu090 Otu091 Otu092 Otu093 Otu094 Otu095 Otu096 
#22     10      5      6      8     11      8      5      4      5      7     12      5      2      3     18 
#Otu097 Otu098 Otu099 Otu100 Otu101 Otu102 Otu103 Otu104 Otu105 Otu106 Otu107 Otu108 Otu109 Otu110 Otu111 Otu112 
#8      1      2      3      2      0      3      1     18      6      9      7      1      3      1      1 
#Otu113 Otu114 Otu115 Otu116 Otu117 Otu118 Otu122 Otu123 Otu124 Otu127 Otu128 Otu129 Otu130 Otu131 Otu132 Otu134 
#10      1     16      3      1      0      7      1     10      0     10      1      0     15      0      2 
#Otu139 Otu141 Otu149 Otu158 Otu168 Otu171 
#1      5     10      2      1      0 

OTU_help_table<-as.data.frame(table(Subs_OTUs_meta_tax_glm_noZero$OTU))
# I have to just try.
# after trial and error, six times of Abundance >0 is needed.
getOTUs_nb<-subset(OTU_help_table$Var1, OTU_help_table$Freq > 6)

getOTUs_nb
getOTUs_nb<-droplevels(getOTUs_nb)

Subs_OTUs89_meta_tax_glm<-Subs_OTUs_meta_tax_glm[Subs_OTUs_meta_tax_glm$OTU %in% getOTUs_nb,]

dim(Subs_OTUs89_meta_tax_glm)
dim(Subs_OTUs_meta_tax_glm)

#Now I will try again....
JohannaIsTheBest_nb<-function(Otux) {
  
  TEMP<-subset(Subs_OTUs89_meta_tax_glm, OTU==Otux)
  fit_nbx <- glm.nb(Abundance ~ treatment, data = TEMP, link = log)
  glht.fit_nbx <- glht(fit_nbx, mcp(treatment = "Tukey"), p.adjust="BH")
  nb.glhtSUMdf<-as.data.frame(summary(glht(glht.fit_nbx))$test[-(1:2)])
  return(nb.glhtSUMdf) 
}
#testing
nb.glhtSUMdf<-as.data.frame(JohannaIsTheBest_nb("Otu001"))
str(nb.glhtSUMdf)

head(unique(Subs_OTUs89_meta_tax_glm$OTU))
nb_tmpOTU <- unique(Subs_OTUs89_meta_tax_glm$OTU)
nb_tmpOTU_test <- lapply(nb_tmpOTU, JohannaIsTheBest_nb)
#There were 50 or more warnings (use warnings() to see the first 50)
# It works!!
head(nb_tmpOTU_test)
names(nb_tmpOTU_test)
names(nb_tmpOTU_test) <- nb_tmpOTU
head(nb_tmpOTU_test)
nb_tmpOTU_test.m<- melt(nb_tmpOTU_test)
head(nb_tmpOTU_test.m)
dim(nb_tmpOTU_test.m)

levels(nb_tmpOTU_test.m$variable)
coefDF<-subset(nb_tmpOTU_test.m, variable=="coefficients")
stderDF<-subset(nb_tmpOTU_test.m, variable=="sigma")
tstatDF<-subset(nb_tmpOTU_test.m, variable=="tstat")
pvalDF<-subset(nb_tmpOTU_test.m, variable=="pvalues")
coefDF<-coefDF[,-c(1,2)]
colnames(coefDF)<-c("coefficients","OTU_c")
stderDF<-stderDF[,-c(1,2)]
colnames(stderDF)<-c("std_error","OTU_s")
tstatDF<-tstatDF[,-c(1,2)]
colnames(tstatDF)<-c("tstat","OTU_t")
pvalDF<-pvalDF[,-c(1,2)]
colnames(pvalDF)<-c("adj_p-value","OTU_p")

eighty9OTUnb_results<-cbind(coefDF,stderDF,tstatDF,pvalDF)
#comparisons
nb_comp_df<-as.data.frame(lapply(nb_tmpOTU, JohannaIsTheBest_nb))
nb_comparisons<-rownames(nb_comp_df)
nb_comparisons

dim(eighty9OTUnb_results)
#[1] 534   8
534/89 
#[1] 6

eighty9OTUnb_results$comparison<-rep(nb_comparisons,89)
eighty9OTUnb_results$OTU<-eighty9OTUnb_results$OTU_c
eighty9OTUnb_results<-eighty9OTUnb_results[,-c(2,4,6,8)]

colnames(eighty9OTUnb_results)<-c("Delta.Estimate", "Std.Error", "z-value", "p.adjusted","Comparison", "OTU")
rownames(eighty9OTUnb_results)<-interaction(eighty9OTUnb_results$Comparison,eighty9OTUnb_results$OTU, sep="_")
str(eighty9OTUnb_results)
eighty9OTUnb_results_sigP<-subset(eighty9OTUnb_results, p.adjusted <= 0.05)

#Fold changes:

eighty9OTUnb_results_sigP$Fold.Change<-exp(eighty9OTUnb_results_sigP$Delta.Estimate)

mean(with (Subs_OTUs89_meta_tax_glm, subset(Abundance, OTU=="Otu002" & treatment=="M")))
#[1] 497.3333
mean(with (Subs_OTUs89_meta_tax_glm, subset(Abundance, OTU=="Otu002" & treatment=="ZnCu")))
#[1] 742.8333
497.3333/742.8333
#[1] 0.6695086
#Corresponding number in eighty9OTUnb_results_sigP$Fold.Change is 0.6695086

unique(eighty9OTUnb_results_sigP$OTU)

#[1] "Otu002" "Otu006" "Otu009" "Otu010" "Otu014" "Otu017" "Otu018" "Otu019" "Otu024" "Otu025" "Otu027"
#[12] "Otu029" "Otu033" "Otu034" "Otu036" "Otu040" "Otu044" "Otu050" "Otu052" "Otu053" "Otu065" "Otu068"
#[23] "Otu077" "Otu085" "Otu086" "Otu096" "Otu115" "Otu131"

#Save NB results:
eighty9OTUnb_results_sigP_forXLSx <- merge(eighty9OTUnb_results_sigP, taxonomy, by = 'OTU')
eighty9OTUnb_results_sigP_forXLSx$genus<- substr(eighty9OTUnb_results_sigP_forXLSx$genus, 1,nchar(eighty9OTUnb_results_sigP_forXLSx$genus)-5)

write.table(eighty9OTUnb_results_sigP_forXLSx, file="eighty9OTUnb_results_sigP_forXLSx.txt",sep="\t")

#Plots
plot_tableNegBinOTU<- Subs_OTUs89_meta_tax_glm[Subs_OTUs89_meta_tax_glm$OTU %in% eighty9OTUnb_results_sigP$OTU,]

unique(plot_tableNegBinOTU$OTU)
unique(eighty9OTUnb_results_sigP$OTU)

plot_tableNegBinOTU$OTU<-droplevels(plot_tableNegBinOTU$OTU)

str(plot_tableNegBinOTU)
levels(plot_tableNegBinOTU$genus)
plot_tableNegBinOTU$genus.f<-as.factor(plot_tableNegBinOTU$genus)
levels(plot_tableNegBinOTU$genus.f)

plot_tableNegBinOTU_o<-plot_tableNegBinOTU[order(-plot_tableNegBinOTU$Abundance),]
unique(plot_tableNegBinOTU_o$genus)

plot_tableNegBinOTU_o$Genus<-factor(plot_tableNegBinOTU_o$genus, levels=unique(plot_tableNegBinOTU_o$genus))

levels(plot_tableNegBinOTU_o$Genus)

str(plot_tableNegBinOTU_o)


dif_abundSUBS<-ggplot(plot_tableNegBinOTU_o, aes(x = Genus, y = Abundance, fill = treatment)) + 
  geom_bar(position="dodge",stat = "identity") +
  #facet_grid(treatment ~ ., scales="free")
  #scale_fill_manual(values = HUE_16subs) +
  # Remove x axis title
  scale_fill_manual(values=c("#aeb646","#924e7d","#8b58c8","#749b9d")) +
  theme_bw()+
  scale_x_discrete(limits = rev(levels(plot_tableNegBinOTU_o$Genus)))+
  theme(axis.title.x = element_blank()) + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1200)) +
  #ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = T, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  coord_flip()

dif_abundSUBS

# Boxplot:
ggplot(plot_tableNegBinOTU_o, aes(x = Genus, y =Abundance)) + 
  geom_boxplot (aes(color=treatment)) +
  #facet_grid(treatment ~ ., scales="free")
  #scale_fill_manual(values = HUE_16subs) +
  # Remove x axis title
  scale_color_manual(values=c("#924e7d","#749b9d","#8b58c8","#aeb646")) +
  theme_bw()+
  #scale_x_discrete(limits = rev(levels(plot_tableGammaOTU_o$Genus)))+
  scale_x_discrete(limits = levels(plot_tableNegBinOTU_o$Genus))+
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 0.95, vjust = 1)) +
  theme(axis.text.x=element_text(colour="black")) +
  #scale_y_log10(expand =c(0, 0.1), limits = c(1e-06, 0.36))+
  scale_y_continuous(expand = c(0.1, 0.1), limits = c(1, 1200)) +
  theme(axis.text.y=element_text(colour="black")) +
  #ylim(c(0,1)) +
  guides(color = guide_legend(title="Treatment", keywidth = .5, keyheight = .5, ncol = 1)) 

#Warning message:
 # Removed 82 rows containing non-finite values (stat_boxplot). 

plot_tableNegBinOTU_o$Genus2<- substr(plot_tableNegBinOTU_o$genus, 1,nchar( plot_tableNegBinOTU_o$genus)-5)
plot_tableNegBinOTU_o$Genus2<-factor(plot_tableNegBinOTU_o$Genus2, levels=unique(plot_tableNegBinOTU_o$Genus2))

plot_tableNegBinOTU_o_abund14<- plot_tableNegBinOTU_o[plot_tableNegBinOTU_o$Genus2 %in% levels(plot_tableNegBinOTU_o$Genus2)[1:14],]
plot_tableNegBinOTU_o_abund14$Genus2<-droplevels(plot_tableNegBinOTU_o_abund14$Genus2)

plot_tableNegBinOTU_o_abund14$treatment.f<-factor(plot_tableNegBinOTU_o_abund14$treatment, 
                                          levels = c("NTC", "AB", "M","ZnCu"))
# 

OTUs14Subs<-ggplot(plot_tableNegBinOTU_o_abund14, aes(x = Genus2, y =Abundance)) + 
  geom_boxplot (aes(color=treatment.f)) +
  #facet_grid(treatment ~ ., scales="free")
  #scale_fill_manual(values = HUE_16subs) +
  # Remove x axis title
  scale_color_manual(values=c("#924e7d","#749b9d","#8b58c8","#aeb646")) +
  theme_bw()+
  #scale_x_discrete(limits = rev(levels(plot_tableGammaOTU_o$Genus)))+
  #scale_x_discrete(limits = levels(plot_tableGammaOTU_abund14$Genus2))+
  scale_x_discrete(limits = rev(levels(plot_tableNegBinOTU_o_abund14$Genus2)))+
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0.95, vjust = 1)) +
  theme(axis.text.x=element_text(colour="black")) +
  scale_y_continuous(expand = c(0, 0.01), limits = c(0, 1200)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 0.36)) +
  theme(axis.text.y=element_text(colour="black")) +
  #ylim(c(0,1)) +
  guides(color = guide_legend(title="Treatment", keywidth = .5, keyheight = .5, ncol = 1)) 
OTUs14Subs

#plot friure 3 a b c:
grid.arrange(OTUs14RA,OTUs14Subs,Dif_abund_ARGs10,nrow=3,ncol=1)

#Now the OTU heatmaps
# Significant OTUs

View(plot_tableNegBinOTU_o)


#And I want to also order the samples
levels(plot_tableNegBinOTU_o$pen_treat)

plot_tableNegBinOTU_o$pen_treatHM<-factor(plot_tableNegBinOTU_o$pen_treat, 
                                         levels = c("7_NTC", "9_NTC","19_NTC","29_NTC","31_NTC",
                                                    "10_AB","13_AB", "21_AB","23_AB","26_AB",  
                                                    "2_M","4_M","15_M", "17_M","22_M","30_M",
                                                    "3_ZnCu","5_ZnCu","14_ZnCu","18_ZnCu","28_ZnCu","32_ZnCu"))

min(plot_tableNegBinOTU_o$Abundance)
#0
plot_tableNegBinOTU_o$valueII<-plot_tableNegBinOTU_o$Abundance
plot_tableNegBinOTU_o$valueII[plot_tableNegBinOTU_o$valueII==0]<-NA

plot_tableNegBinOTU_o$genusII<- substr(plot_tableNegBinOTU_o$genus, 1,nchar(plot_tableNegBinOTU_o$genus)-5)

# Now I can plot
ggplot(plot_tableNegBinOTU_o, aes(x=pen_treatHM, y=genusII, fill=valueII))+
  geom_tile()+
  scale_fill_continuous_sequential(palette =  "Lajolla",
   #I don't want the most darkest color and I want to reverse the colors
    begin = 1, end = 0.1, na.value="white", name="Relative abundance")+
  theme(panel.background=element_rect(fill="black", colour="black")) +
  theme(axis.text.y=element_text(colour= "black",size=10))+
  theme(axis.text.x = element_text(colour="black",angle=90,size=8))+
  theme(panel.border=element_blank())+
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(legend.position="bottom")+
  theme(legend.key.size=unit(0.2, "cm"))+
  theme(legend.key.width=unit(1, "cm"))


#Now the Otus that didnt have significant differences:
plot_table_NB_INSIGF_Otus<- Subs_OTUs_meta_tax_glm[!(Subs_OTUs_meta_tax_glm$OTU %in% eighty9OTUnb_results_sigP$OTU),]
unique(plot_table_NB_INSIGF_Otus$OTU)
#[1] Otu001 Otu003 Otu004 Otu005 Otu007 Otu008 Otu011 Otu012 Otu013 Otu015 Otu016 Otu020 Otu021 Otu022 Otu023 Otu026 Otu028 Otu030 Otu031 Otu032 Otu035
#[22] Otu037 Otu038 Otu039 Otu041 Otu042 Otu043 Otu045 Otu046 Otu047 Otu048 Otu049 Otu051 Otu054 Otu055 Otu056 Otu057 Otu058 Otu059 Otu060 Otu061 Otu062
#[43] Otu063 Otu064 Otu066 Otu067 Otu069 Otu070 Otu071 Otu072 Otu073 Otu074 Otu075 Otu076 Otu078 Otu079 Otu080 Otu081 Otu082 Otu083 Otu084 Otu087 Otu088
#[64] Otu089 Otu090 Otu091 Otu092 Otu093 Otu094 Otu095 Otu097 Otu098 Otu099 Otu100 Otu101 Otu102 Otu103 Otu104 Otu105 Otu106 Otu107 Otu108 Otu109 Otu110
#[85] Otu111 Otu112 Otu113 Otu114 Otu116 Otu117 Otu118 Otu122 Otu123 Otu124 Otu127 Otu128 Otu129 Otu130 Otu132 Otu134 Otu139 Otu141 Otu149 Otu158 Otu168
#[106] Otu171
#134 Levels: Otu001 Otu002 Otu003 Otu004 Otu005 Otu006 Otu007 Otu008 Otu009 Otu010 Otu011 Otu012 Otu013 Otu014 Otu015 Otu016 Otu017 Otu018 Otu019 ... Otu171
#And I want to order the samples

levels(plot_table_NB_INSIGF_Otus$pen_treat)

plot_table_NB_INSIGF_Otus$pen_treatHM<-factor(plot_table_NB_INSIGF_Otus$pen_treat, 
                                          levels = c("7_NTC", "9_NTC","19_NTC","29_NTC","31_NTC",
                                                     "10_AB","13_AB", "21_AB","23_AB","26_AB",  
                                                     "2_M","4_M","15_M", "17_M","22_M","30_M",
                                                     "3_ZnCu","5_ZnCu","14_ZnCu","18_ZnCu","28_ZnCu","32_ZnCu"))
# and for heatmap with ggplot I need to change the minumum value back to NA
plot_table_NB_INSIGF_Otus$valueII<-plot_table_NB_INSIGF_Otus$Abundance
min(plot_table_NB_INSIGF_Otus$Abundance)
plot_table_NB_INSIGF_Otus$valueII[plot_table_NB_INSIGF_Otus$valueII==0]<-NA
plot_table_NB_INSIGF_Otus$genusII<- substr(plot_table_NB_INSIGF_Otus$genus, 1,nchar(plot_table_NB_INSIGF_Otus$genus)-5)
plot_table_NB_INSIGF_OtusII<-subset(plot_table_NB_INSIGF_Otus, genusII!="Weissella")


# Plot:
ggplot(plot_table_NB_INSIGF_OtusII, aes(x=pen_treatHM, y=genusII, fill=valueII))+
  geom_tile()+
  scale_fill_continuous_sequential(palette =  "Light Grays", 
  begin = 0.05, end = 1, na.value="white", name="Relative abundance")+#,
  theme(panel.background=element_rect(fill="black", colour="black")) +
  theme(axis.text.y=element_text(colour= "black",size=7))+
  theme(axis.text.x = element_text(colour="black",angle=90,size=10))+
  theme(panel.border=element_blank())+
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(legend.position="bottom")+
  theme(legend.key.size=unit(0.2, "cm"))+
  theme(legend.key.width=unit(1, "cm"))


###############################################################

   ### Comparison of ARG & MGE and 16S ordinations ###

###############################################################
#        Both ARG and OTU data need to be analyzed            #
#                   for the next section                      #
###############################################################

# We want to compare ordinations from 16S OTUs and ARGs and MGEs
# For this we need to have matching matrices and ordinations.

#I will first compare ARGs and TSS OTUs (relative abundance):
NMDS1$points  #ARGs and MGEs 
mdsIII.f$points #OTUs, TSS normalized

#"2_M " is missing from NMDS1 (ARGs and MGEs)
#  "11_NTC" is missing from mdsIII.f (TSS OTUs)
# "16_AB" is missing from mdsIII.f (TSS OTUs)

# So I neeed to remove some samples from matrices and do again the NMDS runs.

# 1) 16 S OTUs: 
View(noMC_phylOTUtable_RA_NoZero_mat)

#So this is easy: 2_M " is the first row
noMC_phylOTUtable_RA_NoZero_matPRO<-noMC_phylOTUtable_RA_NoZero_mat
noMC_phylOTUtable_RA_NoZero_matPRO<-noMC_phylOTUtable_RA_NoZero_matPRO[-1,]

View(results_relative_mat)
rownames(results_relative_mat)
#"11_NTC" is the 7th row and "16_AB" is the 11th row
results_relative_mat_PRO<-results_relative_mat
results_relative_mat_PRO<-results_relative_mat_PRO[-c(7,11),]
dim(noMC_phylOTUtable_RA_NoZero_matPRO)
dim(results_relative_mat_PRO)

ARG_nmdsPRO<- metaMDS(results_relative_mat_PRO)
OTUrel_nmdsPRO<- metaMDS(noMC_phylOTUtable_RA_NoZero_matPRO)

ARG_nmdsPRO
#Call:
# metaMDS(comm = results_relative_mat_PRO) 

#global Multidimensional Scaling using monoMDS

#Data:     results_relative_mat_PRO 
#Distance: bray 

#Dimensions: 2 
#Stress:      0.1562437 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘results_relative_mat_PRO’

OTUrel_nmdsPRO
#Call:
# metaMDS(comm = noMC_phylOTUtable_RA_NoZero_matPRO) 

#global Multidimensional Scaling using monoMDS

#Data:     noMC_phylOTUtable_RA_NoZero_matPRO 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.153237 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘noMC_phylOTUtable_RA_NoZero_matPRO’ 

ProCruARG_OTUrel_sym<-procrustes(ARG_nmdsPRO, OTUrel_nmdsPRO, symmertric=TRUE)


ProCruARG_OTUrel_sym

#Call:
# procrustes(X = ARG_nmdsPRO, Y = OTUrel_nmdsPRO, symmertric = TRUE) 

#Procrustes sum of squares:
# 1.084 

summary(ProCruARG_OTUrel_sym)

#Call:
# procrustes(X = ARG_nmdsPRO, Y = OTUrel_nmdsPRO, symmertric = TRUE) 

#Number of objects: 21    Number of dimensions: 2 

#Procrustes sum of squares:  
#  1.084353 
#Procrustes root mean squared error: 
# 0.2272352 
#Quantiles of Procrustes errors:
# Min         1Q     Median         3Q        Max 
#0.06696301 0.10151114 0.18288304 0.29352587 0.36511343 

#Rotation matrix:
#  [,1]        [,2]
#[1,] 0.9996394 -0.0268534
#[2,] 0.0268534  0.9996394

#Translation of averages:
# [,1]         [,2]
#[1,] 4.433795e-18 5.094463e-18

#Scaling of target:
#  0.4546066

plot(ProCruARG_OTUrel_sym)
str(ProCruARG_OTUrel_sym)

ProCruARG_OTUrel_sym$Yrot #OTUs
ProCruARG_OTUrel_sym$X #ARGs

plot(0,0, xlim=range(c(ProCruARG_OTUrel_sym$X[,1], ProCruARG_OTUrel_sym$Yrot[,1])),
ylim=range(c(ProCruARG_OTUrel_sym$X[,2], ProCruARG_OTUrel_sym$Yrot[,2])), type="n")
points(ProCruARG_OTUrel_sym$X, pch=17, col="red")
points(ProCruARG_OTUrel_sym$Yrot, pch=19, col="blue")
with(ProCruARG_OTUrel_sym, arrows(Yrot[,1], Yrot[,2], X[,1], X[,2]))

#Lets make a metadata table for nicer colors
meta_forPRO<-data.frame(X1=ProCruARG_OTUrel_sym$X[,1],X2=ProCruARG_OTUrel_sym$X[,2],
      Y1= ProCruARG_OTUrel_sym$Yrot[,1],Y2=ProCruARG_OTUrel_sym$Yrot[,2])
meta_forPRO$treatment<-W_metadata$treatment[-c(7,11)]
meta_forPRO$color<-W_metadata$color[-c(7,11)]
meta_forPRO$pen_treat<-W_metadata$pen_treat[-c(7,11)]

plot(ProCruARG_OTUrel_sym, xlim=range(c(ProCruARG_OTUrel_sym$X[,1], 
    ProCruARG_OTUrel_sym$Yrot[,1])), ylim=range(c(ProCruARG_OTUrel_sym$X[,2], 
    ProCruARG_OTUrel_sym$Yrot[,2])), type="n", main= "Procrustes errors, ARGs vs OTUs_RA")
points(ProCruARG_OTUrel_sym$X, pch=17, col=meta_forPRO$color, cex=1.7)
points(ProCruARG_OTUrel_sym$Yrot, pch=19,col =meta_forPRO$color, cex=1.7)
with(ProCruARG_OTUrel_sym, segments(Yrot[,1], Yrot[,2], X[,1], X[,2],col =meta_forPRO$color))
#residuals plot
plot(ProCruARG_OTUrel_sym, kind=2, col=meta_forPRO$color, lwd=4, main= "Procrustes errors, ARGs vs OTUs_RA")
text(residuals(ProCruARG_OTUrel_sym),labels=meta_forPRO$pen_treat,cex =0.7 )

#statistics:
ProARG_OTUrelTEST<-protest(ARG_nmdsPRO, OTUrel_nmdsPRO, scores="sites", permutations=9999)

ProARG_OTUrelTEST

#Call:
# protest(X = ARG_nmdsPRO, Y = OTUrel_nmdsPRO, scores = "sites",      permutations = 9999) 

#Procrustes Sum of Squares (m12 squared):         0.8341
#Correlation in a symmetric Procrustes rotation:  0.4073 
#Significance: 0.0639 

#Permutation: free
#Number of permutations: 9999

##############################################


#Same with rarefied & subsampled OTU data: 

View(noMC_phyl_OTUsubsTable_mat)

noMC_phyl_OTUsubsTable_matPRO<-noMC_phyl_OTUsubsTable_mat
noMC_phyl_OTUsubsTable_matPRO<-noMC_phyl_OTUsubsTable_matPRO[-1,]

OTUsubs_nmdsPRO<- metaMDS(noMC_phyl_OTUsubsTable_matPRO)

OTUsubs_nmdsPRO

#Call:
# metaMDS(comm = noMC_phyl_OTUsubsTable_matPRO) 

#global Multidimensional Scaling using monoMDS

#Data:     wisconsin(sqrt(noMC_phyl_OTUsubsTable_matPRO)) 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.1569223 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘wisconsin(sqrt(noMC_phyl_OTUsubsTable_matPRO))’ 


ProCruARG_OTUsubs_sym<-procrustes(ARG_nmdsPRO, OTUsubs_nmdsPRO, symmertric=TRUE)

summary(ProCruARG_OTUsubs_sym)

#Call:
# procrustes(X = ARG_nmdsPRO, Y = OTUsubs_nmdsPRO, symmertric = TRUE) 

#Number of objects: 21    Number of dimensions: 2 

#Procrustes sum of squares:  
#  0.8784061  
#Procrustes root mean squared error: 
#  0.2045211
#Quantiles of Procrustes errors:
# Min         1Q     Median         3Q        Max 
#0.01464197 0.09374667 0.15602916 0.22859689 0.38861040 

#Rotation matrix:
# [,1]      [,2]
#[1,]   0.8994692 0.4369841
#[2,] -0.4369841 0.8994692

#Translation of averages:
# [,1]          [,2]
#[1,] 1.626164e-18 3.530026e-18

#Scaling of target:
#  [1] 0.6794088

ProARG_OTUsubsTEST<-protest(ARG_nmdsPRO, OTUsubs_nmdsPRO, scores="sites", permutations=9999)

ProARG_OTUsubsTEST
#Call:
# protest(X = ARG_nmdsPRO, Y = OTUsubs_nmdsPRO, scores = "sites",      permutations = 9999) 

#Procrustes Sum of Squares (m12 squared):        0.6757  
#Correlation in a symmetric Procrustes rotation: 0.5695 
#Significance: 0.0012 

#Permutation: free
#Number of permutations: 9999

# Procrustes ordination
plot(ProCruARG_OTUsubs_sym, xlim=range(c(ProCruARG_OTUsubs_sym$X[,1], 
ProCruARG_OTUsubs_sym$Yrot[,1])), ylim=range(c(ProCruARG_OTUsubs_sym$X[,2],
ProCruARG_OTUsubs_sym$Yrot[,2])), type="n", main= "Procrustes errors, ARGs vs OTUs_subs")
points(ProCruARG_OTUsubs_sym$X, pch=17, col=meta_forPRO$color, cex=1.7)
points(ProCruARG_OTUsubs_sym$Yrot, pch=15,col =meta_forPRO$color, cex=1.7)
with(ProCruARG_OTUsubs_sym, segments(Yrot[,1], Yrot[,2], X[,1], X[,2],col =meta_forPRO$color))

# Residual errors:
plot(ProCruARG_OTUsubs_sym, kind=2, col=meta_forPRO$color, 
     lwd=4, main= "Procrustes errors, ARGs vs OTUs_subs")
text(residuals(ProCruARG_OTUsubs_sym),labels=meta_forPRO$pen_treat,cex =0.7 )

##########################################################################################

# I want to also compare the relative abundace OTUs and rarefied and subsampled OTUs:

ProCruOTUrel_OTUsubs_sym<-procrustes(mdsII, mdsIII.f, symmertric=TRUE)

summary(ProCruOTUrel_OTUsubs_sym)
#Call:
 # procrustes(X = mdsII, Y = mdsIII.f, symmertric = TRUE) 

#Number of objects: 22    Number of dimensions: 2 

#Procrustes sum of squares:  
 #  0.4970521 
#Procrustes root mean squared error: 
 # 0.1503106 
#Quantiles of Procrustes errors:
 # Min         1Q     Median         3Q        Max 
#0.01213190 0.06704475 0.11477361 0.17310909 0.35255725 

#Rotation matrix:
 # [,1]        [,2]
#[1,] -0.99865311 -0.05188416
#[2,] -0.05188416  0.99865311

#Translation of averages:
 # [,1]          [,2]
#[1,] -4.733838e-20 -2.100613e-18

#Scaling of target:
 # [1] 0.6868032

ProOTUrel_OTUsubsTEST<-protest(mdsII, mdsIII.f, scores="sites", permutations=9999)

ProOTUrel_OTUsubsTEST

#Call:
# protest(X = mdsII, Y = mdsIII.f, scores = "sites", permutations = 9999) 

#Procrustes Sum of Squares (m12 squared):        0.5079 
#Correlation in a symmetric Procrustes rotation: 0.7015 
#Significance:  2e-04 

#Permutation: free
#Number of permutations: 9999

#metadata for plots:
meta_forOTU_PRO<-data.frame(X1=ProCruOTUrel_OTUsubs_sym$X[,1],X2=ProCruOTUrel_OTUsubs_sym$X[,2],
  Y1= ProCruOTUrel_OTUsubs_sym$Yrot[,1],Y2=ProCruOTUrel_OTUsubs_sym$Yrot[,2])
meta_forOTU_PRO$pen_treat<-meta_noMCnmdsII$pen_treat.x
meta_forOTU_PRO$color<-meta_noMCnmdsII$color
meta_forOTU_PRO$treatment<-meta_noMCnmdsII$treatment

#Procrustes ordination plot:
plot(ProCruOTUrel_OTUsubs_sym, xlim=range(c(ProCruOTUrel_OTUsubs_sym$X[,1], ProCruOTUrel_OTUsubs_sym$Yrot[,1])), 
     ylim=range(c(ProCruOTUrel_OTUsubs_sym$X[,2], ProCruOTUrel_OTUsubs_sym$Yrot[,2])), type="n",
     main="Procrustes errors, OTUs_subs vs OTUs_rel")
points(ProCruOTUrel_OTUsubs_sym$X, pch=15, col=meta_forOTU_PRO$color, cex=1.7)
points(ProCruOTUrel_OTUsubs_sym$Yrot, pch=19,col =meta_forOTU_PRO$color, cex=1.7)
with(ProCruOTUrel_OTUsubs_sym, segments(Yrot[,1], Yrot[,2], X[,1], X[,2],col =meta_forOTU_PRO$color))

#residual errors:
plot(ProCruOTUrel_OTUsubs_sym, kind=2, col=meta_forOTU_PRO$color, lwd=4,main="Procrustes errors, OTUs_subs vs OTUs_rel")
text(residuals(ProCruOTUrel_OTUsubs_sym),labels=meta_forOTU_PRO$pen_treat,cex =0.7 )

#######################################################
###                  Mantel tests                   ###
######################################################

#The idea is to examine if ecological distance matrices correlate.
# In other words, if the other explains the other.

#Distances (Bray-Curtis):
ARGandMGE_dist<-vegdist(results_relative_mat_PRO)
OTU_RA_dist<-vegdist(noMC_phylOTUtable_RA_NoZero_matPRO)
OTU_subs_dist<-vegdist(noMC_phyl_OTUsubsTable_matPRO)

#lets see what the relationships look like:
plot(ARGandMGE_dist, OTU_RA_dist)
plot(ARGandMGE_dist, OTU_subs_dist)

# I will use Spearman's correlation, TSS OTUs first:
mantel(ARGandMGE_dist, OTU_RA_dist, method = "spearman")

#Mantel statistic based on Spearman's rank correlation rho 

#Call:
#mantel(xdis = ARGandMGE_dist, ydis = OTU_RA_dist, method = "spearman") 

#Mantel statistic r: 0.2462 
#Significance: 0.012 

#Upper quantiles of permutations (null model):
#90%   95% 97.5%   99% 
#0.139 0.177 0.204 0.251 
#Permutation: free
#Number of permutations: 999

#Significant correlation but quite low.

#Same with subsampled OTUs: 

mantel(ARGandMGE_dist, OTU_subs_dist, method = "spearman")

#Mantel statistic based on Spearman's rank correlation rho 

#Call:
#mantel(xdis = ARGandMGE_dist, ydis = OTU_subs_dist, method = "spearman") 

#Mantel statistic r: 0.2343 
 #     Significance:  0.016 

#Upper quantiles of permutations (null model):
 # 90%   95% 97.5%   99% 
#0.141 0.178 0.205 0.236 
#Permutation: free
#Number of permutations: 999

#Significant correlation but quite low.

#I will compare distance matrices of each treatment groups with Mantels test. 

#First I want to have NTC's in separate matrices:
NTC_ARG_MGE_mat<-results_relative_mat_PRO[c("7_NTC","9_NTC","19_NTC","29_NTC","31_NTC"), ]
NTC_OTU_RA_mat<-noMC_phylOTUtable_RA_NoZero_matPRO[c("7_NTC","9_NTC","19_NTC","29_NTC","31_NTC"), ]
NTC_OTU_Subs_mat<-noMC_phyl_OTUsubsTable_matPRO[c("7_NTC","9_NTC","19_NTC","29_NTC","31_NTC"), ]
# distances
NTC_ARG_MGE_dist <- vegdist(NTC_ARG_MGE_mat)
NTC_OTU_RA_dist <- vegdist(NTC_OTU_RA_mat)
NTC_OTU_Subs_dist <- vegdist(NTC_OTU_Subs_mat)

#ARGs and TSS OTUs
mantel(NTC_ARG_MGE_dist, NTC_OTU_RA_dist, method = "spearman")

#Mantel statistic based on Spearman's rank correlation rho 

#Call:
#mantel(xdis = NTC_ARG_MGE_dist, ydis = NTC_OTU_RA_dist, method = "spearman") 

#Mantel statistic r: 0.6606 
#Significance: 0.025 

#Upper quantiles of permutations (null model):
#90%   95% 97.5%   99% 
#0.411 0.518 0.566 0.662 
#Permutation: free
#Number of permutations: 119

mantel(NTC_ARG_MGE_dist, NTC_OTU_Subs_dist, method = "spearman")
#'nperm' >= set of all permutations: complete enumeration.
#Set of permutations < 'minperm'. Generating entire set.


#Mantel statistic based on Spearman's rank correlation rho 

#Call:
#mantel(xdis = NTC_ARG_MGE_dist, ydis = NTC_OTU_Subs_dist, method = "spearman",      permutations = 999) 

#Mantel statistic r: 0.6242 
#Significance: 0.025 

#Upper quantiles of permutations (null model):
#90%   95% 97.5%   99% 
#0.442 0.516 0.588 0.620 
#Permutation: free
#Number of permutations: 119

#These have stronger correlation

#Carbadox group:
AB_ARG_MGE_mat<-results_relative_mat_PRO[c("10_AB","13_AB","21_AB","23_AB","26_AB"), ]
AB_OTU_RA_mat<-noMC_phylOTUtable_RA_NoZero_matPRO[c("10_AB","13_AB","21_AB","23_AB","26_AB"), ]
AB_OTU_Subs_mat<-noMC_phyl_OTUsubsTable_matPRO[c("10_AB","13_AB","21_AB","23_AB","26_AB"), ]

AB_ARG_MGE_dist <- vegdist(AB_ARG_MGE_mat)
AB_OTU_RA_dist <- vegdist(AB_OTU_RA_mat)
AB_OTU_Subs_dist <- vegdist(AB_OTU_Subs_mat)

mantel(AB_ARG_MGE_dist, AB_OTU_RA_dist, method = "spearman")

#Mantel statistic based on Spearman's rank correlation rho 

#Call:
#mantel(xdis = AB_ARG_MGE_dist, ydis = AB_OTU_RA_dist, method = "spearman") 

#Mantel statistic r: -0.3697 
#Significance: 0.925 

#Upper quantiles of permutations (null model):
#90%   95% 97.5%   99% 
#0.321 0.440 0.588 0.620 
#Permutation: free
#Number of permutations: 119

mantel(AB_ARG_MGE_dist, AB_OTU_Subs_dist, method = "spearman")

#Mantel statistic based on Spearman's rank correlation rho 

#Call:
#mantel(xdis = AB_ARG_MGE_dist, ydis = AB_OTU_Subs_dist, method = "spearman") 

#Mantel statistic r: -0.4667 
#Significance: 0.94167 

#Upper quantiles of permutations (null model):
#90%   95% 97.5%   99% 
#0.370 0.467 0.479 0.511 
#Permutation: free
#Number of permutations: 119

#Not significant.

#Mushroom group:
M_ARG_MGE_mat<-results_relative_mat_PRO[c("4_M","15_M","17_M","22_M","30_M"), ]
M_OTU_RA_mat<-noMC_phylOTUtable_RA_NoZero_matPRO[c("4_M","15_M","17_M","22_M","30_M"), ]
M_OTU_Subs_mat<-noMC_phyl_OTUsubsTable_matPRO[c("4_M","15_M","17_M","22_M","30_M"), ]

M_ARG_MGE_dist <- vegdist(M_ARG_MGE_mat)
M_OTU_RA_dist <- vegdist(M_OTU_RA_mat)
M_OTU_Subs_dist <- vegdist(M_OTU_Subs_mat)

mantel(M_ARG_MGE_dist, M_OTU_RA_dist, method = "spearman")

#Mantel statistic based on Spearman's rank correlation rho 

#Call:
#mantel(xdis = M_ARG_MGE_dist, ydis = M_OTU_RA_dist, method = "spearman") 

#Mantel statistic r:  -0.1515 
#Significance: 0.68333 

#Upper quantiles of permutations (null model):
#90%   95% 97.5%   99% 
# 0.527 0.564 0.591 0.688 
#Permutation: free
#Number of permutations: 119

mantel(M_ARG_MGE_dist, M_OTU_Subs_dist, method = "spearman")

#Mantel statistic based on Spearman's rank correlation rho 

#Call:
#mantel(xdis = M_ARG_MGE_dist, ydis = M_OTU_Subs_dist, method = "spearman") 

#Mantel statistic r: -0.1636 
#Significance: 0.675 

#Upper quantiles of permutations (null model):
#90%   95% 97.5%   99% 
#0.491 0.553 0.577 0.689 
#Permutation: free
#Number of permutations: 119

# Not significant
# Zinc and copper:
ZnCu_ARG_MGE_mat<-results_relative_mat_PRO[c("3_ZnCu","5_ZnCu","14_ZnCu","18_ZnCu","28_ZnCu","32_ZnCu"), ]
ZnCu_OTU_RA_mat<-noMC_phylOTUtable_RA_NoZero_matPRO[c("3_ZnCu","5_ZnCu","14_ZnCu","18_ZnCu","28_ZnCu","32_ZnCu"), ]
ZnCu_OTU_Subs_mat<-noMC_phyl_OTUsubsTable_matPRO[c("3_ZnCu","5_ZnCu","14_ZnCu","18_ZnCu","28_ZnCu","32_ZnCu"), ]

ZnCu_ARG_MGE_dist <- vegdist(ZnCu_ARG_MGE_mat)
ZnCu_OTU_RA_dist <- vegdist(ZnCu_OTU_RA_mat)
ZnCu_OTU_Subs_dist <- vegdist(ZnCu_OTU_Subs_mat)

mantel(ZnCu_ARG_MGE_dist, ZnCu_OTU_RA_dist, method = "spearman")

#Mantel statistic based on Spearman's rank correlation rho 

#Call:
#mantel(xdis = ZnCu_ARG_MGE_dist, ydis = ZnCu_OTU_RA_dist, method = "spearman") 

#Mantel statistic r: 0.3571 
#Significance: 0.14167 

#Upper quantiles of permutations (null model):
#90%   95% 97.5%   99% 
#0.429 0.536 0.593 0.686 
#Permutation: free
#Number of permutations: 719

mantel(ZnCu_ARG_MGE_dist, ZnCu_OTU_Subs_dist, method = "spearman")

#Mantel statistic based on Spearman's rank correlation rho 

#Call:
#mantel(xdis = ZnCu_ARG_MGE_dist, ydis = ZnCu_OTU_Subs_dist, method = "spearman") 

#Mantel statistic r: 0.3964 
#Significance: 0.12222 

#Upper quantiles of permutations (null model):
#90%   95% 97.5%   99% 
#0.436 0.532 0.589 0.676 
#Permutation: free
#Number of permutations: 719

#This would have worked as well...
#as.matrix(noMC_phylOTUtable_RA_NoZero_mat[grep("NTC",rownames(noMC_phylOTUtable_RA_NoZero_mat)),])

# Also these are not significant. 

# So OTUs explain ARGs and MGEs in NTC groups but on in growth promoter groups

##################################################################################
#                                                                                #
#    Network analysis within treatment groups (Visualization with Gephi, not R)  #
#                                                                                #
##################################################################################

# The goal is to do correlations between ARGs and MGEs.
# Challeneges are that the correlation matrix will not 
# be symmetrical and I want to have only those assays that
# were positive at least in 3 samples within treatment groups. 

#I will use these data frames 
View(results_relative_mat.mII)
dim(results_relative_mat.mII)
unique(results_relative_mat.mII$Assay)
3128/136 
#[1] 23
dim(results_relative_mat.m)
#[1] 3128    7
str(results_relative_mat.m)
results_relative_mat.mII$treat<-results_relative_mat.m$treat

View(Assay_annot) #ARG annotation
results_relative_mat.mII_annot<-merge(results_relative_mat.mII, Assay_annot,by = "Assay")
str(results_relative_mat.mII_annot)
# unused levels need to be removed
results_relative_mat.mII_annot$Mechanism<-droplevels(results_relative_mat.mII_annot$Mechanism)
results_relative_mat.mII_annot$Classification<-droplevels(results_relative_mat.mII_annot$Classification)
#Checking that everihting is ok:
min(results_relative_mat.mII_annot$value)
unique(results_relative_mat.mII_annot$pen_treat)
#I want to have only those assays that are detected in at least 3 samples / treatment group
ARGcounts<-count(subset(results_relative_mat.mII_annot, value>9.536743e-07), Assay, treat)
str(ARGcounts)
ARGcountsdf<-as.data.frame(ARGcounts)
rownames(ARGcountsdf)<-interaction(ARGcountsdf$Assay, ARGcountsdf$treat, sep="_")
#this should give the idea 
ARGcountsdf$n
#[1] 2 2 2 4 2 3 1 2 4 1 6 5 6 6 2 1 2 5 3 6 6 6 5 6 4 1 3 4 5 6 5 1 5 3 5 3 3 2 1 3 6 5 6 6 3 4 1 1 3 2 2 3 5 5 4 2 5 2 4 6 4 6 6 6 5 6 6 6 5 6 6 2 6 4 6
#[76] 6 1 1 3 3 4 3 2 5 4 5 4 6 5 6 6 2 2 4 2 3 4 1 1 1 1 1 5 1 6 5 6 6 1 1 3 3 2 2 2 1 6 4 4 4 1 2 1 1 1 1 2 2 5 1 3 6 5 6 6 1 6 5 6 6 6 5 6 6 5 5 6 3 6 5
#[151] 6 6 5 5 6 6 5 5 6 5 2 1 1 1 1 1 1 3 3 1 5 3 2 2 2 6 4 6 6 1 1 2 2 2 1 1 2 2 4 5 4 1 1 2 1 6 5 6 6 3 2 3 3 1 1 2 1 3 1 5 3 1 6 5 6 6 4 3 4 4 3 3 4 5 2
#[226] 1 3 2 5 2 3 4 3 3 1 1 4 1 2 1 1 1 2 4 2 3 1 1 3 6 4 6 4 1 1 1 1 1 1 2 1 6 5 6 6 2 6 5 6 6 2 4 1 4 3 1 1 2 5 4 6 5 1 3 1 2 5 2 5 3 1 6 5 6 6 2 1 6 5 6
#[301] 6 2 3 3 4 1 4 3 6 5 6 6 4 1 4 5 3 2 3 3 5 2 4 5 6 5 6 6 5 5 6 6 6 5 6 6 6 5 6 6 1 1 1 2 1 3 2 6 5 6 6 6 5 6 6 1 3 2 2 4 5 5 4 6 5 6 6 5 3 3 4 1 2 3 1
#[376] 1 1 6 5 5 6 3 1 4 2 3 1 5 2 2 1 1 6 5 6 6 6 5 6 6 1 2 1 3 3 2 2 1 5 3 3 6 5 6 6 2 3 6 3

#I want to check that I got what I wanted and not something else
subset(results_relative_mat.mII_annot$value, results_relative_mat.mII_annot$Assay=="aac(3)-iid_iii_iif_iia_iie_410" 
       & results_relative_mat.mII_annot$treat=="ZnCu")
#[1] 1.147449e-05 8.649883e-06 9.536743e-07 1.597981e-05 9.536743e-07 7.680827e-06
#This matches with the ARGcounts data frame
#So the goal is to get only those genes that are detected minimum 3x in a treatment group.
# --> the number needs to be bigger than 2.
#I will subset these assays / tratment group and save them as target vectors.
# and do every treatment group separately

# I will start with NTC group

NTC_target<-subset(ARGcountsdf$Assay, ARGcountsdf$n > 2 & ARGcountsdf$treat=="NTC")
unique(NTC_target)
length(NTC_target)
#[1] 68
#removing unused levels is always a good idea
NTC_target<-droplevels(NTC_target)
#First just NTC sampels
NTC_annot_dat<-subset(results_relative_mat.mII_annot, treat=="NTC")
dim(NTC_annot_dat)
#[1] 816   6
#Then only those assays that are in the target I saved
NTC_annot_dat<-NTC_annot_dat[NTC_annot_dat$Assay %in% levels(NTC_target),]
dim(NTC_annot_dat)
#[1] 408   6
#Check that levels are ok
NTC_annot_dat$Mechanism<-droplevels(NTC_annot_dat$Mechanism)
levels(NTC_annot_dat$Mechanism)
#1] "deactivate" "efflux"     "MGE"        "other"      "protection"
#Lets also check again the parent data frame:
levels(results_relative_mat.mII_annot$Mechanism)
which(results_relative_mat.mII_annot$Mechanism==  "other" )
results_relative_mat.mII_annot[2232:2254,]
#terW_1572 
which(results_relative_mat.mII_annot$Mechanism== "unknown")
results_relative_mat.mII_annot[2577:2599,]
#tetU_69... (Is not a resistance gene actually)
#But everything is ok for now.
#I need to have MGEs and ARGs in separate data frames and for every treatment group.
NTC_MGE_dat<-subset(NTC_annot_dat, Mechanism=="MGE")
NTC_ARG_dat<-subset(NTC_annot_dat, Mechanism!="MGE")
dim(NTC_MGE_dat) #66
dim(NTC_ARG_dat) #342
#Next I need to have these in matrix
#And only Assays and values
NTC_MGE<-NTC_MGE_dat[,-c(4,5,6)]
NTC_MGE<-dcast(NTC_MGE, pen_treat ~ Assay, value.var = "value")
rownames(NTC_MGE)<-NTC_MGE$pen_treat
NTC_MGE$pen_treat<-NULL

NTC_ARG<-NTC_ARG_dat[,-c(4,5,6)]
NTC_ARG<-dcast(NTC_ARG, pen_treat ~ Assay, value.var = "value")
rownames(NTC_ARG)<-NTC_ARG$pen_treat
NTC_ARG$pen_treat<-NULL

#Now I need to chance the 9.536743e-07 to NA's
#so that they can be omitted
#I dont want to create "artificial results"
NTC_MGE[NTC_MGE==9.536743e-07] <- NA
NTC_ARG[NTC_ARG==9.536743e-07] <- NA
#And I will change them to matrix for calculations:
NTC_ARG_mat<-as.matrix(NTC_ARG)
NTC_MGE_mat<-as.matrix(NTC_MGE)
#Because I need to have a correlation between two
#DIFFERENT matrices, I need to do the correlation analysis with 
#two R-packages.
#This is kind of a "dirty trick" and there is probably also another
#way to do this but this way I can also double check my results.

#For getting correlation coefficients and false discovery rate adjusted p-values I will use
#function corr.test from package psych:
ctNTC_ARG_MGEfdr<-corr.test(NTC_ARG_mat, NTC_MGE_mat, use="pairwise", method="spearman",
                            adjust="fdr",alpha=.05,ci=TRUE)
#Spearmans correlation is used because we might have 
#nonlinear dependencies. 

#Warning messages:
# 1: In sqrt(n - 2) : NaNs produced
#2: In pt(abs(t), (n - 2), log.p = TRUE) : NaNs produced
#3: In corr.test(NTC_ARG_mat, NTC_MGE_mat, use = "pairwise", method = "spearman",  :
#                 Number of subjects must be greater than 3 to find confidence intervals.
#              4: In sqrt(n - 3) : NaNs produced

#I get warnings because I have NA's. No need to worry in this case.
str(ctNTC_ARG_MGEfdr)
#I will put the results in a data frame
ctNTC_ARG_MGEfdr_ci<-ctNTC_ARG_MGEfdr$ci
View(ctNTC_ARG_MGEfdr_ci)
#corr.test combines the names and shortens them as it returns the results 
# in a "melted format". I will get the Assay names through the function "cor":
corNTC_ARG_MGE<-cor(NTC_ARG_mat, NTC_MGE_mat,use="pairwise.complete.obs", method="spearman")
dim(corNTC_ARG_MGE)
corNTC_ARG_MGE.m<-melt(corNTC_ARG_MGE)
View(corNTC_ARG_MGE.m)
#Now I have ARG, MGE and their correlation in the data frame. But this is missing the 
#p-values, which I can get from corr.test. 
dim(corNTC_ARG_MGE.m)
dim(ctNTC_ARG_MGEfdr_ci)
#And they have the same length. 
#I will add r and p columns to the df that has assay names:
corNTC_ARG_MGE.m$r<-ctNTC_ARG_MGEfdr_ci$r
corNTC_ARG_MGE.m$p.adj<-ctNTC_ARG_MGEfdr_ci$p
#I will compare the r's and values to make sure they match
View(corNTC_ARG_MGE.m)
#They do so I can remove the column named value
corNTC_ARG_MGE.m$value<-NULL
#Now I want to have only positive correletions >0.8 and p.adj <0.05:
corNTC_ARG_MGE.m_sig<-filter(corNTC_ARG_MGE.m, p.adj<= 0.05)
corNTC_ARG_MGE.m_sig<-filter(corNTC_ARG_MGE.m_sig, r>= 0.8)
#before saving the table, it is possible to get rid of the primer number if one wants to:
#corNTC_ARG_MGE.m_sigF<-corNTC_ARG_MGE.m_sig
#corNTC_ARG_MGE.m_sigF$Var1<-gsub("_.*","",corNTC_ARG_MGE.m_sigF$Var1) 
#corNTC_ARG_MGE.m_sigF$Var2<-gsub("_.*","",corNTC_ARG_MGE.m_sigF$Var2)
#Now these files can be used for getting the edges and nodes for network vizualization with Gephi.
#But I don't want to do that since some assays have the same name and numbers are useful because of that

#Edges are easier:
NTCedges<-corNTC_ARG_MGE.m_sig[,1:3]
write.table(NTCedges,file="NTCedges.txt",quote=FALSE, row.names=FALSE, col.names=FALSE,sep=",")

#Nodes:
NTCnodesI<-as.data.frame(unique(NTCedges$Var1)) #ARGs first
colnames(NTCnodesI)[1]<-"node"
NTCnodesI$type<-(rep("ARG",nrow(NTCnodesI)))
#I need to get their annotation somehow. 
str(NTC_ARG_dat)
#I have it there.
#so I need a target vector with assays in the right order:
NTC_assay_target<-unique(NTCnodesI$node)
#and I can get the "Classification" with the match function:
NTC_annot_dat[match(NTC_assay_target,NTC_annot_dat$Assay),]$Classification
NTCnodesI$class<-NTC_annot_dat[match(NTC_assay_target,NTC_annot_dat$Assay),]$Classification
#The same procedure for MGEs:
NTCnodesII<-as.data.frame(unique(NTCedges$Var2)) 
colnames(NTCnodesII)[1]<-"node"
NTCnodesII$type<-(rep("MGE",nrow(NTCnodesII)))
NTC_assay_targetII<-unique(NTCnodesII$node)
NTCnodesII$class<-NTC_annot_dat[match(NTC_assay_targetII,NTC_annot_dat$Assay),]$Classification
#Now I will combine ARG and MGE nodes:
NTCnodes<-rbind(NTCnodesI,NTCnodesII)
#And write the node table for Gephi
write.table(NTCnodes,file="NTCnodes.txt",quote=FALSE, row.names=FALSE, col.names=FALSE,sep=",")

#Now the same for AB group:

AB_target<-subset(ARGcountsdf$Assay, ARGcountsdf$n > 2 & ARGcountsdf$treat=="AB")
AB_target<-droplevels(AB_target)
unique(AB_target)
length(AB_target)
#[1] 71

AB_annot_dat<-subset(results_relative_mat.mII_annot, treat=="AB")
dim(AB_annot_dat)
#[1] 816   6
AB_annot_dat<-AB_annot_dat[AB_annot_dat$Assay %in% levels(AB_target),]
dim(AB_annot_dat)
#[1] 426   6
AB_annot_dat$Mechanism<-droplevels(AB_annot_dat$Mechanism)
levels(AB_annot_dat$Mechanism)
#1] "deactivate" "efflux"     "MGE"        "protection"
AB_MGE_dat<-subset(AB_annot_dat, Mechanism=="MGE")
AB_ARG_dat<-subset(AB_annot_dat, Mechanism!="MGE")
dim(AB_MGE_dat)
#72 6
dim(AB_ARG_dat)
#354 6
#Next I need to have these in matrix
#And only Assays and values
AB_MGE<-AB_MGE_dat[,-c(4,5,6,7)]
AB_MGE<-dcast(AB_MGE, pen_treat ~ Assay, value.var = "value")
rownames(AB_MGE)<-AB_MGE$pen_treat
AB_MGE$pen_treat<-NULL

AB_ARG<-AB_ARG_dat[,-c(4,5,6,7)]
AB_ARG<-dcast(AB_ARG, pen_treat ~ Assay, value.var = "value")
rownames(AB_ARG)<-AB_ARG$pen_treat
AB_ARG$pen_treat<-NULL

#Now I need to chance the 9.536743e-07 to NA's
#so that they can be omitted
AB_MGE[AB_MGE==9.536743e-07] <- NA
AB_ARG[AB_ARG==9.536743e-07] <- NA
#And I will change them to matrix for calculations:
AB_ARG_mat<-as.matrix(AB_ARG)
AB_MGE_mat<-as.matrix(AB_MGE)
#For getting correlation coefficients and false discovery rate adjusted p-values I will use
#function corr.test from package psych:
ctAB_ARG_MGEfdr<-corr.test(AB_ARG_mat, AB_MGE_mat, use="pairwise", method="spearman",
                           adjust="fdr",alpha=.05,ci=TRUE)
#Spearmans correlation is used because we might have 
#nonlinear dependencies. 

#Warning messages:
# 1: In sqrt(n - 2) : NaNs produced
#2: In pt(abs(t), (n - 2), log.p = TRUE) : NaNs produced
#3: In corr.test(NTC_ARG_mat, NTC_MGE_mat, use = "pairwise", method = "spearman",  :
#                 Number of subjects must be greater than 3 to find confidence intervals.
#              4: In sqrt(n - 3) : NaNs produced

#I get warnings because I have NA's. No need to worry in this case.
#(It is expected to have lot of NA's and NaN's)
str(ctAB_ARG_MGEfdr)
#I will put the results in a data frame
ctAB_ARG_MGEfdr_ci<-ctAB_ARG_MGEfdr$ci
View(ctAB_ARG_MGEfdr_ci)
#corr.test combines the names and shortens them as it returns the results
# in a "melted format". I will get the Assay names through the function "cor":
corAB_ARG_MGE<-cor(AB_ARG_mat, AB_MGE_mat,use="pairwise.complete.obs", method="spearman")
dim(corAB_ARG_MGE)
corAB_ARG_MGE.m<-melt(corAB_ARG_MGE)
View(corAB_ARG_MGE.m)
#Now I have ARG, MGE and their correlation in the data frame. But this is missing the 
#p-values, which I can get from corr.test. 
dim(corAB_ARG_MGE.m)
dim(ctAB_ARG_MGEfdr_ci)
#And they have the same dimensions. 
#I will add r and p columns to to df that has assay names:
corAB_ARG_MGE.m$r<-ctAB_ARG_MGEfdr_ci$r
corAB_ARG_MGE.m$p.adj<-ctAB_ARG_MGEfdr_ci$p
#I will compare the r's and values to make sure they match
View(corAB_ARG_MGE.m)
#They do so I can remove the column named value
corAB_ARG_MGE.m$value<-NULL
#Now I want to have only positive correletions >0.8 and p.adj <0.05:
corAB_ARG_MGE.m_sig<-filter(corAB_ARG_MGE.m, p.adj<= 0.05)
corAB_ARG_MGE.m_sig<-filter(corAB_ARG_MGE.m_sig, r>= 0.8)
#before saving the table, it is possible to get rid of the primer number if one wants to:
#corNTC_ARG_MGE.m_sigF<-corNTC_ARG_MGE.m_sig
#corNTC_ARG_MGE.m_sigF$Var1<-gsub("_.*","",corNTC_ARG_MGE.m_sigF$Var1) 
#corNTC_ARG_MGE.m_sigF$Var2<-gsub("_.*","",corNTC_ARG_MGE.m_sigF$Var2)
#Now these files can be used for getting the edges and nodes for network vizualization with Gephi.
#But I don't want to do that since some assays have the same name and numbers are useful because of that

#Edges:
ABedges<-corAB_ARG_MGE.m_sig[,1:3]
write.table(ABedges,file="ABedges.txt",quote=FALSE, row.names=FALSE, col.names=FALSE,sep=",")

#Nodes:
ABnodesI<-as.data.frame(unique(ABedges$Var1)) #ARGs first
colnames(ABnodesI)[1]<-"node"
ABnodesI$type<-(rep("ARG",nrow(ABnodesI)))
#I need to get their annotation 
#so I need a target vector with assays in the right order:
AB_assay_target<-unique(ABnodesI$node)
#and I can get the "Classification" with the match function:
AB_annot_dat[match(AB_assay_target, AB_annot_dat$Assay),]$Classification
ABnodesI$class<-AB_annot_dat[match(AB_assay_target, AB_annot_dat$Assay),]$Classification
#The same procedure for MGEs:
ABnodesII<-as.data.frame(unique(ABedges$Var2)) 
colnames(ABnodesII)[1]<-"node"
ABnodesII$type<-(rep("MGE",nrow(ABnodesII)))
AB_assay_targetII<-unique(ABnodesII$node)
ABnodesII$class<-AB_annot_dat[match(AB_assay_targetII, AB_annot_dat$Assay),]$Classification
#Now I will combine ARG and MGE nodes:
ABnodes<-rbind(ABnodesI,ABnodesII)
#And write the node table for Gephi
write.table(ABnodes,file="ABnodes.txt",quote=FALSE, row.names=FALSE, col.names=FALSE,sep=",")

#For Mushroom:

M_target<-subset(ARGcountsdf$Assay, ARGcountsdf$n > 2 & ARGcountsdf$treat=="M")
M_target<-droplevels(M_target)
unique(M_target)
length(M_target)
#[1] 62

M_annot_dat<-subset(results_relative_mat.mII_annot, treat=="M")
dim(M_annot_dat)
#[1] 680   6
M_annot_dat<-M_annot_dat[M_annot_dat$Assay %in% levels(M_target),]
dim(M_annot_dat)
#[1] 310   6
M_annot_dat$Mechanism<-droplevels(M_annot_dat$Mechanism)
levels(M_annot_dat$Mechanism)
#[1] "deactivate" "efflux"     "MGE"        "other"      "protection" "unknown"  

M_MGE_dat<-subset(M_annot_dat, Mechanism=="MGE")
M_ARG_dat<-subset(M_annot_dat, Mechanism!="MGE")
dim(M_MGE_dat)
#60 6
dim(M_ARG_dat)
#250 6
#Next I need to have these in matrix
#And only Assays and values
M_MGE<-M_MGE_dat[,-c(4,5,6,7)]
M_MGE<-dcast(M_MGE, pen_treat ~ Assay, value.var = "value")
rownames(M_MGE)<-M_MGE$pen_treat
M_MGE$pen_treat<-NULL

M_ARG<-M_ARG_dat[,-c(4,5,6,7)]
M_ARG<-dcast(M_ARG, pen_treat ~ Assay, value.var = "value")
rownames(M_ARG)<-M_ARG$pen_treat
M_ARG$pen_treat<-NULL

#Now I need to chance the 9.536743e-07 to NA's
#so that they can be omitted
M_MGE[M_MGE==9.536743e-07] <- NA
M_ARG[M_ARG==9.536743e-07] <- NA
#And I will change them to matrix for calculations:
M_ARG_mat<-as.matrix(M_ARG)
M_MGE_mat<-as.matrix(M_MGE)
#Because I need to have a correlation between two
#DIFFERENT matrices, I need to do the correlation analysis with 
#two R-packages.
#This is a kind of "dirty trick" and there is probably also another
#way to do this but this way I can also double check my results.

#For getting correlation coefficients and false discovery rate adjusted p-values I will use
#function corr.test from package psych:
ctM_ARG_MGEfdr<-corr.test(M_ARG_mat, M_MGE_mat, use="pairwise", method="spearman",
                          adjust="fdr",alpha=.05,ci=TRUE)
#Spearmans correlation is used because we might have 
#nonlinear dependencies. 

#Warning messages:
# 1: In sqrt(n - 2) : NaNs produced
#2: In pt(abs(t), (n - 2), log.p = TRUE) : NaNs produced
#3: In corr.test(NTC_ARG_mat, NTC_MGE_mat, use = "pairwise", method = "spearman",  :
#                 Number of subjects must be greater than 3 to find confidence intervals.
#              4: In sqrt(n - 3) : NaNs produced

#I get warnings because I have NA's. No need to worry in this case.
#(It is expected to have lot of NA's and NaN's)
str(ctM_ARG_MGEfdr)
#I will put the results in a data frame
ctM_ARG_MGEfdr_ci<-ctM_ARG_MGEfdr$ci
View(ctM_ARG_MGEfdr_ci)
#corr.test combines the names and shortens them as it returns the results
# in a "melted format". I will get the Assay names through the function "cor":
corM_ARG_MGE<-cor(M_ARG_mat, M_MGE_mat,use="pairwise.complete.obs", method="spearman")
dim(corM_ARG_MGE)
corM_ARG_MGE.m<-melt(corM_ARG_MGE)
View(corM_ARG_MGE.m)
#Now I have ARG, MGE and their correlation in the data frame. But this is missing the 
#p-values, which I can get from corr.test. 
dim(corM_ARG_MGE.m)
dim(ctM_ARG_MGEfdr_ci)
#And they have the same dimensions. 
#I will add r and p columns to to df that has assay names:
corM_ARG_MGE.m$r<-ctM_ARG_MGEfdr_ci$r
corM_ARG_MGE.m$p.adj<-ctM_ARG_MGEfdr_ci$p
#I will compare the r's and values to make sure they match
View(corM_ARG_MGE.m)
#They do so I can remove the column named value
corM_ARG_MGE.m$value<-NULL
#Now I want to have only positive correletions >0.8 and p.adj <0.05:
corM_ARG_MGE.m_sig<-filter(corM_ARG_MGE.m, p.adj<= 0.05)
corM_ARG_MGE.m_sig<-filter(corM_ARG_MGE.m_sig, r>= 0.8)
#before saving the table, it is possible to get rid of the primer number if one wants to:
#corNTC_ARG_MGE.m_sigF<-corNTC_ARG_MGE.m_sig
#corNTC_ARG_MGE.m_sigF$Var1<-gsub("_.*","",corNTC_ARG_MGE.m_sigF$Var1) 
#corNTC_ARG_MGE.m_sigF$Var2<-gsub("_.*","",corNTC_ARG_MGE.m_sigF$Var2)
#Now these files can be used for getting the edges and nodes for network vizualization with Gephi.
#But I don't want to do that since some assays have the same name and numbers are useful because of that

#Edges:
Medges<-corM_ARG_MGE.m_sig[,1:3]
write.table(Medges,file="Medges.txt",quote=FALSE, row.names=FALSE, col.names=FALSE,sep=",")

#Nodes:
MnodesI<-as.data.frame(unique(Medges$Var1)) #ARGs first
colnames(MnodesI)[1]<-"node"
MnodesI$type<-(rep("ARG",nrow(MnodesI)))
#I need to get their annotation 
#so I need a target vector with assays in the right order:
M_assay_target<-unique(MnodesI$node)
#and I can get the "Classification" with the match function:
M_annot_dat[match(M_assay_target, M_annot_dat$Assay),]$Classification
MnodesI$class<-M_annot_dat[match(M_assay_target, M_annot_dat$Assay),]$Classification
#The same procedure for MGEs:
MnodesII<-as.data.frame(unique(Medges$Var2)) 
colnames(MnodesII)[1]<-"node"
MnodesII$type<-(rep("MGE",nrow(MnodesII)))
M_assay_targetII<-unique(MnodesII$node)
MnodesII$class<-M_annot_dat[match(M_assay_targetII, M_annot_dat$Assay),]$Classification
#Now I will combine ARG and MGE nodes:
Mnodes<-rbind(MnodesI,MnodesII)
#And write the node table for Gephi
write.table(Mnodes,file="Mnodes.txt",quote=FALSE, row.names=FALSE, col.names=FALSE,sep=",")


#For ZnCu:

ZnCu_target<-subset(ARGcountsdf$Assay, ARGcountsdf$n > 2 & ARGcountsdf$treat=="ZnCu")
ZnCu_target<-droplevels(ZnCu_target)
unique(ZnCu_target)
length(ZnCu_target)
#[1] 69

ZnCu_annot_dat<-subset(results_relative_mat.mII_annot, treat=="ZnCu")
dim(ZnCu_annot_dat)
#[1] 816   7
ZnCu_annot_dat<-ZnCu_annot_dat[ZnCu_annot_dat$Assay %in% levels(ZnCu_target),]
dim(ZnCu_annot_dat)
#[1] 414   7
ZnCu_annot_dat$Mechanism<-droplevels(ZnCu_annot_dat$Mechanism)
levels(ZnCu_annot_dat$Mechanism)
#[1] "deactivate" "efflux"     "MGE"        "protection"

ZnCu_MGE_dat<-subset(ZnCu_annot_dat, Mechanism=="MGE")
ZnCu_ARG_dat<-subset(ZnCu_annot_dat, Mechanism!="MGE")
dim(ZnCu_MGE_dat)
#72  7
dim(ZnCu_ARG_dat)
# 342   7
#Next I need to have these in matrix
#And only Assays and values
ZnCu_MGE<-ZnCu_MGE_dat[,-c(4,5,6,7)]
ZnCu_MGE<-dcast(ZnCu_MGE, pen_treat ~ Assay, value.var = "value")
rownames(ZnCu_MGE)<-ZnCu_MGE$pen_treat
ZnCu_MGE$pen_treat<-NULL

ZnCu_ARG<-ZnCu_ARG_dat[,-c(4,5,6,7)]
ZnCu_ARG<-dcast(ZnCu_ARG, pen_treat ~ Assay, value.var = "value")
rownames(ZnCu_ARG)<-ZnCu_ARG$pen_treat
ZnCu_ARG$pen_treat<-NULL

#Now I need to chance the 9.536743e-07 to NA's
#so that they can be omitted
ZnCu_MGE[ZnCu_MGE==9.536743e-07] <- NA
ZnCu_ARG[ZnCu_ARG==9.536743e-07] <- NA
#And I will change them to matrix for calculations:
ZnCu_ARG_mat<-as.matrix(ZnCu_ARG)
ZnCu_MGE_mat<-as.matrix(ZnCu_MGE)
#Because I need to have a correlation between two
#DIFFERENT matrices, I need to do the correlation analysis with 
#two R-packages.
#This is a kind of "dirty trick" and there is probably also another
#way to do this but this way I can also double check my results.

#For getting correlation coefficients and false discovery rate adjusted p-values I will use
#function corr.test from package psych:
ctZnCu_ARG_MGEfdr<-corr.test(ZnCu_ARG_mat, ZnCu_MGE_mat, use="pairwise", method="spearman",
                             adjust="fdr",alpha=.05,ci=TRUE)
#Spearmans correlation is used because we might have 
#nonlinear dependencies. 

#Warning messages:
# 1: In sqrt(n - 2) : NaNs produced
#2: In pt(abs(t), (n - 2), log.p = TRUE) : NaNs produced
#3: In corr.test(NTC_ARG_mat, NTC_MGE_mat, use = "pairwise", method = "spearman",  :
#                 Number of subjects must be greater than 3 to find confidence intervals.
#              4: In sqrt(n - 3) : NaNs produced

#So I get warnings because I have NA's. No need to worry in this case.
#(It is expected to have lot of NA's and NaN's)
str(ctZnCu_ARG_MGEfdr)
#I will put the results in a data frame
ctZnCu_ARG_MGEfdr_ci<-ctZnCu_ARG_MGEfdr$ci
View(ctZnCu_ARG_MGEfdr_ci)
#corr.test combines the names and shortens them as it returns the results 
# in a "melted format". I will get the Assay names through the function "cor":
corZnCu_ARG_MGE<-cor(ZnCu_ARG_mat, ZnCu_MGE_mat,use="pairwise.complete.obs", method="spearman")
dim(corZnCu_ARG_MGE)
corZnCu_ARG_MGE.m<-melt(corZnCu_ARG_MGE)
View(corZnCu_ARG_MGE.m)
#Now I have ARG, MGE and their correlation in the data frame. But this is missing the 
#p-values, which I can get from corr.test. 
dim(corZnCu_ARG_MGE.m)
dim(ctZnCu_ARG_MGEfdr_ci)
#And they have the same dimensions. 
#I will add r and p columns to to df that has assay names:
corZnCu_ARG_MGE.m$r<-ctZnCu_ARG_MGEfdr_ci$r
corZnCu_ARG_MGE.m$p.adj<-ctZnCu_ARG_MGEfdr_ci$p
#I will compare the r's and values to make sure they match
View(corZnCu_ARG_MGE.m)
#They do so I can remove the column named value
corZnCu_ARG_MGE.m$value<-NULL
#Now I want to have only positive correletions >0.8 and p.adj <0.05:
corZnCu_ARG_MGE.m_sig<-filter(corZnCu_ARG_MGE.m, p.adj<= 0.05)
corZnCu_ARG_MGE.m_sig<-filter(corZnCu_ARG_MGE.m_sig, r>= 0.8)
#before saving the table, it is possible to get rid of the primer number if one wants to:
#corNTC_ARG_MGE.m_sigF<-corNTC_ARG_MGE.m_sig
#corNTC_ARG_MGE.m_sigF$Var1<-gsub("_.*","",corNTC_ARG_MGE.m_sigF$Var1) 
#corNTC_ARG_MGE.m_sigF$Var2<-gsub("_.*","",corNTC_ARG_MGE.m_sigF$Var2)
#Now these files can be used for getting the edges and nodes for network vizualization with Gephi.
#But I don't want to do that since some assays have the same name and numbers are useful because of that

#Edges:
ZnCuedges<-corZnCu_ARG_MGE.m_sig[,1:3]
write.table(ZnCuedges,file="ZnCuedges.txt",quote=FALSE, row.names=FALSE, col.names=FALSE,sep=",")

#Nodes:
ZnCunodesI<-as.data.frame(unique(ZnCuedges$Var1)) #ARGs first
colnames(ZnCunodesI)[1]<-"node"
ZnCunodesI$type<-(rep("ARG",nrow(ZnCunodesI)))
#I need to get their annotation 
#so I need a target vector with assays in the right order:
ZnCu_assay_target<-unique(ZnCunodesI$node)
#and I can get the "Classification" with the match function:
ZnCu_annot_dat[match(ZnCu_assay_target, ZnCu_annot_dat$Assay),]$Classification
ZnCunodesI$class<-ZnCu_annot_dat[match(ZnCu_assay_target, ZnCu_annot_dat$Assay),]$Classification
#The same procedure for MGEs:
ZnCunodesII<-as.data.frame(unique(ZnCuedges$Var2)) 
colnames(ZnCunodesII)[1]<-"node"
ZnCunodesII$type<-(rep("MGE",nrow(ZnCunodesII)))
ZnCu_assay_targetII<-unique(ZnCunodesII$node)
ZnCunodesII$class<-ZnCu_annot_dat[match(ZnCu_assay_targetII, ZnCu_annot_dat$Assay),]$Classification
#Now I will combine ARG and MGE nodes:
ZnCunodes<-rbind(ZnCunodesI,ZnCunodesII)
#And write the node table for Gephi
write.table(ZnCunodes,file="ZnCunodes.txt",quote=FALSE, row.names=FALSE, col.names=FALSE,sep=",")

########


############################################################################################# 
#                                                                                           #
#          How treatments influenced the microbiome and resistome: machine learning         #
#                                                                                           #
#############################################################################################



library(Rtsne)
library(dbscan)
library(ggrepel)
library(caret)
library(edarf)
library(ranger)

#t-SNE 
#Useful links:
# https://www.machinegurning.com/rstats/tsne/
# https://distill.pub/2016/misread-tsne/

# I want to combine ARGs & MGEs with TSS OTUs; they are both relative abundance

dim(noMC_phylOTUtable_RA_NoZero_matPRO)
str(noMC_phylOTUtable_RA_NoZero_matPRO)
dim(results_relative_mat_PRO)
str(results_relative_mat_PRO)
dim(meta_forPRO)

OTU_RAandARG_tSNE_mat<-cbind(noMC_phylOTUtable_RA_NoZero_matPRO,results_relative_mat_PRO)

rownames(OTU_RAandARG_tSNE_mat)<-rownames(results_relative_mat_PRO)

# using t-SNE for clustering is based on visualization, it's trial and error,
# althogh some rules apply.
#Fitting
set.seed(42)
ARG_OTUtsne500_five <- Rtsne(OTU_RAandARG_tSNE_mat, dims = 2,
                             perplexity = 5,
                             verbose = TRUE, max_iter = 500)
#ploting:
plot(ARG_OTUtsne500_five$Y, type="n",
     ann = FALSE, xaxt = "n")
title(main="OTUs & ARGs perp5 it 500")
with(meta_forPRO,text(ARG_OTUtsne500_five$Y, labels = pen_treat,
                      col = meta_forPRO$color,
                      cex = 1.5))
# I will increase the iterations
set.seed(43)
ARG_OTUtsne5000_five <- Rtsne(OTU_RAandARG_tSNE_mat, dims = 2,
                              perplexity = 5,
                              verbose = TRUE, max_iter = 5000)
plot(ARG_OTUtsne5000_five$Y, type="n",
     ann = FALSE, xaxt = "n")
title(main="OTUs & ARGs perp5 it 5000")
with(meta_forPRO,text(ARG_OTUtsne5000_five$Y, labels = pen_treat,
                      col = meta_forPRO$color,
                      cex = 1.5))

# and try with 1000 iterations
set.seed(44)
ARG_OTUtsne1000_five <- Rtsne(OTU_RAandARG_tSNE_mat, dims = 2,
                              perplexity = 5,
                              verbose = TRUE, max_iter = 1000)

plot(ARG_OTUtsne1000_five$Y, type="n",
     ann = FALSE, xaxt = "n")
title(main="OTUs & ARGs perp5 it 1000")
with(meta_forPRO,text(ARG_OTUtsne1000_five$Y, labels = pen_treat,
                      col = meta_forPRO$color,
                      cex = 1.5))
# 1500
set.seed(44)
ARG_OTUtsne1500_five <- Rtsne(OTU_RAandARG_tSNE_mat, dims = 2,
                              perplexity = 5,
                              verbose = TRUE, max_iter = 1500)


plot(ARG_OTUtsne1500_five$Y, type="n",
     ann = FALSE, xaxt = "n")
title(main="OTUs & ARGs perp5 it 1500")
with(meta_forPRO,text(ARG_OTUtsne1500_five$Y, labels = pen_treat,
                      col = meta_forPRO$color,
                      cex = 1.5))

# What I want to see are "islands". So maybe here?

#2000
ARG_OTUtsne2000_five <- Rtsne(OTU_RAandARG_tSNE_mat, dims = 2,
                              perplexity = 5,
                              verbose = TRUE, max_iter = 2000)

plot(ARG_OTUtsne2000_five$Y, type="n",
     ann = FALSE, xaxt = "n")
title(main="OTUs & ARGs perp5 it 2000")
with(meta_forPRO,text(ARG_OTUtsne2000_five$Y, labels = pen_treat,
                      col = meta_forPRO$color,
                      cex = 1.5))

#also maybe?

# 2500
ARG_OTUtsne2500_five <- Rtsne(OTU_RAandARG_tSNE_mat, dims = 2,
                              perplexity = 5,
                              verbose = TRUE, max_iter = 2500)

plot(ARG_OTUtsne2500_five$Y, type="n",
     ann = FALSE, xaxt = "n")
title(main="OTUs & ARGs perp5 it 2500")
with(meta_forPRO,text(ARG_OTUtsne2500_five$Y, labels = pen_treat,
                      col = meta_forPRO$color,
                      cex = 1.5))

# this looks weird.
set.seed(47)
ARG_OTUtsne3000_five <- Rtsne(OTU_RAandARG_tSNE_mat, dims = 2,
                              perplexity = 5,
                              verbose = TRUE, max_iter = 3000)


plot(ARG_OTUtsne3000_five$Y, type="n",
     ann = FALSE, xaxt = "n")
title(main="OTUs & ARGs perp5 it 3000")
with(meta_forPRO,text(ARG_OTUtsne3000_five$Y, labels = pen_treat,
                      col = meta_forPRO$color,
                      cex = 1.5))

#no. 

set.seed(48)
ARG_OTUtsne3500_five <- Rtsne(OTU_RAandARG_tSNE_mat, dims = 2,
                              perplexity = 5,
                              verbose = TRUE, max_iter = 3500)


plot(ARG_OTUtsne3500_five$Y, type="n",
     ann = FALSE, xaxt = "n")
title(main="OTUs & ARGs perp5 it 3500")
with(meta_forPRO,text(ARG_OTUtsne3500_five$Y, labels = pen_treat,
                      col = meta_forPRO$color,
                      cex = 1.5))

#hmm. smallest error so far, but the plot doesnt look nice

set.seed(49)
ARG_OTUtsne4000_five <- Rtsne(OTU_RAandARG_tSNE_mat, dims = 2,
                              perplexity = 5,
                              verbose = TRUE, max_iter = 4000)


plot(ARG_OTUtsne4000_five$Y, type="n",
     ann = FALSE, xaxt = "n")
title(main="OTUs & ARGs perp5 it 4000")
with(meta_forPRO,text(ARG_OTUtsne4000_five$Y, labels = pen_treat,
                      col = meta_forPRO$color,
                      cex = 1.5))

#slightly bigger error, plot is "more evenly distributed"

set.seed(50)
ARG_OTUtsne4500_five <- Rtsne(OTU_RAandARG_tSNE_mat, dims = 2,
                              perplexity = 5,
                              verbose = TRUE, max_iter = 4500)


plot(ARG_OTUtsne4500_five$Y, type="n",
     ann = FALSE, xaxt = "n")
title(main="OTUs & ARGs perp5 it 4500")
with(meta_forPRO,text(ARG_OTUtsne4500_five$Y, labels = pen_treat,
                      col = meta_forPRO$color,
                      cex = 1.5))

#error gets bigger

set.seed(51)
ARG_OTUtsne10000_five <- Rtsne(OTU_RAandARG_tSNE_mat, dims = 2,
                               perplexity = 5,
                               verbose = TRUE, max_iter = 10000)


plot(ARG_OTUtsne10000_five$Y, type="n",
     ann = FALSE, xaxt = "n")
title(main="OTUs & ARGs perp5 it 10000")
with(meta_forPRO,text(ARG_OTUtsne10000_five$Y, labels = pen_treat,
                      col = meta_forPRO$color,
                      cex = 1.5))

# This is quite ok! But so is 5 & 1500

set.seed(52)
ARG_OTUtsne50000_five <- Rtsne(OTU_RAandARG_tSNE_mat, dims = 2,
                               perplexity = 5,
                               verbose = TRUE, max_iter = 50000)


plot(ARG_OTUtsne50000_five$Y, type="n",
     ann = FALSE, xaxt = "n")
title(main="OTUs & ARGs perp5 it 50000")
with(meta_forPRO,text(ARG_OTUtsne50000_five$Y, labels = pen_treat,
                      col = meta_forPRO$color,
                      cex = 1.5))

#this is my choice!

# next step is to identyfy the clusters

library(dbscan)

ARG_OTUtsne50000_five_mat<- as.matrix(ARG_OTUtsne50000_five$Y)
colnames(ARG_OTUtsne50000_five_mat)<-c("X1","X2")
rownames(ARG_OTUtsne50000_five_mat)<-rownames(meta_forPRO)


cl50000five <- hdbscan(ARG_OTUtsne50000_five_mat, minPts = 3)
cl50000five

#HDBSCAN clustering for 21 objects.
#Parameters: minPts = 3
#The clustering contains 3 cluster(s) and 0 noise points.
#1  2  3 
#4 14  3 
#Available fields: cluster, minPts, cluster_scores, membership_prob, outlier_scores, hc
plot(ARG_OTUtsne50000_five_mat, col=cl50000five$cluster+1, pch=cl50000five$cluster)
# I have 3 clusters
plot(cl50000five)
cl50000five$hc
#Call:
 # hdbscan(x = ARG_OTUtsne50000_five_mat, minPts = 3)
#Cluster method   : robust single 
#Distance         : mutual reachability 
#Number of objects: 21 
#plot(cl1500six$hc, main="HDBSCAN* Hierarchy")

print(cl50000five$cluster_scores)

#1          2          3 
#0.2215518 0.8908191 0.1339982 

head(cl50000five$membership_prob)

plot(ARG_OTUtsne50000_five_mat, col=cl50000five$cluster, pch=21)
colors <- mapply(function(col, i) adjustcolor(col, alpha.f = cl50000five$membership_prob[i]), 
                 palette()[cl50000five$cluster], seq_along(cl50000five$cluster))
points(ARG_OTUtsne50000_five_mat, col=colors, pch=20)


top_outliers <- order(cl50000five$outlier_scores, decreasing = T)[1:10]
colors <- mapply(function(col, i) adjustcolor(col, alpha.f = cl50000five$outlier_scores[i]), 
                 palette()[cl50000five$cluster+1], seq_along(cl50000five$cluster))
plot(ARG_OTUtsne50000_five_mat, col=colors, pch=20)
text(ARG_OTUtsne50000_five_mat[top_outliers, ], labels = top_outliers,pos=4, cex=0.5)

top_outliers

# So the result is not super reliable, but at least the two smaller clusters seem to be 
# actual clusters. 

# I want to make a nice plot with ggplot.

df.ARG_OTUtsne50000_five<-as.data.frame(ARG_OTUtsne50000_five$Y)

rownames(df.ARG_OTUtsne50000_five)<-rownames(meta_forPRO)
df.ARG_OTUtsne50000_five$treatment<-meta_forPRO$treatment
df.ARG_OTUtsne50000_five$color<-meta_forPRO$color
df.ARG_OTUtsne50000_five$cluster<-cl50000five$cluster
df.ARG_OTUtsne50000_five$pen_treat<-meta_forPRO$pen_treat

library(ggrepel)
#library(colorspace)
cl_col<-sequential_hcl(5, palette = "Lajolla")

colnames(df.ARG_OTUtsne50000_five)[1]<-"Dimension1"
colnames(df.ARG_OTUtsne50000_five)[2]<-"Dimension2"

tSNE<-ggplot(df.ARG_OTUtsne50000_five, aes(Dimension1, Dimension2,label=pen_treat)) +
  geom_point(aes(color=factor(cluster)), size=4)+
  geom_text_repel(color=df.ARG_OTUtsne50000_five$color, size=4, fontface="bold",
                  aes(Dimension1, Dimension2, label =pen_treat))+
  #geom_text(aes(color=factor(cluster))) +
  scale_color_manual(name = "Cluster", values = rev(cl_col)) +
  theme_bw(base_size = 16)

#ggplot(df.ARG_OTUtsne10000_five) +
# geom_point(aes(V1, V2, color = Experiment), size = 1) +
#geom_label_repel(aes(V1, V2, color=color, label = df.ARG_OTUtsne10000_five$x),size=1.9, fontface = 'bold',fill="white",
#                box.padding = unit(0.01,"lines"))+
#geom_text_repel(
# aes(V1, V2, color=Experiment, label =cluster),
#size=4,fontface = 'bold', segment.size = 2, segment.alpha=0.15,
#box.padding = unit(0.8,"lines"), point.padding = unit(0.5,"lines")) +
#scale_color_manual(name = "Experiment", values = unique(df.ARG_OTUtsne10000_five$color)) +
#theme_classic(base_size = 14)

# Now the question is what genes or OTUs determine the clustering.
# Random forest answers.

# The script below is modified from: https://github.com/Begia/Hazen16S

library(caret)
library(ranger)
library(edarf)
install.packages('e1071', dependencies=TRUE)
library(e1071)

models <- list(NULL)
models_2 <- list(NULL)
pd_data <- list(NULL)
pd_data_2 <- list(NULL)
plot_pd_data <- list(NULL)
plot_pd_data_2 <- list(NULL)
tax_level_nums <- NA
full_model_errors <- list(NULL)
model_errors_2 <- list(NULL)
all_full_model_errors <- list(NULL)
all_model_errors_2 <- list(NULL)

response <- factor(df.ARG_OTUtsne50000_five$cluster)

rf_data <- data.frame(response, OTU_RAandARG_tSNE_mat)

# Now I will run a random forest model
set.seed(42)
classify <- ranger(response ~ ., data = rf_data, num.trees=10000, importance="impurity")
classify 
#Call:
#ranger(response ~ ., data = rf_data, num.trees = 10000, importance = "impurity") 

#Type:                             Classification 
#Number of trees:                  5000 
#Sample size:                      21 
#Number of independent variables:  268 
#Mtry:                             16 
#Target node size:                 1 
#Variable importance mode:         impurity 
#Splitrule:                        gini 
#OOB prediction error:             23.81 % 

#OOB is the mean prediction error on each training sample xᵢ, 
#using only the trees that did not have xᵢ in their bootstrap sample.

#Next I will calculate model Kappa for the full model 
#(inherently imbalanced data sets so we are using Cohen's Kappa 
#to compare the models)

#install.packages('caret', dependencies = TRUE)
pred1 <- classify$predictions

kappa1 <- postResample(pred1, rf_data$response)[[2]]
kappa1
#[1] 0.3786982 
# This means that the classifier explains ~ 37 % of the data

#I will sort the data by feature importance
importances <- sort(importance(classify), decreasing = T)
#restrict the number of variables to something reasonable to help with memory management
if (length(importances) >= 250) {importances <- importances[1:250]}

# And add one feature at a time to the model in order of importance 
# and compare kappas

# reorder the data sets
rf_data_2 <- rf_data[c("response", names(importances))]
importances #These are the predictors in the order of decreasing importance
# I will run another rf model to get the "effect" of each predictor
comp_classify <- list(NULL)
kappa2 <- NA
for (k in 1:length(importances)) {
  #if all the importances are below the mean break the loop
  if (k == 0) {break}
  new_data <- data.frame(response=rf_data_2$response, rf_data_2[seq(2,k+1)])
  set.seed(42)
  comp_classify[[k]] <- ranger(response ~ ., data = new_data, num.trees=10000, importance="impurity")
  pred2 <- comp_classify[[k]]$predictions
  kappa2[k] <- postResample(pred2, new_data$response)[[2]]
  if (kappa2[k] == 1) {
    break
  }
}
if (length(kappa2) == length(importances)) {
  k <- min(which(kappa2 %in% max(kappa2)))
}
#this will store either values of the first "1" 
#or the highest kappa value of all the models run
all_model_errors_2[[1]] <- kappa2[k]
kappa2
all_model_errors_2
#[1] 0.8082192

#partial dependence plots of the best model using all of the data
nfeatures <- comp_classify[[k]]$num.independent.variables
rf_data2 <- rf_data_2[1:(nfeatures+1)]
pd <- partial_dependence(comp_classify[[k]], vars=colnames(rf_data2)[-1], data=rf_data2, n=c(25,nrow(rf_data2)))

plot_pd(pd)
# of cource I want to make a fancier plot, so I will save the data
plot_pd_data <- plot_pd(pd)$data
# processing
plot_pd_data$variable <- gsub("\\.", " ", plot_pd_data$variable)
unique(plot_pd_data$variable)
plot_pd_data$variable.f<-factor(plot_pd_data$variable, levels=unique(plot_pd_data$variable))

levels(plot_pd_data$variable.f)
# I need the genera names:
subset(taxonomy, OTU=="Otu036" | OTU== "Otu043" | OTU== "Otu028" | 
         OTU=="Otu024",
       select=genus)

#genus
#24               Olsenella(100)
#28         Acidaminococcus(100)
#36 Clostridia_unclassified(100)
#43             Collinsella(100)

plot_pd_data$variableII.f<-plot_pd_data$variable.f

plot_pd_data$variableII.f<-factor(plot_pd_data$variableII.f, 
                                  levels=c(levels(plot_pd_data$variableII.f),
                 "Olsenella","Acidaminococcus","Clostridia_uncl.","Collinsella","tnpA-06/IS1216_206"))

levels(plot_pd_data$variableII.f)

plot_pd_data<- within(plot_pd_data, variableII.f[variable=="Otu024"]<-"Olsenella")
plot_pd_data<- within(plot_pd_data, variableII.f[variable=="Otu028"]<-"Acidaminococcus")
plot_pd_data<- within(plot_pd_data, variableII.f[variable=="Otu036"]<-"Clostridia_uncl.")
plot_pd_data<- within(plot_pd_data, variableII.f[variable=="Otu043"]<-"Collinsella")
plot_pd_data<- within(plot_pd_data, variableII.f[variable=="tnpA 06 IS1216_206"]<-"tnpA-06/IS1216_206")


plot_pd_data$variableII.f<-droplevels(plot_pd_data$variableII.f)
unique(plot_pd_data$variableII.f)
plot_pd_data$variableII.f<-factor(plot_pd_data$variableII.f, levels=unique(plot_pd_data$variableII.f))


PDwrap<-ggplot(data = plot_pd_data, aes(value, prediction*100)) + geom_line(aes(colour=class), size= 1) +
  #scale_x_continuous(trans="log2") + 
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(colour = "black", fill = "white"))+
  #theme(strip.background = element_blank())+
  theme(panel.border = element_rect(colour = "black"))+
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1, size=9))+
  theme(axis.text.y = element_text(size=9))+
  labs(x="Relative abundance", y="Prediction (% chance to be classified)", color="Function") +
  #facet_grid(variableII.f ~., scales="free") + 
 facet_wrap( ~ variableII.f, ncol=3, scales="free_x")+
  theme(strip.text = element_text(size = 10, colour = "black"))+
  # scale_x_continuous(limits=c(0,0.2))+
  scale_y_continuous(limits = c(0, 100))+
  scale_color_manual(name = "Cluster", values = rev(cl_col))#+
  #theme(legend.position="none")
PDwrap

tSNE

grid.arrange(tSNE,PDwrap,nrow=1,ncol=2 )
