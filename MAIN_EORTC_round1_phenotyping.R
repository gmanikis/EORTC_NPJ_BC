
library('haven')
library('tidyr')
library('stringr')
library('SPIAT')
library('plyr')
library('dplyr')
library(data.table)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
path_eortc_txt <- "inform metadata folder"
path <- file.path('working path')
setwd(path)

#help functions
source('return_basic_phenotypes.R')
source('return_combined_phenotypes.R')
source('redundant_phenos.R')
source('return_MERGE_phenotypes.R')
source('custom_create_phenotype_from_thres.R')
source('return_embedded_phenos.R')
source('make_custom_fcs.R')
source('return_MERGE_Phen_For_Estimations_2.R')

delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) < n,]
}

#-------------------------------------------------------------------------------
#load clinical data
clinical_info <- read_dta(paste0(path,' ....dta')) #load clinical data from dta file 

#load InForm meta-data in txt format
merged_data <- read.table('.....txt',sep = '\t', header = T)
#-------------------------------------------------------------------------------


################################################################################
#Step #1 return basic phenotyping using inform classification
merged_data <- return_basic_phenotypes(merged_data)

#filtering: exclude panCK in stroma
merged_data <- merged_data[  !(merged_data$thres__PanCK__cyto_mean_unique == 1 & merged_data$Tissue.Category == 'Stroma')    , ]

#Step #2 return combined phenotyping using inform classification (e.g. T Helpers)
merged_data <- return_combined_phenotypes(merged_data)

#Step #3 do exhaustive filtering due to mixed phenotyping
merged_data <- redundant_phenos(merged_data)

#Step #4 merged phenotypes, again do filtering
merged_data <- return_MERGE_phenotypes(merged_data)

#Step #5 help function
merged_data <- custom_create_phenotype_from_thres(merged_data)

#Step #6 embedded sub-phenotypes
merged_data <- return_embedded_phenos(merged_data)

#multiplex phenotyping ready!
################################################################################
