library('ggplot2')
library('haven')
library('tidyr')
library('plyr')
library('dplyr')
library('data.table')
library('stringr')
library('writexl')
library('openxlsx')
library('raster')
library('splitstackshape')
library('reshape2')

library('readr')

library('imcRtools')
library('spatstat')
library('SPIAT')
library('gtsummary')

library('ComplexHeatmap')
library('fmsb')
library('gridExtra')
#____________________________________________________________________________________________________________________________________________________

path_eortc_txt <- "txt file with phenotyping"
imaging_path_eortc <- "meta images from informs"
tiff_path <- 'folder where tiff images from inform are stored'

folder_path <- 'optional: additional folder'

#predefined pseudocolouring for phenotypes
cols <- c( 
  'CancerCells' = 'white' ,
  'CancerCells_PDL1neg' = 'white',
  'CancerCells_PDL1pos' = 'gray45',
  'T_helpers' =  'violet'  ,  
  'T_helpers_exhausted' =  'violetred4'  ,  
  'Cyto_T_cells' = 'darkgreen'  ,   
  'Cyto_T_cells_exhausted'      = 'aquamarine' ,
  'OTHER_exhausted'      = 'chartreuse2' ,
  "T_regs" = 'blue'  , 
  "PD1" = 'black'  , 
  'PDL1'  = 'red' ,
  'Macrophages'  = 'darkgoldenrod1'  ,
  'PDL1pos_Macrophages_RARE' = 'magenta',
  'SPECIAL_CLASS'   =  'mediumpurple'                )


source(file.path(path_eortc_txt , folder_path , 'GMAN_EORTIC_format_inform_to_spe.R'))
source(file.path(path_eortc_txt , folder_path ,'GMAN_EORTIC_entropy_gradients.R'))
source(file.path(path_eortc_txt , folder_path ,'GMAN_EORTIC_format_inform_to_spe_immune.R'))
source(file.path(path_eortc_txt , folder_path ,'GMAN_EORTIC_mixing_score_summary.R'))
source(file.path(path_eortc_txt , folder_path ,'help_functions_for_spatial_analysis.R'))
source(file.path(path_eortc_txt , folder_path ,'supplementary_figure_s6.R'))
#===============================================================================================================================
#===============================================================================================================================
#===============================================================================================================================
#===============================================================================================================================

#load data
merged_data_celldensinfo <- readxl::read_xlsx(file.path(path_eortc_txt , folder_path,'load cell densities'))


### load clinical info
clinical_info <- read_dta(file.path(path_eortc_txt , folder_path, 'load clinical data from dta' )) 
#----------------------------------------------------------------------------------------------------------------------------------------

#load images full paths
all_files <- list.files(  imaging_path_eortc )
multiplex_raw_images  <- str_subset(all_files, 'composite_image.jpg'   )

tiff_file <- list.files(  tiff_path )
tiff_file  <- str_subset(tiff_file, '\\.tif'   )
#----------------------------------------------------------------------------------------------------------------------------------------

#retrieve image size information
file_info <- get_field_info( file.path( tiff_path  , tiff_file ) )
x_size_ = file_info$image_size[1]* file_info$microns_per_pixel     ############# ----> equal to field_size[2]
y_size_ = file_info$image_size[2]* file_info$microns_per_pixel     ############# ----> equal to field_size[1]


#===============================================================================================================================
#===============================================================================================================================
#===============================================================================================================================
#===============================================================================================================================
#    SPATIAL ANALYSIS START



##### split data per patient
all_samples.split_sample_name <- merged_data_celldensinfo %>% group_by(patient_id)

#split for each patient  NEW with updates
all_samples.split_sample_name.names <- group_keys(all_samples.split_sample_name)
all_samples.split_sample_name.data  <- group_split(all_samples.split_sample_name)

#_________________________________________________________________________________________________________________________________________________________
#1 CancerCells vs T_helpers 
test_sample_cancercelss_thelpers <- lapply(all_samples.split_sample_name.data  , GMAN_EORTIC_format_inform_to_spe,
                                           tissue_ref = 'Entire' , pheno_ref = 'CancerCells' ,
                                           tissue_target = 'Entire' , pheno_target = 'T_helpers' , ref_advanced = FALSE , target_advanced = FALSE)

#___________________________________________________________________________________________________________________________________________________________
#2 CancerCells vs Cyto_T_cells 
test_sample_cancercelss_cytotcells <- lapply(all_samples.split_sample_name.data  , GMAN_EORTIC_format_inform_to_spe,
                                           tissue_ref = 'Entire' , pheno_ref = 'CancerCells' ,
                                           tissue_target = 'Entire' , pheno_target = 'Cyto_T_cells' , ref_advanced = FALSE , target_advanced = FALSE)

#___________________________________________________________________________________________________________________________________________________________
#3 CancerCells vs T_regs 
test_sample_cancercelss_tregs <- lapply(all_samples.split_sample_name.data  , GMAN_EORTIC_format_inform_to_spe,
                                             tissue_ref = 'Entire' , pheno_ref = 'CancerCells' ,
                                             tissue_target = 'Entire' , pheno_target = 'T_regs' , ref_advanced = FALSE , target_advanced = FALSE)

#___________________________________________________________________________________________________________________________________________________________
#4 CancerCells vs Macrophages 
test_sample_cancercelss_macrophages <- lapply(all_samples.split_sample_name.data  , GMAN_EORTIC_format_inform_to_spe,
                                        tissue_ref = 'Entire' , pheno_ref = 'CancerCells' ,
                                        tissue_target = 'Entire' , pheno_target = 'Macrophages' , ref_advanced = FALSE , target_advanced = FALSE)

#___________________________________________________________________________________________________________________________________________________________
#5 CancerCells vs  T_helpers_exhausted
test_sample_cancercelss_thelpersexhausted <- lapply(all_samples.split_sample_name.data  , GMAN_EORTIC_format_inform_to_spe,
                                              tissue_ref = 'Entire' , pheno_ref = 'CancerCells' ,
                                              tissue_target = 'Entire' , pheno_target = 'T_helpers_exhausted' , ref_advanced = FALSE , target_advanced = TRUE)

#___________________________________________________________________________________________________________________________________________________________
#6 CancerCells vs  Cyto_T_cells_exhausted
test_sample_cancercelss_cytotcellsexhausted <- lapply(all_samples.split_sample_name.data  , GMAN_EORTIC_format_inform_to_spe,
                                                    tissue_ref = 'Entire' , pheno_ref = 'CancerCells' ,
                                                    tissue_target = 'Entire' , pheno_target = 'Cyto_T_cells_exhausted' , ref_advanced = FALSE , target_advanced = TRUE)

#___________________________________________________________________________________________________________________________________________________________
#7 CancerCells
test_sample_cancercelss <- lapply(all_samples.split_sample_name.data  , GMAN_EORTIC_format_inform_to_spe,
                                                      tissue_ref = 'Entire' , pheno_ref = 'CancerCells' ,
                                                      tissue_target = 'Entire' , pheno_target = 'CancerCells' , ref_advanced = FALSE , target_advanced = FALSE)

#___________________________________________________________________________________________________________________________________________________________
#8 ALL entire
test_sample_ALL_entire <- lapply(all_samples.split_sample_name.data  , GMAN_EORTIC_format_inform_to_spe,
                                           tissue_ref = 'Entire' , pheno_ref = 'ALL' ,
                                           tissue_target = 'Entire' , pheno_target = 'ALL' , ref_advanced = FALSE , target_advanced = FALSE)

#9 CancerCells VS Immune
test_sample_cancercelss_immune <- lapply(all_samples.split_sample_name.data  , GMAN_EORTIC_format_inform_to_spe_immune,
                                 tissue_ref = 'Entire' , pheno_ref = 'CancerCells' ,
                                 tissue_target = 'Entire' , pheno_target = c('T_helpers','T_helpers_exhausted','Cyto_T_cells',
                                                                             'Cyto_T_cells_exhausted','T_regs','Macrophages','OTHER_exhausted'),
                                 ref_advanced = FALSE , target_advanced = FALSE, immune_phenos = c('T_helpers','T_helpers_exhausted','Cyto_T_cells',
                                                                                                   'Cyto_T_cells_exhausted','T_regs','Macrophages','OTHER_exhausted') , 'ImmuneCells')
#___________________________________________________________________________________________________________________________________________________________


#===============================================================================================================================
#===============================================================================================================================
#===============================================================================================================================
#===============================================================================================================================
#1 ENTROPY GRADIENTS

#a sample below 

gradient_positions <- seq(10,600,10)

entropy_gradient_sample <- GMAN_EORTIC_entropy_gradients(test_sample_cancercelss_immune, gradient_positions , 
                                      cell_types_of_interest =c("Immunes","CancerCells"), feature_colname= "inform_spatial")

rownames(entropy_gradient_sample) <- paste0('Pat_',all_samples.split_sample_name.names$patient_id)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#===============================================================================================================================
#===============================================================================================================================
#===============================================================================================================================
#===============================================================================================================================
#2 MIXING SCORE

mixing_score_sample <- GMAN_EORTIC_mixing_score_summary(test_sample_cancercelss_tregs, gradient_positions ,  "T_regs","CancerCells", "inform_spatial")
rownames(mixing_score_sample) <- paste0('Pat_',all_samples.split_sample_name.names$patient_id)

