
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(meann = median(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  colnames(data_sum)[colnames(data_sum)%in%c('meann')] <- 'value'
  return(data_sum)
}
#____________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________
create_beautiful_radarchart <- function(data, color = "#00AFBB", seg ,
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, 
    seg = seg,
    pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 1,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}
#____________________________________________________________________________________________________________________________________________________

data_for_radar_plot <- function(study_cohort3 , cell_name , cell_name2 , clinic_outcome, pcr_flag){
  
  biomarker_names <- colnames(study_cohort3[,(grepl(paste0('EntropyGrasSlope__' , cell_name ,   '_CancerCells__Radi_'), colnames(study_cohort3)))])
  
  radar_colnames <- c(  paste0('Radius 50μm') ,
                        paste0('Radius 70μm') ,
                        paste0('Radius 100μm') ,
                        paste0('Radius 150μm') ,
                        paste0('Radius 200μm') ,
                        paste0('Radius 250μm') ,
                        
                        paste0('Radius 300μm'),
                        paste0('Radius 350μm'),
                        paste0('Radius 400μm'),
                        paste0('Radius 450μm'),
                        paste0('Radius 500μm'),
                        paste0('Radius 550μm'),
                        paste0('Radius 600μm'),
                        paste0('Radius 650μm'),
                        paste0('Radius 700μm')
  )
  
  
  radar_colnames <- str_replace(biomarker_names, paste0('EntropyGrasSlope__' , cell_name ,   '_CancerCells__Radi_'), "Radius ")
  
  
  
  if(pcr_flag==0)
  {
    biomarker_names <- c(  biomarker_names , 'p53')
    pcr_TN <- study_cohort3[,colnames(study_cohort3)%in%biomarker_names]
    pcr_TN_0 <- pcr_TN[pcr_TN$p53==1,!colnames(pcr_TN)%in%c('p53')]
    pcr_TN_1 <- pcr_TN[pcr_TN$p53==2,!colnames(pcr_TN)%in%c('p53')]
    
  }else{
    biomarker_names <- c(  biomarker_names , 'PCR')
    pcr_TN <- study_cohort3[,colnames(study_cohort3)%in%biomarker_names]
    pcr_TN_0 <- pcr_TN[pcr_TN$PCR==0,!colnames(pcr_TN)%in%c('PCR')]
    pcr_TN_1 <- pcr_TN[pcr_TN$PCR==1,!colnames(pcr_TN)%in%c('PCR')]
  }
  
  
  pcr_TN_0_mean <- as.data.frame(t(apply(pcr_TN_0,2,function(x) median(x, na.rm = TRUE))))  
  pcr_TN_1_mean <- as.data.frame(t(apply(pcr_TN_1,2,function(x) median(x, na.rm = TRUE))))  
  
  pcr_TN_mean <- rbind.data.frame(pcr_TN_0_mean , pcr_TN_1_mean)
  
  if(pcr_flag==0)
  {rownames(pcr_TN_mean) <- c('Wild Type' , 'Mutated')}else{rownames(pcr_TN_mean) <- c('no PCR' , 'PCR')}
  
  colnames(pcr_TN_mean) <- radar_colnames
  pcr_TN_mean <- rbind(rep(1,length(radar_colnames))   , rep(0,length(radar_colnames)) , pcr_TN_mean)
  
  return(pcr_TN_mean)
}

#____________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________
data_for_radar_plot_mix_score <- function(study_cohort3 , cell_name , cell_name2 , clinic_outcome, pcr_flag){
  
  biomarker_names <- c(
    paste0('MixScoreNorm__50' , cell_name ,   '_CancerCells') ,
    paste0('MixScoreNorm__100' , cell_name ,   '_CancerCells') ,
    paste0('MixScoreNorm__150' , cell_name ,   '_CancerCells') ,
    paste0('MixScoreNorm__200' , cell_name ,   '_CancerCells') ,
    paste0('MixScoreNorm__250' , cell_name ,   '_CancerCells'),
    
    paste0('MixScoreNorm__300' , cell_name ,   '_CancerCells'),
    paste0('MixScoreNorm__350' , cell_name ,   '_CancerCells'),
    paste0('MixScoreNorm__400' , cell_name ,   '_CancerCells'),
    paste0('MixScoreNorm__450' , cell_name ,   '_CancerCells'),
    paste0('MixScoreNorm__500' , cell_name ,   '_CancerCells'),
    paste0('MixScoreNorm__550' , cell_name ,   '_CancerCells'),
    paste0('MixScoreNorm__600' , cell_name ,   '_CancerCells')
  )
  
  radar_colnames <- c(
    paste0('Radius 50μm') ,
    paste0('Radius 100μm') ,
    paste0('Radius 150μm') ,
    paste0('Radius 200μm') ,
    paste0('Radius 250μm') ,
    
    paste0('Radius 300μm'),
    paste0('Radius 350μm'),
    paste0('Radius 400μm'),
    paste0('Radius 450μm'),
    paste0('Radius 500μm'),
    paste0('Radius 550μm'),
    paste0('Radius 600μm')
  )
  
  if(pcr_flag==0)
  {
    biomarker_names <- c(  biomarker_names , 'p53')
    pcr_TN <- study_cohort3[,colnames(study_cohort3)%in%biomarker_names]
    pcr_TN_0 <- pcr_TN[pcr_TN$p53==1,!colnames(pcr_TN)%in%c('p53')]
    pcr_TN_1 <- pcr_TN[pcr_TN$p53==2,!colnames(pcr_TN)%in%c('p53')]
    
  }else{
    biomarker_names <- c(  biomarker_names , 'PCR')
    pcr_TN <- study_cohort3[,colnames(study_cohort3)%in%biomarker_names]
    pcr_TN_0 <- pcr_TN[pcr_TN$PCR==0,!colnames(pcr_TN)%in%c('PCR')]
    pcr_TN_1 <- pcr_TN[pcr_TN$PCR==1,!colnames(pcr_TN)%in%c('PCR')]
  }
  
  
  pcr_TN_0_mean <- as.data.frame(t(apply(pcr_TN_0,2,function(x) median(x, na.rm = TRUE))))  
  pcr_TN_1_mean <- as.data.frame(t(apply(pcr_TN_1,2,function(x) median(x, na.rm = TRUE))))  
  
  pcr_TN_mean <- rbind.data.frame(pcr_TN_0_mean , pcr_TN_1_mean)
  
  if(pcr_flag==0)
  {rownames(pcr_TN_mean) <- c('Wild Type' , 'Mutated')}else{rownames(pcr_TN_mean) <- c('no PCR' , 'PCR')}
  
  colnames(pcr_TN_mean) <- radar_colnames
  pcr_TN_mean <- rbind(rep(1,length(radar_colnames))   , rep(0,length(radar_colnames)) , pcr_TN_mean)
  
  return(pcr_TN_mean)
}

#____________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________

get_merged_data_celldensinfo <- function(merged_data_celldensinfo){
  
  
  merged_data_celldensinfo <- separate(data = merged_data_celldensinfo, col = Case, into = c("tmp", "patient_id"), sep = "__")
  merged_data_celldensinfo$patient_id <- str_replace_all(merged_data_celldensinfo$patient_id, 'T','')
  merged_data_celldensinfo$patient_id <- str_replace_all(merged_data_celldensinfo$patient_id, ' ','')
  merged_data_celldensinfo <- merged_data_celldensinfo[ , -which(names(merged_data_celldensinfo) %in% c("tmp"))]
  
  #remove cells that belong to patients with no clinical information
  merged_data_celldensinfo <- merged_data_celldensinfo[merged_data_celldensinfo$patient_id%in%as.character(clinical_info$PATID),]
  
  return(merged_data_celldensinfo)
}
#____________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________

# convert spe object to a data frame with only colData
get_colData <- function(spe_object){
  formatted_data <- data.frame(SummarizedExperiment::colData(spe_object))
  formatted_data <- cbind(formatted_data, 
                          data.frame(SpatialExperiment::spatialCoords(spe_object)))
  if (is.null(formatted_data$Cell.ID)){
    formatted_data <- formatted_data %>% tibble::rownames_to_column("Cell.ID")
  }
  
  # delete column `sample_id`
  formatted_data$sample_id <- NULL
  
  return(formatted_data)
}
#____________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________

get_field_info = function(path) {
  info = readTIFFDirectory(path, all=FALSE)
  
  required_attributes = c('width', 'length', 'x.resolution', 'resolution.unit')
  missing_attributes = setdiff(required_attributes, names(info))
  if (length(missing_attributes) > 0) {
    missing = paste(missing_attributes, collapse=', ')
    stop(paste0('Image file is missing required attributes: ', missing))
  }
  
  if (info$resolution.unit != 'cm')
    stop(paste('Unsupported resolution unit:', info$resolution.unit))
  
  result = list()
  result$image_size = c(info$width, info$length)
  result$microns_per_pixel = as.numeric(10000/info$x.resolution)
  result$pixels_per_micron = 1/result$microns_per_pixel
  result$field_size = result$image_size * result$microns_per_pixel
  
  # Location directly from TIFF info
  result$location = c(info$x.position, info$y.position) * 10000
  result
}
#____________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________
readTIFFDirectory = function(path, all=FALSE) {
  # This function is a shim that calls
  # readTIFF from the latest s-u/tiff using the interface from
  # akoyabio/tiff readTIFFDirectory,
  if ('payload' %in% names(formals(tiff::readTIFF))) {
    # s-u/tiff returns a data.frame
    info = tiff::readTIFF(path, all=all, payload=FALSE) %>%
      purrr::transpose()
    if (!all)
      info = info[[1]]
    info
  } else
    stop('Please install a more recent tiff package.')
}

#____________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________
