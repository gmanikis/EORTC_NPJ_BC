redundant_phenos <- function(merged_data){
  thres_variables <- dplyr::select(merged_data,contains("thres__"))
  thres_variables <- thres_variables %>%mutate_all(as.character)
  all_expressions <- data.frame( var_name = rep(NA, nrow(thres_variables)))
  
  for (i in 1:ncol(thres_variables)){ 
    tmp <- thres_variables[[i]]
    tmp <- data.frame( var_name = ifelse( tmp == '1' , as.character(colnames(thres_variables)[i]) , NA) ) 
    
    tmp <- separate(data = tmp , col = var_name, into = c("redan1", "expression" , 'redan2'), sep = "__")
    tmp = subset(tmp, select = c(expression))
    
    all_expressions <- data.frame( var_name = paste( all_expressions$var_name , tmp$expression , sep = '__'))
  }
  
  all_expressions$var_name <- gsub("NA__", "", all_expressions$var_name)
  all_expressions$var_name <- gsub("__NA", "", all_expressions$var_name)
  all_expressions$var_name <- gsub("NA", "negative", all_expressions$var_name)
  
  merged_data <- cbind.data.frame( Phenotype =  all_expressions$var_name , merged_data )
  
  
  #-------------------------------------------------------------------
  merged_data <- merged_data[!grepl('negative',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD4__PDL1__PanCK__T_helps__PDL1_CD4__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('PDL1__PD1__PanCK',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD4__PDL1__CD8__T_helps__PDL1_CD4__PDL1_CD8__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD68__CD4__PDL1__T_helps__PDL1_CD4__PDL1_CD68__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD4__CD8__PD1__T_helps__PD1_CD4__PD1_CD8__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD4__CD8__T_helps__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD4__PDL1__CD8__PD1__T_helps__PDL1_CD4__PDL1_CD8__PD1_CD4__PD1_CD8__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD68__CD4__PDL1__PD1__T_helps__PDL1_CD4__PDL1_CD68__PD1_CD4__PD1_CD68__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD68__PD1__PD1_CD68__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('PDL1__CD8__PanCK__PDL1_CD8__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD68__CD8__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD68__CD4__PD1__T_helps__PD1_CD4__PD1_CD68__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD68__CD4__T_helps__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD68__PDL1__CD8__PDL1_CD8__PDL1_CD68__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD4__PDL1__FoxP3__PD1__T_regs__PDL1_CD4__PDL1_T_regs__PD1_CD4__PD1_T_regs__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD68__FoxP3__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('PDL1__CD8__PD1__PanCK__PDL1_CD8__PD1_CD8__PDL1_PD1_CD8__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('PDL1__FoxP3__PanCK',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD4__PDL1__PD1__PanCK__T_helps__PDL1_CD4__PD1_CD4__PDL1_PD1_CD4__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD68__CD4__PDL1__CD8__T_helps__PDL1_CD4__PDL1_CD8__PDL1_CD68__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD68__CD8__PD1__PD1_CD8__PD1_CD68__immune_cells',merged_data$Phenotype),]
  merged_data <- merged_data[!grepl('CD4__FoxP3__PD1__T_regs__PD1_CD4__PD1_T_regs__immune_cells',merged_data$Phenotype),]
  
  ### 
  tmp <- data.frame(props = table(merged_data$Phenotype))
  tmp <- tmp[tmp$props.Freq<=24,]
  merged_data <- merged_data[ !merged_data$Phenotype%in%tmp$props.Var1  , ]
  #-------------------------------------------------------------------
  
  return(merged_data)
}