
return_combined_phenotypes <- function(merged_data){
  
  #--------------------------------------------------------------------------------------------------------------------------
  #T_regs
  merged_data$thres__T_regs <- ifelse(( merged_data$thres__CD4__entire_mean == 1 & merged_data$thres__FoxP3__nucleus_mean_unique == 1   ), 1, 0)
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  #T_helpers
  merged_data$thres__T_helps <- ifelse(( merged_data$thres__CD4__entire_mean == 1 & merged_data$thres__FoxP3__nucleus_mean_unique == 0   ), 1, 0)
  --------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  #PDL1_CD4   (Treg + Thelper)
  merged_data$thres__PDL1_CD4 <- ifelse(( merged_data$thres__CD4__entire_mean == 1 & merged_data$thres__PDL1__cyto_mean_unique == 1   ), 1, 0)
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  #PDL1_CD8
  merged_data$thres__PDL1_CD8 <- ifelse(( merged_data$thres__CD8__cyto_mean == 1 & merged_data$thres__PDL1__cyto_mean_unique == 1   ), 1, 0)
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  #PDL1_CD68
  merged_data$thres__PDL1_CD68 <- ifelse(( merged_data$thres__CD68__cyto_mean == 1 & merged_data$thres__PDL1__cyto_mean_unique == 1 ), 1, 0)
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  #PDL1_T_regs
  merged_data$thres__PDL1_T_regs__unique <- ifelse(( merged_data$thres__CD4__entire_mean == 1 
                                                   & merged_data$thres__FoxP3__nucleus_mean_unique == 1   
                                                   & merged_data$thres__PDL1__cyto_mean_unique == 1), 1, 0)
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  #PD1_CD4
  merged_data$thres__PD1_CD4 <- ifelse(( merged_data$thres__CD4__entire_mean == 1 & merged_data$thres__PD1__cyto_mean_unique == 1 ), 1, 0)
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  #PD1_CD8
  merged_data$thres__PD1_CD8 <- ifelse(( merged_data$thres__CD8__cyto_mean == 1 & merged_data$thres__PD1__cyto_mean_unique == 1 ), 1, 0)
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  #PD1_CD68
  merged_data$thres__PD1_CD68 <- ifelse(( merged_data$thres__CD68__cyto_mean == 1 & merged_data$thres__PD1__cyto_mean_unique == 1 ), 1, 0)
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  # rest phenotyping
  
  merged_data$thres__PDL1_PD1_CD4__unique  <- ifelse((merged_data$thres__CD4__entire_mean==1 & merged_data$thres__CD8__cyto_mean==0 & 
                                                      merged_data$thres__PDL1__cyto_mean_unique==1 & merged_data$thres__PD1__cyto_mean_unique==1 & 
                                                      merged_data$thres__FoxP3__nucleus_mean_unique==0 & merged_data$thres__CD68__cyto_mean==0), 1, 0)
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  merged_data$thres__PDL1_PD1_CD8__unique     <- ifelse((merged_data$thres__CD4__entire_mean==0 & merged_data$thres__CD8__cyto_mean==1 & 
                                                         merged_data$thres__PDL1__cyto_mean_unique==1 & merged_data$thres__PD1__cyto_mean_unique==1 & 
                                                         merged_data$thres__FoxP3__nucleus_mean_unique==0 & merged_data$thres__CD68__cyto_mean==0), 1, 0)
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  merged_data$thres__CD8_T_regs__unique <- ifelse((merged_data$thres__CD4__entire_mean==0 & 
                                                   merged_data$thres__CD8__cyto_mean==1 & 
                                                   merged_data$thres__FoxP3__nucleus_mean_unique==1 & 
                                                   merged_data$thres__CD68__cyto_mean==0), 1, 0)
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  merged_data$thres__PDL1_PD1_CD68__unique     <- ifelse((merged_data$thres__CD4__entire_mean==0 & merged_data$thres__CD8__cyto_mean==0 & 
                                                          merged_data$thres__PDL1__cyto_mean_unique==1 & merged_data$thres__PD1__cyto_mean_unique==1 & 
                                                          merged_data$thres__FoxP3__nucleus_mean_unique==0 & merged_data$thres__CD68__cyto_mean==1), 1, 0)
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  merged_data$thres__PanCK_PDL1_neg__unique     <- ifelse((merged_data$thres__CD4__entire_mean==0 & 
                                                           merged_data$thres__CD8__cyto_mean==0 & 
                                                           merged_data$thres__PDL1__cyto_mean_unique==0 & 
                                                           merged_data$thres__PD1__cyto_mean_unique==0 & 
                                                           merged_data$thres__FoxP3__nucleus_mean_unique==0 & 
                                                           merged_data$thres__CD68__cyto_mean==0 &
                                                           merged_data$thres__PanCK__cyto_mean_unique==1), 1, 0)
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  merged_data$thres__PanCK_PDL1_pos__unique     <- ifelse((merged_data$thres__CD4__entire_mean==0 & 
                                                           merged_data$thres__CD8__cyto_mean==0 & 
                                                           merged_data$thres__PDL1__cyto_mean_unique==1 & 
                                                           merged_data$thres__PD1__cyto_mean_unique==0 & 
                                                           merged_data$thres__FoxP3__nucleus_mean_unique==0 & 
                                                           merged_data$thres__CD68__cyto_mean==0 &
                                                           merged_data$thres__PanCK__cyto_mean_unique==1), 1, 0)
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  #immune cells
  merged_data$thres__immune_cells <- ifelse(( merged_data$thres__CD4__entire_mean == 1 | merged_data$thres__CD8__cyto_mean == 1 | merged_data$thres__CD68__cyto_mean == 1), 1, 0)
  
  return(merged_data)
}