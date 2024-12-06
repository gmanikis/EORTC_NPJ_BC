return_embedded_phenos <- function(merged_data){
  
  required_densities <- c( 'Cytoplasm.CD68.Mean' ,  
                           'Entire.Cell.CD4.Mean',
                           'Cytoplasm.PD.L1.Mean',
                           'Cytoplasm.CD8.Mean',
                           'Nucleus.FoXP3.Mean',
                           'Cytoplasm.PD.1.Mean',
                           'Cytoplasm.PanCK.Mean'
                           )
  
  #### keep only necessary infos
  extra_cols <- c('Phen_For_Estimations','Phenotype','patient_id','sample','Tissue.Category','Center.Name','Centre_TMA','Cell.ID','Cell.X.Position.um','Cell.Y.Position.um')
  extra_cols <- c( extra_cols ,  required_densities  )
  
  #find columns containing the required_densities
  columns_need <- unlist(lapply(  extra_cols   ,   grep , colnames(merged_data)))
  
  useful_data  <- merged_data[, columns_need] 
  
  #create 2 different Phenotype Columns ---> Cell.Type.Basic  & Cell.Type.Advanced
  
  #2
  useful_data <- useful_data[!grepl('CD4__PDL1__CD8__PanCK__T_helps__PDL1_CD4__PDL1_CD8__immune_cells',useful_data$Phen_For_Estimations),]
  #3
  useful_data <- useful_data[!grepl('CD4__PDL1__CD8__PD1__PanCK__T_helps__PDL1_CD4__PDL1_CD8__PD1_CD4__PD1_CD8__immune_cells',useful_data$Phen_For_Estimations),]
  #8
  useful_data <- useful_data[!grepl('CD68__CD4__PDL1__CD8__PD1__T_helps__PDL1_CD4__PDL1_CD8__PDL1_CD68__PD1_CD4__PD1_CD8__PD1_CD68__immune_cells',useful_data$Phen_For_Estimations),]
  #10
  useful_data <- useful_data[!grepl('CD8__T_regs',useful_data$Phen_For_Estimations),]
  
  #===========================================================================

  #1
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'CD4__PD1__T_helps__PD1_CD4__immune_cells__PD1_T_helps'] <- 'T_helpers'
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'CD4__PD1__T_helps__PD1_CD4__immune_cells__PD1_T_helps'] <- 'T_helpers_exhausted'
  
  #4
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'CD4__PDL1__PD1__T_helps__PDL1_CD4__PD1_CD4__PDL1_PD1_CD4__immune_cells__PDL1_PD1_T_helps'] <- 'T_helpers'
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'CD4__PDL1__PD1__T_helps__PDL1_CD4__PD1_CD4__PDL1_PD1_CD4__immune_cells__PDL1_PD1_T_helps'] <- 'T_helpers_exhausted'
  
  #5
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'CD4__PDL1__T_helps__PDL1_CD4__immune_cells__PDL1_T_helps'] <- 'T_helpers'
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'CD4__PDL1__T_helps__PDL1_CD4__immune_cells__PDL1_T_helps'] <- 'T_helpers_exhausted'
  
  #6
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'CD4__T_helps__immune_cells__PDL1neg_T_helps'] <- 'T_helpers'
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'CD4__T_helps__immune_cells__PDL1neg_T_helps'] <- 'T_helpers'
  
  #7
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'CD68'] <- 'Macrophages'
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'CD68'] <- 'Macrophages'
  
  #9
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'CD8'] <- 'Cyto_T_cells'
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'CD8'] <- 'Cyto_T_cells'
  
  #11
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'FoxP3'] <- 'T_regs'
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'FoxP3'] <- 'T_regs'
  
  #12
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PanCK__PanCK_PDL1_neg'] <- 'CancerCells'
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PanCK__PanCK_PDL1_neg'] <- 'CancerCells_PDL1neg'
  
  #13
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PD1'] <- 'PD1'          #### SOS check with Vangelis!!!!!!!!!!s
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PD1'] <- 'PD1'
  
  #14
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PD1__CD8'] <- 'Cyto_T_cells_exhausted'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PD1__CD8'] <- 'Cyto_T_cells_exhausted'
  
  #15
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PD1__PDL1__CD68'] <- 'PDL1pos_Macrophages_RARE'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PD1__PDL1__CD68'] <- 'PDL1pos_Macrophages_RARE'
  
  #16
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PD1__T_helps'] <- 'T_helpers'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PD1__T_helps'] <- 'T_helpers_exhausted'
  
  #17
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PD1__T_regs'] <- 'T_regs'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PD1__T_regs'] <- 'T_regs'
  
  #18
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PDL1'] <- 'PDL1'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PDL1'] <- 'PDL1'
  
  #19
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PDL1__CD68'] <- 'PDL1pos_Macrophages_RARE'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PDL1__CD68'] <- 'PDL1pos_Macrophages_RARE'
  
  #20
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PDL1__CD8'] <- 'Cyto_T_cells_exhausted'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PDL1__CD8'] <- 'Cyto_T_cells_exhausted'
  
  #21
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PDL1__CD8__FoxP3__PDL1_CD8__CD8_T_regs__immune_cells'] <- 'SPECIAL_CLASS'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PDL1__CD8__FoxP3__PDL1_CD8__CD8_T_regs__immune_cells'] <- 'SPECIAL_CLASS'
  
  #22
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PDL1__FoxP3__PD1'] <- 'T_regs'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PDL1__FoxP3__PD1'] <- 'T_regs'
  
  #23
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PDL1__FoxP3__PDL1_FoxP3'] <- 'T_regs'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PDL1__FoxP3__PDL1_FoxP3'] <- 'T_regs'
  
  #24
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PDL1__PanCK__PanCK_PDL1_pos'] <- 'CancerCells'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PDL1__PanCK__PanCK_PDL1_pos'] <- 'CancerCells_PDL1pos'
  
  #25
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PDL1__PD1__CD8'] <- 'Cyto_T_cells_exhausted'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PDL1__PD1__CD8'] <- 'Cyto_T_cells_exhausted'
  
  #26
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PDL1__PD1__PDL1_PD1'] <- 'OTHER_exhausted'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PDL1__PD1__PDL1_PD1'] <- 'OTHER_exhausted'
  
  #27
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PDL1__T_regs'] <- 'T_regs'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PDL1__T_regs'] <- 'T_regs'
  
  #28
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'PDL1neg__T_helps'] <- 'T_helpers'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'PDL1neg__T_helps'] <- 'T_helpers'
  
  #29
  useful_data$Cell.Type.Basic[useful_data$Phen_For_Estimations == 'T_regs'] <- 'T_regs'          
  useful_data$Cell.Type.Advanced[useful_data$Phen_For_Estimations == 'T_regs'] <- 'T_regs'
  
  # 1                                                        CD4__PD1__T_helps__PD1_CD4__immune_cells__PD1_T_helps   2039    ---> T_helps      //     T_helps__exhausted
  # 2                                             CD4__PDL1__CD8__PanCK__T_helps__PDL1_CD4__PDL1_CD8__immune_cells     64    ---->   OUT
  # 3                      CD4__PDL1__CD8__PD1__PanCK__T_helps__PDL1_CD4__PDL1_CD8__PD1_CD4__PD1_CD8__immune_cells     27     ---->   OUT
  # 4                     CD4__PDL1__PD1__T_helps__PDL1_CD4__PD1_CD4__PDL1_PD1_CD4__immune_cells__PDL1_PD1_T_helps   1485     ---->  T_helps      //     T_helps__exhausted
  # 5                                                     CD4__PDL1__T_helps__PDL1_CD4__immune_cells__PDL1_T_helps   6993     ---->  T_helps      //     T_helps__exhausted
  # 6                                                                  CD4__T_helps__immune_cells__PDL1neg_T_helps   5617     ---->  T_helps  
  # 7                                                                                                         CD68  11751     ----->  Macrophages
  # 8  CD68__CD4__PDL1__CD8__PD1__T_helps__PDL1_CD4__PDL1_CD8__PDL1_CD68__PD1_CD4__PD1_CD8__PD1_CD68__immune_cells     41    ---->   OUT
  # 9                                                                                                          CD8  22691    -----> Cytotoxic Tcells
  # 10                                                                                                 CD8__T_regs     56    ---->   OUT
  # 11                                                                                                       FoxP3   3958    ---->   Tregs
  # 12                                                                                       PanCK__PanCK_PDL1_neg 296923    ---->   CancerCells_PDL1neg
  # 13                                                                                                         PD1  12432    ---->   
  # 14                                                                                                    PD1__CD8   2016    ---->   Cytotoxic Tcells_exhausted
  # 15                                                                                             PD1__PDL1__CD68    124    ---->   PDL1pos_Macrophages
  # 16                                                                                                PD1__T_helps     98    ---->  T_helps      //     T_helps__exhausted
  # 17                                                                                                 PD1__T_regs    128    ---->   Tregs
  # 18                                                                                                        PDL1 122012    ---->  PDL1
  # 19                                                                                                  PDL1__CD68   2227    ---->   PDL1pos_Macrophages
  # 20                                                                                                   PDL1__CD8   4753    ---->   Cytotoxic Tcells_exhausted
  # 21                                                        PDL1__CD8__FoxP3__PDL1_CD8__CD8_T_regs__immune_cells     25    ---->   SPECIAL_class
  # 22                                                                                            PDL1__FoxP3__PD1     29   ---->   Tregs
  # 23                                                                                     PDL1__FoxP3__PDL1_FoxP3    735   ---->   Tregs
  # 24                                                                                 PDL1__PanCK__PanCK_PDL1_pos  73009       ---->   CancerCells_PDL1pos
  # 25                                                                                              PDL1__PD1__CD8    542    ---->   Cytotoxic Tcells_exhausted
  # 26                                                                                         PDL1__PD1__PDL1_PD1   1743    ---->  OTHER_exhausted
  # 27                                                                                                PDL1__T_regs    247    ---->   Tregs
  # 28                                                                                            PDL1neg__T_helps    838    ---> T_helps 
  # 29                                                                                                      T_regs    208    ---->   Tregs
  
  
  # T_helps__exhausted   ---> pro_tumoral
  # T_helps              ---> anti tumoral
  # Macrophages          ---> pro_tumoral
  # Cytotoxic Tcells     ---> anti tumoral
  # Tregs                ---> pro_tumoral
  #Cytotoxic Tcells_exhausted  ---> pro_tumoral
  
  
  return(useful_data)
}