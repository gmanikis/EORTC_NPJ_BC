custom_create_phenotype_from_thres <- function(input_TMA){
  thres_variables <- dplyr::select(input_TMA,contains("thres__"))
  
  thres_variables <- thres_variables %>%
    mutate_all(as.character)
  
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
  
  input_TMA <- cbind.data.frame( Phen_For_Estimations =  all_expressions$var_name , input_TMA )
}