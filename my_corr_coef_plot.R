plot.my_corr_coef <- function(x,
                           type = "lower",
                           diag = TRUE,
                           reorder = FALSE,
                           signif = c("stars", "pval"),
                           show = c("all", "signif"),
                           p_val = 0.05,
                           caption = TRUE,
                           digits.cor = 2,
                           digits.pval = 3,
                           col.low = "red",
                           col.mid = "white",
                           col.high = "blue",
                           lab.x.position = NULL,
                           lab.y.position = NULL,
                           legend.position = NULL,
                           legend.title = "Pearson's\nCorrelation",
                           size.text.cor = 3,
                           size.text.signif = 3,
                           size.text.lab = 10,
                           ...){
  signif <- rlang::arg_match(signif)
  show <- rlang::arg_match(show)
  
  
  
  
  pval <- x$pval
  correl <- x$cor
  if(reorder == TRUE){
    correl <- reorder_cormat(correl)
    pval <- pval[rownames(correl), colnames(correl)]
  }
  
  nonsign <- which(pval >= p_val)
  if(show == "signif"){
    correl[nonsign] <- NA
  }
  
  if(type == "lower"){
    pval <- my_make_lower_tri(pval)
    print('bike')
    correl <- my_make_lower_tri(correl)
  }
  if(type == "upper"){
    pval <- make_upper_tri(pval)
    correl <- make_upper_tri(correl)
  }
  
  
  
  bind_data <-
    expand.grid(dimnames(correl)) %>%
    mutate(cor = as.vector(correl),
           pval = as.vector(pval),
           pval_star = pval) %>%
    rlang::set_names("v1", "v2", "cor", "pval", "pval_star") %>%
    dplyr::filter(!is.na(pval))
  if(signif == "stars"){
    bind_data <-
      bind_data |>
      mutate(pval_star = case_when(pval < 0.001 ~ "***",
                                   between(pval, 0.001, 0.01) ~ "**",
                                   between(pval, 0.01, 0.05) ~ "*",
                                   pval >= 0.05 ~ "ns")) %>%
      rlang::set_names("v1", "v2", "cor", "pval", "pval_star") %>%
      dplyr::filter(!is.na(pval))
  } else{
    bind_data <-
      expand.grid(dimnames(correl)) %>%
      mutate(cor = as.vector(correl),
             pval = as.vector(pval),
             pval_star <- case_when(
               pval < 0.001 ~ "***",
               between(pval, 0.001, 0.01) ~ "**",
               between(pval, 0.01, 0.05) ~ "*",
               pval >= 0.05 ~ "ns"
             )) %>%
      rlang::set_names("v1", "v2", "cor", "pval", "pval_star") %>%
      dplyr::filter(!is.na(pval))
  }
  
  
  if(show == "signif"){
    bind_data <-
      bind_data |>
      mutate(pval_star = ifelse(pval >= p_val, NA, pval_star),
             pval = ifelse(pval >= p_val, NA, pval))
  }
  
  if(diag == FALSE){
    bind_data <- dplyr::filter(bind_data, !v1 == v2)
  }
  if(type == "lower"){
    lab.y.position <- ifelse(missing(lab.y.position), "right", lab.y.position)
    lab.x.position <- ifelse(missing(lab.x.position), "bottom", lab.x.position)
  }
  if(type == "upper"){
    lab.y.position <- ifelse(missing(lab.y.position), "left", lab.y.position)
    lab.x.position <- ifelse(missing(lab.x.position), "top", lab.x.position)
  }
  axis.text.x.hjust <- ifelse(lab.x.position == "bottom", 1, 0)
  if(missing(legend.position) & type == "lower"){
    legend.position <- c(0.2, 0.85)
  }
  if(missing(legend.position) & type == "upper"){
    legend.position <- c(0.8, 0.17)
  }
  if(!missing(legend.position)){
    legend.position <- legend.position
  }
  
  
  p <-
    ggplot(bind_data, aes(v1, v2, fill = cor)) +
    geom_tile(aes(fill = cor),
              colour = "gray") +
    geom_text(aes(label = ifelse(!is.na(cor), format(round(cor, digits.cor), nsmall = digits.cor), "")),
              vjust = 0,
              size = size.text.cor) +
    {if(signif == "stars")geom_text(aes(label = ifelse(!is.na(pval_star), pval_star, "")),
                                    vjust = 1.5,
                                    size = size.text.signif)} +
    {if(signif == "pval")geom_text(aes(label = ifelse(!is.na(pval), as.vector(format(signif(pval, digits = digits.pval), nsmall = digits.pval)), "")),
                                   vjust = 1.5,
                                   size = size.text.signif)} +
    scale_fill_gradient2(low = col.low,
                         high = col.high,
                         mid = col.mid,
                         midpoint = 0,
                         limit = c(-1, 1),
                         space = "Lab",
                         na.value = "transparent",
                         name = legend.title)+
    theme_bw() +
    xlab(NULL) +
    ylab(NULL) +
    scale_y_discrete(expand = expansion(mult = c(0,0)),
                     position = lab.y.position)+
    scale_x_discrete(expand = expansion(mult = c(0,0)),
                     position = lab.x.position) +
    coord_fixed() +
    theme(axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.position = legend.position,
          legend.direction = "horizontal",
          axis.text = element_text(color = "black", size = size.text.lab),
          axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = axis.text.x.hjust))+
    guides(fill = guide_colorbar(barwidth = 6,
                                 barheight = 1,
                                 title.position = "top",
                                 title.hjust = 0.5))
  if(signif == "stars" & caption == TRUE){
    p <-
      p + labs(caption = c("ns p >= 0.05; * p < 0.05; ** p < 0.01; and *** p < 0.001"))
  }
  suppressWarnings(return(p))
}