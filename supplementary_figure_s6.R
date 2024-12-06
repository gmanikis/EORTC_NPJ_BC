supplementary_figure_s6 <- function(NEW_eortc_fixed_coords, match_reference, multiplex_raw_images_info, TMAs_to_keep, folder_path){


#################################################################################################
######## PRINT THE IMAGES
pat_IDs_eortc <- unique(NEW_eortc_fixed_coords$patient_id)
pheno_to_print <- 'Cell.Type.Basic'

cols <- c( 
  'Tumor Cells' = 'white' ,
  'Macrophages'  = 'darkgoldenrod1'  ,
  'T Helpers' =  'violet'  ,  
  "Regulatory T Cells" = 'blue'  , 
  'Cytotoxic T Cells Exhausted' = 'darkgreen'  ,
  'Cytotoxic T Cells'  = 'red'      )

#
for (ooo in 1:length(pat_IDs_eortc))
{
  
  print(paste0( 'Patient with ID:  ' , pat_IDs_eortc[ooo], '   starts'))
  pat_IDs_eortc[ooo]
  
  patient_match_ref <-  match_reference[match_reference$patient_id==pat_IDs_eortc[ooo],]
  patient_match_ref_codes <- patient_match_ref$code
  patient_match_ref_centre <- unique(patient_match_ref$Centre)
  
  patient_images <- multiplex_raw_images_info[ (multiplex_raw_images_info$patID%in%patient_match_ref_codes) &  (multiplex_raw_images_info$centre%in%paste0(patient_match_ref_centre,"_")   )   ,]
  
  xx <- NEW_eortc_fixed_coords %>% filter(Center.Name ==  patient_match_ref_centre )
  xx <- xx %>% filter(patient_id ==  pat_IDs_eortc[ooo] )
  
  
  ######### remove classes in plots
  xx <- xx[  !xx$Cell.Type.Basic%in%c('OTHER_exhausted' , 'PD1'  ,  'PDL1pos_Macrophages_RARE' ,  'PDL1') ,]
  
  xx$Cell.Type.Basic[xx$Cell.Type.Basic=='CancerCells'] <- 'Tumor Cells'
  xx$Cell.Type.Basic[xx$Cell.Type.Basic=='Cyto_T_cells'] <- 'Cytotoxic T Cells'
  xx$Cell.Type.Basic[xx$Cell.Type.Basic=='Cyto_T_cells_exhausted'] <- 'Cytotoxic T Cells Exhausted'
  xx$Cell.Type.Basic[xx$Cell.Type.Basic=='T_helpers'] <- 'T Helpers'
  xx$Cell.Type.Basic[xx$Cell.Type.Basic=='T_regs'] <- 'Regulatory T Cells'
  
  p <- NULL
  pixels_per_micron_STATS <- data.frame()
  
  if(is.null(p)) {
    p <<-
      ggplot(xx,                                                                           ############# Basic OR Advanced phenotyping plot
             aes(x = Cell.X.Position, y = Cell.Y.Position, colour = Cell.Type.Basic))
  } else{
    p <<-
      p + ggplot(xx,
                 aes(x = Cell.X.Position, y = Cell.Y.Position, colour = Cell.Type.Basic))
  }
  #___________________________________________________________________________
  background <- vector(mode = "list", length = nrow(patient_images))
  
  for (ii in 1:nrow(patient_images))
    #for (ii in 2:2)
  {
    
    zzz <- filter(TMAs_to_keep, patient_ID == pat_IDs_eortc[ooo]  &   patID  == patient_images$patID[ii])
    if(zzz$keep_from_eortic){
      
      background[[ii]] = jpeg::readJPEG(file.path(imaging_path_eortc , patient_images$image_path[ii])) %>% as.raster()
      
      pixels_per_micron_eortc <-
        cbind.data.frame(
          x_axis = patient_images$x_axis[ii]  ,
          y_axis = patient_images$y_axis[ii]           + 1000*ii  ,
          x_size = x_size_   ,
          y_size = y_size_
        )
      pixels_per_micron_STATS <-
        rbind.data.frame(
          pixels_per_micron_STATS ,
          cbind.data.frame(pixels_per_micron_eortc , image_file = patient_images$image_path[ii])
        )
      
      ylim = c((pixels_per_micron_eortc$y_axis)    ,
               (pixels_per_micron_eortc$y_axis) + pixels_per_micron_eortc$y_size
      )
      xlim = c((pixels_per_micron_eortc$x_axis)    ,
               (pixels_per_micron_eortc$x_axis) + pixels_per_micron_eortc$x_size
      )
      
      if (is.null(p)) {
        p <<-
          ggplot2::annotation_raster(
            background[[ii]],
            xmin = xlim[1],
            xmax = xlim[2],
            ymin = -ylim[1],
            ymax = -ylim[2]
          )
      } else{
        p <<-
          p + ggplot2::annotation_raster(
            background[[ii]],
            xmin = xlim[1],
            xmax = xlim[2],
            ymin = -ylim[1],
            ymax = -ylim[2]
          )
      }
      #___________________________________________________________________________
    }
  }
  x_limit_min <- round(min(pixels_per_micron_STATS$x_axis))
  x_limit_max <- round(max(pixels_per_micron_STATS$x_axis) + max(pixels_per_micron_STATS$x_size))
  y_limit_min <- round(min(pixels_per_micron_STATS$y_axis))
  y_limit_max <- round(max(pixels_per_micron_STATS$y_axis) + max(pixels_per_micron_STATS$y_size))
  
  
  mikrometer <- 50
  
  line_coords_x <- c(round(x_limit_min) ,  round(x_limit_min) + mikrometer)
  line_coords_y <- c(round(y_limit_max) ,  round(y_limit_max))
  
  p <-     p +
    theme(
      legend.text = element_text(size = 10),legend.title = element_blank(),
      # legend.title = element_text(size = 10),
      legend.position = "top",
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "black", colour = "black")
    )  + xlab('') + ylab('') +
    
    guides(color = guide_legend(override.aes = list(size = 10)))+
    ggplot2::scale_y_reverse(limits = rev(c(y_limit_min, y_limit_max)))+
    scale_x_continuous(limits = c(x_limit_min, x_limit_max))+
    scale_colour_manual(values = cols) +
    geom_point(size = 0.7) +
    
    ggplot2::geom_segment(
      ggplot2::aes(
        x = line_coords_x[1],
        xend = line_coords_x[2],
        y = line_coords_y[1],
        yend = line_coords_y[2]
      ),
      color = 'white',
      alpha = 1,
      linewidth = 3,
      na.rm = TRUE
    ) +
    
    ggplot2::geom_text(
      ggplot2::aes(
        x = round(x_limit_min) + mikrometer / 2,
        y = round(y_limit_max) - 60,
        label = paste(mikrometer, '~mu*m')
      ),
      size = 10,
      hjust = 0.5,
      vjust = 1,
      color = 'white',
      alpha = 1,
      parse = TRUE
    )

  ggsave(
    file.path(
      path_eortc_txt ,
      folder_path,
      'eortc_image_results' ,
      paste0('EORTC_' ,   pat_IDs_eortc[ooo] , "_all_TMAs.png")
    ),
    plot = print(p),
    scale = 3,
    height = round(y_limit_max - y_limit_min),
    width = round(x_limit_max - x_limit_min),
    units = "px" ,
    dpi = 300,
    limitsize = FALSE,
    bg = 'black'
  )
  
}
}

