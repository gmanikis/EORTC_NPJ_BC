library('metan')
library('ggplot2')
library('table1')

my_make_lower_tri<-function(x, diag = NA){
  x[upper.tri(x)] <- NA
  #diag(x) <- diag
  return(x)
}

pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- kruskal.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- lalonde[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- kruskal.test(y ~ lalonde$age.quartiles)$p.value
    } else {
      p <- chisq.test(table(y, droplevels(lalonde$age.quartiles)))$p.value
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

# help functions, credits to @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
source(file.path(path_eortc_txt , 'my_corr_coef.R'))
source(file.path(path_eortc_txt , 'my_corr_coef_plot.R'))

#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
# FIGURE 3

merged_data_celldensinfo <- get_merged_data_celldensinfo(merged_data_celldensinfo)
merged_data_celldensinfo <- data.frame(merged_data_celldensinfo)

rownames(merged_data_celldensinfo) <- paste0('Pat_',merged_data_celldensinfo$patient_id)

merged_data_celldensinfo  <- merge( clinical_info_to_use, merged_data_celldensinfo, by='row.names') 

##remove NA pCR
merged_data_celldensinfo <- merged_data_celldensinfo %>% drop_na(PCR)

merged_data_celldensinfo <- select( merged_data_celldensinfo , -c("PDL1_in_Tumor_areas" ))

tumor_names <- colnames(merged_data_celldensinfo)[grepl("_T$", colnames(merged_data_celldensinfo))]
stroma_names <- colnames(merged_data_celldensinfo)[grepl("_S$", colnames(merged_data_celldensinfo))]
entire_names <- colnames(merged_data_celldensinfo)[(!colnames(merged_data_celldensinfo)%in%tumor_names) & (!colnames(merged_data_celldensinfo)%in%stroma_names)]

## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   ENTIRE AREA MAIN PHENOS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
main_pheno <- c('CD4'  ,  'CD8'  ,  'CD68'  ,  'T_regs')

entire_names_main_pheno <- merged_data_celldensinfo[  colnames(merged_data_celldensinfo)%in%entire_names  &  colnames(merged_data_celldensinfo)%in%main_pheno   ]
entire_names_main_pheno[] <- sapply( entire_names_main_pheno, as.numeric )
entire_names_main_pheno$ID <- merged_data_celldensinfo$patient_id

table1::label(entire_names_main_pheno$CD4) <- 'T helpers'
table1::label(entire_names_main_pheno$CD8) <- 'Cytotoxic T cells'
table1::label(entire_names_main_pheno$T_regs) <- 'Regulatory T cells'
table1::label(entire_names_main_pheno$CD68) <- 'Macrophages'

tbl1 <- table1::table1(~ CD4 + CD8 + CD68 + T_regs  , data = entire_names_main_pheno)

t1flex(tbl1) %>% 
  save_as_docx(path=file.path(save_figs ,   "Table3Cmain.docx"))

#--------

entire_names_main_pheno2 <- cbind.data.frame(  entire_names_main_pheno  ,  SubType = as.factor(merged_data_celldensinfo$SubType) )
table1::label(entire_names_main_pheno2$CD4) <- 'T helpers'
table1::label(entire_names_main_pheno2$CD8) <- 'Cytotoxic T cells'
table1::label(entire_names_main_pheno2$T_regs) <- 'Regulatory T cells'
table1::label(entire_names_main_pheno2$CD68) <- 'Macrophages'

tbl1 <- table1::table1(~  CD4 + CD8 + CD68 + T_regs | SubType , data = entire_names_main_pheno2,
                       overall = F,extra.col=list(`P-value`=pvalue))

sect_properties <- prop_section(page_size = page_size(orient = "landscape"))

t1flex(tbl1) %>% 
  save_as_docx(path=file.path(save_figs ,   "Table3Cmain_per_subtype.docx"))
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


entire_names_main_pheno_melt <-  melt(entire_names_main_pheno2, id.vars=c("ID",'SubType'))


levels(entire_names_main_pheno_melt$variable)[levels(entire_names_main_pheno_melt$variable)=='CD8'] <- 'Cytotoxic T cells'
levels(entire_names_main_pheno_melt$variable)[levels(entire_names_main_pheno_melt$variable)=='CD4'] <- 'T helpers'
levels(entire_names_main_pheno_melt$variable)[levels(entire_names_main_pheno_melt$variable)=='T_regs'] <- 'Regulatory T cells'
levels(entire_names_main_pheno_melt$variable)[levels(entire_names_main_pheno_melt$variable)=='CD68'] <- 'Macrophages'

cols2 <- c( 
  'Cytotoxic T cells' =  'chocolate3'  ,  
  "T helpers" = 'lightseagreen',
  'Regulatory T cells' = 'skyblue1',
  'Macrophages' = 'lightyellow3'
)

#main phenos
p3 <- ggplot(entire_names_main_pheno_melt, aes(factor(variable), y = value)) + 
  geom_violin(aes(fill = factor(variable)), trim = FALSE) +
  geom_boxplot(width = 0.1) +
  facet_wrap(~factor(variable),scales="free",nrow=1)+
  theme(legend.position="top",legend.title = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),text = element_text(size = 20) ,strip.text = element_text(face = "bold")  )+
  scale_y_continuous( trans = "log2",breaks = trans_breaks("log2", function(x) 2^x),
                      labels = trans_format("log2", math_format(2^.x)))+
  scale_colour_manual(values = cols2)+
  scale_fill_manual(values = cols2)

ggsave(
  file.path(save_figs ,   "Figure3C_main.png"),
  plot = print(p3),
  scale = 3,
  height = 800, width = 1200, units = "px" , dpi = 300,
  limitsize = FALSE,
  bg= 'black'
)

#-------------------------------------------------------------------------------------
cols2 <- c( 
  'HER2+' =  'chocolate3'  ,  
  "HR+HER2-" = 'lightseagreen',
  'TN' = 'skyblue1'
)
#main phenos per subtype
p7 <- ggplot(entire_names_main_pheno_melt, aes(factor(SubType), y = value)) + 
  geom_violin(aes(fill = factor(SubType)), trim = FALSE) +
  geom_boxplot(width = 0.1) +
  facet_wrap(~factor(variable),scales="free",nrow=3)+
  theme(legend.position="top",legend.title = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),text = element_text(size = 20) ,strip.text = element_text(face = "bold")  )+
  scale_y_continuous( trans = "log2",breaks = trans_breaks("log2", function(x) 2^x),
                      labels = trans_format("log2", math_format(2^.x)))+
  scale_colour_manual(values = cols2)+
  scale_fill_manual(values = cols2)

ggsave(
  file.path(save_figs ,   "Figure3C_main_per_subtype.png"),
  plot = print(p7),
  scale = 3,
  height = 800, width = 1000, units = "px" , dpi = 300,
  limitsize = FALSE,
  bg= 'black'
)


#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# FIGURE 4

main_pheno <- c('CD4'  ,  'CD8'  ,  'CD68'  ,  'T_regs')

merged_data_celldensinfo <- data.frame(merged_data_celldensinfo)

rownames(merged_data_celldensinfo) <- paste0('Pat_',as.character(merged_data_celldensinfo$patient_id))

work_tils <- select(work_tils , -c('PatID'))
merged_data_celldensinfo    <- select( merged_data_celldensinfo  , -c('Row.names','PatID'))


de <- merge( work_tils, merged_data_celldensinfo, by='row.names')  # merge by row names (by=0 or by="row.names")

de <-  subset( de, select = -c(Row.names) )

de_plot <- de[  colnames(de)%in%c(    names(work_tils)[names(work_tils)%in%names(cols)] , names(de)[names(de)%in%main_pheno])    ]
de_plot <- data.frame(sapply( de_plot, as.numeric ))


de_second_pheno_all <-     de[  colnames(de)%in%entire_names  &  !colnames(de)%in%main_pheno   ]   
de_second_pheno_all <- select( de_second_pheno_all  ,   -c('PR_TMA','SubType','patient_id')  )
de_second_pheno_all <- data.frame(sapply( de_second_pheno_all, as.numeric ))

colnames(de_second_pheno_all)[colnames(de_second_pheno_all) == 'PDL1_CD4'] <- "CD4+ PD-L1+"
colnames(de_second_pheno_all)[colnames(de_second_pheno_all) == 'PDL1_CD8'] <- "CD8+ PD-L1+"
colnames(de_second_pheno_all)[colnames(de_second_pheno_all) == 'PDL1_CD68'] <- "CD68+ PD-L1+"
colnames(de_second_pheno_all)[colnames(de_second_pheno_all) == 'PDL1_T_regs'] <- "FoxP3+ PD-L1+"
colnames(de_second_pheno_all)[colnames(de_second_pheno_all) == 'PD1_CD4'] <- "CD4+ PD-1+"
colnames(de_second_pheno_all)[colnames(de_second_pheno_all) == 'PD1_CD8'] <- "CD8+ PD-1+"
colnames(de_second_pheno_all)[colnames(de_second_pheno_all) == 'PD1_CD68'] <- "CD68+ PD-1+"
colnames(de_second_pheno_all)[colnames(de_second_pheno_all) == 'PD1_T_regs'] <- "FoxP3+ PD-1+"
de_second_pheno_all <- select(de_second_pheno_all , -c('Row.names'))


#-------------------------------------------------------------------------------
de2 <- de[  colnames(de)%in%c(names(de)[names(de)%in%main_pheno])   ] 
tmp <- colnames(de2)
de2 <- data.frame(sapply( de2, as.numeric ))
colnames(de2) <- tmp


colnames(de2)[colnames(de2) == 'CD4'] <- 'T helpers'
colnames(de2)[colnames(de2) == 'CD8'] <- 'Cytotoxic T cells'
colnames(de2)[colnames(de2) == 'T_regs'] <- 'Regulatory T cells'
colnames(de2)[colnames(de2) == 'CD68'] <- "Macrophages"
#-------------------------------------------------------------------------------
zzz <- paste(c('CD4'  ,  'CD8'  ,  'CD68'  ,  'T_regs'),'S',sep='_')
de2_stroma <- de[  colnames(de)%in%c(names(de)[names(de)%in%zzz   ])   ] 
tmp <- colnames(de2_stroma)
de2_stroma <- data.frame(sapply( de2_stroma, as.numeric ))
colnames(de2_stroma) <- tmp

colnames(de2_stroma)[colnames(de2_stroma) == 'CD4_S'] <- 'T helpers'
colnames(de2_stroma)[colnames(de2_stroma) == 'CD8_S'] <- 'Cytotoxic T cells'
colnames(de2_stroma)[colnames(de2_stroma) == 'T_regs_S'] <- 'Regulatory T cells'
colnames(de2_stroma)[colnames(de2_stroma) == 'CD68_S'] <- "Macrophages"
#-------------------------------------------------------------------------------
zzz <- paste(c('CD4'  ,  'CD8'  ,  'CD68'  ,  'T_regs'),'T',sep='_')
de2_tumor <- de[  colnames(de)%in%c(names(de)[names(de)%in%zzz   ])   ] 
tmp <- colnames(de2_tumor)
de2_tumor <- data.frame(sapply( de2_tumor, as.numeric ))
colnames(de2_tumor) <- tmp

colnames(de2_tumor)[colnames(de2_tumor) == 'CD4_T'] <- 'T helpers'
colnames(de2_tumor)[colnames(de2_tumor) == 'CD8_T'] <- 'Cytotoxic T cells'
colnames(de2_tumor)[colnames(de2_tumor) == 'T_regs_T'] <- 'Regulatory T cells'
colnames(de2_tumor)[colnames(de2_tumor) == 'CD68_T'] <- "Macrophages"


data_ge2 <- cbind.data.frame(  de[  colnames(de)%in%c(names(work_tils)[names(work_tils)%in%names(cols)])   ]   , de_second_pheno_all )

xxxx <- my_corr_coef(data_ge2,method='spearman')

png(file.path(save_figs ,   paste0("Figure4_extra_ALL_new.png")),height = 5000, width = 5000, units = "px" , res = 300)

plot.my_corr_coef(xxxx,diagonal=FALSE,reorder = FALSE,type='lower',
     legend.title = "Spearman's\nCorrelation",
     col.low = "#543005",
     col.mid = "white",
     col.high = "#003C30",
     size.text.cor = 5,
     size.text.signif = 5,
     size.text.lab = 15)

dev.off()

data_ge2 <- cbind.data.frame(  de[  colnames(de)%in%c(names(work_tils)[names(work_tils)%in%names(cols)])   ]   , de2 )

my_cor2 <-corr.test(de[  colnames(de)%in%c(names(work_tils)[names(work_tils)%in%names(cols)])   ]   , de2 , use = "pairwise",method="spearman",adjust="none")

png(file.path(save_figs ,   paste0("Figure4_main_ALL.png")),height = 5000, width = 5000, units = "px" , res = 300)

xxxx <- my_corr_coef(data_ge2,method='spearman')

plot.my_corr_coef(xxxx,diagonal=FALSE,reorder = FALSE,type='lower',
                  legend.title = "Spearman's\nCorrelation",
                  col.low = "#543005",
                  col.mid = "white",
                  col.high = "#003C30",
                  size.text.cor = 5,
                  size.text.signif = 5,
                  size.text.lab = 15)



dev.off()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
my_cor2 <-corr.test(de[  colnames(de)%in%c(names(work_tils)[names(work_tils)%in%names(cols)])   ]   , de2_stroma , use = "pairwise",method="spearman",adjust="none")

data_ge2 <- cbind.data.frame( de[  colnames(de)%in%c(names(work_tils)[names(work_tils)%in%names(cols)])   ]    , de2_stroma )

png(file.path(save_figs ,   paste0("Figure4_main_STROMA.png")),height = 5000, width = 5000, units = "px" , res = 300)

xxxx <- my_corr_coef(data_ge2,method='spearman')

plot.my_corr_coef(xxxx,diagonal=FALSE,reorder = FALSE,type='lower',
                  legend.title = "Spearman's\nCorrelation",
                  col.low = "#543005",
                  col.mid = "white",
                  col.high = "#003C30",
                  size.text.cor = 5,
                  size.text.signif = 5,
                  size.text.lab = 15)

dev.off()
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

my_cor2 <-corr.test(de[  colnames(de)%in%c(names(work_tils)[names(work_tils)%in%names(cols)])   ]   , de2_tumor , use = "pairwise",method="spearman",adjust="none")

png(file.path(save_figs ,   paste0("Figure4_main_TUMOR.png")),height = 5000, width = 5000, units = "px" , res = 300)

data_ge2 <- cbind.data.frame( de[  colnames(de)%in%c(names(work_tils)[names(work_tils)%in%names(cols)])   ]    , de2_tumor )
xxxx <- my_corr_coef(data_ge2,method='spearman')
plot.my_corr_coef(xxxx,diagonal=FALSE,reorder = FALSE,type='lower',
                  legend.title = "Spearman's\nCorrelation",
                  col.low = "#543005",
                  col.mid = "white",
                  col.high = "#003C30",
                  size.text.cor = 5,
                  size.text.signif = 5,
                  size.text.lab = 15)


dev.off()
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


#--------------------------------------------------------------
de_plot_stroma <- de[  colnames(de)%in%c(    names(work_tils)[names(work_tils)%in%names(cols)] , names(de)[names(de)%in%stroma_names])    ]
de_plot_stroma <- data.frame(sapply( de_plot_stroma, as.numeric ))

de2_stroma <- de[  colnames(de)%in%c(names(de)[names(de)%in%stroma_names])   ] 
tmp <- colnames(de2_stroma)
de2_stroma <- data.frame(sapply( de2_stroma, as.numeric ))
colnames(de2_stroma) <- tmp

de2_stroma <- select( de2_stroma  , -c('CD4_S' , 'CD8_S' , 'CD68_S' , 'T_regs_S'))

colnames(de2_stroma)[colnames(de2_stroma) == 'PDL1_CD4_S'] <- "CD4+ PD-L1+"
colnames(de2_stroma)[colnames(de2_stroma) == 'PDL1_CD8_S'] <- "CD8+ PD-L1+"
colnames(de2_stroma)[colnames(de2_stroma) == 'PDL1_CD68_S'] <- "CD68+ PD-L1+"
colnames(de2_stroma)[colnames(de2_stroma) == 'PDL1_T_regs_S'] <- "FoxP3+ PD-L1+"
colnames(de2_stroma)[colnames(de2_stroma) == 'PD1_CD4_S'] <- "CD4+ PD-1+"
colnames(de2_stroma)[colnames(de2_stroma) == 'PD1_CD8_S'] <- "CD8+ PD-1+"
colnames(de2_stroma)[colnames(de2_stroma) == 'PD1_CD68_S'] <- "CD68+ PD-1+"
colnames(de2_stroma)[colnames(de2_stroma) == 'PD1_T_regs_S'] <- "FoxP3+ PD-1+"


    
my_cor2 <-corr.test(de[  colnames(de)%in%c(names(work_tils)[names(work_tils)%in%names(cols)])   ]   , de2_stroma , use = "pairwise",method="spearman",adjust="none")

png(file.path(save_figs ,   paste0("Figure4_extra_STROMA.png")),height = 5000, width = 5000, units = "px" , res = 300)

data_ge2 <- cbind.data.frame( de[  colnames(de)%in%c(names(work_tils)[names(work_tils)%in%names(cols)])   ]    , de2_stroma )
xxxx <- my_corr_coef(data_ge2,method='spearman')
plot.my_corr_coef(xxxx,diagonal=FALSE,reorder = FALSE,type='lower',
                  legend.title = "Spearman's\nCorrelation",
                  col.low = "#543005",
                  col.mid = "white",
                  col.high = "#003C30",
                  size.text.cor = 5,
                  size.text.signif = 5,
                  size.text.lab = 15)


dev.off()


#--------------------------------------------------------------
de_plot_tumor <- de[  colnames(de)%in%c(    names(work_tils)[names(work_tils)%in%names(cols)] , names(de)[names(de)%in%tumor_names])    ]
de_plot_tumor <- data.frame(sapply( de_plot_tumor, as.numeric ))

de2_tumor <- de[  colnames(de)%in%c(names(de)[names(de)%in%tumor_names])   ] 
tmp <- colnames(de2_tumor)
de2_tumor <- data.frame(sapply( de2_tumor, as.numeric ))
colnames(de2_tumor) <- tmp

de2_tumor <- select( de2_tumor  , -c('CD4_T' , 'CD8_T' , 'CD68_T' , 'T_regs_T'))

colnames(de2_tumor)[colnames(de2_tumor) == 'PDL1_CD4_T'] <- "CD4+ PD-L1+"
colnames(de2_tumor)[colnames(de2_tumor) == 'PDL1_CD8_T'] <- "CD8+ PD-L1+"
colnames(de2_tumor)[colnames(de2_tumor) == 'PDL1_CD68_T'] <- "CD68+ PD-L1+"
colnames(de2_tumor)[colnames(de2_tumor) == 'PDL1_T_regs_T'] <- "FoxP3+ PD-L1+"
colnames(de2_tumor)[colnames(de2_tumor) == 'PD1_CD4_T'] <- "CD4+ PD-1+"
colnames(de2_tumor)[colnames(de2_tumor) == 'PD1_CD8_T'] <- "CD8+ PD-1+"
colnames(de2_tumor)[colnames(de2_tumor) == 'PD1_CD68_T'] <- "CD68+ PD-1+"
colnames(de2_tumor)[colnames(de2_tumor) == 'PD1_T_regs_T'] <- "FoxP3+ PD-1+"


my_cor2 <-corr.test(de[  colnames(de)%in%c(names(work_tils)[names(work_tils)%in%names(cols)])   ]   , de2_tumor , use = "pairwise",method="spearman",adjust="none")

png(file.path(save_figs ,   paste0("Figure4_extra_TUMOR.png")),height = 5000, width = 5000, units = "px" , res = 300)

data_ge2 <- cbind.data.frame( de[  colnames(de)%in%c(names(work_tils)[names(work_tils)%in%names(cols)])   ]    , de2_tumor )
xxxx <- my_corr_coef(data_ge2,method='spearman')
plot.my_corr_coef(xxxx,diagonal=FALSE,reorder = FALSE,type='lower',
                  legend.title = "Spearman's\nCorrelation",
                  col.low = "#543005",
                  col.mid = "white",
                  col.high = "#003C30",
                  size.text.cor = 5,
                  size.text.signif = 5,
                  size.text.lab = 15)

dev.off()

