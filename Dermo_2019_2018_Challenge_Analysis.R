# Dermo Challenge Sept. 2019 Data Analysis 

# Code is optimized for performing analysis and starts with analyzing the in vivo Dermo challenge data from 2019 with 1 family and then 
# reanalyzes the summer data from 2018 Dermo in vivo challenge that was not optimized in last years script (Dermo_Apoptosis_Assay_Data_analysis.R)

#### Load Packages ####
library(metafor)
library(car)
library(lsmeans)
library(rgr)
library(multcompView)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(readxl)
library(magrittr)
library(stringr)
library(data.table)
library(ggsignif)
library(scales)
library(export) # graph2ppt() function
library(broom)
library(ggbiplot)
library(ggtext)
library(factoextra)
library(tidyverse)

##### LOAD ASSAY CSV'S and format #####
getwd()
# In excel
  # 1. Copied the first row,column header for the plot name to gates for that plot where the name
     # was not preserved. Copied the second row column headers for counts and % of this plot that were the same
  # 2. Addded in column 1 rows 1,2,3 "ID", "Treat", and "Assay"
  # 3. Saved this new file in the "FORMATTED_CSVS" folder
  
## VIABILITY ASSAY

# Make new header column
VIA_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/ANALYSIS_FILES/FORMATTED_CSVS/DAY7_VIABILITY_ASSAY_OPTIMIZED_WORKSPACE_DATA_header_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Day7_2019_VIA <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/ANALYSIS_FILES/FORMATTED_CSVS/DAY7_VIABILITY_ASSAY_OPTIMIZED_WORKSPACE_DATA_header_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(VIA_nms)

# Split column 1 by space and remove
Day7_2019_VIA <- Day7_2019_VIA %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
Day7_2019_VIA <- Day7_2019_VIA[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
Day7_2019_VIA_percent <- Day7_2019_VIA[,c(1:3,5,7,9,11)]
Day7_2019_VIA_counts <- Day7_2019_VIA[,c(1:4,6,8,10,12:15)]

# Gather count and percent columns separately
Day7_2019_VIA_percent <- Day7_2019_VIA_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:7))
Day7_2019_VIA_counts <-  Day7_2019_VIA_counts %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:11))

Day7_2019_VIA_percent <- Day7_2019_VIA_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Day7_2019_VIA_counts <-  Day7_2019_VIA_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

Day7_2019_VIA_percent <- Day7_2019_VIA_percent[,-c(4,6,8)]
Day7_2019_VIA_counts <-  Day7_2019_VIA_counts[,-c(4,6,8)]

# Full join together so all columns preserved
Day7_2019_VIA_join <- full_join(Day7_2019_VIA_percent, Day7_2019_VIA_counts, by = c("ID","Treat","Assay","Plot_number","Gate"))

# Add in cell type column based on plot number and gate
VIA_cell_type <- data.frame(Plot_number = c("2","3","4","4","9","10","6","6"), Gate= c("M4","P1","E1","E3","This","This","V2-L","V2-R"), Cell_type=c(
                                                               "all_hemocytes","all_hemocytes", "granular","agranular",
                                                               "live_granular", "live_agranular", "all_live_hemocytes","all_dead_hemocytes"))
# Join cell type
Day7_2019_VIA_join <- left_join(Day7_2019_VIA_join, VIA_cell_type, by=c("Plot_number","Gate"))
unique_Day7_2019_VIA_join_gate <- unique(Day7_2019_VIA_join [,c("Plot_number","Gate")])

## APOPTOSIS ASSAY
# Make new header column
Apop_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/ANALYSIS_FILES/FORMATTED_CSVS/DAY7_APOPTOSIS_ASSAY_OPTIMIZED_WORKSPACE_DATA_header_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Day7_2019_APOP <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/ANALYSIS_FILES/FORMATTED_CSVS/DAY7_APOPTOSIS_ASSAY_OPTIMIZED_WORKSPACE_DATA_header_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(Apop_nms)

# Split column 1 by space and remove
Day7_2019_APOP <- Day7_2019_APOP %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
Day7_2019_APOP <- Day7_2019_APOP[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
Day7_2019_APOP_percent <- Day7_2019_APOP[,c(1:3,7,9,11,13,15,17,19,21, 23,25,27,29,31,33,35)]
Day7_2019_APOP_counts <-  Day7_2019_APOP[,c(1:6,8,10,12,14,16,18,20,22,24,26,28,30,32)]

# Gather count and percent columns separately
Day7_2019_APOP_percent <- Day7_2019_APOP_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:18))
Day7_2019_APOP_counts <-  Day7_2019_APOP_counts %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:19))

Day7_2019_APOP_percent <- Day7_2019_APOP_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Day7_2019_APOP_counts <-  Day7_2019_APOP_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

Day7_2019_APOP_percent <- Day7_2019_APOP_percent[,-c(4,6,8)]
Day7_2019_APOP_counts <-  Day7_2019_APOP_counts[,-c(4,6,8)]

# Full join together so all columns preserved
Day7_2019_APOP_join <- full_join(Day7_2019_APOP_percent, Day7_2019_APOP_counts, by = c("ID","Treat","Assay","Plot_number","Gate"))

# Add in cell type column based on plot number and gate
APOP_cell_type <- data.frame(Plot_number = c("13","12","3","8","8","4","4","4","4","7","7","7","7","6","6","6","6"), 
                             Gate= c("This","This","P1","P3","P4","Q7-UL","Q7-UR","Q7-LL","Q7-LR","Q3-UL","Q3-UR","Q3-LL","Q3-LR","Q8-UL","Q8-UR","Q8-LL","Q8-LR"), 
                             Cell_type=c("live_granular", "live_agranular", "all_hemocytes", "agranular", "granular", 
                                         "all_hemocytes_dead", "all_hemocytes_dead_apoptotic", "all_hemocytes_live", 
                                         "all_hemocytes_live_apoptotic", "agranular_dead", "agranular_dead_apoptotic",
                                         "agranular_live", "agranular_live_apoptotic", "granular_dead", "granular_dead_apoptotic", 
                                         "granular_live", "granular_live_apoptotic"))

# Join cell type
Day7_2019_APOP_join <- left_join(Day7_2019_APOP_join, APOP_cell_type, by=c("Plot_number","Gate"))
unique_Day7_2019_APOP_join <- unique(Day7_2019_APOP_join[,c("Plot_number","Gate")])

## CASPASE ASSAY
# Make new header column
Casp_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/ANALYSIS_FILES/FORMATTED_CSVS/DAY7_CASPASE_ASSAY_OPTIMIZED_WORKSPACE_DATA_header_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Day7_2019_CASP <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/ANALYSIS_FILES/FORMATTED_CSVS/DAY7_CASPASE_ASSAY_OPTIMIZED_WORKSPACE_DATA_header_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(Casp_nms)

# Split column 1 by space and remove
Day7_2019_CASP <- Day7_2019_CASP %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
Day7_2019_CASP <- Day7_2019_CASP[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
Day7_2019_CASP_percent <- Day7_2019_CASP[,c(1:3,4,7,9,11,13,15,17,19,21, 23,25,27,29,31,33)]
Day7_2019_CASP_counts <-  Day7_2019_CASP[,c(1:3,5,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,35)]

# Gather count and percent columns separately
Day7_2019_CASP_percent <- Day7_2019_CASP_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:18))
Day7_2019_CASP_counts <-  Day7_2019_CASP_counts %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:20))

Day7_2019_CASP_percent <- Day7_2019_CASP_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Day7_2019_CASP_counts <-  Day7_2019_CASP_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

Day7_2019_CASP_percent <- Day7_2019_CASP_percent[,-c(4,6,8)]
Day7_2019_CASP_counts <-  Day7_2019_CASP_counts[,-c(4,6,8)]

# Full join together so all columns preserved
Day7_2019_CASP_join <- full_join(Day7_2019_CASP_percent, Day7_2019_CASP_counts, by = c("ID","Treat","Assay","Plot_number","Gate"))

# Add in cell type column based on plot number and gate
CASP_cell_type <- data.frame(Plot_number = c("3","8","8","4","4","4","4","7","7","7","7","2","2","2","2","6","9"), 
                             Gate= c("P1","P3","P4","Q2-UL","Q2-UR","Q2-LL","Q2-LR","Q5-UL","Q5-UR","Q5-LL","Q5-LR","Q1-UL","Q1-UR","Q1-LL","Q1-LR","This","This"), 
                             Cell_type=c("all_hemocytes", "agranular", "granular", "agranular_dead", "agranular_dead_caspase_active", 
                                         "agranular_live", "agranular_live_caspase_active", "granular_dead", "granular_dead_caspase_active", 
                                         "granular_live", "granular_live_caspase_active", "all_dead", "all_dead_caspase_active", "all_live", 
                                         "all_live_caspase_active", "live_agranular", "live_granular"))

# Join cell type
Day7_2019_CASP_join <- left_join(Day7_2019_CASP_join, CASP_cell_type, by=c("Plot_number","Gate"))
unique_Day7_2019_CASP_join <- unique(Day7_2019_CASP_join[,c("Plot_number","Gate")])

## LMP ASSAY
# Make new header column
LMP_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/ANALYSIS_FILES/FORMATTED_CSVS/DAY7_LMP_ASSAY_OPTIMIZED_WORKSPACE_DATA_header_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Day7_2019_LMP <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/ANALYSIS_FILES/FORMATTED_CSVS/DAY7_LMP_ASSAY_OPTIMIZED_WORKSPACE_DATA_header_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(LMP_nms)

# Split column 1 by space and remove
Day7_2019_LMP <- Day7_2019_LMP %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
Day7_2019_LMP <- Day7_2019_LMP[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
Day7_2019_LMP_percent <- Day7_2019_LMP[,c(1:3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33)]
Day7_2019_LMP_counts <- Day7_2019_LMP[,c(1:4,6,8,10,12,14,16,18,20,22,24,26,28,30,32)]

# Gather count and percent columns separately
Day7_2019_LMP_percent <- Day7_2019_LMP_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:18))
Day7_2019_LMP_counts <- Day7_2019_LMP_counts %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:18))

Day7_2019_LMP_percent <- Day7_2019_LMP_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Day7_2019_LMP_counts <- Day7_2019_LMP_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

Day7_2019_LMP_percent <- Day7_2019_LMP_percent[,-c(4,6,8)]
Day7_2019_LMP_counts <- Day7_2019_LMP_counts[,-c(4,6,8)]

# Full join together so all columns preserved
Day7_2019_LMP_join <- full_join(Day7_2019_LMP_percent, Day7_2019_LMP_counts, by = c("ID","Treat","Assay","Plot_number","Gate"))

# Add in cell type column based on plot number and gate
LMP_cell_type <- data.frame(Plot_number = c("3","8", "8","4","4","4","4", "7","7","7","7","6","6","6","6"),     
           Gate= c("P1","P3","P4","Q4-UL","Q4-UR","Q4-LL","Q4-LR","Q3-UL","Q3-UR","Q3-LL","Q3-LR","Q2-UL","Q2-UR","Q2-LL","Q2-LR"),
           Cell_type=c("all_hemocytes","agranular","granular","all_hemocytes_total_lysosome_rupture","all_hemocytes_intact_dead_other_means",
                       "all_hemocytes_onset_LMP","all_hemocytes_intact_lysosome_live","agranular_lysosome_rupture",
                       "agranular_dead_other_means","agranular_LMP","agranular_lysosome_live","granular_lysosome_rupture",
                       "granular_dead_other_means","granular_LMP","granular_lysosome_live"))
# Join cell type
Day7_2019_LMP_join <- left_join(Day7_2019_LMP_join, LMP_cell_type, by=c("Plot_number","Gate"))
  
#### Remove bad samples from all ####
#SAMPLES TO REMOVE FROM EACH ANALYSIS = 265 (notched control), 193 (Dermo injected), 218 (Dermo injected), 221 (Dermo injected)

# vector with ID lines for each assay
ID_samples_remove <- c("265","193","218","221")

# List of assay dataframes
assay_list <- list(Day7_2019_VIA_join=Day7_2019_VIA_join, Day7_2019_APOP_join=Day7_2019_APOP_join, 
                   Day7_2019_CASP_join=Day7_2019_CASP_join, Day7_2019_LMP_join=Day7_2019_LMP_join)

# Save lapply output after removing samples
remove_temp <- lapply(assay_list, function(x) x[ ! x$ID %in% ID_samples_remove, ])
# set names of each lapply output object in correct order and save using list2env
names(remove_temp) <- c("Day7_2019_VIA_bad_removed",  "Day7_2019_APOP_bad_removed", 
                        "Day7_2019_CASP_bad_removed", "Day7_2019_LMP_bad_removed")
list2env(remove_temp, envir = .GlobalEnv)

# check output
unique(Day7_2019_VIA_bad_removed$ID) # output is missing correct samples and has correct format 

#### Arcsine Transformation of percentages ####

# Make list of dataframes
bad_removed_assay_list <- list(Day7_2019_VIA_bad_removed=Day7_2019_VIA_bad_removed, Day7_2019_APOP_bad_removed=Day7_2019_APOP_bad_removed, 
                   Day7_2019_CASP_bad_removed=Day7_2019_CASP_bad_removed, Day7_2019_LMP_bad_removed=Day7_2019_LMP_bad_removed)

# Make new column and perform arcsine
Day7_2019_VIA_bad_removed$Percent_of_this_plot_arcsine <- transf.arcsin(Day7_2019_VIA_bad_removed$Percent_of_this_plot)
Day7_2019_APOP_bad_removed$Percent_of_this_plot_arcsine <- transf.arcsin(Day7_2019_APOP_bad_removed$Percent_of_this_plot)
Day7_2019_CASP_bad_removed$Percent_of_this_plot_arcsine <- transf.arcsin(Day7_2019_CASP_bad_removed$Percent_of_this_plot)
Day7_2019_LMP_bad_removed$Percent_of_this_plot_arcsine <- transf.arcsin(Day7_2019_LMP_bad_removed$Percent_of_this_plot)

# Change percent column to percent
Day7_2019_VIA_bad_removed$Percent_of_this_plot <- Day7_2019_VIA_bad_removed$Percent_of_this_plot*100
Day7_2019_APOP_bad_removed$Percent_of_this_plot<- Day7_2019_APOP_bad_removed$Percent_of_this_plot*100
Day7_2019_CASP_bad_removed$Percent_of_this_plot<- Day7_2019_CASP_bad_removed$Percent_of_this_plot*100
Day7_2019_LMP_bad_removed$Percent_of_this_plot <- Day7_2019_LMP_bad_removed$Percent_of_this_plot *100

#### PLOTS AND STATISTICS ####

### VIABILITY ASSAY Statistics ####

# Calculate summary statistics
# make sure group_by not masked by plyr
detach(package:Rmisc)
detach(package:plyr)
Day7_2019_VIA_bad_removed_summary <- 
  Day7_2019_VIA_bad_removed %>% group_by(Treat, Plot_number,Gate) %>% summarize(mean_percent=mean(Percent_of_this_plot), sd_percent =sd(Percent_of_this_plot))

unique(Day7_2019_VIA_bad_removed$Cell_type)
#all_hemocytes      granular           agranular          live_granular      live_agranular     all_live_hemocytes all_dead_hemocytes

# Percent Granular Hemocytes VIA all samples
Day7_2019_VIA_Percent_Granular <- Day7_2019_VIA_bad_removed %>% filter(Gate == "E1")
Day7_2019_VIA_Percent_Granular_Hemocytes <- ggplot(data=Day7_2019_VIA_Percent_Granular,
                              aes(y=Percent_of_this_plot, x=ID, color=Treat)) + 
                              geom_point() + 
                              ggtitle("Percent of Granular Hemocyte Events") + 
                              xlab("ID") +
                              ylab("Percent Granular Hemocytes") + ylim(0,100)

Day7_2019_VIA_Percent_Granular_summary <- Day7_2019_VIA_bad_removed_summary %>% filter(Gate =="E1")

# Percent Granular hemocytes summary (just means + sd) Bar plot
Day7_2019_VIA_Percent_Granular_Hemocytes_summary <- ggplot(data=Day7_2019_VIA_Percent_Granular_summary,
                                                    aes(y=mean_percent, x=Treat, fill=Treat)) + geom_col()+  xlab("ID") +
                                                    ylab("Percent Granular Hemocytes") + 
                                                    ylim(0,100) + ggtitle("Mean Percent of Granular Hemocyte Events") + 
                                                    geom_errorbar(aes(ymin=mean_percent-sd_percent,
                                                    ymax=mean_percent+sd_percent), width=0.2, position=position_dodge(0.9)) +
                                                    theme(panel.background=element_blank(),
                                                    panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
                                                    text=element_text(family="serif",size=16), 
                                                    axis.title.y=element_text(family="serif",size=16),
                                                    axis.title.x=element_text(family="serif",size=16),
                                                    legend.key=element_rect(fill=NA)) + 
                                                    theme(text=element_text(size=16)) + 
                                                    theme(axis.text.x = element_text(size=16)) +
                                                    theme(legend.text = element_text(size=16))

# Percent Agranular Hemocytes VIA All samples 
Day7_2019_VIA_Percent_Agranular <- Day7_2019_VIA_bad_removed %>% filter(Gate=="E3")
Day7_2019_VIA_Percent_Agranular_Hemocytes <- ggplot(data=Day7_2019_VIA_Percent_Agranular,
                                                    aes(y=Percent_of_this_plot, x=ID, color=Treat)) + 
                                                    geom_point() + 
                                                    ggtitle("Percent of Agranular Hemocyte Events") + 
                                                    xlab("ID") +
                                                    ylab("Percent Agranular Hemocytes") + ylim(0,100)



# Agranular and Granular cells Boxplot with significance bars, grouped by treatment color by gate(significant based on ANOVA)
Day7_2019_VIA_Percent_Agranular_Granular <- Day7_2019_VIA_bad_removed %>% filter(Gate =="E3" | Gate=="E1")
Day7_2019_VIA_Percent_Agranular_Granular_Hemocytes <- ggplot(data=Day7_2019_VIA_Percent_Agranular_Granular,
                                                    aes(y=Percent_of_this_plot, x=Treat, fill=Gate)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75))+ xlab("Treatment") +
                                                    ylab("Percent Hemocytes") + 
                                                      ggtitle("Percent of Granular and Agranular Hemocyte Events") + 
                                                    theme(panel.background=element_blank(),
                                                          panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
                                                          text=element_text(family="serif",size=12), 
                                                          axis.title.y=element_text(family="serif",size=12),
                                                          axis.title.x=element_text(family="serif",size=12),
                                                          legend.key=element_rect(fill=NA)) + 
                                                    theme(text=element_text(size=12)) + 
                                                    theme(axis.text.x = element_text(size=12)) +
                                                    theme(legend.text = element_text(size=12)) + 
                                                    scale_x_discrete(labels=c("D"="Dermo Injected", "NC"="Notched Control")) + 
                                                     scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
                                                    scale_fill_manual(name="Cell Type", labels=c("Granular","Agranular"), values=c("#6c81d9","#50b47b")) +
                                                    geom_signif(y_position = c(95,95), xmin = c(0.7,1.7), 
                                                                xmax = c(1.3,2.3), annotation = c("***","***"),
                                                                tip_length = 0.01)

# Agranular and Granular cells Boxplot with significance bars, grouped by GATE colored by treatment (no significance here based on ANOVA)
Day7_2019_VIA_Percent_Agranular_Granular_Hemocytes_cell_type <- ggplot(data=Day7_2019_VIA_Percent_Agranular_Granular,
                                   aes(y=Percent_of_this_plot, x=Gate, fill=Treat)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Treatment") +
  ylab("Percent Agranular and Granular Hemocytes") + 
   ggtitle("Percent of Granular and Agranular Hemocyte Events") +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) + 
  scale_x_discrete(labels=c("E1"="Granular", "E3"="Agranular")) + 
  scale_fill_manual(name="Treatment", labels=c("Dermo Injected","Notched Control"), values=c("#7e78d4",
                                                                                             "#cd4272")) +
  geom_signif(y_position = c(95,95), xmin = c(0.7,1.7), 
              xmax = c(1.3,2.3), annotation = c("NS","NS"),
              tip_length = 0.01)

## Calculate percent of live cells for granular and agranular within each treatment
# Cell type= live_granular and live_agranular 
#Day7_2019_VIA_Percent_Agranular_Granular %>% filter group_by(Treat,Assay,Gate)

# VIA Granular Agranular ANOVA with arcsine transformed
# one way anova of cell type within notched control
Day7_2019_VIA_bad_removed_granular_agranular_NC <- Day7_2019_VIA_bad_removed %>% filter(Gate == "E1" | Gate == "E3") %>% filter(Treat =="NC")
Day7_2019_VIA_bad_removed_cell_NC_aov <- aov(Percent_of_this_plot_arcsine ~ Gate, data= Day7_2019_VIA_bad_removed_granular_agranular_NC)
summary(Day7_2019_VIA_bad_removed_cell_NC_aov)

# Df Sum Sq Mean Sq F value   Pr(>F)    
#  Gate         1  1.728  1.7280   68.16 7.14e-08 ***
#  Residuals   20  0.507  0.0254                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# one way anove of cell type within dermo injected
Day7_2019_VIA_bad_removed_granular_agranular_D <- Day7_2019_VIA_bad_removed %>% filter(Gate == "E1" | Gate == "E3") %>% filter(Treat =="D")
Day7_2019_VIA_bad_removed_cell_D_aov <- aov(Percent_of_this_plot_arcsine ~ Gate, data= Day7_2019_VIA_bad_removed_granular_agranular_D)
summary(Day7_2019_VIA_bad_removed_cell_D_aov)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#Gate         1 0.6536  0.6536   88.33 6.47e-08 ***
#  Residuals   16 0.1184  0.0074                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# one way anova granular hemocytes between treatments
Day7_2019_VIA_Percent_Granular_aov <- aov(Percent_of_this_plot_arcsine ~ Treat , data=Day7_2019_VIA_Percent_Granular)
summary(Day7_2019_VIA_Percent_Granular_aov) 

# Df  Sum Sq Mean Sq F value Pr(>F)
# Treat        1 0.04085 0.04085   2.337  0.144
# Residuals   18 0.31469 0.01748 

# one way anova agranular hemocytes between treatments
Day7_2019_VIA_Percent_Agranular_aov <-  aov(Percent_of_this_plot_arcsine ~ Treat , data=Day7_2019_VIA_Percent_Agranular)
summary(Day7_2019_VIA_Percent_Agranular_aov) 

# Df  Sum Sq Mean Sq F value Pr(>F)
# Treat        1 0.03882 0.03882   2.249  0.151
# Residuals   18 0.31073 0.01726

## Summary
# Significantly more agranular than granular in both groups and treatment doesn't correlate with a difference in number of granular or agranular

##### Apoptosis Assay Statistics and Plotting ####

Day7_2019_APOP_bad_removed_summary <- 
  Day7_2019_APOP_bad_removed %>% group_by(Treat, Plot_number,Gate) %>% summarize(mean_percent=mean(Percent_of_this_plot), sd_percent =sd(Percent_of_this_plot))

# Q3 = agranular, Q8= granular

# Analysis for all Granular and Agranular Combined from the quad plots
# Filter out granular and agranular quad plot
Day7_2019_APOP_ALL_Granular_Agranular <- Day7_2019_APOP_bad_removed %>% filter(Gate == "Q8-UL" | Gate == "Q8-UR" | Gate == "Q8-LL" | Gate == "Q8-LR" | Gate == "Q3-UL" | Gate == "Q3-UR" | Gate == "Q3-LL" | Gate == "Q3-LR")

# Plot all gates barplot with ID on X axis
Day7_2019_APOP_ALL_Granular_Agranular_Hemocytes <- ggplot(data=Day7_2019_APOP_ALL_Granular_Agranular,
          aes(y=Percent_of_this_plot, x=ID, color=Gate)) + 
          geom_col(position = "dodge") + 
          ggtitle("Percent of Granular and Agranular Hemocyte Events") + 
          xlab("ID") + facet_grid(.~Treat, scales="free") +
          ylab("Percent Granular Hemocytes") + ylim(0,100)

# Make APOP combined gate for agranular and granular separately 
Day7_2019_APOP_ALL_Granular_Apop_combined <- Day7_2019_APOP_ALL_Granular_Agranular  %>% filter(Gate == "Q8-LR" | Gate == "Q8-UR") %>% group_by(ID, Treat) %>% summarise(Percent_of_this_plot = sum(Percent_of_this_plot))
Day7_2019_APOP_ALL_Agranular_Apop_combined <- Day7_2019_APOP_ALL_Granular_Agranular  %>% filter(Gate == "Q3-LR" | Gate == "Q3-UR") %>% group_by(ID, Treat) %>% summarise(Percent_of_this_plot = sum(Percent_of_this_plot))

# Add new gate name for each 
Day7_2019_APOP_ALL_Granular_Apop_combined$Gate <- "apop_combined_granular"
Day7_2019_APOP_ALL_Agranular_Apop_combined$Gate <- "apop_combined_agranular"

# Combined data frames for each cell type
Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined <- rbind(Day7_2019_APOP_ALL_Granular_Apop_combined, Day7_2019_APOP_ALL_Agranular_Apop_combined)

# Add arcsine transformed data
Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined$Percent_of_this_plot_arcsine <- transf.arcsin(Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined$Percent_of_this_plot*0.01)

# Merge with original data frame
Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined <- full_join(Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined,  Day7_2019_APOP_ALL_Granular_Agranular , by =c("Gate", "ID","Treat", "Percent_of_this_plot","Percent_of_this_plot_arcsine"))

# Calculate summary statistics that now include the apop_combined
Day7_2019_APOP_Percent_Granular_Agranular_summary <- Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined %>% group_by(Treat, Plot_number,Gate) %>% summarize(mean_percent=mean(Percent_of_this_plot), sd_percent =sd(Percent_of_this_plot)) 

# Boxplot of all gates of agranular and granular cells with faceting by gate, quadrants next to each other
# Change levels of faceting so that box plot quadrants are next to each other with agranular first, starting with live and dead
unique(Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined$Gate) # NOTE  Q3 = agranular, Q8= granular
          #"apop_combined_granular"  "apop_combined_agranular" "Q3-UL"                   "Q3-UL"                  "Q3-LL"                  
          #"Q3-LR"                   "Q8-UL"                   "Q8-UR"                   "Q8-LL"                   "Q8-LR"                  
 
Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined$Gate <- factor(Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined$Gate, 
        levels = c("Q3-UL", "Q8-UL","Q3-LL","Q8-LL","Q3-LR", "Q8-LR","Q3-UR", "Q8-UR","apop_combined_agranular","apop_combined_granular"))
# Make plot 
Day7_2019_APOP_Percent_Granular_Agranular_Hemocytes <- ggplot(data=Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined,
  aes(y=Percent_of_this_plot, x=Treat, fill=Gate)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Agranular and Granular FITC vs. PI Quad Plot") + 
  facet_grid(.~Gate, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  scale_fill_manual(name="Cell Type", labels=c( "Dead Agranular","Dead Granular", "Live Agranular", "Live Granular",
               "Live Agranular Apoptotic",   "Live Granular Apoptotic",  "Dead Agranular Apoptotic",  "Dead Granular Apoptotic", "Combined Agranular Apoptotic","Combined Granular Apoptotic"), 
              values = c("#cc57b4", "#88bf3b", "#aa4dce", "#56b464", "#5b2c90", 
                         "#c89832", "#5a6ee6", "#ca4e33", "#7e78d4", "#cd4272")) 

# Same boxplot as above by facet by gate and not by treatment 
Day7_2019_APOP_Percent_Granular_Agranular_Hemocytes <- ggplot(data=Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined,
  aes(y=Percent_of_this_plot, x=Gate, fill=Treat)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Quadrant") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Agranular and Granular FITC vs. PI Quad Plot") + 
  facet_grid(.~Gate, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(legend.text = element_text(size=12)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_manual(name="Treatment", labels=c("Dermo Injected","Notched Control"), values=c("#7e78d4","#cd4272")) +
  scale_x_discrete(labels=c("Q3-UL"= "Dead Agranular",
                            "Q8-UL"="Dead Granular", 
                            "Q3-LL"="Live Agranular", 
                            "Q8-LL"="Live Granular",
                            "Q3-LR"= "Live Agranular Apoptotic",   
                            "Q8-LR"="Live Granular Apoptotic",  
                            "Q3-UR"= "Dead Agranular Apoptotic",  
                            "Q8-UR"="Dead Granular Apoptotic", 
                            "apop_combined_agranular"="Combined Agranular Apoptotic",
                            "apop_combined_granular"="Combined Granular Apoptotic"))

# BoxPlot of just combined granular apoptotic and agranular combined apoptotic
Day7_2019_APOP_Percent_combined_apop <- Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined %>% filter(Gate== "apop_combined_agranular"| Gate == "apop_combined_granular") 
# change factor order of treatment so that the granular is first
Day7_2019_APOP_Percent_combined_apop$Treat <- factor(Day7_2019_APOP_Percent_combined_apop$Treat, levels=c("NC", "D"))
Day7_2019_APOP_ALL_Agranular_Granular_D_vs_NC_plot <-ggplot(data=Day7_2019_APOP_Percent_combined_apop,
  aes(y=Percent_of_this_plot, x=Gate, fill=Treat)) + 
  geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Cell Type") +
  ylab("Percent Hemocytes") +
  ggtitle("Percent Agranular and Granular Hemocyte Apoptosis") +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=16), 
        axis.title.y=element_text(family="serif",size=16),
        axis.title.x=element_text(family="serif",size=16),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=16))  +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100))+
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text = element_text(size=16)) +
  scale_x_discrete(labels=c("apop_combined_granular"="Granular Apoptotic", "apop_combined_agranular"="Agranular Apoptotic")) +
  scale_fill_manual(name="Cell Type", labels=c("Notched Control", "Dermo Injected"), values = c( "#6778d0","#ba496b")) + geom_signif(y_position = c(95,95), xmin = c(0.7,1.7), 
              xmax = c(1.3,2.3), annotation = c("**","***"),
              tip_length = 0.01)

# export plot 
ggsave(plot = Day7_2019_APOP_ALL_Agranular_Granular_D_vs_NC_plot, path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       filename = "Day7_2019_APOP_ALL_Agranular_Granular_D_vs_NC_plot.tiff", device = "tiff", units = "cm",
       width = 20, height = 15)

# BoxPlot of just combined granular apoptotic
Day7_2019_APOP_Percent_combined_apop_granular <- Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined %>% filter(Gate == "apop_combined_granular") 
Day7_2019_APOP_ALL_Granular_D_vs_NC_plot <-ggplot(data=Day7_2019_APOP_Percent_combined_apop_granular,
                                                            aes(y=Percent_of_this_plot, x=Gate, fill=Treat)) + 
  geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Cell Type") +
  ylab("Percent Hemocytes") +
  ggtitle("Percent Granular FITC vs. PI Quad Plot") +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=16), 
        axis.title.y=element_text(family="serif",size=16),
        axis.title.x=element_text(family="serif",size=16),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=16))  +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100))+
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text = element_text(size=16)) +
  scale_x_discrete(labels=c("apop_combined_granular"="Granular Apoptotic")) +
  scale_fill_manual(name="Cell Type", labels=c("Dermo Injected","Notched Control"), values = c("#ba496b","#6778d0")) + geom_signif(y_position = 95, xmin = 0.7, 
                                                                                                                                    xmax = 1.3, annotation = "***",
                                                                                                                                    tip_length = 0.01)


# ANOVA
# Apop combined granular vs apop combined agranular within each treatment
# Notched control
Day7_2019_APOP_ALL_Granular_Agranular_Apop_NC <- Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined %>% filter(Gate =="apop_combined_agranular" | Gate== "apop_combined_granular") %>% filter(Treat == "NC")
Day7_2019_APOP_ALL_Granular_Agranular_Apop_NC_aov <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_APOP_ALL_Granular_Agranular_Apop_NC)                                                  
summary(Day7_2019_APOP_ALL_Granular_Agranular_Apop_NC_aov)

      #Df Sum Sq Mean Sq F value   Pr(>F)    
      #Gate         1 0.8553  0.8553   39.14 4.14e-06 ***
      #  Residuals   20 0.4371  0.0219                     
      #---
      #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Dermo
Day7_2019_APOP_ALL_Granular_Agranular_Apop_D <- Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined %>% filter(Gate =="apop_combined_agranular" | Gate== "apop_combined_granular") %>% filter(Treat == "D")
Day7_2019_APOP_ALL_Granular_Agranular_Apop_D_aov <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_APOP_ALL_Granular_Agranular_Apop_D)                                                  
summary(Day7_2019_APOP_ALL_Granular_Agranular_Apop_D_aov)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#Gate         1 0.1939 0.19393   29.76 5.29e-05 ***
#Residuals   16 0.1042 0.00652 

# Apop combined granular between treatments
Day7_2019_APOP_ALL_Granular_D_vs_NC <- Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined %>% filter(Gate== "apop_combined_granular") 
Day7_2019_APOP_ALL_Granular_D_vs_NC_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_APOP_ALL_Granular_D_vs_NC)                                                  
summary(Day7_2019_APOP_ALL_Granular_D_vs_NC_aov)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#Treat        1 0.4499  0.4499   18.17 0.000468 ***
# Residuals   18 0.4456  0.0248 

# Apop combined agranular between treatments 
Day7_2019_APOP_ALL_Agranular_D_vs_NC <- Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined %>% filter(Gate== "apop_combined_agranular") 
Day7_2019_APOP_ALL_Agranular_D_vs_NC_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_APOP_ALL_Agranular_D_vs_NC)                                                  
summary(Day7_2019_APOP_ALL_Agranular_D_vs_NC_aov)

# Df  Sum Sq Mean Sq F value  Pr(>F)   
# Treat        1 0.06514 0.06514   12.25 0.00255 **
# Residuals   18 0.09570 0.00532

## SUMMARY
  # Significantly increased granulocytes within treatments, but between treatments notched control had significantly greater agranular and granular hemocyte apoptosis
  # this is indicative of reduced apoptotic immune response following challenge or parasite inhibition of apoptosis

#### Caspase Assay Statistics ####

Day7_2019_CASP_bad_removed
unique(Day7_2019_CASP_bad_removed[,c(2:5)])

# Q2= agranular, Q5= granular

# Analysis for all Granular and Agranular Combined caspase apoptotic from the quad plots
# Filter out granular and agranular quad plot
Day7_2019_CASP_ALL_Granular_Agranular <- Day7_2019_CASP_bad_removed %>% filter(Gate == "Q2-UL" | Gate == "Q2-UR" | Gate == "Q2-LL" | Gate == "Q2-LR" | Gate == "Q5-UL" | Gate == "Q5-UR" | Gate == "Q5-LL" | Gate == "Q5-LR")

# Make casp_apop combined gate for agranular and granular separately 
Day7_2019_CASP_ALL_Granular_Apop_combined <-  Day7_2019_CASP_ALL_Granular_Agranular  %>% filter(Gate == "Q5-LR" | Gate == "Q5-UR") %>% group_by(ID, Treat) %>% summarise(Percent_of_this_plot = sum(Percent_of_this_plot))
Day7_2019_CASP_ALL_Agranular_Apop_combined <- Day7_2019_CASP_ALL_Granular_Agranular  %>% filter(Gate == "Q2-LR" | Gate == "Q2-UR") %>% group_by(ID, Treat) %>% summarise(Percent_of_this_plot = sum(Percent_of_this_plot))

# Add new gate name for each 
Day7_2019_CASP_ALL_Granular_Apop_combined$Gate <- "casp_apop_combined_granular"
Day7_2019_CASP_ALL_Agranular_Apop_combined$Gate <- "casp_apop_combined_agranular"

# Combined data frames for each cell type
Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined <- rbind(Day7_2019_CASP_ALL_Granular_Apop_combined, Day7_2019_CASP_ALL_Agranular_Apop_combined)

# Add arcsine transformed data
Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined$Percent_of_this_plot_arcsine <- transf.arcsin(Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined$Percent_of_this_plot*0.01)

# Merge with original data frame
Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined <- full_join(Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined,  Day7_2019_CASP_ALL_Granular_Agranular , by =c("Gate", "ID","Treat", "Percent_of_this_plot","Percent_of_this_plot_arcsine"))

# Boxplot of all gates of agranular and granular cells with faceting by gate, quadrants next to each other
# Change levels of faceting so that box plot quadrants are next to each other with agranular first, starting with live and dead
unique(Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined$Gate) # NOTE  Q2 = agranular, Q5= granular
#[1] "casp_apop_combined_granular"  "casp_apop_combined_agranular" "Q2-UL"                        "Q2-UR"                       
#[5] "Q2-LL"                        "Q2-LR"                        "Q5-UL"                        "Q5-UR"                       
#[9] "Q5-LL"                        "Q5-LR"                       

Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined$Gate <- factor(Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined$Gate, 
    levels = c("Q2-UL", "Q5-UL","Q2-LL","Q5-LL","Q2-LR", "Q5-LR","Q2-UR", "Q5-UR","casp_apop_combined_agranular","casp_apop_combined_granular"))

# Same boxplot as above by facet by gate and not by treatment 
Day7_2019_CASP_Percent_Granular_Agranular_Hemocytes <- ggplot(data=Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined,
  aes(y=Percent_of_this_plot, x=Gate, fill=Treat)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Quadrant") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Agranular and Granular Caspase Apoptotic FITC vs. APC Quad Plot") + 
  facet_grid(.~Gate, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(legend.text = element_text(size=12)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_manual(name="Treatment", labels=c("Dermo Injected","Notched Control"), values=c("#7e78d4","#cd4272")) +
  scale_x_discrete(labels=c("Q2-UL"= "Dead Agranular",
                            "Q5-UL"="Dead Granular", 
                            "Q2-LL"="Live Agranular", 
                            "Q5-LL"="Live Granular",
                            "Q2-LR"= "Live Agranular Caspase 3/7 Apoptotic",   
                            "Q5-LR"="Live Granular Caspase 3/7 Apoptotic",  
                            "Q2-UR"= "Dead AgranularCaspase 3/7 Apoptotic",  
                            "Q5-UR"="Dead Granular Caspase 3/7 Apoptotic", 
                            "casp_apop_combined_agranular"="Combined Caspase 3/7 Agranular Apoptotic",
                            "casp_apop_combined_granular"="Combined Caspase 3/7 Granular Apoptotic"))

# BoxPlot of just combined granular and agranular apoptotic
Day7_2019_CASP_Granular_Agranular_Apop_combined <- Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined %>% filter(Gate == "casp_apop_combined_granular" | Gate == "casp_apop_combined_agranular") 
# change facet so that control is first on the plot
Day7_2019_CASP_Granular_Agranular_Apop_combined$Treat <- factor(Day7_2019_CASP_Granular_Agranular_Apop_combined$Treat, levels=c("NC","D"))
Day7_2019_CASP_Granular_D_vs_NC_plot <-ggplot(data=Day7_2019_CASP_Granular_Agranular_Apop_combined,
                                                  aes(y=Percent_of_this_plot, x=Gate, fill=Treat)) + 
  geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Cell Type") +
  ylab("Percent Hemocytes") +
  ggtitle("Percent Caspase 3/7 Active Hemocytes") +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=16), 
        axis.title.y=element_text(family="serif",size=16),
        axis.title.x=element_text(family="serif",size=16),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=16))  +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100))+
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text = element_text(size=16)) +
  scale_x_discrete(labels=c("casp_apop_combined_granular"="Granular", "casp_apop_combined_agranular"="Agranular")) +
  scale_fill_manual(name="Cell Type", labels=c("Notched Control", "Dermo Injected"), values = c("#e08c67","#6388ca")) 
                                                                                                                                   
ggsave(plot = Day7_2019_CASP_Granular_D_vs_NC_plot, path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       filename = "Day7_2019_CASP_Granular_D_vs_NC_plot.tiff", device = "tiff", units = "cm",
       width = 20, height = 15)

# ANOVA
# Casp active combined granular vs apop combined agranular within each treatment
# Notched control
Day7_2019_CASP_ALL_Granular_Agranular_Apop_NC <- Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined %>% filter(Gate =="casp_apop_combined_agranular" | Gate== "casp_apop_combined_granular") %>% filter(Treat == "NC")
Day7_2019_CASP_ALL_Granular_Agranular_Apop_NC_aov <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_CASP_ALL_Granular_Agranular_Apop_NC)                                                  
summary(Day7_2019_CASP_ALL_Granular_Agranular_Apop_NC_aov)

            #Df Sum Sq Mean Sq F value   Pr(>F)    
            #Gate         1  3.258   3.258   205.6 5.52e-12 ***
            #  Residuals   20  0.317   0.016                     
            #---
            #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Dermo
Day7_2019_CASP_ALL_Granular_Agranular_Apop_D <- Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined %>% filter(Gate =="casp_apop_combined_agranular" | Gate== "casp_apop_combined_granular") %>% filter(Treat == "D")
Day7_2019_CASP_ALL_Granular_Agranular_Apop_D_aov <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_CASP_ALL_Granular_Agranular_Apop_D)                                                  
summary(Day7_2019_CASP_ALL_Granular_Agranular_Apop_D_aov)

            #Df Sum Sq Mean Sq F value   Pr(>F)    
            #Gate         1 2.5285  2.5285   80.56 1.21e-07 ***
            #  Residuals   16 0.5022  0.0314                     
            #---
            #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Apop combined granular between treatments
Day7_2019_CASP_ALL_Granular_D_vs_NC <- Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined %>% filter(Gate== "casp_apop_combined_granular") 
Day7_2019_CASP_ALL_Granular_D_vs_NC_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_CASP_ALL_Granular_D_vs_NC)                                                  
summary(Day7_2019_CASP_ALL_Granular_D_vs_NC_aov)

            #Df Sum Sq  Mean Sq F value Pr(>F)
            #Treat        1 0.0063 0.006297   0.374  0.549
            #Residuals   18 0.3032 0.016846

# Apop combined agranular between treatments 
Day7_2019_CASP_ALL_Agranular_D_vs_NC <- Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined %>% filter(Gate== "casp_apop_combined_agranular") 
Day7_2019_CASP_ALL_Agranular_D_vs_NC_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_CASP_ALL_Agranular_D_vs_NC)                                                  
summary(Day7_2019_CASP_ALL_Agranular_D_vs_NC_aov)

            #Df Sum Sq Mean Sq F value Pr(>F)
            #Treat        1 0.0154 0.01538   0.537  0.473
            #Residuals   18 0.5160 0.02866 

## SUMMARY
# Significantly increased granulocytes caspase 3/7 activation, but no difference between treatments

#### LMP Assay Statistics ####

Day7_2019_LMP_bad_removed
unique(Day7_2019_LMP_bad_removed[,c(2:5)])

# Q3 = agranular, Q2 = granular

# Analysis for all Granular and Agranular Combined caspase apoptotic from the quad plots
# Filter out granular and agranular quad plot
Day7_2019_LMP_ALL_Granular_Agranular <- Day7_2019_LMP_bad_removed %>% filter(Gate == "Q3-UL" | Gate == "Q3-UR" | Gate == "Q3-LL" | Gate == "Q3-LR" | Gate == "Q2-UL" | Gate == "Q2-UR" | Gate == "Q2-LL" | Gate == "Q2-LR")

# Boxplot of all gates of agranular and granular cells with faceting by gate, quadrants next to each other
# Change levels of faceting so that box plot quadrants are next to each other with agranular first, starting with live and dead

unique(Day7_2019_LMP_ALL_Granular_Agranular$Gate) # NOTE Q3 = agranular, Q2 = granular
# "Q3-UL" "Q3-UR" "Q3-LL" "Q3-LR" "Q2-UL" "Q2-UR" "Q2-LL" "Q2-LR"                     

Day7_2019_LMP_ALL_Granular_Agranular$Gate <- factor(Day7_2019_LMP_ALL_Granular_Agranular$Gate, 
                                                                   levels = c("Q3-UL","Q2-UL", "Q3-LL","Q2-LL","Q3-UR", "Q2-UR", "Q3-LR","Q2-LR"))

# Boxplot as above by facet by gate and not by treatment 
Day7_2019_LMP_Percent_Granular_Agranular_Hemocytes <- ggplot(data=Day7_2019_LMP_ALL_Granular_Agranular,
                                                              aes(y=Percent_of_this_plot, x=Gate, fill=Treat)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Quadrant") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Agranular and Granular Intact Lysosome FITC vs. PI Quad Plot") + 
  facet_grid(.~Gate, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, hjust=1, angle=70)) +
  theme(legend.text = element_text(size=12)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_manual(name="Treatment", labels=c("Dermo Injected","Notched Control"), values=c("#7e78d4","#cd4272")) +
  scale_x_discrete(labels=c("Q3-UL"="agranular_lysosome_rupture",
                            "Q2-UL"="granular_lysosome_rupture", 
                            "Q3-LL"="agranular_LMP_onset",
                            "Q2-LL"="granular_LMP_onset", 
                            "Q3-UR"="agranular_dead_other_means", 
                            "Q2-UR"="granular_dead_other_means",
                            "Q3-LR"="agranular_lysosome_live",
                            "Q2-LR"="granular_lysosome_live"))

  

# ANOVA

## LMP RUPTURE ANOVAS
# LMP rupture granular vs LMP rupture granular within each treatment
# Notched control
Day7_2019_LMP_ALL_Granular_Agranular_rupture_NC <- Day7_2019_LMP_ALL_Granular_Agranular %>% filter(Gate =="Q3-UL" | Gate== "Q2-UL") %>% filter(Treat == "NC")
Day7_2019_LMP_ALL_Granular_Agranular_rupture_NC_aov <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_LMP_ALL_Granular_Agranular_rupture_NC)                                                  
summary(Day7_2019_LMP_ALL_Granular_Agranular_rupture_NC_aov)

        #Df  Sum Sq Mean Sq F value Pr(>F)  
        #Gate         1 0.07806 0.07806   6.252 0.0212 *
        #  Residuals   20 0.24972 0.01249                 
        #---
        #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Dermo
Day7_2019_LMP_ALL_Granular_Agranular_rupture_D <- Day7_2019_LMP_ALL_Granular_Agranular %>% filter(Gate =="Q3-UL" | Gate== "Q2-UL") %>% filter(Treat == "D")
Day7_2019_LMP_ALL_Granular_Agranular_rupture_D_aov <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_LMP_ALL_Granular_Agranular_rupture_D)                                                  
summary(Day7_2019_LMP_ALL_Granular_Agranular_rupture_D_aov)

          #Df  Sum Sq Mean Sq F value Pr(>F)  
          #Gate         1 0.06203 0.06203       8 0.0121 *
          #  Residuals   16 0.12407 0.00775                 
          #---
          #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# LMP rupture granular between treatments
# Q3 = agranular, Q2 = granular
Day7_2019_LMP_ALL_Granular_rupture <- Day7_2019_LMP_ALL_Granular_Agranular %>% filter(Gate== "Q2-UL") 
Day7_2019_LMP_ALL_Granular_rupture_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_LMP_ALL_Granular_rupture )                                                  
summary(Day7_2019_LMP_ALL_Granular_rupture_aov)

          #Df  Sum Sq Mean Sq F value Pr(>F)
          #Treat        1 0.01218 0.01218   0.946  0.344
          #Residuals   18 0.23185 0.01288  

# LMP rupture agranular between treatments
Day7_2019_LMP_ALL_Agranular_rupture <- Day7_2019_LMP_ALL_Granular_Agranular %>% filter(Gate== "Q3-UL") 
Day7_2019_LMP_ALL_Agranular_rupture_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_LMP_ALL_Agranular_rupture )                                                  
summary(Day7_2019_LMP_ALL_Agranular_rupture_aov)

          #Df Sum Sq  Mean Sq F value Pr(>F)
          #Treat        1 0.0115 0.011500   1.441  0.246
          #Residuals   18 0.1437 0.007981   

## LMP ONSET ANOVAs
# LMP onset granular vs LMP rupture granular within each treatment
# Notched control
Day7_2019_LMP_ALL_Granular_Agranular_onset_NC <- Day7_2019_LMP_ALL_Granular_Agranular %>% filter(Gate =="Q3-LL" | Gate== "Q2-LL") %>% filter(Treat == "NC")
Day7_2019_LMP_ALL_Granular_Agranular_onset_NC_aov <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_LMP_ALL_Granular_Agranular_onset_NC )                                                  
summary(Day7_2019_LMP_ALL_Granular_Agranular_onset_NC_aov)

          #Df Sum Sq Mean Sq F value   Pr(>F)    
          #Gate         1 2.8191  2.8191   333.2 6.12e-14 ***
          #  Residuals   20 0.1692  0.0085                     
          #---
          #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Dermo
Day7_2019_LMP_ALL_Granular_Agranular_onset_D <- Day7_2019_LMP_ALL_Granular_Agranular %>% filter(Gate =="Q3-LL" | Gate== "Q2-LL") %>% filter(Treat == "D")
Day7_2019_LMP_ALL_Granular_Agranular_onset_D_aov <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_LMP_ALL_Granular_Agranular_onset_D)                                                  
summary(Day7_2019_LMP_ALL_Granular_Agranular_onset_D_aov)

          #Df Sum Sq Mean Sq F value   Pr(>F)    
          #Gate         1 1.7932  1.7932   200.7 1.79e-10 ***
          #  Residuals   16 0.1429  0.0089                     
          #---
          #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# LMP onset granular between treatments
# Q3 = agranular, Q2 = granular
Day7_2019_LMP_ALL_Granular_onset <- Day7_2019_LMP_ALL_Granular_Agranular %>% filter(Gate== "Q2-LL") 
Day7_2019_LMP_ALL_Granular_onset_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_LMP_ALL_Granular_onset )                                                  
summary(Day7_2019_LMP_ALL_Granular_onset_aov)

          #Df  Sum Sq  Mean Sq F value Pr(>F)
          #Treat        1 0.00003 0.000032   0.009  0.926
          #Residuals   18 0.06448 0.003582  

# LMP onset agranular between treatments
Day7_2019_LMP_ALL_Agranular_onset <- Day7_2019_LMP_ALL_Granular_Agranular %>% filter(Gate== "Q3-LL") 
Day7_2019_LMP_ALL_Agranular_onset_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_LMP_ALL_Agranular_onset )                                                  
summary(Day7_2019_LMP_ALL_Agranular_onset_aov)

          #Df  Sum Sq Mean Sq F value Pr(>F)
          #Treat        1 0.03746 0.03746    2.72  0.116
          #Residuals   18 0.24796 0.01378  

### LMP Intact lysosome live 
# LMP intact granular vs LMP rupture granular within each treatment
# Notched control
Day7_2019_LMP_ALL_Granular_Agranular_intact_NC <- Day7_2019_LMP_ALL_Granular_Agranular %>% filter(Gate =="Q3-LR" | Gate== "Q2-LR") %>% filter(Treat == "NC")
Day7_2019_LMP_ALL_Granular_Agranular_intact_NC_aov <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_LMP_ALL_Granular_Agranular_intact_NC )                                                  
summary(Day7_2019_LMP_ALL_Granular_Agranular_intact_NC_aov)

          #Df Sum Sq Mean Sq F value Pr(>F)  
          #Gate         1 0.2235 0.22351   5.409 0.0307 *
          #  Residuals   20 0.8265 0.04133                 
          #---
          #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Dermo
Day7_2019_LMP_ALL_Granular_Agranular_intact_D <- Day7_2019_LMP_ALL_Granular_Agranular %>% filter(Gate =="Q3-LR" | Gate== "Q2-LR") %>% filter(Treat == "D")
Day7_2019_LMP_ALL_Granular_Agranular_intact_D_aov <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_LMP_ALL_Granular_Agranular_intact_D)                                                  
summary(Day7_2019_LMP_ALL_Granular_Agranular_intact_D_aov)

          #Df Sum Sq Mean Sq F value Pr(>F)  
          #Gate         1 0.1669 0.16689   6.614 0.0205 *
          # Residuals   16 0.4038 0.02523                 
          # ---
          #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# LMP intact granular between treatments
# Q3 = agranular, Q2 = granular
Day7_2019_LMP_ALL_Granular_intact <- Day7_2019_LMP_ALL_Granular_Agranular %>% filter(Gate== "Q2-LR") 
Day7_2019_LMP_ALL_Granular_intact_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_LMP_ALL_Granular_intact )                                                  
summary(Day7_2019_LMP_ALL_Granular_intact_aov)

          #Df Sum Sq Mean Sq F value Pr(>F)
          #Treat        1 0.0736 0.07361   1.469  0.241
          #Residuals   18 0.9022 0.05012   

# LMP intact agranular between treatments
Day7_2019_LMP_ALL_Agranular_intact <- Day7_2019_LMP_ALL_Granular_Agranular %>% filter(Gate== "Q3-LR") 
Day7_2019_LMP_ALL_Agranular_intact_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_LMP_ALL_Agranular_intact )                                                  
summary(Day7_2019_LMP_ALL_Agranular_intact_aov)

          #Df Sum Sq Mean Sq F value Pr(>F)  
          #Treat        1 0.0849 0.08489   4.658 0.0447 *
          #  Residuals   18 0.3281 0.01823                 
          #---
          #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##### PCA of all assays ####
library(ggbiplot)
library(factoextra)

# Create combined data frame with all assays

# check number of columns in each
ncol(Day7_2019_VIA_Percent_Agranular_Granular)
ncol(Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined)
ncol(Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined)
ncol(Day7_2019_LMP_ALL_Granular_Agranular)

# check number of gates in each, figuring out how to eliminate duplicates between assays
unique(Day7_2019_VIA_Percent_Agranular_Granular$Gate) # [1] "E1" "E3"
unique(Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined$Gate) # [1] apop_combined_granular  apop_combined_agranular Q3-UL                  
          # [4] Q3-UR                   Q3-LL                   Q3-LR                  
          # [7] Q8-UL                   Q8-UR                   Q8-LL                  
          # [10] Q8-LR
unique(Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined$Gate) 
          # [1] casp_apop_combined_granular  casp_apop_combined_agranular Q2-UL                       
          # [4] Q2-UR                        Q2-LL                        Q2-LR                       
          # [7] Q5-UL                        Q5-UR                        Q5-LL                       
          # [10] Q5-LR 
unique(Day7_2019_LMP_ALL_Granular_Agranular$Gate)
     
# All have the same columns but all are in a different order
Day_7_All_Assays_Merged <- full_join(Day7_2019_VIA_Percent_Agranular_Granular, Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined, by = c("ID","Treat","Assay",
                                                                                                                "Gate","Plot_number", "Cell_type", "Counts", "Percent_of_this_plot", "Percent_of_this_plot_arcsine"))
Day_7_All_Assays_Merged <- full_join(Day_7_All_Assays_Merged, Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined, by = c("ID","Treat","Assay",
                                                                                                                          "Gate","Plot_number", "Cell_type","Counts","Percent_of_this_plot", "Percent_of_this_plot_arcsine"))
Day_7_All_Assays_Merged <- full_join(Day_7_All_Assays_Merged, Day7_2019_LMP_ALL_Granular_Agranular, by = c("ID","Treat","Assay",
                                                                                                           "Gate","Plot_number", "Cell_type", "Counts", "Percent_of_this_plot", "Percent_of_this_plot_arcsine"))
# Goal: create data frame where each parameter for each sample is side by side
# Add Assay + Gate column to use for spreading
Day_7_All_Assays_Merged <- Day_7_All_Assays_Merged %>% unite(Assay_Gate, c("Assay","Gate"), remove=FALSE)

# reduce columns to facilitate spread using arcsine as the values
Day_7_All_Assays_Merged_reduced <- Day_7_All_Assays_Merged[,c("ID","Assay_Gate","Percent_of_this_plot_arcsine")]

# Spread data to use for PCA calculation
  # Some Gates are shared, so can't spread based on Gate alone, spread using combined Assay_Gate
Day_7_All_Assays_Merged_reduced_spread <- spread(data=Day_7_All_Assays_Merged_reduced, Assay_Gate, Percent_of_this_plot_arcsine)

# create key for samples in each treatment
ID_treat_key <- Day7_2019_VIA_Percent_Agranular_Granular[c(1:20),c(1:2)]
Day_7_All_Assays_Merged_reduced_spread <- full_join(Day_7_All_Assays_Merged_reduced_spread, ID_treat_key , by ="ID" )
colnames(Day_7_All_Assays_Merged_reduced_spread)
Day_7_All_Assays_Merged_reduced_spread_subset <- Day_7_All_Assays_Merged_reduced_spread[,c(1, 18:32)]

#calculate PCA using the dataframe before
All_assays_PCA <- prcomp(Day_7_All_Assays_Merged_reduced_spread_subset[2:15], center=TRUE, scale=TRUE)
summary(All_assays_PCA)

ggbiplot(All_assays_PCA, varname.adjust= 0.1, varname.abbrev = TRUE) + geom_text(vjust="inward",hjust="inward", label=Day_7_All_Assays_Merged_reduced_spread$ID)

#### OUTLIER REMOVAL FROM DAY 7 ####

# Outliers identified in plots to remove = Notched control oysters 244 and 275 

##### Remove bad samples from all
#SAMPLES TO REMOVE FROM EACH ANALYSIS = 265 (notched control), 193 (Dermo injected), 218 (Dermo injected), 221 (Dermo injected), 244 (Nothched control), 275 (notched control)

# vector with ID lines for each assay
ID_samples_remove_outliers <- c("265","193","218","221", "244","275")

# List of assay dataframes
assay_list <- list(Day7_2019_VIA_join=Day7_2019_VIA_join, Day7_2019_APOP_join=Day7_2019_APOP_join, 
                   Day7_2019_CASP_join=Day7_2019_CASP_join, Day7_2019_LMP_join=Day7_2019_LMP_join)

# Save lapply output after removing samples
remove_temp_outliers <- lapply(assay_list, function(x) x[ ! x$ID %in% ID_samples_remove_outliers, ])
# set names of each lapply output object in correct order and save using list2env
names(remove_temp_outliers) <- c("Day7_2019_VIA_bad_removed_no_outlier",  "Day7_2019_APOP_bad_removed_no_outlier", 
                        "Day7_2019_CASP_bad_removed_no_outlier", "Day7_2019_LMP_bad_removed_no_outlier")
list2env(remove_temp_outliers, envir = .GlobalEnv)

# check output
unique(Day7_2019_VIA_bad_removed_no_outlier$ID) # output is missing correct samples and has correct format 

#### Arcsine Transformation of percentages ####

# Make list of dataframes
bad_removed__no_outlier_assay_list <- list(Day7_2019_VIA_bad_removed_no_outlier=Day7_2019_VIA_bad_removed_no_outlier, Day7_2019_APOP_bad_removed_no_outlier=Day7_2019_APOP_bad_removed_no_outlier, 
                               Day7_2019_CASP_bad_removed_no_outlier=Day7_2019_CASP_bad_removed_no_outlier, Day7_2019_LMP_bad_removed_no_outlier=Day7_2019_LMP_bad_removed_no_outlier)

# Make new column and perform arcsine
Day7_2019_VIA_bad_removed_no_outlier$Percent_of_this_plot_arcsine <- transf.arcsin(Day7_2019_VIA_bad_removed_no_outlier$Percent_of_this_plot)
Day7_2019_APOP_bad_removed_no_outlier$Percent_of_this_plot_arcsine <- transf.arcsin(Day7_2019_APOP_bad_removed_no_outlier$Percent_of_this_plot)
Day7_2019_CASP_bad_removed_no_outlier$Percent_of_this_plot_arcsine <- transf.arcsin(Day7_2019_CASP_bad_removed_no_outlier$Percent_of_this_plot)
Day7_2019_LMP_bad_removed_no_outlier$Percent_of_this_plot_arcsine <- transf.arcsin(Day7_2019_LMP_bad_removed_no_outlier$Percent_of_this_plot)

# Change percent column to percent
 Day7_2019_VIA_bad_removed_no_outlier$Percent_of_this_plot <- Day7_2019_VIA_bad_removed_no_outlier$Percent_of_this_plot*100
Day7_2019_APOP_bad_removed_no_outlier$Percent_of_this_plot<- Day7_2019_APOP_bad_removed_no_outlier$Percent_of_this_plot*100
Day7_2019_CASP_bad_removed_no_outlier$Percent_of_this_plot<- Day7_2019_CASP_bad_removed_no_outlier$Percent_of_this_plot*100
 Day7_2019_LMP_bad_removed_no_outlier$Percent_of_this_plot <- Day7_2019_LMP_bad_removed_no_outlier$Percent_of_this_plot *100

#### PLOTS AND STATISTICS ####

### Viability Assay Statistics ####

# Calculate summary statistics
# make sure group_by not masked by plyr
#detach(package:Rmisc)
#detach(package:plyr)
unique(Day7_2019_VIA_bad_removed_no_outlier$Cell_type)
#all_hemocytes      granular           agranular          live_granular      live_agranular     all_live_hemocytes all_dead_hemocytes

Day7_2019_VIA_Percent_Granular_no_outlier <- Day7_2019_VIA_bad_removed_no_outlier %>% filter(Gate =="E1")
Day7_2019_VIA_Percent_Agranular_no_outlier <- Day7_2019_VIA_bad_removed_no_outlier %>% filter(Gate =="E3")

# Agranular and Granular cells Boxplot with significance bars, grouped by GATE colored by treatment (no significance here based on ANOVA)
Day7_2019_VIA_Percent_Agranular_Granular_Hemocytes_cell_type_no_outlier <- ggplot(data=Day7_2019_VIA_Percent_Agranular_Granular_no_outlier,
                                                                       aes(y=Percent_of_this_plot, x=Gate, fill=Treat)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Treatment") +
  ylab("Percent Agranular and Granular Hemocytes") + 
  ggtitle("Percent of Granular and Agranular Hemocyte Events") +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) + 
  scale_x_discrete(labels=c("E1"="Granular", "E3"="Agranular")) + 
  scale_fill_manual(name="Treatment", labels=c("Dermo Injected","Notched Control"), values=c("#7e78d4",
                                                                                             "#cd4272")) +
  geom_signif(y_position = c(95,95), xmin = c(0.7,1.7), 
              xmax = c(1.3,2.3), annotation = c("***","***"),
              tip_length = 0.01)

# VIA Granular Agranular ANOVA with arcsine transformed no outliers
# one way anova of cell type within notched control
Day7_2019_VIA_bad_removed_granular_agranular_NC_no_outlier <- Day7_2019_VIA_bad_removed_no_outlier %>% filter(Gate == "E1" | Gate == "E3") %>% filter(Treat =="NC")
Day7_2019_VIA_bad_removed_cell_NC_no_outlier_aov <- aov(Percent_of_this_plot_arcsine ~ Gate, data= Day7_2019_VIA_bad_removed_granular_agranular_NC_no_outlier)
summary(Day7_2019_VIA_bad_removed_cell_NC_no_outlier_aov)

      #Df Sum Sq Mean Sq F value   Pr(>F)    
      #Gate         1 2.1804  2.1804   680.8 1.53e-14 ***
      #  Residuals   16 0.0512  0.0032                     
      #---
      #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# one way anove of cell type within dermo injected
Day7_2019_VIA_bad_removed_granular_agranular_D_no_outlier <- Day7_2019_VIA_bad_removed_no_outlier %>% filter(Gate == "E1" | Gate == "E3") %>% filter(Treat =="D")
Day7_2019_VIA_bad_removed_cell_D_no_outlier_aov <- aov(Percent_of_this_plot_arcsine ~ Gate, data= Day7_2019_VIA_bad_removed_granular_agranular_D_no_outlier)
summary(Day7_2019_VIA_bad_removed_cell_D_no_outlier_aov)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#Gate         1 0.6536  0.6536   88.33 6.47e-08 ***
#  Residuals   16 0.1184  0.0074                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# one way anova granular hemocytes between treatments
Day7_2019_VIA_Percent_Granular_no_outlier_aov <- aov(Percent_of_this_plot_arcsine ~ Treat , data=Day7_2019_VIA_Percent_Granular_no_outlier)
summary(Day7_2019_VIA_Percent_Granular_no_outlier_aov) 

      #Df  Sum Sq Mean Sq F value   Pr(>F)    
      #Treat        1 0.11372 0.11372   21.54 0.000272 ***
      #  Residuals   16 0.08449 0.00528                     
      #---
      #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# one way anova agranular hemocytes between treatments
Day7_2019_VIA_Percent_Agranular_no_outlier_aov <-  aov(Percent_of_this_plot_arcsine ~ Treat , data=Day7_2019_VIA_Percent_Agranular_no_outlier)
summary(Day7_2019_VIA_Percent_Agranular_no_outlier_aov) 

      #Df  Sum Sq Mean Sq F value   Pr(>F)    
      #Treat        1 0.10953 0.10953   20.58 0.000337 ***
      #  Residuals   16 0.08515 0.00532                     
      #---
      #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Summary
# Significantly more agranular than granular in both groups and treatment does correlate with a difference in number of granular or agranular

##### Apoptosis Assay Statistics and Plotting ####

# Q3 = agranular, Q8= granular

# Analysis for all Granular and Agranular Combined from the quad plots
# Filter out granular and agranular quad plot
Day7_2019_APOP_ALL_Granular_Agranular_no_outlier <- Day7_2019_APOP_bad_removed_no_outlier %>% filter(Gate == "Q8-UL" | Gate == "Q8-UR" | Gate == "Q8-LL" | Gate == "Q8-LR" | Gate == "Q3-UL" | Gate == "Q3-UR" | Gate == "Q3-LL" | Gate == "Q3-LR")


# Make APOP combined gate for agranular and granular separately 
 Day7_2019_APOP_ALL_Granular_Apop_combined_no_outlier <- Day7_2019_APOP_ALL_Granular_Agranular_no_outlier  %>% filter(Gate == "Q8-LR" | Gate == "Q8-UR") %>% group_by(ID, Treat) %>% summarise(Percent_of_this_plot = sum(Percent_of_this_plot))
Day7_2019_APOP_ALL_Agranular_Apop_combined_no_outlier <- Day7_2019_APOP_ALL_Granular_Agranular_no_outlier  %>% filter(Gate == "Q3-LR" | Gate == "Q3-UR") %>% group_by(ID, Treat) %>% summarise(Percent_of_this_plot = sum(Percent_of_this_plot))

# Add new gate name for each 
 Day7_2019_APOP_ALL_Granular_Apop_combined_no_outlier$Gate <- "apop_combined_granular"
Day7_2019_APOP_ALL_Agranular_Apop_combined_no_outlier$Gate <- "apop_combined_agranular"

# Combined data frames for each cell type
Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined_no_outlier <- rbind(Day7_2019_APOP_ALL_Granular_Apop_combined_no_outlier, Day7_2019_APOP_ALL_Agranular_Apop_combined_no_outlier)

# Add arcsine transformed data
Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined_no_outlier$Percent_of_this_plot_arcsine <- transf.arcsin(Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined_no_outlier$Percent_of_this_plot*0.01)

# Merge with original data frame
Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined_no_outlier <- full_join(Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined_no_outlier,  Day7_2019_APOP_ALL_Granular_Agranular_no_outlier , by =c("Gate", "ID","Treat", "Percent_of_this_plot","Percent_of_this_plot_arcsine"))

# Boxplot of all gates of agranular and granular cells with faceting by gate, quadrants next to each other
# Change levels of faceting so that box plot quadrants are next to each other with agranular first, starting with live and dead
unique(Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined_no_outlier$Gate) # NOTE  Q3 = agranular, Q8= granular
#"apop_combined_granular"  "apop_combined_agranular" "Q3-UL"                   "Q3-UL"                  "Q3-LL"                  
#"Q3-LR"                   "Q8-UL"                   "Q8-UR"                   "Q8-LL"                   "Q8-LR"                  

Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined_no_outlier$Gate <- factor(Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined_no_outlier$Gate, 
                                                                   levels = c("Q3-UL", "Q8-UL","Q3-LL","Q8-LL","Q3-LR", "Q8-LR","Q3-UR", "Q8-UR","apop_combined_agranular","apop_combined_granular"))

# Same boxplot as above by facet by gate and not by treatment 
Day7_2019_APOP_Percent_Granular_Agranular_Hemocytes_no_outlier <- ggplot(data=Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined_no_outlier,
                                                              aes(y=Percent_of_this_plot, x=Gate, fill=Treat)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Quadrant") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Agranular and Granular FITC vs. PI Quad Plot") + 
  facet_grid(.~Gate, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(legend.text = element_text(size=12)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_manual(name="Treatment", labels=c("Dermo Injected","Notched Control"), values=c("#7e78d4","#cd4272")) +
  scale_x_discrete(labels=c("Q3-UL"= "Dead Agranular",
                            "Q8-UL"="Dead Granular", 
                            "Q3-LL"="Live Agranular", 
                            "Q8-LL"="Live Granular",
                            "Q3-LR"= "Live Agranular Apoptotic",   
                            "Q8-LR"="Live Granular Apoptotic",  
                            "Q3-UR"= "Dead Agranular Apoptotic",  
                            "Q8-UR"="Dead Granular Apoptotic", 
                            "apop_combined_agranular"="Combined Agranular Apoptotic",
                            "apop_combined_granular"="Combined Granular Apoptotic"))

# BoxPlot of just combined granular apoptotic and agranular combined apoptotic
Day7_2019_APOP_Percent_combined_apop_no_outlier <- Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined_no_outlier %>% filter(Gate== "apop_combined_agranular"| Gate == "apop_combined_granular") 
Day7_2019_APOP_ALL_Agranular_Granular_D_vs_NC_plot_no_outlier <-ggplot(data=Day7_2019_APOP_Percent_combined_apop_no_outlier,
                                                            aes(y=Percent_of_this_plot, x=Gate, fill=Treat)) + 
  geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Cell Type") +
  ylab("Percent Hemocytes") +
  ggtitle("Percent Agranular and Granular FITC vs. PI Quad Plot") +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=16), 
        axis.title.y=element_text(family="serif",size=16),
        axis.title.x=element_text(family="serif",size=16),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=16))  +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100))+
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text = element_text(size=16)) +
  scale_x_discrete(labels=c("apop_combined_granular"="Combined Granular Apoptotic", "apop_combined_agranular"="Combined Agranular Apoptotic")) +
  scale_fill_manual(name="Cell Type", labels=c("Dermo Injected","Notched Control"), values = c("#cc57b4", "#88bf3b")) + geom_signif(y_position = c(95,95), xmin = c(0.7,1.7), 
                                                                                                                                    xmax = c(1.3,2.3), annotation = c("**","***"),
                                                                                                                                    tip_length = 0.01)
# ANOVA
# Apop combined granular vs apop combined agranular within each treatment
# Notched control
Day7_2019_APOP_ALL_Granular_Agranular_Apop_NC_no_outlier <- Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined_no_outlier %>% filter(Gate =="apop_combined_agranular" | Gate== "apop_combined_granular") %>% filter(Treat == "NC")
Day7_2019_APOP_ALL_Granular_Agranular_Apop_NC_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_APOP_ALL_Granular_Agranular_Apop_NC_no_outlier)                                                  
summary(Day7_2019_APOP_ALL_Granular_Agranular_Apop_NC_aov_no_outlier)

        #Df Sum Sq Mean Sq F value   Pr(>F)    
        #Gate         1 0.9115  0.9115   57.28 1.13e-06 ***
        #  Residuals   16 0.2546  0.0159                     
        #---
        #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Dermo
Day7_2019_APOP_ALL_Granular_Agranular_Apop_D_no_outlier <- Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined_no_outlier %>% filter(Gate =="apop_combined_agranular" | Gate== "apop_combined_granular") %>% filter(Treat == "D")
Day7_2019_APOP_ALL_Granular_Agranular_Apop_D_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_APOP_ALL_Granular_Agranular_Apop_D_no_outlier)                                                  
summary(Day7_2019_APOP_ALL_Granular_Agranular_Apop_D_aov_no_outlier)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#Gate         1 0.1939 0.19393   29.76 5.29e-05 ***
#Residuals   16 0.1042 0.00652 

# Apop combined granular between treatments
Day7_2019_APOP_ALL_Granular_D_vs_NC_no_outlier <- Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined_no_outlier %>% filter(Gate== "apop_combined_granular") 
Day7_2019_APOP_ALL_Granular_D_vs_NC_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_APOP_ALL_Granular_D_vs_NC_no_outlier)                                                  
summary(Day7_2019_APOP_ALL_Granular_D_vs_NC_aov_no_outlier)

      #Df Sum Sq Mean Sq F value   Pr(>F)    
      #Treat        1 0.5750  0.5750   33.54 2.76e-05 ***
      #  Residuals   16 0.2743  0.0171                     
      #---
      #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Apop combined agranular between treatments 
Day7_2019_APOP_ALL_Agranular_D_vs_NC_no_outlier <- Day7_2019_APOP_ALL_Granular_Agranular_Apop_combined_no_outlier %>% filter(Gate== "apop_combined_agranular") 
Day7_2019_APOP_ALL_Agranular_D_vs_NC_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_APOP_ALL_Agranular_D_vs_NC_no_outlier)                                                  
summary(Day7_2019_APOP_ALL_Agranular_D_vs_NC_aov_no_outlier)

      #Df  Sum Sq Mean Sq F value  Pr(>F)   
      #Treat        1 0.05950 0.05950   11.26 0.00402 **
      #  Residuals   16 0.08455 0.00528                   
      #---
      #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## SUMMARY
# Significantly increased granulocytes within treatments, but between treatments notched control had significantly greater agranular and granular hemocyte apoptosis
# this is indicative of reduced apoptotic immune response following challenge or parasite inhibition of apoptosis

#### CASPASE ASSAY

Day7_2019_CASP_bad_removed_no_outlier
unique(Day7_2019_CASP_bad_removed_no_outlier[,c(2:5)])

# Q2= agranular, Q5= granular

# Analysis for all Granular and Agranular Combined caspase apoptotic from the quad plots
# Filter out granular and agranular quad plot
Day7_2019_CASP_ALL_Granular_Agranular_no_outlier <- Day7_2019_CASP_bad_removed_no_outlier %>% filter(Gate == "Q2-UL" | Gate == "Q2-UR" | Gate == "Q2-LL" | Gate == "Q2-LR" | Gate == "Q5-UL" | Gate == "Q5-UR" | Gate == "Q5-LL" | Gate == "Q5-LR")

# Make casp_apop combined gate for agranular and granular separately 
Day7_2019_CASP_ALL_Granular_Apop_combined_no_outlier <-  Day7_2019_CASP_ALL_Granular_Agranular_no_outlier  %>% filter(Gate == "Q5-LR" | Gate == "Q5-UR") %>% group_by(ID, Treat) %>% summarise(Percent_of_this_plot = sum(Percent_of_this_plot))
Day7_2019_CASP_ALL_Agranular_Apop_combined_no_outlier <- Day7_2019_CASP_ALL_Granular_Agranular_no_outlier  %>% filter(Gate == "Q2-LR" | Gate == "Q2-UR") %>% group_by(ID, Treat) %>% summarise(Percent_of_this_plot = sum(Percent_of_this_plot))

# Add new gate name for each 
 Day7_2019_CASP_ALL_Granular_Apop_combined_no_outlier$Gate <- "casp_apop_combined_granular"
Day7_2019_CASP_ALL_Agranular_Apop_combined_no_outlier$Gate <- "casp_apop_combined_agranular"

# Combined data frames for each cell type
Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined_no_outlier <- rbind(Day7_2019_CASP_ALL_Granular_Apop_combined_no_outlier, Day7_2019_CASP_ALL_Agranular_Apop_combined_no_outlier)

# Add arcsine transformed data
Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined_no_outlier$Percent_of_this_plot_arcsine <- transf.arcsin(Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined_no_outlier$Percent_of_this_plot*0.01)

# Merge with original data frame
Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined_no_outlier <- full_join(Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined_no_outlier,  Day7_2019_CASP_ALL_Granular_Agranular_no_outlier , by =c("Gate", "ID","Treat", "Percent_of_this_plot","Percent_of_this_plot_arcsine"))

# Boxplot of all gates of agranular and granular cells with faceting by gate, quadrants next to each other
# Change levels of faceting so that box plot quadrants are next to each other with agranular first, starting with live and dead
unique(Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined_no_outlier$Gate) # NOTE  Q2 = agranular, Q5= granular
#[1] "casp_apop_combined_granular"  "casp_apop_combined_agranular" "Q2-UL"                        "Q2-UR"                       
#[5] "Q2-LL"                        "Q2-LR"                        "Q5-UL"                        "Q5-UR"                       
#[9] "Q5-LL"                        "Q5-LR"                       

Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined_no_outlier$Gate <- factor(Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined_no_outlier$Gate, 
                                                                   levels = c("Q2-UL", "Q5-UL","Q2-LL","Q5-LL","Q2-LR", "Q5-LR","Q2-UR", "Q5-UR","casp_apop_combined_agranular","casp_apop_combined_granular"))

# Same boxplot as above by facet by gate and not by treatment 
Day7_2019_CASP_Percent_Granular_Agranular_Hemocytes_no_outlier <- ggplot(data=Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined_no_outlier,
                                                              aes(y=Percent_of_this_plot, x=Gate, fill=Treat)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Quadrant") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Agranular and Granular Caspase Apoptotic FITC vs. APC Quad Plot") + 
  facet_grid(.~Gate, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(legend.text = element_text(size=12)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_manual(name="Treatment", labels=c("Dermo Injected","Notched Control"), values=c("#7e78d4","#cd4272")) +
  scale_x_discrete(labels=c("Q2-UL"= "Dead Agranular",
                            "Q5-UL"="Dead Granular", 
                            "Q2-LL"="Live Agranular", 
                            "Q5-LL"="Live Granular",
                            "Q2-LR"= "Live Agranular Caspase 3/7 Apoptotic",   
                            "Q5-LR"="Live Granular Caspase 3/7 Apoptotic",  
                            "Q2-UR"= "Dead AgranularCaspase 3/7 Apoptotic",  
                            "Q5-UR"="Dead Granular Caspase 3/7 Apoptotic", 
                            "casp_apop_combined_agranular"="Combined Caspase 3/7 Agranular Apoptotic",
                            "casp_apop_combined_granular"="Combined Caspase 3/7 Granular Apoptotic"))

# ANOVA
# Casp active combined granular vs apop combined agranular within each treatment
# Notched control
Day7_2019_CASP_ALL_Granular_Agranular_Apop_NC_no_outlier <- Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined_no_outlier %>% filter(Gate =="casp_apop_combined_agranular" | Gate== "casp_apop_combined_granular") %>% filter(Treat == "NC")
Day7_2019_CASP_ALL_Granular_Agranular_Apop_NC_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_CASP_ALL_Granular_Agranular_Apop_NC_no_outlier)                                                  
summary(Day7_2019_CASP_ALL_Granular_Agranular_Apop_NC_aov_no_outlier)

      #Df Sum Sq Mean Sq F value   Pr(>F)    
      #Gate         1 2.7535  2.7535   332.4 3.97e-12 ***
      #  Residuals   16 0.1325  0.0083                     
      #---
      #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Dermo
Day7_2019_CASP_ALL_Granular_Agranular_Apop_D_no_outlier <- Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined_no_outlier %>% filter(Gate =="casp_apop_combined_agranular" | Gate== "casp_apop_combined_granular") %>% filter(Treat == "D")
Day7_2019_CASP_ALL_Granular_Agranular_Apop_D_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_CASP_ALL_Granular_Agranular_Apop_D_no_outlier)                                                  
summary(Day7_2019_CASP_ALL_Granular_Agranular_Apop_D_aov_no_outlier)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#Gate         1 2.5285  2.5285   80.56 1.21e-07 ***
#  Residuals   16 0.5022  0.0314                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Apop combined granular between treatments
Day7_2019_CASP_ALL_Granular_D_vs_NC_no_outlier <- Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined_no_outlier %>% filter(Gate== "casp_apop_combined_granular") 
Day7_2019_CASP_ALL_Granular_D_vs_NC_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_CASP_ALL_Granular_D_vs_NC_no_outlier)                                                  
summary(Day7_2019_CASP_ALL_Granular_D_vs_NC_aov_no_outlier)

      #Df  Sum Sq Mean Sq F value Pr(>F)
      #Treat        1 0.01413 0.01413   1.003  0.331
      #Residuals   16 0.22537 0.01409 

# Apop combined agranular between treatments 
Day7_2019_CASP_ALL_Agranular_D_vs_NC_no_outlier <- Day7_2019_CASP_ALL_Granular_Agranular_Apop_combined_no_outlier %>% filter(Gate== "casp_apop_combined_agranular") 
Day7_2019_CASP_ALL_Agranular_D_vs_NC_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_CASP_ALL_Agranular_D_vs_NC_no_outlier)                                                  
summary(Day7_2019_CASP_ALL_Agranular_D_vs_NC_aov_no_outlier)

      #Df Sum Sq Mean Sq F value Pr(>F)
      #Treat        1 0.0354 0.03540   1.384  0.257
      #Residuals   16 0.4094 0.02558  

## SUMMARY
# Significantly increased granulocytes caspase 3/7 activation, but no difference between treatments

#### LMP Assay Statistics ####

Day7_2019_LMP_bad_removed_no_outlier
unique(Day7_2019_LMP_bad_removed_no_outlier[,c(2:5)])

# Q3 = agranular, Q2 = granular

# Analysis for all Granular and Agranular Combined caspase apoptotic from the quad plots
# Filter out granular and agranular quad plot
Day7_2019_LMP_ALL_Granular_Agranular_no_outlier <- Day7_2019_LMP_bad_removed_no_outlier %>% filter(Gate == "Q3-UL" | Gate == "Q3-UR" | Gate == "Q3-LL" | Gate == "Q3-LR" | Gate == "Q2-UL" | Gate == "Q2-UR" | Gate == "Q2-LL" | Gate == "Q2-LR")

# Boxplot of all gates of agranular and granular cells with faceting by gate, quadrants next to each other
# Change levels of faceting so that box plot quadrants are next to each other with agranular first, starting with live and dead

unique(Day7_2019_LMP_ALL_Granular_Agranular_no_outlier$Gate) # NOTE Q3 = agranular, Q2 = granular
# "Q3-UL" "Q3-UR" "Q3-LL" "Q3-LR" "Q2-UL" "Q2-UR" "Q2-LL" "Q2-LR"                     

Day7_2019_LMP_ALL_Granular_Agranular_no_outlier$Gate <- factor(Day7_2019_LMP_ALL_Granular_Agranular_no_outlier$Gate, 
                                                    levels = c("Q3-UL","Q2-UL", "Q3-LL","Q2-LL","Q3-UR", "Q2-UR", "Q3-LR","Q2-LR"))

# Boxplot as above by facet by gate and not by treatment 
Day7_2019_LMP_Percent_Granular_Agranular_Hemocytes_no_outlier <- ggplot(data=Day7_2019_LMP_ALL_Granular_Agranular_no_outlier,
                                                             aes(y=Percent_of_this_plot, x=Gate, fill=Treat)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Quadrant") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Agranular and Granular Intact Lysosome FITC vs. PI Quad Plot") + 
  facet_grid(.~Gate, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, hjust=1, angle=70)) +
  theme(legend.text = element_text(size=12)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_manual(name="Treatment", labels=c("Dermo Injected","Notched Control"), values=c("#7e78d4","#cd4272")) +
  scale_x_discrete(labels=c("Q3-UL"="agranular_lysosome_rupture",
                            "Q2-UL"="granular_lysosome_rupture", 
                            "Q3-LL"="agranular_LMP_onset",
                            "Q2-LL"="granular_LMP_onset", 
                            "Q3-UR"="agranular_dead_other_means", 
                            "Q2-UR"="granular_dead_other_means",
                            "Q3-LR"="agranular_lysosome_live",
                            "Q2-LR"="granular_lysosome_live"))



# ANOVA

## LMP RUPTURE ANOVAS
# LMP rupture granular vs LMP rupture granular within each treatment
# Notched control
Day7_2019_LMP_ALL_Granular_Agranular_rupture_NC_no_outlier <- Day7_2019_LMP_ALL_Granular_Agranular_no_outlier %>% filter(Gate =="Q3-UL" | Gate== "Q2-UL") %>% filter(Treat == "NC")
Day7_2019_LMP_ALL_Granular_Agranular_rupture_NC_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_LMP_ALL_Granular_Agranular_rupture_NC_no_outlier)                                                  
summary(Day7_2019_LMP_ALL_Granular_Agranular_rupture_NC_aov_no_outlier)

        #Df Sum Sq Mean Sq F value Pr(>F)  
        #Gate         1 0.0837 0.08370   5.931  0.027 *
        #  Residuals   16 0.2258 0.01411                 
        #---
        #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Dermo
Day7_2019_LMP_ALL_Granular_Agranular_rupture_D_no_outlier <- Day7_2019_LMP_ALL_Granular_Agranular_no_outlier %>% filter(Gate =="Q3-UL" | Gate== "Q2-UL") %>% filter(Treat == "D")
Day7_2019_LMP_ALL_Granular_Agranular_rupture_D_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_LMP_ALL_Granular_Agranular_rupture_D_no_outlier)                                                  
summary(Day7_2019_LMP_ALL_Granular_Agranular_rupture_D_aov_no_outlier)

#Df  Sum Sq Mean Sq F value Pr(>F)  
#Gate         1 0.06203 0.06203       8 0.0121 *
#  Residuals   16 0.12407 0.00775                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# LMP rupture granular between treatments
# Q3 = agranular, Q2 = granular
Day7_2019_LMP_ALL_Granular_rupture_no_outlier <- Day7_2019_LMP_ALL_Granular_Agranular_no_outlier %>% filter(Gate== "Q2-UL") 
Day7_2019_LMP_ALL_Granular_rupture_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_LMP_ALL_Granular_rupture_no_outlier )                                                  
summary(Day7_2019_LMP_ALL_Granular_rupture_aov_no_outlier)

      #Df  Sum Sq Mean Sq F value Pr(>F)
      #Treat        1 0.01744 0.01744   1.335  0.265
      #Residuals   16 0.20901 0.01306  

# LMP rupture agranular between treatments
Day7_2019_LMP_ALL_Agranular_rupture_no_outlier <- Day7_2019_LMP_ALL_Granular_Agranular %>% filter(Gate== "Q3-UL") 
Day7_2019_LMP_ALL_Agranular_rupture_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_LMP_ALL_Agranular_rupture_no_outlier )                                                  
summary(Day7_2019_LMP_ALL_Agranular_rupture_aov_no_outlier)

#Df Sum Sq  Mean Sq F value Pr(>F)
#Treat        1 0.0115 0.011500   1.441  0.246
#Residuals   18 0.1437 0.007981   

## LMP ONSET ANOVAs
# LMP onset granular vs LMP rupture granular within each treatment
# Notched control
Day7_2019_LMP_ALL_Granular_Agranular_onset_NC_no_outlier <- Day7_2019_LMP_ALL_Granular_Agranular %>% filter(Gate =="Q3-LL" | Gate== "Q2-LL") %>% filter(Treat == "NC")
Day7_2019_LMP_ALL_Granular_Agranular_onset_NC_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_LMP_ALL_Granular_Agranular_onset_NC_no_outlier )                                                  
summary(Day7_2019_LMP_ALL_Granular_Agranular_onset_NC_aov_no_outlier)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#Gate         1 2.8191  2.8191   333.2 6.12e-14 ***
#  Residuals   20 0.1692  0.0085                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Dermo
Day7_2019_LMP_ALL_Granular_Agranular_onset_D_no_outlier <- Day7_2019_LMP_ALL_Granular_Agranular_no_outlier %>% filter(Gate =="Q3-LL" | Gate== "Q2-LL") %>% filter(Treat == "D")
Day7_2019_LMP_ALL_Granular_Agranular_onset_D_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_LMP_ALL_Granular_Agranular_onset_D_no_outlier)                                                  
summary(Day7_2019_LMP_ALL_Granular_Agranular_onset_D_aov_no_outlier)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#Gate         1 1.7932  1.7932   200.7 1.79e-10 ***
#  Residuals   16 0.1429  0.0089                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# LMP onset granular between treatments
# Q3 = agranular, Q2 = granular
Day7_2019_LMP_ALL_Granular_onset_no_outlier <- Day7_2019_LMP_ALL_Granular_Agranular_no_outlier %>% filter(Gate== "Q2-LL") 
Day7_2019_LMP_ALL_Granular_onset_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_LMP_ALL_Granular_onset_no_outlier )                                                  
summary(Day7_2019_LMP_ALL_Granular_onset_aov_no_outlier)

      #Df  Sum Sq  Mean Sq F value Pr(>F)
      #Treat        1 0.00076 0.000755   0.204  0.657
      #Residuals   16 0.05909 0.003693  

# LMP onset agranular between treatments
Day7_2019_LMP_ALL_Agranular_onset_no_outlier <- Day7_2019_LMP_ALL_Granular_Agranular_no_outlier %>% filter(Gate== "Q3-LL") 
Day7_2019_LMP_ALL_Agranular_onset_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_LMP_ALL_Agranular_onset_no_outlier )                                                  
summary(Day7_2019_LMP_ALL_Agranular_onset_aov_no_outlier)

      #Df Sum Sq Mean Sq F value Pr(>F)  
      #Treat        1 0.0572 0.05720   5.624 0.0306 *
      #  Residuals   16 0.1627 0.01017                 
      #---
      #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

### LMP Intact lysosome live 
# LMP intact granular vs LMP rupture granular within each treatment
# Notched control
Day7_2019_LMP_ALL_Granular_Agranular_intact_NC_no_outlier <- Day7_2019_LMP_ALL_Granular_Agranular_no_outlier %>% filter(Gate =="Q3-LR" | Gate== "Q2-LR") %>% filter(Treat == "NC")
Day7_2019_LMP_ALL_Granular_Agranular_intact_NC_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_LMP_ALL_Granular_Agranular_intact_NC_no_outlier )                                                  
summary(Day7_2019_LMP_ALL_Granular_Agranular_intact_NC_aov_no_outlier)

        #Df Sum Sq Mean Sq F value Pr(>F)  
        #Gate         1 0.1972 0.19716   4.517 0.0495 *
        #  Residuals   16 0.6984 0.04365                 
        #---
        #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Dermo
Day7_2019_LMP_ALL_Granular_Agranular_intact_D_no_outlier <- Day7_2019_LMP_ALL_Granular_Agranular_no_outlier %>% filter(Gate =="Q3-LR" | Gate== "Q2-LR") %>% filter(Treat == "D")
Day7_2019_LMP_ALL_Granular_Agranular_intact_D_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Gate, data=Day7_2019_LMP_ALL_Granular_Agranular_intact_D_no_outlier)                                                  
summary(Day7_2019_LMP_ALL_Granular_Agranular_intact_D_aov_no_outlier)

#Df Sum Sq Mean Sq F value Pr(>F)  
#Gate         1 0.1669 0.16689   6.614 0.0205 *
# Residuals   16 0.4038 0.02523                 
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# LMP intact granular between treatments
# Q3 = agranular, Q2 = granular
Day7_2019_LMP_ALL_Granular_intact_no_outlier <- Day7_2019_LMP_ALL_Granular_Agranular_no_outlier %>% filter(Gate== "Q2-LR") 
Day7_2019_LMP_ALL_Granular_intact_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_LMP_ALL_Granular_intact_no_outlier )                                                  
summary(Day7_2019_LMP_ALL_Granular_intact_aov_no_outlier)

    #Df Sum Sq Mean Sq F value Pr(>F)
    #Treat        1 0.0841 0.08408   1.565  0.229
    #Residuals   16 0.8595 0.05372    

# LMP intact agranular between treatments
Day7_2019_LMP_ALL_Agranular_intact_no_outlier <- Day7_2019_LMP_ALL_Granular_Agranular_no_outlier %>% filter(Gate== "Q3-LR") 
Day7_2019_LMP_ALL_Agranular_intact_aov_no_outlier <- aov(Percent_of_this_plot_arcsine ~ Treat, data=Day7_2019_LMP_ALL_Agranular_intact_no_outlier )                                                  
summary(Day7_2019_LMP_ALL_Agranular_intact_aov_no_outlier)

    #Df Sum Sq Mean Sq F value Pr(>F)  
    #Treat        1 0.1059 0.10593   6.986 0.0177 *
    #  Residuals   16 0.2426 0.01516                 
    #---
    #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


######## RE-ANALYSIS OF LAST YEAR'S DATA ########

##### LOAD ASSAY CSV'S and format #####

# Make new header column
VIA_nms_Day7_2018 <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2018_Dermo_Challenge/Data/ANALYSIS FILES/FORMATTED_CSVS/DAY7/VIA_DAY7_UPDATED_2018_DATA_header_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

VIA_nms_Day50_2018 <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2018_Dermo_Challenge/Data/ANALYSIS FILES/FORMATTED_CSVS/DAY50/VIA_DAY50_UPDATED_2018_DATA_header_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Day7_2018_VIA <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2018_Dermo_Challenge/Data/ANALYSIS FILES/FORMATTED_CSVS/DAY7/VIA_DAY7_UPDATED_2018_DATA_header_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(VIA_nms_Day7_2018)

Day50_2018_VIA <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2018_Dermo_Challenge/Data/ANALYSIS FILES/FORMATTED_CSVS/DAY50/VIA_DAY50_UPDATED_2018_DATA_header_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(VIA_nms_Day50_2018)

# check that they are the same
colnames(Day7_2018_VIA) %in% colnames(Day50_2018_VIA) # all return TRUE
colnames(Day50_2018_VIA) %in% colnames(Day7_2018_VIA) # all return TRUE

# compare columns with the 2018 data as well because they should be the same
colnames(Day7_2018_VIA) %in% colnames(Day7_2019_VIA) # all return TRUE plot 6 is missing percent of this plot in 2019 plots, the
  # the counts of dead cells to subtract from total is better than subtracting percentage so its fine that this is missing
colnames(Day7_2019_VIA) %in% colnames(Day7_2018_VIA) # all return TRUE

# Split column 1 by space and remove
head(Day7_2018_VIA)
Day7_2018_VIA <- Day7_2018_VIA %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID"))

head(Day50_2018_VIA)
Day50_2018_VIA <- Day50_2018_VIA %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID"))

# Separate new column 1 by dash, remove spaces from column names
Day7_2018_VIA  <- Day7_2018_VIA [,-1] %>% separate(ID, sep="-", into=c("Family","ID","Assay"))
colnames(Day7_2018_VIA)
head(Day7_2018_VIA)
Day7_2018_VIA_percent <- Day7_2018_VIA [,c(1:3,5,7,9,11,15,17)] # keep header and percent columns 
Day7_2018_VIA_counts <- Day7_2018_VIA [,c(1:4,6,8,10,12:14,16)] # keep header and counts columns 

Day50_2018_VIA  <- Day50_2018_VIA [,-1] %>% separate(ID, sep="-", into=c("Family","ID","Assay"))
colnames(Day50_2018_VIA)
Day50_2018_VIA_percent <- Day50_2018_VIA [,c(1:3,5,7,9,11,15,17)] # keep header and percent columns 
Day50_2018_VIA_counts <-  Day50_2018_VIA [,c(1:4,6,8,10,12:14,16)] # keep header and counts columns 

# Gather count and percent columns separately
colnames(Day7_2018_VIA_percent)
colnames(Day7_2018_VIA_counts)
Day7_2018_VIA_percent <- Day7_2018_VIA_percent %>% group_by(Family,ID,Assay) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:9))
Day7_2018_VIA_counts <-  Day7_2018_VIA_counts %>% group_by(Family,ID,Assay) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:11))

Day7_2018_VIA_percent <- Day7_2018_VIA_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Day7_2018_VIA_counts <-  Day7_2018_VIA_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

# remove columns for plot_name, channel, %
colnames(Day7_2018_VIA_percent)
colnames(Day7_2018_VIA_counts)
Day7_2018_VIA_percent <- Day7_2018_VIA_percent[,-c(4,6,8)]
Day7_2018_VIA_counts <-  Day7_2018_VIA_counts[,-c(4,6,8)]

colnames(Day50_2018_VIA_percent)
colnames(Day50_2018_VIA_counts)
Day50_2018_VIA_percent <- Day50_2018_VIA_percent %>% group_by(Family,ID,Assay) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:9))
Day50_2018_VIA_counts <-  Day50_2018_VIA_counts %>% group_by(Family,ID,Assay) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:11))

Day50_2018_VIA_percent <- Day50_2018_VIA_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Day50_2018_VIA_counts <-  Day50_2018_VIA_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

# remove columns for plot_name, channel, %
colnames(Day50_2018_VIA_percent)
colnames(Day50_2018_VIA_counts)
Day50_2018_VIA_percent <- Day50_2018_VIA_percent[,-c(4,6,8)]
Day50_2018_VIA_counts <-  Day50_2018_VIA_counts[,-c(4,6,8)]

# Full join together so all columns preserved
Day7_2018_VIA_join <- full_join(Day7_2018_VIA_percent, Day7_2018_VIA_counts, by = c("Family","ID","Assay","Plot_number","Gate"))
Day50_2018_VIA_join <- full_join(Day50_2018_VIA_percent, Day50_2018_VIA_counts, by = c("Family","ID","Assay","Plot_number","Gate"))

# Add Day column for each
Day7_2018_VIA_join$Day <- "7"
Day50_2018_VIA_join$Day <- "50"

# Combine days into one dataframe
Day7_Day50_2018_VIA_join <- full_join(Day7_2018_VIA_join , Day50_2018_VIA_join , by = c("Family","ID","Assay","Plot_number","Gate","Counts","Percent_of_this_plot","Day"))
unique_Day7_Day50_2018_VIA_join <- unique(Day7_Day50_2018_VIA_join[,c("Plot_number", "Gate")])
#  Plot_number Gate 
# <chr>       <chr>
#  2           M4   
#  3           P1   
#  4           E1   
#  4           E3   
#  6           V2-L 
#  6           V2-R 
#  9           This 
#  10          This

unique_Day7_2019_VIA_join_gate$Gate %in% unique_Day7_Day50_2018_VIA_join$Gate # all return TRUE
unique_Day7_Day50_2018_VIA_join$Gate %in% unique_Day7_2019_VIA_join_gate$Gate # all return TRUE

# Add in cell type column based on plot number and gate
VIA_cell_type <- data.frame(Plot_number = c("2","3","4","4","9","10","6","6"), Gate= c("M4","P1","E1","E3","This","This","V2-L","V2-R"), Cell_type=c(
  "all_hemocytes","all_hemocytes", "granular","agranular",
  "live_granular", "live_agranular", "all_live_hemocytes","all_dead_hemocytes"))
# Join cell type
Day7_Day50_2018_VIA_join  <- left_join(Day7_Day50_2018_VIA_join , VIA_cell_type, by=c("Plot_number","Gate"))

## APOPTOSIS ASSAY
# Make new header column
Apop_nms_Day7_2018 <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2018_Dermo_Challenge/Data/ANALYSIS FILES/FORMATTED_CSVS/DAY7/APOP_DAY7_UPDATED_2018_DATA_formatted_header.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

Apop_nms_Day50_2018 <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2018_Dermo_Challenge/Data/ANALYSIS FILES/FORMATTED_CSVS/DAY50/APOP_DAY50_UPDATED_2018_DATA_header_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Day7_2018_APOP <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2018_Dermo_Challenge/Data/ANALYSIS FILES/FORMATTED_CSVS/DAY7/APOP_DAY7_UPDATED_2018_DATA_formatted_header.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(Apop_nms_Day7_2018)
Day50_2018_APOP <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2018_Dermo_Challenge/Data/ANALYSIS FILES/FORMATTED_CSVS/DAY50/APOP_DAY50_UPDATED_2018_DATA_header_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(Apop_nms_Day50_2018 )

# check that they are the same
colnames(Day7_2018_APOP) %in% colnames(Day50_2018_APOP) # all return TRUE
colnames(Day50_2018_APOP) %in% colnames(Day7_2018_APOP) # all return TRUE

colnames(Day7_2018_APOP) %in% colnames(Day7_2019_APOP) # the 2019 has just already been processed to have the ID treatment and assay columns combined
colnames(Day7_2019_APOP) %in% colnames(Day7_2018_APOP)

# Split column 1 by space and remove
Day7_2018_APOP <- Day7_2018_APOP %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID"))

Day50_2018_APOP <- Day50_2018_APOP %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID"))

# Separate new column 1 by dash, remove spaces from column names
Day7_2018_APOP <- Day7_2018_APOP[,-1] %>% separate(ID, sep="-", into=c("Family","ID","Assay"))
colnames(Day7_2018_APOP)
head(Day7_2018_APOP)
Day7_2018_APOP_percent <- Day7_2018_APOP[,c(1:3,7,9,11,13,15,17,19,21, 23,25,27,29,31,33,35)]
Day7_2018_APOP_counts <-  Day7_2018_APOP[,c(1:6,8,10,12,14,16,18,20,22,24,26,28,30,32,34)]

Day50_2018_APOP <- Day50_2018_APOP[,-1] %>% separate(ID, sep="-", into=c("Family","ID","Assay"))
colnames(Day50_2018_APOP)
head(Day50_2018_APOP)
Day50_2018_APOP_percent <- Day50_2018_APOP[,c(1:3,7,9,11,13,15,17,19,21, 23,25,27,29,31,33,35)]
Day50_2018_APOP_counts <-  Day50_2018_APOP[,c(1:6,8,10,12,14,16,18,20,22,24,26,28,30,32,34)]

# Gather count and percent columns separately
colnames(Day7_2018_APOP_percent)
colnames(Day7_2018_APOP_counts)
Day7_2018_APOP_percent <- Day7_2018_APOP_percent %>% group_by(Family,ID,Assay) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:18))
Day7_2018_APOP_counts <-  Day7_2018_APOP_counts %>% group_by(Family,ID,Assay) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:20))

Day7_2018_APOP_percent <- Day7_2018_APOP_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Day7_2018_APOP_counts <-  Day7_2018_APOP_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

colnames(Day7_2018_APOP_percent)
colnames(Day7_2018_APOP_counts)
Day7_2018_APOP_percent <- Day7_2018_APOP_percent[,-c(4,6,8)]
Day7_2018_APOP_counts <-  Day7_2018_APOP_counts[,-c(4,6,8)]

colnames(Day50_2018_APOP_percent)
colnames(Day50_2018_APOP_counts)
Day50_2018_APOP_percent <- Day50_2018_APOP_percent %>% group_by(Family,ID,Assay) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:18))
Day50_2018_APOP_counts <-  Day50_2018_APOP_counts %>% group_by(Family,ID,Assay) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:20))

Day50_2018_APOP_percent <- Day50_2018_APOP_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Day50_2018_APOP_counts <-  Day50_2018_APOP_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

Day50_2018_APOP_percent <- Day50_2018_APOP_percent[,-c(4,6,8)]
Day50_2018_APOP_counts <-  Day50_2018_APOP_counts[,-c(4,6,8)]

# Full join together so all columns preserved
Day7_2018_APOP_join <- full_join(Day7_2018_APOP_percent, Day7_2018_APOP_counts, by = c("Family","ID","Assay","Plot_number","Gate"))
Day50_2018_APOP_join <- full_join(Day50_2018_APOP_percent, Day50_2018_APOP_counts, by = c("Family","ID","Assay","Plot_number","Gate"))
unique_Day50_2018_APOP_join <- unique(Day50_2018_APOP_join[,c("Plot_number", "Gate")])
unique_Day7_2018_APOP_join <- unique(Day7_2018_APOP_join[,c("Plot_number", "Gate")])

# check that the same columns are present between days and years
unique_Day50_2018_APOP_join$Gate %in% unique_Day7_2018_APOP_join$Gate # all return TRUE
unique_Day7_2018_APOP_join$Gate %in% unique_Day50_2018_APOP_join$Gate # all return TRUE
unique_Day7_2019_APOP_join$Gate %in% unique_Day7_2018_APOP_join$Gate # 8 columns don't match

#   Plot_number Gate  
#<chr>       <chr> 
# 3           P1    
# 8           P3    
# 8           P4    
# 4           Q7-UL 
# 4           Q7-UR 
# 4           Q7-LL 
# 4           Q7-LR 
# 7           Q12-UL
# 7           Q12-UR
# 7           Q12-LL
# 7           Q12-LR
# 6           Q13-UR
# 6           Q13-UL
# 6           Q13-LL
# 6           Q13-LR
# 12          This  
# 13          This  

# Add in cell type column based on plot number and gate, this is different than the 2019 gates

APOP_cell_type_2018 <- data.frame(Plot_number = c("3",  "8",  "8",  "4",  "4",  "4",  "4",  
                "7",  "7",  "7",  "7", "6",  "6",  "6",  "6",  "12", "13" ), Gate= c("P1",  "P3",    "P4",  
                                                              "Q7-UL", "Q7-UR", "Q7-LL", "Q7-LR", 
                                                              "Q12-UL","Q12-UR","Q12-LL","Q12-LR","Q13-UR",
                                                              "Q13-UL","Q13-LL","Q13-LR","This",  "This"), 
                Cell_type = c("all_hemocytes", "agranular", "granular", 
                              "all_hemocytes_dead", "all_hemocytes_dead_apoptotic", "all_hemocytes_live", "all_hemocytes_live_apoptotic",
                              "agranular_dead", "agranular_dead_apoptotic",
                              "agranular_live", "agranular_live_apoptotic", "granular_dead", "granular_dead_apoptotic", 
                              "granular_live", "granular_live_apoptotic","live_agranular","live_granular"))

# Join cell type
Day7_2018_APOP_join <- left_join(Day7_2018_APOP_join, APOP_cell_type_2018, by=c("Plot_number","Gate"))
Day50_2018_APOP_join <- left_join(Day50_2018_APOP_join, APOP_cell_type_2018, by=c("Plot_number","Gate"))

# Add Day column for each
Day7_2018_APOP_join$Day <- "7"
Day50_2018_APOP_join$Day <- "50"

# Combine days into one dataframe
Day7_Day50_2018_APOP_join <- full_join(Day7_2018_APOP_join , Day50_2018_APOP_join , by = c("Family","ID","Assay","Plot_number","Gate","Counts","Percent_of_this_plot","Day","Cell_type"))

## CASPASE ASSAY
# Make new header column
Casp_nms_Day7_2018 <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2018_Dermo_Challenge/Data/ANALYSIS FILES/FORMATTED_CSVS/DAY7/CASP_DAY7_UPDATED_2018_DATA_header_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

Casp_nms_Day50_2018 <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2018_Dermo_Challenge/Data/ANALYSIS FILES/FORMATTED_CSVS/DAY50/CASP_DAY50_UPDATED_2018_DATA_header_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Day7_2018_CASP <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2018_Dermo_Challenge/Data/ANALYSIS FILES/FORMATTED_CSVS/DAY7/CASP_DAY7_UPDATED_2018_DATA_header_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(Casp_nms_Day7_2018)

Day50_2018_CASP <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2018_Dermo_Challenge/Data/ANALYSIS FILES/FORMATTED_CSVS/DAY50/CASP_DAY50_UPDATED_2018_DATA_header_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(Casp_nms_Day50_2018)

# check that they are the same
colnames(Day7_2018_CASP) %in% colnames(Day50_2018_CASP) # all return TRUE
colnames(Day50_2018_CASP) %in% colnames(Day7_2018_CASP) # all return TRUE

# Split column 1 by space and remove
Day7_2018_CASP <- Day7_2018_CASP %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 
Day50_2018_CASP <- Day50_2018_CASP %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
Day7_2018_CASP <- Day7_2018_CASP[,-1] %>% separate(ID, sep="-", into=c("Family","ID","Assay"))
colnames(Day7_2018_CASP)
head(Day7_2018_CASP)
Day7_2018_CASP_percent <- Day7_2018_CASP[,c(1:3,4,7,9,11,13,15,17,21,23,25,27,29,31,33,35)]
Day7_2018_CASP_counts <-  Day7_2018_CASP[,c(1:3,5,6,8,10,12,14,16,18,19,20,22,24,26,28,30,32,34)]

Day50_2018_CASP <- Day50_2018_CASP[,-1] %>% separate(ID, sep="-", into=c("Family","ID","Assay"))
colnames(Day50_2018_CASP)
Day50_2018_CASP_percent <- Day50_2018_CASP[,c(1:3,4,7,9,11,13,15,17,21,23,25,27,29,31,33,35)]
Day50_2018_CASP_counts <-  Day50_2018_CASP[,c(1:3,5,6,8,10,12,14,16,18,19,20,22,24,26,28,30,32,34)]

# Gather count and percent columns separately
colnames(Day7_2018_CASP_percent)
colnames(Day7_2018_CASP_counts)
Day7_2018_CASP_percent_gather <- Day7_2018_CASP_percent %>% group_by(Family,ID,Assay) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:18))
Day7_2018_CASP_counts_gather <-  Day7_2018_CASP_counts %>% group_by(Family,ID,Assay) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:20))

Day7_2018_CASP_percent <- Day7_2018_CASP_percent_gather %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Day7_2018_CASP_counts <-  Day7_2018_CASP_counts_gather %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

head(Day7_2018_CASP_percent)
colnames(Day7_2018_CASP_counts)
Day7_2018_CASP_percent <- Day7_2018_CASP_percent[,-c(4,6,8)]
Day7_2018_CASP_counts <-  Day7_2018_CASP_counts[,-c(4,6,8)]

Day50_2018_CASP_percent_gather <- Day50_2018_CASP_percent %>% group_by(Family,ID,Assay) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:18))
Day50_2018_CASP_counts_gather <-  Day50_2018_CASP_counts %>% group_by(Family,ID,Assay) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:20))

Day50_2018_CASP_percent <- Day50_2018_CASP_percent_gather %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Day50_2018_CASP_counts <-  Day50_2018_CASP_counts_gather %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

colnames(Day50_2018_CASP_percent)
colnames(Day50_2018_CASP_counts)
Day50_2018_CASP_percent <- Day50_2018_CASP_percent[,-c(4,6,8)]
Day50_2018_CASP_counts <-  Day50_2018_CASP_counts[,-c(4,6,8)]

# Full join together so all columns preserved
Day7_2018_CASP_join <- full_join(Day7_2018_CASP_percent, Day7_2018_CASP_counts, by = c("Family","ID","Assay","Plot_number","Gate"))
Day50_2018_CASP_join <- full_join(Day50_2018_CASP_percent, Day50_2018_CASP_counts, by = c("Family","ID","Assay","Plot_number","Gate"))

unique_plot_gate_Casp_2018_day7 <- unique(Day7_2018_CASP_join[,c("Plot_number", "Gate")])
unique_plot_gate_Casp_2018_day50 <- unique(Day50_2018_CASP_join[,c("Plot_number", "Gate")])

unique_plot_gate_Casp_2018_day7$Gate %in% unique_plot_gate_Casp_2018_day50$Gate # All return true
unique_plot_gate_Casp_2018_day50$Gate %in% unique_plot_gate_Casp_2018_day7$Gate
unique_Day7_2019_CASP_join$Gate %in% unique_plot_gate_Casp_2018_day50$Gate # 8 are different between last year and this year

unique_plot_gate_Casp_2018_day7
#Plot_number Gate 
#<chr>       <chr>
# 3           P1   
# 8           P3   
# 8           P4   
# 4           Q2-UL
# 4           Q2-UR
# 4           Q2-LL
# 4           Q2-LR
# 7           Q6-UL
# 7           Q6-UR
# 7           Q6-LL
# 7           Q6-LR
# 2           Q7-UL
# 2           Q7-UR
# 2           Q7-LL
# 2           Q7-LR
# 6           This 
# 9           This 

# Add in cell type column based on plot number and gate
CASP_cell_type_2018 <- data.frame(Plot_number = c("3","8","8","4","4","4","4","7","7","7","7","2","2","2","2","6","9"), 
                             Gate= c("P1","P3","P4","Q2-UL","Q2-UR","Q2-LL","Q2-LR","Q6-UL","Q6-UR","Q6-LL","Q6-LR","Q7-UL","Q7-UR","Q7-LL","Q7-LR","This","This"), 
                             Cell_type=c("all_hemocytes", "agranular", "granular", "agranular_dead", "agranular_dead_caspase_active", 
                                         "agranular_live", "agranular_live_caspase_active", "granular_dead", "granular_dead_caspase_active", 
                                         "granular_live", "granular_live_caspase_active", "all_dead", "all_dead_caspase_active", "all_live", 
                                         "all_live_caspase_active", "live_agranular", "live_granular"))
# Join cell type
Day7_2018_CASP_join <- left_join(Day7_2018_CASP_join, CASP_cell_type_2018, by=c("Plot_number","Gate"))
Day50_2018_CASP_join <- left_join(Day50_2018_CASP_join, CASP_cell_type_2018, by=c("Plot_number","Gate"))

# Add Day column for each
Day7_2018_CASP_join$Day <- "7"
Day50_2018_CASP_join$Day <- "50"

# Combine days into one dataframe
Day7_Day50_2018_CASP_join <- full_join(Day7_2018_CASP_join , Day50_2018_CASP_join , by = c("Family","ID","Assay","Plot_number","Gate","Counts","Percent_of_this_plot","Day","Cell_type"))

#### Combine all 2018 Data frames ####

Day7_Day50_2018_all_assays <- full_join(Day7_Day50_2018_VIA_join, Day7_Day50_2018_APOP_join, by =c("Family","ID","Assay",
                                                                                               "Plot_number", "Gate", "Cell_type", "Day", "Counts", "Percent_of_this_plot"))
Day7_Day50_2018_all_assays <- full_join(Day7_Day50_2018_all_assays, Day7_Day50_2018_CASP_join, by =c("Family","ID","Assay",
                                                                                                     "Plot_number", "Gate", "Cell_type", "Day", "Counts", "Percent_of_this_plot"))
# Create combined ID column 
Day7_Day50_2018_all_assays$Flow_Code <- paste(Day7_Day50_2018_all_assays$Family, Day7_Day50_2018_all_assays$ID, sep="-")
head(Day7_Day50_2018_all_assays)

#### Add oyster ID's for Day7 #####

Day7_Day_50_ID_key <- read_csv("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2018_Dermo_Challenge/Data/ANALYSIS FILES/FORMATTED_CSVS/Day_7_Day_50_ID_Key_R.csv")
class(Day7_Day_50_ID_key$Day)
class(Day7_Day50_2018_all_assays$Day)
Day7_Day_50_ID_key$Day <- as.character(Day7_Day_50_ID_key$Day)

Day7_Day50_2018_all_assays <- left_join(Day7_Day50_2018_all_assays, Day7_Day_50_ID_key, by =c("Family","Flow_Code","Day"))
head(Day7_Day50_2018_all_assays)

#### Remove bad samples from all 2018 analysis ####
#SAMPLES TO REMOVE FROM DAY 7 analysis = A-09, J-09 ,J-12 ,D-04 ,E-05 ,E-09 ,L-07,L-09,L-13
#SAMPLES TO REMOVE FROM DAY 50 analysis = J-25, J-40, J-166, E-30  (was E-36), D-23  , D-142  ,D-187  ,E-18  (was E-12), B-154   

# vector with ID lines for each assay
ID_samples_remove_2018 <- c("A-09", "J-09" ,"J-12" ,"D-04" ,"E-05" ,"E-09" ,"L-07","L-09","L-13",
                            "J-25", "J-40", "J-166", "E-30", "D-23", "D-142"  ,"D-187"  ,"E-18", "B-154")

#Remove all of those lines with those ID's from the data
Day7_Day50_2018_all_assays_bad_removed <-  Day7_Day50_2018_all_assays[!(Day7_Day50_2018_all_assays$Flow_Code %in% ID_samples_remove_2018),] 

# check output
unique(Day7_Day50_2018_all_assays_bad_removed$Flow_Code) # output is missing correct samples and has correct format 

#### Arcsine Transformation of percentages ####
Day7_Day50_2018_all_assays_bad_removed$Percent_of_this_plot_arcsine <- transf.arcsin(Day7_Day50_2018_all_assays_bad_removed$Percent_of_this_plot)

# Change percent column to percent
Day7_Day50_2018_all_assays_bad_removed$Percent_of_this_plot <- Day7_Day50_2018_all_assays_bad_removed$Percent_of_this_plot*100

write.csv(file="Day7_Day50_2018_all_assays_bad_removed.csv", Day7_Day50_2018_all_assays_bad_removed)

#### PICK SAMPLES FOR SEQUENCING ####

### Create spread sheet with combined load data, apoptosis data, caspase, viability assay data side by side for comparison to pick samples for sequence
## GOAL: Select samples with extreme high and low apoptosis phenotype - sort by high and low 

# Load in qPCR data 
QPCRDataQAed <- readxl::read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2018_Dermo_Challenge/Data/ANALYSIS FILES/ANALYSIS_CSVs_DAY7_BAD_2018_ANALYSIS/MASTER_DATA/MasterData_10.15.2018.XLSX",
                                                           col_names = c("qPCR_plate", "ID" , "Pconc_rep1", "log_pconc_rep_1", "Pconc_rep2", "log_pconc_rep2", "ave_log_pconc", "notes"),
                                                           sheet = "QPCRDataQAed", skip = 1)

# keep avg_log pconc and combine with assay data by first putting back together the family and oyster ID info that matches the load data 
Day7_Day50_2018_all_assays_bad_removed_pconc <- Day7_Day50_2018_all_assays_bad_removed %>% mutate(ID = paste(Family, Oyster_ID, sep ="-")) %>% 
                                                    left_join(., QPCRDataQAed[,c("ID", "ave_log_pconc")], by = "ID")

# Assess levels and which ones I should be viewing for comparison for which samples to select
levels(factor(Day7_Day50_2018_all_assays_bad_removed_pconc$Cell_type))
# levels of interest: "granular_live_apoptotic"  , "granular_dead_apoptotic", "granular_live_caspase_active" , "granular_dead_caspase_active"

# Subset for these levels 
subset_list <- c("granular_live_apoptotic"  , "granular_dead_apoptotic", "granular_live_caspase_active" , "granular_dead_caspase_active")
Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp <- Day7_Day50_2018_all_assays_bad_removed_pconc %>% 
  filter(Cell_type %in% subset_list, # subset for only the four types above
          Day == "7") # subset for only day 7
levels(factor(Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp$Cell_type))

Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_1 <- Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp %>% ungroup() %>% filter(Cell_type == "granular_live_apoptotic") %>% rename(granular_live_apoptotic_percent = Percent_of_this_plot) %>% 
  dplyr::select(Family, ID, Day, Treat, ave_log_pconc, granular_live_apoptotic_percent)
Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_2 <- Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp %>%  ungroup() %>% filter(Cell_type == "granular_dead_apoptotic") %>% rename(granular_dead_apoptotic_percent = Percent_of_this_plot) %>% 
  dplyr::select(Family, ID, Day, Treat, ave_log_pconc, granular_dead_apoptotic_percent)
Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_3 <- Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp %>%  ungroup() %>% filter(Cell_type == "granular_live_caspase_active") %>% rename(granular_live_caspase_percent = Percent_of_this_plot) %>% 
  dplyr::select(Family, ID, Day, Treat, ave_log_pconc, granular_live_caspase_percent)
Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_4 <- Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp %>%  ungroup() %>% filter(Cell_type == "granular_dead_caspase_active") %>% rename(granular_dead_caspase_percent = Percent_of_this_plot) %>% 
  dplyr::select(Family, ID, Day, Treat, ave_log_pconc, granular_dead_caspase_percent)

Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table <- left_join(Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_1, Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_2)
Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table <- left_join(Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table, Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_3)
Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table <- left_join(Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table, Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_4)

# create combined caspase and apoptosis columns
Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb <- Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table %>% group_by(ID) %>% mutate(apoptosis_combined = granular_live_apoptotic_percent +granular_dead_apoptotic_percent) %>%
  mutate(caspase_combined = granular_live_caspase_percent + granular_dead_caspase_percent) %>% dplyr::select(Family, ID, Day, Treat, ave_log_pconc, apoptosis_combined, caspase_combined) %>% 
  # sort by apoptosis levels to get a feeling of high and low 
  arrange(Family, Treat, desc(apoptosis_combined))
  
# Export this sheet to file to go through with Marta and Dina
write.csv(Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb, "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/Day7_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb.csv")

#### PLOTS AND STATISTICS 2018 ####

#Calculate mean and sd for all assays by family
Day7_Day50_2018_all_assays_bad_removed_summary <- Day7_Day50_2018_all_assays_bad_removed %>% group_by(Family, Treat, Day, Assay, Plot_number, Gate) %>% summarize(mean_percent=mean(Percent_of_this_plot), sd_percent =sd(Percent_of_this_plot))

### VIABILITY ASSAY Statistics and Plotting ####

# make sure group_by not masked by plyr
#detach(package:Rmisc)
#detach(package:plyr)

Day7_Day50_2018_all_assays_bad_removed_VI <- Day7_Day50_2018_all_assays_bad_removed %>% filter(Assay == "VI") 
unique(Day7_Day50_2018_all_assays_bad_removed_VI$Cell_type)
# "all_hemocytes"      "granular"           "agranular"          "all_live_hemocytes" "all_dead_hemocytes" "live_granular"      "live_agranular"

# Agranular and Granular cells Boxplot with significance bars, grouped by treatment color by gate(significant based on ANOVA)
Day7_Day50_2018_VIA_Percent_Agranular_Granular <- Day7_Day50_2018_all_assays_bad_removed_VI %>% filter(Gate =="E3" | Gate=="E1")

Day7_Day50_2018_VIA_Percent_Agranular_Granular_plot <- ggplot(data=Day7_Day50_2018_VIA_Percent_Agranular_Granular,
                                                             aes(y=Percent_of_this_plot, x=Treat, fill=Gate)) + geom_boxplot()+ 
  geom_point(position=position_dodge(width=0.75))+ 
  xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent of Granular and Agranular Viability Assay 2018") + 
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) + 
  facet_grid(Day~Family) + 
  scale_x_discrete(labels=c("D"="Dermo Injected", "NC"="Notched Control")) + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  scale_fill_manual(name="Cell Type", labels=c("Granular","Agranular"), values=c("#6c81d9","#50b47b")) 

# Agranular and Granular cells Boxplot with significance bars, grouped by GATE colored by treatment (no significance here based on ANOVA)
Day7_Day50_2018_VIA_Percent_Agranular_Granular_cell_type_plot <- ggplot(data=Day7_Day50_2018_VIA_Percent_Agranular_Granular,
                                                                       aes(y=Percent_of_this_plot, x=Gate, fill=Treat)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent of Granular and Agranular Hemocyte Events") +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12),
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) + 
  facet_grid(Day~Family) + 
  scale_x_discrete(labels=c("E1"="Granular", "E3"="Agranular")) + 
  scale_fill_manual(name="Cell Type", labels=c("Dermo Injected","Notched Control"), values=c("#7e78d4",
                                                                                             "#cd4272")) 

# VIA Granular vs. Agranular ANOVA with arcsine transformed

# Perform one wway AOV for each family and treatment and day E1 vs E3 in loop with TukeyHSD
Day7_Day50_2018_VIA_Percent_Agranular_Granular_AOV <- Day7_Day50_2018_VIA_Percent_Agranular_Granular %>%
  group_by(Family, Treat, Day) %>%
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Gate, data = .)))  %>%
  ungroup
  # only family B and J on Day 50 not significant 

Day7_Day50_2018_VIA_Percent_Agranular_Granular_AOV_TUKEY <-Day7_Day50_2018_VIA_Percent_Agranular_Granular %>%
  group_by(Family, Treat, Day) %>%
  do(broom::tidy(TukeyHSD(aov(Percent_of_this_plot_arcsine ~ Gate, data = .)))) %>%
  ungroup

# Anova looking at Day 
Day7_Day50_2018_VIA_Percent_Agranular_Granular_AOV_DAY <-Day7_Day50_2018_VIA_Percent_Agranular_Granular %>%
  group_by(Family, Treat, Gate) %>%
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Day, data = .))) %>%
  ungroup

# two way ANOVA
Day7_Day50_2018_VIA_Percent_Agranular_Granular_AOV_DAY_TREAT <-Day7_Day50_2018_VIA_Percent_Agranular_Granular %>%
  group_by(Family) %>%
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Day+Treat, data = .))) %>%
  ungroup

# one way anova granular and agranular hemocytes between treatments
Day7_Day50_2018_VIA_Percent_Agranular_Granular_TREAT_AOV <- Day7_Day50_2018_VIA_Percent_Agranular_Granular %>%
  group_by(Family, Day, Gate) %>%
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Treat, data = .)))  %>%
  ungroup
# Family L day 50 significant

#### Apoptosis Assay Statistics and Plotting ####

Day7_Day50_2018_APOP <- Day7_Day50_2018_all_assays_bad_removed %>% filter(Assay=="A")
unique(Day7_Day50_2018_APOP[,c("Gate","Cell_type")])
# Q12 = agranular, Q13 = granular
Day7_Day50_2018_APOP_Agranular_Granular <- Day7_Day50_2018_APOP %>% filter(Gate =="Q12-UL" | Gate =="Q12-UR" | Gate =="Q12-LL" | Gate =="Q12-LR"|
                                                                             Gate =="Q13-UL" | Gate =="Q13-UR" | Gate =="Q13-LL" | Gate =="Q13-LR")


# Analysis for all Granular and Agranular Combined from the quad plots
# Make APOP combined gate for agranular and granular separately 
Day7_Day50_2018_APOP_Granular_Apop_combined <- Day7_Day50_2018_APOP_Agranular_Granular  %>% filter(Gate == "Q13-LR" | Gate == "Q13-UR") %>% group_by(Oyster_ID, Treat,Family,Day) %>% summarise(Percent_of_this_plot = sum(Percent_of_this_plot))
Day7_Day50_2018_APOP_Agranular_Apop_combined <- Day7_Day50_2018_APOP_Agranular_Granular %>% filter(Gate == "Q12-LR" | Gate == "Q12-UR") %>% group_by(Oyster_ID, Treat,Family, Day) %>% summarise(Percent_of_this_plot = sum(Percent_of_this_plot))

# Add new gate name for each and cell type column 
Day7_Day50_2018_APOP_Granular_Apop_combined$Gate <- "apop_combined_granular"
Day7_Day50_2018_APOP_Agranular_Apop_combined$Gate <- "apop_combined_agranular"
Day7_Day50_2018_APOP_Granular_Apop_combined$Cell_type <- "apop_combined_granular"
Day7_Day50_2018_APOP_Agranular_Apop_combined$Cell_type <- "apop_combined_agranular"

# Combined data frames for each cell type
Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined <- rbind(Day7_Day50_2018_APOP_Granular_Apop_combined, Day7_Day50_2018_APOP_Agranular_Apop_combined)

# Add arcsine transformed data
Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined$Percent_of_this_plot_arcsine <- transf.arcsin(Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined$Percent_of_this_plot*0.01)

# Merge with original data frame
Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_full <- full_join(Day7_Day50_2018_APOP_Agranular_Granular, Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined, by =c("Gate", "Oyster_ID","Treat", "Family","Day","Percent_of_this_plot","Percent_of_this_plot_arcsine", "Cell_type"))

# Boxplot of all gates of agranular and granular cells with faceting by gate, quadrants next to each other
# Change levels of faceting so that box plot quadrants are next to each other with agranular first, starting with live and dead
unique(Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_full$Gate) # NOTE  
#  "apop_combined_granular"  "apop_combined_agranular" "Q12-UL"                  "Q12-UR"                  "Q12-LL"                 
#  "Q12-LR"                  "Q13-UR"                  "Q13-UL"                  "Q13-LL"                  "Q13-LR" 

Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_full$Gate <- factor(Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_full$Gate, 
                                                                   levels = c("Q12-UL", "Q13-UL","Q12-LL","Q13-LL","Q12-LR", "Q13-LR","Q12-UR", "Q13-UR","apop_combined_agranular","apop_combined_granular"))
# Make plot 
Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_plot <- ggplot(data=Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_full,
  aes(y=Percent_of_this_plot, x=Treat, fill=Gate)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Agranular and Granular Apoptosis Assay") + 
  facet_grid(Family~Day, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  scale_fill_manual(name="Cell Type", labels=c( "Dead Agranular","Dead Granular", "Live Agranular", "Live Granular",
                                                "Live Agranular Apoptotic",   "Live Granular Apoptotic",  "Dead Agranular Apoptotic",  "Dead Granular Apoptotic", "Combined Agranular Apoptotic","Combined Granular Apoptotic"), 
                    values = c("#cc57b4", "#88bf3b", "#aa4dce", "#56b464", "#5b2c90", 
                               "#c89832", "#5a6ee6", "#ca4e33", "#7e78d4", "#cd4272")) 
# plot without live and dead granular agranular
Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_full_subset <- Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_full %>% filter(Gate == "Q12-LR"| Gate ==  "Q13-LR"| Gate == "Q12-UR" |
                                                                         Gate == "Q13-UR")
Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_subset_plot <- ggplot(data=Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_full_subset,
  aes(y=Percent_of_this_plot, x=Treat, fill=Gate)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Agranular and Granular Apoptosis Assay") + 
  facet_grid(Family~Day, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  scale_fill_manual(name="Cell Type", labels=c("Live Agranular Apoptotic",   "Live Granular Apoptotic",  "Dead Agranular Apoptotic",  "Dead Granular Apoptotic"), 
                    values = c("#6b62b9","#8da54f","#b55d97","#bb5542")) 


# Plot only the apop combined 
Day7_Day50_2018_APOP_Granular_Agranular_Apop_just_combined_plot <- ggplot(data=Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined,
                                                                     aes(y=Percent_of_this_plot, x=Treat, fill=Gate)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Agranular and Granular Apoptotic Combined") + 
  facet_grid(Day~Family, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  scale_fill_manual(name="Cell Type", labels=c("Combined Agranular Apoptotic","Combined Granular Apoptotic"), 
                    values = c("#7e78d4", "#cd4272")) 

# Plot of just combined apoptotic granular in control and treated 
Day7_Day50_2018_APOP_Granular_Apop_just_combined_plot <- Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined %>% filter(Gate == "apop_combined_granular") %>%
        ggplot(aes(y=Percent_of_this_plot, x=Treat, fill=Treat)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Granular Apoptotic Combined") + 
  facet_grid(Day~Family, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  scale_fill_manual(name="Cell Type", labels=c("Control","Dermo Challenge"), 
                    values = c("#7e78d4", "#cd4272")) 
ggsave(plot = Day7_Day50_2018_APOP_Granular_Apop_just_combined_plot, filename  = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2018_Dermo_Challenge/Data/ANALYSIS FILES/APOPTOSIS_ANOVA/Figures/Day7_Day50_2018_APOP_Granular_Apop_just_combined_plot.tiff",
       device = "tiff")

# ANOVA
# Apop combined granular vs apop combined agranular within each treatment
Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_AOV <- Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined %>%
  group_by(Family, Treat, Day) %>%
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Gate, data = .)))  %>%
  ungroup

Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_AOV_TUKEY <- Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined %>%
  group_by(Family, Treat, Day) %>%
  do(broom::tidy(TukeyHSD(aov(Percent_of_this_plot_arcsine ~ Gate, data = .)))) %>%
  ungroup

# Apop combined granular vs apop combined agranular within each treatment between Days
Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_DAY_AOV <- Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined %>%
  group_by(Family, Treat, Gate) %>%
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Day, data = .)))  %>%
  ungroup
# B dermo granular, D granular control and dermo, E Dermo Granular (not control), J dermo granular (not control)
  
# one way anova granular and agranular hemocytes between treatments
Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_TREAT_AOV <-Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined %>%
  group_by(Family, Day, Gate) %>%
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Treat, data = .)))  %>%
  ungroup
# A day 7 aganular significant
Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_TREAT_AOV %>% filter(p.value <= 0.05)
#Family Day   Gate                    term     df  sumsq meansq statistic p.value
#<chr>  <chr> <chr>                   <chr> <dbl>  <dbl>  <dbl>     <dbl>   <dbl>
#  1 A      7     apop_combined_agranular Treat     1 0.0293 0.0293      7.42  0.0234

Day7_Day50_2018_APOP_Granular_Agranular_all_TREAT_AOV <-Day7_Day50_2018_APOP_Granular_Agranular_Apop_combined_full %>%
  group_by(Family, Day, Gate) %>%
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Treat, data = .)))  %>%
  ungroup
# significant include:
  #D	50	Q13-LL	Treat	1	6.11E-02	6.11E-02	7.99E+00	0.01984361
  #A	7	apop_combined_agranular	Treat	1	2.93E-02	2.93E-02	7.42E+00	0.02344723
  #A	7	Q12-LL	Treat	1	2.72E-02	2.72E-02	7.15E+00	0.02542913  #agranular
  #A	7	Q12-LR	Treat	1	2.47E-02	2.47E-02	6.96E+00	0.02700916 #agranular
  #E	7	Q12-UL	Treat	1	2.99E-02	2.99E-02	6.71E+00	0.03592898 #agranular
  #D	7	Q12-UR	Treat	1	5.15E-03	5.15E-03	5.79E+00	0.05289778 #agranular
  #B	7	Q13-UL	Treat	1	7.32E-02	7.32E-02	4.94E+00	0.05696606 # granular


## Combined granular combined control vs. treated for each family


#### Caspase Assay Statistics ####

Day7_Day50_2018_CASP <- Day7_Day50_2018_all_assays_bad_removed %>% filter(Assay=="C")
unique(Day7_Day50_2018_CASP[,c("Gate","Cell_type")])
# Q2= agranular, Q6= granular
Day7_Day50_2018_CASP_Agranular_Granular <- Day7_Day50_2018_CASP %>% filter(Gate =="Q2-UL" | Gate =="Q2-UR" | Gate =="Q2-LL" | Gate =="Q2-LR"|
                                                                             Gate =="Q6-UL" | Gate =="Q6-UR" | Gate =="Q6-LL" | Gate =="Q6-LR")

# Analysis for all Granular and Agranular Combined from the quad plots
# Make APOP combined gate for agranular and granular separately 
Day7_Day50_2018_CASP_Granular_casp_combined <-  Day7_Day50_2018_CASP_Agranular_Granular  %>% filter(Gate == "Q6-LR" | Gate == "Q6-UR") %>% group_by(Oyster_ID, Treat,Family,Day) %>% summarise(Percent_of_this_plot = sum(Percent_of_this_plot))
Day7_Day50_2018_CASP_Agranular_casp_combined <- Day7_Day50_2018_CASP_Agranular_Granular %>% filter(Gate == "Q2-LR" | Gate == "Q2-UR") %>% group_by(Oyster_ID, Treat,Family, Day) %>% summarise(Percent_of_this_plot = sum(Percent_of_this_plot))

# Add new gate name for each 
Day7_Day50_2018_CASP_Granular_casp_combined$Gate <- "casp_active_combined_granular"
Day7_Day50_2018_CASP_Agranular_casp_combined$Gate <- "casp_active_combined_agranular"

# Combined data frames for each cell type
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined <- rbind(Day7_Day50_2018_CASP_Granular_casp_combined, Day7_Day50_2018_CASP_Agranular_casp_combined)

# Add arcsine transformed data
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined$Percent_of_this_plot_arcsine <- transf.arcsin(Day7_Day50_2018_CASP_Granular_Agranular_casp_combined$Percent_of_this_plot*0.01)

# Merge with original data frame
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_full <- full_join(Day7_Day50_2018_CASP_Granular_Agranular_casp_combined, Day7_Day50_2018_CASP_Agranular_Granular, by =c("Gate", "Oyster_ID","Treat", "Family","Day","Percent_of_this_plot","Percent_of_this_plot_arcsine"))

# Boxplot of all gates of agranular and granular cells with faceting by gate, quadrants next to each other
# Change levels of faceting so that box plot quadrants are next to each other with agranular first, starting with live and dead
unique(Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_full$Gate) # NOTE  
# "casp_active_combined_granular"  "casp_active_combined_agranular" "Q2-UL"                          "Q2-UR"                         
# "Q2-LL"                          "Q2-LR"                          "Q6-UL"                          "Q6-UR"                         
# "Q6-LL"                          "Q6-LR" 

Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_full$Gate <- factor(Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_full$Gate, 
                                                                          levels = c("Q2-UL", "Q6-UL","Q2-LL","Q6-LL","Q2-LR", "Q6-LR","Q2-UR", "Q6-UR","casp_active_combined_agranular" ,"casp_active_combined_granular" ))
# Make plot 
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_full_plot <- ggplot(data=Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_full,
                                                                     aes(y=Percent_of_this_plot, x=Treat, fill=Gate)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Agranular and Granular Caspase Assay") + 
  facet_grid(Family~Day, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  scale_fill_manual(name="Cell Type", labels=c( "Dead Agranular","Dead Granular", "Live Agranular", "Live Granular",
                                                "Live Agranular Casp. 3/7 Active",   "Live Granular Casp. 3/7 Active",  "Dead Agranular Casp. 3/7 Active",  "Dead Granular Casp. 3/7 Active", "Combined Agranular Casp. 3/7 Active","Combined Granular Casp. 3/7 Active"), 
                    values = c("#cc57b4", "#88bf3b", "#aa4dce", "#56b464", "#5b2c90", 
                               "#c89832", "#5a6ee6", "#ca4e33", "#7e78d4", "#cd4272")) 

# Plot only the apop combined 
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_just_combined_plot <- ggplot(data=Day7_Day50_2018_CASP_Granular_Agranular_casp_combined,
                                                                          aes(y=Percent_of_this_plot, x=Treat, fill=Gate)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Agranular and Granular Caspase 3/7 Active Combined") + 
  facet_grid(Day~Family, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  scale_fill_manual(name="Cell Type", labels=c("Combined Agranular Casp. 3/7 Active","Combined Granular Casp. 3/7 Active"), 
                    values = c("#7e78d4", "#cd4272")) 

# Plot only granular at Day 7, families A and J removed because no control
Day7_2018_CASP_Granular_casp_combined <- Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>% filter(Gate =="casp_active_combined_granular" & Day == "7" & Family != "A" & Family !="J")
Day7_2018_CASP_Granular_casp_combined_just_combined_plot <- ggplot(data=Day7_2018_CASP_Granular_casp_combined,
                                                                                   aes(y=Percent_of_this_plot, x=Treat, fill=Treat)) + geom_boxplot()+ geom_point(position=position_dodge(width=0.75)) + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Granular Caspase 3/7 Active Combined 2018") + 
  facet_grid(Day~Family, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  scale_fill_manual(name="Treatment", labels=c("Control","Dermo Injected"), 
                    values = c("#7e78d4", "#cd4272")) 

#export this plot
ggsave("Day7_2018_CASP_Granular_casp_combined_just_combined_plot.tiff", plot = Day7_2018_CASP_Granular_casp_combined_just_combined_plot ,
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES", device = "tiff")

# ANOVA
# Casp active combined granular vs apop combined agranular within each treatment
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_AOV <- Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>%
  group_by(Family, Treat, Day) %>%
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Gate, data = .)))  %>%
  ungroup
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_AOV %>% filter(p.value <= 0.05)

Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_AOV_TUKEY <- Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>%
  group_by(Family, Treat, Day) %>%
  do(broom::tidy(TukeyHSD(aov(Percent_of_this_plot_arcsine ~ Gate, data = .)))) %>%
  ungroup

# one way anova granular and agranular hemocytes between treatments
# check levels
lapply(Day7_Day50_2018_CASP_Granular_Agranular_casp_combined[c("Family", "Gate", "Day")], unique)
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>% filter(Family == "A") %>% View()
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>% filter(Family == "B") %>% View() # only has two controls
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>% filter(Family == "D") %>% View() # has one control
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>% filter(Family == "E") %>% View()
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>% filter(Family == "J") %>% View() # has no controls Day 7 
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>% filter(Family == "L") %>% View() # has one control

Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_TREAT_A_AOV <- Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>% # won't work 
  group_by(Family, Gate, Day) %>% filter(Family == "A") %>% 
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Treat, data = .)))  %>%
  ungroup
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_TREAT_B_AOV <- Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>%
    group_by(Family, Gate, Day) %>% filter(Family == "B") %>% 
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Treat, data = .)))  %>%
  ungroup 
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_TREAT_B_AOV %>% filter(p.value <= 0.05) # NONE significant

Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_TREAT_D_AOV <- Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>%
  group_by(Family, Gate, Day) %>% filter(Family == "D") %>% 
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Treat, data = .)))  %>%
  ungroup
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_TREAT_D_AOV %>% filter(p.value <= 0.05) # none significant

Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_TREAT_E_AOV <- Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>%
  group_by(Family, Gate, Day) %>% filter(Family == "E") %>% 
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Treat, data = .)))  %>%
  ungroup
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_TREAT_E_AOV %>% filter(p.value <= 0.05) # none significant

Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_TREAT_J_AOV <- Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>% # not enough samples
  group_by(Family, Gate, Day) %>% filter(Family == "J") %>% 
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Treat, data = .)))  %>%
  ungroup


Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_TREAT_L_AOV <- Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>% # not enough samples
  group_by(Family, Gate, Day) %>% filter(Family == "L") %>% 
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Treat, data = .)))  %>%
  ungroup
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_TREAT_L_AOV %>% filter(p.value <= 0.05) # none significant

# none are significant

# Anova of changes in gate between days
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_DAY_B_AOV <- Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>%
  group_by(Family, Gate, Treat) %>% filter(Family == "B") %>% 
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Day, data = .)))  %>%
  ungroup
  # control Dermo for both agranular and granular are significantly different

Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_DAY_D_AOV <- Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>%
  group_by(Family, Gate, Treat) %>% filter(Family == "D") %>% 
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Day, data = .)))  %>%
  ungroup
# granular casp active control different by day 

Day7_Day50_2018_CASP_Granular_Agranular_casp_combined_DAY_E_AOV <- Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>%
  group_by(Family, Gate, Treat) %>% filter(Family == "E") %>% 
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Day, data = .)))  %>%
  ungroup

# Anova changes in granular hemocyte apoptosis between families in the treated group
Day7_Day50_2018_CASP_Granular_Agranular_casp_combined %>%
  filter(Gate == "casp_active_combined_granular" & Treat == "Dermo") %>%
 aov(Percent_of_this_plot_arcsine ~ Family, data = .)
#Call:
#  aov(formula = Percent_of_this_plot_arcsine ~ Family, data = .)
#
#Terms:
#  Family Residuals
#Sum of Squares  0.223390  3.471664
#Deg. of Freedom        5        67
#
#Residual standard error: 0.227631
#Estimated effects may be unbalanced


#### COMBINED 2018 2019 PLOTS ####

# GOAL: PLOT ALL SAMPLES FOR EACH PARAMETER IN EACH YEAR, BOTH CONTROL AND CHALLENGE TO LOOK FOR TRENDS 

# Combine together all dataframes with phenotype data for both years
colnames(Day7_2019_VIA_bad_removed_no_outlier)
colnames(Day7_2019_APOP_bad_removed_no_outlier)
colnames(Day7_2019_CASP_bad_removed_no_outlier)
colnames(Day7_2019_LMP_bad_removed_no_outlier)
  # all the four above are the same columns
colnames(Day7_Day50_2018_all_assays_bad_removed) # 
      # [1] "Family"                       "ID"                           "Assay"                        "Plot_number"                  "Gate"                        
      # [6] "Percent_of_this_plot"         "Counts"                       "Day"                          "Cell_type"                    "Flow_Code"                   
      # [11] "Treat"                        "Oyster_ID"                    "Percent_of_this_plot_arcsine"
levels(factor(Day7_Day50_2018_all_assays_bad_removed$Assay)) #  "A"  "C"  "VI"

# join 2019 dataframes and add in family column - call it S for susceptible
Day7_2019_VIA_bad_removed_no_outlier_all_assays <- rbind(Day7_2019_VIA_bad_removed_no_outlier, Day7_2019_APOP_bad_removed_no_outlier, 
                                                        Day7_2019_CASP_bad_removed_no_outlier,Day7_2019_LMP_bad_removed_no_outlier) %>%
  mutate(Family = "S", Day = "7")
colnames(Day7_2019_VIA_bad_removed_no_outlier_all_assays) #  
  # [1] "ID"                           "Treat"                        "Assay"                        "Plot_number"                  "Gate"                        
  # [6] "Percent_of_this_plot"         "Counts"                       "Cell_type"                    "Percent_of_this_plot_arcsine" "Family"                      
  # [11] "Day"

levels(factor(Day7_2019_VIA_bad_removed_no_outlier_all_assays$Assay)) # "AP"  "C"   "LMP" "VIA"
Day7_2019_VIA_bad_removed_no_outlier_all_assays <- Day7_2019_VIA_bad_removed_no_outlier_all_assays %>% mutate(Assay = case_when(
  Assay == "AP" ~ "A",
  Assay == "VIA" ~ "VI",
  TRUE ~ Assay),
  Treat = case_when(
    Treat == "D" ~ "Dermo",
    Treat == "NC" ~ "control",
    TRUE ~ Treat
  ))
levels(factor(Day7_2019_VIA_bad_removed_no_outlier_all_assays$Assay))

# Join all assays
Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays <- rbind(Day7_Day50_2018_all_assays_bad_removed,Day7_2019_VIA_bad_removed_no_outlier_all_assays)

## 2018 2019 VIA plot
Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_VIA_Fam <- Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays %>% 
  filter(Assay == "VI" & Cell_type == "granular") %>% ungroup() %>%
  ggplot(aes(y=Percent_of_this_plot, x=Treat)) + geom_point() + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Total (live and dead) Granular Hemocytes") + 
  facet_grid(Day~Family, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12))

ggsave(plot = Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_VIA_Fam, 
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES/",
       filename = "Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_VIA_Fam.tiff",
       width = 5, height = 4,
       device = "tiff")

Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_Day_allfam <- Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays %>% 
  filter(Assay == "VI" & Cell_type == "granular") %>% ungroup() %>%
  ggplot(aes(y=Percent_of_this_plot, x=Treat)) + geom_point() + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Total (live and dead) Granular Hemocytes") + 
  facet_grid(Day~., scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12))

ggsave(plot = Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_Day_allfam, 
       width = 5, height = 4,
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES/",
       filename = "Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_Day_allfam.tiff",
       device = "tiff")

## 2018 2019 APOP PLOT 
Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_APOP_Fam <- Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays %>% 
  filter(Assay == "A" & 
           # Cell_type == "granular_live_apoptotic" | 
           Cell_type == "granular_dead_apoptotic") %>% ungroup() %>% 
  ggplot(aes(y=Percent_of_this_plot, x=Treat, color = Cell_type)) + geom_point() + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent live and dead Granular apoptotic Hemocytes") + 
  facet_grid(Day~Family, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=10))

ggsave(plot = Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_APOP_Fam, 
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES/",
       filename = "Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_APOP_Fam.tiff",
       width = 6, height = 4,
       device = "tiff")

Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_Day_APOP_allfam <- Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays %>% 
  filter(Assay == "A" & 
           #Cell_type == "granular_live_apoptotic" | 
           Cell_type == "granular_dead_apoptotic") %>% ungroup() %>% 
  ggplot(aes(y=Percent_of_this_plot, x=Treat, color = Cell_type)) + geom_point() + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent live and dead Granular apoptotic Hemocytes") + 
  facet_grid(Day~., scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=10))

ggsave(plot = Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_Day_APOP_allfam, 
       width = 5, height = 4,
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES/",
       filename = "Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_Day_APOP_allfam.tiff",
       device = "tiff")

# add together the granular cell percentages for live and dead 
Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_gran_dead <- Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays %>% 
  filter(Assay == "A" & 
            # Cell_type == "granular_live_apoptotic" | 
           Cell_type == "granular_dead_apoptotic") %>% rename(Percent_of_this_plot_dead = Percent_of_this_plot)
Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_gran_live <- Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays %>% 
  filter(Assay == "A" & 
            Cell_type == "granular_live_apoptotic" ) %>% rename(Percent_of_this_plot_live = Percent_of_this_plot)

Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_gran_all <- left_join(Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_gran_dead[,c("Family","ID","Treat","Oyster_ID","Plot_number","Percent_of_this_plot_dead", "Day")],
                                                                                      Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_gran_live[,c("Family","ID","Treat","Oyster_ID","Plot_number","Percent_of_this_plot_live", "Day")] ,
                                                                                      by = c("Family","ID","Treat","Oyster_ID","Plot_number", "Day")) %>% mutate(Percent_of_this_plot = Percent_of_this_plot_dead + Percent_of_this_plot_live)

Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_APOP_Fam_gran_all <- Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_gran_all %>%
  ggplot(aes(y=Percent_of_this_plot, x=Treat, color = Treat)) + geom_point()  + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Total Apoptotic Granular Hemocytes") + 
  facet_grid(Day~Family, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=10))

ggsave(plot = Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_APOP_Fam_gran_all, 
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES/",
       filename = "Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_APOP_Fam_gran_all.tiff",
       width = 6, height = 4,
       device = "tiff")


## 2018 2019 CASPASE PLOT 
Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_CASP_Fam <- Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays %>% 
  filter(Assay == "C" & Cell_type == "granular_live_caspase_active" | Cell_type == "granular_dead_caspase_active") %>% ungroup() %>% 
  ggplot(aes(y=Percent_of_this_plot, x=Treat, color = Cell_type)) + geom_point() + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent live and dead Granular caspase active Hemocytes") + 
  facet_grid(Day~Family, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=10))

ggsave(plot = Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_CASP_Fam, 
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES/",
       filename = "Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_CASP_Fam.tiff",
       width = 6, height = 4,
       device = "tiff")

Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_Day_CASP_allfam <- Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays %>% 
  filter(Assay == "C" & Cell_type == "granular_live_caspase_active" | Cell_type == "granular_dead_caspase_active") %>% ungroup() %>% 
  ggplot(aes(y=Percent_of_this_plot, x=Treat, color = Cell_type)) + geom_point() + xlab("Treatment") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent live and dead Granular caspase active Hemocytes") + 
  facet_grid(Day~., scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(axis.text.x = element_text(size=12, hjust=1, angle=50)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=10))

ggsave(plot = Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_Day_CASP_allfam, 
       width = 5, height = 4,
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES/",
       filename = "Day7_Day50_2018_Day7_2019_VIA_bad_removed_no_outlier_all_assays_Day_CASP_allfam.tiff",
       device = "tiff")


## 2019 LYSOSOME PLOT 
Day7_2019_LMP_Percent_Granular_Agranular_Hemocytes_no_outlier_points <- Day7_2019_LMP_ALL_Granular_Agranular_no_outlier %>% 
  ggplot(aes(y=Percent_of_this_plot, x=Gate, color=Treat)) + geom_point(position=position_dodge(width=0.75)) + xlab("Quadrant") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Agranular and Granular Intact Lysosome FITC vs. PI Quad Plot") + 
  facet_grid(.~Gate, scales="free") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, hjust=1, angle=70)) +
  theme(legend.text = element_text(size=12)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_color_manual(name="Treatment", labels=c("Dermo Injected","Notched Control"), values=c("#7e78d4","#cd4272")) +
  scale_x_discrete(labels=c("Q3-UL"="agranular_lysosome_rupture",
                            "Q2-UL"="granular_lysosome_rupture", 
                            "Q3-LL"="agranular_LMP_onset",
                            "Q2-LL"="granular_LMP_onset", 
                            "Q3-UR"="agranular_dead_other_means", 
                            "Q2-UR"="granular_dead_other_means",
                            "Q3-LR"="agranular_lysosome_live",
                            "Q2-LR"="granular_lysosome_live"))

ggsave(plot = Day7_2019_LMP_Percent_Granular_Agranular_Hemocytes_no_outlier_points, 
       width = 7, height = 6,
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES/",
       filename = "Day7_2019_LMP_Percent_Granular_Agranular_Hemocytes_no_outlier_points.tiff",
       device = "tiff")


##### 2019 IN VIVO 1HR HEMOCYTE ASSAY ANALYSIS #####

## LOAD ANALYSIS CSVs

## VIABILITY ASSAY

# Make new header column
MOI_1hr_2019_VIA_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/2019_Hemocyte_Assays/1HR_MOI_ASSAYS/ANALYSIS_FILES/FORMATTED_CSVS/MOI_1hr_12_09_19_VIA_header_formatted2.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
MOI_1hr_2019_VIA <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/2019_Hemocyte_Assays/1HR_MOI_ASSAYS/ANALYSIS_FILES/FORMATTED_CSVS/MOI_1hr_12_09_19_VIA_header_formatted2.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(MOI_1hr_2019_VIA_nms)

# Split column 1 by space and remove
MOI_1hr_2019_VIA <- MOI_1hr_2019_VIA %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
MOI_1hr_2019_VIA <- MOI_1hr_2019_VIA[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
MOI_1hr_2019_VIA_percent <- MOI_1hr_2019_VIA[,c(1:3,5,7,9,11,13,15,17,19,21,23)]
MOI_1hr_2019_VIA_counts <- MOI_1hr_2019_VIA[,c(1:4,6,8,10,12,14,16,18,20,22)]

# Gather count and percent columns separately
MOI_1hr_2019_VIA_percent <- MOI_1hr_2019_VIA_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:13))
MOI_1hr_2019_VIA_counts <-  MOI_1hr_2019_VIA_counts %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:13))

MOI_1hr_2019_VIA_percent <- MOI_1hr_2019_VIA_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
MOI_1hr_2019_VIA_counts <-  MOI_1hr_2019_VIA_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

MOI_1hr_2019_VIA_percent <- MOI_1hr_2019_VIA_percent[,-c(4,8)]
MOI_1hr_2019_VIA_counts <-  MOI_1hr_2019_VIA_counts[,-c(4,8)]

# Full join together so all columns preserved
MOI_1hr_2019_VIA_join <- full_join(MOI_1hr_2019_VIA_percent, MOI_1hr_2019_VIA_counts, by = c("ID","Treat","Assay","Plot_number","Gate", "Channel"))

# Add in cell type column based on plot number and gate
MOI_1hr_2019_VIA_cell_type <- data.frame(Plot_number = c("8","8", "19","19", "20","20", "21","21", "22","22"), Gate= c("P4","P3","V6-L","V6-R", "V7-L","V7-R", "V10-L", "V10-R","V11-L","V11-R"), Cell_type=c(
  "granular","agranular", "agranular","live_agranular", "granular","live_granular", 
  "agranular","dead_agranular","granular","dead_granular"))
# Join cell type
MOI_1hr_2019_VIA_join <- left_join(MOI_1hr_2019_VIA_join, MOI_1hr_2019_VIA_cell_type, by=c("Plot_number","Gate"))
unique_MOI_1hr_2019_VIA_join_gate <- unique(MOI_1hr_2019_VIA_join [,c("Plot_number","Gate")])

## APOPTOSIS ASSAY
## HK and FSW samples correctly switched to be the right labels in excel 
# Make new header column
MOI_1hr_2019_Apop_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/2019_Hemocyte_Assays/1HR_MOI_ASSAYS/ANALYSIS_FILES/FORMATTED_CSVS/MOI_1hr_12_09_19_APOP_header_formatted2.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
MOI_1hr_2019_APOP <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/2019_Hemocyte_Assays/1HR_MOI_ASSAYS/ANALYSIS_FILES/FORMATTED_CSVS/MOI_1hr_12_09_19_APOP_header_formatted2.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(MOI_1hr_2019_Apop_nms)

# Split column 1 by space and remove
MOI_1hr_2019_APOP <- MOI_1hr_2019_APOP %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
MOI_1hr_2019_APOP <- MOI_1hr_2019_APOP[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
MOI_1hr_2019_APOP_percent <- MOI_1hr_2019_APOP[,c(1:3,5,7,9,11,13,15,17,19,21, 23,25,27)]
MOI_1hr_2019_APOP_counts <-  MOI_1hr_2019_APOP[,c(1:4,6,8,10,12,14,16,18,20,22,24,26)]

# Gather count and percent columns separately
MOI_1hr_2019_APOP_percent <- MOI_1hr_2019_APOP_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:15))
MOI_1hr_2019_APOP_counts <-  MOI_1hr_2019_APOP_counts %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:15))

MOI_1hr_2019_APOP_percent <- MOI_1hr_2019_APOP_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
MOI_1hr_2019_APOP_counts <-  MOI_1hr_2019_APOP_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

MOI_1hr_2019_APOP_percent <- MOI_1hr_2019_APOP_percent[,-c(4,8)]
MOI_1hr_2019_APOP_counts <-  MOI_1hr_2019_APOP_counts[,-c(4,8)]

# Full join together so all columns preserved
MOI_1hr_2019_APOP_join <- full_join(MOI_1hr_2019_APOP_percent, MOI_1hr_2019_APOP_counts, by = c("ID","Treat","Assay","Plot_number","Gate","Channel"))

# Add in cell type column based on plot number and gate
MOI_1hr_2019_APOP_cell_type <- data.frame(Plot_number = c("2","2","5","5","8","8","4","4","4","4","7","7","7","7"), 
                             Gate= c("V2-L","V2-R","V1-L","V1-R","P3","P4","Q1-UL","Q1-UR","Q1-LL","Q1-LR","Q2-UL","Q2-UR","Q2-LL","Q2-LR"), 
                             Cell_type=c("non_apop_cells", "apop_cells", "dead_non_stained_perkinsus", "live_perkinsus", "agranular", 
                                         "granular", "apoptotic_granular_no_parasite", "apoptotic_granular_parasite", "unstained_granular", "live_perkinsus", 
                                         "apoptotic_agranular_no_parasite", "apoptotic_agranular_parasite", "unstained_agranular", "live_perkinsus"))

# Join cell type
MOI_1hr_2019_APOP_join <- left_join(MOI_1hr_2019_APOP_join, MOI_1hr_2019_APOP_cell_type, by=c("Plot_number","Gate"))
unique_MOI_1hr_2019_APOP_join <- unique(MOI_1hr_2019_APOP_join[,c("Plot_number","Gate")])

## PI ASSAY
# Make new header column
MOI_1hr_2019_PI_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/2019_Hemocyte_Assays/1HR_MOI_ASSAYS/ANALYSIS_FILES/FORMATTED_CSVS/MOI_1hr_12_09_19_PI_header_formatted2.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
MOI_1hr_2019_PI <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/2019_Hemocyte_Assays/1HR_MOI_ASSAYS/ANALYSIS_FILES/FORMATTED_CSVS/MOI_1hr_12_09_19_PI_header_formatted2.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(MOI_1hr_2019_PI_nms)

# Split column 1 by space and remove
MOI_1hr_2019_PI <- MOI_1hr_2019_PI %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
MOI_1hr_2019_PI <- MOI_1hr_2019_PI[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
MOI_1hr_2019_PI_percent <- MOI_1hr_2019_PI[,c(1:3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)]
MOI_1hr_2019_PI_counts <-  MOI_1hr_2019_PI[,c(1:4,6,8,10,12,14,16,18,20,22,24,26,28,30)]

# Gather count and percent columns separately
MOI_1hr_2019_PI_percent <- MOI_1hr_2019_PI_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:17))
MOI_1hr_2019_PI_counts <-  MOI_1hr_2019_PI_counts %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:17))

MOI_1hr_2019_PI_percent <- MOI_1hr_2019_PI_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
MOI_1hr_2019_PI_counts <-  MOI_1hr_2019_PI_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

MOI_1hr_2019_PI_percent <- MOI_1hr_2019_PI_percent[,-c(4,8)]
MOI_1hr_2019_PI_counts <-  MOI_1hr_2019_PI_counts[,-c(4,8)]

# Full join together so all columns preserved
MOI_1hr_2019_PI_join <- full_join(MOI_1hr_2019_PI_percent, MOI_1hr_2019_PI_counts, by = c("ID","Treat","Assay","Plot_number","Gate","Channel"))

# Add in cell type column based on plot number and gate
MOI_1hr_2019_PI_cell_type <- data.frame(Plot_number = c("2","2","5","5","8","8","4","4","4","4","7","7","7","7"), 
                             Gate= c("V1-L","V1-R","V3-L","V3-R","P3","P4","Q2-UL","Q2-UR", "Q2-LL", "Q2-LR", "Q4-UL","Q4-UR","Q4-LL","Q4-LR"), 
                             Cell_type=c("dead_non_stained_perkinsus", "live_perkinsus", "live_cells", "dead_cells", "agranular", 
                                         "granular", "dead_granular_no_parasite", "dead_granular_parasite", "unstained_granular", "live_perkinsus", 
                                         "dead_agranular_no_parasite", "dead_agranular_parasite", "unstained_agranular", "live_perkinsus"))

# Join cell type
MOI_1hr_2019_PI_join <- left_join(MOI_1hr_2019_PI_join, MOI_1hr_2019_PI_cell_type, by=c("Plot_number","Gate"))
unique_MOI_1hr_2019_PI_join <- unique(MOI_1hr_2019_PI_join[,c("Plot_number","Gate")])

## JC1 ASSAY
# Make new header column
MOI_1hr_2019_JC1_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/2019_Hemocyte_Assays/1HR_MOI_ASSAYS/ANALYSIS_FILES/FORMATTED_CSVS/MOI_1hr_12_09_19_JC1_header_formatted2.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
MOI_1hr_2019_JC1 <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2019 Dermo Challenge Experiments/DATA/2019_Hemocyte_Assays/1HR_MOI_ASSAYS/ANALYSIS_FILES/FORMATTED_CSVS/MOI_1hr_12_09_19_JC1_header_formatted2.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(MOI_1hr_2019_JC1_nms)

# Split column 1 by space and remove
MOI_1hr_2019_JC1 <- MOI_1hr_2019_JC1 %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
MOI_1hr_2019_JC1 <- MOI_1hr_2019_JC1[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
MOI_1hr_2019_JC1_percent <- MOI_1hr_2019_JC1[,c(1:3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67)]
MOI_1hr_2019_JC1_counts <-  MOI_1hr_2019_JC1[,c(1:4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66)]

# Gather count and percent columns separately
MOI_1hr_2019_JC1_percent <- MOI_1hr_2019_JC1_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:35))
MOI_1hr_2019_JC1_counts <-  MOI_1hr_2019_JC1_counts %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:35))

MOI_1hr_2019_JC1_percent <- MOI_1hr_2019_JC1_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
MOI_1hr_2019_JC1_counts <-  MOI_1hr_2019_JC1_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

MOI_1hr_2019_JC1_percent <- MOI_1hr_2019_JC1_percent[,-c(4,8)]
MOI_1hr_2019_JC1_counts <-  MOI_1hr_2019_JC1_counts[,-c(4,8)]

# Full join together so all columns preserved
MOI_1hr_2019_JC1_join <- full_join(MOI_1hr_2019_JC1_percent, MOI_1hr_2019_JC1_counts, by = c("ID","Treat","Assay","Plot_number","Gate", "Channel"))

# Add in cell type column based on plot number and gate
MOI_1hr_2019_JC1_cell_type <- data.frame(Plot_number = c("2","2","5", "5","8","8","9","9", "14","14","14","14","16","16","16","16","17","17","17","17","18","18","18","18",
                                                         "19","19","19","19","20","20","20","20"),     
                            Gate= c("V3-L","V3-R","V1-L","V1-R", "P3","P4","V2-L","V2-R", "Q17-UL","Q17-UR","Q17-LL","Q17-LR","Q11-UL","Q11-UR","Q11-LL","Q11-LR",
                                    "Q22-UL","Q22-UR","Q22-LL","Q22-LR","Q28-UL","Q28-UR","Q28-LL","Q28-LR","Q25-UL","Q25-UR","Q25-LL","Q25-LR",
                                    "Q26-UL","Q26-UR","Q26-LL","Q26-LR"),
                            Cell_type=c("normal_hemocytes", "potentially_apoptotic_hemocytes","dead_non_stained_perkinsus","live_perkinsus", "agranular","granular", "potentially_apoptotic_cells_non_stained", "normal_hemocytes",
                                        "normal_agranular","normal_agranular", "apoptotic_agranular_or_parasite", "unstained_agranular_parasite",
                                        "apoptotic_agranular_no_parasite", "apoptotic_agranular_live_parasite","non_apoptotic_live_parasite", "dead_parasite_unstained_agranular",
                                        "normal_granular","normal_granular", "apoptotic_granular_or_parasite", "unstained_granular_parasite",
                                        "apoptotic_granular_no_parasite", "apoptotic_granular_live_parasite","non_apoptotic_live_parasite", "dead_parasite_unstained_granular",
                                        "normal_agranular_no_parasite", "normal_agranular_live_parasite","normal_live_parasite", "dead_parasite_unstained_agranular",
                                        "normal_granular_no_parasite", "normal_granular_live_parasite","normal_live_parasite", "dead_parasite_unstained_granular"))
# Join cell type
MOI_1hr_2019_JC1_join <- left_join(MOI_1hr_2019_JC1_join, MOI_1hr_2019_JC1_cell_type, by=c("Plot_number","Gate"))

#### MOI 1hr 2019 Arcsine Transformation of percentages ####

# Make list of dataframes
MOI_1hr_all_list <- list(MOI_1hr_2019_VIA_join=MOI_1hr_2019_VIA_join, 
                         MOI_1hr_2019_APOP_join=MOI_1hr_2019_APOP_join, 
                         MOI_1hr_2019_PI_join=MOI_1hr_2019_PI_join, 
                         MOI_1hr_2019_JC1_join=MOI_1hr_2019_JC1_join)

MOI_1hr_2019_VIA_join
MOI_1hr_2019_APOP_join
MOI_1hr_2019_PI_join
MOI_1hr_2019_JC1_join

# Make new column and perform arcsine
MOI_1hr_2019_VIA_join$Percent_of_this_plot_arcsine <-  transf.arcsin(MOI_1hr_2019_VIA_join$Percent_of_this_plot)
MOI_1hr_2019_APOP_join$Percent_of_this_plot_arcsine <- transf.arcsin(MOI_1hr_2019_APOP_join$Percent_of_this_plot)
MOI_1hr_2019_PI_join$Percent_of_this_plot_arcsine <-   transf.arcsin(MOI_1hr_2019_PI_join$Percent_of_this_plot)
MOI_1hr_2019_JC1_join$Percent_of_this_plot_arcsine <-  transf.arcsin(MOI_1hr_2019_JC1_join$Percent_of_this_plot)

# Change percent column to percent
MOI_1hr_2019_VIA_join$Percent_of_this_plot <- MOI_1hr_2019_VIA_join$Percent_of_this_plot*100
MOI_1hr_2019_APOP_join$Percent_of_this_plot<- MOI_1hr_2019_APOP_join$Percent_of_this_plot*100
MOI_1hr_2019_PI_join$Percent_of_this_plot<-   MOI_1hr_2019_PI_join$Percent_of_this_plot*100
MOI_1hr_2019_JC1_join$Percent_of_this_plot <- MOI_1hr_2019_JC1_join$Percent_of_this_plot *100

#### MOI 1HR 2019 PLOTS AND STATISTICS ####

### VIABILITY ASSAY ####

MOI_1hr_2019_VIA_join

# Agranular and Granular cells Boxplot with significance bars, grouped by treatment color by gate(significant based on ANOVA)
MOI_1hr_2019_VIA_join_Percent_Agranular_Granular <- MOI_1hr_2019_VIA_join %>% filter(Gate =="P3" | Gate=="P4")
levels(factor(MOI_1hr_2019_VIA_join_Percent_Agranular_Granular$Plot_number)) # 8 correct

MOI_1hr_2019_VIA_join_Percent_Agranular_Granular_plot <- ggplot(data=MOI_1hr_2019_VIA_join_Percent_Agranular_Granular,
                                                              aes(y=Percent_of_this_plot, x=ID, fill=Gate)) + geom_col()  + 
  xlab("Individual") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent of Granular and Agranular MOI 1hr 2019") + 
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) + 
#  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  scale_fill_manual(name="Cell Type", labels=c("Agranular","Granular"), values=c("#6c81d9","#50b47b")) 

#save
ggsave(plot = MOI_1hr_2019_VIA_join_Percent_Agranular_Granular_plot, device = "tiff", filename = "MOI_1hr_2019_VIA_join_Percent_Agranular_Granular_plot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 10, width = 12)

# Plot of live and dead cells
MOI_1hr_2019_VIA_join_Percent_Agranular_Granular_LIVE_DEAD <- MOI_1hr_2019_VIA_join %>% filter(Gate =="V6-R" | Gate=="V7-R" |Gate=="V7-L"|Gate=="V6-L")
levels(factor(MOI_1hr_2019_VIA_join_Percent_Agranular_Granular_LIVE$Plot_number)) # 19, 20 correct

MOI_1hr_2019_VIA_join_Percent_Granular_LIVE_DEAD_plot <- MOI_1hr_2019_VIA_join_Percent_Agranular_Granular_LIVE_DEAD %>%
  # filter for granular, V7 from plot 20
  filter(Gate=="V7-R" |Gate=="V7-L") %>%
  ggplot(aes(y=Percent_of_this_plot, x=ID, fill=Gate)) + geom_col()  + 
  xlab("Individual") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Live and Dead Granular Hemocytes") + 
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) + 
  #  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  scale_fill_manual(name="Cell Type", labels=c("Dead","Live"), values=c("#5b2c90", "#88bf3b")) 

  #"#cc57b4", "#88bf3b", "#aa4dce", "#56b464", "#5b2c90"
  #"#c89832", "#5a6ee6", "#ca4e33", "#7e78d4", "#cd4272"

#save
ggsave(plot = MOI_1hr_2019_VIA_join_Percent_Granular_LIVE_DEAD_plot , device = "tiff", filename = "MOI_1hr_2019_VIA_join_Percent_Granular_LIVE_DEAD_plot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 10, width = 12)

MOI_1hr_2019_VIA_join_Percent_Agranular_LIVE_DEAD_plot <- MOI_1hr_2019_VIA_join_Percent_Agranular_Granular_LIVE_DEAD %>%
  # filter for agranular, V7 from plot 20
  filter(Gate=="V6-R" |Gate=="V6-L") %>%
  ggplot(aes(y=Percent_of_this_plot, x=ID, fill=Gate)) + geom_col()  + 
  xlab("Individual") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Live and Dead Agranular Hemocytes") + 
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) + 
  #  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  scale_fill_manual(name="Cell Type", labels=c("Dead","Live"), values=c("#5b2c90", "#88bf3b")) 

#save
ggsave(plot = MOI_1hr_2019_VIA_join_Percent_Agranular_LIVE_DEAD_plot, device = "tiff", filename = "MOI_1hr_2019_VIA_join_Percent_Agranular_LIVE_DEAD_plot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 10, width = 12)

# Viability assay results: the oysters have similar levels of granular and agranular hemocytes 


### Apoptosis Assay Statistics and Plotting ####

  # code in the sections below is for exploration of how to remove contaminating free apoptotic parasite from my Q1-UR quadrant

MOI_1hr_2019_APOP_join_Granular <- MOI_1hr_2019_APOP_join %>% filter(Plot_number == "4")

## Plot both the apoptotic quadrants together 

# combine the UL and UR quadrants to get combined apoptosis levels 
MOI_1hr_2019_APOP_join_Granular_combined_apoptotic <- MOI_1hr_2019_APOP_join_Granular %>% filter(Gate ==  "Q1-UR" | Gate ==  "Q1-UL") %>% 
  filter(Treat != "PERK") %>% 
  mutate(ID_full = paste(ID, Treat, Assay, sep = "_")) %>%
  ungroup() %>% dplyr::group_by(ID_full) %>%
  mutate(Percent_of_this_plot_combined = sum(Percent_of_this_plot)) %>% distinct(ID_full, .keep_all = TRUE)

# calculate the arcsine transformed percentages
MOI_1hr_2019_APOP_join_Granular_combined_apoptotic$Percent_of_this_plot_combined_arcsine <- transf.arcsin(MOI_1hr_2019_APOP_join_Granular_combined_apoptotic$Percent_of_this_plot_combined*0.01)

##Plot apoptosis granulocytes with BOTH parasite and non-parasite combined in format for multipanel figure with multiple comparisons run 
MOI_1hr_2019_APOP_join_Granular_combined_apoptotic_sd <-   MOI_1hr_2019_APOP_join_Granular_combined_apoptotic %>% ungroup() %>%
 group_by(Treat) %>% mutate(mean = mean(Percent_of_this_plot_combined), sd = sd(Percent_of_this_plot_combined))

MOI_1hr_2019_APOP_join_Granular_combined_apoptotic_sd$Treat <- factor(MOI_1hr_2019_APOP_join_Granular_combined_apoptotic_sd$Treat,
                                                                       levels = c("Beads","FSW","HK","P11","P51","P101","P251","PERK"))

MOI_1hr_2019_APOP_join_Granular_apoptotic_sd_multipanel <- 
  ggplot(data=MOI_1hr_2019_APOP_join_Granular_combined_apoptotic_sd,
         aes(y=Percent_of_this_plot_combined, x=Treat)) + 
  geom_bar(aes(fill=Treat), position="dodge", stat = "summary", fill = "#6d8dd7")  + 
  geom_point(aes(x= Treat, shape = ID), size = 3) +
  labs(x = NULL , y ="% Granular Apoptotic") + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold", angle = 90, hjust = 1),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) +
  #scale_shape_manual(values = c(15,16,17)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  scale_x_discrete(labels = c("Beads"="Beads",
                              "FSW"="FSW",
                              "HK" = "Heat Killed *P. mar.*",
                              "P101"="*P. mar.* 10:1",
                              "P11"="*P. mar.* 1:1",
                              "P51"="*P. mar.* 5:1",
                              "P251"="*P. mar.* 25:1"))
                              #"PERK"="*P. mar.* alone")) 

MOI_1hr_2019_APOP_join_Granular_apoptotic_sd_multipanel <- 
  MOI_1hr_2019_APOP_join_Granular_apoptotic_sd_multipanel  + 
  theme(axis.text.x=ggtext::element_markdown(),
        legend.text = ggtext::element_markdown()) 

# Perform anova with Tukey test and generate stats dataframe
MOI_1hr_2019_APOP_join_Granular_combined_apoptotic_sd_AOV <- aov(Percent_of_this_plot_combined_arcsine ~ Treat, MOI_1hr_2019_APOP_join_Granular_combined_apoptotic_sd)
summary(MOI_1hr_2019_APOP_join_Granular_combined_apoptotic_sd_AOV)
stat_test_tukey <- tukey_hsd(MOI_1hr_2019_APOP_join_Granular_combined_apoptotic_sd_AOV) %>%
  add_significance(p.col = "p.adj")

# take only the significant columns
stat_test_tukey <- stat_test_tukey %>% filter(p.adj <= 0.05)

MOI_1hr_2019_APOP_join_Granular_apoptotic_sd_multipanel_sig <- 
  MOI_1hr_2019_APOP_join_Granular_apoptotic_sd_multipanel + stat_pvalue_manual(
    stat_test_tukey, label = "{p.adj} {p.adj.signif}",  tip.length = 0.01, y.position = c(70, 73, 76, 79,82,85), size = 3) +
  # add overall anova values 
  #stat_compare_means(method= "anova") +
  labs(subtitle = "Tukey HSD, Arcsine Percent ~ Treat")

# export plot 
ggsave(plot = MOI_1hr_2019_APOP_join_Granular_apoptotic_sd_multipanel_sig, device = "tiff", filename = "MOI_1hr_2019_APOP_join_Granular_apoptotic_sd_multipanel_sig.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 8, width = 5)

## Plot data as counts 
# combine the UL and UR quadrants to get combined apoptosis levels 
MOI_1hr_2019_APOP_join_Granular_combined_apoptotic_counts <- MOI_1hr_2019_APOP_join_Granular %>% filter(Gate ==  "Q1-UR" | Gate ==  "Q1-UL") %>% 
  #filter(Treat != "PERK") %>% 
  mutate(ID_full = paste(ID, Treat, Assay, sep = "_")) %>%
  ungroup() %>% dplyr::group_by(ID_full) %>%
  mutate(counts_combined = sum(Counts)) %>% distinct(ID_full, .keep_all = TRUE)

#Plot apoptosis granulocytes with BOTH parasite and non-parasite combined in format for multipanel figure with multiple comparisons run 
MOI_1hr_2019_APOP_join_Granular_combined_apoptotic_counts_sd <-   MOI_1hr_2019_APOP_join_Granular_combined_apoptotic_counts %>% ungroup() %>%
  group_by(Treat) %>% mutate(mean = mean(counts_combined), sd = sd(counts_combined))

MOI_1hr_2019_APOP_join_Granular_combined_apoptotic_counts_sd$Treat <- factor(MOI_1hr_2019_APOP_join_Granular_combined_apoptotic_counts_sd$Treat,
                                                                      levels = c("Beads","FSW","HK","P11","P51","P101","P251","PERK"))

MOI_1hr_2019_APOP_join_Granular_apoptotic_counts_sd_multipanel <- 
  ggplot(data=MOI_1hr_2019_APOP_join_Granular_combined_apoptotic_counts_sd,
         aes(y=counts_combined, x=Treat)) + 
  geom_bar(aes(fill=Treat), position="dodge", 
          stat = "summary", 
          fill = "#6d8dd7")  + 
  geom_point(aes(x= Treat, shape = ID), size = 3) +
  #geom_boxplot(aes(fill=Treat),fill = "#6d8dd7") + 
  labs(x = NULL , y ="Granular Apoptotic Counts") + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold", angle = 90, hjust = 1),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) +
  scale_shape_manual(values = c(15,16,17,0,1,2,3)) +
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  #scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  scale_x_discrete(labels = c("Beads"="Beads",
                              "FSW"="FSW",
                              "HK" = "Heat Killed *P. mar.*",
                              "P101"="*P. mar.* 10:1",
                              "P11"="*P. mar.* 1:1",
                              "P51"="*P. mar.* 5:1",
                              "P251"="*P. mar.* 25:1",
                              "PERK"="*P. mar.* alone")) 

MOI_1hr_2019_APOP_join_Granular_apoptotic_counts_sd_multipanel <- 
  MOI_1hr_2019_APOP_join_Granular_apoptotic_counts_sd_multipanel  + 
  theme(axis.text.x=ggtext::element_markdown(),
        legend.text = ggtext::element_markdown()) 

# export plot 
ggsave(plot = MOI_1hr_2019_APOP_join_Granular_apoptotic_counts_sd_multipanel, device = "tiff", filename = "MOI_1hr_2019_APOP_join_Granular_apoptotic_counts_sd_multipanel.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 8, width = 5)

## Plot apoptotic parasite and compare between the different agranular and granular plots
MOI_1hr_2019_APOP_join_total_PERK_apop_UR_plot <- MOI_1hr_2019_APOP_join %>% filter(Treat =="PERK") %>% filter(Plot_number == "4" | Plot_number == "7") %>% 
  filter(Gate == "Q1-UR" | Gate == "Q2-UR") %>%  
  ggplot(data=.,
         aes(y=Percent_of_this_plot, x=Cell_type)) + 
  geom_bar(aes(fill=Cell_type), position="dodge", 
           stat = "summary")  +
  geom_point(aes(x = Cell_type), size = 3)+
  #geom_boxplot(aes(fill=Treat),fill = "#6d8dd7") + 
  labs(x = NULL , y ="Percent Apoptotic") + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold", angle = 90, hjust = 1),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) +
  scale_fill_manual(values=c("#6778d0","#b0923b"))
  
  ggsave(plot = MOI_1hr_2019_APOP_join_total_PERK_apop_UR_plot, device = "tiff", filename = "MOI_1hr_2019_APOP_join_total_PERK_apop_UR_plot.tiff",
         path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
         height = 5, width = 5)

## Plot total apoptotic parasite!
View(MOI_1hr_2019_APOP_join)

# The parasite counts were split up into "granular" and "Agranular" though these are really artificial separations. 
# Going to add together the counts for both UR and LR quadrants and recalculate the percent apoptosis
MOI_1hr_2019_APOP_join_total_PERK_apop <- MOI_1hr_2019_APOP_join %>% filter(Treat =="PERK") %>% filter(Plot_number == "4" | Plot_number == "7") %>% 
  filter(Gate == "Q1-UR" | Gate == "Q1-LR" | Gate == "Q2-UR" | Gate == "Q2-LR") %>%
  mutate(ID_full = paste(ID,Treat,Assay), group = case_when(
    grepl("UR", Gate) ~ "apoptotic",
    grepl("LR",Gate) ~ "live")) %>% 
  ungroup() %>% group_by(ID_full, group) %>% 
  mutate(Counts_group_sum = sum(Counts)) %>% distinct(ID_full, group, .keep_all = TRUE) %>%
  ungroup() %>%
  group_by(ID_full) %>%
  mutate(Counts_total = sum(Counts_group_sum), Percent_of_this_plot = Counts_group_sum / Counts_total * 100)

MOI_1hr_2019_APOP_join_total_PERK_apop$ID <- factor(MOI_1hr_2019_APOP_join_total_PERK_apop$ID, levels = c("P11","P5","P101","P251"))

MOI_1hr_2019_APOP_join_total_PERK_apop %>% filter(Gate == "Q1-UR") %>% group_by(Assay) %>% summarize(mean=mean(Percent_of_this_plot))

# Make plot of the parasite PERK alone treatment counts and percentage when not separated out by "granular" or "agranular"
MOI_1hr_2019_APOP_join_total_PERK_apop_counts <-  MOI_1hr_2019_APOP_join_total_PERK_apop %>% filter(group == "apoptotic") %>%
  ggplot(data=.,
         aes(y=Counts_group_sum, x=ID)) + 
  geom_bar(aes(fill=group), position="dodge", 
           stat = "summary")  +
  #geom_boxplot(aes(fill=Treat),fill = "#6d8dd7") + 
  labs(x = NULL , y ="Combined Granular Agranular Counts") + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold", angle = 90, hjust = 1),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) +
  scale_fill_manual(values=c("#6778d0","#b0923b"))

MOI_1hr_2019_APOP_join_total_PERK_apop_perc <- MOI_1hr_2019_APOP_join_total_PERK_apop %>% filter(group == "apoptotic") %>%
  ggplot(data=.,
         aes(y=Percent_of_this_plot, x=ID)) + 
  geom_bar(aes(fill=group), position="dodge", 
           stat = "summary")  + 
  #geom_boxplot(aes(fill=Treat),fill = "#6d8dd7") + 
  labs(x = NULL , y ="Combined Granular Agranular Percent Apoptotic") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold", angle = 90, hjust = 1),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,25)) +
  scale_fill_manual(values="#6778d0")

# export plots
PERK_alone_apop <- cowplot::plot_grid(MOI_1hr_2019_APOP_join_total_PERK_apop_counts, MOI_1hr_2019_APOP_join_total_PERK_apop_perc)
ggsave(plot = PERK_alone_apop, device = "tiff", filename = "MOI_1hr_2019_APOP_PERK_alone_apop.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 10)

## Plotting Apoptotic parasite in agranular plots with hemocytes added to see changes in response to exposure to hemocytes
MOI_1hr_2019_APOP_join_agraular_PERK <- MOI_1hr_2019_APOP_join %>% filter(Plot_number == "7") %>% filter(Gate ==  "Q2-UR" | Gate ==  "Q2-LR") %>% 
  filter( Treat == "P11" | Treat == "P51" | Treat == "P101" | Treat == "P251" | Treat =="FSW" ) %>% 
  mutate(ID_full= paste(ID, Treat,Assay)) %>% group_by(ID_full) %>%
  mutate(Combined_parasite = sum(Counts)) %>% ungroup() %>% 
  mutate(Percent_of_combined_parasite = Counts/Combined_parasite*100) %>%
  distinct(ID_full, Gate , .keep_all = TRUE)

MOI_1hr_2019_APOP_join_agraular_PERK $Treat <- factor(MOI_1hr_2019_APOP_join_agraular_PERK $Treat,
                                                                      levels = c("FSW","P11","P51",
                                                                                 "P101","P251"))

# Make plot of the parasite PERK alone treatment counts and percentage when not separated out by "granular" or "agranular"
MOI_1hr_2019_APOP_join_agraular_PERK_counts <- MOI_1hr_2019_APOP_join_agraular_PERK %>% filter(Gate == "Q2-UR") %>% 
  ggplot(data=.,
         aes(y=Counts, x=Treat)) + 
  geom_bar(aes(x = Treat), position="dodge", 
           stat = "summary", fill = "#6d8dd7")  +
  geom_point(aes(x = Treat, shape = ID), size = 3) +
  #geom_boxplot(aes(fill=Treat),fill = "#6d8dd7") + 
  labs(x = NULL , y ="Agranular Parasite Counts in Hemocyte Treatments") + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold", angle = 90, hjust = 1),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) 

MOI_1hr_2019_APOP_join_agraular_PERK_perc <- MOI_1hr_2019_APOP_join_agraular_PERK %>% filter(Gate == "Q2-UR") %>%
  ggplot(data=.,
         aes(y=Percent_of_combined_parasite, x=Treat)) + 
  geom_bar(aes(fill=Treat), position="dodge", 
           stat = "summary", fill = "#6d8dd7")  + 
  geom_point(aes(x = Treat, shape = ID), size = 3) +
  #geom_boxplot(aes(fill=Treat),fill = "#6d8dd7") + 
  labs(x = NULL , y ="Agranular Percent Apoptotic in Hemocyte Treatments") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold", angle = 90, hjust = 1),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,25)) 

# Plot agranular hemocyte percent apoptosis from PERK alone 
MOI_1hr_2019_APOP_join_Agrnular_PERK_alone <- MOI_1hr_2019_APOP_join %>% filter(Plot_number == "7") %>%  
  filter(Treat == "PERK" ) %>% filter(Gate ==  "Q2-UR" | Gate == "Q2-LR") %>%
  mutate(ID_full= paste(ID, Treat,Assay)) %>% group_by(ID_full) %>%
  mutate(Combined_parasite = sum(Counts)) %>% ungroup() %>% 
  mutate(Percent_of_combined_parasite = Counts/Combined_parasite*100) %>%
  distinct(ID_full, Gate , .keep_all = TRUE)

MOI_1hr_2019_APOP_join_Agrnular_PERK_alone$ID <- factor(MOI_1hr_2019_APOP_join_Agrnular_PERK_alone$ID , levels= c("P11","P5","P101","P251"))

MOI_1hr_2019_APOP_join_agraular_PERK_agranular_perc <-
  MOI_1hr_2019_APOP_join_Agrnular_PERK_alone %>% filter(Gate == "Q2-UR") %>%
  ggplot(data=.,
         aes(y=Percent_of_this_plot, x=ID)) + 
  geom_bar(aes(fill=Gate), position="dodge", 
           stat = "summary")  + 
  #geom_boxplot(aes(fill=Treat),fill = "#6d8dd7") + 
  labs(x = NULL , y ="Agranular Percent Apoptotic Perkinsus Control") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold", angle = 90, hjust = 1),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,25)) +
  scale_fill_manual(values="#6778d0")

MOI_1hr_2019_APOP_join_agraular_PERK_agranular_count <-
  MOI_1hr_2019_APOP_join_Agrnular_PERK_alone %>% filter(Gate == "Q2-UR") %>%
  ggplot(data=.,
         aes(y=Counts, x=ID)) + 
  geom_bar(aes(fill=Gate), position="dodge", 
           stat = "summary")  + 
  #geom_boxplot(aes(fill=Treat),fill = "#6d8dd7") + 
  labs(x = NULL , y ="Agranular Perkinsus Counts Control") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold", angle = 90, hjust = 1),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) + 
  #scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,25)) +
  scale_fill_manual(values="#6778d0")

# export plots
PERK_agranular_apop <- cowplot::plot_grid(MOI_1hr_2019_APOP_join_agraular_PERK_counts, MOI_1hr_2019_APOP_join_agraular_PERK_agranular_count, MOI_1hr_2019_APOP_join_agraular_PERK_perc, MOI_1hr_2019_APOP_join_agraular_PERK_agranular_perc,
                                          nrow = 1, ncol = 4)
ggsave(plot = PERK_agranular_apop, device = "tiff", filename = "MOI_1hr_2019_APOP_PERK_agranular_apop.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 15)


### Finally, compare Q2-UR in the total hemocytes rather than those in just the agranular or granular quadrants and then compare the percentage
MOI_1hr_2019_APOP_join_apop_granular_agranular_hemocyte <- MOI_1hr_2019_APOP_join %>% filter(Treat !="PERK") %>% filter(Treat !="HK") %>%
  filter(Plot_number == "4" | Plot_number == "7") %>% 
  filter(Gate == "Q1-UR" | Gate == "Q2-UR" ) %>%
  mutate(ID_full = paste(ID,Treat,Assay)) %>% 
  mutate(Counts_sum = sum(Counts)) %>% distinct(ID_full, .keep_all = TRUE) 

MOI_1hr_2019_APOP_join_apop_granular_agranular_hemocyte$Treat <- factor(MOI_1hr_2019_APOP_join_apop_granular_agranular_hemocyte$Treat, 
                                                                        levels = c("FSW","Beads","P11","P51","P101","P251"))

MOI_1hr_2019_APOP_join_apop_granular_agranular_hemocyte_count <-
  ggplot(data=MOI_1hr_2019_APOP_join_apop_granular_agranular_hemocyte,
         aes(y=Counts, x=Treat)) + 
  geom_bar(aes(fill=Treat), position="dodge", 
           stat = "summary", fill = "#6d8dd7")  + 
  geom_point(aes(x = Treat, shape = ID), size = 3) + 
  #geom_boxplot(aes(fill=Treat),fill = "#6d8dd7") + 
  labs(x = NULL , y ="Combined Granular Agranular Counts") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold", angle = 90, hjust = 1),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) + 
  scale_y_continuous(limits=c(0,15000), breaks = c(0,1000,2000,3000,4000,5000,10000,15000)) +
  scale_fill_manual(values="#6778d0")

# compare perkinsus apoptotic counts and hemocyte + perkinsus apoptotic counts 

perk_hemo_total_compare <- cowplot::plot_grid(MOI_1hr_2019_APOP_join_apop_granular_agranular_hemocyte_count, MOI_1hr_2019_APOP_join_total_PERK_apop_counts)
ggsave(plot = perk_hemo_total_compare, device = "tiff", filename = "MOI_1hr_2019_APOP_perk_hemo_total_compare.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 10)




### CORRECTION OF DATA FOR APOPTOSIS OF FREE PARASITE ###

### STRATEGY TO REMOVE CONTAMINATING FREE PARASITE APOPTOSIS FROM HEMOCYTE + LIVE PARASITE TREATMENTS 
#1. Add together the total number of cells in quadrants across Q1 and Q2 
#2. Get count of how many parasite cells were measured by adding together UR and LR 
#3. To get the approximate amount of apoptotic parasite in granular P4 (Q1-UR) - multiply the total apoptotic cells by the average proportion of Dermo
# cells that are usually in the granular cell plot (from the P1 plot P3 and P4 gate proportions for the perkinsus only control samples), 
# then multiply that by the average proportion of granular cell apoptosis in the perkinsus only control samples. 
#4. Multiply that number by the approximate amount of phagocytosis of granular sized cells from my 2020 assay. Subtract the amount of phagocytosed apoptotic granular from 
    # the estimated perkinsus apoptosis granular.
#5. Subtract this final number of estimated non phagocytosed apoptotic granular perkinsus from the original Q1-UR quadrant for just the parasite samples 

MOI_1hr_2019_APOP_join_total <- MOI_1hr_2019_APOP_join %>% filter(Gate == "Q1-UR" | Gate == "Q1-UL" | Gate == "Q1-LL" | Gate == "Q1-LR" | 
                                                                    Gate == "Q2-UR" | Gate == "Q2-UL" | Gate == "Q2-LL" | Gate == "Q2-LR") %>%
  separate(Gate, c("plot","quadrant"), sep = "-") %>%
  mutate(ID_full = paste(ID, Treat, Assay)) %>% ungroup() %>% 
  group_by(ID_full, quadrant) %>% 
  mutate(Counts_total = sum(Counts)) %>% distinct(ID_full, quadrant, .keep_all = TRUE) %>% select(ID,Treat,Assay, ID_full, Channel,quadrant,Cell_type,Counts_total)

MOI_1hr_2019_APOP_join_total_Perkinsus <- MOI_1hr_2019_APOP_join_total %>% filter(quadrant == "UR" | quadrant == "LR") %>% ungroup() %>%
  group_by(ID_full) %>% mutate(Total_parasite_counts = sum(Counts_total)) %>% distinct(ID_full, .keep_all = TRUE) %>% select(ID, Treat, Assay, ID_full, Channel, Total_parasite_counts)

# find average percent of perkinsus cells in perkinsus controls that are in the granular area of the plot 
MOI_1hr_2019_APOP_join_granular_percent <-  MOI_1hr_2019_APOP_join %>% filter(Gate == "P4") %>% filter(Treat == "PERK") %>% ungroup() %>%
  summarize(mean = mean(Percent_of_this_plot))
    # mean = 10.3 percent 

# find average percent of perkinsus cell apoptosis in the Q1-UR granular apoptotic quadrant 
MOI_1hr_2019_APOP_join_granular_percent_perk_apop <-  MOI_1hr_2019_APOP_join %>% filter(Gate == "Q1-UR") %>% filter(Treat == "PERK") %>% ungroup() %>%
  summarize(mean = mean(Percent_of_this_plot))
  # mean = 23.1 

# calculate approximate apoptotic granular perkinsus cells 
MOI_1hr_2019_APOP_join_total_Perkinsus_apoptotic <- MOI_1hr_2019_APOP_join_total_Perkinsus %>% filter(Treat == "P101" | Treat ==  "P11" | Treat ==  "P251" | Treat == "P51") %>%
  mutate(Parasite_granular = Total_parasite_counts * 0.103, Parasite_granular_apoptosis = Parasite_granular * 0.231, Parasite_gran_apop_phago = Parasite_granular_apoptosis *0.13219,
         Parasite_granular_apop_minus_phago = Parasite_granular_apoptosis - Parasite_gran_apop_phago) %>%
  # ronud to whole numbers
  mutate(across(6:9, ceiling))

# Subtract granular apoptotic parasite from Q1-UR for hemocyte treatments by joining to original plot
MOI_1hr_2019_APOP_join_Q1UR_hemo <- MOI_1hr_2019_APOP_join %>% filter(Gate == "Q1-UR") %>% filter(Treat == "P101" | Treat ==  "P11" | Treat ==  "P251" | Treat == "P51") %>%
  mutate(ID_full = paste(ID,Treat,Assay), sep = "_")

MOI_1hr_2019_APOP_join_total_Perkinsus_apoptotic_Q1UR_join <- left_join(MOI_1hr_2019_APOP_join_total_Perkinsus_apoptotic, 
                                                                        MOI_1hr_2019_APOP_join_Q1UR_hemo[,c("ID_full","Counts")], by = "ID_full") %>%
  mutate(Counts_minus_free_PM_apop = Counts - Parasite_granular_apop_minus_phago) %>% select(ID,Treat,Assay, ID_full, Counts_minus_free_PM_apop) %>%
  #create column that says gate for joining
  mutate(Gate = "Q1-UR") %>%
  # remove Counts_minus_free_PM_apop to counts for joining back with rest of hemocyte treatment data
  rename(Counts = Counts_minus_free_PM_apop)

# Separate dataframes so that I can bring them back together..
MOI_1hr_2019_APOP_join_Q1_other_hemo <- MOI_1hr_2019_APOP_join %>% filter(Gate == "Q1-UL" | Gate == "Q1-LL" | Gate == "Q1-LR") %>% 
  filter(Treat ==  "P11" ) %>%
  mutate(ID_full = paste(ID,Treat,Assay), sep = "_") %>%
  # remove percent of this plot and percent of this plot arcsine since they will be recalculated 
  select(ID, Treat, Assay, ID_full, Gate, Counts)

# Join full hemo treatment quadrants and recalculate percent and percent of this plot arcsine, and focus on Q1 granular plot
MOI_1hr_2019_APOP_join_hemo_recalc <- rbind(MOI_1hr_2019_APOP_join_Q1_other_hemo, MOI_1hr_2019_APOP_join_total_Perkinsus_apoptotic_Q1UR_join) %>%
  ungroup() %>%  filter(Treat == "P11")  %>% group_by(ID_full) %>%
  mutate(Total_counts = sum(Counts)) %>% ungroup() %>%
  mutate(Percent_of_this_plot = Counts/Total_counts * 100) %>% select(ID, Treat,Assay,ID_full, Gate,Counts, Percent_of_this_plot)

MOI_1hr_2019_APOP_join_hemo_recalc$Percent_of_this_plot_arcsine <- transf.arcsin(MOI_1hr_2019_APOP_join_hemo_recalc$Percent_of_this_plot*0.01)

# Subset other treatments, with only granular plots that I want to compare and subset columns
MOI_1hr_2019_APOP_join_granular_other <- MOI_1hr_2019_APOP_join %>% filter(Treat == "Beads" |Treat == "FSW"  |Treat == "HK") %>% 
  filter(Gate == "Q1-UL" | Gate == "Q1-LL" | Gate == "Q1-LR" | Gate == "Q1-UR")  %>% mutate(ID_full = paste(ID, Treat, Assay, sep = "_")) %>%
  select(ID, Treat,Assay,ID_full, Gate,Counts, Percent_of_this_plot, Percent_of_this_plot_arcsine)

# Join all together into one data frame of just the adjusted granular, the treatments, and no Perkinsus only
MOI_1hr_2019_APOP_join_granular_recalc_all_treat <- rbind(as.data.frame(MOI_1hr_2019_APOP_join_hemo_recalc),as.data.frame(MOI_1hr_2019_APOP_join_granular_other))

### Plotting percent combined apototic from all treatmetns (Q1-UR + Q1-UL)
# combine the UL and UR quadrants to get combined apoptosis levels 
MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic <- MOI_1hr_2019_APOP_join_granular_recalc_all_treat %>% filter(Gate ==  "Q1-UR" | Gate ==  "Q1-UL") %>% 
  ungroup() %>% dplyr::group_by(ID_full) %>%
  mutate(Percent_of_this_plot_combined = sum(Percent_of_this_plot)) %>% distinct(ID_full, .keep_all = TRUE)

# calculate the arcsine transformed percentages
MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic$Percent_of_this_plot_combined_arcsine <- transf.arcsin(MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic$Percent_of_this_plot_combined*0.01)

##Plot apoptosis granulocytes with BOTH parasite and non-parasite combined in format for multipanel figure with multiple comparisons run 
MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd <-   MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic %>% ungroup() %>%
  group_by(Treat) %>% mutate(mean = mean(Percent_of_this_plot_combined), sd = sd(Percent_of_this_plot_combined))

MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd$Treat <- factor(MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd$Treat,
                                                                      levels = c("Beads","FSW","HK","P11"))

MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd_multipanel <- 
  ggplot(data=MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd,
         aes(y=Percent_of_this_plot_combined, x=Treat)) + 
  geom_bar(aes(fill=Treat), position="dodge", stat = "summary", fill = "#6d8dd7")  + 
  geom_point(aes(x= Treat, shape = ID), size = 3) +
  labs(x = NULL , y ="% Granular Apoptotic") + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold", angle = 90, hjust = 1),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) +
  #scale_shape_manual(values = c(15,16,17)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  scale_x_discrete(labels = c("Beads"="Beads",
                              "FSW"="FSW",
                              "HK" = "Heat Killed *P. mar.*",
                              "P11"="*P. mar.* 1:1"))

MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd_multipanel <- 
  MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd_multipanel  + 
  theme(axis.text.x=ggtext::element_markdown(),
        legend.text = ggtext::element_markdown()) 

# Perform anova with Tukey test and generate stats dataframe
MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd_AOV <- aov(Percent_of_this_plot_combined_arcsine ~ Treat, MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd)
summary(MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd_AOV)
stat_test_tukey <- tukey_hsd(MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd_AOV) %>%
  add_significance(p.col = "p.adj")

# take only the significant columns
stat_test_tukey <- stat_test_tukey %>% filter(p.adj <= 0.05)

MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd_multipanel_sig <- 
  MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd_multipanel + stat_pvalue_manual(
    stat_test_tukey, label = "{p.adj} {p.adj.signif}",  tip.length = 0.01, y.position = c(75,80,85), size = 3) +
  # add overall anova values 
  #stat_compare_means(method= "anova") +
  labs(subtitle = "Tukey HSD, Arcsine Percent ~ Treat")

# export plot 
ggsave(plot = MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd_multipanel_sig, device = "tiff", filename = "MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd_multipanel_sig.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 8, width = 5)


## T-test

MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd_ttest <- MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd %>% filter(Treat == "FSW" | Treat == "P11")
t.test(Percent_of_this_plot_combined_arcsine ~Treat,  MOI_1hr_2019_APOP_join_granular_recalc_all_treat_combined_apoptotic_sd_ttest) # pvalue = p-value = 0.2476

### PI Assay Statistics and Plotting ####

MOI_1hr_2019_PI_join
# plot 4 and plot 7 are the ones I again care about , but now the gates are Q2 and Q4

MOI_1hr_2019_PI_join_Agranular_Granular <- MOI_1hr_2019_PI_join %>% filter(Gate =="Q4-UL" | Gate =="Q4-UR" | Gate =="Q4-LL" | Gate =="Q4-LR"|
                                                                                 Gate =="Q2-UL" | Gate =="Q2-UR" | Gate =="Q2-LL" | Gate =="Q2-LR")

MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic <- MOI_1hr_2019_PI_join_Agranular_Granular %>% filter(Gate ==  "Q4-UR" | Gate =="Q2-UR" | Gate ==  "Q4-UL" | Gate =="Q2-UL")
levels(factor(MOI_1hr_2019_PI_join_Agranular_Granular$Plot_number)) # 4, 7
levels(factor(MOI_1hr_2019_PI_join$Treat)) # "Beads" "FSW"   "HK"    "P101"  "P11"   "P251"  "P51"   "PERK" 

# Analysis for all Granular and Agranular Combined from the quad plots
MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic$Treat <- factor(MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic$Treat, 
                                                                    levels = c( "FSW","Beads",   "HK",    "P11" , "P51",  "P101",   "P251",   "PERK" ),
                                                                    labels = c( "FSW and\n Hemocytes","Beads to\n Hemocytes\n 1:1",   "Heat-Killed\n P. mar to\n Hemocytes\n 1:1",    "P. mar to\n Hemocytes\n 1:1" , 
                                                                                "P. mar to\n Hemocytes\n 5:1",  
                                                                                "P. mar to\n Hemocytes\n 10:1",   "P. mar to\n Hemocytes\n 25:1",   "Perkinsus\n Only\n Control" ))
# Make plot 
MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_all_samplesplot <- ggplot(data=MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic,
                                                                              aes(y=Percent_of_this_plot, x=Treat, color=Gate)) + geom_point(position=position_dodge(width=0.75)) + 
  geom_boxplot() +
  xlab("Treatment") +
  ylab("Percent of Each Cell Type") + 
  ggtitle("Percent of Each Cell Type Positive\n for Cell Death and Parasite, Beads, or HK P. marinus") + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100), breaks = c(0,5,10,20,30,40,50,60,70,80,90,100)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  scale_color_manual(name="Cell Type", labels=c("Dead \nGranular", "Dead Granular,\nPerkinsus, or Beads", "Dead\nAgranular", "Dead Agranular,\nPerkinsus or Beads"), 
                     values = c("#cc57b4", "#7e78d4", "#56b464", "#5b2c90")) 

#save
ggsave(plot = MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_all_samplesplot, device = "tiff", filename = "MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_all_samplesplot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 10)

### Do cell types differ in cell death levels? - 
# Need to run the FSW control test as a separate t.test because I'm comparing the UL quadrants instead of UR
MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_cell_type_AOV <- MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic %>%
  filter(Treat == "Beads to\n Hemocytes\n 1:1" | Treat ==  "Heat-Killed\n P. mar to\n Hemocytes\n 1:1"| Treat ==   "P. mar to\n Hemocytes\n 1:1" | 
           Treat == "P. mar to\n Hemocytes\n 5:1" | Treat == "P. mar to\n Hemocytes\n 10:1"| Treat == "P. mar to\n Hemocytes\n 25:1"|
           Treat ==  "Perkinsus\n Only\n Control" ) %>%
  filter(Gate == "Q4-UR" | Gate == "Q2-UR") %>%
  group_by(Treat) %>%
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Gate, data = .)))  %>%
  ungroup

MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_cell_type_AOV %>% filter(p.value <= 0.05)
# granular have more apoptosis in all categories

# FSW alone t.test
MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_FSW <- MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic %>%
  filter(Treat == "FSW and\n Hemocytes" ) %>% filter(Gate == "Q4-UL" | Gate == "Q2-UL")
t.test(Percent_of_this_plot_arcsine ~ Gate , data =MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_FSW) 
# p-value = 0.056

### Does Perkinsus cause significantly greater granular cell death than beads? ### - not correcting for amount of apoptosis caused by Perkinsus alone
# comparing only the 1:1
MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_beads_vs_hemo <- MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic %>% 
  filter(Treat == "Beads to\n Hemocytes\n 1:1"  | Treat == "P. mar to\n Hemocytes\n 1:1") %>% 
  # filter out just the granulocytes
  filter(Gate == "Q2-UR")

MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_beads_vs_hemo_ttest <- t.test(Percent_of_this_plot_arcsine ~ Treat, data = MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_beads_vs_hemo)
# p-value = 0.00703

# test all concentrations with AOV
MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_beads_vs_hemo_aov <- MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic %>% 
  filter(Treat == "Beads to\n Hemocytes\n 1:1"  | Treat == "P. mar to\n Hemocytes\n 1:1" | Treat == "P. mar to\n Hemocytes\n 5:1" |
           Treat == "P. mar to\n Hemocytes\n 10:1" | Treat == "P. mar to\n Hemocytes\n 25:1") %>% 
  # filter out just the granulocytes
  filter(Gate == "Q2-UR")

MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_beads_vs_hemo_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data = MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_beads_vs_hemo_aov)
TukeyHSD(MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_beads_vs_hemo_aov)
## All concentrations have significantly different granular cell death as compared to the beads
#P. mar to Hemocytes 1:1-Beads to Hemocytes 1:1    0.0054883
#P. mar to Hemocytes 5:1-Beads to Hemocytes 1:1    0.0005208
#P. mar to Hemocytes 10:1-Beads to Hemocytes 1:1   0.0007265
#P. mar to Hemocytes 25:1-Beads to Hemocytes 1:1   0.0009005

### Does Perkinsus cause significantly greater cell death than FSW? - not correcting for amount of apoptosis caused by Perkinsus alone ###
MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_FSW_granular <- MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic %>% 
  filter(Treat == "FSW and\n Hemocytes" & Gate == "Q2-UL")
MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_hemo_granular <- MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic %>% 
  filter(Treat == "P. mar to\n Hemocytes\n 1:1" | Treat == "P. mar to\n Hemocytes\n 5:1" |
           Treat == "P. mar to\n Hemocytes\n 10:1" | Treat == "P. mar to\n Hemocytes\n 25:1") %>% filter(Gate == "Q2-UR")        

MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_FSW_hemo_granular <- rbind(MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_FSW_granular, MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_hemo_granular)

# test with AOV
MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_FSW_vs_hemo_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data = MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_FSW_hemo_granular)
TukeyHSD(MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_FSW_vs_hemo_aov )  
# P. marinus not significantly different from 

# Does granular apoptosis differ between the different MOIs?
MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_MOI <- MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic %>% 
  filter(Treat == "P. mar to\n Hemocytes\n 1:1"  | Treat == "P. mar to\n Hemocytes\n 5:1" | Treat == "P. mar to\n Hemocytes\n 10:1" |
           Treat == "P. mar to\n Hemocytes\n 25:1") %>% 
  # filter out just the granulocytes
  filter(Gate == "Q2-UR")

MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_MOI_aov <- aov(Percent_of_this_plot_arcsine ~ Treat,
                                                                   data = MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_MOI)
TukeyHSD(MOI_1hr_2019_PI_join_Agranular_Granular_apoptotic_MOI_aov)
# no MOIs are different from one another 
#diff        lwr       upr     p adj
#P. mar to Hemocytes 5:1-P. mar to Hemocytes\n 1:1    0.4234763
#P. mar to Hemocytes 10:1-P. mar to Hemocytes\n 1:1   0.5496332
#P. mar to Hemocytes 25:1-P. mar to Hemocytes\n 1:1   0.6353521
#P. mar to Hemocytes 10:1-P. mar to Hemocytes\n 5:1   0.9944681
#P. mar to Hemocytes 25:1-P. mar to Hemocytes\n 5:1   0.9771173
#P. mar to Hemocytes 25:1-P. mar to Hemocytes\n 10:1  0.9985764

### JC-1 Assay Statistics and Plotting ####
MOI_1hr_2019_JC1_join
# plot 14, 16, 18 and plot 16 are the ones I again care about
# Q11, Q17 = agranulocytes, Q22, Q28 are the granulocytes

MOI_1hr_2019_JC1_join_Agranular_Granular <- MOI_1hr_2019_JC1_join %>% 
  filter(Gate =="Q11-UL" | Gate =="Q11-UR" | Gate =="Q11-LL" | Gate =="Q11-LR"|
         Gate =="Q28-UL" | Gate =="Q28-UR" | Gate =="Q28-LL" | Gate =="Q28-LR" |
           Gate =="Q17-UL" | Gate =="Q17-UR" | Gate =="Q17-LL" | Gate =="Q17-LR"|
           Gate =="Q22-UL" | Gate =="Q22-UR" | Gate =="Q22-LL" | Gate =="Q22-LR")

# Analysis for all Granular and Agranular Combined from the quad plots
MOI_1hr_2019_JC1_join_Agranular_Granular$Treat <- factor(MOI_1hr_2019_JC1_join_Agranular_Granular$Treat, 
                                                                  levels = c( "FSW","Beads",   "HK",    "P11" , "P51",  "P101",   "P251",   "PERK", "hemocytes_CCP","PERK_CCP" ),
                                                                  labels = c( "FSW and\n Hemocytes","Beads to\n Hemocytes\n 1:1",   "Heat-Killed\n P. mar to\n Hemocytes\n 1:1",    "P. mar to\n Hemocytes\n 1:1" , 
                                                                              "P. mar to\n Hemocytes\n 5:1",  
                                                                              "P. mar to\n Hemocytes\n 10:1",   "P. mar to\n Hemocytes\n 25:1",   "Perkinsus\n Only\n Control", "Hemocytes\nInhibitor",
                                                                              "Perkinsus\nInhibitor"))
## Plot the percent of mitochondrial permeabilized hemocytes vs. all hemocytes of each type

# First need to add together the quadrants for the normal cells in the UR and UL categories for each cell type 
MOI_1hr_2019_JC1_join_Agranular_normal <- MOI_1hr_2019_JC1_join_Agranular_Granular %>%
  # subset for the quadrants in plot 14 and plot 17 - LR are the mitochondria permeabilized cells while UL and UR are the red normal 
  # need to add together the red normal cells 
  filter( Gate =="Q17-UL" | Gate =="Q17-UR") %>% 
  group_by(Treat, ID) %>%
  mutate(normal_percent_combined = sum(Percent_of_this_plot)) %>%
# also add together the counts 
mutate(normal_counts_combined = sum(Counts)) %>% 
# change the gate designation
  mutate(Gate = "Q17-UL_Q17-UR")

MOI_1hr_2019_JC1_join_Granular_normal <- MOI_1hr_2019_JC1_join_Agranular_Granular %>%
  # subset for the quadrants in plot 14 and plot 17 - LR are the mitochondria permeabilized cells while UL and UR are the red normal 
  # need to add together the red normal cells 
  filter( Gate =="Q22-UL" | Gate =="Q22-UR") %>% 
  group_by(Treat, ID) %>%
  mutate(normal_percent_combined = sum(Percent_of_this_plot))  %>%
# also add together the counts 
  mutate(normal_counts_combined = sum(Counts)) %>%
  # change the gate designation
  mutate(Gate = "Q22-UL_Q22-UR")

#rbind two cell types after addition
MOI_1hr_2019_JC1_join_Agranular_Granular_normal <- rbind(MOI_1hr_2019_JC1_join_Agranular_normal, MOI_1hr_2019_JC1_join_Granular_normal) %>%
  # remove the old percent of this plot column and rename
  dplyr::select(ID, Treat, Assay, Plot_number, Channel, Gate, Cell_type, normal_percent_combined, normal_counts_combined) %>%
  rename(Percent_of_this_plot = normal_percent_combined) %>% rename(Counts = normal_counts_combined)

# extract just the mitochondrial permeabilization positive
MOI_1hr_2019_JC1_join_Agranular_Granular_mito_perm <- MOI_1hr_2019_JC1_join_Agranular_Granular %>%
  # subset for the quadrants in plot 14 and plot 17 - LR are the mitochondria permeabilized cells while UL and UR are the red normal 
  # need to add together the red normal cells 
  filter( Gate =="Q17-LR" | Gate =="Q22-LR")

# Join mito perm quadrants back with the normal cell percentages
MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm <- rbind(MOI_1hr_2019_JC1_join_Agranular_Granular_normal, MOI_1hr_2019_JC1_join_Agranular_Granular_mito_perm)

# transform the newly combined percentages
MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm$Percent_of_this_plot_arcsine <- transf.arcsin(MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm$Percent_of_this_plot*0.01)

## Plot the normal combined quadrant and the mitochondrial permeabilized quadrants
MOI_1hr_2019_JC1_join_Agranular_Granular_all_samplesplot <- MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm %>%
  ggplot(aes(y=Percent_of_this_plot, x=Treat, color=Gate)) + geom_point(position=position_dodge(width=0.75)) + 
  geom_boxplot() +
  xlab("Treatment") +
  ylab("Percent of Each Cell Type") + 
  ggtitle("Percent of Each Cell Type Positive\n for Mitochondrial Permeabilization") + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100), breaks = c(0,5,10,20,30,40,50,60,70,80,90,100)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  scale_color_manual(name="Cell Type", labels=c("Apoptotic\nAgranular", "Normal\nAgranular", "Apoptotic\nGranular", "Normal\nGranular"), 
                     values = c("#5b2c90", "#56b464", "#7e78d4", "#cc57b4")) 

#save
ggsave(plot = MOI_1hr_2019_JC1_join_Agranular_Granular_all_samplesplot, device = "tiff", filename = "MOI_1hr_2019_JC1_join_Agranular_Granular_all_samplesplot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 12)


### ANAlYSIS OF ALL MITOCHONDRIAL PERMEABILIZED VS. NORMAL CELLS 
### Do cell types differ in mitochondrial permeabilization levels? - 
# Need to run the FSW control test as a separate t.test because I'm comparing the UL quadrants instead of UR
MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm_cell_type_AOV <- MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm %>%
  filter(Treat == "FSW and\n Hemocytes" | Treat == "Beads to\n Hemocytes\n 1:1" | Treat ==  "Heat-Killed\n P. mar to\n Hemocytes\n 1:1"| Treat ==   "P. mar to\n Hemocytes\n 1:1" | 
           Treat == "P. mar to\n Hemocytes\n 5:1" | Treat == "P. mar to\n Hemocytes\n 10:1"| Treat == "P. mar to\n Hemocytes\n 25:1"|
           Treat ==  "Perkinsus\n Only\n Control" ) %>%
  filter( Gate =="Q17-LR" | Gate =="Q22-LR") %>%
  group_by(Treat) %>%
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Gate, data = .)))  %>%
  ungroup

MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm_cell_type_AOV %>% filter(p.value <= 0.05)
#Treat                                       term     df  sumsq meansq statistic  p.value
#<fct>                                       <chr> <dbl>  <dbl>  <dbl>     <dbl>    <dbl>
#  1 "FSW and\n Hemocytes"                       Gate      1 0.0871 0.0871     101.  0.000552
#2 "Beads to\n Hemocytes\n 1:1"                Gate      1 0.0724 0.0724      63.6 0.00134 
#3 "Heat-Killed\n P. mar to\n Hemocytes\n 1:1" Gate      1 0.0353 0.0353      70.9 0.00109 


### Does Perkinsus cause significantly different GRANULAR mitochondrial permeabilization than beads? ### - not correcting for amount of apoptosis caused by Perkinsus alone
# test all concentrations with AOV
MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm_beads_vs_hemo <- MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm %>% 
  filter(Treat == "Beads to\n Hemocytes\n 1:1"  | Treat == "P. mar to\n Hemocytes\n 1:1" | Treat == "P. mar to\n Hemocytes\n 5:1" |
           Treat == "P. mar to\n Hemocytes\n 10:1" | Treat == "P. mar to\n Hemocytes\n 25:1") %>% 
  # filter out just the granulocytes
  filter(Gate == "Q22-LR")

MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm_beads_vs_hemo_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data = MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm_beads_vs_hemo)
TukeyHSD(MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm_beads_vs_hemo_aov)
# no significant difference

## Compare different between FSW and each concentration
MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm_FSW_vs_hemo <- MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm %>% 
  filter(Treat == "FSW and\n Hemocytes" | Treat == "P. mar to\n Hemocytes\n 1:1" | Treat == "P. mar to\n Hemocytes\n 5:1" |
           Treat == "P. mar to\n Hemocytes\n 10:1" | Treat == "P. mar to\n Hemocytes\n 25:1") %>% 
  # filter out just the granulocytes
  filter(Gate == "Q22-LR")

MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm_FSW_vs_hemo_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data = MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm_FSW_vs_hemo)
TukeyHSD(MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm_FSW_vs_hemo_aov)

# Does granular mitochondrial permeabilization differ between the different MOIs?
MOI_1hr_2019_JC1_join_Agranular_Granular_MOI <-  MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm %>% 
  filter(Treat == "P. mar to\n Hemocytes\n 1:1"  | Treat == "P. mar to\n Hemocytes\n 5:1" | Treat == "P. mar to\n Hemocytes\n 10:1" |
           Treat == "P. mar to\n Hemocytes\n 25:1") %>% 
  # filter out just the granulocytes
  filter(Gate == "Q22-LR")

MOI_1hr_2019_JC1_join_Agranular_Granular_MOI_aov <- aov(Percent_of_this_plot_arcsine ~ Treat,
                                                                  data = MOI_1hr_2019_JC1_join_Agranular_Granular_MOI )
TukeyHSD(MOI_1hr_2019_JC1_join_Agranular_Granular_MOI_aov)

### Does Perkinsus cause significantly different AGRANULAR mitochondrial permeabilization than beads? ### - not correcting for amount of apoptosis caused by Perkinsus alone
# test all concentrations with AOV
MOI_1hr_2019_JC1_join_Agranular_normal_mito_perm_beads_vs_hemo <- MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm %>% 
  filter(Treat == "Beads to\n Hemocytes\n 1:1"  | Treat == "P. mar to\n Hemocytes\n 1:1" | Treat == "P. mar to\n Hemocytes\n 5:1" |
           Treat == "P. mar to\n Hemocytes\n 10:1" | Treat == "P. mar to\n Hemocytes\n 25:1") %>% 
  # filter out just the granulocytes
  filter(Gate == "Q17-LR")

MOI_1hr_2019_JC1_join_Agranular_normal_mito_perm_beads_vs_hemo_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data = MOI_1hr_2019_JC1_join_Agranular_normal_mito_perm_beads_vs_hemo)
TukeyHSD(MOI_1hr_2019_JC1_join_Agranular_normal_mito_perm_beads_vs_hemo_aov)
# significant differences with all perkinsus concentrations 
#$Treat
#diff         lwr         upr     p adj
#P. mar to\n Hemocytes\n 1:1-Beads to\n Hemocytes\n 1:1    -0.13156472 -0.22046150 -0.04266794 0.0045700
#P. mar to\n Hemocytes\n 5:1-Beads to\n Hemocytes\n 1:1    -0.19664564 -0.28554242 -0.10774886 0.0002002
#P. mar to\n Hemocytes\n 10:1-Beads to\n Hemocytes\n 1:1   -0.16422171 -0.25311849 -0.07532493 0.0008719
#P. mar to\n Hemocytes\n 25:1-Beads to\n Hemocytes\n 1:1   -0.14801267 -0.23690944 -0.05911589 0.0019403

## Compare different between FSW and each concentration
MOI_1hr_2019_JC1_join_Agranular_normal_mito_perm_FSW_vs_hemo <- MOI_1hr_2019_JC1_join_Agranular_Granular_normal_mito_perm %>% 
  filter(Treat == "FSW and\n Hemocytes" | Treat == "P. mar to\n Hemocytes\n 1:1" | Treat == "P. mar to\n Hemocytes\n 5:1" |
           Treat == "P. mar to\n Hemocytes\n 10:1" | Treat == "P. mar to\n Hemocytes\n 25:1") %>% 
  # filter out just the granulocytes
  filter(Gate == "Q17-LR")

MOI_1hr_2019_JC1_join_Agranular_normal_mito_perm_FSW_vs_hemo_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data = MOI_1hr_2019_JC1_join_Agranular_normal_mito_perm_FSW_vs_hemo)
TukeyHSD(MOI_1hr_2019_JC1_join_Agranular_normal_mito_perm_FSW_vs_hemo_aov)
#$Treat
#diff         lwr         upr     p adj
#P. mar to\n Hemocytes\n 1:1-FSW and\n Hemocytes           -0.12821949 -0.21784254 -0.03859644 0.0057848
#P. mar to\n Hemocytes\n 5:1-FSW and\n Hemocytes           -0.19330041 -0.28292346 -0.10367736 0.0002476
#P. mar to\n Hemocytes\n 10:1-FSW and\n Hemocytes          -0.16087648 -0.25049953 -0.07125343 0.0010918
#P. mar to\n Hemocytes\n 25:1-FSW and\n Hemocytes          -0.14466744 -0.23429049 -0.05504438 0.0024441

# no difference between MOIs however



### COMPARISON OF PERCENTAGE OF PARASITE POSITIVE VS. PARASITE NEGATIVE CELLS 
## Plot the percent of parasite positive vs. parasite negative for the mitochondrial permeabilized 
# UR are for the parasite and apoptosis positive, while UL are for FSW hemocytes
MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic <- MOI_1hr_2019_JC1_join_Agranular_Granular %>% filter(Gate ==  "Q11-UR" | Gate =="Q28-UR" | Gate ==  "Q11-UL" | Gate =="Q28-UL")
levels(factor(MOI_1hr_2019_JC1_join_Agranular_Granular$Plot_number)) # 16, 18
levels(factor(MOI_1hr_2019_JC1_join$Treat)) # "Beads" "FSW" "hemocytes_CCP" "HK" "P101" "P11" "P251" "P51" "PERK" "PERK_CCP"

# Analysis for all Granular and Agranular Combined from the quad plots
MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic$Treat <- factor(MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic$Treat, 
                                                                   levels = c( "FSW","Beads",   "HK",    "P11" , "P51",  "P101",   "P251",   "PERK", "hemocytes_CCP","PERK_CCP" ),
                                                                   labels = c( "FSW and\n Hemocytes","Beads to\n Hemocytes\n 1:1",   "Heat-Killed\n P. mar to\n Hemocytes\n 1:1",    "P. mar to\n Hemocytes\n 1:1" , 
                                                                               "P. mar to\n Hemocytes\n 5:1",  
                                                                               "P. mar to\n Hemocytes\n 10:1",   "P. mar to\n Hemocytes\n 25:1",   "Perkinsus\n Only\n Control", "Hemocytes\nInhibitor",
                                                                               "Perkinsus\nInhibitor"))

# Make plot 
MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_all_samplesplot <- ggplot(data=MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic,
                                                                   aes(y=Percent_of_this_plot, x=Treat, color=Gate)) + geom_point(position=position_dodge(width=0.75)) + 
  geom_boxplot() +
  xlab("Treatment") +
  ylab("Percent of Each Cell Type") + 
  ggtitle("Percent of Each Cell Type Positive\n for Mitochondrial Permeabilization\nand Parasite, Beads, or HK P. marinus") + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100), breaks = c(0,5,10,20,30,40,50,60,70,80,90,100)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  scale_color_manual(name="Cell Type", labels=c("Apoptotic\nAgranular", "Apoptotic Agranular,\nPerkinsus, or Beads", "Apoptotic\nGranular", "Apoptotic Granular,\nPerkinsus or Beads"), 
                     values = c("#56b464", "#5b2c90","#cc57b4", "#7e78d4")) 

#save
ggsave(plot = MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_all_samplesplot, device = "tiff", filename = "MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_all_samplesplot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 12)


### Do cell types differ in mitochondrial permeabilization levels? - 
# Need to run the FSW control test as a separate t.test because I'm comparing the UL quadrants instead of UR
MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_cell_type_AOV <- MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic %>%
  filter(Treat == "Beads to\n Hemocytes\n 1:1" | Treat ==  "Heat-Killed\n P. mar to\n Hemocytes\n 1:1"| Treat ==   "P. mar to\n Hemocytes\n 1:1" | 
           Treat == "P. mar to\n Hemocytes\n 5:1" | Treat == "P. mar to\n Hemocytes\n 10:1"| Treat == "P. mar to\n Hemocytes\n 25:1"|
           Treat ==  "Perkinsus\n Only\n Control" ) %>%
  filter(Gate == "Q28-UR" | Gate == "Q11-UR") %>%
  group_by(Treat) %>%
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Gate, data = .)))  %>%
  ungroup

MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_cell_type_AOV %>% filter(p.value <= 0.05)
#  Treat         p.value
#1 "Beads to Hemocytes 1:1" 0.00199
#2 "P. mar to Hemocytes 5:1"   0.0495 
#3 "P. mar to Hemocytes 10:1"  0.0127 

# FSW alone t.test - not running because so few cells were in these quadrants its not vaild to test for this

### Does Perkinsus cause significantly greater granular mitochondrial permeabilization than beads? ### - not correcting for amount of apoptosis caused by Perkinsus alone
# comparing only the 1:1
MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_beads_vs_hemo <- MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic %>% 
  filter(Treat == "Beads to\n Hemocytes\n 1:1"  | Treat == "P. mar to\n Hemocytes\n 1:1") %>% 
  # filter out just the granulocytes
  filter(Gate == "Q28-UR")

MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_beads_vs_hemo_ttest <- t.test(Percent_of_this_plot_arcsine ~ Treat, data = MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_beads_vs_hemo)
#  p-value = 0.003408

# test all concentrations with AOV
MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_beads_vs_hemo_aov <- MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic %>% 
  filter(Treat == "Beads to\n Hemocytes\n 1:1"  | Treat == "P. mar to\n Hemocytes\n 1:1" | Treat == "P. mar to\n Hemocytes\n 5:1" |
           Treat == "P. mar to\n Hemocytes\n 10:1" | Treat == "P. mar to\n Hemocytes\n 25:1") %>% 
  # filter out just the granulocytes
  filter(Gate == "Q28-UR")

MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_beads_vs_hemo_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data = MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_beads_vs_hemo_aov)
TukeyHSD(MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_beads_vs_hemo_aov)
#diff        lwr       upr     p adj
#P. mar to\n Hemocytes\n 1:1-Beads to\n Hemocytes\n 1:1     0.696106214  0.4598227 0.9323898 0.0000163
#P. mar to\n Hemocytes\n 5:1-Beads to\n Hemocytes\n 1:1     1.172514007  0.9362305 1.4087975 0.0000001
#P. mar to\n Hemocytes\n 10:1-Beads to\n Hemocytes\n 1:1    1.162857775  0.9265742 1.3991413 0.0000001
#P. mar to\n Hemocytes\n 25:1-Beads to\n Hemocytes\n 1:1    1.155823731  0.9195402 1.3921073 0.0000001
#P. mar to\n Hemocytes\n 5:1-P. mar to\n Hemocytes\n 1:1    0.476407793  0.2401243 0.7126913 0.0004321
#P. mar to\n Hemocytes\n 10:1-P. mar to\n Hemocytes\n 1:1   0.466751562  0.2304680 0.7030351 0.0005103
#P. mar to\n Hemocytes\n 25:1-P. mar to\n Hemocytes\n 1:1   0.459717517  0.2234340 0.6960011 0.0005769

# can't compare the the FSW because so few

# Does granular mitochondrial permeabilization differ between the different MOIs?
MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_MOI <- MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic %>% 
  filter(Treat == "P. mar to\n Hemocytes\n 1:1"  | Treat == "P. mar to\n Hemocytes\n 5:1" | Treat == "P. mar to\n Hemocytes\n 10:1" |
           Treat == "P. mar to\n Hemocytes\n 25:1") %>% 
  # filter out just the granulocytes
  filter(Gate == "Q28-UR")

MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_MOI_aov <- aov(Percent_of_this_plot_arcsine ~ Treat,
                                                                 data = MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_MOI)
TukeyHSD(MOI_1hr_2019_JC1_join_Agranular_Granular_apoptotic_MOI_aov)
#diff        lwr       upr     p adj
#P. mar to\n Hemocytes\n 5:1-P. mar to\n Hemocytes\n 1:1    0.476407793  0.2387422 0.7140733 0.0009278
#P. mar to\n Hemocytes\n 10:1-P. mar to\n Hemocytes\n 1:1   0.466751562  0.2290860 0.7044171 0.0010644
#P. mar to\n Hemocytes\n 25:1-P. mar to\n Hemocytes\n 1:1   0.459717517  0.2220520 0.6973831 0.0011779

#### PCA ANALYSIS to help decide 2018 sequencing samples ####

Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb

# arcsine transform the percentages for PCA analysis
Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb$arcsine_apop <- transf.arcsin(Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb$apoptosis_combined*0.01)

# subset for families A, B and L
Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb_subset <-  Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb %>%
  filter(Family =="A" | Family == "B" | Family == "L")

## Add information about which samples are going to be sequenced or not 
# only 6 samples are not going to be sequenced
Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb_subset <- Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb_subset %>%
  mutate(sequencing = case_when(
    ID == "A-189" | ID == "A-166"| ID == "B-141"| ID == "B-159"| ID == "L-168"| ID == "L-187" ~ "no",
    TRUE ~ "yes"))

#calculate PCA using, exclusing the caspase assay because samples are missing 
PCA_data <- prcomp(Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb_subset[c(5,8)], center=TRUE, scale=TRUE)
summary(PCA_data)

# plot by treatment
ggbiplot(PCA_data, varname.adjust= 0.1, varname.abbrev = TRUE) +
  geom_text(vjust="inward",hjust="inward", 
            label=Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb_subset$Treat)

# export plot 
ggsave(plot = last_plot(), path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       filename = "Day7_Day50_2018_pconc_Granular_apop_comb_PCA_TREAT.tiff", device = "tiff", units = "cm",
       width = 20, height = 20)

# plot by sequencing
ggbiplot(PCA_data, varname.adjust= 0.1, varname.abbrev = TRUE) +
  geom_text(vjust="inward",hjust="inward", 
            label=Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb_subset$sequencing)

# export plot 
ggsave(plot = last_plot(), path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       filename = "Day7_Day50_2018_pconc_Granular_apop_comb_PCA_SEQ.tiff", device = "tiff", units = "cm",
       width = 20, height = 20)

# plot by ave_log_pconc
ggbiplot(PCA_data, varname.adjust= 0.1, varname.abbrev = TRUE) +
  geom_text(vjust="inward",hjust="inward", 
            label=Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb_subset$ave_log_pconc)

# export plot 
ggsave(plot = last_plot(), path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       filename = "Day7_Day50_2018_pconc_Granular_apop_comb_PCA_PCONC.tiff", device = "tiff", units = "cm",
       width = 20, height = 20)

# plot by Family
ggbiplot(PCA_data, varname.adjust= 0.1, varname.abbrev = TRUE) +
  geom_text(vjust="inward",hjust="inward", 
            label=Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb_subset$Family)

# export plot 
ggsave(plot = last_plot(), path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       filename = "Day7_Day50_2018_pconc_Granular_apop_comb_PCA_FAM.tiff", device = "tiff", units = "cm",
       width = 20, height = 20)

# plot by apop
ggbiplot(PCA_data, varname.adjust= 0.1, varname.abbrev = TRUE) +
  geom_text(vjust="inward",hjust="inward", 
            label=Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb_subset$apoptosis_combined)

# export plot
ggsave(plot = last_plot(), path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       filename = "Day7_Day50_2018_pconc_Granular_apop_comb_PCA_APOP.tiff", device = "tiff", units = "cm",
       width = 20, height = 20)

## PCA for family B only, where they all have caspase samples

# subset for family B 
Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb_B <-  Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb %>%
  filter(Family == "B" )

#calculate PCA using, exclusing the caspase assay because samples are missing 
PCA_data_B <- prcomp(Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb_B[c(5,7,8)], center=TRUE, scale=TRUE)
summary(PCA_data_B)

# plot by treatment
ggbiplot(PCA_data_B, varname.adjust= 0.1, varname.abbrev = TRUE) +
  geom_text(vjust="inward",hjust="inward", 
            label=Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb_B$Treat)

# export plot 
ggsave(plot = last_plot(), path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       filename = "Day7_Day50_2018_pconc_Granular_apop_comb_PCA_B.tiff", device = "tiff", units = "cm",
       width = 20, height = 20)

# plot by caspase
ggbiplot(PCA_data_B, varname.adjust= 1, varname.abbrev = TRUE) +
  geom_text(vjust="inward",hjust="inward", 
            label=Day7_Day50_2018_all_assays_bad_removed_pconc_Granular_apop_casp_table_comb_B$caspase_combined)

# export plot 
ggsave(plot = last_plot(), path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       filename = "Day7_Day50_2018_pconc_Granular_apop_comb_PCA_casp.tiff", device = "tiff", units = "cm",
       width = 20, height = 20)

##### 2020 INHIBITOR CONCENTRATION ASSAYS #####

## LOAD ANALYSIS CSVs

## VIABILITY ASSAY

# Make new header column
Inhibitor_2020_VIA_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Inhibitor_Concentration_Assay/ANALYSIS_FILES/FORMATTED_CSVs/2020_Inhibitor_Concentration_VIA_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Inhibitor_2020_VIA <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Inhibitor_Concentration_Assay/ANALYSIS_FILES/FORMATTED_CSVs/2020_Inhibitor_Concentration_VIA_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(Inhibitor_2020_VIA_nms)

# Split column 1 by space and remove
Inhibitor_2020_VIA <- Inhibitor_2020_VIA %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
Inhibitor_2020_VIA <- Inhibitor_2020_VIA[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
Inhibitor_2020_VIA_percent <- Inhibitor_2020_VIA[,c(1:3,5,7,9,11,13,15,17,19,21,23)]
Inhibitor_2020_VIA_counts <- Inhibitor_2020_VIA[,c(1:4,6,8,10,12,14,16,18,20,22)]

# Gather count and percent columns separately
Inhibitor_2020_VIA_percent <- Inhibitor_2020_VIA_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:13))
Inhibitor_2020_VIA_counts <-  Inhibitor_2020_VIA_counts %>%  group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:13))

Inhibitor_2020_VIA_percent <- Inhibitor_2020_VIA_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Inhibitor_2020_VIA_counts <-  Inhibitor_2020_VIA_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

Inhibitor_2020_VIA_percent <- Inhibitor_2020_VIA_percent[,-c(4,8)]
Inhibitor_2020_VIA_counts <-  Inhibitor_2020_VIA_counts[,-c(4,8)]

# Full join together so all columns preserved
Inhibitor_2020_VIA_join <- full_join(Inhibitor_2020_VIA_percent, Inhibitor_2020_VIA_counts, by = c("ID","Treat","Assay","Plot_number","Gate", "Channel"))

# Add in cell type column based on plot number and gate
Inhibitor_2020_VIA_cell_type <- data.frame(Plot_number = c("8","8", "19","19", "20","20", "21","21", "22","22"), Gate= c("P4","P3","V9-L","V9-R", "V12-L","V12-R", "V5-L", "V5-R","V8-L","V8-R"), Cell_type=c(
  "granular","agranular", "agranular","live_agranular", "granular","live_granular", 
  "agranular","dead_agranular","granular","dead_granular"))
# Join cell type
Inhibitor_2020_VIA_join <- left_join(Inhibitor_2020_VIA_join, Inhibitor_2020_VIA_cell_type, by=c("Plot_number","Gate"))
unique_Inhibitor_2020_VIA_join_gate <- unique(Inhibitor_2020_VIA_join[,c("Plot_number","Gate")])

## APOPTOSIS ASSAY
# Make new header column
Inhibitor_2020_APOP_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Inhibitor_Concentration_Assay/ANALYSIS_FILES/FORMATTED_CSVs/2020_Inhibitor_concentration_APOP_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Inhibitor_2020_APOP <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Inhibitor_Concentration_Assay/ANALYSIS_FILES/FORMATTED_CSVs/2020_Inhibitor_concentration_APOP_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(Inhibitor_2020_APOP_nms)

# Split column 1 by space and remove
Inhibitor_2020_APOP <- Inhibitor_2020_APOP %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
Inhibitor_2020_APOP <- Inhibitor_2020_APOP[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
Inhibitor_2020_APOP_percent <- Inhibitor_2020_APOP[,c(1:3,5,7,9,11,13,15,17,19,21, 23,25,27,29,31)]
Inhibitor_2020_APOP_counts <-  Inhibitor_2020_APOP[,c(1:4,6,8,10,12,14,16,18,20,22,24,26,28,30)]

# Gather count and percent columns separately
Inhibitor_2020_APOP_percent <- Inhibitor_2020_APOP_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:17))
Inhibitor_2020_APOP_counts <-  Inhibitor_2020_APOP_counts %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:17))

Inhibitor_2020_APOP_percent <- Inhibitor_2020_APOP_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Inhibitor_2020_APOP_counts <-  Inhibitor_2020_APOP_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

Inhibitor_2020_APOP_percent <- Inhibitor_2020_APOP_percent[,-c(4,8)]
Inhibitor_2020_APOP_counts <-  Inhibitor_2020_APOP_counts[,-c(4,8)]

# Full join together so all columns preserved
Inhibitor_2020_APOP_join <- full_join(Inhibitor_2020_APOP_percent, Inhibitor_2020_APOP_counts, by = c("ID","Treat","Assay","Plot_number","Gate","Channel"))

# Add in cell type column based on plot number and gate
Inhibitor_2020_APOP_cell_type <- data.frame(Plot_number = c("2","2","6","6","8","8","9","9","9","9","11","11","11","11"), 
                                          Gate= c("V2-L","V2-R","V3-L","V3-R","P3","P4","Q11-UL","Q11-UR","Q11-LL","Q11-LR","Q8-UL","Q8-UR","Q8-LL","Q8-LR"), 
                                          Cell_type=c("non_apop_cells", "apop_cells", "live_cells", "dead_cells", "agranular", 
                                                      "granular", "dead_agranular", "dead_apoptotic_agranular", "unstained_agranular", "live_apoptotic_agranular", 
                                                      "dead_granular", "dead_apoptotic_granular", "unstained_granular", "live_apoptotic_granular"))

# Join cell type
Inhibitor_2020_APOP_join <- left_join(Inhibitor_2020_APOP_join, Inhibitor_2020_APOP_cell_type, by=c("Plot_number","Gate"))
Inhibitor_2020_APOP_join_unique <- unique(Inhibitor_2020_APOP_join[,c("Plot_number","Gate")])

## CASPASE ASSAY
# Make new header column
Inhibitor_2020_CASP_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Inhibitor_Concentration_Assay/ANALYSIS_FILES/FORMATTED_CSVs/2020_Inhibitor_concentration_CASP_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Inhibitor_2020_CASP <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Inhibitor_Concentration_Assay/ANALYSIS_FILES/FORMATTED_CSVs/2020_Inhibitor_concentration_CASP_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(Inhibitor_2020_CASP_nms)

# Split column 1 by space and remove
Inhibitor_2020_CASP <- Inhibitor_2020_CASP %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
Inhibitor_2020_CASP <- Inhibitor_2020_CASP[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
Inhibitor_2020_CASP_percent <- Inhibitor_2020_CASP[,c(1:3,5,7,9,11,13,15,17,19,21,23)]
Inhibitor_2020_CASP_counts <-  Inhibitor_2020_CASP[,c(1:4,6,8,10,12,14,16,18,20,22)]

# Gather count and percent columns separately
Inhibitor_2020_CASP_percent <- Inhibitor_2020_CASP_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:13))
Inhibitor_2020_CASP_counts <-  Inhibitor_2020_CASP_counts %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:13))

Inhibitor_2020_CASP_percent <- Inhibitor_2020_CASP_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Inhibitor_2020_CASP_counts <-  Inhibitor_2020_CASP_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

Inhibitor_2020_CASP_percent <- Inhibitor_2020_CASP_percent[,-c(4,8)]
Inhibitor_2020_CASP_counts <-  Inhibitor_2020_CASP_counts[,-c(4,8)]

# Full join together so all columns preserved
Inhibitor_2020_CASP_join <- full_join(Inhibitor_2020_CASP_percent, Inhibitor_2020_CASP_counts, by = c("ID","Treat","Assay","Plot_number","Gate","Channel"))

# Add in cell type column based on plot number and gate
Inhibitor_2020_CASP_cell_type <- data.frame(Plot_number = c("8","8","4","4","4","4","7","7","7","7"), 
                                        Gate= c("P3","P4","Q7-UL","Q7-UR", "Q7-LL", "Q7-LR", "Q3-UL","Q3-UR","Q3-LL","Q3-LR"), 
                                        Cell_type=c("agranular", 
                                                    "granular", "dead_agranular", "caspase_active_agranular", "unstained_agranular", "caspase_active_agranular", 
                                                     "dead_granular", "caspase_active_granular", "unstained_granular", "caspase_active_granular"))

# Join cell type
Inhibitor_2020_CASP_join <- left_join(Inhibitor_2020_CASP_join, Inhibitor_2020_CASP_cell_type, by=c("Plot_number","Gate"))
unique_Inhibitor_2020_CASP_join <- unique(Inhibitor_2020_CASP_join[,c("Plot_number","Gate")])

#### 2020 Inhibitor Concentration Arcsine Transformation of percentages ####

# Make list of dataframes
Inhibitor_2020_list <- list(Inhibitor_2020_VIA_join=Inhibitor_2020_VIA_join, 
                            Inhibitor_2020_APOP_join=Inhibitor_2020_APOP_join, 
                            Inhibitor_2020_CASP_join=Inhibitor_2020_CASP_join)

Inhibitor_2020_VIA_join
Inhibitor_2020_APOP_join
Inhibitor_2020_CASP_join

# Make new column and perform arcsine
Inhibitor_2020_VIA_join$Percent_of_this_plot_arcsine <-  transf.arcsin(Inhibitor_2020_VIA_join$Percent_of_this_plot)
Inhibitor_2020_APOP_join$Percent_of_this_plot_arcsine <- transf.arcsin(Inhibitor_2020_APOP_join$Percent_of_this_plot)
Inhibitor_2020_CASP_join$Percent_of_this_plot_arcsine <- transf.arcsin(Inhibitor_2020_CASP_join$Percent_of_this_plot)

# Change percent column to percent
Inhibitor_2020_VIA_join$Percent_of_this_plot <- Inhibitor_2020_VIA_join$Percent_of_this_plot*100
Inhibitor_2020_APOP_join$Percent_of_this_plot<- Inhibitor_2020_APOP_join$Percent_of_this_plot*100
Inhibitor_2020_CASP_join$Percent_of_this_plot<- Inhibitor_2020_CASP_join$Percent_of_this_plot*100

#### 2020 INHIBITOR CONCENTRATION PLOTS AND STATISTICS ####

### VIABILITY ASSAY ####

Inhibitor_2020_VIA_join

# set levels in correct order
Inhibitor_2020_VIA_join$Treat <- factor(Inhibitor_2020_VIA_join$Treat, levels= c("FSW_Control", "GDC_10_1HR","GDC_10_2HR","GDC_50_1HR","GDC_50_2HR",
                                                                                 "GDC_100_1HR", "GDC_100_2HR","ZVAD_10_1HR","ZVAD_10_2HR","ZVAD_50_1HR","ZVAD_50_2HR", 
                                                                                 "ZVAD_100_1HR", "ZVAD_100_2HR" ))

# Make combined ID_Treat column 
Inhibitor_2020_VIA_join_ID_treat <- Inhibitor_2020_VIA_join %>% mutate(ID_Treat = paste(ID,Treat, sep = "_"))

# Agranular and Granular cells Boxplot with significance bars, grouped by treatment color by gate(significant based on ANOVA)
Inhibitor_2020_VIA_join_Percent_Agranular_Granular <- Inhibitor_2020_VIA_join_ID_treat %>% filter(Gate =="P3" | Gate=="P4")
levels(factor(Inhibitor_2020_VIA_join_Percent_Agranular_Granular$Plot_number)) # 8 correct

Inhibitor_2020_VIA_join_Percent_Agranular_Granular_plot <- ggplot(data=Inhibitor_2020_VIA_join_Percent_Agranular_Granular,
                                                                aes(y=Percent_of_this_plot, x=Treat, fill=Gate)) + geom_col()  + 
  xlab("Individual") +
  ylab("Percent Hemocytes") +
  ggtitle("Percent of Granular and Agranular") + 
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, hjust = 1, angle = 90)) +
  theme(legend.text = element_text(size=12)) + 
  #  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  scale_fill_manual(name="Cell Type", labels=c("Agranular","Granular"), values=c("#6c81d9","#50b47b")) 

#save
ggsave(plot = Inhibitor_2020_VIA_join_Percent_Agranular_Granular_plot, device = "tiff", filename = "Inhibitor_2020_VIA_join_Percent_Agranular_Granular_plot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 10, width = 12)

# Granular Agranular plot formatted for multipanel figure 
# calculate sd and mean 
Inhibitor_2020_VIA_join_Percent_Agranular_Granular_sd <- Inhibitor_2020_VIA_join_Percent_Agranular_Granular %>%
  # get mean and sd
  dplyr::group_by(Gate) %>% mutate(mean = mean(Percent_of_this_plot), sd = sd(Percent_of_this_plot)) %>%
  # add inhibitor column for faceting
  mutate(inhibitor = case_when(grepl("GDC",Treat)~ "GDC",
                               grepl("Z",Treat)~ "ZVAD",
                               grepl("FSW", Treat)~"Control"))

Inhibitor_2020_VIA_join_Percent_Agranular_Granular_plot_multipanel <- 
  ggplot(data=Inhibitor_2020_VIA_join_Percent_Agranular_Granular_sd,
         aes(y=Percent_of_this_plot, x=Gate)) + geom_bar(aes(fill= Gate) ,position = position_dodge(), stat = "summary")  + 
  geom_point(aes(fill = Gate), size = 3, position=position_dodge(width = 0.9)) + 
  #facet_grid(.~inhibitor, scales = "free") +
  xlab(NULL) +
  ylab("% Hemocytes") + 
  theme_classic() +
  theme(text=element_text(size=12, face = "bold"), 
        axis.title.y=element_text(size=12, face = "bold"),
        axis.title.x=element_text(size=12,face = "bold"),
        axis.text.x=element_text(size=10,face = "bold", angle = 90, hjust = 1 ),
        axis.text.y=element_text(size=12,face = "bold")) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position=position_dodge(width = 0.9), width = 0.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  scale_x_discrete(labels = c("P3" = "Agranular","P4"="Granular")) + 
  scale_fill_manual(name="Cell Type", labels=c("Agranular","Granular"), values=c("#6c81d9","#50b47b")) 

Inhibitor_2020_VIA_join_Percent_Agranular_Granular_plot_multipanel <- Inhibitor_2020_VIA_join_Percent_Agranular_Granular_plot_multipanel +
    theme(axis.text.x = ggtext::element_markdown())
  
# Perform t.test and plot results onto boxplot 
stat.test <- as.data.frame(Inhibitor_2020_VIA_join_Percent_Agranular_Granular_sd) %>%
  # make sure to use the arcsine values 
  t_test(Percent_of_this_plot_arcsine ~ Gate) %>%
  add_significance(p.col = "p") %>% 
  add_xy_position(x = "Gate")

Inhibitor_2020_VIA_join_Percent_Agranular_Granular_plot_multipanel_sig <- 
  Inhibitor_2020_VIA_join_Percent_Agranular_Granular_plot_multipanel + stat_pvalue_manual(
    stat.test, label = "{p} {p.signif}",  tip.length = 0.02, y.position = 75)

# Compare the percent of dead cells between treatments
Inhibitor_2020_VIA_join_Percent_Agranular_Granular_DEAD_plot <- Inhibitor_2020_VIA_join_ID_treat %>%
  # filter for granular, V7 from plot 20
  filter(Gate=="V5-R" |Gate=="V8-R") %>%
  ggplot(aes(y=Percent_of_this_plot, x=Treat, fill=Gate))  + geom_boxplot() +
  xlab("Individual") +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Live and Dead Granular Hemocytes") + 
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, hjust = 1, angle = 90)) +
  theme(legend.text = element_text(size=12)) + 
  #  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  scale_fill_manual(name="Cell Type", labels=c("Percent Dead Agranular","Percent Dead Granular"), values=c("#5b2c90", "#88bf3b")) 

#"#cc57b4", "#88bf3b", "#aa4dce", "#56b464", "#5b2c90"
#"#c89832", "#5a6ee6", "#ca4e33", "#7e78d4", "#cd4272"

#save
ggsave(plot = Inhibitor_2020_VIA_join_Percent_Agranular_Granular_DEAD_plot , device = "tiff", filename = "Inhibitor_2020_VIA_join_Percent_Agranular_Granular_DEAD_plot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 10, width = 12)

## Calculate live cells using dead cell values and plot 

# Plot Dead granular hemocytes as live and prepare for multi-panel figure 
Inhibitor_2020_VIA_join_Percent_Agranular_Granular_LIVE_sd <- Inhibitor_2020_VIA_join_ID_treat %>%
  # filter for granular
  filter(Gate=="V8-R") %>%
  dplyr::group_by(Gate,Treat) %>%
  mutate(Percent_of_this_plot_live = (100 - Percent_of_this_plot)) %>%
  # get mean and sd
  mutate(mean = mean(Percent_of_this_plot_live), sd = sd(Percent_of_this_plot_live)) %>%
  mutate(inhibitor = case_when(grepl("GDC",Treat)~ "GDC",
                               grepl("Z",Treat)~ "ZVAD",
                               grepl("FSW", Treat)~"Control"))

# add arcine transformed percentages 
Inhibitor_2020_VIA_join_Percent_Agranular_Granular_LIVE_sd$Percent_of_this_plot_live_arcsine <- transf.arcsin(Inhibitor_2020_VIA_join_Percent_Agranular_Granular_LIVE_sd$Percent_of_this_plot_live*0.01)

# Make plot of live hemocytes for multipanel figure
Inhibitor_2020_VIA_join_Percent_Agranular_Granular_LIVE_sd_multipanel <- 
  ggplot(data=  Inhibitor_2020_VIA_join_Percent_Agranular_Granular_LIVE_sd,
         aes(y=Percent_of_this_plot_live, x=Treat)) + geom_bar(aes(fill= inhibitor),position = "dodge", stat = "summary")  + 
  geom_point(shape = 15, aes(fill = inhibitor)) + 
  #facet_grid(.~Gate) +
  xlab(NULL) +
  ylab("% Live Granular Hemocytes") + 
  theme_classic() +
  theme(text=element_text(size=12, face = "bold"), 
        axis.title.y=element_text(size=12, face = "bold"),
        axis.title.x=element_text(size=12,face = "bold"),
        axis.text.x=element_text(size=10,face = "bold", angle = 90, hjust = 1),
        axis.text.y=element_text(size=12,face = "bold")) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,105), breaks = c(0,25,50,75,100)) +
  scale_x_discrete(labels = c( "FSW_Control" = "Control", 
                               "GDC_10_1HR"  = "GDC-0152<br>10um 1hr" ,
                               "GDC_10_2HR"  = "GDC-0152<br>10um 2hr" ,
                               "GDC_50_1HR"  = "GDC-0152<br>50um 1hr" ,
                               "GDC_50_2HR"  = "GDC-0152<br>50um 2hr" ,
                               "GDC_100_1HR" = "GDC-0152<br>100um 1hr", 
                               "GDC_100_2HR" = "GDC-0152<br>100um 2hr", 
                               "ZVAD_10_1HR" = "Z-VAD-fmk<br>10um 1hr", 
                               "ZVAD_10_2HR" = "Z-VAD-fmk<br>10um 2hr", 
                               "ZVAD_50_1HR" = "Z-VAD-fmk<br>50um 1hr", 
                               "ZVAD_50_2HR" = "Z-VAD-fmk<br>50um 2hr",
                               "ZVAD_100_1HR"= "Z-VAD-fmk<br>100um 1hr", 
                               "ZVAD_100_2HR"= "Z-VAD-fmk<br>100um 2hr")) + 
  scale_fill_manual(name="Treatment", labels=c("Control","GDC-0152","Z-VAD-fmk"), values=c("#b94973", "#45c097","#9fac3a")) 

Inhibitor_2020_VIA_join_Percent_Agranular_Granular_LIVE_sd_multipanel <- Inhibitor_2020_VIA_join_Percent_Agranular_Granular_LIVE_sd_multipanel + 
  theme(axis.text.x = ggtext::element_markdown())

# Perform aov and plot results onto boxplot #
Inhibitor_2020_VIA_join_Percent_Agranular_Granular_LIVE_sd_AOV <-  aov(Percent_of_this_plot_live_arcsine ~ Treat, Inhibitor_2020_VIA_join_Percent_Agranular_Granular_LIVE_sd)

stat_test_tukey <- tukey_hsd(Inhibitor_2020_VIA_join_Percent_Agranular_Granular_LIVE_sd_AOV) %>%
  add_significance(p.col = "p.adj")  

# take only the significant columns
stat_test_tukey <- stat_test_tukey %>% filter(p.adj <= 0.05) %>% filter(group1 == "FSW_Control")

Inhibitor_2020_VIA_join_Percent_Agranular_Granular_LIVE_sd_multipanel_sig <- 
  Inhibitor_2020_VIA_join_Percent_Agranular_Granular_LIVE_sd_multipanel + stat_pvalue_manual(
    stat_test_tukey, label = "{p.adj} {p.adj.signif}",  tip.length = 0.02, y.position = c(97,101),
    size = 4, vjust = 0) + labs(subtitle = "Tukey HSD, Arcsine Percent ~ Treat")

## ANALYSIS

# Do treatments increase cell death as compared to the control 
Inhibitor_2020_VIA_join_ID_treat_gran_dead <- Inhibitor_2020_VIA_join_ID_treat %>%
  filter(Gate=="V8-R")

Inhibitor_2020_VIA_join_ID_treat_gran_dead_aov <- aov(Percent_of_this_plot_arcsine ~ Treat, data = Inhibitor_2020_VIA_join_ID_treat_gran_dead)
Inhibitor_2020_VIA_join_ID_treat_gran_dead_aov_tukey <- TukeyHSD(Inhibitor_2020_VIA_join_ID_treat_gran_dead_aov)

  # GDC_10_1HR-FSW_Control    0.9997797
  # GDC_10_2HR-FSW_Control    0.7479487
  # GDC_50_1HR-FSW_Control    0.2079920
  # GDC_50_2HR-FSW_Control    0.3509844
  # GDC_100_1HR-FSW_Control   0.0210438
  # GDC_100_2HR-FSW_Control   0.0036492
  # ZVAD_10_1HR-FSW_Control   0.9893511
  # ZVAD_10_2HR-FSW_Control   0.7547903
  # ZVAD_50_1HR-FSW_Control   0.9236884
  # ZVAD_50_2HR-FSW_Control   0.6614097
  # ZVAD_100_1HR-FSW_Control  0.9699563
  # ZVAD_100_2HR-FSW_Control  0.4637972


### Apoptosis Assay Statistics and Plotting ####

Inhibitor_2020_APOP_join
# set levels in correct order
Inhibitor_2020_APOP_join$Treat <- factor(Inhibitor_2020_APOP_join$Treat, levels= c("FSW_Control", "GDC_10_1HR","GDC_10_2HR","GDC_50_1HR","GDC_50_2HR",
                                                                                 "GDC_100_1HR", "GDC_100_2HR","ZVAD_10_1HR","ZVAD_10_2HR","ZVAD_50_1HR","ZVAD_50_2HR", 
                                                                                 "ZVAD_100_1HR", "ZVAD_100_2HR" ))
# Make combined ID_Treat column 
Inhibitor_2020_APOP_join_ID_treat <- Inhibitor_2020_APOP_join %>% mutate(ID_Treat = paste(ID,Treat, sep = "_"))

# Combine live and dead apoptotic cells 
# Make APOP combined gate for agranular and granular separately 
Inhibitor_2020_APOP_join_ID_treat_Granular_Apop_combined <- Inhibitor_2020_APOP_join_ID_treat  %>% filter(Gate == "Q8-LR" | Gate == "Q8-UR") %>% group_by(ID_Treat, Treat) %>% dplyr::summarise(Percent_of_this_plot = sum(Percent_of_this_plot))
Inhibitor_2020_APOP_join_ID_treat_Agranular_Apop_combined <- Inhibitor_2020_APOP_join_ID_treat  %>% filter(Gate == "Q11-LR" | Gate == "Q11-UR") %>% group_by(ID_Treat, Treat) %>% dplyr::summarise(Percent_of_this_plot = sum(Percent_of_this_plot))

# Add new gate name for each 
Inhibitor_2020_APOP_join_ID_treat_Granular_Apop_combined$Gate <- "apop_combined_granular"
Inhibitor_2020_APOP_join_ID_treat_Agranular_Apop_combined$Gate <- "apop_combined_agranular"

# Combined data frames for each cell type
Inhibitor_2020_APOP_join_ID_treat_Granular_Agranular_Apop_combined <- rbind(Inhibitor_2020_APOP_join_ID_treat_Granular_Apop_combined, Inhibitor_2020_APOP_join_ID_treat_Agranular_Apop_combined)

# Add arcsine transformed data
Inhibitor_2020_APOP_join_ID_treat_Granular_Agranular_Apop_combined$Percent_of_this_plot_arcsine <- transf.arcsin(Inhibitor_2020_APOP_join_ID_treat_Granular_Agranular_Apop_combined$Percent_of_this_plot*0.01)

# Merge with original data frame
Inhibitor_2020_APOP_join_ID_treat_combined <- full_join(Inhibitor_2020_APOP_join_ID_treat,Inhibitor_2020_APOP_join_ID_treat_Granular_Agranular_Apop_combined , by =c("Gate", "ID_Treat","Treat", "Percent_of_this_plot","Percent_of_this_plot_arcsine"))


# Plot granular and agranluar combined apoptosis (both live and dead apoptotic cells)
# Make plot 
Inhibitor_2020_APOP_join_ID_treat_combined_apoptotic <- Inhibitor_2020_APOP_join_ID_treat_combined %>%
  filter(Gate == "apop_combined_granular" | Gate == "apop_combined_agranular") %>%
    ggplot(aes(y=Percent_of_this_plot, x=Treat, color=Gate)) + 
  geom_point(position=position_dodge(width=0.75)) + 
  geom_boxplot() +
  xlab("Treatment") +
  ylab("Percent of Each Cell Type") + 
  ggtitle("Percent of Granular and Agranular Apoptotic Cells") + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,20), breaks = c(0,5,10,15,20)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, hjust = 1, angle = 90)) +
  theme(legend.text = element_text(size=12)) +
  scale_color_manual(name="Cell Type", labels=c("Apoptotic Agranular", "Apoptotic Granular"), 
                     values = c("#cc57b4", "#7e78d4")) 

# color options
#"#cc57b4", "#88bf3b", "#aa4dce", "#56b464", "#5b2c90"
#"#c89832", "#5a6ee6", "#ca4e33", "#7e78d4", "#cd4272"

#save
ggsave(plot = Inhibitor_2020_APOP_join_ID_treat_combined_apoptotic, device = "tiff", filename = "Inhibitor_2020_APOP_join_ID_treat_combined_apoptotic.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 8)

## Plot just the granular apoptosis for multi-panel figure 
Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd <- Inhibitor_2020_APOP_join_ID_treat_combined %>%
  filter(Gate == "apop_combined_granular" ) %>%
  group_by(Treat) %>% 
  # get mean and sd
  mutate(mean = mean(Percent_of_this_plot), sd = sd(Percent_of_this_plot)) %>%
  mutate(inhibitor = case_when(grepl("GDC",Treat)~ "GDC_0152",
                               grepl("Z",Treat)~ "Z_VAD_fmk",
                               grepl("FSW", Treat)~"Control")) %>%
  mutate(time = case_when(
    grepl("1HR",Treat) ~ "1hr",
        grepl("2HR", Treat)~ "2hr",
    grepl("FSW",Treat)~"Control")) %>%
  mutate(level = case_when(
    grepl("10_",Treat) ~ "10",
    grepl("50_",Treat) ~ "50",
    grepl("100_",Treat) ~ "100",
    grepl("FSW",Treat )~ "Control"
  ))

# Make plot 
Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_multipanel <- 
  ggplot(data=Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd ,
         aes(y=Percent_of_this_plot, x=Treat)) + geom_bar(aes(fill= inhibitor),position = "dodge", stat = "summary")  + 
  geom_point(shape = 15, aes(fill = inhibitor)) + 
  #facet_grid(.~Gate) +
  xlab(NULL) +
  ylab("% Apoptotic Granular Hemocytes") + 
  theme_classic() +
  theme(text=element_text(size=12, face = "bold"), 
        axis.title.y=element_text(size=12, face = "bold"),
        axis.title.x=element_text(size=12,face = "bold"),
        axis.text.x=element_text(size=10,face = "bold", angle = 90, hjust = 1),
        axis.text.y=element_text(size=12,face = "bold")) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,25)) +
  scale_x_discrete(labels = c( "FSW_Control" = "Control", 
                               "GDC_10_1HR"  = "GDC-0152<br>10um 1hr" ,
                               "GDC_10_2HR"  = "GDC-0152<br>10um 2hr" ,
                               "GDC_50_1HR"  = "GDC-0152<br>50um 1hr" ,
                               "GDC_50_2HR"  = "GDC-0152<br>50um 2hr" ,
                               "GDC_100_1HR" = "GDC-0152<br>100um 1hr", 
                               "GDC_100_2HR" = "GDC-0152<br>100um 2hr", 
                               "ZVAD_10_1HR" = "Z-VAD-fmk<br>10um 1hr", 
                               "ZVAD_10_2HR" = "Z-VAD-fmk<br>10um 2hr", 
                               "ZVAD_50_1HR" = "Z-VAD-fmk<br>50um 1hr", 
                               "ZVAD_50_2HR" = "Z-VAD-fmk<br>50um 2hr",
                               "ZVAD_100_1HR"= "Z-VAD-fmk<br>100um 1hr", 
                               "ZVAD_100_2HR"= "Z-VAD-fmk<br>100um 2hr")) + 
  scale_fill_manual(name="Treatment", labels=c("Control","GDC-0152","Z-VAD-fmk"), values=c("#b94973", "#45c097","#9fac3a")) 

Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_multipanel <- Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_multipanel + 
  theme(axis.text.x = ggtext::element_markdown())

# Perform aov and plot results onto barplot 
Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_AOV <- aov(Percent_of_this_plot_arcsine ~ Treat, Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd)

stat.test_tukey <- 
  tukey_hsd(Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_AOV) %>%
  add_significance(p.col = "p.adj")

# take only significant columns
stat.test_tukey <- stat.test_tukey %>% filter(p.adj <= 0.05)

Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_multipanel_sig <- 
  Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_multipanel + stat_pvalue_manual(
    stat.test_tukey, label = "{p.adj} {p.adj.signif}",  tip.length = 0.02, y.position = c(24),
    size = 4, vjust = 0) + labs(subtitle = "Tukey HSD, Arcsine Percent ~ Treat")

## Plot with treatments combined and bars by inhibitor
Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_inhibitor <- Inhibitor_2020_APOP_join_ID_treat_combined %>%
  filter(Gate == "apop_combined_granular" ) %>%
  mutate(inhibitor = case_when(grepl("GDC",Treat)~ "GDC_0152",
                               grepl("Z",Treat)~ "Z_VAD_fmk",
                               grepl("FSW", Treat)~"Control")) %>% 
  ungroup() %>%
  group_by(inhibitor) %>% 
  # get mean and sd
  mutate(mean = mean(Percent_of_this_plot), sd = sd(Percent_of_this_plot))

Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_multipanel_inhibitor <- 
  ggplot(data=Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_inhibitor,
         aes(y=Percent_of_this_plot, x=inhibitor)) + geom_bar(aes(fill= inhibitor),position = "dodge", stat = "summary")  + 
  geom_point(shape = 15, aes(fill = inhibitor)) + 
  #facet_grid(.~Gate) +
  xlab(NULL) +
  ylab("% Apoptotic Granular Hemocytes") + 
  theme_classic() +
  theme(text=element_text(size=12, face = "bold"), 
        axis.title.y=element_text(size=12, face = "bold"),
        axis.title.x=element_text(size=12,face = "bold"),
        axis.text.x=element_text(size=10,face = "bold", angle = 90, hjust = 1),
        axis.text.y=element_text(size=12,face = "bold")) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,25)) +
  scale_x_discrete(labels = c( "GDC_0152" = "GDC-0152","Control"="Control","Z_VAD_fmk"="Z-VAD-fmk")) + 
  scale_fill_manual(name="Treatment", labels=c("Control","GDC-0152","Z-VAD-fmk"), values=c("#b94973", "#45c097","#9fac3a")) 

# Perform aov and plot results onto barplot 
Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_AOV <- aov(Percent_of_this_plot_arcsine ~ inhibitor, Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_inhibitor)

stat.test_tukey <- 
  tukey_hsd(Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_AOV) %>%
  add_significance(p.col = "p.adj")

# take only significant columns
stat.test_tukey <- stat.test_tukey %>% filter(p.adj <= 0.05)

Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_multipanel_inhibitor_sig <- 
  Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_multipanel_inhibitor + stat_pvalue_manual(
    stat.test_tukey, label = "{p.adj} {p.adj.signif}",  tip.length = 0.02, y.position = c(24,25),
    size = 4, vjust = 0) + labs(subtitle = "Tukey HSD, Arcsine Percent ~ inhibitor")


### ANALYSIS QUESTIONS ###

### Do cell types differ in apoptosis levels? 
Inhibitor_2020_APOP_join_ID_treat_combined_apoptotic_AOV <- Inhibitor_2020_APOP_join_ID_treat_combined %>%
  filter(Gate == "apop_combined_granular" | Gate == "apop_combined_agranular") %>%
  group_by(Treat) %>%
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Gate, data = .)))  %>%
  ungroup

Inhibitor_2020_APOP_join_ID_treat_combined_apoptotic_AOV %>% filter(p.value <= 0.05)
# granular have more apoptosis in all categories

## Is there a difference in granular between treatments?
Inhibitor_2020_APOP_join_ID_treat_combined_granular <- Inhibitor_2020_APOP_join_ID_treat_combined %>% filter(Gate == "apop_combined_granular")

Inhibitor_2020_APOP_join_ID_treat_combined_apoptotic_TREAT_AOV <-  aov(Percent_of_this_plot_arcsine ~ Treat, data = Inhibitor_2020_APOP_join_ID_treat_combined_granular)
TukeyHSD(Inhibitor_2020_APOP_join_ID_treat_combined_apoptotic_TREAT_AOV)

#GDC_10_1HR-FSW_Control   0.1643456
#GDC_10_2HR-FSW_Control   0.9068076
#GDC_50_1HR-FSW_Control   0.3184684
#GDC_50_2HR-FSW_Control   0.0759855
#GDC_100_1HR-FSW_Control  0.1818360
#GDC_100_2HR-FSW_Control  0.0454713
#ZVAD_10_1HR-FSW_Control  0.8078644
#ZVAD_10_2HR-FSW_Control  0.2705107
#ZVAD_50_1HR-FSW_Control  0.4192132
#ZVAD_50_2HR-FSW_Control  0.2406204
#ZVAD_100_1HR-FSW_Control 0.5096372
#ZVAD_100_2HR-FSW_Control 0.2655081

### Caspase Assay Statistics and Plotting ####

Inhibitor_2020_CASP_join
# set levels in correct order
Inhibitor_2020_CASP_join$Treat <- factor(Inhibitor_2020_CASP_join$Treat, levels= c("FSW_Control", "GDC_10_1HR","GDC_10_2HR","GDC_50_1HR","GDC_50_2HR",
                                                                                   "GDC_100_1HR", "GDC_100_2HR","ZVAD_10_1HR","ZVAD_10_2HR","ZVAD_50_1HR","ZVAD_50_2HR", 
                                                                                   "ZVAD_100_1HR", "ZVAD_100_2HR" ))
# Make combined ID_Treat column 
Inhibitor_2020_CASP_join_ID_treat <- Inhibitor_2020_CASP_join %>% mutate(ID_Treat = paste(ID,Treat, sep = "_"))

# Combine live and dead apoptotic cells 
# Make APOP combined gate for agranular and granular separately 
Inhibitor_2020_CASP_join_ID_treat_Granular_combined <-  Inhibitor_2020_CASP_join_ID_treat  %>% filter(Gate == "Q3-LR" | Gate == "Q3-UR") %>% group_by(ID_Treat, Treat) %>% dplyr::summarise(Percent_of_this_plot = sum(Percent_of_this_plot))
Inhibitor_2020_CASP_join_ID_treat_Agranular_combined <- Inhibitor_2020_CASP_join_ID_treat  %>% filter(Gate == "Q7-LR" | Gate == "Q7-UR") %>% group_by(ID_Treat, Treat) %>% dplyr::summarise(Percent_of_this_plot = sum(Percent_of_this_plot))

# Add new gate name for each 
Inhibitor_2020_CASP_join_ID_treat_Granular_combined$Gate <- "casp_active_combined_granular"
Inhibitor_2020_CASP_join_ID_treat_Agranular_combined$Gate <- "casp_active_combined_agranular"

# Combined data frames for each cell type
Inhibitor_2020_CASP_join_ID_treat_Granular_Agranular_combined <- rbind(Inhibitor_2020_CASP_join_ID_treat_Granular_combined, Inhibitor_2020_CASP_join_ID_treat_Agranular_combined)

# Add arcsine transformed data
Inhibitor_2020_CASP_join_ID_treat_Granular_Agranular_combined$Percent_of_this_plot_arcsine <- transf.arcsin(Inhibitor_2020_CASP_join_ID_treat_Granular_Agranular_combined$Percent_of_this_plot*0.01)

# Merge with original data frame
Inhibitor_2020_CASP_join_ID_treat_combined <- full_join(Inhibitor_2020_CASP_join_ID_treat,Inhibitor_2020_CASP_join_ID_treat_Granular_Agranular_combined, by =c("Gate", "ID_Treat","Treat", "Percent_of_this_plot","Percent_of_this_plot_arcsine"))


# Plot granular and agranluar combined apoptosis (both live and dead apoptotic cells)

# Make plot 
Inhibitor_2020_CASP_join_ID_treat_combined_active <- Inhibitor_2020_CASP_join_ID_treat_combined %>%
  filter(Gate == "casp_active_combined_granular" | Gate == "casp_active_combined_agranular") %>%
  ggplot(aes(y=Percent_of_this_plot, x=Treat, color=Gate)) + 
  geom_point(position=position_dodge(width=0.75)) + 
  geom_boxplot() +
  xlab("Treatment") +
  ylab("Percent of Each Cell Type") + 
  ggtitle("Percent of Granular and Agranular Caspase 3/7 Active Cells") + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,20), breaks = c(0,5,10,15,20)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, hjust = 1, angle = 90)) +
  theme(legend.text = element_text(size=12)) +
  scale_color_manual(name="Cell Type", labels=c("Caspase 3/7 Active Agranular", "Caspase 3/7 Active Granular"), 
                     values = c("#cc57b4", "#7e78d4")) 

# color options
#"#cc57b4", "#88bf3b", "#aa4dce", "#56b464", "#5b2c90"
#"#c89832", "#5a6ee6", "#ca4e33", "#7e78d4", "#cd4272"

#save
ggsave(plot = Inhibitor_2020_CASP_join_ID_treat_combined_active, device = "tiff", filename = "Inhibitor_2020_CASP_join_ID_treat_combined_active.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 8)

##Plot Caspase active granulocytes in format for multipanel figure with multiple comparisons run 
## Plot just the granular apoptosis for multi-panel figure 
Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd <- Inhibitor_2020_CASP_join_ID_treat_combined %>%
  filter(Gate == "casp_active_combined_granular" ) %>%
  group_by(Treat) %>% 
  # get mean and sd
  mutate(mean = mean(Percent_of_this_plot), sd = sd(Percent_of_this_plot)) %>%
  mutate(inhibitor = case_when(grepl("GDC",Treat)~ "GDC_0152",
                               grepl("Z",Treat)~ "Z_VAD_fmk",
                               grepl("FSW", Treat)~"Control")) %>%
  mutate(time = case_when(
    grepl("1HR",Treat) ~ "1hr",
    grepl("2HR", Treat)~ "2hr",
    grepl("FSW",Treat)~"Control")) %>%
  mutate(level = case_when(
    grepl("10_",Treat) ~ "10",
    grepl("50_",Treat) ~ "50",
    grepl("100_",Treat) ~ "100",
    grepl("FSW",Treat )~ "Control"
  ))

# Make plot 
Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_multipanel <- 
  ggplot(data=Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd ,
         aes(y=Percent_of_this_plot, x=Treat)) + geom_bar(aes(fill= inhibitor),position = "dodge", stat = "summary")  + 
  geom_point(shape = 15, aes(fill = inhibitor)) + 
  #facet_grid(.~Gate) +
  xlab(NULL) +
  ylab("% Caspase 3/7 Active Granular Hemocytes") + 
  theme_classic() +
  theme(text=element_text(size=12, face = "bold"), 
        axis.title.y=element_text(size=12, face = "bold"),
        axis.title.x=element_text(size=12,face = "bold"),
        axis.text.x=element_text(size=10,face = "bold", angle = 90, hjust = 1),
        axis.text.y=element_text(size=12,face = "bold")) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,30), breaks = c(0,5,15,20,25,30)) +
  scale_x_discrete(labels = c( "FSW_Control" = "Control", 
                               "GDC_10_1HR"  = "GDC-0152<br>10um 1hr" ,
                               "GDC_10_2HR"  = "GDC-0152<br>10um 2hr" ,
                               "GDC_50_1HR"  = "GDC-0152<br>50um 1hr" ,
                               "GDC_50_2HR"  = "GDC-0152<br>50um 2hr" ,
                               "GDC_100_1HR" = "GDC-0152<br>100um 1hr", 
                               "GDC_100_2HR" = "GDC-0152<br>100um 2hr", 
                               "ZVAD_10_1HR" = "Z-VAD-fmk<br>10um 1hr", 
                               "ZVAD_10_2HR" = "Z-VAD-fmk<br>10um 2hr", 
                               "ZVAD_50_1HR" = "Z-VAD-fmk<br>50um 1hr", 
                               "ZVAD_50_2HR" = "Z-VAD-fmk<br>50um 2hr",
                               "ZVAD_100_1HR"= "Z-VAD-fmk<br>100um 1hr", 
                               "ZVAD_100_2HR"= "Z-VAD-fmk<br>100um 2hr")) + 
  scale_fill_manual(name="Treatment", labels=c("Control","GDC-0152","Z-VAD-fmk"), values=c("#b94973", "#45c097","#9fac3a")) 

Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_multipanel <- Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_multipanel + 
  theme(axis.text.x = ggtext::element_markdown())

# Perform aov and plot results onto barplot 
Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_AOV <- aov(Percent_of_this_plot_arcsine ~ Treat, Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd)

stat.test_tukey <- 
  tukey_hsd(Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_AOV) %>%
  add_significance(p.col = "p.adj")

# take only significant columns
stat.test_tukey <- stat.test_tukey %>% filter(p.adj <= 0.05)

Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_multipanel_sig <- 
  Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_multipanel + stat_pvalue_manual(
    stat.test_tukey, label = "{p.adj} {p.adj.signif}",  tip.length = 0.02, y.position = c(21,22,23,24,25,26,27),
    size = 4, vjust = 0) + labs(subtitle = "Tukey HSD, Arcsine Percent ~ Treat")

## Plot with treatments combined and bars by inhibitor
Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_inhibitor <- Inhibitor_2020_CASP_join_ID_treat_combined %>%
  filter(Gate == "casp_active_combined_granular" ) %>%
  mutate(inhibitor = case_when(grepl("GDC",Treat)~ "GDC_0152",
                               grepl("Z",Treat)~ "Z_VAD_fmk",
                               grepl("FSW", Treat)~"Control")) %>% 
  ungroup() %>%
  group_by(inhibitor) %>% 
  # get mean and sd
  mutate(mean = mean(Percent_of_this_plot), sd = sd(Percent_of_this_plot))

Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_multipanel_inhibitor <- 
  ggplot(data=Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_inhibitor,
         aes(y=Percent_of_this_plot, x=inhibitor)) + geom_bar(aes(fill= inhibitor),position = "dodge", stat = "summary")  + 
  geom_point(shape = 15, aes(fill = inhibitor)) + 
  #facet_grid(.~Gate) +
  xlab(NULL) +
  ylab("% Caspase 3/7 Active Granular Hemocytes") + 
  theme_classic() +
  theme(text=element_text(size=12, face = "bold"), 
        axis.title.y=element_text(size=12, face = "bold"),
        axis.title.x=element_text(size=12,face = "bold"),
        axis.text.x=element_text(size=10,face = "bold", angle = 90, hjust = 1),
        axis.text.y=element_text(size=12,face = "bold")) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,25)) +
  scale_x_discrete(labels = c("Control"="Control","GDC_0152"="GDC-0152","Z_VAD_fmk"="Z-VAD-fmk")) + 
scale_fill_manual(name="Treatment", labels=c("Control","GDC-0152","Z-VAD-fmk"), values=c("#b94973", "#45c097","#9fac3a")) 

# Perform aov and plot results onto barplot 
Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_AOV <- aov(Percent_of_this_plot_arcsine ~ inhibitor, Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_inhibitor)

stat.test_tukey <- 
  tukey_hsd(Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_AOV) %>%
  add_significance(p.col = "p.adj")

# take only significant columns
stat.test_tukey <- stat.test_tukey %>% filter(p.adj <= 0.05)

Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_multipanel_inhibitor_sig <- 
  Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_multipanel_inhibitor + stat_pvalue_manual(
    stat.test_tukey, label = "{p.adj} {p.adj.signif}",  tip.length = 0.02, y.position = c(24),
    size = 4, vjust = 0) + labs(subtitle = "Tukey HSD, Arcsine Percent ~ inhibitor")


### ANALYSIS QUESTIONS ###

### Do cell types differ in apoptosis levels? 
Inhibitor_2020_CASP_join_ID_treat_combined_AOV <- Inhibitor_2020_CASP_join_ID_treat_combined %>%
  filter(Gate == "casp_active_combined_granular" | Gate == "casp_active_combined_agranular") %>%
  group_by(Treat) %>%
  do(broom::tidy(aov(Percent_of_this_plot_arcsine ~ Gate, data = .)))  %>%
  ungroup

Inhibitor_2020_CASP_join_ID_treat_combined_AOV %>% filter(p.value <= 0.05)
# granular have more caspase activation in all categories

## Is there a difference in granular between treatments?
Inhibitor_2020_CASP_join_ID_treat_combined_granular <- Inhibitor_2020_CASP_join_ID_treat_combined %>% filter(Gate == "casp_active_combined_granular")

Inhibitor_2020_CASP_join_ID_treat_combined_granular_TREAT_AOV <-  aov(Percent_of_this_plot_arcsine ~ Treat, data = Inhibitor_2020_CASP_join_ID_treat_combined_granular)
TukeyHSD(Inhibitor_2020_CASP_join_ID_treat_combined_granular_TREAT_AOV)

#### 2020 Inhibitor Assay Multipanel figure ####

Inhibitor_2020_multipanel <- cowplot::plot_grid(Inhibitor_2020_VIA_join_Percent_Agranular_Granular_plot_multipanel_sig,Inhibitor_2020_VIA_join_Percent_Agranular_Granular_LIVE_sd_multipanel_sig, 
                                                 Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_multipanel_inhibitor_sig, Inhibitor_2020_APOP_join_ID_treat_combined_granular_sd_multipanel_sig, 
                                                 Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_multipanel_inhibitor_sig, Inhibitor_2020_CASP_join_ID_treat_combined_granular_sd_multipanel_sig,
                                                 ncol = 2, nrow = 3, labels = "AUTO", label_size = 16, label_fontface = "bold", 
                                                rel_widths = c(0.5,1))

ggsave(Inhibitor_2020_multipanel, device = "tiff", filename = "Inhibitor_2020_multipanel.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 17, width = 12 ) 


#### 2020 Dermo and Inhibitor Experiment Load Data ####

## LOAD ANALYSIS CSVs

## VIABILITY ASSAY
# First load the Viability data from the hemocytes 
# Make new header column
Dermo_Inhibitor_hemo_2020_VIA_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/FORMATTED_CSVs/2020_Dermo_inhibitors_hemocytes_VIA_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Dermo_Inhibitor_hemo_2020_VIA <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/FORMATTED_CSVs/2020_Dermo_inhibitors_hemocytes_VIA_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(Dermo_Inhibitor_hemo_2020_VIA_nms)

# Split column 1 by space and remove
Dermo_Inhibitor_hemo_2020_VIA <- Dermo_Inhibitor_hemo_2020_VIA %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
Dermo_Inhibitor_hemo_2020_VIA <- Dermo_Inhibitor_hemo_2020_VIA[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
Dermo_Inhibitor_hemo_2020_VIA_percent <- Dermo_Inhibitor_hemo_2020_VIA[,c(1:3,5,7,9,11,13,15,17,19,21,23)]
Dermo_Inhibitor_hemo_2020_VIA_counts <- Dermo_Inhibitor_hemo_2020_VIA[,c(1:4,6,8,10,12,14,16,18,20,22)]

# Gather count and percent columns separately
Dermo_Inhibitor_hemo_2020_VIA_percent <- Dermo_Inhibitor_hemo_2020_VIA_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:13))
Dermo_Inhibitor_hemo_2020_VIA_counts <-  Dermo_Inhibitor_hemo_2020_VIA_counts %>%  group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:13))

Dermo_Inhibitor_hemo_2020_VIA_percent <- Dermo_Inhibitor_hemo_2020_VIA_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Dermo_Inhibitor_hemo_2020_VIA_counts <-  Dermo_Inhibitor_hemo_2020_VIA_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

Dermo_Inhibitor_hemo_2020_VIA_percent <- Dermo_Inhibitor_hemo_2020_VIA_percent[,-c(4,8)]
Dermo_Inhibitor_hemo_2020_VIA_counts <-  Dermo_Inhibitor_hemo_2020_VIA_counts[,-c(4,8)]

# Full join together so all columns preserved
Dermo_Inhibitor_hemo_2020_VIA_join <- full_join(Dermo_Inhibitor_hemo_2020_VIA_percent, Dermo_Inhibitor_hemo_2020_VIA_counts, by = c("ID","Treat","Assay","Plot_number","Gate", "Channel"))

# Add in cell type column based on plot number and gate
Dermo_Inhibitor_hemo_2020_VIA_cell_type <- data.frame(Plot_number = c("8","8", "19","19", "20","20", "21","21", "22","22"), Gate= c("P4","P3","V9-L","V9-R", "V12-L","V12-R", "V1-L", "V1-R","V4-L","V4-R"), Cell_type=c(
  "granular","agranular", "agranular","live_agranular", "granular","live_granular", 
  "agranular","dead_agranular","granular","dead_granular"))
# Join cell type
Dermo_Inhibitor_hemo_2020_VIA_join <- left_join(Dermo_Inhibitor_hemo_2020_VIA_join, Dermo_Inhibitor_hemo_2020_VIA_cell_type, by=c("Plot_number","Gate"))
unique_Dermo_Inhibitor_hemo_2020_VIA_join_gate <- unique(Dermo_Inhibitor_hemo_2020_VIA_join[,c("Plot_number","Gate")])

## Load the Viability data from the parasite alone 
# Make new header column
Dermo_Inhibitor_perk_2020_VIA_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/FORMATTED_CSVs/2020_Dermo_inhibitors_PERK_VIA_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Dermo_Inhibitor_perk_2020_VIA <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/FORMATTED_CSVs/2020_Dermo_inhibitors_PERK_VIA_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(Dermo_Inhibitor_perk_2020_VIA_nms)

# Split column 1 by space and remove
Dermo_Inhibitor_perk_2020_VIA <- Dermo_Inhibitor_perk_2020_VIA %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
Dermo_Inhibitor_perk_2020_VIA <- Dermo_Inhibitor_perk_2020_VIA[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
Dermo_Inhibitor_perk_2020_VIA_percent <- Dermo_Inhibitor_perk_2020_VIA[,c(1:3,5,7,9,11,13)]
Dermo_Inhibitor_perk_2020_VIA_counts <- Dermo_Inhibitor_perk_2020_VIA[,c(1:4,6,8,10,12)]

# Gather count and percent columns separately
Dermo_Inhibitor_perk_2020_VIA_percent <- Dermo_Inhibitor_perk_2020_VIA_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:8))
Dermo_Inhibitor_perk_2020_VIA_counts <-  Dermo_Inhibitor_perk_2020_VIA_counts %>%  group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:8))

Dermo_Inhibitor_perk_2020_VIA_percent <- Dermo_Inhibitor_perk_2020_VIA_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Dermo_Inhibitor_perk_2020_VIA_counts <-  Dermo_Inhibitor_perk_2020_VIA_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

Dermo_Inhibitor_perk_2020_VIA_percent <- Dermo_Inhibitor_perk_2020_VIA_percent[,-c(4,8)]
Dermo_Inhibitor_perk_2020_VIA_counts <-  Dermo_Inhibitor_perk_2020_VIA_counts[,-c(4,8)]

# Full join together so all columns preserved
Dermo_Inhibitor_perk_2020_VIA_join <- full_join(Dermo_Inhibitor_perk_2020_VIA_percent, Dermo_Inhibitor_perk_2020_VIA_counts, by = c("ID","Treat","Assay","Plot_number","Gate", "Channel"))

# Add in cell type column based on plot number and gate
Dermo_Inhibitor_perk_2020_VIA_cell_type <- data.frame(Plot_number = c("3", "19","19", "21","21"), Gate= c("P1","V9-L","V9-R","V4-L","V4-R"), Cell_type=c(
  "all_perk", "perkinsus","live_perkinsus", "perkinsus","dead_perkinsus"))
# Join cell type
Dermo_Inhibitor_perk_2020_VIA_join <- left_join(Dermo_Inhibitor_perk_2020_VIA_join, Dermo_Inhibitor_perk_2020_VIA_cell_type, by=c("Plot_number","Gate"))
unique_Dermo_Inhibitor_perk_2020_VIA_join_gate <- unique(Dermo_Inhibitor_perk_2020_VIA_join[,c("Plot_number","Gate")])

## APOPTOSIS ASSAY
# Make new header column
Dermo_Inhibitor_2020_APOP_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/FORMATTED_CSVs/2020_Dermo_inhibitors_dermo_hemocyte_APOP_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Dermo_Inhibitor_2020_APOP <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/FORMATTED_CSVs/2020_Dermo_inhibitors_dermo_hemocyte_APOP_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(Dermo_Inhibitor_2020_APOP_nms)

# Split column 1 by space and remove
Dermo_Inhibitor_2020_APOP <- Dermo_Inhibitor_2020_APOP %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
Dermo_Inhibitor_2020_APOP <- Dermo_Inhibitor_2020_APOP[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
Dermo_Inhibitor_2020_APOP_percent <- Dermo_Inhibitor_2020_APOP[,c(1:3,5,7,9,11,13,15,17,19,21, 23,25,27)]
Dermo_Inhibitor_2020_APOP_counts <-  Dermo_Inhibitor_2020_APOP[,c(1:4,6,8,10,12,14,16,18,20,22,24,26)]

# Gather count and percent columns separately
Dermo_Inhibitor_2020_APOP_percent <- Dermo_Inhibitor_2020_APOP_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:15))
Dermo_Inhibitor_2020_APOP_counts <-  Dermo_Inhibitor_2020_APOP_counts %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:15))

Dermo_Inhibitor_2020_APOP_percent <- Dermo_Inhibitor_2020_APOP_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Dermo_Inhibitor_2020_APOP_counts <-  Dermo_Inhibitor_2020_APOP_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

Dermo_Inhibitor_2020_APOP_percent <-Dermo_Inhibitor_2020_APOP_percent[,-c(4,8)]
Dermo_Inhibitor_2020_APOP_counts <- Dermo_Inhibitor_2020_APOP_counts[,-c(4,8)]

# Full join together so all columns preserved
Dermo_Inhibitor_2020_APOP_join <- full_join(Dermo_Inhibitor_2020_APOP_percent, Dermo_Inhibitor_2020_APOP_counts, by = c("ID","Treat","Assay","Plot_number","Gate","Channel"))

# Add in cell type column based on plot number and gate
Dermo_Inhibitor_2020_APOP_join_cell_type <- data.frame(Plot_number = c("2","2","8","8","4","4","4","4","7","7","7","7","9","9","9","9"), 
                                            Gate= c("V2-L","V2-R","P3","P4","Q16-UL","Q16-UR","Q16-LL","Q16-LR","Q15-UL","Q15-UR","Q15-LL","Q15-LR","Q10-UL","Q10-UR","Q10-LL","Q10-LR"), 
                                            Cell_type=c("non_apop_cells", "apop_cells", "agranular", "granular", "apop_granular_no_perk", "apop_granular_perk", "unstained_granular", "parasite_granular", 
                                                        "apop_agranular_no_perk", "apop_agranular_perk", "unstained_agranular", "parasite_agranular","apop_all_hemo_no_perk", "apop_all_hemo_perk", "unstained_all_hemo", "parasite"))

# Join cell type
Dermo_Inhibitor_2020_APOP_join <- left_join(Dermo_Inhibitor_2020_APOP_join, Dermo_Inhibitor_2020_APOP_join_cell_type, by=c("Plot_number","Gate"))
Dermo_Inhibitor_2020_APOP_join_unique <- unique(Dermo_Inhibitor_2020_APOP_join[,c("Plot_number","Gate")])

## CASPASE ASSAY
# Make new header column
Dermo_Inhibitor_2020_CASP_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/FORMATTED_CSVs/2020_Dermo_inhibitors_dermo_hemocyte_CASP_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Dermo_Inhibitor_2020_CASP <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/FORMATTED_CSVs/2020_Dermo_inhibitors_dermo_hemocyte_CASP_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(Dermo_Inhibitor_2020_CASP_nms)

# Split column 1 by space and remove
Dermo_Inhibitor_2020_CASP <- Dermo_Inhibitor_2020_CASP %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
Dermo_Inhibitor_2020_CASP <- Dermo_Inhibitor_2020_CASP[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
Dermo_Inhibitor_2020_CASP_percent <- Dermo_Inhibitor_2020_CASP[,c(1:3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35)]
Dermo_Inhibitor_2020_CASP_counts <-  Dermo_Inhibitor_2020_CASP[,c(1:4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34)]

# Gather count and percent columns separately
Dermo_Inhibitor_2020_CASP_percent <- Dermo_Inhibitor_2020_CASP_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:19))
Dermo_Inhibitor_2020_CASP_counts <-  Dermo_Inhibitor_2020_CASP_counts %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:19))

Dermo_Inhibitor_2020_CASP_percent <- Dermo_Inhibitor_2020_CASP_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Dermo_Inhibitor_2020_CASP_counts <-  Dermo_Inhibitor_2020_CASP_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

Dermo_Inhibitor_2020_CASP_percent <- Dermo_Inhibitor_2020_CASP_percent[,-c(4,8)]
Dermo_Inhibitor_2020_CASP_counts <-  Dermo_Inhibitor_2020_CASP_counts[,-c(4,8)]

# Full join together so all columns preserved
Dermo_Inhibitor_2020_CASP_join <- full_join(Dermo_Inhibitor_2020_CASP_percent, Dermo_Inhibitor_2020_CASP_counts, by = c("ID","Treat","Assay","Plot_number","Gate","Channel"))

# Remove data for plots 3, 6, and 9 which are not needed
Dermo_Inhibitor_2020_CASP_join <- Dermo_Inhibitor_2020_CASP_join %>% filter(Plot_number != "3" & Plot_number != "6" & Plot_number != "9")

# Add in cell type column based on plot number and gate
Dermo_Inhibitor_2020_CASP_cell_type <- data.frame(Plot_number = c("8","8","4","4","4","4","7","7","7","7","2","2","2","2"), 
                                            Gate= c("P3","P4","Q11-UL","Q11-UR", "Q11-LL", "Q11-LR", "Q10-UL","Q10-UR","Q10-LL","Q10-LR","Q5-UL","Q5-UR","Q5-LL","Q5-LR"), 
                                            Cell_type=c("agranular", 
                                                        "granular",  "caspase_active_agranular", "caspase_active_agranular_perk", "unstained_agranular", "caspase_active_agranular_perk",
                                                        "caspase_active_granular", "caspase_active_granular_perk", "unstained_granular", "caspase_active_granular_perk",
                                                        "caspase_active_all_hemo", "caspase_active_all_perk", "unstained_all_hemo", "caspase_active_all_perk"))

# Join cell type
Dermo_Inhibitor_2020_CASP_join <- left_join(Dermo_Inhibitor_2020_CASP_join, Dermo_Inhibitor_2020_CASP_cell_type, by=c("Plot_number","Gate"))
unique_Dermo_Inhibitor_2020_CASP_join <- unique(Dermo_Inhibitor_2020_CASP_join[,c("Plot_number","Gate")])

## JC1 ASSAY
# Make new header column
Dermo_Inhibitor_2020_JC1_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/FORMATTED_CSVs/2020_Dermo_inhibitors_dermo_hemocyte_JC1_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Dermo_Inhibitor_2020_JC1 <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/FORMATTED_CSVs/2020_Dermo_inhibitors_dermo_hemocyte_JC1_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(Dermo_Inhibitor_2020_JC1_nms)

# Split column 1 by space and remove
Dermo_Inhibitor_2020_JC1 <- Dermo_Inhibitor_2020_JC1 %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
Dermo_Inhibitor_2020_JC1 <- Dermo_Inhibitor_2020_JC1[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
Dermo_Inhibitor_2020_JC1_percent <- Dermo_Inhibitor_2020_JC1[,c(1:3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43)]
Dermo_Inhibitor_2020_JC1_counts <-  Dermo_Inhibitor_2020_JC1[,c(1:4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42)]

# Gather count and percent columns separately
Dermo_Inhibitor_2020_JC1_percent <- Dermo_Inhibitor_2020_JC1_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:23))
Dermo_Inhibitor_2020_JC1_counts <-  Dermo_Inhibitor_2020_JC1_counts %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:23))

Dermo_Inhibitor_2020_JC1_percent <- Dermo_Inhibitor_2020_JC1_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Dermo_Inhibitor_2020_JC1_counts <-  Dermo_Inhibitor_2020_JC1_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

Dermo_Inhibitor_2020_JC1_percent <- Dermo_Inhibitor_2020_JC1_percent[,-c(4,8)]
Dermo_Inhibitor_2020_JC1_counts <-  Dermo_Inhibitor_2020_JC1_counts[,-c(4,8)]

# Full join together so all columns preserved
Dermo_Inhibitor_2020_JC1_join <- full_join(Dermo_Inhibitor_2020_JC1_percent, Dermo_Inhibitor_2020_JC1_counts, by = c("ID","Treat","Assay","Plot_number","Gate", "Channel"))

# remove rows I don't need: Q11-LL, Q11-LR
Dermo_Inhibitor_2020_JC1_join <- Dermo_Inhibitor_2020_JC1_join %>% filter(Gate !="Q11-LL" & Gate !="Q11-LR")

# Add in cell type column based on plot number and gate
Dermo_Inhibitor_2020_JC1_cell_type <- data.frame(Plot_number = c("5", "5","8","8", "14","14","14","14","16","16","17","17","17","17","18","18",
                                                         "19","19","20","20"),     
                                         Gate= c("V1-L","V1-R", "P3","P4","H11-1","H11-2","H11-3","H11-4","Q11-UL","Q11-UR",
                                                 "H10-1","H10-2","H10-3","H10-4","Q28-UL","Q28-UR","Q6-UL","Q6-UR",
                                                 "Q7-UL","Q7-UR"),
                                         Cell_type=c("dead_non_stained_perkinsus","live_perkinsus", "agranular","granular", 
                                                     "normal_agranular", "normal_agranular", "apoptotic_agranular","unstained_agranular", 
                                                     "apoptotic_agranular_no_parasite", "apoptotic_agranular_parasite",
                                                    "normal_granular", "normal_granular", "apoptotic_granular","unstained_granular", 
                                                    "apoptotic_granular_no_parasite", "apoptotic_granular_parasite",
                                                    "normal_agranular_no_parasite", "normal_agranular_parasite",
                                                    "normal_agranular_no_parasite", "normal_agranular_parasite"))
# Join cell type
Dermo_Inhibitor_2020_JC1_join <- left_join(Dermo_Inhibitor_2020_JC1_join, Dermo_Inhibitor_2020_JC1_cell_type, by=c("Plot_number","Gate"))

## PHAGOCYTOSIS ASSAY
# Make new header column
Dermo_Inhibitor_2020_PHAGO_nms <-                                                                      
  read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/FORMATTED_CSVs/2020_Dermo_inhibitors_dermo_hemocyte_PHAGOCYTOSIS_formatted.xlsx", range = cell_rows(1:3), col_names = F) %>%       
  summarise_all(funs(paste(na.omit(.), collapse = "_"))) %>%                  
  unlist()

# Set the new column names
Dermo_Inhibitor_2020_PHAGO <- read_excel("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/FORMATTED_CSVs/2020_Dermo_inhibitors_dermo_hemocyte_PHAGOCYTOSIS_formatted.xlsx", skip = 2) %>%                                    
  magrittr::set_colnames(Dermo_Inhibitor_2020_PHAGO_nms)

# Split column 1 by space and remove
Dermo_Inhibitor_2020_PHAGO <- Dermo_Inhibitor_2020_PHAGO %>% 
  separate(ID_Treat_Assay, sep=" ", into = c("remove","ID")) 

# Separate new column 1 by dash, remove spaces from column names
Dermo_Inhibitor_2020_PHAGO <- Dermo_Inhibitor_2020_PHAGO[,-1] %>% separate(ID, sep="-", into=c("ID","Treat","Assay"))
Dermo_Inhibitor_2020_PHAGO_percent <- Dermo_Inhibitor_2020_PHAGO[,c(1:3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37)]
Dermo_Inhibitor_2020_PHAGO_counts <-  Dermo_Inhibitor_2020_PHAGO[,c(1:4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36)]

# Gather count and percent columns separately
Dermo_Inhibitor_2020_PHAGO_percent <- Dermo_Inhibitor_2020_PHAGO_percent %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_percent", value="Percent_of_this_plot", c(4:20))
Dermo_Inhibitor_2020_PHAGO_counts <-  Dermo_Inhibitor_2020_PHAGO_counts %>% group_by(ID,Assay,Treat) %>%  gather(key = "Plot_name_counts", value="Counts", c(4:20))

Dermo_Inhibitor_2020_PHAGO_percent <- Dermo_Inhibitor_2020_PHAGO_percent %>% separate(Plot_name_percent, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","%"), sep="_")
Dermo_Inhibitor_2020_PHAGO_counts <-  Dermo_Inhibitor_2020_PHAGO_counts %>% separate(Plot_name_counts, into = c("Plot_name","Plot_number", "Gate"), sep=' ') %>% separate(Gate, into=c("Channel","Gate","Count"), sep="_")

Dermo_Inhibitor_2020_PHAGO_percent <- Dermo_Inhibitor_2020_PHAGO_percent[,-c(4,8)]
Dermo_Inhibitor_2020_PHAGO_counts <-  Dermo_Inhibitor_2020_PHAGO_counts[,-c(4,8)]

# Full join together so all columns preserved
Dermo_Inhibitor_2020_PHAGO_join <- full_join(Dermo_Inhibitor_2020_PHAGO_percent, Dermo_Inhibitor_2020_PHAGO_counts, by = c("ID","Treat","Assay","Plot_number","Gate","Channel"))

# remove plots 6 and 9
Dermo_Inhibitor_2020_PHAGO_join <- Dermo_Inhibitor_2020_PHAGO_join %>% filter(Plot_number !="6" & Plot_number !="9")

# Add in cell type column based on plot number and gate
Dermo_Inhibitor_2020_PHAGO_cell_type <- data.frame(Plot_number = c("10","10","8","8","4","4","4","4","7","7","7","7","2","2","2","2"), 
                                                  Gate= c("V2-L","V2-R","P3","P4","Q2-UL","Q2-UR", "Q2-LL", "Q2-LR", "Q1-UL","Q1-UR","Q1-LL","Q1-LR","Q5-UL","Q5-UR","Q5-LL","Q5-LR"), 
                                                  Cell_type=c("unstained_parasite","stained_parasite","agranular",  "granular",  
                                                              "live_agranular", "agranular_phagocytosis", "unstained_agranular", "agranular_parasite_bead_alone",
                                                              "live_granular", "granular_phagocytosis", "unstained_agranular", "granular_parasite_bead_alone",
                                                              "live_hemo", "hemo_phagocytosis", "unstained_hemo", "hemo_parasite_bead_alone"))

# Join cell type
Dermo_Inhibitor_2020_PHAGO_join <- left_join(Dermo_Inhibitor_2020_PHAGO_join, Dermo_Inhibitor_2020_PHAGO_cell_type, by=c("Plot_number","Gate"))
unique_Dermo_Inhibitor_2020_PHAGO_join <- unique(Dermo_Inhibitor_2020_PHAGO_join[,c("Plot_number","Gate")])

#### 2020 Dermo and Inhibitors Arcsine Transformation of percentages ####

# Make list of dataframes
Dermo_Inhibitor_2020_list <- list(Dermo_Inhibitor_hemo_2020_VIA_join=Dermo_Inhibitor_hemo_2020_VIA_join, 
                                  Dermo_Inhibitor_perk_2020_VIA_join=Dermo_Inhibitor_perk_2020_VIA_join, 
                                  Dermo_Inhibitor_2020_APOP_join=Dermo_Inhibitor_2020_APOP_join, 
                                  Dermo_Inhibitor_2020_CASP_join=Dermo_Inhibitor_2020_CASP_join,
                                  Dermo_Inhibitor_2020_JC1_join=Dermo_Inhibitor_2020_JC1_join,
                                  Dermo_Inhibitor_2020_PHAGO_join=Dermo_Inhibitor_2020_PHAGO_join)

Dermo_Inhibitor_hemo_2020_VIA_join
Dermo_Inhibitor_perk_2020_VIA_join
Dermo_Inhibitor_2020_APOP_join
Dermo_Inhibitor_2020_CASP_join
Dermo_Inhibitor_2020_JC1_join
Dermo_Inhibitor_2020_PHAGO_join

# Make new column and perform arcsine
Dermo_Inhibitor_hemo_2020_VIA_join$Percent_of_this_plot_arcsine <-  transf.arcsin(Dermo_Inhibitor_hemo_2020_VIA_join$Percent_of_this_plot)
Dermo_Inhibitor_perk_2020_VIA_join$Percent_of_this_plot_arcsine <-  transf.arcsin(Dermo_Inhibitor_perk_2020_VIA_join$Percent_of_this_plot)
Dermo_Inhibitor_2020_APOP_join$Percent_of_this_plot_arcsine <- transf.arcsin(Dermo_Inhibitor_2020_APOP_join$Percent_of_this_plot)
Dermo_Inhibitor_2020_CASP_join$Percent_of_this_plot_arcsine <- transf.arcsin(Dermo_Inhibitor_2020_CASP_join$Percent_of_this_plot)
Dermo_Inhibitor_2020_JC1_join$Percent_of_this_plot_arcsine <- transf.arcsin(Dermo_Inhibitor_2020_JC1_join$Percent_of_this_plot)
Dermo_Inhibitor_2020_PHAGO_join$Percent_of_this_plot_arcsine <- transf.arcsin(Dermo_Inhibitor_2020_PHAGO_join$Percent_of_this_plot)

# Change percent column to percent
Dermo_Inhibitor_hemo_2020_VIA_join$Percent_of_this_plot <- Dermo_Inhibitor_hemo_2020_VIA_join$Percent_of_this_plot*100
Dermo_Inhibitor_perk_2020_VIA_join$Percent_of_this_plot <- Dermo_Inhibitor_perk_2020_VIA_join$Percent_of_this_plot*100
Dermo_Inhibitor_2020_APOP_join$Percent_of_this_plot<- Dermo_Inhibitor_2020_APOP_join$Percent_of_this_plot*100
Dermo_Inhibitor_2020_CASP_join$Percent_of_this_plot<- Dermo_Inhibitor_2020_CASP_join$Percent_of_this_plot*100
Dermo_Inhibitor_2020_JC1_join$Percent_of_this_plot <-  Dermo_Inhibitor_2020_JC1_join$Percent_of_this_plot*100
Dermo_Inhibitor_2020_PHAGO_join$Percent_of_this_plot <-  Dermo_Inhibitor_2020_PHAGO_join$Percent_of_this_plot*100

### 2020 Dermo Inhibitors VIABILITY ASSAY Statistics and Plotting ####

## Analyze Viability from the hemocyte data  
Dermo_Inhibitor_hemo_2020_VIA_join

# Agranular and Granular cells Boxplot with significance bars, grouped by treatment color by gate(significant based on ANOVA)
Dermo_Inhibitor_hemo_2020_VIA_join_Percent_Agranular_Granular <- Dermo_Inhibitor_hemo_2020_VIA_join %>% filter(Gate =="P3" | Gate=="P4")

Dermo_Inhibitor_hemo_2020_VIA_join_Percent_Agranular_Granular_plot <- ggplot(data=Dermo_Inhibitor_hemo_2020_VIA_join_Percent_Agranular_Granular,
                                                                  aes(y=Percent_of_this_plot, x=ID, fill=Gate)) + geom_col()  + 
  xlab("Individual") +
  ylab("Percent Hemocytes") +
  ggtitle("Percent of Granular and Agranular") + 
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, hjust = 1, angle = 90)) +
  theme(legend.text = element_text(size=12)) + 
  #  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  scale_fill_manual(name="Cell Type", labels=c("Agranular","Granular"), values=c("#6c81d9","#50b47b")) 

#save
ggsave(plot = Dermo_Inhibitor_hemo_2020_VIA_join_Percent_Agranular_Granular_plot, device = "tiff", filename = "Dermo_Inhibitor_hemo_2020_VIA_join_Percent_Agranular_Granular_plot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 6)

# Granular Agranular plot formatted for multipanel figure 
# calculate sd and mean 
Dermo_Inhibitor_hemo_2020_VIA_join_Percent_Agranular_Granular_sd <- Dermo_Inhibitor_hemo_2020_VIA_join_Percent_Agranular_Granular %>%
  # get mean and sd
  dplyr::group_by(Gate) %>% mutate(mean = mean(Percent_of_this_plot), sd = sd(Percent_of_this_plot)) 

Dermo_Inhibitor_hemo_2020_VIA_join_Percent_Agranular_Granular_plot_multipanel <- 
  ggplot(data=Dermo_Inhibitor_hemo_2020_VIA_join_Percent_Agranular_Granular_sd,
  aes(y=Percent_of_this_plot, x=Gate)) + geom_bar(aes(fill= Gate),position = "dodge", stat = "summary")  + 
  geom_point(shape = 1, aes(fill = Gate)) + 
  xlab(NULL) +
  ylab("% Hemocytes") + 
  theme_classic() +
  theme(text=element_text(size=12, face = "bold"), 
                 axis.title.y=element_text(size=12, face = "bold"),
                 axis.title.x=element_text(size=12,face = "bold"),
        axis.text.x=element_text(size=12,face = "bold"),
        axis.text.y=element_text(size=12,face = "bold")) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  scale_x_discrete(labels = c("P3" = "Agranular", "P4"= "Granular")) + 
  scale_fill_manual(name="Cell Type", labels=c("Agranular","Granular"), values=c("#6c81d9","#50b47b")) 

# Perform t.test and plot results onto boxplot 
stat.test <- as.data.frame(Dermo_Inhibitor_hemo_2020_VIA_join_Percent_Agranular_Granular_sd) %>%
  # make sure to use the arcsine values 
  t_test(Percent_of_this_plot_arcsine ~ Gate) %>%
  add_significance(p.col = "p") %>% 
  add_xy_position(x = "Gate")

Dermo_Inhibitor_hemo_2020_VIA_join_Percent_Agranular_Granular_plot_multipanel_sig <- 
  Dermo_Inhibitor_hemo_2020_VIA_join_Percent_Agranular_Granular_plot_multipanel + stat_pvalue_manual(
  stat.test, label = "{p} {p.signif}",  tip.length = 0.02, y.position = 75)

# Average granular and agranular cells
Dermo_Inhibitor_hemo_2020_VIA_join_Percent_Agranular_Granular %>% dplyr::group_by(Gate) %>% summarise_at(vars(Percent_of_this_plot), list(name = mean))
  #Gate   name
  #<chr> <dbl>
  #1 P3     59.8
  #2 P4     40.4

# Compare the percent of dead cells
Dermo_Inhibitor_2020_VIA_join_Percent_Agranular_Granular_DEAD_plot <- Dermo_Inhibitor_hemo_2020_VIA_join %>%
  # filter for granular and agranular dead hemocytes
  filter(Gate=="V1-R" |Gate=="V4-R") %>%
  ggplot(aes(y=Percent_of_this_plot, x=Treat, fill=Gate))  + geom_boxplot() +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Dead Granular and Agranular Hemocytes") + 
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, hjust = 1, angle = 90)) +
  theme(legend.text = element_text(size=12)) + 
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,20)) +
  scale_fill_manual(name="Cell Type", labels=c("Percent Dead Agranular","Percent Dead Granular"), values=c("#5b2c90", "#88bf3b")) 

#save
ggsave(plot = Dermo_Inhibitor_2020_VIA_join_Percent_Agranular_Granular_DEAD_plot, device = "tiff", filename = "Dermo_Inhibitor_2020_VIA_join_Percent_Agranular_Granular_DEAD_plot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 6)

# Plot Dead hemocytes as live and prepare for multi-panel figure 
Dermo_Inhibitor_hemo_2020_VIA_join_LIVE_sd <- Dermo_Inhibitor_hemo_2020_VIA_join %>%
  # filter for granular and agranular dead hemocytes
  filter(Gate=="V1-R" |Gate=="V4-R") %>%
  mutate(Percent_of_this_plot_live = (100 - Percent_of_this_plot)) %>%
  # get mean and sd
  dplyr::group_by(Gate) %>% mutate(mean = mean(Percent_of_this_plot_live), sd = sd(Percent_of_this_plot_live)) 

# add arcine transformed percentages 
Dermo_Inhibitor_hemo_2020_VIA_join_LIVE_sd$Percent_of_this_plot_live_arcsine <- transf.arcsin(Dermo_Inhibitor_hemo_2020_VIA_join_LIVE_sd$Percent_of_this_plot_live*0.01)

# Make plot of live hemocytes 
Dermo_Inhibitor_hemo_2020_VIA_join_LIVE_sd_multipanel <- 
  ggplot(data=Dermo_Inhibitor_hemo_2020_VIA_join_LIVE_sd,
         aes(y=Percent_of_this_plot_live, x=Gate)) + geom_bar(aes(fill= Gate),position = "dodge", stat = "summary")  + 
  geom_point(shape = 1, aes(fill = Gate)) + 
  xlab(NULL) +
  ylab("% Live Hemocytes") + 
  theme_classic() +
  theme(text=element_text(size=12, face = "bold"), 
        axis.title.y=element_text(size=12, face = "bold"),
        axis.title.x=element_text(size=12,face = "bold"),
        axis.text.x=element_text(size=12,face = "bold"),
        axis.text.y=element_text(size=12,face = "bold")) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  scale_x_discrete(labels = c("V1-R" = "Agranular", "V4-R"= "Granular")) + 
  scale_fill_manual(name="Cell Type", labels=c("Agranular","Granular"), values=c("#6c81d9","#50b47b")) 

# Perform t.test and plot results onto boxplot 
stat.test <- as.data.frame(Dermo_Inhibitor_hemo_2020_VIA_join_LIVE_sd) %>%
  # make sure to use the arcsine values 
  t_test(Percent_of_this_plot_live_arcsine ~ Gate) %>%
  add_significance(p.col = "p") %>% 
  add_xy_position(x = "Gate")

Dermo_Inhibitor_hemo_2020_VIA_join_LIVE_sd_multipanel_sig <- 
  Dermo_Inhibitor_hemo_2020_VIA_join_LIVE_sd_multipanel + stat_pvalue_manual(
    stat.test, label = "{p} {p.signif}",  tip.length = 0.02, y.position = 100)

## Analysis of dead hemocytes
Dermo_Inhibitor_hemo_2020_VIA_join %>%
  # filter for granular and agranular dead hemocytes
  filter(Gate=="V1-R" |Gate=="V4-R") %>% dplyr::group_by(Gate) %>% summarise_at(vars(Percent_of_this_plot), list(name = mean))

#Gate   name
#<chr> <dbl>
#1 V1-R   1.33
#2 V4-R   7.94

### Viability Analysis for Parasite

Dermo_Inhibitor_2020_VIA_join_Percent_PERK_DEAD_plot <- Dermo_Inhibitor_perk_2020_VIA_join %>%
  # filter for granular and agranular dead hemocytes
  filter(Gate=="V4-R") %>%
  ggplot(aes(y=Percent_of_this_plot, x=Treat, fill=Gate))  + geom_boxplot() +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Dead Perkinsus") + 
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, hjust = 1, angle = 90)) +
  theme(legend.text = element_text(size=12)) + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,20)) +
  scale_fill_manual(name="Cell Type", labels=c("Percent Dead Perkinsus"), values=c("#5b2c90")) 

#save
ggsave(plot = Dermo_Inhibitor_2020_VIA_join_Percent_PERK_DEAD_plot, device = "tiff", filename = "Dermo_Inhibitor_2020_VIA_join_Percent_PERK_DEAD_plot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 6)

## Reformat plot to match others for multi-panel
Dermo_Inhibitor_perk_2020_VIA_join_LIVE_multipanel <- 
  Dermo_Inhibitor_perk_2020_VIA_join %>%
  # filter for granular dead P. marinus
  filter(Gate=="V4-R") %>% 
  # plot LIVE
  mutate(Percent_of_this_plot_live = 100 - Percent_of_this_plot) %>%
  dplyr::group_by(Gate) %>% mutate(mean = mean(Percent_of_this_plot_live), sd = sd(Percent_of_this_plot_live)) %>%
  ggplot(data=.,
         aes(y=Percent_of_this_plot_live, x=Gate)) + geom_bar(aes(fill= Gate),position = "dodge", stat = "summary")  + 
  geom_point(shape = 1, aes(fill = Gate)) +
  labs(x = NULL , y ="% Live P. marinus") + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 12, face= "bold"),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  scale_x_discrete(labels = c("V4-R" = "*P. marinus*")) + 
  scale_fill_manual(name="Cell Type", labels=c("*P. marinus*"), values=c("#5b2c90")) 

Dermo_Inhibitor_perk_2020_VIA_join_LIVE_multipanel_sig <- 
  Dermo_Inhibitor_perk_2020_VIA_join_LIVE_multipanel + 
  theme(axis.text.x=ggtext::element_markdown(),
   legend.text = ggtext::element_markdown()) 

  # Analysis of dead P. marinus
Dermo_Inhibitor_perk_2020_VIA_join %>%
  # filter for granular and agranular dead hemocytes
  filter(Gate=="V4-R") %>% dplyr::group_by(Gate) %>% summarise_at(vars(Percent_of_this_plot), list(name = mean))

#Gate   name
#<chr> <dbl>
#  1 V4-R   4.98

#### 2020 Dermo and Inhibitors PHAGOCYTOSIS ASSAY ####

Dermo_Inhibitor_2020_PHAGO_join
# REMINDER: LPS was not actually added to this treatment (I had accidentally left this step out of the protocol)

# Plot the percent of hemocytes that have phagocytosed a bead 
# to do this I need to get rid of the quandrant that has only beads and no hemocytes so I can get a more accurate percentage

# first perform this for the agranular cells
Dermo_Inhibitor_2020_PHAGO_join_agranular_phago <- Dermo_Inhibitor_2020_PHAGO_join %>% filter(Gate == "Q2-UL" | Gate == "Q2-UR" | Gate == "Q2-LL") %>% ungroup() %>%
  dplyr::group_by(ID, Treat, Plot_number, Assay) %>%
  dplyr::mutate(agranular_sum_no_beads_perk_alone = sum(Counts)) %>%
  dplyr::mutate(revised_percent_of_this_plot = Counts/agranular_sum_no_beads_perk_alone) %>%
  mutate(revised_percent_of_this_plot_arcsine = revised_percent_of_this_plot) %>%
  mutate(revised_percent_of_this_plot = revised_percent_of_this_plot*100)
Dermo_Inhibitor_2020_PHAGO_join_agranular_phago$revised_percent_of_this_plot_arcsine <- transf.arcsin(Dermo_Inhibitor_2020_PHAGO_join_agranular_phago$revised_percent_of_this_plot_arcsine)

# repeat for the granular cells
Dermo_Inhibitor_2020_PHAGO_join_granular_phago <- Dermo_Inhibitor_2020_PHAGO_join %>% filter(Gate == "Q1-UL" | Gate == "Q1-UR" | Gate == "Q1-LL") %>% ungroup() %>%
  dplyr::group_by(ID, Treat, Plot_number, Assay) %>%
  dplyr::mutate(granular_sum_no_beads_perk_alone = sum(Counts)) %>%
  dplyr::mutate(revised_percent_of_this_plot = Counts/granular_sum_no_beads_perk_alone) %>%
  mutate(revised_percent_of_this_plot_arcsine = revised_percent_of_this_plot) %>%
  mutate(revised_percent_of_this_plot = revised_percent_of_this_plot*100)
Dermo_Inhibitor_2020_PHAGO_join_granular_phago$revised_percent_of_this_plot_arcsine <- transf.arcsin(Dermo_Inhibitor_2020_PHAGO_join_granular_phago$revised_percent_of_this_plot_arcsine)

# Combine together so I can plot the revised percentage of phagocytosis
Dermo_Inhibitor_2020_PHAGO_join_phago_combined <- 
  rbind(Dermo_Inhibitor_2020_PHAGO_join_agranular_phago[,c("ID","Treat","Assay","Gate","revised_percent_of_this_plot", "revised_percent_of_this_plot_arcsine")],
  Dermo_Inhibitor_2020_PHAGO_join_granular_phago[,c("ID","Treat","Assay", "Gate","revised_percent_of_this_plot", "revised_percent_of_this_plot_arcsine")])

# Plot the Percent of cells in each quadrant
Dermo_Inhibitor_2020_PHAGO_join_phago_combined_plot <- ggplot(Dermo_Inhibitor_2020_PHAGO_join_phago_combined, aes(x=Gate, y = revised_percent_of_this_plot, fill = Treat)) + 
         geom_boxplot() +
           ylab("Percent Hemocytes") + 
           ggtitle("Percent Phagocytic Agranular and Granular Hemocytes") + 
           theme(panel.background=element_blank(),
                 panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
                 text=element_text(family="serif",size=12), 
                 axis.title.y=element_text(family="serif",size=12),
                 axis.title.x=element_text(family="serif",size=12),
                 legend.key=element_rect(fill=NA)) + 
           theme(text=element_text(size=12)) + 
           theme(legend.text = element_text(size=12)) + 
           scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
           scale_fill_manual(name="Treatment", labels=c("Ruby Red Microspheres", "Perkinsus marinus"), values=c("#aa4dce","#c89832")) +
           scale_x_discrete(breaks= c("Q1-LL","Q1-UL","Q1-UR","Q2-LL","Q2-UL","Q2-UR"),
                            labels = c("Unstained\nGranular","Live\nGranular","Live\nGranular\nPhagocytosis","Unstained\nAgranular","Live\nAgranular","Live\nAgranular\nPhagocytosis"))

#"#cc57b4", "#88bf3b", "#aa4dce", "#56b464", "#5b2c90"
#"#c89832", "#5a6ee6", "#ca4e33", "#7e78d4", "#cd4272"         
         
#save
ggsave(plot = Dermo_Inhibitor_2020_PHAGO_join_phago_combined_plot , device = "tiff", filename = "Dermo_Inhibitor_2020_PHAGO_join_phago_combined_plot.tiff",
      path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
      height = 5, width = 6.5)

# Plot just the Q1-UR granular phagocytic quadrant of interest for simplicity
Dermo_Inhibitor_2020_PHAGO_join_phago_granular_plot <- Dermo_Inhibitor_2020_PHAGO_join_phago_combined %>%
  filter(Gate == "Q1-UR") %>%
  ggplot(., aes(x=Gate, y = revised_percent_of_this_plot, fill = Treat)) + 
  geom_boxplot() +
  ylab("Percent Hemocytes") + 
  ggtitle("Percent Phagocytic Agranular and Granular Hemocytes") + 
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(legend.text = element_text(size=12)) + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,20), breaks = c(0,5,10,15,20)) +
  scale_fill_manual(name="Treatment", labels=c("Ruby Red Microspheres", "Perkinsus marinus"), values=c("#aa4dce","#c89832")) +
  scale_x_discrete(breaks= c("Q1-UR"),
                   labels = c("Live Granular Phagocytosis"))

#"#cc57b4", "#88bf3b", "#aa4dce", "#56b464", "#5b2c90"
#"#c89832", "#5a6ee6", "#ca4e33", "#7e78d4", "#cd4272"         

#save
ggsave(plot = Dermo_Inhibitor_2020_PHAGO_join_phago_granular_plot , device = "tiff", filename = "Dermo_Inhibitor_2020_PHAGO_join_phago_granular_plot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 6.5)

## Format plot for multi-panel figure 
Dermo_Inhibitor_2020_PHAGO_join_phago_combined_Q1_UR <- Dermo_Inhibitor_2020_PHAGO_join_phago_combined %>%
  filter(Gate == "Q1-UR") %>% 
  dplyr::group_by(Treat) %>% mutate(mean = mean(revised_percent_of_this_plot), sd = sd(revised_percent_of_this_plot)) 

# average percent phagocytosis of parasites in the granular range is 13.219

Dermo_Inhibitor_2020_PHAGO_join_phago_combined_Q1_UR_multipanel <- 
  Dermo_Inhibitor_2020_PHAGO_join_phago_combined_Q1_UR %>%
  ggplot(data=.,
         aes(y=revised_percent_of_this_plot, x=Treat)) + geom_bar(aes(fill= Treat),position = "dodge", stat = "summary")  + 
  geom_point(shape = 1, aes(fill = Treat)) +
  labs(x = NULL , y ="% Granular Phagocytosis") + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 12, face= "bold"),
       legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100)) +
  scale_x_discrete(labels = c("PERK_hemo" = "*P. marinus*", "BEADS_LPS"="Beads")) + 
  scale_fill_manual(name="Treatment", labels=c("Beads","*P. marinus*"), values=c("#88bf3b","#5b2c90")) 

Dermo_Inhibitor_2020_PHAGO_join_phago_combined_Q1_UR_multipanel_sig <- 
  Dermo_Inhibitor_2020_PHAGO_join_phago_combined_Q1_UR_multipanel + 
  theme(axis.text.x=ggtext::element_markdown(),
        legend.text = ggtext::element_markdown()) 
  
# Perform t.test and plot results onto boxplot 
stat.test <- as.data.frame(Dermo_Inhibitor_2020_PHAGO_join_phago_combined_Q1_UR) %>%
  # make sure to use the arcsine values 
  t_test(revised_percent_of_this_plot_arcsine ~ Treat) %>%
  add_significance(p.col = "p") %>% 
  add_xy_position(x = "Treat")

Dermo_Inhibitor_2020_PHAGO_join_phago_combined_Q1_UR_multipanel_sig <- 
  Dermo_Inhibitor_2020_PHAGO_join_phago_combined_Q1_UR_multipanel_sig + stat_pvalue_manual(
    stat.test, label = "{p} {p.signif}",  tip.length = 0.02, y.position = 25)

# Analysis of Phagocytic cells percentages
Dermo_Inhibitor_2020_PHAGO_join_phago_combined %>%
  filter(Gate == "Q1-UR")  %>% dplyr::group_by(Treat) %>% summarise_at(vars(revised_percent_of_this_plot), list(name = mean))

### Statistical analysis

# goal is to compare the percent phagocytosis of beads vs parasite in agranulocytes
Dermo_Inhibitor_2020_PHAGO_join_phago_combined_phago_only <- Dermo_Inhibitor_2020_PHAGO_join_phago_combined %>% filter(Gate == "Q1-UR")

summary(aov(revised_percent_of_this_plot_arcsine~Treat, data =Dermo_Inhibitor_2020_PHAGO_join_phago_combined_phago_only))
#Df  Sum Sq Mean Sq F value   Pr(>F)    
#Treat        1 0.07318 0.07318   106.8 0.000495 ***
#  Residuals    4 0.00274 0.00069                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Perkinsus is phagocytosed significantly more than the BEADS treatment

#### 2020 Dermo and Inhibitors Multi-panel VIA and PHAGO figure ####

# include two spaces to insert flow cytometry plots 
VIA_PHAGO_multipanel <- cowplot::plot_grid( Dermo_Inhibitor_hemo_2020_VIA_join_Percent_Agranular_Granular_plot_multipanel_sig, 
# insert white space for hemocyte gating 
          NULL, Dermo_Inhibitor_hemo_2020_VIA_join_LIVE_sd_multipanel_sig, 
Dermo_Inhibitor_perk_2020_VIA_join_LIVE_multipanel_sig,
Dermo_Inhibitor_2020_PHAGO_join_phago_combined_Q1_UR_multipanel_sig, NULL,
ncol = 2, nrow = 3, labels = c("A","B","C","D","E","F"),
label_size = 16, label_fontface = "bold", align = "hv")

ggsave(VIA_PHAGO_multipanel, device = "tiff", filename = "VIA_PHAGO_multipanel_plot.tiff",
        path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
        height = 10, width = 8 )

#### 2020 Dermo and Inhibitors APOPTOSIS ASSAY Statistics and Plotting ####

Dermo_Inhibitor_2020_APOP_join

## Start by comparing the number of parasite cells in Q16-LR in the stained dermo plots to the Q16-LL in the unstained Dermo plots
  # this will help assure my self that similar levels of non-apoptotic P. marinus are present in both
Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk <- Dermo_Inhibitor_2020_APOP_join %>% filter(Gate == "Q16-LL" & Treat == "PERK" ) %>% 
  mutate(sample_ID = paste(ID, Treat, sep = "-"))
Dermo_Inhibitor_2020_APOP_join_non_apop_granular_treat <- Dermo_Inhibitor_2020_APOP_join %>% filter(Gate == "Q16-LR") %>%
  filter(Treat == "Dermo" | Treat == "Dermo_GDC" | Treat == "Dermo_ZVAD") %>% 
  mutate(sample_ID = paste(ID, Treat, sep = "-"))
Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk_treat <- rbind(Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk, 
                                                                     Dermo_Inhibitor_2020_APOP_join_non_apop_granular_treat)
#Plot 
Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk_treat_plot <- 
  ggplot(Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk_treat, aes(x= sample_ID, y = Counts, fill = Gate)) + 
  geom_col(position = "dodge") +
  xlab("Treatment") +
  ylab("Cell Counts") + 
  ggtitle("Counts of Non-Apoptotic Perkinsus") + 
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, angle= 90, hjust=1.0)) +
  theme(legend.text = element_text(size=12)) +
  scale_fill_manual(name="Cell Type", labels=c("Parasite Control", "Treatment Groups"), 
                     values = c("#50b47b",
                                "#7f63b8")) 
#save
ggsave(plot = Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk_treat_plot, device = "tiff", filename = "Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk_treat_plot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 10)
# Remember, each pool got a different amount of Dermo cells..even the ones for the same pool had different amounts
# each of the Perk pools has less Perk though than their pools. Pretty standard amount of non-apoptotic perkinsus in each one though

## What is the counts ratio in the PERK only pool of apoptotic to non apoptotic cells
Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk_all <- Dermo_Inhibitor_2020_APOP_join %>% 
  filter(Gate == "Q16-LL" | Gate == "Q16-UL") %>% filter(Treat == "PERK" ) %>% 
  mutate(sample_ID = paste(ID, Treat, sep = "-"))
# calculate ratio
Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk_all %>% dplyr::group_by(ID) %>%
  summarise(apop_ratio = Counts[Gate == "Q16-UL"] / Counts[Gate == "Q16-LL"]) %>% ungroup() %>% summarise(mean = mean(apop_ratio))
      #apop_ratio: average is 0.39
  # apop ratio is between 35% and 46%..so about 40% of the perkinsus cells are apoptotic...so I can roughly calculate what the ratio would
    # be for the other assay

#Plot apoptotic perkinsus
Dermo_Inhibitor_2020_APOP_join_non_apop_granular_apoptotic_perkinsus_plot <- 
  ggplot(Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk_all, aes(x= sample_ID, y = Counts, fill = Gate)) + 
  geom_col(position = "dodge") +
  xlab("Treatment") +
  ylab("Cell Counts") + 
  ggtitle("Counts of Apoptotic Perkinsus") + 
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, angle= 90, hjust=1.0)) +
  theme(legend.text = element_text(size=12)) +
  scale_fill_manual(name="Cell Type", labels=c("Non-apoptotic P. marinus", "Apoptotic P. marinus"), 
                    values = c("#50b47b",
                               "#7f63b8")) 
#save
ggsave(plot = Dermo_Inhibitor_2020_APOP_join_non_apop_granular_apoptotic_perkinsus_plot, device = "tiff", filename = "Dermo_Inhibitor_2020_APOP_join_non_apop_granular_apoptotic_perkinsus_plot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 10)

## Plot percent of apoptotic perkinsus and format for multi-panel figure
Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk_all_sd <- Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk_all %>%
  group_by(Gate) %>% mutate(mean = mean(Percent_of_this_plot), sd = sd(Percent_of_this_plot))

Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk_all_multipanel <- 
ggplot(data=Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk_all_sd,
         aes(fill = factor(Gate,levels = c("Q16-UL","Q16-LL")), y=Percent_of_this_plot, x=Assay)) + geom_bar(position="fill", stat="identity")  + 
  labs(x = NULL , y ="% P. marinus") + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold"),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold"),
        legend.position = "bottom") +
  guides(fill=guide_legend(nrow=2, title.position = "top")) +
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = scales::percent) +
 scale_x_discrete(labels = c("APOP" = "*P. mar.*")) + 
  scale_fill_manual(name="Apoptosis<br>Assay", labels=c("Apoptotic", "Non-Apoptotic"),
                   values=c("#a44f9a","#5b2c90")) 

Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk_all_multipanel_sig <- 
  Dermo_Inhibitor_2020_APOP_join_non_apop_granular_perk_all_multipanel + 
  theme(axis.text.x=ggtext::element_markdown(),
        legend.text = ggtext::element_markdown(),
        legend.title = ggtext::element_markdown()) 

## Calculate rough counts of apoptotic parasite based on Q16-LR apoptotic parasite ratio from the perkinsus only assay
Dermo_Inhibitor_2020_APOP_join_non_apop_granular_treat_approx_perk_apop <- Dermo_Inhibitor_2020_APOP_join_non_apop_granular_treat %>% 
  mutate(apop_parasite_alone = Counts * 0.396) %>%
  mutate(Percent_of_this_plot = NA, Percent_of_this_plot_arcsine = NA) %>%
  # gather this column
  gather("Gate", "Counts", apop_parasite_alone)

# now add this gate and counts back in to the treatments table and subtract the apop_parasite_counts from the Q16-UR
Dermo_Inhibitor_2020_APOP_join_granular_approx_perk_apop_treat <- 
  # subset the Q16-UR for all the assays
  Dermo_Inhibitor_2020_APOP_join %>% filter(Gate == "Q16-UR") %>%
  filter(Treat == "Dermo" | Treat == "Dermo_GDC" | Treat == "Dermo_ZVAD") %>% 
  mutate(sample_ID = paste(ID, Treat, sep = "-")) %>%
  # join the calculated number of counts that need to be subtracted
  rbind(.,Dermo_Inhibitor_2020_APOP_join_non_apop_granular_treat_approx_perk_apop) %>%
  group_by(sample_ID) %>%
    # subtract the apop_parasite_alone counts
  summarise(counts_minus_apop_parasite = Counts[Gate == "Q16-UR"] - Counts[Gate == "apop_parasite_alone"])

##### this produces negative numbers! Not going to use this normalization technique
## Decided to remove the correction for P. marinus apoptosis since the apoptosis levels were similar in all pools, 
  ## and parasite engulfment was similar across all pools, so any correction for this would apply equally
  ## to all groups and not affect statistical significance 

### Plotting granular cell data

# change treatment labels, remember that we are not plotting here parasite alone or the UV
Dermo_Inhibitor_2020_APOP_join_granular_apoptotic$Treat <- factor(Dermo_Inhibitor_2020_APOP_join_granular_apoptotic$Treat, 
                                                                    levels = c( "BEADS_LPS","Control_hemo", "Dermo","Dermo_GDC","Dermo_ZVAD"),
                                                                    labels = c( "Beads and\n LPS", "Control hemocytes",    "P. mar to\n Hemocytes\n" , 
                                                                                "P. mar to\n Hemocytes\n GDC",  
                                                                                "P. mar to\n Hemocytes\n ZVAD"))
# Make plot of percentage of apoptotic hemocytes with parasite and without
Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_plot <- Dermo_Inhibitor_2020_APOP_join %>%
  # plot percentage of apoptotic hemocytes with and without parasite
  filter(Gate == "Q16-UR" | Gate == "Q16-UL") %>%
  ggplot(data=., aes(y=Percent_of_this_plot, x=Treat, color=Gate)) + geom_point(position=position_dodge(width=0.75)) + 
  geom_boxplot() +
  xlab("Treatment") +
  ylab("Percent of Granulocytes") + 
  ggtitle("Percent of Apoptotic Granulocytes") + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,40), breaks = c(0,5,10,20,30,40)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  scale_color_manual(name="Cell Type", labels=c("Apoptotic\nGranular", "Apoptotic Granular,\nwith P. mar., or Beads"), 
                     values = c("#56b464", "#5b2c90")) 

# color options
#"#cc57b4", "#88bf3b", "#aa4dce", "#56b464", "#5b2c90"
#"#c89832", "#5a6ee6", "#ca4e33", "#7e78d4", "#cd4272"

#save
ggsave(plot = Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_plot, device = "tiff", filename = "Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_plot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 10)


## Plot just the 16-UR data 

Dermo_Inhibitor_2020_APOP_join_granular_apoptotic <- Dermo_Inhibitor_2020_APOP_join %>% filter(Gate == "Q16-UR") %>% filter(Treat !="UV") %>% filter(Treat !="PERK")

##Plot apoptosis granulocytes in format for multipanel figure with multiple comparisons run 
Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_sd <-   Dermo_Inhibitor_2020_APOP_join_granular_apoptotic %>%
  filter(Gate == "Q16-UR") %>% group_by(Treat) %>% mutate(mean = mean(Percent_of_this_plot), sd = sd(Percent_of_this_plot))

Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_sd$Treat <- factor(Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_sd$Treat,
                                                                         levels = c("Control_hemo","BEADS_LPS","Dermo","Dermo_GDC","Dermo_ZVAD"))

Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_sd_multipanel <- 
ggplot(data=Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_sd,
       aes(y=Percent_of_this_plot, x=Treat)) + 
  geom_bar(aes(fill=Treat), position="dodge", stat = "summary", fill = "#6d8dd7")  + 
  geom_point(aes(x= Treat, shape = ID), size = 3) +
  labs(x = NULL , y ="% Granular Apoptotic") + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold", angle = 90, hjust = 1),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) +
  scale_shape_manual(values = c(15,16,17)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,20)) +
  scale_x_discrete(labels = c("BEADS_LPS"="Beads,<br> LPS",
                              "Control_hemo"="Control",
                              "Dermo"="*P. mar.*",
                              "Dermo_GDC" ="*P. mar.*,<br>GDC-0152",
                              "Dermo_ZVAD"= "*P. mar.*,<br>Z-VAD-fmk")) 
 # scale_fill_manual(name="Treatment", labels=c("Beads and LPS","Control","*P. mar*",
 #                                              "*P. mar*, GDC-0152",
 #                                              "*P. mar*, Z-VAD-fmk"), values=c("#88bf3b", "#6c81d9","#5b2c90",
 #                                                                                              "#ba4b41","#bfac3e")) 

Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_sd_multipanel <- 
  Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_sd_multipanel  + 
  theme(axis.text.x=ggtext::element_markdown(),
        legend.text = ggtext::element_markdown()) 

# Perform anova with Tukey test instead and generate stats dataframe
Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_sd_AOV <- aov(Percent_of_this_plot_arcsine ~ Treat + ID, Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_sd)
summary(Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_sd_AOV)
stat_test_tukey <- tukey_hsd(Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_sd_AOV) %>%
  add_significance(p.col = "p.adj")

# take only the significant columns
stat_test_tukey <- stat_test_tukey[c(1,2,5,8,9,10),]

Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_sd_multipanel_sig <- 
  Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_sd_multipanel + stat_pvalue_manual(
    stat_test_tukey, label = "{p.adj} {p.adj.signif}",  tip.length = 0.01, y.position = c(8, 9, 10,11,12,14), size = 3) +
  # add overall anova values 
  #stat_compare_means(method= "anova") +
  labs(subtitle = "Tukey HSD, Arcsine Percent ~ Treat + Pool")

## Plot pool since my statistical analysis revealed Pool is significant 
Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_pool_plot <- 
  Dermo_Inhibitor_2020_APOP_join_granular_apoptotic %>%
  # plot only UR and remove the control hemocyte
  filter(Gate == "Q16-UR") %>% 
  ggplot(., aes(y=Percent_of_this_plot, x=Treat, fill=ID)) + geom_col(position = "dodge") + 
  xlab("Treatment") +
  ylab("Percent of Granulocytes") + 
  ggtitle("Percent of Apoptotic Granulocytes") + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,10), breaks = c(0,2,4,6,8,10)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  scale_fill_manual(name="Cell Type", labels=c("Pool 1", "Pool 2", "Pool3"), 
                    values = c("#7f63b8", "#50b47b", "#ba583b")) 

#save
ggsave(plot = Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_pool_plot, device = "tiff", filename = "Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_pool_plot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 10)

#### 2020 Dermo and Inhibitors CASPASE ASSAY Statistics and Plotting ####

Dermo_Inhibitor_2020_CASP_join

## Start by comparing the number of parasite cells in Q10-LR in the stained dermo plots to the Q10-LL in the unstained Dermo plots

Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_perk <- Dermo_Inhibitor_2020_CASP_join %>% filter(Gate == "Q10-LL" & Treat == "PERK" ) %>% 
  mutate(sample_ID = paste(ID, Treat, sep = "-"))
Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_treat <- Dermo_Inhibitor_2020_CASP_join %>% filter(Gate == "Q10-LR") %>%
  filter(Treat == "Dermo" | Treat == "Dermo_GDC" | Treat == "Dermo_ZVAD") %>% 
  mutate(sample_ID = paste(ID, Treat, sep = "-"))
Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_perk_treat <- rbind(Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_perk, 
                                                                     Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_treat)
#Plot 
Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_perk_treat_plot <- 
  ggplot(Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_perk_treat, aes(x= sample_ID, y = Counts, fill = Gate)) + 
  geom_col(position = "dodge") +
  xlab("Treatment") +
  ylab("Cell Counts") + 
  ggtitle("Counts of Non-Caspase Active Perkinsus") + 
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12, angle= 90, hjust=1.0)) +
  theme(legend.text = element_text(size=12)) +
  scale_fill_manual(name="Cell Type", labels=c("Parasite Control", "Treatment Groups"), 
                    values = c("#50b47b",
                               "#7f63b8")) 
#save
ggsave(plot = Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_perk_treat_plot, device = "tiff", filename = "Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_perk_treat_plot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 10)
# Remember, each pool got a different amount of Dermo cells..even the ones for the same pool had slightly different amounts

## What is the counts ratio in the PERK only pool of caspase 3/7active to non caspase 3/7 active cells
Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_perk_all <- Dermo_Inhibitor_2020_CASP_join %>% 
  filter(Gate == "Q10-LL" | Gate == "Q10-UL") %>% filter(Treat == "PERK" ) %>% 
  mutate(sample_ID = paste(ID, Treat, sep = "-"))
# calculate ratio
Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_perk_all %>% group_by(ID) %>%
  summarize(CASP_ratio = Counts[Gate == "Q10-UL"] / Counts[Gate == "Q10-LL"]) %>% ungroup() %>% summarise(mean = mean(CASP_ratio))
#CASP_ratio: average is 40.9% caspase active
# CASP ratio is between  of the perkinsus cells are Caspase 3/7 active...so I can roughly calculate what the ratio would
# be for the other assay

## Plot percent of Caspase active perkinsus and format for multi-panel figure
Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_perk_all_sd <- Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_perk_all %>%
  group_by(Gate) %>% mutate(mean = mean(Percent_of_this_plot), sd = sd(Percent_of_this_plot))

Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_perk_all_multipanel <- 
  ggplot(data=Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_perk_all_sd,
         aes(fill = factor(Gate,levels = c("Q10-UL","Q10-LL")), y=Percent_of_this_plot, x=Assay)) + geom_bar(position="fill", stat="identity")  + 
  labs(x = NULL , y ="% P. marinus") + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold"),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold"),
        legend.position = "bottom") +
  guides(fill=guide_legend(nrow=2, title.position = "top")) +
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels = c("CASP" = "*P. mar.*")) + 
  scale_fill_manual(name="Caspase<br>Assay", labels=c("Caspase 3/7 Active", "Non-Caspase 3/7 Active"),
                    values=c("#a44f9a","#5b2c90")) 

Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_perk_all_multipanel_sig <- 
  Dermo_Inhibitor_2020_CASP_join_non_CASP_granular_perk_all_multipanel + 
  theme(axis.text.x=ggtext::element_markdown(),
        legend.text = ggtext::element_markdown(),
        legend.title = ggtext::element_markdown()) 

### Plotting granular cell data

# Make plot of percentage of casp hemocytes with parasite and without
Dermo_Inhibitor_2020_CASP_join_granular_casp_plot <- Dermo_Inhibitor_2020_CASP_join %>%
  # plot percentage of casp hemocytes with and without parasite
  filter(Gate == "Q10-UR" | Gate == "Q10-UL") %>%
  ggplot(data=., aes(y=Percent_of_this_plot, x=Treat, color=Gate)) + geom_point(position=position_dodge(width=0.75)) + 
  geom_boxplot() +
  xlab("Treatment") +
  ylab("Percent of Granulocytes") + 
  ggtitle("Percent of casp Granulocytes") + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,40), breaks = c(0,5,10,20,30,40)) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=12), 
        axis.title.y=element_text(family="serif",size=12),
        axis.title.x=element_text(family="serif",size=12),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=12)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  scale_color_manual(name="Cell Type", labels=c("casp\nGranular", "casp Granular,\nwith P. mar., or Beads"), 
                     values = c("#56b464", "#5b2c90")) 

# color options
#"#cc57b4", "#88bf3b", "#aa4dce", "#56b464", "#5b2c90"
#"#c89832", "#5a6ee6", "#ca4e33", "#7e78d4", "#cd4272"

#save
ggsave(plot = Dermo_Inhibitor_2020_CASP_join_granular_casp_plot, device = "tiff", filename = "Dermo_Inhibitor_2020_CASP_join_granular_casp_plot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 5, width = 10)


## Plot just the 10-UR data 

Dermo_Inhibitor_2020_CASP_join_granular_casp <- Dermo_Inhibitor_2020_CASP_join %>% filter(Gate == "Q10-UR") %>% filter(Treat !="UV") %>% filter(Treat !="PERK")

##Plot Caspase active granulocytes in format for multipanel figure with multiple comparisons run 
Dermo_Inhibitor_2020_CASP_join_granular_casp_sd <-   Dermo_Inhibitor_2020_CASP_join_granular_casp %>%
  filter(Gate == "Q10-UR") %>% group_by(Treat) %>% mutate(mean = mean(Percent_of_this_plot), sd = sd(Percent_of_this_plot))

Dermo_Inhibitor_2020_CASP_join_granular_casp_sd$Treat <- factor(Dermo_Inhibitor_2020_CASP_join_granular_casp_sd$Treat,
                                                                     levels = c("Control_hemo","BEADS_LPS","Dermo","Dermo_GDC","Dermo_ZVAD"))

Dermo_Inhibitor_2020_CASP_join_granular_casp_sd_multipanel <- 
  ggplot(data=Dermo_Inhibitor_2020_CASP_join_granular_casp_sd,
         aes(y=Percent_of_this_plot, x=Treat)) + 
  geom_bar(aes(fill=Treat), position="dodge", stat = "summary", fill = "#b84c3f")  + 
  geom_point(aes(x= Treat, shape = ID), size = 3) +
  labs(x = NULL , y ="% Granular Caspase 3/7 Active") + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold", angle = 90, hjust = 1),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) +
  scale_shape_manual(values = c(15,16,17)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,20)) +
  scale_x_discrete(labels = c("BEADS_LPS"="Beads,<br> LPS",
                              "Control_hemo"="Control",
                              "Dermo"="*P. mar.*",
                              "Dermo_GDC" ="*P. mar.*,<br>GDC-0152",
                              "Dermo_ZVAD"= "*P. mar.*,<br>Z-VAD-fmk")) 
  #scale_fill_manual(name="Treatment", labels=c("Beads and LPS","Control","*P. mar*",
  #                                             "*P. mar*, GDC-0152",
  #                                             "*P. mar*, Z-VAD-fmk"), values=c("black", "black","black",
  #                                                                              "black","black")) 
    # colors option: values=c("#88bf3b", "#6c81d9","#5b2c90","#ba4b41","#bfac3e")

Dermo_Inhibitor_2020_CASP_join_granular_casp_sd_multipanel <- 
  Dermo_Inhibitor_2020_CASP_join_granular_casp_sd_multipanel  + 
  theme(axis.text.x=ggtext::element_markdown(),
        legend.text = ggtext::element_markdown()) 

# Perform anova with Tukey test instead and generate stats dataframe
Dermo_Inhibitor_2020_CASP_join_granular_casp_sd_AOV <- aov(Percent_of_this_plot_arcsine ~ Treat + ID, Dermo_Inhibitor_2020_CASP_join_granular_casp_sd)
summary(Dermo_Inhibitor_2020_CASP_join_granular_casp_sd_AOV)
stat_test_tukey <- tukey_hsd(Dermo_Inhibitor_2020_CASP_join_granular_casp_sd_AOV) %>%
  add_significance(p.col = "p.adj")

# take only the significant columns
stat_test_tukey <- stat_test_tukey[c(1,2,5,8,9,10),]

Dermo_Inhibitor_2020_CASP_join_granular_casp_sd_multipanel_sig <- 
  Dermo_Inhibitor_2020_CASP_join_granular_casp_sd_multipanel + stat_pvalue_manual(
    stat_test_tukey, label = "{p.adj} {p.adj.signif}",  tip.length = 0.01, y.position = c(8, 9, 10,11,12,14), size = 3) +
  # add overall anova values 
  #stat_compare_means(method= "anova") +
  labs(subtitle = "Tukey HSD, Arcsine Percent ~ Treat + Pool")

#### 2020 Dermo and Inhibitors JC1 ASSAY statistics and plotting ####

Dermo_Inhibitor_2020_JC1_join

# Plot 18 is what we care about for comparing between treatments(Q28-UL is mitochondria permeabilization positive hemocytes alone) (Q28-UR is mitochondria permeabilization positive hemocytes with parasite) 
# however, to compare these treatments to control will be more difficult because this plot no longer applies
# The CCCP data also has to be plotted separately

### Plotting granular cell data due to the parasite across parasite treatment groups

# isolate just the parasite or bead treated samples where plot 18 applies
Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated <- Dermo_Inhibitor_2020_JC1_join %>% 
  filter(Gate == "Q28-UR") %>% filter(Treat !="UV_JC1") %>%  filter(Treat !="UV") %>% filter(Treat !="PERK") %>% filter(Treat != "Control_hemo") %>% 
  filter(Treat != "CCCP")

# Plot mitochondria permeabilized granulocytes in format for multipanel figure with multiple comparisons run 
Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_sd <- Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated %>%
  group_by(Treat) %>% mutate(mean = mean(Percent_of_this_plot), sd = sd(Percent_of_this_plot))

Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_sd$Treat <- factor(Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_sd$Treat,
                                                                levels = c("BEADS_LPS","Dermo","Dermo_GDC","Dermo_ZVAD"))
class(Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_sd$Percent_of_this_plot)

Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_sd_multipanel <- 
  ggplot(data=Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_sd,
         aes(y=Percent_of_this_plot, x=Treat)) + 
  geom_bar(aes(fill=Treat), position="dodge", stat = "summary", fill = "#52b77f")  + 
  geom_point(aes(x= Treat, shape = ID), size = 3) +
  labs(x = NULL , y ="% Granular Mitochondria Permeabilized") + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, face= "bold"),
        axis.title.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10, face= "bold", angle = 90, hjust = 1),
        legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 12, face= "bold")) +
  scale_shape_manual(values = c(15,16,17)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0,25)) +
  scale_x_discrete(labels = c("BEADS_LPS"="Beads,<br> LPS",
                              "Dermo"="*P. mar.*",
                              "Dermo_GDC" ="*P. mar.*,<br>GDC-0152",
                              "Dermo_ZVAD"= "*P. mar.*,<br>Z-VAD-fmk")) 

Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_sd_multipanel <- 
  Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_sd_multipanel  + 
  theme(axis.text.x=ggtext::element_markdown(),
        legend.text = ggtext::element_markdown()) 

# Perform anova with Tukey test instead and generate stats dataframe
Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_sd_AOV <- aov(Percent_of_this_plot_arcsine ~ Treat + ID, Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_sd)
summary(Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_sd_AOV)
stat_test_tukey <- tukey_hsd(Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_sd_AOV) %>%
  add_significance(p.col = "p.adj")



# take only the significant columns
stat_test_tukey <- stat_test_tukey %>% filter(p.adj <= 0.05)

Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_sd_multipanel_sig <- 
  Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_sd_multipanel + stat_pvalue_manual(
    stat_test_tukey, label = "{p.adj} {p.adj.signif}",  tip.length = 0.01, y.position = c(23,25), size = 3) +
  # add overall anova values 
  #stat_compare_means(method= "anova") +
  labs(subtitle = "Tukey HSD, Arcsine Percent ~ Treat + Pool")

## does T test make a difference?

Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_dermo_GDC <- Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated %>% filter(Treat == "Dermo" | Treat == "Dermo_GDC") %>%
  filter(ID == "Pool2" | ID == "Pool1") 

# No the high variance in this data set makes there be not that much difference 

#### 2020 Dermo and Inhibitors Multi-panel APOP, CASP, JC-1, figure ####

Flow_2020_multipanel <- cowplot::plot_grid(Dermo_Inhibitor_2020_APOP_join_granular_apoptotic_sd_multipanel_sig, NULL,
                                           Dermo_Inhibitor_2020_CASP_join_granular_casp_sd_multipanel_sig, NULL,
                                           Dermo_Inhibitor_2020_JC1_join_granular_JC1_treated_sd_multipanel_sig, NULL,
                                           ncol = 2, nrow = 3, labels = "AUTO", label_size = 16, label_fontface = "bold")

ggsave(Flow_2020_multipanel, device = "tiff", filename = "Flow_2020_multipanel.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/FIGURES",
       height = 15, width = 11.2 ) 

