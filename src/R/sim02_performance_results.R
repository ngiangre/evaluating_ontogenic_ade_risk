#' Title: "Evaluating risk detection methods to uncover ontogenic-mediated adverse drug effect mechanisms in children" study
#' 
#' Author details: Nicholas Giangreco
#' 
#' This script calculates the data description and performance 
#' statistics and figures (minus the real world validation results)

# PURPOSE -----------------------------------------------------------------

#' Derive simulation data for evaluating detection methods
#' 
#' deriving population adverse drug event (ADE or DE) reporting rate from 
#' the all pediatric FAERS (population) to generate a representative
#' distribution for random ADEs
#' 
#' these random ADEs will be "causal" 
#' we can then compare methods for assessing which have higher performance
#' giving evidence for hypothesized ADEs


# Load libraries and set variables and functions ----------------------------------------------------------

#' General loading of libraries, setting of global variables, and 
#' creaing simulaion data output directory

#install.packages("pacman")
pacman::p_load(tidyverse,data.table)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
out_dir <- paste0(data_dir,"ade_simulation_data/")
basename <- "sim_"

seed = 0
set.seed(seed)

library(doParallel)
cores=4
registerDoParallel(cores=cores)


theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))

score_colors <- 
  c("PRR" = RColorBrewer::brewer.pal(n=3,name = "Set1")[2],
    "GAM" = RColorBrewer::brewer.pal(n=3,name = "Set1")[1])

gam_score_colors <- 
  c("GAM" = RColorBrewer::brewer.pal(n=3,name = "Set1")[1],
    "full GAM" = RColorBrewer::brewer.pal(n=3,name = "Set2")[2],
    "tensor GAM" = RColorBrewer::brewer.pal(n=3,name = "Set2")[3])

direction_colors <- 
  c("0" = RColorBrewer::brewer.pal(n=3,name = "Set2")[1],
    "1" = RColorBrewer::brewer.pal(n=3,name = "Set2")[2])

dynamics_colors <- 
  c("uniform"=RColorBrewer::brewer.pal(n=5,name = "Dark2")[1],
    "increase"=RColorBrewer::brewer.pal(n=5,name = "Dark2")[2],
    "decrease"=RColorBrewer::brewer.pal(n=5,name = "Dark2")[3],
    "plateau"=RColorBrewer::brewer.pal(n=5,name = "Dark2")[4],
    "inverse_plateau"=RColorBrewer::brewer.pal(n=5,name = "Dark2")[5]
  )
dynamics_colors_no_uniform <- 
  c("increase"=RColorBrewer::brewer.pal(n=5,name = "Dark2")[2],
    "decrease"=RColorBrewer::brewer.pal(n=5,name = "Dark2")[3],
    "plateau"=RColorBrewer::brewer.pal(n=5,name = "Dark2")[4],
    "inverse_plateau"=RColorBrewer::brewer.pal(n=5,name = "Dark2")[5]
  )
dynamics_colors_with_na <- 
  c("uniform"=RColorBrewer::brewer.pal(n=5,name = "Dark2")[1],
    "increase"=RColorBrewer::brewer.pal(n=5,name = "Dark2")[2],
    "decrease"=RColorBrewer::brewer.pal(n=5,name = "Dark2")[3],
    "plateau"=RColorBrewer::brewer.pal(n=5,name = "Dark2")[4],
    "inverse_plateau"=RColorBrewer::brewer.pal(n=5,name = "Dark2")[5],
    "NA" = "darkgray"
  )

auc_ci <- function(y_pred,y_true,stat="auc",nboot=100){
  sapply(1:nboot,function(x){
    set.seed(x)
    inds = sample(1:length(y_pred),length(y_pred),replace=T)
    if(length(unique(y_true[inds]))==2){
    ROCR::performance(
      ROCR::prediction(
        y_pred[inds],
        y_true[inds]),
      stat)@y.values[[1]]
    }else{
      c(NA,NA,NA)
    }
  }) %>% quantile(probs=c(0.025,0.5,0.975),na.rm = T)
}

auc_probability <- function(labels, scores, N=1e7,split=F){
  #https://www.r-bloggers.com/2016/11/calculating-auc-the-area-under-a-roc-curve/
  p <- sample(scores[labels], N, replace=TRUE)
  n <- sample(scores[!labels], N, replace=TRUE)
  tp = sum(p > n)
  tie = sum(p == n)/2
  fp = length(p) - tp
  # sum( (1 + sign(pos - neg))/2)/N # does the same thing
  if(split==T){
    c(tp/N,tie/N)
  }else{
    (tp + tie) / N # give partial credit for ties
  }
}

#https://www.r-bloggers.com/2016/11/calculating-auc-the-area-under-a-roc-curve/
simple_auc <- function(FPR, TPR){
  # inputs already sorted, best scores first 
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}

# Load pediatric data --------------------------------------------

#' Loading of FAERS pediatric population dataset
#' Extracted by code at https://zenodo.org/record/4464544#.YEExxV1Kg6A

drug_col = "atc_concept_id"
drug_col_name = "atc_concept_name"
rx_col = "meddra_concept_id"
rx_col_name = "meddra_concept_name"
stage_col = "nichd"
age_col = "age"
id_col="safetyreportid"
raw_data <- 
  fread(paste0(data_dir,"sim_pediatric_FAERS.csv.gz"))

raw_data$years <- raw_data[,ceiling(age)]
raw_data <- raw_data[order(age)]

# Define category levels -------------------------------------------------------------

category_levels = list()

category_levels[["years"]] =
  seq(1,21,1)

category_levels[["nichd"]] =
  c("term_neonatal","infancy",
    "toddler","early_childhood",
    "middle_childhood","early_adolescence",
    "late_adolescence")

category_levels[["reporter_qualification"]] =
  c("Consumer or non-health professional",
    "Lawyer","Other health professional",
    "Pharmacist","Physician")

# Load positive ADEs ----------------------------------------------

positive_ade_pairs <- 
  fread(paste0(data_dir,basename,"positive_ade_pairs.csv"))


# Load negative ADEs ----------------------------------------------------

negative_ade_pairs <- 
  fread(paste0(data_dir,basename,"negative_ade_pairs.csv"))

# Process positive control data -------------------------------------------------------
l=0.75
stage_col="nichd"
spikein_cols=c("uniform","increase","decrease","plateau","inverse_plateau")
for(spikein_col in spikein_cols){
  base = paste0("*",stage_col,"_simulation_positive_l",l,"_",spikein_col,"_rate.rds")
  gam_data <- lapply(
    list.files(out_dir,base),
    function(x){
      readRDS(paste0(out_dir,x))$gam
      }) %>% 
    bind_rows() %>% 
    data.table() %>% 
    merge(positive_ade_pairs,all.y=T)
  
  glm_data <- lapply(
    list.files(out_dir,base),
    function(x){
      readRDS(paste0(out_dir,x))$glm
    }) %>% 
    bind_rows() %>% 
    data.table() %>% 
    merge(positive_ade_pairs,all.y=T)
  
  poly2_data <- lapply(
    list.files(out_dir,base),
    function(x){
      readRDS(paste0(out_dir,x))$poly2
    }) %>% 
    bind_rows() %>% 
    data.table() %>% 
    merge(positive_ade_pairs,all.y=T)
  colnames(poly2_data)[grepl("poly",colnames(poly2_data))] <- 
    gsub("poly","poly2",colnames(poly2_data)[grepl("poly",colnames(poly2_data))])
  
  poly3_data <- lapply(
    list.files(out_dir,base),
    function(x){
      readRDS(paste0(out_dir,x))$poly3
    }) %>% 
    bind_rows() %>% 
    data.table() %>% 
    merge(positive_ade_pairs,all.y=T)
  colnames(poly3_data)[grepl("poly",colnames(poly3_data))] <- 
    gsub("poly","poly3",colnames(poly3_data)[grepl("poly",colnames(poly3_data))])
  
  poly4_data <- lapply(
    list.files(out_dir,base),
    function(x){
      readRDS(paste0(out_dir,x))$poly4
    }) %>% 
    bind_rows() %>% 
    data.table() %>% 
    merge(positive_ade_pairs,all.y=T)
  colnames(poly4_data)[grepl("poly",colnames(poly4_data))] <- 
    gsub("poly","poly4",colnames(poly4_data)[grepl("poly",colnames(poly4_data))])
  
  prr_data <- lapply(
    list.files(out_dir,base),
    function(x){
      readRDS(paste0(out_dir,x))$prr
      }) %>% 
    bind_rows() %>% 
    data.table() %>% 
    merge(positive_ade_pairs,all.y=T)
  
  positive_data <- 
    merge(prr_data[,
                   c(drug_col,rx_col,stage_col,"a","b","c","d",
                     colnames(prr_data)[grepl("PRR",colnames(prr_data))]
                   ),with=F],
          gam_data[,
               c(drug_col,rx_col,stage_col,"D","E","DE","Ej","Fij","Tij",
                 colnames(gam_data)[grepl("gam",colnames(gam_data))]
                 ),with=F],
      by=c(drug_col,rx_col,stage_col),
      suffixes = c("","_gam"),
      allow.cartesian = T) %>% 
    unique()  %>% 
    merge(
      poly2_data[,
                 c(drug_col,rx_col,stage_col,
                   colnames(poly2_data)[grepl("poly2",colnames(poly2_data))]
                 ),with=F],
      by=c(drug_col,rx_col,stage_col),
      allow.cartesian = T) %>% 
    merge(
      poly3_data[,
                 c(drug_col,rx_col,stage_col,
                   colnames(poly3_data)[grepl("poly3",colnames(poly3_data))]
                 ),with=F],
      by=c(drug_col,rx_col,stage_col),
      allow.cartesian = T) %>% 
    merge(
      poly4_data[,
                 c(drug_col,rx_col,stage_col,
                   colnames(poly4_data)[grepl("poly4",colnames(poly4_data))]
                 ),with=F],
      by=c(drug_col,rx_col,stage_col),
      allow.cartesian = T) %>% 
    merge(
      glm_data[,
                 c(drug_col,rx_col,stage_col,
                   colnames(glm_data)[grepl("glm",colnames(glm_data))]
                 ),with=F],
      by=c(drug_col,rx_col,stage_col),
      allow.cartesian = T)

  positive_data <- 
    merge(positive_data,
          prr_data[,.(atc_concept_id,meddra_concept_id,nichd,time_prr = time, aic_prr = NA)],
          by=c(drug_col,rx_col,stage_col)
          )
  positive_data <- 
    merge(positive_data,
          gam_data[,.(atc_concept_id,meddra_concept_id,nichd,time_gam = time, aic_gam = AIC)],
          by=c(drug_col,rx_col,stage_col)
    )
  positive_data <- 
    merge(positive_data,
          poly2_data[,
                     .(atc_concept_id,meddra_concept_id,nichd,time_poly2 = time, aic_poly2 = AIC)
          ],
          by=c(drug_col,rx_col,stage_col)
    )
  
  positive_data <- 
    merge(positive_data,
          poly3_data[,
                     .(atc_concept_id,meddra_concept_id,nichd,time_poly3 = time, aic_poly3 = AIC)
          ],
          by=c(drug_col,rx_col,stage_col)
    )
  
  positive_data <- 
    merge(positive_data,
          poly4_data[,
                     .(atc_concept_id,meddra_concept_id,nichd,time_poly4 = time, aic_poly4 = AIC)
          ],
          by=c(drug_col,rx_col,stage_col)
    )
  
  positive_data <- 
    merge(positive_data,
          glm_data[,
                     .(atc_concept_id,meddra_concept_id,nichd,time_glm = time, aic_glm = AIC)
          ],
          by=c(drug_col,rx_col,stage_col)
    )
  
  positive_data$spikein <- spikein_col
  positive_data$ade <- paste0(positive_data$atc_concept_id,"_",positive_data$meddra_concept_id)
  
  positive_data[nichd=="all"] %>% fwrite(paste0(data_dir,basename,"positive_",spikein_col,"_data.csv"))
    
  positive_data[nichd!="all"] %>% fwrite(paste0(data_dir,basename,"positive_",stage_col,"_",spikein_col,"_data.csv"))
    
}

spikein_col <- "original"
base = paste0("*",stage_col,"_simulation_",spikein_col,"_positive.rds")
gam_data <- lapply(
  list.files(out_dir,base),
  function(x){
    readRDS(paste0(out_dir,x))$gam
  }) %>% 
  bind_rows() %>% 
  data.table() %>% 
  merge(positive_ade_pairs,all.y=T)

glm_data <- lapply(
  list.files(out_dir,base),
  function(x){
    readRDS(paste0(out_dir,x))$glm
  }) %>% 
  bind_rows() %>% 
  data.table() %>% 
  merge(positive_ade_pairs,all.y=T)

poly2_data <- lapply(
  list.files(out_dir,base),
  function(x){
    readRDS(paste0(out_dir,x))$poly2
  }) %>% 
  bind_rows() %>% 
  data.table() %>% 
  merge(positive_ade_pairs,all.y=T)
colnames(poly2_data)[grepl("poly",colnames(poly2_data))] <- 
  gsub("poly","poly2",colnames(poly2_data)[grepl("poly",colnames(poly2_data))])

poly3_data <- lapply(
  list.files(out_dir,base),
  function(x){
    readRDS(paste0(out_dir,x))$poly3
  }) %>% 
  bind_rows() %>% 
  data.table() %>% 
  merge(positive_ade_pairs,all.y=T)
colnames(poly3_data)[grepl("poly",colnames(poly3_data))] <- 
  gsub("poly","poly3",colnames(poly3_data)[grepl("poly",colnames(poly3_data))])

poly4_data <- lapply(
  list.files(out_dir,base),
  function(x){
    readRDS(paste0(out_dir,x))$poly4
  }) %>% 
  bind_rows() %>% 
  data.table() %>% 
  merge(positive_ade_pairs,all.y=T)
colnames(poly4_data)[grepl("poly",colnames(poly4_data))] <- 
  gsub("poly","poly4",colnames(poly4_data)[grepl("poly",colnames(poly4_data))])

prr_data <- lapply(
  list.files(out_dir,base),
  function(x){
    readRDS(paste0(out_dir,x))$prr
  }) %>% 
  bind_rows() %>% 
  data.table() %>% 
  merge(positive_ade_pairs,all.y=T)

positive_data <- 
  merge(prr_data[,
                 c(drug_col,rx_col,stage_col,"a","b","c","d",
                   colnames(prr_data)[grepl("PRR",colnames(prr_data))]
                 ),with=F],
        gam_data[,
             c(drug_col,rx_col,stage_col,"D","E","DE","Ej","Fij","Tij",
               colnames(gam_data)[grepl("gam",colnames(gam_data))]
             ),with=F],
    by=c(drug_col,rx_col,stage_col),
    suffixes = c("","_gam"),
    allow.cartesian = T) %>% 
  unique() %>% 
  merge(
    poly2_data[,
               c(drug_col,rx_col,stage_col,
                 colnames(poly2_data)[grepl("poly2",colnames(poly2_data))]
               ),with=F],
    by=c(drug_col,rx_col,stage_col),
    allow.cartesian = T) %>% 
  merge(
    poly3_data[,
               c(drug_col,rx_col,stage_col,
                 colnames(poly3_data)[grepl("poly3",colnames(poly3_data))]
               ),with=F],
    by=c(drug_col,rx_col,stage_col),
    allow.cartesian = T) %>% 
  merge(
    poly4_data[,
               c(drug_col,rx_col,stage_col,
                 colnames(poly4_data)[grepl("poly4",colnames(poly4_data))]
               ),with=F],
    by=c(drug_col,rx_col,stage_col),
    allow.cartesian = T) %>% 
  merge(
    glm_data[,
             c(drug_col,rx_col,stage_col,
               colnames(glm_data)[grepl("glm",colnames(glm_data))]
             ),with=F],
    by=c(drug_col,rx_col,stage_col),
    allow.cartesian = T)

positive_data <- 
  merge(positive_data,
        prr_data[,.(atc_concept_id,meddra_concept_id,nichd,time_prr = time, aic_prr = NA)],
        by=c(drug_col,rx_col,stage_col)
  )
positive_data <- 
  merge(positive_data,
        gam_data[,.(atc_concept_id,meddra_concept_id,nichd,time_gam = time, aic_gam = AIC)],
        by=c(drug_col,rx_col,stage_col)
  )

positive_data <- 
  merge(positive_data,
        poly2_data[,
                   .(atc_concept_id,meddra_concept_id,nichd,time_poly2 = time, aic_poly2 = AIC)
        ],
        by=c(drug_col,rx_col,stage_col)
  )

positive_data <- 
  merge(positive_data,
        poly3_data[,
                   .(atc_concept_id,meddra_concept_id,nichd,time_poly3 = time, aic_poly3 = AIC)
        ],
        by=c(drug_col,rx_col,stage_col)
  )

positive_data <- 
  merge(positive_data,
        poly4_data[,
                   .(atc_concept_id,meddra_concept_id,nichd,time_poly4 = time, aic_poly4 = AIC)
        ],
        by=c(drug_col,rx_col,stage_col)
  )

positive_data <- 
  merge(positive_data,
        glm_data[,
                 .(atc_concept_id,meddra_concept_id,nichd,time_glm = time, aic_glm = AIC)
        ],
        by=c(drug_col,rx_col,stage_col)
  )

positive_data$spikein <- spikein_col
positive_data$ade <- paste0(positive_data$atc_concept_id,"_",positive_data$meddra_concept_id)

positive_data[nichd!="all"] %>% 
  fwrite(paste0(data_dir,basename,"positive_",stage_col,"_",spikein_col,"_data.csv"))

# Process negative control data -------------------------------------------------------

base = paste0("*",stage_col,"_simulation_negative.rds")
gam_data <- lapply(
  list.files(out_dir,base),
  function(x){
    readRDS(paste0(out_dir,x))$gam
  }) %>% 
  bind_rows() %>% 
  data.table() %>% 
  merge(negative_ade_pairs,all.y=T)

glm_data <- lapply(
  list.files(out_dir,base),
  function(x){
    readRDS(paste0(out_dir,x))$glm
  }) %>% 
  bind_rows() %>% 
  data.table() %>% 
  merge(negative_ade_pairs,all.y=T)

poly2_data <- lapply(
  list.files(out_dir,base),
  function(x){
    readRDS(paste0(out_dir,x))$poly2
  }) %>% 
  bind_rows() %>% 
  data.table() %>% 
  merge(negative_ade_pairs,all.y=T)
colnames(poly2_data)[grepl("poly",colnames(poly2_data))] <- 
  gsub("poly","poly2",colnames(poly2_data)[grepl("poly",colnames(poly2_data))])

poly3_data <- lapply(
  list.files(out_dir,base),
  function(x){
    readRDS(paste0(out_dir,x))$poly3
  }) %>% 
  bind_rows() %>% 
  data.table() %>% 
  merge(negative_ade_pairs,all.y=T)
colnames(poly3_data)[grepl("poly",colnames(poly3_data))] <- 
  gsub("poly","poly3",colnames(poly3_data)[grepl("poly",colnames(poly3_data))])

poly4_data <- lapply(
  list.files(out_dir,base),
  function(x){
    readRDS(paste0(out_dir,x))$poly4
  }) %>% 
  bind_rows() %>% 
  data.table() %>% 
  merge(negative_ade_pairs,all.y=T)
colnames(poly4_data)[grepl("poly",colnames(poly4_data))] <- 
  gsub("poly","poly4",colnames(poly4_data)[grepl("poly",colnames(poly4_data))])

prr_data <- lapply(
  list.files(out_dir,base),
  function(x){
    readRDS(paste0(out_dir,x))$prr
  }) %>% 
  bind_rows() %>% 
  data.table() %>% 
  merge(negative_ade_pairs,all.y=T)

negative_data <- 
  merge(prr_data[,
                 c(drug_col,rx_col,stage_col,"a","b","c","d",
                   colnames(prr_data)[grepl("PRR",colnames(prr_data))]
                 ),with=F],
        gam_data[,
             c(drug_col,rx_col,stage_col,"D","E","DE","Ej","Fij","Tij",
               colnames(gam_data)[grepl("gam",colnames(gam_data))]
             ),with=F],
    by=c(drug_col,rx_col,stage_col),
    suffixes = c("","_gam"),
    allow.cartesian = T) %>% 
  unique() %>% 
  merge(
    poly2_data[,
               c(drug_col,rx_col,stage_col,
                 colnames(poly2_data)[grepl("poly2",colnames(poly2_data))]
               ),with=F],
    by=c(drug_col,rx_col,stage_col),
    allow.cartesian = T) %>% 
  merge(
    poly3_data[,
               c(drug_col,rx_col,stage_col,
                 colnames(poly3_data)[grepl("poly3",colnames(poly3_data))]
               ),with=F],
    by=c(drug_col,rx_col,stage_col),
    allow.cartesian = T) %>% 
  merge(
    poly4_data[,
               c(drug_col,rx_col,stage_col,
                 colnames(poly4_data)[grepl("poly4",colnames(poly4_data))]
               ),with=F],
    by=c(drug_col,rx_col,stage_col),
    allow.cartesian = T) %>% 
  merge(
    glm_data[,
               c(drug_col,rx_col,stage_col,
                 colnames(glm_data)[grepl("glm",colnames(glm_data))]
               ),with=F],
    by=c(drug_col,rx_col,stage_col),
    allow.cartesian = T)

negative_data <- 
  merge(negative_data,
        prr_data[,.(atc_concept_id,meddra_concept_id,nichd,time_prr = time, aic_prr = NA)],
        by=c(drug_col,rx_col,stage_col)
  )
negative_data <- 
  merge(negative_data,
        gam_data[,.(atc_concept_id,meddra_concept_id,nichd,time_gam = time, aic_gam = AIC)],
        by=c(drug_col,rx_col,stage_col)
  )
negative_data <- 
  merge(negative_data,
        poly2_data[,
                   .(atc_concept_id,meddra_concept_id,nichd,time_poly2 = time, aic_poly2 = AIC)
        ],
        by=c(drug_col,rx_col,stage_col)
  )
negative_data <- 
  merge(negative_data,
        poly3_data[,
                   .(atc_concept_id,meddra_concept_id,nichd,time_poly3 = time, aic_poly3 = AIC)
        ],
        by=c(drug_col,rx_col,stage_col)
  )
negative_data <- 
  merge(negative_data,
        poly4_data[,
                   .(atc_concept_id,meddra_concept_id,nichd,time_poly4 = time, aic_poly4 = AIC)
        ],
        by=c(drug_col,rx_col,stage_col)
  )
negative_data <- 
  merge(negative_data,
        glm_data[,
                   .(atc_concept_id,meddra_concept_id,nichd,time_glm = time, aic_glm = AIC)
        ],
        by=c(drug_col,rx_col,stage_col)
  )
negative_data$ade <- paste0(negative_data$atc_concept_id,"_",negative_data$meddra_concept_id)

negative_data[nichd=="all"] %>% fwrite(paste0(data_dir,basename,"negative_data.csv"))
negative_data[nichd!="all"] %>% fwrite(paste0(data_dir,basename,"negative_",stage_col,"_data.csv"))


# Load processed data ------------------------------------------------

negative_data <- fread(paste0(data_dir,basename,"negative_data.csv"))
stage_negative_data <- fread(paste0(data_dir,basename,"negative_",stage_col,"_data.csv"))

spikein_cols=c("increase","decrease","plateau","uniform","inverse_plateau","original")
positive_data <- NULL
for(spikein_col in spikein_cols){
  tmp <- fread(paste0(data_dir,basename,"positive_",spikein_col,"_data.csv"))
  tmp$spikein <- spikein_col
  tmp$gam_score <- tmp$gam_score %>% as.numeric()
  positive_data <- bind_rows(
    positive_data,
    tmp
  )
}

stage_col="nichd"
spikein_cols=c("increase","decrease","plateau","uniform","inverse_plateau","original")
stage_positive_data <- NULL
for(spikein_col in spikein_cols){
  tmp <- fread(paste0(data_dir,basename,"positive_",stage_col,"_",spikein_col,"_data.csv"))
  tmp$spikein <- spikein_col
  tmp$gam_score <- tmp$gam_score %>% as.numeric()
  stage_positive_data <- bind_rows(
    stage_positive_data,
    tmp
  )
}

stage_positive_data[,stage_col] <- 
  factor(stage_positive_data[,get(stage_col)],levels=category_levels[[stage_col]])

# Define high/low reporting rates at stages per spikein -------------------

stage_spikein_class <- 
  bind_rows(
    data.table(nichd=category_levels[[stage_col]],
               spikein="decrease",class=c(1,1,NA,NA,NA,0,0)),
    data.table(nichd=category_levels[[stage_col]],
               spikein="increase",class=c(0,0,NA,NA,NA,1,1)),
    data.table(nichd=category_levels[[stage_col]],
               spikein="plateau",class=c(0,NA,1,1,1,NA,0)),
    data.table(nichd=category_levels[[stage_col]],
               spikein="inverse_plateau",class=c(1,NA,0,0,0,NA,1))
  )

stage_spikein_class %>% 
  na.omit() %>% 
  dcast(spikein ~ factor(nichd,levels=category_levels[[stage_col]]))

# Stage score summary ------------------------------------------------

tmp <- 
  bind_rows(
    stage_positive_data[
      spikein %in% stage_spikein_class[,unique(spikein)]
    ][,.(ade,nichd,PRR,gam_score,spikein,control="positive ADEs")],
    stage_negative_data[,.(ade,nichd,PRR,gam_score,spikein="NA",control="negative ADEs")] 
  )%>%
  pivot_longer(cols=c("PRR","gam_score")) %>% 
  data.table() %>% 
  .[,.('Score==0' = sum(value==0,na.rm = T),
       'Score>0' = sum(value!=0,na.rm = T),
       'Score==NaN' = sum(is.na(value))),
    .(spikein,control,name,nichd)
  ] %>% 
  pivot_longer(cols=c('Score==0','Score>0','Score==NaN'),names_to="summary",values_to="summary_value") %>% 
  data.table()

tmp$summary_value_percent <- 
  tmp[,ifelse(control=='positive ADEs',round(summary_value/500,2),round(summary_value/10000,2))]

tmp[name=="PRR"] %>% 
  .[,quantile(summary_value_percent,c(0.025),na.rm=T),.(summary,control)]
tmp[name=="PRR"] %>% 
  .[,mean(summary_value_percent),.(summary,control)]
tmp[name=="PRR"] %>% 
  .[,quantile(summary_value_percent,c(0.975),na.rm=T),.(summary,control)]

g <- 
tmp[name=="PRR"] %>% 
  ggplot(aes(
    factor(nichd,levels=c(category_levels[[stage_col]])),
    summary_value_percent,fill=spikein)) +
  geom_bar(stat="identity",position="dodge") +
  facet_grid(control ~ summary,scales="free") +
  guides(fill=guide_legend(title="Drug event dynamics class",title.position = "top")) +
  scale_fill_manual(values=dynamics_colors_with_na) +
  scale_y_continuous(labels=scales::percent) +
  xlab("") +
  ylab("Percent of Scores") +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle=45,hjust=1,vjust=1)
  )
ggsave(paste0(img_dir,basename,"prr_score_summary_faers_controls_by_spikein_stages.png"),g,width=12,height=6)


tmp <- bind_rows(
stage_positive_data[spikein %in% names(dynamics_colors_no_uniform),
                    .(method="GAM",type="score",ade,nichd,
                      score = gam_score,spikein,control="positive")],
stage_positive_data[spikein %in% names(dynamics_colors_no_uniform),
                    .(method="GAM",type="90mse",ade,nichd,
                      score = gam_score_90mse,spikein,control="positive")],
stage_positive_data[spikein %in% names(dynamics_colors_no_uniform),
                    .(method="PRR",type="score",ade,nichd,
                      score = PRR,spikein,control="positive")],
stage_positive_data[spikein %in% names(dynamics_colors_no_uniform),
                    .(method="PRR",type="90mse",ade,nichd,
                      score = PRR_90mse,spikein,control="positive")],
lapply(names(dynamics_colors_no_uniform),
       function(x)stage_negative_data[,
                                      .(method="GAM",type="score",ade,nichd,
                                        score = gam_score,spikein=x,control="negative")
                                      ]) %>%
  bind_rows(),
lapply(names(dynamics_colors_no_uniform),
       function(x)stage_negative_data[,
                                      .(method="GAM",type="90mse",ade,nichd,
                                        score = gam_score_90mse,spikein=x,control="negative")
       ]) %>%
  bind_rows(),
lapply(names(dynamics_colors_no_uniform),
       function(x)stage_negative_data[,
                                      .(method="PRR",type="score",ade,nichd,
                                        score = PRR,spikein=x,control="negative")
       ]) %>%
  bind_rows(),
lapply(names(dynamics_colors_no_uniform),
       function(x)stage_negative_data[,
                                      .(method="PRR",type="90mse",ade,nichd,
                                        score = PRR_90mse,spikein=x,control="negative")
       ]) %>%
  bind_rows()
) 

g <- tmp %>% 
  .[type=="score"] %>% 
  ggplot(aes(score,fill=control)) +
  geom_density(alpha=0.5) +
  scale_x_continuous(trans="log10",labels=scales::comma) +
  facet_grid(method~spikein,scales="free")
ggsave(paste0(img_dir,basename,"score_spikein_by_control_distributions.png"),g,width=15,height=7)

g <- tmp %>% 
  .[type=="90mse"] %>% 
  ggplot(aes(score,fill=control)) +
  geom_density(alpha=0.5) +
  scale_x_continuous(trans="log10",labels=scales::comma) +
  facet_grid(method~spikein,scales="free")
ggsave(paste0(img_dir,basename,"90mse_spikein_by_control_distributions.png"),g,width=15,height=7)

tmp[name=="PRR"]
tmp[name=="gam_score"]

g <- stage_positive_data[
  DE==0 & spikein=="uniform",
  .(ade,nichd,gam_score,spikein)
] %>% 
  ggplot(aes(gam_score)) +
  geom_histogram() +
  scale_y_continuous(trans="log10") +
  facet_wrap(~factor(nichd,levels=category_levels$nichd),scales="free") +
  xlab("GAM score") +
  ylab("Number of ADEs")
ggsave(paste0(img_dir,basename,"gam_score_distribution_0reports_per_stage.png"),
       g,width=8,height=6)

g <- stage_positive_data[
  DE==0 & spikein=="uniform",
  .(ade,nichd,gam_score_90mse,spikein)
] %>% 
  ggplot(aes(gam_score_90mse)) +
  geom_histogram() +
  scale_y_continuous(trans="log10") +
  facet_wrap(~factor(nichd,levels=category_levels$nichd),scales="free") +
  xlab("90% lower bound GAM score") +
  ylab("Number of ADEs")
ggsave(paste0(img_dir,basename,"gam_score_90mse_distribution_0reports_per_stage.png"),
       g,width=8,height=6)

g <- stage_positive_data[
  DE<=5 & spikein=="uniform",
  .(ade,DE,nichd,gam_score,spikein)
] %>% 
  ggplot(aes(gam_score,fill=factor(DE))) +
  geom_density(alpha=0.5) +
  facet_wrap(.~factor(nichd,levels=category_levels$nichd),scales="free") +
  guides(fill=guide_legend(title="Number of reports",title.position = "top",nrow=1)) +
  xlab("GAM score") +
  ylab("Number of ADEs") +
  theme(
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"gam_score_distribution_<=5reports_per_stage.png"),
       g,width=12,height=7)

# Stage drug event summary ------------------------------

tmp <- 
  bind_rows(
    raw_data[,
             .(DE = .N),.(atc_concept_id,meddra_concept_id)][
               ,
               .(DE,atc_concept_id,meddra_concept_id,control="FAERS")
             ],
    positive_data[,.(DE,atc_concept_id,meddra_concept_id,control="Positive ADEs")] %>% unique(),
    negative_data[,.(DE,atc_concept_id,meddra_concept_id,control="Negative ADEs")] %>% unique()
  )

len=1e8
set.seed(seed)
pop = sample(1:len,tmp[control=="FAERS",DE],replace = T)
set.seed(seed)
neg = sample(1:len,tmp[control=="Negative ADEs",DE],replace=T)
set.seed(seed)
pos = sample(1:len,tmp[control=="Positive ADEs",DE],replace=T)

t.test(
  pop,
  neg,
  paired=F
)$p.value

set.seed(seed)
t.test(
  pop,
  pos,
  paired=F
)$p.value

g <- tmp %>% 
  ggplot(aes(control,DE)) +
  scale_y_continuous(trans='log10',breaks=c(1,10,100,300,500,1000,10000),labels=scales::label_number(accuracy=1,big.mark=",")) +
  geom_boxplot() +
  ggbeeswarm::geom_quasirandom(method="tukeyDense",shape=1) +
  xlab("") +
  ylab("Number of drug, event reports")
ggsave(paste0(img_dir,basename,"drug_event_reporting_faers_controls.png"),g,width=8,height=5)

g <- 
  bind_rows(
    stage_positive_data[,.(DE,atc_concept_id,meddra_concept_id,nichd,spikein,control="Positive ADEs")] %>% unique(),
    stage_negative_data[,.(DE,atc_concept_id,meddra_concept_id,nichd,spikein="NA",control="Negative ADEs")] %>% unique()
  ) %>% 
  .[spikein!="original"] %>% 
  ggplot(aes(spikein,DE,color=spikein)) +
  geom_boxplot() +
  scale_y_continuous(trans="log10") +
  facet_grid(~control,scales="free") +
  guides(color=guide_legend(title="Drug event dynamics class",title.position = "top")) +
  scale_color_manual(values=dynamics_colors_with_na) +
  xlab("") +
  ylab("Number of drug event reports") +
  theme(
    axis.text.x = element_blank(),
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"drug_event_reporting_faers_controls_by_spikein.png"),g,width=10,height=6)

dts <- NULL
for(sp in stage_positive_data[,unique(spikein)]){
  a = stage_positive_data[order(atc_concept_id,meddra_concept_id,nichd)][spikein==sp,DE]
  b = stage_negative_data[order(atc_concept_id,meddra_concept_id,nichd)][,DE]
  tt = t.test(a,b)
  dt <- data.table(test = "ttest",pvalue = tt$p.value,spikein = sp)
  dts <- bind_rows(dts,dt)
}
dts

g <- 
  bind_rows(
    stage_positive_data[
      spikein %in% c("uniform",stage_spikein_class[,unique(spikein)]),
      .(DE,atc_concept_id,meddra_concept_id,nichd,spikein,control="Positive ADEs")] %>% unique(),
    stage_negative_data[,.(DE,atc_concept_id,meddra_concept_id,nichd,spikein="NA",control="Negative ADEs")] %>% unique()
  ) %>% 
  ggplot(aes(factor(nichd,levels=category_levels[["nichd"]]),DE,color=spikein)) +
  geom_boxplot() +
  scale_y_continuous(trans="log10",labels=scales::comma) +
  facet_grid(control~.,scales="free") +
  scale_color_manual(values=dynamics_colors_with_na) +
  xlab("") +
  ylab("Number of drug event reports") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,hjust=1,vjust=1)
  )
ggsave(paste0(img_dir,basename,"drug_event_reporting_faers_controls_by_spikein_stages.png"),g,width=10,height=6)

# Stage score distribution and normality ------------------------------------------------

#' 

stage_score_correlation <- function(dat=stage_positive_data,method="PRR",score="PRR_95mse",type="95mse"){
  sub <- dat[,
             c(stage_col,"spikein",score,"ade"),
             with=F
             ]
  dts <- NULL
  for(spikein_col in spikein_cols){
    for(s in as.integer(sub[,unique(get(stage_col))])){
      sn <- levels(sub[,unique(get(stage_col))])[s]
      subsub = sub[nichd==sn] %>% na.omit()
      a <- subsub[spikein=="uniform"][is.finite(get(score)),get(score)] %>% as.numeric()
      b <- subsub[spikein==spikein_col][is.finite(get(score)),get(score)] %>% as.numeric()
      if(length(a) == 0 | length(b)==0){next}
      asamps = sapply(1:100,function(x){set.seed(x);sample(a,size=100) %>% mean(na.rm=T)})
      bsamps = sapply(1:100,function(x){set.seed(x);sample(b,size=100) %>% mean(na.rm=T)})
      tt <- t.test(b,a,conf.level = .95,conf.int = T,paired=F)
      dts <- bind_rows(
        dts,
        data.table(nichd = sn,level = s,
                   spikein = spikein_col,
                   test = "ttest",
                   statistic = tt$statistic %>% as.numeric(), 
                   estimate = tt$estimate,
                   pvalue = tt$p.value,
                   method = method,
                   type = type,
                   a = "uniform",
                   b = spikein_col,
                   lwr = quantile(bsamps-asamps,c(0.025))[1],
                   diff = mean(bsamps) - mean(asamps),
                   upr = quantile(bsamps-asamps,c(0.975))[1],
                   a_len = length(a),
                   b_len = length(b))
      )
    }
  }
  
  
  dts[spikein!="original" & spikein!="uniform"]
  
}

tmp = 
bind_rows(
  stage_score_correlation(stage_positive_data,"PRR","PRR","score"),
  stage_score_correlation(stage_positive_data,"PRR","PRR_90mse","90mse"),
  stage_score_correlation(stage_positive_data,"GAM","gam_score","score"),
  stage_score_correlation(stage_positive_data,"GAM","gam_score_90mse","90mse"),
  stage_score_correlation(stage_positive_data,"GLM","glm_score","score"),
  stage_score_correlation(stage_positive_data,"GLM","glm_score_90mse","90mse"),
  stage_score_correlation(stage_positive_data,"LowOrder2","poly2_score","score"),
  stage_score_correlation(stage_positive_data,"LowOrder2","poly2_score_90mse","90mse"),
  stage_score_correlation(stage_positive_data,"LowOrder3","poly3_score","score"),
  stage_score_correlation(stage_positive_data,"LowOrder3","poly3_score_90mse","90mse"),
  stage_score_correlation(stage_positive_data,"LowOrder4","poly4_score","score"),
  stage_score_correlation(stage_positive_data,"LowOrder4","poly4_score_90mse","90mse")
) 

g <- tmp[method=="GAM"] %>% 
  ggplot(aes(forcats::fct_reorder(nichd,level),diff,color=spikein)) +
  geom_point(position = position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1,position = position_dodge(width=0.4)) +
  geom_line(aes(group=spikein),position = position_dodge(width=0.4)) +
  facet_grid(type~.,scales="free") +
  guides(color=guide_legend(title="Drug event dynamics class",title.position = "top")) +
  scale_color_manual(values=dynamics_colors_no_uniform) +
  xlab("") +
  ylab("Mean score difference from random (uniform class scores)") +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  ) 
ggsave(paste0(img_dir,basename,"gam_score_summary_by_spikein_stages.png"),g,width=10,height=8)

g <- tmp[method=="PRR" & type=="90mse"] %>% 
  ggplot(aes(forcats::fct_reorder(nichd,level),diff,color=spikein)) +
  geom_point(position = position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1,position = position_dodge(width=0.4)) +
  geom_line(aes(group=spikein),position = position_dodge(width=0.4)) +
  facet_grid(type~.,scales="free") +
  scale_color_manual(values=dynamics_colors_no_uniform) +
  xlab("") +
  ylab("Mean score difference from random") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank()
  )
ggsave(paste0(img_dir,basename,"prr_score_summary_90mse_by_spikein_stages.png"),g,width=10,height=6)

g <- tmp[method=="PRR" & type=="score"] %>% 
  ggplot(aes(forcats::fct_reorder(nichd,level),diff,color=spikein)) +
  geom_point(position = position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1,position = position_dodge(width=0.4)) +
  geom_line(aes(group=spikein),position = position_dodge(width=0.4)) +
  facet_grid(type~.,scales="free") +
  coord_cartesian(ylim=c(-20,20)) +
  scale_color_manual(values=dynamics_colors_no_uniform) +
  xlab("") +
  ylab("Mean score difference from random") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"prr_score_summary_score_by_spikein_stages.png"),g,width=10,height=6)

g <- tmp[grepl("LowOrder",method)] %>% 
  ggplot(aes(forcats::fct_reorder(nichd,level),diff,color=spikein)) +
  geom_point(position = position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1,position = position_dodge(width=0.4)) +
  geom_line(aes(group=spikein),position = position_dodge(width=0.4)) +
  facet_grid(method~type,scales="free") +
  scale_color_manual(values=dynamics_colors_no_uniform) +
  xlab("") +
  ylab("Mean score difference from random") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank()
  )
ggsave(paste0(img_dir,basename,"loworder_score_summary_by_spikein_stages.png"),g,width=10,height=10)

g <- tmp[grepl("GLM",method)] %>% 
  ggplot(aes(forcats::fct_reorder(nichd,level),diff,color=spikein)) +
  geom_point(position = position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1,position = position_dodge(width=0.4)) +
  geom_line(aes(group=spikein),position = position_dodge(width=0.4)) +
  facet_grid(type~method,scales="free") +
  scale_color_manual(values=dynamics_colors_no_uniform) +
  xlab("") +
  ylab("Mean score difference from random") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank()
  )
ggsave(paste0(img_dir,basename,"glm_score_summary_by_spikein_stages.png"),g,width=10,height=10)

stage_score_dist_test <- function(dat=stage_positive_data,method="PRR",score="PRR_90mse",type="90mse"){
  sub <- dat[spikein!="original" & spikein!="uniform",
             c(stage_col,"spikein",score,"ade"),
             with=F
             ]
  dts <- NULL
  for(spikein_col in sub[,unique(spikein)]){
    for(st in sub[,unique(get(stage_col))]){
      subsub = sub[nichd==st] %>% na.omit()
      b <- subsub[spikein==spikein_col][is.finite(get(score)),get(score)] %>% as.numeric()
      bsamps = sapply(1:100,function(x){
        set.seed(x);sample(b,size=100,replace=T) %>% mean(na.rm=T)
        })
      if(method=="PRR"){
        bsamps <- log10(bsamps)
      }
      tt <- shapiro.test(bsamps)
      dts <- bind_rows(
        dts,
        data.table(nichd = st,
                   spikein = spikein_col,
                   test = tt$method,
                   statistic = tt$statistic %>% as.numeric(), 
                   estimate = tt$estimate,
                   pvalue = tt$p.value,
                   method = method,
                   type = type,
                   len = length(b))
      )
    }
  }
  
  
  dts
  
}

tmp = 
  bind_rows(
    stage_score_dist_test(stage_positive_data,"PRR","PRR","score"),
    stage_score_dist_test(stage_positive_data,"PRR","PRR_90mse","90mse"),
    stage_score_dist_test(stage_positive_data,"GAM","gam_score","score"),
    stage_score_dist_test(stage_positive_data,"GAM","gam_score_90mse","90mse")
  ) 

tmp[,
    quantile(pvalue,c(0.025),na.rm=T),
    .(method,type)
    ] %>% dcast(method ~ type)
tmp[,
    mean(pvalue,na.rm=T),
    .(method,type)
    ] %>% dcast(method ~ type)
tmp[,
    quantile(pvalue,c(0.975),na.rm=T),
    .(method,type)
    ] %>% dcast(method ~ type)

g <- 
tmp %>% 
  ggplot(aes(factor(nichd,levels=category_levels$nichd),-log10(pvalue),color=spikein)) +
  geom_point(position = position_dodge(width=0.4)) +
  geom_line(aes(group=spikein),position = position_dodge(width=0.4)) +
  facet_grid(type~method,scales="free") +
  guides(color=guide_legend(title="Drug event dynamics class",title.position = "top")) +
  scale_color_manual(values=dynamics_colors_no_uniform) +
  xlab("") +
  ylab("-log10(pvalue)")+
  ggtitle("") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"detection_method_score_normality_tests.png"),g,width=12,height=6)

# Power by effect size per method ------------------------------------------------

stage_score_power <- function(dat=stage_positive_data,method="PRR",score = "PRR",score_threshold=1,type='overall'){
  
  esize='Tij'
  dat = merge(stage_spikein_class,
              dat[,c(stage_col,"spikein",score,"ade",esize,"Ej","Fij"),with=F],
              all.x=T) %>% 
    .[is.finite(get(score))] %>% 
    na.omit()
  
  es <- dat[,get(esize)] %>% unique()
  es <- c(es,seq(min(es),max(es),length=10))
  es <- es[order(es)]
  es <- expand.grid(es) %>% data.table() %>% unique()
  colnames(es) <- c("esize")
  tmp <- foreach(i=1:nrow(es),.combine = "rbind") %dopar% {
    esthresh <- es[i][,esize]
    sub <- dat[get(esize)>=esthresh]
    pred_true <- sub[,get(score)]>score_threshold
    pred_false <- !pred_true
    cond_true <- sub[,class==1]
    cond_false <- !cond_true
    c(
      ((as.integer(pred_true)+as.integer(cond_true))==2) %>% which() %>% length(),
      ((as.integer(pred_false)+as.integer(cond_true))==2) %>% which() %>% length(),
      ((as.integer(pred_true)+as.integer(cond_false))==2) %>% which() %>% length(),
      ((as.integer(pred_false)+as.integer(cond_false))==2) %>% which() %>% length(),
      length(pred_true)
    )
  } %>% data.table()
  colnames(tmp) <- c("tp","fn","fp","tn","len")
  ts_es_perf <- cbind(es,tmp)
  ts_es_perf[,
             c("tpr","fnr","tnr","power") :=
               list(
                 tp/(tp+fn),
                 fn/(tp+fn),
                 tn/(tn+fp),
                 tp/(tp+fn)
               )
             ]
  ts_es_perf$size <- esize
  ts_es_perf$type <- type
  ts_es_perf$type_name <- score
  ts_es_perf$method <-  method
  ts_es_perf
}

tmp <- 
  bind_rows(
    stage_score_power(method="PRR",score = "PRR",type='score'),
    stage_score_power(method="PRR",score = "PRR_90mse",type='90mse'),
    stage_score_power(method="GAM",score = "gam_score",type='score'),
    stage_score_power(method="GAM",score = "gam_score_90mse",type='90mse')
  )

g <- tmp %>% 
  ggplot(aes(esize,power,color=method)) +
  geom_point() +
  geom_path() +
  facet_grid(type~.) +
  scale_color_manual(values=score_colors) +
  xlab("Effect size") +
  ylab("Power")
ggsave(paste0(img_dir,basename,"detection_method_power.png"),g,width=6,height=6)

# ADE dynamic detection performance: without power filter ------------------------------------------------

understanding_auc <- function(method="gam_score",score="gam_score",type="overall",N=1e6,seed=0){
  set.seed(seed)
  dts <- 
    foreach(spikein_col=stage_spikein_class[,unique(spikein)],.combine="rbind") %dopar% {
      p = 
        stage_positive_data[spikein==spikein_col] 
      neg=stage_negative_data
      y_true <- c(rep(1,nrow(p)),rep(0,nrow(neg)))
      y_pred <- c(p[,get(score)],neg[,get(score)])
      y_true <- y_true[is.finite(y_pred)]
      y_pred <- y_pred[is.finite(y_pred)]
      if(length(unique(y_true))==2){
        #https://www.r-bloggers.com/2016/11/calculating-auc-the-area-under-a-roc-curve/
        p <- sample(y_pred[as.logical(y_true)], N, replace=TRUE)
        n <- sample(y_pred[!as.logical(y_true)], N, replace=TRUE)
        kst = ks.test(n,p,alternative = "greater")
        tp = sum(p > n)
        tie = sum(p == n)/2
        fp = length(p) - tp
        dt <- 
          bind_rows(
            data.table(value=p,sample="positive control"),
            data.table(value=n,sample="negative control")
          )
        dt$ptp = tp
        dt$pfp = fp
        dt$pfpr = fp/(tp)
        dt$ptie = tie
        dt$pauc=(tp + tie) / N
        dt$spikein=spikein_col
        dt$method=method
        dt$score=score
        dt$type=type
        dt$test_statistic <- kst$statistic %>% unname
        dt$test_pvalue <- kst$p.value
        dt$test_method <- kst$method
        dt
      }
  }
  dts
}

tmp <- bind_rows(
  understanding_auc(method="PRR",score="PRR",type="score"),
  understanding_auc(method="PRR",score="PRR_90mse",type="90mse"),
  understanding_auc(method="GAM",score="gam_score",type="score"),
  understanding_auc(method="GAM",score="gam_score_90mse",type="90mse")
)

g <- tmp %>% 
  .[,.(spikein,method,type,sample,test_statistic)] %>% 
  unique() %>% 
  ggplot(aes(spikein,test_statistic,color=method)) +
  geom_point() +
  facet_grid(~type) +
  scale_color_manual(values=score_colors) +
  xlab("") +
  ylab("KS test statistic") 
ggsave(paste0(img_dir,basename,"detection_before_performance_kstest_control_distributions.png"),g,width=12,height=5)

g <- tmp %>% 
  .[type=="score"] %>% 
  ggplot(aes(value,fill=sample)) +
  geom_density(alpha=0.5) +
  scale_x_continuous(trans="log10",labels=scales::comma) +
  facet_grid(method~spikein) +
  guides(fill=guide_legend(title="Control")) +
  xlab("") +
  ylab("Score density") +
  ggtitle("Distribution of resampled scores used for calculating probabilistic AUC")
ggsave(paste0(img_dir,basename,"detection_before_performance_bootstrap_score_by_control_pauc.png"),g,width=15,height=7)

g <- tmp %>% 
  .[type=="90mse"] %>% 
  ggplot(aes(value,fill=sample)) +
  geom_density(alpha=0.5) +
  scale_x_continuous(trans="log10",labels=scales::comma) +
  facet_grid(method~spikein) +
  guides(fill=guide_legend(title="Control")) +
  xlab("") +
  ylab("Score density") +
  ggtitle("Distribution of resampled scores used for calculating probabilistic AUC")
ggsave(paste0(img_dir,basename,"detection_before_performance_bootstrap_90mse_score_by_control_pauc.png"),g,width=15,height=7)

g <- tmp %>% 
  ggplot(aes(spikein,value,color=sample)) +
  geom_boxplot(outlier.shape=NA) +
  scale_y_continuous(trans="log10",labels=scales::comma) +
  facet_grid(method~type,scales="free") +
  xlab("") +
  ylab("Score")
ggsave(paste0(img_dir,basename,"detection_before_performance_bootstrap_score_by_control_pauc_boxplot.png"),g,width=15,height=7)

g <- tmp %>% 
  ggplot(aes(spikein,value,color=sample)) +
  geom_boxplot(outlier.shape=NA) +
  geom_violin() +
  scale_y_continuous(trans="log10",labels=scales::comma) +
  facet_grid(method~type,scales="free") +
  xlab("") +
  ylab("Score")
ggsave(paste0(img_dir,basename,"detection_before_performance_bootstrap_score_by_control_pauc_boxplot_overlap.png"),g,width=15,height=7)

generate_filtered_perf_curve <- function(method="gam_score",score="gam_score",type="overall",x="prec",y="rec",stat="auc",thresh=0,high=T){
  
  if(high){rates <- c(1)}else{rates <- c(0,1)}
  
  dts <- 
    foreach(spikein_col=stage_spikein_class[,unique(spikein)],.combine="rbind") %dopar% {
      p = 
        stage_positive_data[spikein==spikein_col] %>% 
        merge(stage_spikein_class[class %in% rates],by=c("nichd","spikein"))
      neg <- stage_negative_data
      y_true <- c(rep(1,nrow(p)),rep(0,nrow(neg)))
      y_pred <- c(p[,get(score)],neg[,get(score)])
      y_true <- y_true[is.finite(y_pred)]
      y_pred <- y_pred[is.finite(y_pred)]
      if(length(unique(y_true))==2){
        score_thresh_perf_dt = 
          lapply(1:100,
                 function(i){
                   set.seed(i)
                   sinds = sample(1:length(y_pred),length(y_pred),replace=T)
                   pos_class_mean = mean(y_pred[sinds][y_true[sinds]==1])
                   neg_class_mean = mean(y_pred[sinds][y_true[sinds]==0])
                   cond_true <- y_true[sinds]==1
                   cond_false <- !cond_true
                   pred_true <- y_pred[sinds]>=thresh
                   pred_false <- !pred_true
                   tp <- ((as.integer(pred_true)+as.integer(cond_true))==2) %>% which() %>% length()
                   fn <- ((as.integer(pred_false)+as.integer(cond_true))==2) %>% which() %>% length()
                   fp <- ((as.integer(pred_true)+as.integer(cond_false))==2) %>% which() %>% length()
                   tn <- ((as.integer(pred_false)+as.integer(cond_false))==2) %>% which() %>% length()
                   pauc = auc_probability(as.logical(y_true),y_pred,N=1e5,split=T)
                   data.table(score_threshold=thresh,tp=tp,tn=tn,fp=fp,fn=fn,
                              tpr=(tp/(tp+fn)),fpr=(fp/(tn+fp)),tnr=(tn/(tn+fp)),fnr=(fn/(fn+tp)),
                              ppv=(fp/(fp+tp)),npv=(tn/(tn+fn)),
                              pauc = sum(pauc),ptp = pauc[1],ptie=pauc[2],
                              pos_class_mean = pos_class_mean,neg_class_mean=neg_class_mean)
                 }) %>% 
          bind_rows()
        pred <- ROCR::prediction(y_pred,y_true)
        perf <- ROCR::performance(pred,x,y)
        ci = auc_ci(y_pred,y_true,stat=stat)
        dt <- data.table(perf@y.values[[1]],perf@x.values[[1]],perf@alpha.values[[1]])
        colnames(dt) <- c(x,y,"threshold")
        dt$method <- method
        dt$type <- type
        dt$spikein <- spikein_col
        dt$N <- length(y_true)
        dt$auc_lwr <- ci[1]
        dt$auc <- ci[2]
        dt$auc_upr <- ci[3]
        dt$pauc_lwr <- score_thresh_perf_dt[,quantile(pauc,c(0.025))]
        dt$pauc <- score_thresh_perf_dt[,mean(pauc)]
        dt$pauc_upr <- score_thresh_perf_dt[,quantile(pauc,c(0.975))]
        dt$ptp_lwr <- score_thresh_perf_dt[,quantile(ptp,c(0.025))]
        dt$ptp <- score_thresh_perf_dt[,mean(ptp)]
        dt$ptp_upr <- score_thresh_perf_dt[,quantile(ptp,c(0.975))]
        dt$ptie_lwr <- score_thresh_perf_dt[,quantile(ptie,c(0.025))]
        dt$ptie <- score_thresh_perf_dt[,mean(ptie)]
        dt$ptie_upr <- score_thresh_perf_dt[,quantile(ptie,c(0.975))]
        dt$tpr_lwr <- score_thresh_perf_dt[,quantile(tpr,c(0.025),na.rm=T)] %>% unname
        dt$tpr_upr <- score_thresh_perf_dt[,quantile(tpr,c(0.975),na.rm=T)] %>% unname
        dt$tpr_mean <- score_thresh_perf_dt[,mean(tpr)]
        dt$tnr_lwr <- score_thresh_perf_dt[,quantile(tnr,c(0.025),na.rm=T)] %>% unname
        dt$tnr_upr <- score_thresh_perf_dt[,quantile(tnr,c(0.975),na.rm=T)] %>% unname
        dt$tnr_mean <- score_thresh_perf_dt[,mean(tnr)]
        dt$ppv_lwr <- score_thresh_perf_dt[,quantile(ppv,c(0.025),na.rm=T)] %>% unname
        dt$ppv_upr <- score_thresh_perf_dt[,quantile(ppv,c(0.975),na.rm=T)] %>% unname
        dt$ppv_mean <- score_thresh_perf_dt[,mean(ppv)]
        dt$npv_lwr <- score_thresh_perf_dt[,quantile(npv,c(0.025),na.rm=T)] %>% unname
        dt$npv_upr <- score_thresh_perf_dt[,quantile(npv,c(0.975),na.rm=T)] %>% unname
        dt$npv_mean <- score_thresh_perf_dt[,mean(npv)]
        dt$fpr_lwr <- score_thresh_perf_dt[,quantile(fpr,c(0.025),na.rm=T)] %>% unname
        dt$fpr_upr <- score_thresh_perf_dt[,quantile(fpr,c(0.975),na.rm=T)] %>% unname
        dt$fpr_mean <- score_thresh_perf_dt[,mean(fpr)]
        dt$fnr_lwr <- score_thresh_perf_dt[,quantile(fnr,c(0.025),na.rm=T)] %>% unname
        dt$fnr_upr <- score_thresh_perf_dt[,quantile(fnr,c(0.975),na.rm=T)] %>% unname
        dt$fnr_mean <- score_thresh_perf_dt[,mean(fnr)]
        dt$pos_class_score_lwr = score_thresh_perf_dt[,quantile(pos_class_mean,c(0.025),na.rm=T)]
        dt$pos_class_score_mean = score_thresh_perf_dt[,mean(pos_class_mean)]
        dt$pos_class_score_upr = score_thresh_perf_dt[,quantile(pos_class_mean,c(0.975),na.rm=T)]
        dt$neg_class_score_lwr = score_thresh_perf_dt[,quantile(neg_class_mean,c(0.025),na.rm=T)]
        dt$neg_class_score_mean = score_thresh_perf_dt[,mean(neg_class_mean)]
        dt$neg_class_score_upr = score_thresh_perf_dt[,quantile(neg_class_mean,c(0.975),na.rm=T)]
        dt
    }
  }
  dts
}

tmp <- bind_rows(
  generate_filtered_perf_curve(method="PRR",score="PRR",type="score",x="tpr",y="fpr",thresh=1),
  generate_filtered_perf_curve(method="PRR",score="PRR_90mse",type="90mse",x="tpr",y="fpr",thresh=1),
  generate_filtered_perf_curve(method="GAM",score="gam_score",type="score",x="tpr",y="fpr",thresh=0),
  generate_filtered_perf_curve(method="GAM",score="gam_score_90mse",type="90mse",x="tpr",y="fpr",thresh=0),
  generate_filtered_perf_curve(method="LowOrder2",score="poly2_score",type="score",x="tpr",y="fpr",thresh=0),
  generate_filtered_perf_curve(method="LowOrder2",score="poly2_score_90mse",type="90mse",x="tpr",y="fpr",thresh=0),
  generate_filtered_perf_curve(method="LowOrder3",score="poly3_score",type="score",x="tpr",y="fpr",thresh=0),
  generate_filtered_perf_curve(method="LowOrder3",score="poly3_score_90mse",type="90mse",x="tpr",y="fpr",thresh=0),
  generate_filtered_perf_curve(method="LowOrder4",score="poly4_score",type="score",x="tpr",y="fpr",thresh=0),
  generate_filtered_perf_curve(method="LowOrder4",score="poly4_score_90mse",type="90mse",x="tpr",y="fpr",thresh=0)
)

g <- tmp %>% 
  ggplot(aes(fpr,tpr,color=method)) +
  geom_line() +
  facet_grid(type~spikein) +
  geom_abline(intercept=0,slope=1,linetype=2,color="red") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  xlim(0,1) +
  ylim(0,1)
ggsave(paste0(img_dir,basename,"detection_withlowordermodels_before_performance_roc_curves.png"),g,width=12,height=5)

tmp <- bind_rows(
  generate_filtered_perf_curve(method="PRR",score="PRR",type="score",x="tpr",y="fpr",thresh=1),
  generate_filtered_perf_curve(method="PRR",score="PRR_90mse",type="90mse",x="tpr",y="fpr",thresh=1),
  generate_filtered_perf_curve(method="GAM",score="gam_score",type="score",x="tpr",y="fpr",thresh=0),
  generate_filtered_perf_curve(method="GAM",score="gam_score_90mse",type="90mse",x="tpr",y="fpr",thresh=0)
)

g <- tmp %>% 
  ggplot(aes(fpr,tpr,color=method)) +
  geom_line() +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein) +
  geom_abline(intercept=0,slope=1,linetype=2,color="red") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  xlim(0,1) +
  ylim(0,1)
ggsave(paste0(img_dir,basename,"detection_before_performance_roc_curves.png"),g,width=12,height=5)


g <- tmp %>% 
  ggplot(aes(fpr,tpr,color=method)) +
  geom_line() +
  facet_grid(type~spikein) +
  geom_abline(intercept=0,slope=1,linetype=2,color="red") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  xlim(0,1) +
  ylim(0,1)
ggsave(paste0(img_dir,basename,"detection_withlowordermodels_before_performance_roc_curves.png"),g,width=12,height=5)

tmp[,.(auc = simple_auc(tpr,fpr)),.(method,type,spikein)] %>% dcast(method + type ~ spikein)

tmp %>% 
  .[fpr<=0.25 & spikein=="decrease"] %>% 
  ggplot(aes(fpr,tpr,color=method)) +
  geom_line() +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein)

tmp[fpr<=.25,.(auc = simple_auc(tpr,fpr)),.(method,type,spikein)] %>% dcast(method + type ~ spikein)

tmp2 <- 
  tmp[,.(spikein,method,type,
         pos_class_score_lwr,
         pos_class_score_mean,
         pos_class_score_upr,
         neg_class_score_lwr,
         neg_class_score_mean,
         neg_class_score_upr)] %>% 
  unique() %>% 
  pivot_longer(cols=c("pos_class_score_lwr",
                      "pos_class_score_mean",
                      "pos_class_score_upr",
                      "neg_class_score_lwr",
                      "neg_class_score_mean",
                      "neg_class_score_upr"))
tmp2$name2 <- str_split_fixed(tmp2$name,"_",2)[,2]
tmp2$class <- str_split_fixed(tmp2$name,"_",2)[,1]
tmp2$class <- ifelse(tmp2$class=="pos","positive control","negative control")
tmp2$summary <- str_split_fixed(tmp2$name2,"_",3)[,3]
g <- 
  tmp2 %>% 
  data.table() %>% 
  dcast(spikein + method + type + class ~ summary,value.var = "value") %>% 
  ggplot(aes(spikein,mean,color=class)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),
                width=0.1,position = position_dodge(width=.5)) +
  scale_color_brewer(palette="Dark2") +
  facet_grid(method~type,scales="free") +
  xlab("") +
  ylab("Bootstrap score in class")
ggsave(paste0(img_dir,basename,"detection_before_performance_bootstrap_score_by_control.png"),g,width=15,height=7)

g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,tpr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=tpr_lwr,ymax=tpr_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("TPR")
ggsave(paste0(img_dir,basename,"detection_before_performance_tpr_at_score_threshold.png"),g,width=10,height=7)

g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,fpr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=fpr_lwr,ymax=fpr_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("FPR")
ggsave(paste0(img_dir,basename,"detection_before_performance_fpr_at_score_threshold.png"),g,width=10,height=7)


g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,fnr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=fnr_lwr,ymax=fnr_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("FNR")
ggsave(paste0(img_dir,basename,"detection_before_performance_fnr_at_score_threshold.png"),g,width=10,height=7)

g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,tnr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=tnr_lwr,ymax=tnr_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("TNR")
ggsave(paste0(img_dir,basename,"detection_before_performance_tnr_at_score_threshold.png"),g,width=10,height=7)

g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,ppv_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=ppv_lwr,ymax=ppv_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("PPV\nTrue positives out of all predicted positive")
ggsave(paste0(img_dir,basename,"detection_before_performance_ppv_at_score_threshold.png"),g,width=10,height=7)

g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,npv_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=npv_lwr,ymax=npv_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("NPV\nTrue negatives out of all predicted negative")
ggsave(paste0(img_dir,basename,"detection_before_performance_npv_at_score_threshold.png"),g,width=10,height=7)

g <- tmp[,.(spikein,method,type,auc_lwr,auc,auc_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,auc,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=auc_lwr,ymax=auc_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("AUROC")
ggsave(paste0(img_dir,basename,"detection_before_performance_auc.png"),g,width=10,height=7)

tmp[,.(spikein,method,type,auc_lwr)] %>% 
  unique() %>% 
  dcast(method + type ~ spikein) %>% .[order(method,type)]


tmp$fpr_cut <- cut(tmp$fpr,breaks=seq(0,1,0.05))
g <- tmp %>% 
  na.omit %>% 
  .[,
    .(
      tpr_lwr = quantile(tpr,c(0.025),na.rm=T),
      tpr_mean = mean(tpr),
      tpr_upr = quantile(tpr,c(0.975),na.rm=T)
    ),
    .(method,type,fpr_cut,spikein)] %>% 
  ggplot(aes(fpr_cut,tpr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=tpr_lwr,ymax=tpr_upr),width=0.3,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein) +
  xlab("False Positive Rate categories") +
  ylab("Sensitivity") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_before_performance_auc_by_fpr.png"),g,width=21,height=7)

tmp <- bind_rows(
  generate_filtered_perf_curve(method="PRR",score="PRR",type="score",x="prec",y="rec",stat="aucpr",thresh=1),
  generate_filtered_perf_curve(method="PRR",score="PRR_90mse",type="90mse",x="prec",y="rec",stat="aucpr",thresh=1),
  generate_filtered_perf_curve(method="GAM",score="gam_score",type="score",x="prec",y="rec",stat="aucpr",thresh=1),
  generate_filtered_perf_curve(method="GAM",score="gam_score_90mse",type="90mse",x="prec",y="rec",stat="aucpr",thresh=1)
)

g <- tmp[,.(spikein,method,type,auc_lwr,auc,auc_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,auc,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=auc_lwr,ymax=auc_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("AUPRC")
ggsave(paste0(img_dir,basename,"detection_before_performance_auprc.png"),g,width=10,height=7)

g <- tmp %>% 
  ggplot(aes(rec,prec,color=method)) +
  geom_line() +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein) +
  xlab("Recall") +
  ylab("Precision") 
ggsave(paste0(img_dir,basename,"detection_before_performance_pr_curves.png"),g,width=12,height=5)

# ADE dynamic detection power analysis ---------------------------------

stage_score_power <- function(dat=stage_positive_data,method="PRR",score = "PRR",score_threshold=1,type='overall',high=T){
  
  if(high){rates <- c(1)}else{rates <- c(0,1)}
  
  esize='Tij'
  tsize="D"
  dat = merge(stage_spikein_class[class %in% rates],
              dat[,c(stage_col,"spikein",score,"ade",esize,tsize),with=F],
              by=c(stage_col,"spikein"),
              all.x=T) %>% 
    .[is.finite(get(score))] %>% 
    na.omit()
  
  es <- dat[is.finite(get(esize)),get(esize)] %>% unique()
  es <- c(es,seq(min(es),max(es),length=10))
  es <- es[order(es)]
  ts <- dat[is.finite(get(tsize)),get(tsize)] %>% unique()
  ts <- c(ts,seq(min(ts),max(ts),length=10))
  ts <- ts[order(ts)]
  es=quantile(es,probs=seq(0,1,0.1)) %>% unname
  ts=quantile(ts,probs=seq(0,1,0.1)) %>% unname
  ts_es <- expand.grid(seq(1,1000,length.out = 50),seq(0,0.2,length.out = 50)) %>% data.table()
  colnames(ts_es) <- c("tsize","esize")
  tmp <- foreach(i=1:nrow(ts_es),.combine = "rbind") %dopar% {
    esthresh <- ts_es[i][,esize]
    tsthresh <- ts_es[i][,tsize]
    sub <- dat[is.finite(get(score)) & get(esize)>=esthresh & get(tsize)>=tsthresh]
    pred_true <- sub[,get(score)]>score_threshold
    pred_false <- !pred_true
    cond_true <- sub[,class==1]
    c(
      ((as.integer(pred_true)+as.integer(cond_true))==2) %>% which() %>% length(),
      ((as.integer(pred_false)+as.integer(cond_true))==2) %>% which() %>% length(),
      length(pred_true),
      esthresh,
      tsthresh
    )
  } %>% data.table()
  colnames(tmp) <- c("tp","fn","len","Tij","D")
  ts_es_perf <- tmp
  ts_es_perf[,
             c("tpr","fnr","power") :=
               list(
                 tp/(tp+fn),
                 fn/(tp+fn),
                 tp/(tp+fn)
               )
  ]
  ts_es_perf$type <- type
  ts_es_perf$score <- score
  ts_es_perf$method <-  method
  ts_es_perf 
}

tmp <- 
  bind_rows(
    stage_score_power(method="PRR",score="PRR",type="score",score_threshold=1),
    stage_score_power(method="PRR",score="PRR_90mse",type="90mse",score_threshold=1),
    stage_score_power(method="GAM",score="gam_score",type="score",score_threshold=0),
    stage_score_power(method="GAM",score="gam_score_90mse",type="90mse",score_threshold=0)
  )

tmp %>% 
  fwrite(paste0(data_dir,basename,"power_analysis_results.csv"))

g <- tmp %>% 
  ggplot(aes(Tij,D,fill=tpr)) +
  geom_tile() +
  scale_fill_continuous(type="viridis",breaks=seq(0,1,0.1)) +
  guides(fill=guide_legend(title="Power")) +
  facet_grid(method~type) +
  xlab("Effect size") +
  ylab("Number of drug reports")
ggsave(paste0(img_dir,basename,"power_analysis_results.png"),g,width=10,height=10)

tmp <- fread(paste0(data_dir,basename,"power_analysis_results.csv"))

tab <- tmp[tpr>=0.8][order(Tij,D,method)][,.SD,.(score,method,type),.SDcols=c("Tij","D")]

powered_score_ades_prr <- stage_positive_data[nichd!='all'] %>% 
  .[spikein!="original" & spikein!="uniform"] %>% 
  pivot_longer(cols=c("PRR","PRR_90mse"),names_to="score") %>% 
  data.table() %>% 
  .[,.(spikein,nichd,ade,score,value,D,Tij)] %>% 
  na.omit() %>% 
  merge(
    tab,
    by=c("score"),
    allow.cartesian = T
  ) %>% 
  .[Tij.x>=Tij.y & D.x>=D.y,.(score,spikein,nichd,ade,Tij.x,D.x)] %>% unique()

powered_score_ades_gam <- stage_positive_data[nichd!='all'] %>% 
  .[spikein!="original" & spikein!="uniform"] %>% 
  pivot_longer(cols=c("gam_score","gam_score_90mse"),names_to="score") %>% 
  data.table() %>% 
  .[,.(spikein,nichd,ade,score,value,D,Tij)] %>% 
  na.omit() %>% 
  merge(
    tab,
    by=c("score"),
    allow.cartesian = T
  ) %>% 
  .[Tij.x>=Tij.y & D.x>=D.y,.(score,spikein,nichd,ade,Tij.x,D.x)] %>% unique()

powered_score_ades <- 
  bind_rows(
    powered_score_ades_prr,
    powered_score_ades_gam
  ) 

powered_score_ades %>% 
  fwrite(paste0(data_dir,basename,"power_analysis_powered_ades.csv"))

powered_score_ades[,.(score,spikein,ade)] %>% unique() %>% dcast(score ~ spikein)

score_ades_gam <- stage_positive_data[nichd!='all'] %>% 
  .[spikein!="original" & spikein!="uniform"] %>% 
  pivot_longer(cols=c("gam_score","gam_score_90mse"),names_to="score") %>% 
  data.table() %>% 
  .[,.(spikein,nichd,ade,score,value,D,Tij)] %>% 
  na.omit() %>% 
  merge(
    tmp[
      order(Tij,D,method)
    ][,
      .SD,
      .(score,method,type),
      .SDcols=c("Tij","D")
    ] ,
    by=c("score"),
    allow.cartesian = T
  ) %>% 
  .[Tij.x>=Tij.y & D.x>=D.y,.(score,spikein,nichd,ade,Tij.x,D.x)] %>% unique()

score_ades_prr <- stage_positive_data[nichd!='all'] %>% 
  .[spikein!="original" & spikein!="uniform"] %>% 
  pivot_longer(cols=c("PRR","PRR_90mse"),names_to="score") %>% 
  data.table() %>% 
  .[,.(spikein,nichd,ade,score,value,D,Tij)] %>% 
  na.omit() %>% 
  merge(
    tmp[
      order(Tij,D,method)
    ][,
      .SD,
      .(score,method,type),
      .SDcols=c("Tij","D")
    ] ,
    by=c("score"),
    allow.cartesian = T
  ) %>% 
  .[Tij.x>=Tij.y & D.x>=D.y,.(score,spikein,nichd,ade,Tij.x,D.x)] %>% unique()

score_ades <- 
  bind_rows(
    score_ades_prr,
    score_ades_gam
  ) 

score_ades %>% 
  fwrite(paste0(data_dir,basename,"power_analysis_all_ades.csv"))

score_ades[,.(score,spikein,ade)] %>% unique() %>% dcast(score ~ spikein)


# Powered reference ade score comparison ------------------------------

powered_score_ade <- fread(paste0(data_dir,basename,"power_analysis_powered_ades.csv"))

ades <- 
  powered_score_ade[
    score=="gam_score_90mse" & 
      spikein=="decrease" &
      Tij.x>0.05,unique(ade)]

for(pair in ades){
  g <- 
    stage_positive_data[spikein=="decrease" & ade==pair] %>% 
    ggplot(aes(factor(nichd,levels=category_levels[[stage_col]]))) +
    geom_point(aes(y=gam_score),
               position = position_nudge(x = -0.2),color=score_colors[["GAM"]],size=2) +
    geom_errorbar(aes(ymin=gam_score_90mse,ymax=gam_score_90pse),
                  position = position_nudge(x = -0.2),color=score_colors[["GAM"]],width=0.1,size=1) +
    geom_point(aes(y=log10(PRR)),
               position = position_nudge(x = 0.2),color=score_colors[["PRR"]],size=2) +
    geom_errorbar(aes(ymin=log10(PRR_90mse),ymax=log10(PRR_90pse)),
                  position = position_nudge(x = 0.2),color=score_colors[["PRR"]],width=0.1,size=1) +
    scale_y_continuous(sec.axis=sec_axis(~.,name="")) +
    guides(color=guide_legend(title="Detection method",title.position = "top")) +
    xlab("") +
    ylab("") +
    ggtitle("") +
    theme(
      #strip.background = element_blank(),
      #strip.text = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(angle=45,hjust=1,vjust=1)
    )
  ggsave(paste0(img_dir,basename,"decrease_class_powered_ade_",pair,"_prr_vs_gam.png"),g,
         width = 6,height=4)
  
}


# ADE dynamic detection performance: with power filter ------------------------------------------------

powered_score_ade <- fread(paste0(data_dir,basename,"power_analysis_powered_ades.csv"))
setnames(powered_score_ade,"score","score_")

understanding_auc <- function(powered_score_ade,method="gam_score",score="gam_score",type="overall",N=1e6,seed=0){
  set.seed(seed)
  dts <- 
    foreach(spikein_col=stage_spikein_class[,unique(spikein)],.combine="rbind") %dopar% {
      p = 
        stage_positive_data[spikein==spikein_col] %>% 
        merge(powered_score_ade,by=c("ade","nichd","spikein"))
      neg=stage_negative_data
      y_true <- c(rep(1,nrow(p)),rep(0,nrow(neg)))
      y_pred <- c(p[,get(score)],neg[,get(score)])
      y_true <- y_true[is.finite(y_pred)]
      y_pred <- y_pred[is.finite(y_pred)]
      if(length(unique(y_true))==2){
        #https://www.r-bloggers.com/2016/11/calculating-auc-the-area-under-a-roc-curve/
        p <- sample(y_pred[as.logical(y_true)], N, replace=TRUE)
        n <- sample(y_pred[!as.logical(y_true)], N, replace=TRUE)
        kst = ks.test(n,p,alternative = "greater")
        tp = sum(p > n)
        tie = sum(p == n)/2
        fp = length(p) - tp
        dt <- 
          bind_rows(
            data.table(value=p,sample="positive control"),
            data.table(value=n,sample="negative control")
          )
        dt$ptp = tp
        dt$pfp = fp
        dt$pfpr = fp/(tp)
        dt$ptie = tie
        dt$pauc=(tp + tie) / N
        dt$spikein=spikein_col
        dt$method=method
        dt$score=score
        dt$type=type
        dt$test_statistic <- kst$statistic %>% unname
        dt$test_pvalue <- kst$p.value
        dt$test_method <- kst$method
        dt
    }
  }
  dts
}

tmp <- bind_rows(
  understanding_auc(powered_score_ade[score_=="PRR"],method="PRR",score="PRR",type="score"),
  understanding_auc(powered_score_ade[score_=="PRR_90mse"],method="PRR",score="PRR_90mse",type="90mse"),
  understanding_auc(powered_score_ade[score_=="gam_score"],method="GAM",score="gam_score",type="score"),
  understanding_auc(powered_score_ade[score_=="gam_score_90mse"],method="GAM",score="gam_score_90mse",type="90mse")
)

g <- tmp %>% 
  .[,.(spikein,method,type,sample,test_statistic)] %>% 
  unique() %>% 
  ggplot(aes(spikein,test_statistic,color=method)) +
  geom_point() +
  facet_grid(~type) +
  scale_color_manual(values=score_colors) +
  xlab("") +
  ylab("KS test statistic") 
ggsave(paste0(img_dir,basename,"detection_after_performance_kstest_control_distributions.png"),g,width=12,height=5)

g <- tmp %>% 
  .[type=="score"] %>% 
  ggplot(aes(value,fill=sample)) +
  geom_density(alpha=0.5) +
  scale_x_continuous(trans="log10",labels=scales::comma) +
  facet_grid(method~spikein) +
  guides(fill=guide_legend(title="Control")) +
  xlab("") +
  ylab("Score density") +
  ggtitle("Distribution of resampled scores used for calculating probabilistic AUC")
ggsave(paste0(img_dir,basename,"detection_after_performance_bootstrap_score_by_control_pauc.png"),g,width=15,height=7)

g <- tmp %>% 
  .[type=="90mse"] %>% 
  ggplot(aes(value,fill=sample)) +
  geom_density(alpha=0.5) +
  scale_x_continuous(trans="log10",labels=scales::comma) +
  facet_grid(method~spikein) +
  guides(fill=guide_legend(title="Control")) +
  xlab("") +
  ylab("Score density") +
  ggtitle("Distribution of resampled scores used for calculating probabilistic AUC")
ggsave(paste0(img_dir,basename,"detection_after_performance_bootstrap_90mse_score_by_control_pauc.png"),g,width=15,height=7)

g <- tmp %>% 
  ggplot(aes(spikein,value,color=sample)) +
  geom_boxplot(outlier.shape=NA) +
  scale_y_continuous(trans="log10",labels=scales::comma) +
  facet_grid(method~type,scales="free") +
  xlab("") +
  ylab("Score")
ggsave(paste0(img_dir,basename,"detection_after_performance_bootstrap_score_by_control_pauc_boxplot.png"),g,width=15,height=7)

g <- tmp %>% 
  ggplot(aes(spikein,value,color=sample)) +
  geom_boxplot(outlier.shape=NA) +
  geom_violin() +
  scale_y_continuous(trans="log10",labels=scales::comma) +
  facet_grid(method~type,scales="free") +
  xlab("") +
  ylab("Score")
ggsave(paste0(img_dir,basename,"detection_after_performance_bootstrap_score_by_control_pauc_boxplot_overlap.png"),g,width=15,height=7)
  
powered_score_ade <- fread(paste0(data_dir,basename,"power_analysis_powered_ades.csv"))
setnames(powered_score_ade,"score","score_")

generate_filtered_perf_curve <- function(method="gam_score",score="gam_score",type="overall",x="tpr",y="fpr",stat="auc",thresh=0,high=T){
  
  if(high){rates <- c(1)}else{rates <- c(0,1)}
  
  dts <- 
    foreach(spikein_col=stage_spikein_class[,unique(spikein)],.combine="rbind") %dopar% {
      pos = 
        stage_positive_data[spikein==spikein_col] %>% 
        merge(powered_score_ade[score_==score],by=c("ade","nichd","spikein")) %>% 
        merge(stage_spikein_class[class %in% rates],by=c("nichd","spikein"))
      neg=stage_negative_data
      y_true <- c(rep(1,nrow(pos)),rep(0,nrow(neg)))
      y_pred <- c(pos[,get(score)],neg[,get(score)])
      y_true <- y_true[is.finite(y_pred)]
      y_pred <- y_pred[is.finite(y_pred)]
      if(length(unique(y_true))==2){
        score_thresh_perf_dt = 
        lapply(1:100,
               function(i){
                 set.seed(i)
                 sinds = sample(1:length(y_pred),length(y_pred),replace=T)
                 pos_class_mean = mean(y_pred[sinds][y_true[sinds]==1])
                 neg_class_mean = mean(y_pred[sinds][y_true[sinds]==0])
                 cond_true <- y_true[sinds]==1
                 cond_false <- !cond_true
                 pred_true <- y_pred[sinds]>=thresh
                 pred_false <- !pred_true
                 tp <- ((as.integer(pred_true)+as.integer(cond_true))==2) %>% which() %>% length()
                 fn <- ((as.integer(pred_false)+as.integer(cond_true))==2) %>% which() %>% length()
                 fp <- ((as.integer(pred_true)+as.integer(cond_false))==2) %>% which() %>% length()
                 tn <- ((as.integer(pred_false)+as.integer(cond_false))==2) %>% which() %>% length()
                 pauc = auc_probability(as.logical(y_true),y_pred,N=1e5,split=T)
                 data.table(score_threshold=thresh,tp=tp,tn=tn,fp=fp,fn=fn,
                            tpr=(tp/(tp+fn)),fpr=(fp/(tn+fp)),tnr=(tn/(tn+fp)),fnr=(fn/(fn+tp)),
                            ppv=(tp/(fp+tp)),npv=(tn/(tn+fn)),
                            pauc = sum(pauc),ptp = pauc[1],ptie=pauc[2],
                            pos_class_mean = pos_class_mean,neg_class_mean=neg_class_mean)
               }) %>% 
          bind_rows()
        pred <- ROCR::prediction(y_pred,y_true)
        perf <- ROCR::performance(pred,x,y)
        ci = auc_ci(y_pred,y_true,stat=stat)
        dt <- data.table(perf@y.values[[1]],perf@x.values[[1]],perf@alpha.values[[1]])
        colnames(dt) <- c(x,y,"threshold")
        aucs_upto25fpr = 
          sapply(1:100,
                 function(i){
                   set.seed(i)
                   sinds = sample(1:nrow(dt),nrow(dt),replace=T)
                   dt[sinds[order(sinds)]][get(x)<=.25,simple_auc(get(x),get(y))]
                 })      
        dt$method <- method
        dt$type <- type
        dt$spikein <- spikein_col
        dt$N <- length(y_true)
        dt$auc_lwr <- ci[1]
        dt$auc <- ci[2]
        dt$auc_upr <- ci[3]
        dt$pauc_lwr <- score_thresh_perf_dt[,quantile(pauc,c(0.025))]
        dt$pauc <- score_thresh_perf_dt[,mean(pauc)]
        dt$pauc_upr <- score_thresh_perf_dt[,quantile(pauc,c(0.975))]
        dt$auc_lwr_upto25fpr = quantile(aucs_upto25fpr,c(0.025))
        dt$auc_upto25fpr = mean(aucs_upto25fpr)
        dt$auc_upr_upto25fpr = quantile(aucs_upto25fpr,c(0.975))
        dt$ptp_lwr <- score_thresh_perf_dt[,quantile(ptp,c(0.025))]
        dt$ptp <- score_thresh_perf_dt[,mean(ptp)]
        dt$ptp_upr <- score_thresh_perf_dt[,quantile(ptp,c(0.975))]
        dt$ptie_lwr <- score_thresh_perf_dt[,quantile(ptie,c(0.025))]
        dt$ptie <- score_thresh_perf_dt[,mean(ptie)]
        dt$ptie_upr <- score_thresh_perf_dt[,quantile(ptie,c(0.975))]
        dt$tpr_lwr <- score_thresh_perf_dt[,quantile(tpr,c(0.025),na.rm=T)] %>% unname
        dt$tpr_upr <- score_thresh_perf_dt[,quantile(tpr,c(0.975),na.rm=T)] %>% unname
        dt$tpr_mean <- score_thresh_perf_dt[,mean(tpr)]
        dt$tnr_lwr <- score_thresh_perf_dt[,quantile(tnr,c(0.025),na.rm=T)] %>% unname
        dt$tnr_upr <- score_thresh_perf_dt[,quantile(tnr,c(0.975),na.rm=T)] %>% unname
        dt$tnr_mean <- score_thresh_perf_dt[,mean(tnr)]
        dt$ppv_lwr <- score_thresh_perf_dt[,quantile(ppv,c(0.025),na.rm=T)] %>% unname
        dt$ppv_upr <- score_thresh_perf_dt[,quantile(ppv,c(0.975),na.rm=T)] %>% unname
        dt$ppv_mean <- score_thresh_perf_dt[,mean(ppv)]
        dt$npv_lwr <- score_thresh_perf_dt[,quantile(npv,c(0.025),na.rm=T)] %>% unname
        dt$npv_upr <- score_thresh_perf_dt[,quantile(npv,c(0.975),na.rm=T)] %>% unname
        dt$npv_mean <- score_thresh_perf_dt[,mean(npv)]
        dt$fpr_lwr <- score_thresh_perf_dt[,quantile(fpr,c(0.025),na.rm=T)] %>% unname
        dt$fpr_upr <- score_thresh_perf_dt[,quantile(fpr,c(0.975),na.rm=T)] %>% unname
        dt$fpr_mean <- score_thresh_perf_dt[,mean(fpr)]
        dt$fnr_lwr <- score_thresh_perf_dt[,quantile(fnr,c(0.025),na.rm=T)] %>% unname
        dt$fnr_upr <- score_thresh_perf_dt[,quantile(fnr,c(0.975),na.rm=T)] %>% unname
        dt$fnr_mean <- score_thresh_perf_dt[,mean(fnr)]
        dt$pos_class_score_lwr = score_thresh_perf_dt[,quantile(pos_class_mean,c(0.025),na.rm=T)]
        dt$pos_class_score_mean = score_thresh_perf_dt[,mean(pos_class_mean)]
        dt$pos_class_score_upr = score_thresh_perf_dt[,quantile(pos_class_mean,c(0.975),na.rm=T)]
        dt$neg_class_score_lwr = score_thresh_perf_dt[,quantile(neg_class_mean,c(0.025),na.rm=T)]
        dt$neg_class_score_mean = score_thresh_perf_dt[,mean(neg_class_mean)]
        dt$neg_class_score_upr = score_thresh_perf_dt[,quantile(neg_class_mean,c(0.975),na.rm=T)]
        dt
      }
  }
  dts
}

tmp <- bind_rows(
  generate_filtered_perf_curve(method="PRR",score="PRR",type="score",thresh=1,high=T),
  generate_filtered_perf_curve(method="PRR",score="PRR_90mse",type="90mse",thresh=1,high=T),
  generate_filtered_perf_curve(method="GAM",score="gam_score",type="score",thresh=0,high=T),
  generate_filtered_perf_curve(method="GAM",score="gam_score_90mse",type="90mse",thresh=0,high=T)
)


g <- tmp %>% 
  ggplot(aes(fpr,tpr,color=method)) +
  geom_line() +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein) +
  geom_abline(intercept=0,slope=1,linetype=2,color="red") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  xlim(0,1) +
  ylim(0,1)
ggsave(paste0(img_dir,basename,"detection_after_performance_roc_curves.png"),g,width=12,height=5)

g <- tmp[,.(method,type,spikein,auc,auc_upto25fpr)] %>% 
  unique() %>% 
  pivot_longer(cols=c("auc","auc_upto25fpr")) %>% 
  ggplot(aes(spikein,value,color=method)) +
  geom_point() +
  geom_line(aes(group=method)) +
  scale_color_manual(values=score_colors) +
  facet_grid(name~type,scales="free") +
  xlab("") +
  ylab("AUROC")
ggsave(paste0(img_dir,basename,"detection_after_performance_auc_to_auc25percentfpr.png"),g,width=12,height=5)

tmp[,.(auc = simple_auc(tpr,fpr)),.(method,type,spikein)] %>% dcast(method + type ~ spikein)

tmp %>% 
  .[fpr<=0.25] %>% 
  ggplot(aes(fpr,tpr,color=method)) +
  geom_line() +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein)

tmp[fpr<=.25,.(auc = simple_auc(tpr,fpr)),.(method,type,spikein)] %>% dcast(method + type ~ spikein)

tmp2 <- 
  tmp[,.(spikein,method,type,
         pos_class_score_lwr,
         pos_class_score_mean,
         pos_class_score_upr,
         neg_class_score_lwr,
         neg_class_score_mean,
         neg_class_score_upr)] %>% 
  unique() %>% 
  pivot_longer(cols=c("pos_class_score_lwr",
                      "pos_class_score_mean",
                      "pos_class_score_upr",
                      "neg_class_score_lwr",
                      "neg_class_score_mean",
                      "neg_class_score_upr"))
tmp2$name2 <- str_split_fixed(tmp2$name,"_",2)[,2]
tmp2$class <- str_split_fixed(tmp2$name,"_",2)[,1]
tmp2$class <- ifelse(tmp2$class=="pos","positive control","negative control")
tmp2$summary <- str_split_fixed(tmp2$name2,"_",3)[,3]
g <- 
tmp2 %>% 
  data.table() %>% 
  dcast(spikein + method + type + class ~ summary,value.var = "value") %>% 
  ggplot(aes(spikein,mean,color=class)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),
                width=0.1,position = position_dodge(width=.5)) +
  scale_color_brewer(palette="Dark2") +
  facet_grid(method~type,scales="free") +
  xlab("") +
  ylab("Bootstrap score in class")
ggsave(paste0(img_dir,basename,"detection_after_performance_bootstrap_score_by_control.png"),g,width=15,height=7)

g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,tpr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=tpr_lwr,ymax=tpr_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("TPR")
ggsave(paste0(img_dir,basename,"detection_after_performance_tpr_at_score_threshold.png"),g,width=10,height=7)

g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,fpr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=fpr_lwr,ymax=fpr_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("FPR")
ggsave(paste0(img_dir,basename,"detection_after_performance_fpr_at_score_threshold.png"),g,width=10,height=7)


g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,fnr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=fnr_lwr,ymax=fnr_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("FNR")
ggsave(paste0(img_dir,basename,"detection_after_performance_fnr_at_score_threshold.png"),g,width=10,height=7)

g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,tnr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=tnr_lwr,ymax=tnr_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("TNR")
ggsave(paste0(img_dir,basename,"detection_after_performance_tnr_at_score_threshold.png"),g,width=10,height=7)

g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,ppv_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=ppv_lwr,ymax=ppv_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("PPV\nTrue positives out of all predicted positive")
ggsave(paste0(img_dir,basename,"detection_after_performance_ppv_at_score_threshold.png"),g,width=10,height=7)

g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,npv_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=npv_lwr,ymax=npv_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("NPV\nTrue negatives out of all predicted negative")
ggsave(paste0(img_dir,basename,"detection_after_performance_npv_at_score_threshold.png"),g,width=10,height=7)

g <- tmp[,.(spikein,method,type,auc_lwr,auc,auc_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,auc,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=auc_lwr,ymax=auc_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("AUROC")
ggsave(paste0(img_dir,basename,"detection_after_performance_auc.png"),g,width=10,height=7)

tmp[,.(spikein,method,type,auc_lwr)] %>% 
  unique() %>% 
  dcast(method + type ~ spikein) %>% .[order(method,type)]


tmp$fpr_cut <- cut(tmp$fpr,breaks=seq(0,1,0.05))
g <- tmp %>% 
  na.omit %>% 
  .[,
    .(
      tpr_lwr = quantile(tpr,c(0.025),na.rm=T),
      tpr_mean = mean(tpr),
      tpr_upr = quantile(tpr,c(0.975),na.rm=T)
    ),
    .(method,type,fpr_cut,spikein)] %>% 
  ggplot(aes(fpr_cut,tpr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=tpr_lwr,ymax=tpr_upr),width=0.3,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein) +
  xlab("False Positive Rate categories") +
  ylab("Sensitivity") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_performance_auc_by_fpr.png"),g,width=21,height=7)

tmp <- bind_rows(
  generate_filtered_perf_curve(pos,method="PRR",score="PRR",type="score",x="prec",y="rec",stat="aucpr",thresh=1),
  generate_filtered_perf_curve(pos,method="PRR",score="PRR_90mse",type="90mse",x="prec",y="rec",stat="aucpr",thresh=1),
  generate_filtered_perf_curve(pos,method="GAM",score="gam_score",type="score",x="prec",y="rec",stat="aucpr",thresh=1),
  generate_filtered_perf_curve(pos,method="GAM",score="gam_score_90mse",type="90mse",x="prec",y="rec",stat="aucpr",thresh=1)
)

g <- tmp[,.(spikein,method,type,auc_lwr,auc,auc_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,auc,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=auc_lwr,ymax=auc_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("AUPRC")
ggsave(paste0(img_dir,basename,"detection_after_performance_auprc.png"),g,width=10,height=7)

g <- tmp %>% 
  ggplot(aes(rec,prec,color=method)) +
  geom_line() +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein) +
  xlab("Recall") +
  ylab("Precision") 
ggsave(paste0(img_dir,basename,"detection_after_performance_pr_curves.png"),g,width=12,height=5)


generate_filtered_perf_curve_wn_stage <- function(powered_score_ade,method="gam_score",score="gam_score",type="overall",x="prec",y="rec",stat="auc",thresh=0,high=T){
  
  if(high){rates <- c(1)}else{rates <- c(0,1)}
  
  dts <- 
    foreach(spikein_col=stage_spikein_class[,unique(spikein)],.combine="rbind") %dopar% {
      dts <- NULL
      for(st in stage_spikein_class[,unique(get(stage_col))]){
        p = 
          stage_positive_data[nichd==st & spikein==spikein_col] %>% 
          merge(powered_score_ade,by=c("ade","nichd","spikein")) %>% 
          merge(stage_spikein_class[class %in% rates],by=c("nichd","spikein"))
        n = stage_negative_data[nichd==st]
        y_true <- c(rep(1,nrow(p)),rep(0,nrow(n)))
        y_pred <- c(p[,get(score)],n[,get(score)])
        y_true <- y_true[is.finite(y_pred)]
        y_pred <- y_pred[is.finite(y_pred)]
        if(length(unique(y_true))==2){
          score_thresh_perf_dt = 
            lapply(1:100,
                   function(i){
                     set.seed(i)
                     sinds = sample(1:length(y_pred),length(y_pred),replace=T)
                     pos_class_mean = mean(y_pred[sinds][y_true[sinds]==1])
                     neg_class_mean = mean(y_pred[sinds][y_true[sinds]==0])
                     cond_true <- y_true[sinds]==1
                     cond_false <- !cond_true
                     pred_true <- y_pred[sinds]>=thresh
                     pred_false <- !pred_true
                     tp <- ((as.integer(pred_true)+as.integer(cond_true))==2) %>% which() %>% length()
                     fn <- ((as.integer(pred_false)+as.integer(cond_true))==2) %>% which() %>% length()
                     fp <- ((as.integer(pred_true)+as.integer(cond_false))==2) %>% which() %>% length()
                     tn <- ((as.integer(pred_false)+as.integer(cond_false))==2) %>% which() %>% length()
                     pauc = auc_probability(as.logical(y_true),y_pred,N=1e5,split=T)
                     data.table(score_threshold=thresh,tp=tp,tn=tn,fp=fp,fn=fn,
                                tpr=(tp/(tp+fn)),fpr=(fp/(tn+fp)),tnr=(tn/(tn+fp)),fnr=(fn/(fn+tp)),
                                ppv=(fp/(fp+tp)),npv=(tn/(tn+fn)),
                                pauc = sum(pauc),ptp = pauc[1],ptie=pauc[2],
                                pos_class_mean = pos_class_mean,neg_class_mean=neg_class_mean)
                   }) %>% 
            bind_rows()
          pred <- ROCR::prediction(y_pred,y_true)
          perf <- ROCR::performance(pred,x,y)
          ci = auc_ci(y_pred,y_true,stat=stat)
          dt <- data.table(perf@y.values[[1]],perf@x.values[[1]],perf@alpha.values[[1]])
          colnames(dt) <- c(x,y,"threshold")
          aucs_upto25fpr = 
            sapply(1:100,
                 function(i){
                   set.seed(i)
                   sinds = sample(1:nrow(dt),nrow(dt),replace=T)
                   dt[sinds[order(sinds)]][fpr<=.25,simple_auc(fpr,tpr)]
                 })
          dt$method <- method
          dt$type <- type
          dt$spikein <- spikein_col
          dt$nichd <- st
          dt$N <- length(y_true)
          dt$auc_lwr <- ci[1]
          dt$auc <- ci[2]
          dt$auc_upr <- ci[3]
          dt$pauc_lwr <- score_thresh_perf_dt[,quantile(pauc,c(0.025))]
          dt$pauc <- score_thresh_perf_dt[,mean(pauc)]
          dt$pauc_upr <- score_thresh_perf_dt[,quantile(pauc,c(0.975))]
          dt$auc_lwr_upto25fpr = quantile(aucs_upto25fpr,c(0.025))
          dt$auc_upto25fpr = mean(aucs_upto25fpr)
          dt$auc_upr_upto25fpr = quantile(aucs_upto25fpr,c(0.975))
          dt$ptp_lwr <- score_thresh_perf_dt[,quantile(ptp,c(0.025))]
          dt$ptp <- score_thresh_perf_dt[,mean(ptp)]
          dt$ptp_upr <- score_thresh_perf_dt[,quantile(ptp,c(0.975))]
          dt$ptie_lwr <- score_thresh_perf_dt[,quantile(ptie,c(0.025))]
          dt$ptie <- score_thresh_perf_dt[,mean(ptie)]
          dt$ptie_upr <- score_thresh_perf_dt[,quantile(ptie,c(0.975))]
          dt$tpr_lwr <- score_thresh_perf_dt[,quantile(tpr,c(0.025),na.rm=T)] %>% unname
          dt$tpr_upr <- score_thresh_perf_dt[,quantile(tpr,c(0.975),na.rm=T)] %>% unname
          dt$tpr_mean <- score_thresh_perf_dt[,mean(tpr)]
          dt$tnr_lwr <- score_thresh_perf_dt[,quantile(tnr,c(0.025),na.rm=T)] %>% unname
          dt$tnr_upr <- score_thresh_perf_dt[,quantile(tnr,c(0.975),na.rm=T)] %>% unname
          dt$tnr_mean <- score_thresh_perf_dt[,mean(tnr)]
          dt$ppv_lwr <- score_thresh_perf_dt[,quantile(ppv,c(0.025),na.rm=T)] %>% unname
          dt$ppv_upr <- score_thresh_perf_dt[,quantile(ppv,c(0.975),na.rm=T)] %>% unname
          dt$ppv_mean <- score_thresh_perf_dt[,mean(ppv)]
          dt$npv_lwr <- score_thresh_perf_dt[,quantile(npv,c(0.025),na.rm=T)] %>% unname
          dt$npv_upr <- score_thresh_perf_dt[,quantile(npv,c(0.975),na.rm=T)] %>% unname
          dt$npv_mean <- score_thresh_perf_dt[,mean(npv)]
          dt$fpr_lwr <- score_thresh_perf_dt[,quantile(fpr,c(0.025),na.rm=T)] %>% unname
          dt$fpr_upr <- score_thresh_perf_dt[,quantile(fpr,c(0.975),na.rm=T)] %>% unname
          dt$fpr_mean <- score_thresh_perf_dt[,mean(fpr)]
          dt$fnr_lwr <- score_thresh_perf_dt[,quantile(fnr,c(0.025),na.rm=T)] %>% unname
          dt$fnr_upr <- score_thresh_perf_dt[,quantile(fnr,c(0.975),na.rm=T)] %>% unname
          dt$fnr_mean <- score_thresh_perf_dt[,mean(fnr)]
          dt$pos_class_score_lwr = score_thresh_perf_dt[,quantile(pos_class_mean,c(0.025),na.rm=T)]
          dt$pos_class_score_mean = score_thresh_perf_dt[,mean(pos_class_mean)]
          dt$pos_class_score_upr = score_thresh_perf_dt[,quantile(pos_class_mean,c(0.975),na.rm=T)]
          dt$neg_class_score_lwr = score_thresh_perf_dt[,quantile(neg_class_mean,c(0.025),na.rm=T)]
          dt$neg_class_score_mean = score_thresh_perf_dt[,mean(neg_class_mean)]
          dt$neg_class_score_upr = score_thresh_perf_dt[,quantile(neg_class_mean,c(0.975),na.rm=T)]
          dts <- bind_rows(dts,dt)
        }
      }
      dts
  }
  dts
}

tmp <- bind_rows(
  generate_filtered_perf_curve_wn_stage(powered_score_ade[score_=="PRR"],method="PRR",score="PRR",type="score",x="fpr",y="tpr",stat="auc",thresh=1),
  generate_filtered_perf_curve_wn_stage(powered_score_ade[score_=="PRR_90mse"],method="PRR",score="PRR_90mse",type="90mse",x="fpr",y="tpr",stat="auc",thresh=1),
  generate_filtered_perf_curve_wn_stage(powered_score_ade[score_=="gam_score"],method="GAM",score="gam_score",type="score",x="fpr",y="tpr",stat="auc",thresh=1),
  generate_filtered_perf_curve_wn_stage(powered_score_ade[score_=="gam_score_90mse"],method="GAM",score="gam_score_90mse",type="90mse",x="fpr",y="tpr",stat="auc",thresh=1)
)

g <- tmp[,.(method,type,nichd,spikein,auc_upto25fpr)] %>% 
  unique() %>% 
  ggplot(aes(factor(nichd,levels=category_levels$nichd),auc_upto25fpr,color=method)) +
  geom_point() +
  geom_line(aes(group=method)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein,scales="free") +
  xlab("") +
  ylab("AUROC") +
  theme(
    axis.text.x = element_text(angle=45,hjust=1,vjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_performance_auc25percentfpr_wn_nichd.png"),g,width=12,height=5)

g <- tmp[,
    .(method,type,spikein,nichd,auc_lwr,auc,auc_upr)
    ] %>% 
  unique() %>% 
  ggplot(aes(factor(nichd,levels=category_levels$nichd),auc,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=auc_lwr,ymax=auc_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein,scales="free") +
  xlab("") +
  ylab("AUROC") +
  theme(
    axis.text.x = element_text(angle=45,hjust=1,vjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_performance_auc_wn_nichd.png"),g,width=15,height=7)

# ADE dynamic detection performance: with power filter and intersecting powered ades ------------------------------------------------

powered_score_ade <- fread(paste0(data_dir,basename,"power_analysis_powered_ades.csv"))
setnames(powered_score_ade,"score","score_")

generate_filtered_perf_curve_inter <- function(method="gam_score",score="gam_score",type="score",x="tpr",y="fpr",stat="auc",thresh=0,high=T){
  
  if(high){rates <- c(1)}else{rates <- c(0,1)}
  
  dts <- 
    foreach(spikein_col=stage_spikein_class[,unique(spikein)],.combine="rbind") %dopar% {
    if(type=="90mse"){
      powered <- powered_score_ade[grepl("90mse",score_) &
                                     spikein==spikein_col] %>% 
        .[,.(score_,ade)] %>% 
        unique() %>% 
        .[,.N,ade] %>% 
        .[N==2,ade]
    }else{
      powered <- powered_score_ade[!grepl("90mse",score_) &
                                     spikein==spikein_col] %>% 
        .[,.(score_,ade)] %>% 
        unique() %>% 
        .[,.N,ade] %>% 
        .[N==2,ade]
    }
    pos = 
      stage_positive_data[spikein==spikein_col & ade %in% powered] %>% 
      merge(stage_spikein_class[class %in% rates],by=c("nichd","spikein"))
    neg=stage_negative_data
    y_true <- c(rep(1,nrow(pos)),rep(0,nrow(neg)))
    y_pred <- c(pos[,get(score)],neg[,get(score)])
    y_true <- y_true[is.finite(y_pred)]
    y_pred <- y_pred[is.finite(y_pred)]
    if(length(unique(y_true))==2){
      score_thresh_perf_dt = 
        lapply(1:100,
               function(i){
                 set.seed(i)
                 sinds = sample(1:length(y_pred),length(y_pred),replace=T)
                 pos_class_mean = mean(y_pred[sinds][y_true[sinds]==1])
                 neg_class_mean = mean(y_pred[sinds][y_true[sinds]==0])
                 cond_true <- y_true[sinds]==1
                 cond_false <- !cond_true
                 pred_true <- y_pred[sinds]>=thresh
                 pred_false <- !pred_true
                 tp <- ((as.integer(pred_true)+as.integer(cond_true))==2) %>% which() %>% length()
                 fn <- ((as.integer(pred_false)+as.integer(cond_true))==2) %>% which() %>% length()
                 fp <- ((as.integer(pred_true)+as.integer(cond_false))==2) %>% which() %>% length()
                 tn <- ((as.integer(pred_false)+as.integer(cond_false))==2) %>% which() %>% length()
                 pauc = auc_probability(as.logical(y_true),y_pred,N=1e5,split=T)
                 data.table(score_threshold=thresh,tp=tp,tn=tn,fp=fp,fn=fn,
                            tpr=(tp/(tp+fn)),fpr=(fp/(tn+fp)),tnr=(tn/(tn+fp)),fnr=(fn/(fn+tp)),
                            ppv=(tp/(fp+tp)),npv=(tn/(tn+fn)),
                            pauc = sum(pauc),ptp = pauc[1],ptie=pauc[2],
                            pos_class_mean = pos_class_mean,neg_class_mean=neg_class_mean)
               }) %>% 
        bind_rows()
      pred <- ROCR::prediction(y_pred,y_true)
      perf <- ROCR::performance(pred,x,y)
      ci = auc_ci(y_pred,y_true,stat=stat)
      dt <- data.table(perf@y.values[[1]],perf@x.values[[1]],perf@alpha.values[[1]])
      colnames(dt) <- c(x,y,"threshold")
      aucs_upto25fpr = 
        sapply(1:100,
               function(i){
                 set.seed(i)
                 sinds = sample(1:nrow(dt),nrow(dt),replace=T)
                 dt[sinds[order(sinds)]][get(x)<=.25,simple_auc(get(x),get(y))]
               })      
      dt$method <- method
      dt$type <- type
      dt$spikein <- spikein_col
      dt$N <- length(y_true)
      dt$auc_lwr <- ci[1]
      dt$auc <- ci[2]
      dt$auc_upr <- ci[3]
      dt$pauc_lwr <- score_thresh_perf_dt[,quantile(pauc,c(0.025))]
      dt$pauc <- score_thresh_perf_dt[,mean(pauc)]
      dt$pauc_upr <- score_thresh_perf_dt[,quantile(pauc,c(0.975))]
      dt$auc_lwr_upto25fpr = quantile(aucs_upto25fpr,c(0.025))
      dt$auc_upto25fpr = mean(aucs_upto25fpr)
      dt$auc_upr_upto25fpr = quantile(aucs_upto25fpr,c(0.975))
      dt$ptp_lwr <- score_thresh_perf_dt[,quantile(ptp,c(0.025))]
      dt$ptp <- score_thresh_perf_dt[,mean(ptp)]
      dt$ptp_upr <- score_thresh_perf_dt[,quantile(ptp,c(0.975))]
      dt$ptie_lwr <- score_thresh_perf_dt[,quantile(ptie,c(0.025))]
      dt$ptie <- score_thresh_perf_dt[,mean(ptie)]
      dt$ptie_upr <- score_thresh_perf_dt[,quantile(ptie,c(0.975))]
      dt$tpr_lwr <- score_thresh_perf_dt[,quantile(tpr,c(0.025),na.rm=T)] %>% unname
      dt$tpr_upr <- score_thresh_perf_dt[,quantile(tpr,c(0.975),na.rm=T)] %>% unname
      dt$tpr_mean <- score_thresh_perf_dt[,mean(tpr)]
      dt$tnr_lwr <- score_thresh_perf_dt[,quantile(tnr,c(0.025),na.rm=T)] %>% unname
      dt$tnr_upr <- score_thresh_perf_dt[,quantile(tnr,c(0.975),na.rm=T)] %>% unname
      dt$tnr_mean <- score_thresh_perf_dt[,mean(tnr)]
      dt$ppv_lwr <- score_thresh_perf_dt[,quantile(ppv,c(0.025),na.rm=T)] %>% unname
      dt$ppv_upr <- score_thresh_perf_dt[,quantile(ppv,c(0.975),na.rm=T)] %>% unname
      dt$ppv_mean <- score_thresh_perf_dt[,mean(ppv)]
      dt$npv_lwr <- score_thresh_perf_dt[,quantile(npv,c(0.025),na.rm=T)] %>% unname
      dt$npv_upr <- score_thresh_perf_dt[,quantile(npv,c(0.975),na.rm=T)] %>% unname
      dt$npv_mean <- score_thresh_perf_dt[,mean(npv)]
      dt$fpr_lwr <- score_thresh_perf_dt[,quantile(fpr,c(0.025),na.rm=T)] %>% unname
      dt$fpr_upr <- score_thresh_perf_dt[,quantile(fpr,c(0.975),na.rm=T)] %>% unname
      dt$fpr_mean <- score_thresh_perf_dt[,mean(fpr)]
      dt$fnr_lwr <- score_thresh_perf_dt[,quantile(fnr,c(0.025),na.rm=T)] %>% unname
      dt$fnr_upr <- score_thresh_perf_dt[,quantile(fnr,c(0.975),na.rm=T)] %>% unname
      dt$fnr_mean <- score_thresh_perf_dt[,mean(fnr)]
      dt$pos_class_score_lwr = score_thresh_perf_dt[,quantile(pos_class_mean,c(0.025),na.rm=T)]
      dt$pos_class_score_mean = score_thresh_perf_dt[,mean(pos_class_mean)]
      dt$pos_class_score_upr = score_thresh_perf_dt[,quantile(pos_class_mean,c(0.975),na.rm=T)]
      dt$neg_class_score_lwr = score_thresh_perf_dt[,quantile(neg_class_mean,c(0.025),na.rm=T)]
      dt$neg_class_score_mean = score_thresh_perf_dt[,mean(neg_class_mean)]
      dt$neg_class_score_upr = score_thresh_perf_dt[,quantile(neg_class_mean,c(0.975),na.rm=T)]
      dt
    }
  }
  dts
}

tmp <- bind_rows(
  generate_filtered_perf_curve_inter(method="PRR",score="PRR",type="score",thresh=1,high=T),
  generate_filtered_perf_curve_inter(method="PRR",score="PRR_90mse",type="90mse",thresh=1,high=T),
  generate_filtered_perf_curve_inter(method="GAM",score="gam_score",type="score",thresh=0,high=T),
  generate_filtered_perf_curve_inter(method="GAM",score="gam_score_90mse",type="90mse",thresh=0,high=T)
)

tmp %>% 
  fwrite(paste0(data_dir,basename,"ade_dynamic_detection_performance_results.csv"))

tmp <- 
  fread(paste0(data_dir,basename,"ade_dynamic_detection_performance_results.csv"))

g <- tmp[type=="score"] %>% 
  ggplot(aes(fpr,tpr,color=method)) +
  geom_line() +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein) +
  geom_abline(intercept=0,slope=1,linetype=2,color="red") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  xlim(0,1) +
  ylim(0,1)
ggsave(paste0(img_dir,basename,"detection_after_inter_ades_performance_roc_curves_score.png"),g,width=12,height=3)
g <- tmp[type=="90mse"] %>% 
  ggplot(aes(fpr,tpr,color=method)) +
  geom_line() +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein) +
  geom_abline(intercept=0,slope=1,linetype=2,color="red") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  xlim(0,1) +
  ylim(0,1)
ggsave(paste0(img_dir,basename,"detection_after_inter_ades_performance_roc_curves_90mse.png"),g,width=12,height=3)

tmp[,.(auc = simple_auc(tpr,fpr)),.(method,type,spikein)] %>% dcast(method + type ~ spikein)

tmp %>% 
  .[fpr<=0.25 & spikein=="decrease"] %>% 
  ggplot(aes(fpr,tpr,color=method)) +
  geom_line() +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein)

tmp[fpr<=.25,.(auc = simple_auc(tpr,fpr)),.(method,type,spikein)] %>% dcast(method + type ~ spikein)

tmp2 <- 
  tmp[,.(spikein,method,type,
         pos_class_score_lwr,
         pos_class_score_mean,
         pos_class_score_upr,
         neg_class_score_lwr,
         neg_class_score_mean,
         neg_class_score_upr)] %>% 
  unique() %>% 
  pivot_longer(cols=c("pos_class_score_lwr",
                      "pos_class_score_mean",
                      "pos_class_score_upr",
                      "neg_class_score_lwr",
                      "neg_class_score_mean",
                      "neg_class_score_upr"))
tmp2$name2 <- str_split_fixed(tmp2$name,"_",2)[,2]
tmp2$class <- str_split_fixed(tmp2$name,"_",2)[,1]
tmp2$class <- ifelse(tmp2$class=="pos","positive control","negative control")
tmp2$summary <- str_split_fixed(tmp2$name2,"_",3)[,3]
g <- 
  tmp2 %>% 
  data.table() %>% 
  dcast(spikein + method + type + class ~ summary,value.var = "value") %>% 
  ggplot(aes(spikein,mean,color=class)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),
                width=0.1,position = position_dodge(width=.5)) +
  scale_color_brewer(palette="Dark2") +
  facet_grid(method~type,scales="free") +
  xlab("") +
  ylab("Bootstrap score in class")
ggsave(paste0(img_dir,basename,"detection_after_inter_ades_performance_bootstrap_score_by_control.png"),g,width=15,height=7)

sub <- tmp[,.(spikein,method,type,
       tpr_lwr,tpr_mean,tpr_upr,
       fpr_lwr,fpr_mean,fpr_upr,
       tnr_lwr,tnr_mean,tnr_upr,
       fnr_lwr,fnr_mean,fnr_upr,
       ppv_lwr,ppv_mean,ppv_upr,
       npv_lwr,npv_mean,npv_upr,
       auc_lwr,auc_mean = auc,auc_upr)] %>% 
  unique() %>% 
  pivot_longer(
    cols=c("tpr_lwr","tpr_mean","tpr_upr",
           "fpr_lwr","fpr_mean","fpr_upr",
           "tnr_lwr","tnr_mean","tnr_upr",
           "fnr_lwr","fnr_mean","fnr_upr",
           "ppv_lwr","ppv_mean","ppv_upr",
           "npv_lwr","npv_mean","npv_upr",
           'auc_lwr',"auc_mean","auc_upr")
    ) %>% data.table()
sub$metric <- 
  sapply(str_split(sub$name,"_"),function(x){x[1]})
sub$statistic <- 
  sapply(str_split(sub$name,"_"),function(x){x[2]})
metric_names <- 
  list(
    "auc" = "AUROC",
    "tpr" = "Power",
    "ppv" = "Positive predictive value",
    "npv" = "Negative predictive value",
    "fpr" = "False positive rate",
    "tnr" = "True negative rate",
    "fnr" = "False negative rate"
    )
sub$metric_name <- 
  sapply(sub$metric,function(x){metric_names[[x]]}) %>% unname
g <- sub %>% 
  .[order(metric)] %>% 
  dcast(
    spikein + method + type + metric_name ~ statistic,
    value.var="value"
    ) %>% 
  .[type=="score" &
      metric_name %in% 
      sapply(c("auc","tpr","ppv","npv"),function(x){as.factor(metric_names[[x]])})
    ] %>% 
  ggplot(aes(spikein,mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_wrap(~factor(metric_name,levels= paste0(unname(metric_names))),scales="free_y") +
  xlab("") +
  ylab("Performance value") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_ades_performance_metrics_score.png"),g,width=7,height=6)
g <- sub %>% 
  .[order(metric)] %>% 
  dcast(
    spikein + method + type + metric_name ~ statistic,
    value.var="value"
  ) %>% 
  .[type=="90mse" &
      metric_name %in% 
      sapply(c("auc","tpr","ppv","npv"),function(x){as.factor(metric_names[[x]])})
  ] %>% 
  ggplot(aes(spikein,mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_wrap(~factor(metric_name,levels= paste0(unname(metric_names))),scales="free_y") +
  xlab("") +
  ylab("Performance value") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_ades_performance_metrics_90mse.png"),g,width=7,height=6)

g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,tpr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=tpr_lwr,ymax=tpr_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("Power") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_ades_performance_power_at_score_threshold.png"),g,width=4,height=5)

g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,fpr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=fpr_lwr,ymax=fpr_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("False positive rate") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_ades_performance_fpr_at_score_threshold.png"),g,width=6,height=5)


g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,fnr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=fnr_lwr,ymax=fnr_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("False negative rate") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_ades_performance_fnr_at_score_threshold.png"),g,width=6,height=5)

g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,tnr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=tnr_lwr,ymax=tnr_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("True negative rate") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_ades_performance_tnr_at_score_threshold.png"),g,width=6,height=5)

g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,ppv_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=ppv_lwr,ymax=ppv_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("Positive predictive value") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_ades_performance_ppv_at_score_threshold.png"),g,width=6,height=5)

g <- tmp[,.(spikein,method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,npv_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=npv_lwr,ymax=npv_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("Negative predictive value") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_ades_performance_npv_at_score_threshold.png"),g,width=6,height=5)

g <- tmp[,.(spikein,method,type,auc_lwr,auc,auc_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,auc,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=auc_lwr,ymax=auc_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("AUROC") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_ades_performance_auc.png"),g,width=6,height=5)

tmp[,.(spikein,method,type,auc_lwr)] %>% 
  unique() %>% 
  dcast(method + type ~ spikein) %>% .[order(method,type)]


tmp$fpr_cut <- cut(tmp$fpr,breaks=seq(0,1,0.05))
g <- tmp %>% 
  na.omit %>% 
  .[,
    .(
      tpr_lwr = quantile(tpr,c(0.025),na.rm=T),
      tpr_mean = mean(tpr),
      tpr_upr = quantile(tpr,c(0.975),na.rm=T)
    ),
    .(method,type,fpr_cut,spikein)] %>% 
  ggplot(aes(fpr_cut,tpr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=tpr_lwr,ymax=tpr_upr),width=0.3,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein) +
  xlab("False Positive Rate categories") +
  ylab("Sensitivity") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_ades_performance_auc_by_fpr.png"),g,width=21,height=7)

tmp <- bind_rows(
  generate_filtered_perf_curve_inter(method="PRR",score="PRR",type="score",x="prec",y="rec",stat="aucpr",thresh=1),
  generate_filtered_perf_curve_inter(method="PRR",score="PRR_90mse",type="90mse",x="prec",y="rec",stat="aucpr",thresh=1),
  generate_filtered_perf_curve_inter(method="GAM",score="gam_score",type="score",x="prec",y="rec",stat="aucpr",thresh=1),
  generate_filtered_perf_curve_inter(method="GAM",score="gam_score_90mse",type="90mse",x="prec",y="rec",stat="aucpr",thresh=1)
)

g <- tmp[,.(spikein,method,type,auc_lwr,auc,auc_upr)] %>% 
  unique() %>% 
  ggplot(aes(spikein,auc,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=auc_lwr,ymax=auc_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("AUPRC")
ggsave(paste0(img_dir,basename,"detection_after_inter_ades_performance_auprc.png"),g,width=10,height=7)

g <- tmp %>% 
  ggplot(aes(rec,prec,color=method)) +
  geom_line() +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein) +
  xlab("Recall") +
  ylab("Precision") 
ggsave(paste0(img_dir,basename,"detection_after_inter_ades_performance_pr_curves.png"),g,width=12,height=5)

generate_filtered_perf_curve_wn_stage_inter <- function(method="gam_score",score="gam_score",type="overall",x="prec",y="rec",stat="auc",thresh=0,high=T){
  
  if(high){rates <- c(1)}else{rates <- c(0,1)}
  
  dts <- 
    foreach(spikein_col=stage_spikein_class[,unique(spikein)],.combine="rbind") %dopar% {
      dts <- NULL
      for(st in stage_spikein_class[,unique(get(stage_col))]){
        if(type=="90mse"){
          powered <- powered_score_ade[grepl("90mse",score_) &
                                         spikein==spikein_col] %>% 
            .[,.(score_,ade)] %>% 
            unique() %>% 
            .[,.N,ade] %>% 
            .[N==2,ade]
        }else{
          powered <- powered_score_ade[!grepl("90mse",score_) &
                                         spikein==spikein_col] %>% 
            .[,.(score_,ade)] %>% 
            unique() %>% 
            .[,.N,ade] %>% 
            .[N==2,ade]
        }
        p = 
          stage_positive_data[nichd==st & ade %in% powered & spikein==spikein_col] %>% 
          merge(stage_spikein_class[class %in% rates],by=c("nichd","spikein"))
        n = stage_negative_data[nichd==st]
        y_true <- c(rep(1,nrow(p)),rep(0,nrow(n)))
        y_pred <- c(p[,get(score)],n[,get(score)])
        y_true <- y_true[is.finite(y_pred)]
        y_pred <- y_pred[is.finite(y_pred)]
        if(length(unique(y_true))==2){
          score_thresh_perf_dt = 
            lapply(1:100,
                   function(i){
                     set.seed(i)
                     sinds = sample(1:length(y_pred),length(y_pred),replace=T)
                     pos_class_mean = mean(y_pred[sinds][y_true[sinds]==1])
                     neg_class_mean = mean(y_pred[sinds][y_true[sinds]==0])
                     cond_true <- y_true[sinds]==1
                     cond_false <- !cond_true
                     pred_true <- y_pred[sinds]>=thresh
                     pred_false <- !pred_true
                     tp <- ((as.integer(pred_true)+as.integer(cond_true))==2) %>% which() %>% length()
                     fn <- ((as.integer(pred_false)+as.integer(cond_true))==2) %>% which() %>% length()
                     fp <- ((as.integer(pred_true)+as.integer(cond_false))==2) %>% which() %>% length()
                     tn <- ((as.integer(pred_false)+as.integer(cond_false))==2) %>% which() %>% length()
                     pauc = auc_probability(as.logical(y_true),y_pred,N=1e5,split=T)
                     data.table(score_threshold=thresh,tp=tp,tn=tn,fp=fp,fn=fn,
                                tpr=(tp/(tp+fn)),fpr=(fp/(tn+fp)),tnr=(tn/(tn+fp)),fnr=(fn/(fn+tp)),
                                ppv=(fp/(fp+tp)),npv=(tn/(tn+fn)),
                                pauc = sum(pauc),ptp = pauc[1],ptie=pauc[2],
                                pos_class_mean = pos_class_mean,neg_class_mean=neg_class_mean)
                   }) %>% 
            bind_rows()
          pred <- ROCR::prediction(y_pred,y_true)
          perf <- ROCR::performance(pred,x,y)
          ci = auc_ci(y_pred,y_true,stat=stat)
          dt <- data.table(perf@y.values[[1]],perf@x.values[[1]],perf@alpha.values[[1]])
          colnames(dt) <- c(x,y,"threshold")
          aucs_upto25fpr = 
            sapply(1:100,
                   function(i){
                     set.seed(i)
                     sinds = sample(1:nrow(dt),nrow(dt),replace=T)
                     dt[sinds[order(sinds)]][fpr<=.25,simple_auc(fpr,tpr)]
                   })
          dt$method <- method
          dt$type <- type
          dt$spikein <- spikein_col
          dt$nichd <- st
          dt$N <- length(y_true)
          dt$auc_lwr <- ci[1]
          dt$auc <- ci[2]
          dt$auc_upr <- ci[3]
          dt$pauc_lwr <- score_thresh_perf_dt[,quantile(pauc,c(0.025))]
          dt$pauc <- score_thresh_perf_dt[,mean(pauc)]
          dt$pauc_upr <- score_thresh_perf_dt[,quantile(pauc,c(0.975))]
          dt$auc_lwr_upto25fpr = quantile(aucs_upto25fpr,c(0.025))
          dt$auc_upto25fpr = mean(aucs_upto25fpr)
          dt$auc_upr_upto25fpr = quantile(aucs_upto25fpr,c(0.975))
          dt$ptp_lwr <- score_thresh_perf_dt[,quantile(ptp,c(0.025))]
          dt$ptp <- score_thresh_perf_dt[,mean(ptp)]
          dt$ptp_upr <- score_thresh_perf_dt[,quantile(ptp,c(0.975))]
          dt$ptie_lwr <- score_thresh_perf_dt[,quantile(ptie,c(0.025))]
          dt$ptie <- score_thresh_perf_dt[,mean(ptie)]
          dt$ptie_upr <- score_thresh_perf_dt[,quantile(ptie,c(0.975))]
          dt$tpr_lwr <- score_thresh_perf_dt[,quantile(tpr,c(0.025),na.rm=T)] %>% unname
          dt$tpr_upr <- score_thresh_perf_dt[,quantile(tpr,c(0.975),na.rm=T)] %>% unname
          dt$tpr_mean <- score_thresh_perf_dt[,mean(tpr)]
          dt$tnr_lwr <- score_thresh_perf_dt[,quantile(tnr,c(0.025),na.rm=T)] %>% unname
          dt$tnr_upr <- score_thresh_perf_dt[,quantile(tnr,c(0.975),na.rm=T)] %>% unname
          dt$tnr_mean <- score_thresh_perf_dt[,mean(tnr)]
          dt$ppv_lwr <- score_thresh_perf_dt[,quantile(ppv,c(0.025),na.rm=T)] %>% unname
          dt$ppv_upr <- score_thresh_perf_dt[,quantile(ppv,c(0.975),na.rm=T)] %>% unname
          dt$ppv_mean <- score_thresh_perf_dt[,mean(ppv)]
          dt$npv_lwr <- score_thresh_perf_dt[,quantile(npv,c(0.025),na.rm=T)] %>% unname
          dt$npv_upr <- score_thresh_perf_dt[,quantile(npv,c(0.975),na.rm=T)] %>% unname
          dt$npv_mean <- score_thresh_perf_dt[,mean(npv)]
          dt$fpr_lwr <- score_thresh_perf_dt[,quantile(fpr,c(0.025),na.rm=T)] %>% unname
          dt$fpr_upr <- score_thresh_perf_dt[,quantile(fpr,c(0.975),na.rm=T)] %>% unname
          dt$fpr_mean <- score_thresh_perf_dt[,mean(fpr)]
          dt$fnr_lwr <- score_thresh_perf_dt[,quantile(fnr,c(0.025),na.rm=T)] %>% unname
          dt$fnr_upr <- score_thresh_perf_dt[,quantile(fnr,c(0.975),na.rm=T)] %>% unname
          dt$fnr_mean <- score_thresh_perf_dt[,mean(fnr)]
          dt$pos_class_score_lwr = score_thresh_perf_dt[,quantile(pos_class_mean,c(0.025),na.rm=T)]
          dt$pos_class_score_mean = score_thresh_perf_dt[,mean(pos_class_mean)]
          dt$pos_class_score_upr = score_thresh_perf_dt[,quantile(pos_class_mean,c(0.975),na.rm=T)]
          dt$neg_class_score_lwr = score_thresh_perf_dt[,quantile(neg_class_mean,c(0.025),na.rm=T)]
          dt$neg_class_score_mean = score_thresh_perf_dt[,mean(neg_class_mean)]
          dt$neg_class_score_upr = score_thresh_perf_dt[,quantile(neg_class_mean,c(0.975),na.rm=T)]
          dts <- bind_rows(dts,dt)
        }
      }
      dts
    }
  dts
}

tmp <- bind_rows(
  generate_filtered_perf_curve_wn_stage_inter(method="PRR",score="PRR",type="score",x="fpr",y="tpr",stat="auc",thresh=1),
  generate_filtered_perf_curve_wn_stage_inter(method="PRR",score="PRR_90mse",type="90mse",x="fpr",y="tpr",stat="auc",thresh=1),
  generate_filtered_perf_curve_wn_stage_inter(method="GAM",score="gam_score",type="score",x="fpr",y="tpr",stat="auc",thresh=1),
  generate_filtered_perf_curve_wn_stage_inter(method="GAM",score="gam_score_90mse",type="90mse",x="fpr",y="tpr",stat="auc",thresh=1)
)

sub <- tmp[,.(spikein,method,type,nichd,
              tpr_lwr,tpr_mean,tpr_upr,
              fpr_lwr,fpr_mean,fpr_upr,
              tnr_lwr,tnr_mean,tnr_upr,
              fnr_lwr,fnr_mean,fnr_upr,
              ppv_lwr,ppv_mean,ppv_upr,
              npv_lwr,npv_mean,npv_upr,
              auc_lwr,auc_mean = auc,auc_upr)] %>% 
  unique() %>% 
  pivot_longer(
    cols=c("tpr_lwr","tpr_mean","tpr_upr",
           "fpr_lwr","fpr_mean","fpr_upr",
           "tnr_lwr","tnr_mean","tnr_upr",
           "fnr_lwr","fnr_mean","fnr_upr",
           "ppv_lwr","ppv_mean","ppv_upr",
           "npv_lwr","npv_mean","npv_upr",
           'auc_lwr',"auc_mean","auc_upr")
  ) %>% data.table()
sub$metric <- 
  sapply(str_split(sub$name,"_"),function(x){x[1]})
sub$statistic <- 
  sapply(str_split(sub$name,"_"),function(x){x[2]})
metric_names <- 
  list(
    "auc" = "AUROC",
    "tpr" = "Power",
    "ppv" = "Positive predictive value",
    "npv" = "Negative predictive value",
    "fpr" = "False positive rate",
    "tnr" = "True negative rate",
    "fnr" = "False negative rate"
  )
sub$metric_name <- 
  sapply(sub$metric,function(x){metric_names[[x]]}) %>% unname
g <- sub %>% 
  .[order(metric)] %>% 
  dcast(
    spikein + method + type + nichd + metric_name ~ statistic,
    value.var="value"
  ) %>% 
  .[type=="score" &
      metric_name %in% 
      sapply(c("auc","tpr","ppv","npv"),function(x){as.factor(metric_names[[x]])})
  ] %>% 
  ggplot(aes(factor(nichd,levels=category_levels$nichd),mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(
    factor(metric_name,levels= paste0(unname(metric_names)))~spikein,
    scales="free",labeller = label_wrap_gen(width = 15)) +
  xlab("") +
  ylab("Performance value") +
  theme(
    axis.text.x = element_text(angle=45,hjust=1,vjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_performance_metrics_wn_nichd_score.png"),g,width=10,height=7)

g <- sub %>% 
  .[order(metric)] %>% 
  dcast(
    spikein + method + type + nichd + metric_name ~ statistic,
    value.var="value"
  ) %>% 
  .[type=="90mse" &
      metric_name %in% 
      sapply(c("auc","tpr","ppv","npv"),function(x){as.factor(metric_names[[x]])})
  ] %>% 
  ggplot(aes(factor(nichd,levels=category_levels$nichd),mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(
    factor(metric_name,levels= paste0(unname(metric_names)))~spikein,
    scales="free",labeller = label_wrap_gen(width = 15)) +
  xlab("") +
  ylab("Performance value") +
  theme(
    axis.text.x = element_text(angle=45,hjust=1,vjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_performance_metrics_wn_nichd_90mse.png"),g,width=10,height=7)

g <- tmp[,
         .(method,type,spikein,nichd,lwr = tnr_lwr,m = tnr_mean,upr = tnr_upr)
] %>% 
  unique() %>% 
  ggplot(aes(factor(nichd,levels=category_levels$nichd),m,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein,scales="free") +
  xlab("") +
  ylab("True negative rate") +
  theme(
    axis.text.x = element_text(angle=45,hjust=1,vjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_performance_tnr_wn_nichd.png"),g,width=15,height=7)

g <- tmp[,
         .(method,type,spikein,nichd,lwr = tpr_lwr,m = tpr_mean,upr = tpr_upr)
] %>% 
  unique() %>% 
  ggplot(aes(factor(nichd,levels=category_levels$nichd),m,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein,scales="free") +
  xlab("") +
  ylab("True positive rate") +
  theme(
    axis.text.x = element_text(angle=45,hjust=1,vjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_performance_power_wn_nichd.png"),g,width=15,height=7)

g <- tmp[,
         .(method,type,spikein,nichd,lwr = npv_lwr,m = npv_mean,upr = npv_upr)
] %>% 
  unique() %>% 
  ggplot(aes(factor(nichd,levels=category_levels$nichd),m,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein,scales="free") +
  xlab("") +
  ylab("Negative predictive value") +
  theme(
    axis.text.x = element_text(angle=45,hjust=1,vjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_performance_npv_wn_nichd.png"),g,width=15,height=7)

g <- tmp[,
         .(method,type,spikein,nichd,lwr = ppv_lwr,m = ppv_mean,upr = ppv_upr)
] %>% 
  unique() %>% 
  ggplot(aes(factor(nichd,levels=category_levels$nichd),m,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein,scales="free") +
  xlab("") +
  ylab("Positive predictive value") +
  theme(
    axis.text.x = element_text(angle=45,hjust=1,vjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_performance_ppv_wn_nichd.png"),g,width=15,height=7)

g <- tmp[,
         .(method,type,spikein,nichd,auc_lwr,auc,auc_upr)
] %>% 
  unique() %>% 
  ggplot(aes(factor(nichd,levels=category_levels$nichd),auc,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=auc_lwr,ymax=auc_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~spikein,scales="free") +
  xlab("") +
  ylab("AUROC") +
  theme(
    axis.text.x = element_text(angle=45,hjust=1,vjust=1)
  )
ggsave(paste0(img_dir,basename,"detection_after_inter_performance_auc_wn_nichd.png"),g,width=15,height=7)

# Model time summary ----------------------------

tmp <- 
  bind_rows(
    positive_data[,.(ade,time_gam = time_gam,time_prr,control="positive ADEs")],
    negative_data[,.(ade,time_gam = time_gam,time_prr,control="negative ADEs")]
  ) %>% 
  unique()

setnames(tmp,c("time_gam","time_prr"),c("GAM","PRR"))

g <- tmp %>% 
  melt(id.vars=c("ade","control")) %>% 
  ggplot() +
  #geom_violin(aes(variable,value,fill=variable)) +
  #scale_fill_manual(values=score_colors) +
  ggbeeswarm::geom_quasirandom(aes(variable,value,color=variable)) +
  scale_color_manual(values=score_colors) +
  scale_y_continuous(trans="log10",labels = scales::number,breaks=c(0,0.01,1,10,100,300)) +
  facet_grid(~control) +
  xlab("") +
  ylab("Computation time\n(seconds)") +
  theme(
    legend.position = "none"
  )
ggsave(paste0(img_dir,"sim_gam_vs_prr_computation_time.png"),g,width=7,height=7)

# Low reporting performance ------------

powered_score_ade <- fread(paste0(data_dir,basename,"power_analysis_powered_ades.csv"))
setnames(powered_score_ade,"score","score_")

inc_reporting_perf <- function(method="gam_score",score="gam_score",type="score",thresh=0,nreports=c(1,2,3,4,5),stat="auc",high=T){
  
  if(high){rates <- c(1)}else{rates <- c(0,1)}
  
  dts <- NULL
  for(spikein_col in stage_spikein_class[,unique(spikein)]){
    if(type=="90mse"){
      powered <- powered_score_ade[grepl("90mse",score_) &
                                     spikein==spikein_col] %>% 
        .[,.(score_,ade)] %>% 
        unique() %>% 
        .[,.N,ade] %>% 
        .[N==2,ade]
    }else{
      powered <- powered_score_ade[!grepl("90mse",score_) &
                                     spikein==spikein_col] %>% 
        .[,.(score_,ade)] %>% 
        unique() %>% 
        .[,.N,ade] %>% 
        .[N==2,ade]
    }
    pos = 
      stage_positive_data[spikein==spikein_col & ade %in% powered] %>% 
      merge(stage_spikein_class[class %in% rates],by=c("nichd","spikein"))
    neg=stage_negative_data
    dt <- 
      foreach(nreport=nreports,.combine="rbind") %dopar% {
        y_true <- c(rep(1,nrow(pos[a<=nreport])),rep(0,nrow(neg[a<=nreport])))
        y_pred <- c(pos[a<=nreport,get(score)],neg[a<=nreport,get(score)])
        y_true <- y_true[is.finite(y_pred)]
        y_pred <- y_pred[is.finite(y_pred)]
        score_thresh_perf_dt = 
          lapply(1:100,
                 function(i){
                   set.seed(i)
                   sinds = sample(1:length(y_pred),length(y_pred),replace=T)
                   pos_class_mean = mean(y_pred[sinds][y_true[sinds]==1])
                   neg_class_mean = mean(y_pred[sinds][y_true[sinds]==0])
                   cond_true <- y_true[sinds]==1
                   cond_false <- !cond_true
                   pred_true <- y_pred[sinds]>=thresh
                   pred_false <- !pred_true
                   tp <- ((as.integer(pred_true)+as.integer(cond_true))==2) %>% which() %>% length()
                   fn <- ((as.integer(pred_false)+as.integer(cond_true))==2) %>% which() %>% length()
                   fp <- ((as.integer(pred_true)+as.integer(cond_false))==2) %>% which() %>% length()
                   tn <- ((as.integer(pred_false)+as.integer(cond_false))==2) %>% which() %>% length()
                   pauc = auc_probability(as.logical(y_true),y_pred,N=1e5,split=T)
                   data.table(score_threshold=thresh,tp=tp,tn=tn,fp=fp,fn=fn,
                              tpr=(tp/(tp+fn)),fpr=(fp/(tn+fp)),tnr=(tn/(tn+fp)),fnr=(fn/(fn+tp)),
                              ppv=(tp/(fp+tp)),npv=(tn/(tn+fn)),
                              pauc = sum(pauc),ptp = pauc[1],ptie=pauc[2],
                              pos_class_mean = pos_class_mean,neg_class_mean=neg_class_mean)
                 }) %>% 
          bind_rows()
        ci = auc_ci(y_pred,y_true,stat=stat)
        dt <- data.table()
        dt$nreport <- nreport
        dt$method <- method
        dt$type <- type
        dt$spikein <- spikein_col
        dt$N <- length(y_true)
        dt$auc_lwr <- ci[1]
        dt$auc <- ci[2]
        dt$auc_upr <- ci[3]
        dt$pauc_lwr <- score_thresh_perf_dt[,quantile(pauc,c(0.025))]
        dt$pauc <- score_thresh_perf_dt[,mean(pauc)]
        dt$pauc_upr <- score_thresh_perf_dt[,quantile(pauc,c(0.975))]
        dt$ptp_lwr <- score_thresh_perf_dt[,quantile(ptp,c(0.025))]
        dt$ptp <- score_thresh_perf_dt[,mean(ptp)]
        dt$ptp_upr <- score_thresh_perf_dt[,quantile(ptp,c(0.975))]
        dt$ptie_lwr <- score_thresh_perf_dt[,quantile(ptie,c(0.025))]
        dt$ptie <- score_thresh_perf_dt[,mean(ptie)]
        dt$ptie_upr <- score_thresh_perf_dt[,quantile(ptie,c(0.975))]
        dt$tpr_lwr <- score_thresh_perf_dt[,quantile(tpr,c(0.025),na.rm=T)] %>% unname
        dt$tpr_upr <- score_thresh_perf_dt[,quantile(tpr,c(0.975),na.rm=T)] %>% unname
        dt$tpr_mean <- score_thresh_perf_dt[,mean(tpr)]
        dt$tnr_lwr <- score_thresh_perf_dt[,quantile(tnr,c(0.025),na.rm=T)] %>% unname
        dt$tnr_upr <- score_thresh_perf_dt[,quantile(tnr,c(0.975),na.rm=T)] %>% unname
        dt$tnr_mean <- score_thresh_perf_dt[,mean(tnr)]
        dt$ppv_lwr <- score_thresh_perf_dt[,quantile(ppv,c(0.025),na.rm=T)] %>% unname
        dt$ppv_upr <- score_thresh_perf_dt[,quantile(ppv,c(0.975),na.rm=T)] %>% unname
        dt$ppv_mean <- score_thresh_perf_dt[,mean(ppv)]
        dt$npv_lwr <- score_thresh_perf_dt[,quantile(npv,c(0.025),na.rm=T)] %>% unname
        dt$npv_upr <- score_thresh_perf_dt[,quantile(npv,c(0.975),na.rm=T)] %>% unname
        dt$npv_mean <- score_thresh_perf_dt[,mean(npv)]
        dt$fpr_lwr <- score_thresh_perf_dt[,quantile(fpr,c(0.025),na.rm=T)] %>% unname
        dt$fpr_upr <- score_thresh_perf_dt[,quantile(fpr,c(0.975),na.rm=T)] %>% unname
        dt$fpr_mean <- score_thresh_perf_dt[,mean(fpr)]
        dt$fnr_lwr <- score_thresh_perf_dt[,quantile(fnr,c(0.025),na.rm=T)] %>% unname
        dt$fnr_upr <- score_thresh_perf_dt[,quantile(fnr,c(0.975),na.rm=T)] %>% unname
        dt$fnr_mean <- score_thresh_perf_dt[,mean(fnr)]
        dt$pos_class_score_lwr = score_thresh_perf_dt[,quantile(pos_class_mean,c(0.025),na.rm=T)]
        dt$pos_class_score_mean = score_thresh_perf_dt[,mean(pos_class_mean)]
        dt$pos_class_score_upr = score_thresh_perf_dt[,quantile(pos_class_mean,c(0.975),na.rm=T)]
        dt$neg_class_score_lwr = score_thresh_perf_dt[,quantile(neg_class_mean,c(0.025),na.rm=T)]
        dt$neg_class_score_mean = score_thresh_perf_dt[,mean(neg_class_mean)]
        dt$neg_class_score_upr = score_thresh_perf_dt[,quantile(neg_class_mean,c(0.975),na.rm=T)]
        dt
        
      }
    dts <- bind_rows(dts,dt)
    }
  dts
}

lst = list()

nreports=round(seq(0,100,10),1)
tmp <- 
  bind_rows(
    inc_reporting_perf(method="GAM",score="gam_score",type="score",thresh=0,nreports=nreports),
    inc_reporting_perf(method="GAM",score="gam_score_90mse",type="90mse",thresh=0,nreports=nreports),
    inc_reporting_perf(method="PRR",score="PRR",type="score",thresh=1,nreports=nreports),
    inc_reporting_perf(method="PRR",score="PRR_90mse",type="90mse",thresh=1,nreports=nreports)
  )

range="0to100by10"
lst[[range]] = tmp

range="0to100by10"
tmp <- lst[[range]]

sub <- tmp[,.(spikein,method,type,nreport,
              tpr_lwr,tpr_mean,tpr_upr,
              fpr_lwr,fpr_mean,fpr_upr,
              tnr_lwr,tnr_mean,tnr_upr,
              fnr_lwr,fnr_mean,fnr_upr,
              ppv_lwr,ppv_mean,ppv_upr,
              npv_lwr,npv_mean,npv_upr,
              auc_lwr,auc_mean = auc,auc_upr)] %>% 
  unique() %>% 
  pivot_longer(
    cols=c("tpr_lwr","tpr_mean","tpr_upr",
           "fpr_lwr","fpr_mean","fpr_upr",
           "tnr_lwr","tnr_mean","tnr_upr",
           "fnr_lwr","fnr_mean","fnr_upr",
           "ppv_lwr","ppv_mean","ppv_upr",
           "npv_lwr","npv_mean","npv_upr",
           'auc_lwr',"auc_mean","auc_upr")
  ) %>% data.table()
sub$metric <- 
  sapply(str_split(sub$name,"_"),function(x){x[1]})
sub$statistic <- 
  sapply(str_split(sub$name,"_"),function(x){x[2]})
metric_names <- 
  list(
    "auc" = "AUROC",
    "tpr" = "Power",
    "ppv" = "Positive predictive value",
    "npv" = "Negative predictive value",
    "fpr" = "False positive rate",
    "tnr" = "True negative rate",
    "fnr" = "False negative rate"
  )
sub$metric_name <- 
  sapply(sub$metric,function(x){metric_names[[x]]}) %>% unname

g <- sub %>% 
  .[order(metric)] %>% 
  dcast(
    spikein + method + type + metric_name + nreport ~ statistic,
    value.var="value"
  ) %>% 
  .[type=="score" &
      metric_name %in% 
      sapply(c("auc","tpr","ppv","npv"),function(x){as.factor(metric_names[[x]])})
  ] %>% 
  ggplot(aes(factor(nreport),mean,color=method)) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(position=position_dodge(0.5),aes(ymin=lwr,ymax=upr),width=0.1) +
  facet_grid(
    factor(metric_name,levels= paste0(unname(metric_names)))~spikein,
    scales="free",
    labeller = label_wrap_gen(width=15)) +
  scale_color_manual(values=score_colors) +
  xlab("Number of drug, event reports") +
  ylab("Performance value") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_by_nreports_",range,"_metrics_score.png"),g,width=10,height=7)
  
g <- sub %>% 
  .[order(metric)] %>% 
  dcast(
    spikein + method + type + metric_name + nreport ~ statistic,
    value.var="value"
  ) %>% 
  .[type=="90mse" &
      metric_name %in% 
      sapply(c("auc","tpr","ppv","npv"),function(x){as.factor(metric_names[[x]])})
  ] %>% 
  ggplot(aes(factor(nreport),mean,color=method)) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(position=position_dodge(0.5),aes(ymin=lwr,ymax=upr),width=0.1) +
  facet_grid(
    factor(metric_name,levels= paste0(unname(metric_names)))~spikein,
    scales="free",
    labeller = label_wrap_gen(width=15)) +
  scale_color_manual(values=score_colors) +
  xlab("Number of drug, event reports") +
  ylab("Performance value") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_by_nreports_",range,"_metrics_90mse.png"),g,width=10,height=7)

g <- 
  bind_rows(
    tmp[,
        .(nreport,method,type,spikein,
          lwr = tpr_lwr,median = tpr_mean,upr = tpr_upr)
    ]
  ) %>%
  ggplot(aes(factor(nreport),median,color=method)) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(position=position_dodge(0.5),aes(ymin=lwr,ymax=upr),width=0.1) +
  facet_grid(type~spikein) +
  scale_color_manual(values=score_colors) +
  scale_y_continuous(limits=c(0,1)) +
  xlab("Number of drug, event reports") +
  ylab("Power")
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_by_nreports_",range,"_power.png"),g,width=12,height=5)

g <- 
  bind_rows(
    tmp[,
        .(nreport,method,type,spikein,
          lwr = tnr_lwr,median = tnr_mean,upr = tnr_upr)
    ]
  ) %>%
  ggplot(aes(factor(nreport),median,color=method)) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(position=position_dodge(0.5),aes(ymin=lwr,ymax=upr),width=0.1) +
  facet_grid(type~spikein) +
  scale_color_manual(values=score_colors) +
  scale_y_continuous(limits=c(0,1)) +
  xlab("Number of drug, event reports") +
  ylab("True negative rate")
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_by_nreports_",range,"_tnr.png"),g,width=12,height=5)

g <- 
  bind_rows(
    tmp[,
        .(nreport,method,type,spikein,
          lwr = fpr_lwr,median = fpr_mean,upr = fpr_upr)
    ]
  ) %>%
  ggplot(aes(factor(nreport),median,color=method)) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(position=position_dodge(0.5),aes(ymin=lwr,ymax=upr),width=0.1) +
  facet_grid(type~spikein) +
  scale_color_manual(values=score_colors) +
  xlab("Number of drug, event reports") +
  ylab("False positive rate")
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_by_nreports_",range,"_fpr.png"),g,width=12,height=5)

g <- 
  bind_rows(
    tmp[,
        .(nreport,method,type,spikein,
          lwr = ppv_lwr,median = ppv_mean,upr = ppv_upr)
    ]
  ) %>%
  ggplot(aes(factor(nreport),median,color=method)) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(position=position_dodge(0.5),aes(ymin=lwr,ymax=upr),width=0.1) +
  facet_grid(type~spikein) +
  scale_color_manual(values=score_colors) +
  xlab("Number of drug, event reports") +
  ylab("Positive predictive value")
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_by_nreports_",range,"_ppv.png"),g,width=12,height=5)

g <- 
  bind_rows(
    tmp[,
        .(nreport,method,type,spikein,
          lwr = npv_lwr,median = npv_mean,upr = npv_upr)
    ]
  ) %>%
  ggplot(aes(factor(nreport),median,color=method)) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(position=position_dodge(0.5),aes(ymin=lwr,ymax=upr),width=0.1) +
  facet_grid(type~spikein) +
  scale_color_manual(values=score_colors) +
  xlab("Number of drug, event reports") +
  ylab("Negative predictive value")
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_by_nreports_",range,"_npv.png"),g,width=12,height=5)

g <- 
  bind_rows(
    tmp[,
        .(nreport,method,type,spikein,
          lwr = auc_lwr,median = auc,upr = auc_upr)
    ]
  ) %>%
  ggplot(aes(factor(nreport),median,color=method)) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(position=position_dodge(0.5),aes(ymin=lwr,ymax=upr),width=0.1) +
  facet_grid(type~spikein) +
  scale_color_manual(values=score_colors) +
  scale_y_continuous(limits=c(0,1)) +
  xlab("Number of drug, event reports") +
  ylab("AUROC")
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_by_nreports_",range,"_auc.png"),g,width=12,height=5)

# ADE dynamics drug report sensitivity analysis: single ADE -------------------------

gam_dt <-  
  fread(paste0(data_dir,"sim_generate_data_single_gam_drug_report_reduction_data.csv"))

prr_dt <-  
  fread(paste0(data_dir,"sim_generate_data_single_prr_drug_report_reduction_data.csv"))

sp="uniform"
g <- 
  gam_dt[
    spikein==sp & stage_reduced=="late_adolescence" &
      stage_reduced==nichd,
    .(nichd,stage_reduced,percent_drug_report_reduction,Drug = D, Event = E, ADE = DE)
  ] %>% 
  unique() %>% 
  pivot_longer(cols=c("Drug","Event","ADE")) %>% 
  ggplot(aes(percent_drug_report_reduction,value,color=name)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(trans="log10",labels=scales::comma) +
  guides(color=guide_legend(title="Reporting type")) +
  scale_color_brewer(palette="Set1") +
  xlab("Percent drug report reduction") +
  ylab("Number of reports")
ggsave(paste0(img_dir,"sensitivity_analysis_detection_D_E_DE_single_ade_one_stage_drug_report_reduction_summary.png"),g,width=5,height=3)


g <- 
  gam_dt[
    spikein==sp 
  ] %>% 
  unique() %>% 
  pivot_longer(cols=c("D","E","DE")) %>% 
  ggplot(aes(percent_drug_report_reduction,value,color=name)) +
  geom_point(size=0.5,position=position_jitterdodge(dodge.width = 0.5)) +
  geom_line() +
  scale_y_continuous(trans="log10",labels=scales::comma) +
  facet_grid(
    factor(stage_reduced,levels=category_levels[[stage_col]])~
      factor(nichd,levels=category_levels[[stage_col]])) +
  guides(color=guide_legend(title="Reporting type")) +
  scale_color_brewer(palette="Set1") +
  xlab("Percent drug report reduction") +
  ylab("Number of reports")
ggsave(paste0(img_dir,"sensitivity_analysis_detection_D_E_DE_single_ade_drug_report_reduction_summary.png"),g,width=25,height=25)

g <- 
  prr_dt[
    spikein==sp 
  ] %>% 
  unique() %>% 
  pivot_longer(cols=c("a","b","c","d")) %>% 
  ggplot(aes(percent_drug_report_reduction,value,color=name)) +
  geom_point(size=0.5,position=position_jitterdodge(dodge.width = 0.5)) +
  geom_line() +
  scale_y_continuous(trans="log10",labels=scales::comma) +
  facet_grid(
    factor(stage_reduced,levels=category_levels[[stage_col]])~
      factor(nichd,levels=category_levels[[stage_col]])) +
  guides(color=guide_legend(title="Reporting type")) +
  scale_color_brewer(palette="Set1")
ggsave(paste0(img_dir,"drug_report_sensitivity_analysis_detection_a_b_c_d_single_ade_drug_report_reduction_summary.png"),g,width=25,height=25)


g <- bind_rows(
  gam_dt[spikein==sp,.(nichd,score = gam_score,method="GAM",type="score",percent_drug_report_reduction,stage_reduced)],
  prr_dt[spikein==sp,.(nichd,score = log10(PRR),type="score",method="log10(PRR)",percent_drug_report_reduction,stage_reduced)]
) %>% 
  ggplot(aes(factor(nichd,levels=category_levels[[stage_col]]),score,color=percent_drug_report_reduction)) +
  geom_point() +
  geom_line(aes(group=percent_drug_report_reduction)) +
  facet_grid(method~factor(stage_reduced,levels=category_levels[[stage_col]]),scales="free") +
  scale_color_continuous(low="blue",high="red",
                         guide=guide_colorbar(title="Percent drug report reduction")) +
  xlab("") +
  ylab("Score") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,"drug_report_sensitivity_analysis_detection_score_single_ade_drug_report_reduction.png"),g,width=25,height=10)

g <- 
  bind_rows(
    prr_dt[spikein==sp,.(nichd,score = log10(PRR),method="log10(PRR)",type="score",percent_drug_report_reduction,stage_reduced)],
    gam_dt[spikein==sp,.(nichd,score = gam_score,method="GAM",type="score",percent_drug_report_reduction,stage_reduced)]
  ) %>% 
  ggplot(aes(percent_drug_report_reduction,score,color=method)) +
  geom_point() +
  geom_path() +
  scale_color_brewer(palette="Set1") +
  facet_grid(factor(stage_reduced,levels=category_levels[[stage_col]])~factor(nichd,levels=category_levels[[stage_col]]),scales="free") +
  xlab("") +
  ylab("Score") 
ggsave(paste0(img_dir,"drug_report_sensitivity_analysis_detection_score_single_ade_drug_report_reduction_comparison.png"),g,width=25,height=25)


# ADE dynamics drug report sensitivity analysis: all ADEs, all stages --------------------------------

dat <- fread(paste0(data_dir,"sim_generate_data_drug_report_reduction_data.csv"))
dat$percent_drug_report_reduction <- as.integer(100-dat$percent_drug_report_reduction)

colnames(dat)
dat[,.N,.(nichd,percent_drug_report_reduction,reduced_stage,spikein)][,unique(N)]

g <- dat[,
         .(PRR_N = sum(is.na(PRR))/500,
           GAM_N = sum(is.na(gam_score_90mse))/500
         ),
         .(spikein,reduced_stage,percent_drug_report_reduction,nichd)
][
  reduced_stage==nichd
] %>% 
  ggplot(aes(percent_drug_report_reduction,PRR_N,color=spikein)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(labels=scales::percent) +
  facet_grid(.~factor(nichd,levels=category_levels$nichd)) +
  scale_color_manual(values=dynamics_colors_no_uniform) +
  guides(color=guide_legend(title="Drug, event dynamics class",title.position="top")) +
  xlab("Percent drug report reduction") +
  ylab("Percent of NaN scores") +
  geom_abline(slope=0,intercept=500,color="red")
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_number_NaN_prr_drug_report_reduction_summary.png"),g,width=18,height=5)

g <- dat[,
         .(PRR_N = sum(is.finite(PRR))/500,
           GAM_N = sum(is.finite(gam_score_90mse))/500
         ),
         .(spikein,reduced_stage,percent_drug_report_reduction,nichd)
][
  reduced_stage==nichd
] %>% 
  ggplot(aes(percent_drug_report_reduction,PRR_N,color=spikein)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(labels=scales::percent) +
  facet_grid(.~factor(nichd,levels=category_levels$nichd)) +
  scale_color_manual(values=dynamics_colors_no_uniform) +
  xlab("Percent drug report reduction") +
  ylab("Percent of finite scores") +
  geom_abline(slope=0,intercept=500,color="red")
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_number_finite_prr_drug_report_reduction_summary.png"),g,width=18,height=5)

g <- dat[spikein=="increase",
         .(
           Event = mean(E),
           Drug = mean(D),
           ADE = mean(DE)
         )
         ,
         .(nichd,percent_drug_report_reduction,reduced_stage,spikein)
] %>% 
  pivot_longer(cols=c("Event","Drug","ADE")) %>% 
  ggplot(aes(percent_drug_report_reduction,value,color=name)) +
  geom_point(size=0.5,position=position_jitterdodge(dodge.width = 0.5)) +
  geom_line() +
  facet_grid(factor(reduced_stage,levels=category_levels[[stage_col]])~factor(nichd,levels=category_levels[[stage_col]]),
             scales="free") +
  scale_y_continuous(trans="log10",labels=scales::comma)+
  guides(color=guide_legend(title="Reporting type")) +
  scale_color_brewer(palette="Set1") +
  xlab("Percent drug report reduction") +
  ylab("Average number of reports")
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_D_E_DE_ade_drug_report_reduction_summary.png"),g,width=25,height=25)

g <- dat[spikein=='plateau' & is.finite(PRR),
         .(nichd,percent_drug_report_reduction,reduced_stage,PRR,GAM = gam_score)] %>% 
  pivot_longer(cols=c("PRR","GAM")) %>% 
  data.table() %>% 
  .[,
    .(
      score = mean(value,na.rm=T)
    ),
    .(nichd,method = name,percent_drug_report_reduction,reduced_stage)] %>% 
  ggplot(aes(factor(nichd,levels=category_levels[[stage_col]]),score,color=percent_drug_report_reduction)) +
  geom_point() +
  geom_line(aes(group=percent_drug_report_reduction)) +
  facet_grid(method~factor(reduced_stage,levels=category_levels[[stage_col]]),scales="free") +
  scale_color_continuous(low="blue",high="red",
                         guide=guide_colorbar(title="Percent drug report reduction")) +
  xlab("") +
  ylab("Average detection score") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_score_ade_drug_report_reduction.png"),g,width=25,height=10)

g <- dat[!is.na(PRR_90mse) & !is.na(PRR_90pse) & reduced_stage==nichd
  ,.(percent_drug_report_reduction,reduced_stage,nichd,ade,
     PRR = PRR_90pse - PRR_90mse,GAM = gam_score_90pse - gam_score_90mse)
] %>% 
  pivot_longer(
    cols=c("PRR","GAM"),
    names_to = "method",
    values_to = "error"
    ) %>% 
  data.table() %>% 
  .[,
    .(
      lwr = mean(error,na.rm=T) - sd(error,na.rm=T),
      err = mean(error,na.rm=T),
      upr = mean(error,na.rm=T) + sd(error,na.rm=T)
      ),
    .(percent_drug_report_reduction,reduced_stage,nichd,method)
    ] %>% 
  ggplot(aes(percent_drug_report_reduction,err,color=method)) +
  geom_point() +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1) +
  scale_color_manual(values=score_colors) +
  xlab("Percent drug report reduction") +
  ylab("Detection score error") +
  facet_grid(
    method~
      factor(nichd,levels=category_levels$nichd),
    scales="free")
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_score_error_ade_drug_report_reduction.png"),g,width=15,height=6)

powered_score_ades <- fread(paste0(data_dir,basename,"power_analysis_powered_ades.csv"))
setnames(powered_score_ades,"score","score_")

reporting_drug_report_reduction_perf <- function(dat,method="PRR",score="PRR",type="score",thresh=1,stat="auc",len=100,high=T){
  
  if(high){rates <- c(1)}else{rates <- c(0,1)}
  
  dts <- NULL
  for(st in dat[,unique(reduced_stage)]){
    dt <- foreach(q=dat[,unique(percent_drug_report_reduction)],.combine="rbind",.errorhandling = "remove") %dopar% {
      if(type=="90mse"){
        powered <- powered_score_ade[grepl("90mse",score_) &
                                       spikein==spikein_col] %>% 
          .[,.(score_,ade)] %>% 
          unique() %>% 
          .[,.N,ade] %>% 
          .[N==2,ade]
      }else{
        powered <- powered_score_ade[!grepl("90mse",score_) &
                                       spikein==spikein_col] %>% 
          .[,.(score_,ade)] %>% 
          unique() %>% 
          .[,.N,ade] %>% 
          .[N==2,ade]
      }
      pos = 
        dat[reduced_stage==st & percent_drug_report_reduction==q &
              nichd==st & ade %in% powered] %>% 
        merge(stage_spikein_class[class %in% rates],by=c("nichd","spikein"))
      neg=stage_negative_data[nichd==st]
      y_true <- c(rep(1,nrow(pos)),rep(0,nrow(neg)))
      y_pred <- c(pos[,get(score)],neg[,get(score)])
      y_true <- y_true[is.finite(y_pred)]
      y_pred <- y_pred[is.finite(y_pred)]
      a <- y_pred[y_true==1]
      b <- y_pred[y_true==0]
      diff <- mean(a)-mean(b)
      diffsamps = sapply(1:len,function(x){
        set.seed(x)
        ( mean(sample(a,size=min(c(length(a),length(b))),replace=T)) - mean(sample(b,size=min(c(length(a),length(b))),replace=T)) )
      })
      diffsamp_lwr = quantile(diffsamps,c(0.025))[1] %>% unname
      diffsamp_mean = mean(diffsamps)
      diffsamp_upr = quantile(diffsamps,c(0.975))[1] %>% unname
      
      ci = auc_ci(y_pred,y_true,stat=stat)
      perf_dt <- 
        lapply(1:len,function(i){
          set.seed(i)
          sinds = sample(1:length(y_pred),length(y_pred),replace=T)
          cond_true <- y_true[sinds]==1
          cond_false <- !cond_true
          pred_true <- y_pred[sinds]>thresh
          pred_false <- !pred_true
          tp <- sum((as.integer(cond_true) + as.integer(pred_true))==2)
          fn <- sum((as.integer(cond_true) + as.integer(pred_false))==2)
          fp <- sum((as.integer(cond_false) + as.integer(pred_true))==2)
          tn <- sum((as.integer(cond_false) + as.integer(pred_false))==2)
          data.table(power = tp/(tp+fn),fpr = (fp/(fp+tn)),
                     tnr=(tn/(tn+fp)),
                     ppv=(tp/(tp+fp)),npv=(tn/(tn+fn)))
        }) %>% bind_rows()
      power_lwr = quantile(perf_dt$power,c(0.025))[1] %>% unname
      power_mean = mean(perf_dt$power)
      power_upr = quantile(perf_dt$power,c(0.975))[1] %>% unname
      fpr_lwr = quantile(perf_dt$fpr,c(0.025))[1] %>% unname
      fpr_mean = mean(perf_dt$fpr)
      fpr_upr = quantile(perf_dt$fpr,c(0.975))[1] %>% unname
      ppv_lwr = quantile(perf_dt$ppv,c(0.025))[1] %>% unname
      ppv_mean = mean(perf_dt$ppv)
      ppv_upr = quantile(perf_dt$ppv,c(0.975))[1] %>% unname
      npv_lwr = quantile(perf_dt$npv,c(0.025))[1] %>% unname
      npv_mean = mean(perf_dt$npv)
      npv_upr = quantile(perf_dt$npv,c(0.975))[1] %>% unname
      dt <- data.table("percent_drug_report_reduction"=q,
                       "power_lwr" = power_lwr,"Power"=power_mean,
                       "power_upr" = power_upr,
                       "fpr_lwr" = fpr_lwr,"FPR"=fpr_mean,"fpr_upr"=fpr_upr,
                       "ppv_lwr" = ppv_lwr,"PPV"=ppv_mean,"ppv_upr"=ppv_upr,
                       "npv_lwr" = npv_lwr,"NPV"=npv_mean,"npv_upr"=npv_upr,
                       "auc_lwr" = ci[1],"AUC" = ci[2],"auc_upr" = ci[3],
                       "high_mean" = mean(a),"low_mean" = mean(b),"diff_mean" = diff,
                       "diffsamp" = diffsamp_mean,"diffsamp_lwr"=diffsamp_lwr,
                       "diffsamp_upr" = diffsamp_upr
      )
      dt$reduced_stage <- st
      dt$threshold <- thresh
      dt$method <- method
      dt$type <- type
      dt
    }
    dts <- rbind(dts,dt)
  }
  dts
}

tmp <- 
  bind_rows(
    reporting_drug_report_reduction_perf(dat,method="GAM",score="gam_score",type="score",thresh=0),
    reporting_drug_report_reduction_perf(dat,method="GAM",score="gam_score_90mse",type="90mse",thresh=0),
    reporting_drug_report_reduction_perf(dat,method="PRR",score="PRR",type="score",thresh=1),
    reporting_drug_report_reduction_perf(dat,method="PRR",score="PRR_90mse",type="90mse",thresh=1)
  )

tmp %>% 
  fwrite(paste0(data_dir,basename,"dynamics_sensitivity_analysis_drug_report_results.csv"))

tmp <- 
  fread(paste0(data_dir,basename,"dynamics_sensitivity_analysis_drug_report_results.csv"))

g <- 
  tmp %>% 
  .[type=="score"] %>% 
  ggplot(aes(percent_drug_report_reduction,diffsamp,color=method)) +
  geom_point(position=position_dodge(width=5)) +
  geom_errorbar(position=position_dodge(width=5),aes(ymin=diffsamp_lwr,ymax=diffsamp_upr),width=1) +
  facet_grid(method~factor(reduced_stage,levels=category_levels$nichd),scales="free") +
  scale_color_manual(values=score_colors) +
  xlab("Percent drug report reduction\nFrom all reported drugs to no drug reports") +
  ylab("High and low dynamic reporting score difference") +
  theme(legend.position = "none")
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_score_drug_report_reduction_score_difference.png"),g,width=20,height=7)

g <- 
  tmp %>% 
  .[type=="90mse"] %>% 
  ggplot(aes(percent_drug_report_reduction,diffsamp,color=method)) +
  geom_point(position=position_dodge(width=5)) +
  geom_errorbar(position=position_dodge(width=5),aes(ymin=diffsamp_lwr,ymax=diffsamp_upr),width=1) +
  facet_grid(method~factor(reduced_stage,levels=category_levels$nichd),scales="free") +
  scale_color_manual(values=score_colors) +
  xlab("Percent drug report reduction\nFrom all reported drugs to no drug reports") +
  ylab("High and low dynamic reporting score difference") +
  theme(legend.position = "none")
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_score_drug_report_reduction_90mse_difference.png"),g,width=20,height=7)

sub <- tmp[,.(percent_drug_report_reduction,
              method,type,reduced_stage,
              tpr_lwr = power_lwr,tpr_mean = Power,tpr_upr = power_upr,
              fpr_lwr,fpr_mean = FPR,fpr_upr,
              ppv_lwr,ppv_mean = PPV,ppv_upr,
              npv_lwr,npv_mean = NPV,npv_upr,
              auc_lwr,auc_mean = AUC,auc_upr)] %>% 
  unique() %>% 
  pivot_longer(
    cols=c("tpr_lwr","tpr_mean","tpr_upr",
           "fpr_lwr","fpr_mean","fpr_upr",
           "ppv_lwr","ppv_mean","ppv_upr",
           "npv_lwr","npv_mean","npv_upr",
           'auc_lwr',"auc_mean","auc_upr")
  ) %>% data.table()
sub$metric <- 
  sapply(str_split(sub$name,"_"),function(x){x[1]})
sub$statistic <- 
  sapply(str_split(sub$name,"_"),function(x){x[2]})
metric_names <- 
  list(
    "auc" = "AUROC",
    "tpr" = "Power",
    "ppv" = "Positive predictive value",
    "npv" = "Negative predictive value",
    "fpr" = "False positive rate"
  )
sub$metric_name <- 
  sapply(sub$metric,function(x){metric_names[[x]]}) %>% unname

g <- sub %>% 
  .[order(metric)] %>% 
  dcast(
    percent_drug_report_reduction + reduced_stage + method + type + metric_name ~ statistic,
    value.var="value"
  ) %>% 
  .[type=="score" &
      metric_name %in% 
      sapply(c("auc","tpr","ppv","npv"),function(x){as.factor(metric_names[[x]])})
  ] %>% 
  ggplot(aes(factor(percent_drug_report_reduction),mean,color=method)) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(position=position_dodge(0.5),aes(ymin=lwr,ymax=upr),width=0.1) +
  facet_grid(
    factor(metric_name,levels= paste0(unname(metric_names)))~
      factor(reduced_stage,levels=category_levels[[stage_col]]),
    scales="free",
    labeller = label_wrap_gen(width=15)) +
  scale_color_manual(values=score_colors) +
  xlab("Percent drug report reduction") +
  ylab("Performance value") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_score_drug_report_reduction_performance_metrics_score.png"),
       g,width=18,height=7)
g <- sub %>% 
  .[order(metric)] %>% 
  dcast(
    percent_drug_report_reduction + reduced_stage + method + type + metric_name ~ statistic,
    value.var="value"
  ) %>% 
  .[type=="90mse" &
      metric_name %in% 
      sapply(c("auc","tpr","ppv","npv"),function(x){as.factor(metric_names[[x]])})
  ] %>% 
  ggplot(aes(factor(percent_drug_report_reduction),mean,color=method)) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(position=position_dodge(0.5),aes(ymin=lwr,ymax=upr),width=0.1) +
  facet_grid(
    factor(metric_name,levels= paste0(unname(metric_names)))~
      factor(reduced_stage,levels=category_levels[[stage_col]]),
    scales="free",
    labeller = label_wrap_gen(width=15)) +
  scale_color_manual(values=score_colors) +
  xlab("Percent drug report reduction") +
  ylab("Performance value") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"sensitivity_analysis_detection_score_drug_report_reduction_performance_metrics_90mse.png"),
       g,width=18,height=7)
