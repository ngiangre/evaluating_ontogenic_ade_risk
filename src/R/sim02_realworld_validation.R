#' Title: "Evaluating risk detection methods to uncover ontogenic-mediated adverse drug effect mechanisms in children" study
#' 
#' Author details: Nicholas Giangreco
#' 
#' This script calculates the data description and performance 
#' statistics and figures for the real world validation results

# PURPOSE -----------------------------------------------------------------

#' To compare the performance of the GAM and PRR for detecting drug-event pairs in a real-world pediatric reference set 

# Load libraries and set variables -------------------------------------------------------------------


#install.packages("pacman")
pacman::p_load(tidyverse,data.table,mgcv,Rcpp,doParallel)

seed = 0
set.seed(seed)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "sim_realworld_validation_"

cores=5
registerDoParallel(cores=cores)

theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))


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

# Load reference sets ---------------------------------------------------


pos_neg_set <- 
  fread(paste0(data_dir,"GRiP_Pediatric_ADE_with_generated_negatives_meddrapt.csv"))

reference_set <- 
  fread(paste0(data_dir,"GRiP_Pediatric_ADE_Gold_Standard_List_minimal_joined.csv"))

ref_ades <- 
  reference_set[
    Population=="C" & 
      MedDRA_concept_class_id=="PT"
  ][,.(ATC_concept_name,MedDRA_concept_name,ATC_concept_id,MedDRA_concept_id)
  ] %>% 
  unique()

nrow(ref_ades)

setnames(ref_ades,colnames(ref_ades),c(drug_col_name,rx_col_name,drug_col,rx_col))

ref_ades$class <- 1

pos_neg_set_merge <- 
  merge(
    ref_ades,
    pos_neg_set[Control=="N"],
    by.x=c("atc_concept_id","meddra_concept_id"),
    by.y=c("ATC_concept_id","MedDRA_concept_id"),
    all=T
  ) 
pos_neg_set_merge[is.na(class),"class"] <- 0

ref_candidate_ades <- 
  inner_join(pos_neg_set_merge[,c(drug_col,rx_col,"class"),with=F],
             raw_data[,.N,c(drug_col,rx_col,drug_col_name,rx_col_name)],
             by=c(drug_col,rx_col)
  ) 

ref_candidate_ades$ade <- 
  paste0(ref_candidate_ades$atc_concept_id,"_",ref_candidate_ades$meddra_concept_id)

ref_candidate_ades[,.N,.(class)]


candidate_ades <- 
  inner_join(ref_ades,
             raw_data[,.N,c(drug_col,rx_col)]
  ) 

candidate_ades

# Define functions ---------------------------------------------------------------

#'
#'
make_ade_data <- function(){
  # Extract IDs
  ids <- raw_data[,unique(get(id_col))]
  dt <- data.table()
  dt[,id_col] <- ids
  dt[,e_col] <- 0
  dt[,d_col] <- 0
  dt[,de_col] <- 0
  
  # Map stages to reports
  group_cols <- c(id_col,x_cols)
  reports_case_control_dat <-
    data.table::merge.data.table(dt,
                                 raw_data[,group_cols,with=F] %>% unique(),
                                 by=id_col
    ) %>% unique
  #Order reports by age
  reports_case_control_dat <- reports_case_control_dat[match(raw_data[,unique(safetyreportid)],get(id_col))]
  gam_data <- reports_case_control_dat 
  # Convert stage categories to integers
  group_cols <- c(d_col,e_col,de_col,x_cols)
  for(x_col in group_cols){
    x <- gam_data[,x_col,with=F] %>%
      unlist %>% unname
    if(x_col %in% names(category_levels)){
      gam_data[,paste0(x_col,"_name")] <- gam_data[,c(x_col),with=F]
      gam_data[,x_col] <- match(x,category_levels[[x_col]])
    }
  }
  
  # Output
  return(gam_data)
  
}

#'
#'
set_ade_pair <- function(ade_data,drug_ids,rx_ids,drug_name,rx_name){
  #drug_ids=21604757;rx_ids=36919145;drug_name="methylphenidate";rx_name="Acute psychosis"
  # Extract event, drug, and not drug reports
  drug_reports <- raw_data[get(drug_col) %in% drug_ids,unique(get(id_col))]
  rx_reports <- raw_data[get(rx_col) %in% rx_ids,unique(get(id_col))]
  
  ade_data[,(e_col) := as.integer(get(id_col) %in% rx_reports)]
  ade_data[,(d_col) := as.integer(get(id_col) %in% drug_reports)]
  ade_data[,(de_col) := get(d_col) * get(e_col)]
  
  ade_data[,drug_name := drug_name]
  ade_data[,rx_name := rx_name]
  
  # Output
  return(ade_data)
  
}


#' Compute proportional reporting ratio of ADE across childhood
#' 
compute_prr <- function(drug_id,rx_id,drug_name,rx_name,dat_new){
  t0=Sys.time()
  # Set GAM data columns as factors
  dat_new[,d_col] <- factor(dat_new[,get(d_col)],levels=c(0,1))
  dat_new[,e_col] <- factor(dat_new[,get(e_col)],levels=c(0,1))
  # Make 2x2 ADE tables at stages
  tab <- dat_new[,
                 .(tab = list(table(D,E))
                 ),by=c(paste0(stage_col,"_name"))
  ]
  # Extract 2x2 parameters
  #https://www.ema.europa.eu/en/documents/regulatory-procedural-guideline/draft-guideline-use-statistical-signal-detection-methods-eudravigilance-data-analysis-system_en.pdf
  # stage order
  # a => drug, event
  # b => drug, no event
  # c => no drug, event
  # d => no drug, no event
  stage <- tab[,paste0(stage_col,"_name"),with=F] %>% unlist %>% unname
  a <- sapply(tab$tab,function(x){x["1","1"]})
  b <- sapply(tab$tab,function(x){x["1","0"]})
  # check:
  # dat_new[which(dat_new$D==1 & dat_new$E==0)][,table(get(stage_col))]
  # names(b) <- stage; b
  c <- sapply(tab$tab,function(x){x["0","1"]})
  d <- sapply(tab$tab,function(x){x["0","0"]})
  
  dt <- data.table()
  dt[,"a":=a]
  dt[,"b":=b]
  dt[,"c":=c]
  dt[,"d":=d]
  # (drug, events / all drug) / (no drug, events / all other drugs)
  dt[,"PRR":=(a/(a+b))/(c/(c+d))]
  dt[,"PRR_se":=sqrt( (1 / a) + (1 / c) - (1 / (a + b)) - (1/(c + d)) )]
  dt[,(stage_col):=stage]
  dt[,(drug_col_name):=drug_name]
  dt[,(rx_col_name):=rx_name]
  dt[,"ade_name"] = paste0(drug_name," and ",rx_name)
  dt[,(drug_col):=drug_id]
  dt[,(rx_col):=rx_id]
  dt[,"ade"] = paste0(drug_id,"_",rx_id)
  dt[,"PRR_90mse":=dt[,"PRR"]/exp(1.645*dt[,"PRR_se"])]
  dt[,"PRR_95mse":=dt[,"PRR"]/exp(1.96*dt[,"PRR_se"])]
  dt[,"PRR_99mse":=dt[,"PRR"]/exp(2.567*dt[,"PRR_se"])]
  dt[,"PRR_999mse":=dt[,"PRR"]/exp(3.291*dt[,"PRR_se"])]
  dt[,"PRR_90pse":=dt[,"PRR"]*exp(1.645*dt[,"PRR_se"])]
  dt[,"PRR_95pse":=dt[,"PRR"]*exp(1.96*dt[,"PRR_se"])]
  dt[,"PRR_99pse":=dt[,"PRR"]*exp(2.567*dt[,"PRR_se"])]
  dt[,"PRR_999pse":=dt[,"PRR"]*exp(3.291*dt[,"PRR_se"])]
  dt[,"time"] = difftime(Sys.time(),t0,"secs") %>% as.numeric()
  
  dt
  
}

#' Compute stage by drug reporting (spline) interaction 
#' 
compute_bam_base <- function(dat_new,bs_="cs"){
  
  #bs="cs" gives lower AIC and df - a better fit model
  bam_mod <-
    bam(E ~ s(get(stage_col),
              by=D,
              bs=bs_,
              k=length(stage_knots)),
        data=dat_new,
        knots = list(x=stage_knots),
        family=binomial(link="logit"),
        discrete=T
    )
  
  bam_mod
  
}

#' Process GAM model data
#' 
process_gam_model <- function(bam_mod,dat_new,drug_id,rx_id,drug_name,rx_name,knots,coef_str="s(get(stage_col)):D.",time_=7,pad=""){
  
  # Extract GAM coefficients
  # mapping spline knots to category names
  # add overall GAM score - weighted average across stages
  #https://www.andrewheiss.com/blog/2016/04/25/convert-logistic-regression-standard-errors-to-odds-ratios-with-r/
  #https://stats.stackexchange.com/questions/364568/extracting-significance-of-gam-smoothing-terms
  
  cnames=names(coef(bam_mod))
  inds = grepl(coef_str,cnames,fixed=T)
  cknots = str_replace(cnames[inds],fixed(coef_str),"") %>% as.integer()
  scores_gam_dt <- data.table(stage = knots)
  scores_gam_dt[match(cknots,scores_gam_dt$stage),"gam_score"] <- coef(bam_mod,complete=T)[inds]
  scores_gam_dt[match(cknots,scores_gam_dt$stage),"gam_score_se"] <- summary(bam_mod,freq=F)$se[inds]
  scores_gam_dt$N <- dat_new[,.N,by=stage_col]$N
  setnames(scores_gam_dt,c("stage"),c(stage_col))
  scores_gam_dt[[stage_col]] <- as.factor(category_levels[[stage_col]])
  
  # Mapping GAM data to categories
  # for summarizing column counts
  dat_new[,(stage_col)] <-
    dat_new[,lapply(.SD,function(x){category_levels[[stage_col]][x]}),.SDcols=stage_col]
  dat_agg <-
    dat_new[,
            lapply(.SD,function(x){sum(x==1)}),
            by=stage_col,
            .SDcols=c(e_col,d_col,de_col)
    ] %>%
    pivot_longer(-!!as.name(stage_col)) %>%
    data.table()
  dat_agg[[stage_col]] <- factor(dat_agg[[stage_col]] %>% unlist %>% unname,
                                 levels=category_levels[[stage_col]])
  
  # Computing confidence (credible) interval scores
  scores_gam_dt$gam_score_90mse <-
    scores_gam_dt[,.(gam_score_90mse = gam_score - (1.645*(gam_score_se)))]$gam_score_90mse
  scores_gam_dt$gam_score_95mse <-
    scores_gam_dt[,.(gam_score_95mse = gam_score - (1.96*(gam_score_se)))]$gam_score_95mse
  scores_gam_dt$gam_score_99mse <-
    scores_gam_dt[,.(gam_score_99mse = gam_score - (2.576*(gam_score_se)))]$gam_score_99mse
  scores_gam_dt$gam_score_999mse <-
    scores_gam_dt[,.(gam_score_999mse = gam_score - (3.291*(gam_score_se)))]$gam_score_999mse
  scores_gam_dt$gam_score_90pse <-
    scores_gam_dt[,.(gam_score_90pse = gam_score + (1.645*(gam_score_se)))]$gam_score_90pse
  scores_gam_dt$gam_score_95pse <-
    scores_gam_dt[,.(gam_score_95pse = gam_score + (1.96*(gam_score_se)))]$gam_score_95pse
  scores_gam_dt$gam_score_99pse <-
    scores_gam_dt[,.(gam_score_99pse = gam_score + (2.576*(gam_score_se)))]$gam_score_99pse
  scores_gam_dt$gam_score_999pse <-
    scores_gam_dt[,.(gam_score_999pse = gam_score + (3.291*(gam_score_se)))]$gam_score_999pse
  
  stage_dat_agg <- dat_agg[,.(value=sum(value,na.rm=T)),by=c(stage_col,"name")] %>% dcast(get(stage_col) ~ name)
  colnames(stage_dat_agg)[1] <- stage_col
  all_dat_agg <- dat_agg[,.(value=sum(value,na.rm=T)),by=c("name")]
  all_dat_agg[[stage_col]] <- "all"
  all_dat_agg <- all_dat_agg %>% dcast(get(stage_col) ~ name)
  colnames(all_dat_agg)[1] <- stage_col
  dat_agg <- bind_rows(stage_dat_agg,all_dat_agg)
  # Joining mapped GAM data
  # to GAM results
  joined <- inner_join(scores_gam_dt,dat_agg,by=stage_col) %>% data.table()
  summ <- summary(bam_mod)
  joined[,drug_col] = drug_id
  joined[,rx_col] = rx_id
  joined[,"ade"] = paste0(drug_id,"_",rx_id)
  joined[,drug_col_name] = drug_name
  joined[,rx_col_name] = rx_name
  joined[,"ade_name"] = paste0(drug_name," and ",rx_name)
  joined[,"formula"] = paste0(as.character(deparse(bam_mod$formula)),collapse="")
  joined[,"AIC"] = bam_mod$aic
  joined[,"gcv_ubre"] = bam_mod$gcv.ubre %>% unname
  
  y_true = bam_mod$y
  y_pred = bam_mod$fitted.values
  joined[,"mean_true_probability"] = mean(y_pred[y_true==1])
  joined[,"mean_false_probability"] = mean(y_pred[y_true==0])
  
  if(length(unique(y_true))==2){
    pred = ROCR::prediction(y_pred,y_true)
    pred = ROCR::performance(prediction.obj = pred,"auc")
    auc = pred@y.values[[1]]
    joined[,"auc"] = auc        
  }else{
    joined[,"auc"] = NA
  }
  
  joined[,"time"] = time_
  
  stringi::stri_sub(colnames(joined)[grepl("gam",colnames(joined))],4,3) <- pad
  joined
}

# Perform GAM modeling ------------------------------------------------------------


stage_knots=1:(category_levels[[stage_col]] %>% length())
e_col="E"
d_col="D"
de_col="DE"
x_cols=c(stage_col,"age")
ade_data <- make_ade_data()

res_prr <- 
  foreach(i=1:(nrow(ref_candidate_ades)),.combine = "rbind") %dopar% {
    p <- ref_candidate_ades[i]
    drug_id <- p[[drug_col]]
    rx_id <- p[[rx_col]]
    drug_name <- p[[drug_col_name]]
    rx_name <- p[[rx_col_name]]
    
    dat <- set_ade_pair(ade_data,drug_id,rx_id,drug_name,rx_name)
    
    compute_prr(drug_id,rx_id,drug_name,rx_name,dat)
    
  }

res_base <- 
foreach(i=1:(nrow(ref_candidate_ades)),.combine = "rbind") %dopar% {
    p <- ref_candidate_ades[i]
    drug_id <- p[[drug_col]]
    rx_id <- p[[rx_col]]
    drug_name <- p[[drug_col_name]]
    rx_name <- p[[rx_col_name]]
    
    dat <- set_ade_pair(ade_data,drug_id,rx_id,drug_name,rx_name)
    
    time_ <- system.time(mod <- compute_bam_base(dat))
    
    process_gam_model(mod,dat,drug_id,rx_id,drug_name,rx_name,knots=stage_knots,
                      coef_str="s(get(stage_col)):D.",time_=time_[1],pad="")

}

dfs <- 
  list("PRR" = res_prr,
       "base nichd spline interaction with drug" = res_base)

dfs %>% 
  write_rds(paste0(data_dir,basename,"method_data.rds"))



# Visualize ADE scores ---------------------------------------------------

dfs <- read_rds(paste0(data_dir,basename,"method_data.rds"))

score_colors <- 
  c("PRR" = RColorBrewer::brewer.pal(n=3,name = "Set1")[2],
    "GAM" = RColorBrewer::brewer.pal(n=3,name = "Set1")[1])

tmp <- 
  merge(
    merge(
      dfs[["PRR"]],
      candidate_ades[,c(drug_col,rx_col),with=F],
      by=c(drug_col,rx_col)
    ) %>% .[,.(ade_name,nichd,PRR,PRR_90mse,PRR_90pse)],
    merge(
      dfs[["base nichd spline interaction with drug"]],
      candidate_ades[,c(drug_col,rx_col),with=F],
      by=c(drug_col,rx_col)
    ) %>% .[,.(ade_name,nichd,gam_score,gam_score_90mse,gam_score_90pse)]
  )

g <- 
  tmp %>% 
  ggplot(aes(factor(nichd,levels=category_levels[[stage_col]]))) +
  geom_point(aes(y=gam_score),position = position_nudge(x = -0.2),color=score_colors[["GAM"]],size=2) +
  geom_errorbar(aes(ymin=gam_score_90mse,ymax=gam_score_90pse),
                position = position_nudge(x = -0.2),color=score_colors[["GAM"]],width=0.1,size=1) +
  geom_point(aes(y=log10(PRR)),position = position_nudge(x = 0.2),color=score_colors[["PRR"]],size=2) +
  geom_errorbar(aes(ymin=log10(PRR_90mse),ymax=log10(PRR_90pse)),
                position = position_nudge(x = 0.2),color=score_colors[["PRR"]],width=0.1,size=1) +
  scale_y_continuous(sec.axis=sec_axis(~.,name="")) +
  facet_wrap(~ade_name,ncol=5,scales="free") +
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
ggsave(paste0(img_dir,basename,"pediatric_reference_pt_prr_vs_gam.png"),g,
       width = 25,height=25)

pairs <- tmp[,unique(ade_name)]
for(pair in pairs){
  g <- 
    tmp[ade_name==pair] %>% 
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
    facet_wrap(~ade_name,ncol=5,scales="free") +
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
  ggsave(paste0(img_dir,basename,"pediatric_reference_",
                str_replace_all(pair," ","_"),"_prr_vs_gam.png"),g,
         width = 6,height=4)
  
}

# Visualize performance ---------------------------------------------------------------

dfs <- read_rds(paste0(data_dir,basename,"method_data.rds"))

score_colors <- 
  c("PRR" = RColorBrewer::brewer.pal(n=3,name = "Set1")[2],
    "GAM" = RColorBrewer::brewer.pal(n=3,name = "Set1")[1])


tmp1 <- dfs[["PRR"]]
tmp1$method <- "PRR"
tmp1 <- 
tmp1[,.(nichd,PRR,PRR_90mse,ade)]
tmp2 <- dfs[["base nichd spline interaction with drug"]]
tmp2$method <- "GAM"
tmp2 <- 
tmp2[,.(nichd,gam_score,gam_score_90mse,
        ade,ade_name)]

val_data <- 
  merge(
    tmp1,
    tmp2,
    by=c("ade","nichd")
    ) %>% 
    merge(
      ref_candidate_ades,
      by="ade"
    )

val_data[,.N,class]

val_data %>% 
  fwrite(paste0(data_dir,basename,"positive_and_negative_drug_event_pairs.csv"))

auc_ci <- function(y_pred,y_true,stat="auc",nboot=100){
  foreach(x=1:nboot,.combine = "c") %dopar% {
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
  } %>% quantile(probs=c(0.025,0.5,0.975),na.rm = T)
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

perf <- function(method="PRR",type="score",score="PRR",thresh=1,x="fpr",y="tpr",stat="auc",max=T){
  dts <- NULL
  sub = val_data[,.(max = max(get(score))),.(ade,class)]
  if(max){
    dat = sub
  }else{
    dat <- val_data
    dat$max <- dat[,get(score)]
  }
  y_pred <- dat[,max]
  y_true <- dat[,class]
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
    dt$method <- method
    dt$type <- type
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
    dts <- 
      bind_rows(
        dts,
        dt
      )
  }
  dts

}

tmp <- 
  bind_rows(
    perf(method="PRR",type="score",score="PRR",thresh=1,x="fpr",y="tpr",stat="auc",max=T),
    perf(method="GAM",type="score",score="gam_score",thresh=0,x="fpr",y="tpr",stat="auc",max=T),
    perf(method="PRR",type="90mse",score="PRR_90mse",thresh=1,x="fpr",y="tpr",stat="auc",max=T),
    perf(method="GAM",type="90mse",score="gam_score_90mse",thresh=0,x="fpr",y="tpr",stat="auc",max=T)
  )

g <- tmp %>% 
  ggplot(aes(fpr,tpr,color=method)) +
  geom_line() +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  geom_abline(intercept=0,slope=1,linetype=2,color="red") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  xlim(0,1) +
  ylim(0,1)
ggsave(paste0(img_dir,basename,"roc_curves.png"),g,height=6,width=6)

g <- tmp[,.(method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(method,tpr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=tpr_lwr,ymax=tpr_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  ylim(0,1) +
  xlab("") +
  ylab("TPR") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank()
  )
ggsave(paste0(img_dir,basename,"performance_tpr_at_score_threshold.png"),g,width=10,height=7)

tmp[,.(method,type,auc_lwr,auc,auc_upr)] %>% 
  unique()

g <- tmp[,.(method,type,auc_lwr,auc,auc_upr)] %>% 
  unique() %>% 
  ggplot(aes(method,auc,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=auc_lwr,ymax=auc_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("AUROC") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank()
  )
ggsave(paste0(img_dir,basename,"performance_auc.png"),g,width=10,height=7)



tmp <- 
  bind_rows(
    perf(method="PRR",type="score",score="PRR",thresh=1,x="fpr",y="tpr",stat="auc",max=F),
    perf(method="GAM",type="score",score="gam_score",thresh=0,x="fpr",y="tpr",stat="auc",max=F),
    perf(method="PRR",type="90mse",score="PRR_90mse",thresh=1,x="fpr",y="tpr",stat="auc",max=F),
    perf(method="GAM",type="90mse",score="gam_score_90mse",thresh=0,x="fpr",y="tpr",stat="auc",max=F)
  )

g <- tmp %>% 
  ggplot(aes(fpr,tpr,color=method)) +
  geom_line() +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  geom_abline(intercept=0,slope=1,linetype=2,color="red") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  xlim(0,1) +
  ylim(0,1)
ggsave(paste0(img_dir,basename,"all_stages_roc_curves.png"),g,height=6,width=6)

tmp[,.(method,type,
       tpr_lwr,tpr_mean,tpr_upr)] %>% 
  unique()

g <- tmp[,.(method,type,
            tpr_lwr,tpr_mean,tpr_upr,
            fpr_lwr,fpr_mean,fpr_upr,
            tnr_lwr,tnr_mean,tnr_upr,
            fnr_lwr,fnr_mean,fnr_upr,
            ppv_lwr,ppv_mean,ppv_upr,
            npv_lwr,npv_mean,npv_upr)] %>% 
  unique() %>% 
  ggplot(aes(method,tpr_mean,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=tpr_lwr,ymax=tpr_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  ylim(0,1) +
  xlab("") +
  ylab("TPR") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank()
  )
ggsave(paste0(img_dir,basename,"all_stages_performance_tpr_at_score_threshold.png"),g,width=10,height=7)

tmp[,.(method,type,auc_lwr,auc,auc_upr)] %>% 
  unique()

g <- tmp[,.(method,type,auc_lwr,auc,auc_upr)] %>% 
  unique() %>% 
  ggplot(aes(method,auc,color=method)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(aes(ymin=auc_lwr,ymax=auc_upr),width=0.1,position = position_dodge(width=.5)) +
  scale_color_manual(values=score_colors) +
  facet_grid(type~.) +
  xlab("") +
  ylab("AUROC") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank()
  )
ggsave(paste0(img_dir,basename,"all_stages_performance_auc.png"),g,width=10,height=7)

