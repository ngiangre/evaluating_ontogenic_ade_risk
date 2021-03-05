#' Title: "Evaluating risk detection methods to uncover ontogenic-mediated adverse drug effect mechanisms in children" study
#' 
#' Author details: Nicholas Giangreco
#' 
#' This script generates the positive and negative control data
#' and generates the sensitivity analysis data

# Purpose -----------------------------------------------------------------

#' Loads adverse drug event (ADE) data source,
#' Defines positive and negative drug-event pair datasets
#' Simulates drug event reporting dynamics,
#' Augments ADE data with the simulated reporting dynamics,
#' Applies the ADE detection methods,
#' And generates the sensitivity analysis data
#' 


# Load libraries and set variables ----------------------------------------------------------

#install.packages("pacman")
pacman::p_load(tidyverse,data.table,mgcv,Rcpp,doParallel)

seed = 0
set.seed(seed)

cores=50
registerDoParallel(cores=cores)

data_dir <- "../../data/"
out_dir <- paste0(data_dir,"ade_simulation_data/")
basename <- "sim_"

ifelse(
    !dir.exists(file.path(out_dir)),
    dir.create(file.path(out_dir)),
    FALSE)

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

# Define positive ADEs ----------------------------------------------

N = 500
all_ades <- unique(raw_data[,c(drug_col,rx_col),with=F])
rand <- sample(1:nrow(all_ades),size = N)
random_ade_pairs <- all_ades[rand]
positive_ade_pairs <- random_ade_pairs

positive_ade_pairs %>% 
  fwrite(paste0(data_dir,basename,"positive_ade_pairs.csv"))

# Define negative ADEs ----------------------------------------------------

drugs <- unique(random_ade_pairs[,get(drug_col)])
rxs <- unique(random_ade_pairs[,get(rx_col)])
negative_ade_pairs <-
  expand.grid(
    drugs,
    rxs
    ) %>% data.table
colnames(negative_ade_pairs) <- c(drug_col,rx_col)
negative_ade_pairs <-
    anti_join(
        negative_ade_pairs,
        positive_ade_pairs
    ) %>% data.table() %>%
    merge(all_ades) %>%
   sample_n(10000,replace=F)

negative_ade_pairs %>% 
  fwrite(paste0(data_dir,basename,"negative_ade_pairs.csv"))

# Define functions ---------------------------------------------------------------

# Generate fold change for producing ADE effect size
random_att_assign <- function(pairs,l=0.75,i=seed){
  set.seed(i)
  att=rexp(nrow(pairs),l) + 1
  pairs$fold_change <- att
  pairs
}

# Produce reporting rate according to dynamics shape (also called class or spikein)
define_report_probs <- function(type,N,tij){
    bprobs <- rep(tij,N)
    if(type=="uniform"){
        dy <- rep(0,N)
        rprobs <- bprobs+dy
    }
    if(type=="increase"){
        dy <- (tanh(seq(-pi,pi,length.out=N))*tij)
        rprobs <- bprobs+dy
    }
    if(type=="decrease"){
        dy <- -(tanh(seq(-pi,pi,length.out=N))*tij)
        rprobs <- bprobs+dy
    }
    if(type=="plateau"){
        dy <- c((tanh(seq(-pi,pi,length.out=floor(N/2)))*tij),(tanh(seq(pi,-pi,length.out=ceiling(N/2)))*tij))
        rprobs <- bprobs+dy
    }
    if(type=="inverse_plateau"){
        dy <- c((tanh(seq(pi,-pi,length.out=floor(N/2)))*tij),(tanh(seq(-pi,pi,length.out=ceiling(N/2)))*tij))
        rprobs <- bprobs+dy
    }
    rprobs
}

# Make and setup model data 
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

# Define drug, event indicators for model data
set_ade_pair <- function(ade_data,pair){
  # Extract ADE IDs
  drug_id <- pair[1]
  rx_id <- pair[2]
  drug_name <- all_data[atc_concept_id==drug_id,unique(atc_concept_name)]
  rx_name <- all_data[meddra_concept_id==rx_id,unique(meddra_concept_name)]
  # Extract event, drug, and not drug reports
  drug_reports <- raw_data[get(drug_col) %in% drug_id,unique(get(id_col))]
  rx_reports <- raw_data[get(rx_col) %in% rx_id,unique(get(id_col))]
  
  ade_data[,(e_col) := as.integer(get(id_col) %in% rx_reports)]
  ade_data[,(d_col) := as.integer(get(id_col) %in% drug_reports)]
  ade_data[,(de_col) := get(d_col) * get(e_col)]
  
  ade_data[,drug_name := drug_name]
  ade_data[,rx_name := rx_name]
  
  # Output
  return(ade_data)
  
}

# Format ADE data for drug, event pair
make_pair_data <- function(pair){
    # Extract ADE IDs
    drug_id <- pair[1]
    rx_id <- pair[2]
    drug_name <- all_data[atc_concept_id==drug_id,unique(atc_concept_name)]
    rx_name <- all_data[meddra_concept_id==rx_id,unique(meddra_concept_name)]
    # Extract event, drug, and not drug reports
    E_reports <-
        all_data[
            get(rx_col)==rx_id,
            unique(id_col),
            with=F] %>%
        unique %>%
        unlist %>%
        unname
    D_reports <-
        all_data[
            get(drug_col)==drug_id,
            unique(id_col),
            with=F] %>%
        unique %>%
        unlist %>%
        unname
    notD_reports <-
        all_data[
            get(drug_col)!=drug_id,
            unique(id_col),
            with=F] %>%
        unique %>%
        unlist %>%
        unname
    D_notD_reports <- union(D_reports,notD_reports)
    # Make GAM dataset
    # one hot encode events an drugs in reports
    dat <- data.table()
    dat[,(id_col) := D_notD_reports]
    dat[,(e_col) := as.integer(D_notD_reports %in% E_reports)]
    dat[,(d_col) := as.integer(D_notD_reports %in% D_reports)]
    dat[,(de_col) := as.integer(D_notD_reports %in% D_reports) * as.integer(D_notD_reports %in% E_reports)]
    # Map stages to reports
    group_cols <- c(id_col,x_cols)
    reports_case_control_dat <-
        data.table::merge.data.table(dat,
                                     all_data[,group_cols,with=F] %>% unique(),
                                     by=id_col
        ) %>% unique
    #Order reports by age
    reports_case_control_dat <- reports_case_control_dat[match(all_data[,unique(safetyreportid)],get(id_col))]
    reports_case_control_dat[,(id_col):=NULL]
    gam_data <- reports_case_control_dat %>% na.omit()
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
    gam_data

}

# Compute relative risk of ADE across childhood
compute_rrr <- function(pair,dat_new){
  t0 = Sys.time()
  # Extract ADE IDs
  drug_id <- pair[1]
  rx_id <- pair[2]
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
  # (drug, event / drug, no event) / (no drug, event / no drug, no event)
  dt[,"RRR":=(a/(b))/(c/(d))]
  dt[,"RRR_se":=sqrt( (1 / a) + (1 / c) - (1 / (c)) - (1/(d)) )]

  RRR_wm = dt[,sum(RRR * (a+b+c+d),na.rm=T)/sum( (a+b+c+d),na.rm=T)]
  si2 = dt[,(1/(a+b+c+d-1))*sum( (RRR - mean(RRR,na.rm=T) )^2,na.rm=T)]
  sp2 = dt[,sqrt( sum( (a+b+c+d-1) * si2, na.rm=T) )]

  dt <- rbind(dt,dt[,lapply(.SD,sum)])
  dt[nrow(dt),"RRR"] = RRR_wm
  dt[nrow(dt),"RRR_se"] = sqrt(sp2)

  dt[,(stage_col):=c(stage,"all")]
  dt[,(drug_col):=drug_id]
  dt[,(rx_col):=rx_id]
  dt[,"RRR_90mse":=dt[,"RRR"]/exp(1.645*dt[,"RRR_se"])]
  dt[,"RRR_95mse":=dt[,"RRR"]/exp(1.96*dt[,"RRR_se"])]
  dt[,"RRR_99mse":=dt[,"RRR"]/exp(2.567*dt[,"RRR_se"])]
  dt[,"RRR_999mse":=dt[,"RRR"]/exp(3.291*dt[,"RRR_se"])]
  dt[,"RRR_90pse":=dt[,"RRR"]*exp(1.645*dt[,"RRR_se"])]
  dt[,"RRR_95pse":=dt[,"RRR"]*exp(1.96*dt[,"RRR_se"])]
  dt[,"RRR_99pse":=dt[,"RRR"]*exp(2.567*dt[,"RRR_se"])]
  dt[,"RRR_999pse":=dt[,"RRR"]*exp(3.291*dt[,"RRR_se"])]
  dt[,"time"] = difftime(Sys.time(),t0,"secs") %>% as.numeric()

  dt

}

# Compute proportional reporting ratio of ADE across childhood
compute_prr <- function(pair,dat_new){
  t0=Sys.time()
    # Extract ADE IDs
    drug_id <- pair[1]
    rx_id <- pair[2]
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

    PRR_wm = dt[,sum(PRR * (a+b+c+d),na.rm=T)/sum( (a+b+c+d),na.rm=T)]
    si2 = dt[,(1/(a+b+c+d-1))*sum( (PRR - mean(PRR,na.rm=T) )^2,na.rm=T)]
    sp2 = dt[,sqrt( sum( (a+b+c+d-1) * si2, na.rm=T) )]

    dt <- rbind(dt,dt[,lapply(.SD,sum)])
    dt[nrow(dt),"PRR"] = PRR_wm
    dt[nrow(dt),"PRR_se"] = sqrt(sp2)

    dt[,(stage_col):=c(stage,"all")]
    dt[,(drug_col):=drug_id]
    dt[,(rx_col):=rx_id]
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

# Compute stage by drug reporting interaction within GLM 
compute_glm_f <- function(dat_new){

  dat_new[,stage_col] <- factor(dat_new[,stage_col,with=F] %>% unlist %>% unname,levels=1:length(stages))

  glm_mod <- glm(E ~ nichd:D,data=dat_new,family=binomial(link="logit"))

  glm_mod

}

# Compute stage by drug reporting (spline) interaction within GAM
compute_gam <- function(dat_new,bs_="cs"){

  nknots <- c(
    length(unique(dat_new[,get(stage_col)]))
  )

  #bs="cs" gives lower AIC and df - a better fit model
  bam_mod <-
    gam(E ~ s(get(stage_col),
              by=D,
              bs=bs_,
              k=nknots[1]),
        data=dat_new,
        knots = list(x=unique(dat_new[,stage_col,with=F]) %>% unlist %>% unname),
        family=binomial(link="logit"),
        nthreads=1
    )

  bam_mod

}

# Compute stage by drug reporting (spline) interaction within discretized GAM
compute_bam <- function(dat_new,bs_="cs"){
  
  nknots <- c(
    length(unique(dat_new[,get(stage_col)]))
  )
  
  #bs="cs" gives lower AIC and df - a better fit model
  bam_mod <-
    bam(E ~ s(get(stage_col),
              by=D,
              bs=bs_,
              k=nknots[1]),
        data=dat_new,
        knots = list(x=unique(dat_new[,stage_col,with=F]) %>% unlist %>% unname),
        family=binomial(link="logit"),
        discrete=T
    )
  
  bam_mod
  
}

# Compute stage by drug reporting interaction within GAM
compute_bam_te <- function(dat_new,bs_="cs"){
  
  nknots <- c(
    length(unique(dat_new[,get(stage_col)]))
  )
  
  bam_mod <-
    bam(E ~ te(get(stage_col),by=D,bs=bs_,k=nknots[1]),
        data=dat_new,
        knots = list(x=unique(dat_new[,stage_col,with=F]) %>% unlist %>% unname),
        family=binomial(link="logit"),
        discrete=T,nthreads=1
    )
  
  bam_mod
  
}

# Process GLM model data
process_glm_model <- function(glm_mod,dat_new,drug_id,rx_id,drug_name,rx_name,ej,fij,knots,coef_str="factor(nichd)",time_=3){

  scores_glm_dt <- data.table(stage = knots,
                              glm_score =
                                coef(glm_mod,complete=T)[
                                  grepl(coef_str,
                                        names(coef(glm_mod,complete=T)),fixed=T)
                                  ],
                              glm_score_se =
                                sqrt(diag(vcov(glm_mod)))[
                                  grepl(coef_str,
                                        names(sqrt(diag(vcov(glm_mod)))),fixed=T)
                                  ],
                              N = dat_new[,.N,by=stage_col]$N
  )
  setnames(scores_glm_dt,c("stage"),c(stage_col))
  scores_glm_dt[[stage_col]] <- as.factor(category_levels[[stage_col]])

  # Mapping GLM data to categories
  # for summarizing column counts
  dat_new[,(x_cols[1])] <-
    dat_new[,lapply(.SD,function(x){category_levels[[x_cols[1]]][x]}),.SDcols=x_cols[1]]
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

  tmp <- scores_glm_dt
  tmp$si2 <- scores_glm_dt[,.(si2 =  sum( (glm_score - mean(glm_score) )^2,na.rm=T ) * (1/(N-1)))]$si2
  sp2 <- tmp[,sum( (N-1)*si2,na.rm=T ) / sum( (N-1),na.rm=T )]

  all_score_glm_dt <- data.table("all", scores_glm_dt[,sum(glm_score*N,na.rm=T)/sum(N,na.rm=T)],sqrt(sp2),nrow(dat_new))
  colnames(all_score_glm_dt) <- c(stage_col,"glm_score","glm_score_se","N")
  scores_glm_dt <- bind_rows(scores_glm_dt,all_score_glm_dt)
  #scores_glm_dt$glm_score <- scores_glm_dt[,exp(glm_score)]
  scores_glm_dt$glm_score_90mse <-
    scores_glm_dt[,.(glm_score_90mse = glm_score - (1.645*(glm_score_se)))]$glm_score_90mse
  scores_glm_dt$glm_score_95mse <-
    scores_glm_dt[,.(glm_score_95mse = glm_score - (1.96*(glm_score_se)))]$glm_score_95mse
  scores_glm_dt$glm_score_99mse <-
    scores_glm_dt[,.(glm_score_99mse = glm_score - (2.576*(glm_score_se)))]$glm_score_99mse
  scores_glm_dt$glm_score_999mse <-
    scores_glm_dt[,.(glm_score_999mse = glm_score - (3.291*(glm_score_se)))]$glm_score_999mse
  scores_glm_dt$glm_score_90pse <-
    scores_glm_dt[,.(glm_score_90pse = glm_score + (1.645*(glm_score_se)))]$glm_score_90pse
  scores_glm_dt$glm_score_95pse <-
    scores_glm_dt[,.(glm_score_95pse = glm_score + (1.96*(glm_score_se)))]$glm_score_95pse
  scores_glm_dt$glm_score_99pse <-
    scores_glm_dt[,.(glm_score_99pse = glm_score + (2.576*(glm_score_se)))]$glm_score_99pse
  scores_glm_dt$glm_score_999pse <-
    scores_glm_dt[,.(glm_score_999pse = glm_score + (3.291*(glm_score_se)))]$glm_score_999pse

  stage_dat_agg <- dat_agg[,.(value=sum(value,na.rm=T)),by=c(stage_col,"name")] %>% dcast(get(stage_col) ~ name)
  colnames(stage_dat_agg)[1] <- stage_col
  all_dat_agg <- dat_agg[,.(value=sum(value,na.rm=T)),by=c("name")]
  all_dat_agg[[stage_col]] <- "all"
  all_dat_agg <- all_dat_agg %>% dcast(get(stage_col) ~ name)
  colnames(all_dat_agg)[1] <- stage_col
  dat_agg <- bind_rows(stage_dat_agg,all_dat_agg)
  # Joining mapped glm data
  # to glm results
  joined <- inner_join(scores_glm_dt,dat_agg,by=stage_col) %>% data.table()
  joined[,drug_col] = drug_id
  joined[,rx_col] = rx_id
  joined[,drug_col_name] = drug_name
  joined[,rx_col_name] = rx_name
  joined[,"Ej"] = ej
  joined[,"Fij"] = fij
  joined[,"Tij"] = fij * ej
  joined[,"AIC"] = glm_mod$aic
  joined[,"time"] = time_

  joined
}

# Process GAM model data
process_gam_model <- function(bam_mod,dat_new,drug_id,rx_id,drug_name,rx_name,ej,fij,knots,coef_str="s(get(stage_col)):D",time_=7,pad=""){

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
  dat_new[,(x_cols[1])] <-
    dat_new[,lapply(.SD,function(x){category_levels[[x_cols[1]]][x]}),.SDcols=x_cols[1]]
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

  # Computing weighted average score with pooled variance
  tmp <- scores_gam_dt
  tmp$si2 <- scores_gam_dt[,.(si2 =  sum( (gam_score - mean(gam_score,na.rm=T) )^2,na.rm=T ) * (1/(N-1)))]$si2
  sp2 <- tmp[,sum( (N-1)*si2,na.rm=T ) / sum( (N-1),na.rm=T )]

  all_score_gam_dt <- data.table("all", scores_gam_dt[,sum(gam_score*N,na.rm=T)/sum(N,na.rm=T)],sqrt(sp2),nrow(dat_new))
  colnames(all_score_gam_dt) <- c(stage_col,"gam_score","gam_score_se","N")
  scores_gam_dt <- bind_rows(scores_gam_dt,all_score_gam_dt)
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
  joined[,drug_col] = drug_id
  joined[,rx_col] = rx_id
  joined[,drug_col_name] = drug_name
  joined[,rx_col_name] = rx_name
  joined[,"Ej"] = ej
  joined[,"Fij"] = fij
  joined[,"Tij"] = fij * ej
  joined[,"AIC"] = bam_mod$aic
  joined[,"time"] = time_

  stringi::stri_sub(colnames(joined)[grepl("gam",colnames(joined))],4,3) <- pad
  colnames(joined)[grepl("gam",colnames(joined))]
  joined
}

# Extract GAM model specs
extract_gam_glance <- function(bam_mod){

  glance_gam <- broom::glance(bam_mod) %>% data.table()
  summary_bam <- summary.gam(bam_mod)
  glance_gam$rsq = summary_bam$r.sq
  glance_gam$fREML <- as.numeric(summary_bam$sp.criterion)
  glance_gam$deviance_explained <- summary_bam$dev.expl
  glance_gam[,drug_col] = drug_id
  glance_gam[,rx_col] = rx_id
  glance_gam[,drug_col_name] = drug_name
  glance_gam[,rx_col_name] = rx_name

  glance_gam

}

# Extract GAM parameter specs 
extract_gam_tidy <- function(bam_mod){

  tidy_gam <- broom::tidy(bam_mod,parametric=F) %>%
    bind_rows(broom::tidy(bam_mod,parametric=T)) %>%
    data.table()
  tidy_gam[,drug_col] = drug_id
  tidy_gam[,rx_col] = rx_id
  tidy_gam[,drug_col_name] = drug_name
  tidy_gam[,rx_col_name] = rx_name

  tidy_gam

}

# Main function that collects applied method results in list
compute_risk_and_process_data <- function(pair,dat_new,ej=NA,fij=NA){
    # Extract ADE IDs and names
    drug_id <- pair[1]
    rx_id <- pair[2]
    drug_name <- all_data[atc_concept_id==drug_id,unique(atc_concept_name)]
    rx_name <- all_data[meddra_concept_id==rx_id,unique(meddra_concept_name)]

    time_ <- system.time(gam_mod <- compute_gam(dat_new))
    gam_cs <- process_gam_model(
      gam_mod,
      dat_new,
      drug_id,rx_id,drug_name,rx_name,ej,fij,
      unique(dat_new[,stage_col,with=F]) %>% unlist %>% unname,
      coef_str="s(get(stage_col)):D.",time_[3],pad="_g"
    )
    time_ <- system.time(bam_mod <- compute_bam(dat_new))
    bam_cs <- process_gam_model(
      bam_mod,
      dat_new,
      drug_id,rx_id,drug_name,rx_name,ej,fij,
      unique(dat_new[,stage_col,with=F]) %>% unlist %>% unname,
      coef_str="s(get(stage_col)):D.",time_[3]
    )
    time_ <- system.time(gam_mod <- compute_bam_te(dat_new))
    bam_te_cs <- process_gam_model(
      gam_mod,
      dat_new,
      drug_id,rx_id,drug_name,rx_name,ej,fij,
      unique(dat_new[,stage_col,with=F]) %>% unlist %>% unname,
      coef_str="te(get(stage_col)):D.",time_[3],pad="_te"
    )
    # time_ <- system.time(glm_mod <- compute_glm_f(dat_new))
    # glm_f <- process_glm_model(
    #   glm_mod,
    #   dat_new,
    #   drug_id,rx_id,drug_name,rx_name,ej,fij,
    #   unique(dat_new[,stage_col,with=F]) %>% unlist %>% unname,
    #   coef_str="nichd",time_[3]
    # )

    # Output list of summarized datta objects
    lst = list(
      #"gam_glance" = extract_gam_glance(compute_gam(pair,dat_new)),
      #"gam_tidy" = extract_gam_tidy(compute_gam(pair,dat_new)),
      #"glm" = glm_f,
      "gam" = bam_cs,
      "gam_g" = gam_cs,
      "gam_te" = bam_te_cs, 
      "prr" = compute_prr(pair,dat_new)
      #"rrr" = compute_rrr(pair,dat_new)
    )

}

# Fast way to get max of drug, event reporting - for not removing drug event 
# reports if they were there before simulation
cppFunction("
  NumericVector max_vec(NumericVector vec1, NumericVector vec2) {
    int n = vec1.size();
    if(n != vec2.size()) return 0;
    else {
      NumericVector out(n);
      for(int i = 0; i < n; i++) {
        out[i] = std::max(vec1[i], vec2[i]);
      }
      return out;
    }
  }
")
# Generate positive control data -----------------------------------------------------------
#
# Applying detection methods to augmented ADE data
# Positive ADE pairs
#

# Define global variables
d_col="D"
e_col="E"
de_col="DE"
stage_col="nichd"
x_cols=c(stage_col,"age")
all_data = raw_data
stages <- all_data[,unique(get(stage_col))]
stage_age <- raw_data[,.(nichd,age)] %>% unique()
ade_data <- make_ade_data()

# For each of the positive ADE pairs,
#### Format ADE data
#### Apply detection methods
#### Output

# Get algorithm start time
t0 = Sys.time()
# Assign fold changes, "effect size", so ADE pairs
l=0.75
positive_ade_pairs_fc <- random_att_assign(positive_ade_pairs,l=l)
# Choose dynamic shapes or spikeins for drug event simulated reporting rate generation
spikein_cols = c("uniform","increase","decrease","plateau","inverse_plateau")
# Go through each dynamics class
for(spikein_col in spikein_cols){
  cat("\nSpike in: ",spikein_col,"\n")
  # Parallize through each ADE pair
  lst = foreach(i=1:nrow(positive_ade_pairs)) %dopar% {
    # Set random number generator seed
    set.seed(i)
    # Get pair ids
    pair <- positive_ade_pairs[i] %>% unlist %>% unname
    # Format ADE data
    dat = set_ade_pair(ade_data,pair)
    # Define ADE simulated rate
    e_old <- dat[,get(e_col)]
    ej = mean(e_old)
    fij = positive_ade_pairs_fc[i]$fold_change
    ej_fij = fij * ej
    #define data container for augmentation
    dat_new <- dat
    # Generate event dynamic data
    stage_age$probs = define_report_probs(spikein_col,nrow(stage_age),ej_fij)
    rs <- merge(dat_new,stage_age,by="age")$probs
    # Insert event dynamic data at drug reports for simulated drug event rate
    e_new <- rbinom(nrow(dat_new),1,rs)
    dat_new[D==1,"E"] <- max_vec(e_old,e_new)[dat_new$D==1]
    dat_new[E==1 & D==1,"DE"] <- 1
    # Format simulated ADE data
    dat_new$Di = dat_new[,mean(D)]
    dat_new$Ej = dat_new[,mean(E)]
    dat_new$Fij = fij
    dat_new$Tij = ej_fij
    dat_new$random_state = rs
    dat_new$random_seed = i
    dat_new[,c(drug_col) := list(pair[1])]
    dat_new[,c(rx_col) := list(pair[2])]
    # Apply ADE detection methods
    lst = compute_risk_and_process_data(pair,dat_new,ej,fij)
    # Make output names
    base =
      paste0(out_dir,pair[1],"_",pair[2],"_",stage_col,
             "_simulation_positive_l",l,"_",spikein_col,"_rate")
    out_file <- paste0(base,".rds")
    # Output
    saveRDS(lst,out_file)
  }
  # Get algorithm end time
  t1 = Sys.time()
  cat("\n",round(as.numeric(difftime(t1,t0,units="mins")),2)," minutes\n")
}


# Generate negative control data -----------------------------------------------------------
#
# Applying detection methods to negative ADE data
# Negative ADE pairs
#

# Define global variables
d_col="D"
e_col="E"
de_col="DE"
stage_col="nichd"
x_cols=c(stage_col)
all_data = raw_data
stages <- all_data[,unique(get(stage_col))]
ade_data <- make_ade_data()

# For each of the negative ADE pairs,
#### Format ADE data
#### Apply detection methods
#### Output
#
# Get algorithm start time
t0 = Sys.time()
# Parallize through each ADE pair
lst = foreach(i=1:nrow(negative_ade_pairs)) %dopar% {
  # Set random number generator seed
  set.seed(i)
  # Extract ADE pair IDs
  pair <- negative_ade_pairs[i] %>% unique() %>% unlist %>% unname
  # Format ADE data
  dat = set_ade_pair(ade_data,pair)
  # Apply ADE detection methods
  lst = compute_risk_and_process_data(pair,dat)
  # Make output names
  base <- paste0(out_dir,pair[1],"_",pair[2],"_",stage_col)
  out_file <- paste0(base,"_simulation_negative.rds")
  # Output
  saveRDS(lst,out_file)
}
# Get algorithm end time
t1 = Sys.time()
cat("\n",round(as.numeric(difftime(t1,t0,units="mins")),2)," minutes\n")

# Generate original positive data -----------------------------------------------------------
#
# Applying detection methods to positive (no augmentation) ADE data
# Positive ADE pairs
# Not used in paper but generated for QC
#

# Define global variables
d_col="D"
e_col="E"
de_col="DE"
stage_col="nichd"
x_cols=c(stage_col)
all_data = raw_data
stages <- all_data[,unique(get(stage_col))]
ade_data <- make_ade_data()

# For each of the positive ADE pairs,
#### Format ADE data
#### Apply detection methods
#### Output
#
# Get algorithm start time
t0 = Sys.time()
# Parallize through each ADE pair
lst = foreach(i=1:nrow(positive_ade_pairs)) %dopar% {
  # Set random number generator seed
  set.seed(i)
  # Extract ADE pair IDs
  pair <- positive_ade_pairs[i] %>% unique() %>% unlist %>% unname
  # Format ADE data
  dat = set_ade_pair(ade_data,pair)
  # Apply ADE detection methods
  lst = compute_risk_and_process_data(pair,dat)
  # Make output names
  base <- paste0(out_dir,pair[1],"_",pair[2],"_",stage_col)
  out_file <- paste0(base,"_simulation_original_positive.rds")
  # Output
  saveRDS(lst,out_file)
}
# Get algorithm end time
t1 = Sys.time()
cat("\n",round(as.numeric(difftime(t1,t0,units="mins")),2)," minutes\n")



# Calculate sensitivity analysis for single ADE -----------------------------------------------------
cat("\nEvent reporting reduction, single ADE\n")
#
# Applying detection methods to augmented ADE data at reduced event rates
# A single positive ADE pairs
#

# Define global variables
d_col="D"
e_col="E"
de_col="DE"
stage_col="nichd"
x_cols=c(stage_col,"age")
all_data = raw_data
stages <- all_data[,unique(get(stage_col))]
stage_age <- raw_data[,.(nichd,age)] %>% unique()
ade_data <- make_ade_data()

# Get algorithm start time
t0 = Sys.time()
# Assign fold changes, "effect size", to ADE pairs
l=0.75
positive_ade_pairs_fc <- random_att_assign(positive_ade_pairs,l=l)

# Choose ADE pair
i=403
# Set random number generator seed
set.seed(i)
# Define stage to reduce and pick pair with many reports in that stage
stages_to_reduce <- category_levels[[stage_col]]
prr_dts <- NULL
gam_dts <- NULL
for(stage_to_reduce in stages_to_reduce){

  # Get pair ids
  pair <- positive_ade_pairs[i] %>% unlist %>% unname
  # Format ADE data - get "old" event reporting
  dat = set_ade_pair(ade_data,pair)
  e_old <- dat[,get(e_col)]
  # Define reporting rate
  ej = mean(e_old)
  fij = positive_ade_pairs_fc[i]$fold_change
  ej_fij = fij * ej
  # Go through dynamics class
  spikein_cols = c("uniform")
  for(spikein_col in spikein_cols){
    cat("\nSpike in: ",spikein_col,"\n")
    #define data container for augmentation
    dat_new <- dat
    # Generate event dynamic data
    stage_age$probs = define_report_probs(spikein_col,nrow(stage_age),ej_fij)
    rs <- merge(dat_new,stage_age,by="age")$probs
    # Insert event dynamic data at drug reports for simulated drug event rate
    e_new <- rbinom(nrow(dat_new),1,rs)
    dat_new[D==1,e_col] <- max_vec(e_old,e_new)[dat_new$D==1]
    dat_new[E==1 & D==1,de_col] <- 1
    # Define event reporting after data simulation/augmentation
    e_new <- dat_new[,get(e_col)]

    lsts <- list()
    # Percentiles for data reduction
    step=0.1
    qs <- seq(0,1,step)
    # Reduce event reporting only at stage and apply detecction methods
    for(q in qs){
      cat("\nPercent reduction at ",stage_to_reduce,": ",q,"\n")
      # define reduced data container
      dat_q <- dat_new

      # Define event indices to sample
      stage_to_reduce_indices <- which(dat_q[,get(stage_col)==which(category_levels[[stage_col]]==stage_to_reduce)])
      stages_not_to_reduce_indices <- which(dat_q[,get(stage_col)!=which(category_levels[[stage_col]]==stage_to_reduce)])
      event_indices <- intersect(which(e_new==1),stage_to_reduce_indices)
      event_at_other_stages_indices <- intersect(which(e_new==1),stages_not_to_reduce_indices)
      # Create new event reporting vector
      e_q <- rep(0,length(e_new))
      e_q[event_at_other_stages_indices] <- 1
      # Sample reports with events per percentile rediction
      sampled_event_indices <-
        sample(
          event_indices,
          floor(q*length(event_indices)),
          replace = F
        )
      # Insert sampled event reporting into event vector
      if(length(sampled_event_indices)>0)e_q[sampled_event_indices] <- 1
      dat_q[,e_col] <- e_q
      dat_q[E==1 & D==1,de_col] <- 1
      dat_q[E==0,de_col] <- 0

      # Format simulated ADE data
      dat_q$Di = dat_q[,mean(D)]
      dat_q$Ej = dat_q[,mean(E)]
      dat_q$Fij = fij
      dat_q$Tij = ej_fij
      dat_q$random_state = rs
      dat_q$random_seed = i
      dat_q[,c(drug_col) := list(pair[1])]
      dat_q[,c(rx_col) := list(pair[2])]
      # Apply ADE detection methods and attribute data reduction
      lst = compute_risk_and_process_data(pair,dat_q,ej,fij)
      for(x in names(lst)){
        lst[[x]]$percent_event_reduction = (1-q)*100
      }

      lsts[[as.character(q)]] <- lst
    }

    gam_dt <- lapply(lsts,
                     function(x){x$gam}
    ) %>%
      bind_rows() %>%
      .[nichd!="all"]
    gam_dt$spikein <- spikein_col
    gam_dt$stage_reduced <- stage_to_reduce

    prr_dt <- lapply(lsts,
                     function(x){x$prr}
    ) %>%
      bind_rows() %>%
      .[nichd!="all"]
    prr_dt$spikein <- spikein_col
    prr_dt$stage_reduced <- stage_to_reduce

    gam_dts <- bind_rows(gam_dts,gam_dt)
    prr_dts <- bind_rows(prr_dts,prr_dt)
  }
  # Get algorithm end time
  t1 = Sys.time()
  cat("\n",round(as.numeric(difftime(t1,t0,units="mins")),2)," minutes\n")

}
# Get algorithm end time
t1 = Sys.time()
cat("\n",round(as.numeric(difftime(t1,t0,units="mins")),2)," minutes\n")


gam_dts %>%
  fwrite(paste0(data_dir,basename,"generate_data_single_gam_event_reduction_data.csv"))

prr_dts %>%
  fwrite(paste0(data_dir,basename,"generate_data_single_prr_event_reduction_data.csv"))

# Calculate sensitivity analysis for all ADEs n all stages -----------------------------------------------------
cat("\nEvent reporting reduction, all stages, sample of positive ades\n")

# Define global variables
d_col="D"
e_col="E"
de_col="DE"
stage_col="nichd"
x_cols=c(stage_col,"age")
all_data = raw_data
stages <- all_data[,unique(get(stage_col))]
stage_age <- raw_data[,.(nichd,age)] %>% unique()
ade_data <- make_ade_data()

# Get algorithm start time
t0 = Sys.time()
set.seed(0)
# Assign fold changes, "effect size", so ADE pairs
l=0.75
positive_ade_pairs_fc <- random_att_assign(positive_ade_pairs,l=l)
# Define stages to reduce and pick pair with many reports in that stage
stages_to_reduce <- category_levels[[stage_col]]
# Choose dynamic shapes or spikeins for drug event simulated reporting rate generation
spikein_cols = c("increase","decrease","plateau","inverse_plateau")
stage_dts <- NULL
# define indices for ADE pairs - used for testing but in use should be all pairs
inds <- sample(1:nrow(positive_ade_pairs),nrow(positive_ade_pairs),replace=F)
for(stage_to_reduce in stages_to_reduce){
  dts <- NULL
  for(spikein_col in spikein_cols){
    cat("\nPercent reduction at ",stage_to_reduce," for ",spikein_col," dynamics class\n")
    # Parallelize through each ADE pairs
    dt <- foreach(i=inds,.combine="rbind") %dopar% {
      # Set random number generator seed
      set.seed(i)
      # Get pair ids
      pair <- positive_ade_pairs[i] %>% unlist %>% unname
      # Format ADE data - get "old" event reporting vector
      dat = set_ade_pair(ade_data,pair)
      e_old <- dat[,get(e_col)]
      # Define reporting rate
      ej = mean(e_old)
      fij = positive_ade_pairs_fc[i]$fold_change
      ej_fij = fij * ej
      #define data container for augmentation
      dat_new <- dat
      # Generate event dynamic data
      stage_age$probs = define_report_probs(spikein_col,nrow(stage_age),ej_fij)
      rs <- merge(dat_new,stage_age,by="age")$probs
      # Insert event dynamic data at drug reports for simulated drug event rate
      e_new <- rbinom(nrow(dat_new),1,rs)
      dat_new[D==1,e_col] <- max_vec(e_old,e_new)[dat_new$D==1]
      dat_new[E==1 & D==1,de_col] <- 1
      e_new <- dat_new[,get(e_col)]

      lsts <- list()
      # Percentiles for data reduction
      step=0.1
      qs <- seq(0,1,step)
      # Reduce event reporting only at stage and apply detecction methods
      for(q in qs){
        # define reduced data container
        dat_q <- dat_new
        # Define event indices to sample
        stage_to_reduce_indices <- which(dat_q[,get(stage_col)==which(category_levels[[stage_col]]==stage_to_reduce)])
        stages_not_to_reduce_indices <- which(dat_q[,get(stage_col)!=which(category_levels[[stage_col]]==stage_to_reduce)])
        event_indices <- intersect(which(e_new==1),stage_to_reduce_indices)
        event_at_other_stages_indices <- intersect(which(e_new==1),stages_not_to_reduce_indices)
        # Create new event reporting vector
        e_q <- rep(0,length(e_old))
        e_q[event_at_other_stages_indices] <- 1
        # Sample reports with events per percentile rediction
        sampled_event_indices <- 
          sample(
            event_indices,
            floor(q*length(event_indices)),
            replace = F
          )
        # Insert sampled event reporting into event vector
        if(length(sampled_event_indices)>0)e_q[sampled_event_indices] <- 1
        dat_q[,e_col] <- e_q
        dat_q[E==1 & D==1,de_col] <- 1
        dat_q[E==0,de_col] <- 0
        
        # Format simulated ADE data
        dat_q$Di = dat_q[,mean(D)]
        dat_q$Ej = dat_q[,mean(E)]
        dat_q$Fij = fij
        dat_q$Tij = ej_fij
        dat_q$random_state = rs
        dat_q$random_seed = i
        dat_q[,c(drug_col) := list(pair[1])]
        dat_q[,c(rx_col) := list(pair[2])]
        # Apply ADE detection methods and attribute data reduction
        lst = compute_risk_and_process_data(pair,dat_q,ej,fij)
        for(x in names(lst)){
          lst[[x]]$percent_event_reduction = (1-q)*100
        }
        
        lsts[[as.character(q)]] <- lst
      }
      
      gam_dt <- lapply(lsts,
                       function(x){x$gam}
      ) %>%
        bind_rows() %>%
        .[nichd!="all"]

      prr_dt <- lapply(lsts,
                       function(x){x$prr}
      ) %>%
        bind_rows() %>%
        .[nichd!="all"]

      dt <-
        merge(
          prr_dt[,c(drug_col,rx_col,stage_col,"a","b","c","d","PRR","PRR_90mse","percent_event_reduction"),with=F],
          gam_dt[,c(drug_col,rx_col,stage_col,"D","E","DE","gam_score","gam_score_90mse","percent_event_reduction"),with=F],
          by=c(drug_col,rx_col,stage_col,"percent_event_reduction")
        )
      dt$ade <- paste0(dt$atc_concept_id,"_",dt$meddra_concept_id)
      dt$spikein <- spikein_col
      dt$reduced_stage <- stage_to_reduce
      dt
    }
    dts <- bind_rows(dts,dt)
    # Get algorithm end time
    t1 = Sys.time()
    cat(round(as.numeric(difftime(t1,t0,units="mins")),2)," minutes\n")
  }
  stage_dts <- bind_rows(stage_dts,dts)
  
  dts %>% 
    fwrite(paste0(data_dir,basename,"generate_",stage_to_reduce,"_data_event_reduction_data.csv"))
  
  # Get algorithm end time
  t1 = Sys.time()
  cat("\n",round(as.numeric(difftime(t1,t0,units="mins")),2)," minutes\n")
}
# Get algorithm end time
t1 = Sys.time()
cat("\n",round(as.numeric(difftime(t1,t0,units="mins")),2)," minutes\n")

stage_dts %>%
  fwrite(paste0(data_dir,basename,"generate_data_event_reduction_data.csv"))

# DEPRECATED DE reporting after fold change --------------------------
# 
# d_col="D"
# e_col="E"
# de_col="DE"
# stage_col="nichd"
# x_cols=c(stage_col,"age")
# all_data = raw_data
# stage_age = raw_data[,.(nichd,age)] %>% unique()
# seed=1
# set.seed(seed)
# spikein_col="increase"
# fcs = c(1,2)
# pair <- positive_ade_pairs[12] %>% unlist %>% unname
# 
# simulated_data <- function(dat_new,fij,spikein_col,stage_age=raw_data[,.(nichd,age)] %>% unique() ){
#   d_old <- dat_new[,get(d_col)]
#   e_old <- dat_new[,get(e_col)]
#   ej <-  mean(e_old)
#   ej_fij <- ej * fij
#   stage_age$probs <-  define_report_probs(spikein_col,nrow(stage_age),ej_fij)
#   rs <- merge(dat_new,stage_age,by="age")$probs
#   e_new <- rbinom(nrow(dat_new),1,rs)
#   dat_new[D==1,"E"] <- max_vec(e_old,e_new)[dat_new$D==1]
#   dat_new[E==1 & D==1,"DE"] <- 1
#   dat_new
# }
# 
# test <- function(pairs,fcs,spikein_col,all_data = raw_data,stage_age = raw_data[,.(nichd,age)] %>% unique(),x_cols=c(stage_col,"age"),d_col="D",e_col="E",de_col="DE",stage_col="nichd",seed=1){
#   set.seed(seed)
#   tmps <- foreach(i=1:nrow(pairs),.combine = "rbind") %dopar% {
#     pair <- pairs[i] %>% unlist %>% unname
#     dat = make_pair_data(pair)
#     tmps <- NULL
#     for(fij in fcs){
#       dat_new <- simulated_data(dat,fij,spikein_col=spikein_col)
#       tmp = dat_new[,.(DE = sum(D==1 & E==1),
#                        pD = mean(D==1),
#                        pE = mean(E==1),
#                        pDE = mean(D==1 & E==1),
#                        D = sum(D==1),
#                        E = sum(E==1)),
#                     nichd]
#       tmp$pDE_pDpE <- tmp$pDE/(tmp$pD*tmp$pE)
#       tmp$fold_change <- fij
#       tmp[,drug_col] <- pair[1]
#       tmp[,rx_col] <- pair[2]
#       tmp$dynamic <- spikein_col
#       tmps <- bind_rows(tmps,tmp)
#     }
#     tmps
#   }
#   tmps
# }
# 
# dat <- test(positive_ade_pairs[14],c(1,2,3,4,5,6,7,8,9,10),"increase")
# dat$ade <- paste0(dat$atc_concept_id,"_",dat$meddra_concept_id)
# dat$nichd <- factor(dat$nichd,labels=c("term_neonatal","infancy","toddler","early_childhood","middle_childhood","early_adolescence","late_adolescence"))
# 
# dat %>%
#   pivot_longer(cols=c("D","E","DE")) %>%
#   ggplot(aes(nichd,value,color=fold_change)) +
#   geom_point() +
#   geom_path(aes(group=fold_change)) +
#   facet_grid(name~.,scales="free") +
#   scale_colour_gradient2(low = "blue", mid="white",midpoint=5, high = "red") +
#   ylab("Number of reports\n(D==drug, E==event, DE=drug, event") +
#   xlab("")
# 
# dat <-
#   bind_rows(
#     test(positive_ade_pairs[14],c(1,2,3,4,5,6,7,8,9,10),"uniform"),
#     test(positive_ade_pairs[14],c(1,2,3,4,5,6,7,8,9,10),"increase"),
#     test(positive_ade_pairs[14],c(1,2,3,4,5,6,7,8,9,10),"decrease"),
#     test(positive_ade_pairs[14],c(1,2,3,4,5,6,7,8,9,10),"plateau"),
#     test(positive_ade_pairs[14],c(1,2,3,4,5,6,7,8,9,10),"inverse_plateau")
#   )
# dat$ade <- paste0(dat$atc_concept_id,"_",dat$meddra_concept_id)
# dat$nichd <- factor(dat$nichd,labels=c("term_neonatal","infancy","toddler","early_childhood","middle_childhood","early_adolescence","late_adolescence"))
# 
# dat[,lapply(.SD,median),.(nichd,fold_change,dynamic),.SDcols=c("D","E","DE")] %>%
#   pivot_longer(cols=c("D","E","DE")) %>%
#   ggplot(aes(nichd,value,color=fold_change)) +
#   geom_point() +
#   geom_path(aes(group=fold_change)) +
#   facet_grid(name~dynamic,scales="free") +
#   scale_colour_gradient2(low = "blue", mid="white",midpoint=5, high = "red") +
#   ylab("Median of reports across ADEs\n(D==drug, E==event, DE=drug, event") +
#   xlab("")
