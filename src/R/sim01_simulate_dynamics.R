#' Title: "Evaluating risk detection methods to uncover ontogenic-mediated adverse drug effect mechanisms in children" study
#' 
#' Author details: Nicholas Giangreco
#' 
#' This script generates simulated drug event reporting dynamics
#' and figures per simulation

# PURPOSE -----------------------------------------------------------------

#' Simulates ADE dydnamics


# Load libraries and set variables ----------------------------------------------------------

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
cores=5
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


# Define dynamics --------------------------------------------

N=10000
tij=0.03

# uniform
bprobs <- rep(tij,N)
dy <- rep(0,N)
rprobs <- bprobs+dy
plot((1:N)/N,rprobs)
mean(rprobs+dy)
plot((1:N)/N,rprobs)

# low flat then increase to flat
bprobs <- rep(tij,N)
dy <- tanh(seq(-pi,pi,length.out=N))*tij
rprobs <- bprobs+dy
plot((1:N)/N,dy)
plot((1:N)/N,rprobs)
mean(rprobs+dy)

# high flat then decrease to flat
bprobs <- rep(tij,N)
dy <- -(tanh(seq(-pi,pi,length.out=N))*tij)
rprobs <- bprobs+dy
plot((1:N)/N,dy)
plot((1:N)/N,rprobs)
mean(rprobs+dy)

# low flat then plateau then decrease to flat
bprobs <- rep(tij,N)
dy <- c((tanh(seq(-pi,pi,length.out=N/2))*tij),(tanh(seq(pi,-pi,length.out=N/2))*tij))
rprobs <- bprobs+dy
plot((1:N)/N,dy)
plot((1:N)/N,rprobs)
mean(rprobs+dy)

# high flat then plateau then increase to flat
bprobs <- rep(tij,N)
dy <- c((tanh(seq(pi,-pi,length.out=N/2))*tij),(tanh(seq(-pi,pi,length.out=N/2))*tij))
rprobs <- bprobs+dy
plot((1:N)/N,dy)
plot((1:N)/N,rprobs)
mean(rprobs+dy)


# Define dynamics function ---------------------------------------------------------


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

# Define fold change distribution ------------------------------------------------

dt <- data.table(id = 1:500,value = rexp(500,0.75)+1)

g <- dt %>% 
    ggplot(aes(value)) +
    geom_histogram(aes(y=..density..)) +
    geom_density(color="red") +
    scale_x_continuous(breaks=seq(1,14,1)) +
    scale_y_continuous(labels=scales::percent) +
    xlab("Fold change") +
    ylab("Percent of drug-event pairs")
ggsave(paste0(img_dir,basename,"fold_change_distribution.png"),g,width=5,height=4)

# Define report probabilities across childhood --------------------------

dt <- raw_data[,.(nichd)] %>% unique() 
dt$nichd <- factor(dt$nichd,levels=dt$nichd)
dt$lwr_age <- c(0,28/365,1,2,5,11,18)
dt$upr_age <- c(28/365,1,2,5,11,18,21)

dt$uniform <- define_report_probs("uniform",nrow(dt),0.3)
dt$increase <- define_report_probs("increase",nrow(dt),0.3)
dt$decrease <- define_report_probs("decrease",nrow(dt),0.3)
dt$plateau <- define_report_probs("plateau",nrow(dt),0.3)
dt$inverse_plateau <- define_report_probs("inverse_plateau",nrow(dt),0.3)

dt

g <- dt %>% 
    pivot_longer(cols=c("uniform","increase","decrease","plateau","inverse_plateau")) %>% 
    ggplot(aes(nichd,value,color=name,group=name)) +
    guides(color=guide_legend(title="Drug event dynamics class",title.position = "top")) +
    geom_line() +
    scale_color_manual(values=dynamics_colors) +
    xlab("") +
    ylab("Reporting probability") +
    theme(
        axis.text.x = element_text(angle=30,hjust=1,vjust=1),
        legend.position = "bottom"
    )
ggsave(paste0(img_dir,basename,"drug_event_reporting_across_nichd_stages.png"),g,width=10,height=7)

merged_dt <- raw_data[,.(nichd,age)] %>% unique() 
merged_dt[,.N,nichd]
merged_dt$log_age <- merged_dt[,log(age)]
merged_dt$uniform <- define_report_probs("uniform",nrow(merged_dt),0.3)
merged_dt$increase <- define_report_probs("increase",nrow(merged_dt),0.3)
merged_dt$decrease <- define_report_probs("decrease",nrow(merged_dt),0.3)
merged_dt$plateau <- define_report_probs("plateau",nrow(merged_dt),0.3)
merged_dt$inverse_plateau <- define_report_probs("inverse_plateau",nrow(merged_dt),0.3)

merged_dt
merged_dt[,lapply(.SD,mean),.SDcols=c("uniform","increase","decrease","plateau","inverse_plateau")]

g <- 
    merged_dt %>% 
    pivot_longer(cols=c("uniform","increase","decrease","plateau","inverse_plateau")) %>% 
    ggplot(aes(log_age,value,color=name)) +
    geom_line() +
    scale_x_continuous(
        label=dt$nichd,
        breaks=dt[,log(upr_age)],
        limits=c(min(dt[,log(upr_age)]),max(dt[,log(upr_age)]))
    ) +
    ylab("Report probability") +
    xlab("Age") +
    guides(color=guide_legend(title="Drug event dynamics class",title.position = "top")) +
    scale_color_manual(values=dynamics_colors) +
    geom_vline(xintercept=log(27/365)) + 
    geom_vline(xintercept=log(1)) + 
    geom_vline(xintercept=log(2)) + 
    geom_vline(xintercept=log(5)) + 
    geom_vline(xintercept=log(11)) + 
    geom_vline(xintercept=log(18)) + 
    geom_vline(xintercept=log(21)) +
    theme(
        axis.text.x = element_text(angle=30,hjust=1,vjust=1),
        legend.position = "bottom"
    ) 
ggsave(paste0(img_dir,basename,"drug_event_reporting_across_childhood.png"),g,width=15,height=7)

merged_dt[,lapply(.SD,mean),.SDcols=c("uniform","increase","decrease","plateau","inverse_plateau")]
g <- merged_dt %>% 
    pivot_longer(cols=c("uniform","increase","decrease","plateau","inverse_plateau")) %>% 
    ggplot(aes(age,value,color=name,group=name)) +
    geom_line() +
    ylab("Report probability") +
    xlab("Age") +
    scale_color_manual(values=dynamics_colors) +
    guides(color=guide_legend(title="Drug event dynamics class",title.position = "top")) +
    theme(
        legend.position = "bottom"
    ) 
ggsave(paste0(img_dir,basename,"drug_event_reporting_across_untransformed_childhood.png"),g,width=15,height=7)

dt_range <- dt[,.(seq = list(seq(lwr_age,upr_age,length.out=100))),.(nichd,lwr_age,upr_age)]
dt_range <- dt_range[,.(age = unlist(seq)),.(nichd,lwr_age,upr_age)]
dt_range$log_age <- dt_range[,log(age)]
dt_range$uniform <- define_report_probs("uniform",nrow(dt_range),0.3)
dt_range$increase <- define_report_probs("increase",nrow(dt_range),0.3)
dt_range$decrease <- define_report_probs("decrease",nrow(dt_range),0.3)
dt_range$plateau <- define_report_probs("plateau",nrow(dt_range),0.3)
dt_range$inverse_plateau <- define_report_probs("inverse_plateau",nrow(dt_range),0.3)

g <- 
    dt_range %>% 
    pivot_longer(cols=c("uniform","increase","decrease","plateau","inverse_plateau")) %>% 
    ggplot(aes(log_age,value,color=name)) +
    geom_line() +
    scale_x_continuous(
        label=dt$nichd,
        breaks=dt[,log(upr_age)],
        limits=c(min(dt[,log(upr_age)]),max(dt[,log(upr_age)]))
    ) +
    ylab("Report probability") +
    xlab("Age") +
    scale_color_manual(values=dynamics_colors) +
    guides(color=guide_legend(title="Drug event dynamics class",title.position = "top")) +
    geom_vline(xintercept=log(0)) + 
    geom_vline(xintercept=log(27/365)) + 
    geom_vline(xintercept=log(1)) + 
    geom_vline(xintercept=log(2)) + 
    geom_vline(xintercept=log(5)) + 
    geom_vline(xintercept=log(11)) + 
    geom_vline(xintercept=log(18)) + 
    geom_vline(xintercept=log(21)) +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle=30,vjust=1,hjust=1)
    )
ggsave(paste0(img_dir,basename,"drug_event_reporting_across_childhood_smooth.png"),g,width=10,height=7)

g <- dt_range %>% 
    pivot_longer(cols=c("uniform","increase","decrease","plateau","inverse_plateau")) %>% 
    ggplot(aes(age,value,color=name,group=name)) +
    guides(color=guide_legend(title="dynamics class")) +
    geom_line() +
    ylab("Report probability") +
    xlab("Age") +
    scale_color_manual(values=dynamics_colors) +
    guides(color=guide_legend(title="Drug event dynamics class",title.position = "top")) +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle=30,vjust=1,hjust=1)
    )
ggsave(paste0(img_dir,basename,"drug_event_reporting_across_untransformed_childhood_smooth.png"),g,width=10,height=7)

