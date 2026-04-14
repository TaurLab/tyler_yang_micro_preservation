# LIBRARIES ------------------------
library(phyloseq)
library(yingtools2)
library(phytools)
library(tidyverse)
library(ggtree)
library(lattice)
library(grid)
library(gridExtra)
library(ggfittext)
library(ggpubr)
library(ggnewscale)
library(vegan)
library(readxl)
library(glue)
library(vegan)
library(forcats)
library(ggh4x)

# DATA IMPORT ---------------
rm(list = ls())
setwd("/Users/Tyler/Documents/Research/IntegrityProject/tyler_yang_micro_preservation")  # change path name if working in a different directory

load("data/phy.tyler.RData")
load("data/phy.tyler.additional.RData")

# Reset the samp
phy.tyler = phy.tyler %>% 
  select(otu,Superkingdom,Phylum,Class,Order,Family,Genus,Species)

s.tyler = phy.tyler %>% get.samp() %>%
  left_join(trace_tbl.tyler,by="sample") %>%
  mutate(temp=ifelse(sample %in% c("1A","1B"),"-80C",temp),
         temp=factor(temp,levels=c("-80C", "-20C", "4C", "room temp")),
         time=factor(time,levels=c("day 0", "day 3", "day 6", "day 8", "day 9", "day 11")),
         treatment=factor(treatment,levels=c("none", "75C", "UV", "75C+UV", "autoclave", "autoclave+UV", "UV DNA")),
         sample_number=order(order(treatment,temp,time)),
         letter=str_extract(sample,"(?<=^[0-9]{1,2})[AB]"),
         sample2=paste(treatment,temp,storage,time,sample_number,sep="|"),
         sample2=fct_reordern(sample2,sample_number),
         sample.comparator=ifelse(experiment==1,"1A","TY.1_D0_NT"))
sample_data(phy.tyler) <- s.tyler %>% set.samp()

# Functions that really just condense ugly code for dataframe cleaning
add_clean_names = function(phy_samp_df) {
  cleaned = phy_samp_df %>%
    mutate(sampname = sample,
           samplabel = sample) %>%
    # Manual renaming of special controls
    mutate(sampname = recode(sampname,
                             "1A" = "Control 1",
                             "1B" = "Control 2",
                             "TY.1_D0_NT" = "Control 3",
                             "TY.19_D0_DNA_UV" = "D0 UV post-extraction",
                             "TY.20_D6_DNA_UV" = "D6 UV post-extraction",
                             "TY.21_D9_DNA_UV" = "D9 UV post-extraction",
                             "TY.8_D6_AC" = "D6 autoclave",
                             "TY.9_D9_AC" = "D9 autoclave",
                             "TY.7_D0_AC" = "D0 autoclave",
                             "TY.16_D0_AC_UV" = "D0 autoclave + UV pre-ex",
                             "TY.17_D6_AC_UV" = "D6 autoclave + UV pre-ex",
                             "TY.18_D9_AC_UV" = "D9 autoclave + UV pre-ex",
                             "TY.2_D6_Dry" = "D6 dried",
                             "TY.3_D9_Dry" = "D9 dried",
                             "TY.4_D0_75C" = "D0 75C heat",
                             "TY.5_D6_75C" = "D6 75C heat",
                             "TY.6_D9_75C" = "D9 75C heat",
                             "TY.10_D0_UV" = "D0 UV pre-extraction",
                             "TY.11_D6_UV" = "D6 UV pre-extraction",
                             "TY.12_D9_UV" = "D9 UV pre-extraction",
                             "TY.13_D0_75C_UV" = "D0 75C + UV pre-ex",
                             "TY.14_D6_75C_UV" = "D6 75C + UV pre-ex",
                             "TY.15_D9_75C_UV" = "D9 75C + UV pre-ex",
                             "2A.RT.D3" = "D3 at 20C A",
                             "2B.RT.D3" = "D3 at 20C B",
                             "3A.4C.D3" = "D3 at 4C A",
                             "3B.4C.D3" = "D3 at 4C B",
                             "4A.-20.D3" = "D3 at -20C A",
                             "4B.-20.D3" = "D3 at -20C B",
                             "5A.-80.D3" = "D3 at -80C A",
                             "5B.-80.D3" = "D3 at -80C B",
                             "6A.RT.D8" = "D8 at 20C A",
                             "6B.RT.D8" = "D8 at 20C B",
                             "7A.4C.D8" = "D8 at 4C A",
                             "7B.4C.D8" = "D8 at 4C B",
                             "8A.-20.D8" = "D8 at -20C A",
                             "8B.-20.D8" = "D8 at -20C B",
                             "9A.-80.D8" = "D8 at -80C A",
                             "9B.-80.D8" = "D8 at -80C B",
                             "10A.RT.D11" = "D11 at 20C A",
                             "10B.RT.D11" = "D11 at 20C B",
                             "11A.4C.D11" = "D11 at 4C A",
                             "11B.4C.D11" = "D11 at 4C B",
                             "12A.-20.D11" = "D11 at -20C A",
                             "12B.-20.D11" = "D11 at -20C B",
                             "13A.-80.D11" = "D11 at -80C A",
                             "13B.-80.D11" = "D11 at -80C B"
    ),
    samplabel = recode(samplabel,
                       "1A" = "T1",
                       "1B" = "T2",
                       "TY.1_D0_NT" = "TY2",
                       "TY.19_D0_DNA_UV" = "TY19",
                       "TY.20_D6_DNA_UV" = "TY20",
                       "TY.21_D9_DNA_UV" = "TY21",
                       "TY.8_D6_AC" = "TY8",
                       "TY.9_D9_AC" = "TY9",
                       "TY.7_D0_AC" = "TY7",
                       "TY.16_D0_AC_UV" = "TY16",
                       "TY.17_D6_AC_UV" = "TY17",
                       "TY.18_D9_AC_UV" = "TY18",
                       "TY.2_D6_Dry" = "TY2",
                       "TY.3_D9_Dry" = "TY3",
                       "TY.4_D0_75C" = "TY4",
                       "TY.5_D6_75C" = "TY5",
                       "TY.6_D9_75C" = "TY6",
                       "TY.10_D0_UV" = "TY10",
                       "TY.11_D6_UV" = "TY11",
                       "TY.12_D9_UV" = "TY12",
                       "TY.13_D0_75C_UV" = "TY13",
                       "TY.14_D6_75C_UV" = "TY14",
                       "TY.15_D9_75C_UV" = "TY15",
                       "2A.RT.D3" = "T3",
                       "2B.RT.D3" = "T4",
                       "3A.4C.D3" = "T5",
                       "3B.4C.D3" = "T6",
                       "4A.-20.D3" = "T7",
                       "4B.-20.D3" = "T8",
                       "5A.-80.D3" = "T9",
                       "5B.-80.D3" = "T10",
                       "6A.RT.D8" = "T11",
                       "6B.RT.D8" = "T12",
                       "7A.4C.D8" = "T13",
                       "7B.4C.D8" = "T14",
                       "8A.-20.D8" = "T15",
                       "8B.-20.D8" = "T16",
                       "9A.-80.D8" = "T17",
                       "9B.-80.D8" = "T18",
                       "10A.RT.D11" = "T19",
                       "10B.RT.D11" = "T20",
                       "11A.4C.D11" = "T21",
                       "11B.4C.D11" = "T22",
                       "12A.-20.D11" = "T23",
                       "12B.-20.D11" = "T24",
                       "13A.-80.D11" = "T25",
                       "13B.-80.D11" = "T26"))
  return(cleaned)
}
create_samp_df = function(phy_object) {
  qPCR = read.csv("/Users/Tyler/Documents/Research/IntegrityProject/TY-qPCR.csv", stringsAsFactors = FALSE) %>%
    filter(Dilution.Number != "#N/A") %>%
    select(Sample.ID:Average.CT.Value) %>%
    slice(1:44) %>%
    rename('samplabel' = 'Sample.ID', 
           'seqs_total' = 'X16S.Adjusted') %>%
    mutate(seqs_total = as.numeric(seqs_total))
  
  allsamps = phy_object %>%
    get.samp(stats = TRUE) %>%
    add_clean_names() %>%
    mutate(
      label = if_else(experiment < 1.5, "one", "two"),
      sample = fct_relevel(sample,
                           "1A", "1B", "2A.RT.D3", "2B.RT.D3", "3A.4C.D3", "3B.4C.D3",
                           "4A.-20.D3", "4B.-20.D3", "5A.-80.D3", "5B.-80.D3", "6A.RT.D8",
                           "6B.RT.D8", "7A.4C.D8", "7B.4C.D8", "8A.-20.D8", "8B.-20.D8",
                           "9A.-80.D8", "9B.-80.D8", "10A.RT.D11", "10B.RT.D11", "11A.4C.D11",
                           "11B.4C.D11", "12A.-20.D11", "12B.-20.D11", "13A.-80.D11", "13B.-80.D11",
                           "TY.1_D0_NT", "TY.2_D6_Dry", "TY.3_D9_Dry", "TY.4_D0_75C", "TY.5_D6_75C",
                           "TY.6_D9_75C", "TY.10_D0_UV", "TY.11_D6_UV", "TY.12_D9_UV", "TY.13_D0_75C_UV",
                           "TY.14_D6_75C_UV", "TY.15_D9_75C_UV", "TY.7_D0_AC", "TY.8_D6_AC", "TY.9_D9_AC",
                           "TY.16_D0_AC_UV", "TY.17_D6_AC_UV", "TY.18_D9_AC_UV", "TY.19_D0_DNA_UV",
                           "TY.20_D6_DNA_UV", "TY.21_D9_DNA_UV"),
      sampname = fct_relevel(sampname,
                             "Control 1", "Control 2", "D3 at 20C A", "D3 at 20C B", "D3 at 4C A", "D3 at 4C B",
                             "D3 at -20C A", "D3 at -20C B", "D3 at -80C A", "D3 at -80C B", "D8 at 20C A",
                             "D8 at 20C B", "D8 at 4C A", "D8 at 4C B", "D8 at -20C A", "D8 at -20C B",
                             "D8 at -80C A", "D8 at -80C B", "D11 at 20C A", "D11 at 20C B", "D11 at 4C A",
                             "D11 at 4C B", "D11 at -20C A", "D11 at -20C B", "D11 at -80C A", "D11 at -80C B",
                             "Control 3", "D6 dried", "D9 dried", "D0 UV pre-extraction", "D6 UV pre-extraction", "D9 UV pre-extraction",
                             "D0 75C heat", "D6 75C heat", "D9 75C heat", "D0 75C + UV pre-ex", "D6 75C + UV pre-ex", "D9 75C + UV pre-ex", 
                             "D0 autoclave", "D6 autoclave", "D9 autoclave", "D0 autoclave + UV pre-ex", "D6 autoclave + UV pre-ex",
                             "D9 autoclave + UV pre-ex", "D0 UV post-extraction", "D6 UV post-extraction", "D9 UV post-extraction"),
      heat = if_else(experiment == 1, "no heat", heat),
      treatment = if_else(experiment == 1, "none", treatment),
      origin = case_when(
        is.na(experiment) ~ "Other",
        experiment == 1 ~ "Year 1 Samples",
        heat == "autoclave" ~ "Autoclave",
        uv == "UV DNA" ~ "UV DNA",
        TRUE ~ "Year 2 Samples")) %>%
    left_join(qPCR, by = "samplabel") 
  return(allsamps)
}
create_otu_df = function(phy_object, phy_samp_df) {
  physub = phy_object %>% 
    phy.collapse()  # combine OTUs of the same taxonomic classification
  otusub = physub %>% 
    get.otu.melt(sample_data = FALSE) %>%
    tax.plot(data = TRUE) %>%  
    left_join(phy_samp_df) %>%
    mutate(absseqs = pctseqs * seqs_total)
}


# Samp and OTU dataframes 
allsamps = create_samp_df(phy.tyler)  # get SAMP dataframe
otusub = create_otu_df(phy.tyler, allsamps)  # get OTU dataframe (with mapping/samp data)

firstsamps = allsamps %>%
  dplyr::filter(experiment == 1) %>%
  mutate(replicate = if_else(grepl("A", sample), 1, 2))  # for Figure 1 pairing
secondsamps = allsamps %>%
  dplyr::filter(experiment == 2)

firstotu = otusub %>% 
  dplyr::filter(experiment == 1) %>%
  mutate(replicate = if_else(grepl("A", sample), 1, 2)) 
secondotu = otusub %>%
  dplyr::filter(experiment == 2)

# Linear Regression (Exp 1) ------

add_dist1 = function(sampdata,method,sample0="1A",phy=phy1) {
  varname <- paste0("dist_",method)
  dist <- phy.tyler %>%
    calc.distance(method=method)  
  pw <- dist %>% get.pairwise() %>%
    filter(sample1==sample0|sample2==sample0) %>%
    mutate(sample=coalesce(na_if(sample1,sample0),sample2)) %>%
    select(sample,!!varname:=dist)  
  sampdata %>% left_join(pw,by="sample")
}
testit1 = function(sampdata,var) {
  var <- ensym(var)
  findp <- function(models,text) {
    models %>% map_lgl(~{
      pvals <- .x %>% filter(term %ilike% text) %>% pull(p.value)
      if (length(pvals)==0) {
        NA
      } else {
        any(pvals<=0.05)
      }
    })
  }
  tests <- rlang::inject(
    list(
      !!var ~ time,
      !!var ~ time.num,
      !!var ~ time.rank,
      !!var ~ temp,
      !!var ~ temp.num,
      !!var ~ temp.rank,
      !!var ~ time+temp,
      !!var ~ time.num+temp.num,
      !!var ~ time.rank+temp.rank
    )
  )
  tbl <- tibble(formula=tests) %>%
    mutate(test=map_chr(formula,deparse1),
           model=map(formula,~{
             broom::tidy(lm(.x,data=sampdata))
           }),
           temp=findp(model,"temp"),
           time=findp(model,"time"))
  return(tbl)
}
plotit1 = function(sampdata,var,eps=0,ymax=NULL) {
  if (!is.null(ymax)) {
    ymax <- expand_limits(y=ymax)
  }
  var <- enquo(var)
  ggplot(sampdata) +
    ymax +
    geom_col(aes(x=letter,y=!!var,
                 color=baseline),linetype="longdash",linewidth=0.75,
             show.legend=FALSE) +
    scale_color_manual(values=c("TRUE"="blue","FALSE"=NA)) +
    scale_y_continuous(trans=log_epsilon_trans(epsilon = eps)) +
    facet_grid(temp ~ time) +
    theme(aspect.ratio=2)
}
export_stats_one = function(samp, variable_list) {
  temp_time = 7
  result = tibble()
  for (variable in variable_list) {
    print(variable)
    test_output = testit1(samp, !!sym(variable))
    model_output = test_output$model[[temp_time]] %>%
      mutate(metric = variable) %>%
      select(metric, everything())
    result = rbind(result, model_output)
  }
  return(result)
}

phy1 = phy.tyler %>% filter(experiment==1) %>%
  mutate(baseline=sample=="1A") %>%
  add.abundance(pct.anaerobe=Class=="Clostridia" |
                  Class=="Negativicutes" |
                  Phylum=="Fusobacteria" |
                  Phylum=="Bacteroidetes",denom=TRUE)

s1 = phy1 %>% 
  get.samp(stats=TRUE) %>%
  add_dist1("mean.bray") %>%
  add_dist1("horn") %>%
  add_dist1("mean.horn") %>%
  add_dist1("unfold.horn") %>%
  mutate(days=as.numeric(str_replace(time,"day ","")),
         celsius=ifelse(temp=="room temp",20,as.numeric(str_replace(temp,"C",""))),
         time.num=days,
         time.rank=dense_rank(days),
         temp.num=celsius,
         temp.rank=dense_rank(temp.num))
s1_test = s1 %>%
  mutate(time = factor(time, levels = c("day 0", "day 3", "day 8", "day 11")),
         temp = factor(temp, levels = c("-80C", "-20C", "4C", "room temp")))

variable_list = list('InvSimpson', 'Shannon',
                     "dist_mean.bray", "dist_horn", "dist_mean.horn", "dist_unfold.horn",
                     'qpcr.totalseqs', 'pct.anaerobe')

one_linreg_stats = export_stats_one(s1_test, variable_list) %>%
  mutate(significant = if_else((p.value < 0.05), TRUE, FALSE))
one_linreg_sig = one_linreg_stats %>%
  filter(significant == TRUE) %>%
  filter(term != "(Intercept)") %>%
  select(metric, term, p.value, significant)

one_linreg_sig %>% dt

one_linreg_horn = one_linreg_stats %>%
  filter(metric == "dist_mean.horn")
  
  
write.csv(one_linreg_stats, '/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/Linear Regression Tables/experiment_one_stats.csv')
write.csv(one_linreg_sig, '/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/Linear Regression Tables/experiment_one_significant.csv')

# more under the hood
one_output = testit1(s1_test, dist_mean.bray)
one_output$model[7]

plotit1(s1,InvSimpson)
plotit1(s1,qpcr.totalseqs,eps=100)
plotit1(s1,dist_horn,ymax=1)
plotit1(s1,dist_mean.horn,ymax=1)
plotit1(s1,dist_unfold.horn,ymax=1)

# Linear Regression Combined Condition (Exp ---------
s1_combcond = s1 %>%
  mutate(time = factor(time, levels = c("day 0", "day 3", "day 8", "day 11")),
         temp = factor(temp, levels = c("-80C", "-20C", "4C", "room temp"))) %>%
  mutate(condition = factor(paste(temp, time, sep = "_"))) %>%
  # Force the -80C_day 0 to be the reference (Intercept)
  mutate(condition = fct_relevel(condition, "-80C_day 0"))

# 2. Update testit1 to include the condition formula
testit1_combcond = function(sampdata,var) {
  var <- ensym(var)
  findp <- function(models,text) {
    models %>% map_lgl(~{
      pvals <- .x %>% filter(term %ilike% text) %>% pull(p.value)
      if (length(pvals)==0) {
        NA
      } else {
        any(pvals<=0.05)
      }
    })
  }
  tests <- rlang::inject(list(!!var ~ time,
      !!var ~ time.num,
      !!var ~ time.rank,
      !!var ~ temp,
      !!var ~ temp.num,
      !!var ~ temp.rank,
      !!var ~ time+temp,
      !!var ~ time.num+temp.num,
      !!var ~ time.rank+temp.rank,
      !!var ~ condition # NEW: This is the 10th test in the list
    ))
  tbl <- tibble(formula=tests) %>%
    mutate(test=map_chr(formula,deparse1),
           model=map(formula,~{
             broom::tidy(lm(.x,data=sampdata))
           }),
           temp=findp(model,"temp"),
           time=findp(model,"time"))
  return(tbl)
}

# 3. Update export_stats_one to extract the 10th model (condition)
export_stats_combcond = function(samp, variable_list) {
  temp_time = 10 # Extracts the !!var ~ condition model
  result = tibble()
  for (variable in variable_list) {
    print(variable)
    test_output = testit1_combcond(samp, !!sym(variable))
    model_output = test_output$model[[temp_time]] %>%
      mutate(metric = variable) %>%
      select(metric, everything())
    result = rbind(result, model_output)
  }
  return(result)
}

variable_list = list('InvSimpson', 'Shannon',
                     "dist_mean.bray", "dist_horn", "dist_mean.horn", "dist_unfold.horn",
                     'qpcr.totalseqs', 'pct.anaerobe')

# 4. Run the stats, add p.adj, and filter for significance
combcond_linreg_stats = export_stats_combcond(s1_combcond, variable_list) %>%
  # Grouping by metric ensures the FDR correction is applied per metric type
  group_by(metric) %>% 
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  # Determine significance based on the adjusted p-value instead of raw p-value
  mutate(significant = if_else((p.adj < 0.05), TRUE, FALSE))

combcond_linreg_sig = combcond_linreg_stats %>%
  filter(significant == TRUE) %>%
  filter(term != "(Intercept)") %>%
  select(metric, term, p.value, p.adj, significant)

combcond_linreg_sig %>% dt

write.csv(combcond_linreg_stats, '/Users/Tyler/Documents/Research/IntegrityProject/tyler_yang_micro_preservation/stats/experiment_one_comb_stats.csv')
write.csv(combcond_linreg_sig, '/Users/Tyler/Documents/Research/IntegrityProject/tyler_yang_micro_preservation/stats/experiment_one_comb_significant.csv')

# Linear Regression (Exp 2) ------

add_dist2 = function(sampdata,method,sample0="TY.1_D0_NT",phy=phy2) {
  varname <- paste0("dist_",method)
  dist <- phy.tyler %>%
    calc.distance(method=method)  
  pw <- dist %>% get.pairwise() %>%
    filter(sample1==sample0|sample2==sample0) %>%
    mutate(sample=coalesce(na_if(sample1,sample0),sample2)) %>%
    select(sample,!!varname:=dist)  
  sampdata %>% left_join(pw,by="sample")
}
testit2 = function(sampdata,var) {
  var <- ensym(var)
  tests <- rlang::inject(
    list(
      !!var ~ time,
      !!var ~ time.num,
      !!var ~ time.rank,
      !!var ~ uv,
      !!var ~ uv.dna,
      !!var ~ heat,
      !!var ~ heat.autoclave,
      !!var ~ time+uv+heat,
      !!var ~ time.num+uv+heat,
      !!var ~ time.rank+uv+heat,
      !!var ~ time.rank+uv.dna+heat.autoclave
    )
  )
  
  findp <- function(models,text) {
    models %>% map_lgl(~{
      pvals <- .x %>% filter(term %ilike% text) %>% pull(p.value)
      if (length(pvals)==0) {
        NA
      } else {
        any(pvals<=0.05)
      }
    })
  }
  tbl <- tibble(formula=tests) %>%
    mutate(test=map_chr(formula,deparse1),
           model=map(formula,~{broom::tidy(lm(.x,data=sampdata))}),
           time=findp(model,"time"),
           heat.75c=findp(model,"heat75C"),
           heat.autoclave=findp(model,"autoclave"),
           uv.regular=findp(model,"uvUV$"),
           uv.dna=findp(model,"dna"))
  return(tbl)
}
plotit2 = function(sampdata,var,eps=0,ymax=NULL) {
  var <- enquo(var)
  
  if (!is.null(ymax)) {
    ymax <- expand_limits(y=ymax)
  }
  ggplot(sampdata) + 
    ymax +
    geom_col(aes(x=factor(1),y=!!var,
                 color=baseline),linetype="longdash",linewidth=0.75,
             show.legend=FALSE) +
    scale_color_manual(values=c("TRUE"="blue","FALSE"=NA)) +
    facet_grid(treatment~time) + 
    scale_y_continuous(trans=log_epsilon_trans(eps)) +
    theme(aspect.ratio=2)
}


phy2 = phy.tyler %>% filter(experiment==2) %>%
  mutate(baseline=sample=="TY.1_D0_NT")
s2 = phy2 %>% get.samp(stats=TRUE) %>%
  add_dist2("mean.bray") %>%
  add_dist2("horn") %>%
  add_dist2("mean.horn") %>%
  add_dist2("unfold.horn") %>%
  mutate(time.num=as.numeric(str_replace(time,"day ","")),
         time.rank=dense_rank(time.num),
         uv=fct_relevel(uv,"no UV"),
         heat=fct_relevel(heat,"no heat"),
         heat.autoclave=heat=="autoclave",
         uv.dna=uv=="UV DNA",
         qpcr.totalseqs=coalesce(qpcr.totalseqs,100),
         log.qpcr.totalseqs=log(qpcr.totalseqs))

export_stats_two = function(samp, variable_list) {
  all_treatments = 8
  result = tibble()
  for (variable in variable_list) {
    print(variable)
    test_output = testit2(samp, !!sym(variable))
    model_output = test_output$model[[all_treatments]] %>%
      mutate(metric = variable) %>%
      select(metric, everything())
    result = rbind(result, model_output)
  }
  return(result)
}

variable_list_two = list('InvSimpson', 'Shannon',
                         "dist_mean.bray", "dist_horn", "dist_mean.horn", "dist_unfold.horn",
                         'qpcr.totalseqs')
two_linreg_stats = export_stats_two(s2, variable_list_two) %>%
  mutate(significant = if_else((p.value < 0.05), TRUE, FALSE)) %>%
  filter(term != "(Intercept)") %>%
  arrange(metric, term) %>%
  select(metric, term, everything())
two_linreg_sig = two_linreg_stats %>%
  filter(significant == TRUE) 
two_linreg_sig %>% dt

write.csv(two_linreg_stats, '/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/Linear Regression Tables/experiment_two_stats.csv')
write.csv(two_linreg_sig, '/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/Linear Regression Tables/experiment_two_significant.csv')


# PERMANOVA (Exp 1) ----------------

comm_df = firstotu %>%
  select(sample, otu, pctseqs) %>%   # or pctseqs if relative
  pivot_wider(
    names_from = otu,
    values_from = pctseqs,
    values_fill = 0)

comm_mat = as.data.frame(comm_df)

rownames(comm_mat) = comm_mat$sample
comm_mat$sample = NULL

comm_mat = as.matrix(comm_mat)
metadata = firstotu %>%
  select(sample, temp, time) %>%
  distinct()
metadata = metadata[match(rownames(comm_mat), metadata$sample), ]

dist_horn = vegdist(comm_mat, method = "horn")
dist_bray = vegdist(comm_mat, method = "bray")

# HORN
# Test for a significant interaction between temp and time
# ex. room temperature after 11 days matters 
adonis2(dist_horn ~ temp*time,
  data = metadata,
  permutations = 999,
  by = "margin")
  # Test for differences across temp and time groups separately 
adonis2(dist_horn ~ temp + time,
  data = metadata,
  permutations = 999,
  by = "margin")

# Test dispersion by temperature — differences in variability within samples of a given temperature group  
bd_horn = betadisper(dist_horn, metadata$temp)
anova(bd_horn)
permutest(bd_horn)
plot(bd_horn)

# BRAY
adonis2(dist_bray ~ temp*time,
        data = metadata,
        permutations = 999,
        by = "margin")

# Test for differences across temp and time groups separately 
adonis2(dist_bray ~ temp + time,
        data = metadata,
        permutations = 999,
        by = "margin")

# Test dispersion by temperature — differences in variability within samples of a given temperature group  
bd_bray = betadisper(dist_bray, metadata$temp)
anova(bd_bray)
permutest(bd_bray)
plot(bd_bray)


# PERMANOVA Pairwise Time (Exp 1) ------------

metadata = metadata[match(rownames(comm_mat), metadata$sample), , drop = FALSE]

# make time a factor in the order you want
metadata$time = factor(metadata$time, levels = c("day 0", "day 3", "day 8", "day 11"))

# all pairwise time comparisons
time_pairs = combn(levels(metadata$time), 2, simplify = FALSE)

pairwise_results = lapply(time_pairs, function(tp) {
  keep <- metadata$time %in% tp
  
  md_sub <- metadata[keep, , drop = FALSE]
  comm_sub <- comm_mat[keep, , drop = FALSE]
  md_sub$time <- droplevels(md_sub$time)
  
  dist_sub <- vegdist(comm_sub, method = "horn")
  fit <- adonis2(dist_sub ~ time, data = md_sub, permutations = 999)
  data.frame(time1 = tp[1],
             time2 = tp[2],
             F = fit$F[1],
             R2 = fit$R2[1],
             p = fit$`Pr(>F)`[1])
})

pairwise_time = bind_rows(pairwise_results) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(p_adj)

pairwise_time  # differences day 3 v day 8 | day 8 v day 11


# PRE-FIGURE Reset -----------------

rm(list=ls())

load("data/phy.tyler.RData")
load("data/phy.tyler.additional.RData")
phy.others = read_rds("data/other.samps.rds")

add_clean_names = function(phy_samp_df) {
  cleaned = phy_samp_df %>%
    mutate(sampname = sample,
           samplabel = sample) %>%
    # Manual renaming of special controls
    mutate(sampname = recode(sampname,
                             "1A" = "Control 1",
                             "1B" = "Control 2",
                             "TY.1_D0_NT" = "Control 3",
                             "TY.19_D0_DNA_UV" = "D0 UV post-extraction",
                             "TY.20_D6_DNA_UV" = "D6 UV post-extraction",
                             "TY.21_D9_DNA_UV" = "D9 UV post-extraction",
                             "TY.8_D6_AC" = "D6 autoclave",
                             "TY.9_D9_AC" = "D9 autoclave",
                             "TY.7_D0_AC" = "D0 autoclave",
                             "TY.16_D0_AC_UV" = "D0 autoclave + UV pre-ex",
                             "TY.17_D6_AC_UV" = "D6 autoclave + UV pre-ex",
                             "TY.18_D9_AC_UV" = "D9 autoclave + UV pre-ex",
                             "TY.2_D6_Dry" = "D6 dried",
                             "TY.3_D9_Dry" = "D9 dried",
                             "TY.4_D0_75C" = "D0 75C heat",
                             "TY.5_D6_75C" = "D6 75C heat",
                             "TY.6_D9_75C" = "D9 75C heat",
                             "TY.10_D0_UV" = "D0 UV pre-extraction",
                             "TY.11_D6_UV" = "D6 UV pre-extraction",
                             "TY.12_D9_UV" = "D9 UV pre-extraction",
                             "TY.13_D0_75C_UV" = "D0 75C + UV pre-ex",
                             "TY.14_D6_75C_UV" = "D6 75C + UV pre-ex",
                             "TY.15_D9_75C_UV" = "D9 75C + UV pre-ex",
                             "2A.RT.D3" = "D3 at 20C A",
                             "2B.RT.D3" = "D3 at 20C B",
                             "3A.4C.D3" = "D3 at 4C A",
                             "3B.4C.D3" = "D3 at 4C B",
                             "4A.-20.D3" = "D3 at -20C A",
                             "4B.-20.D3" = "D3 at -20C B",
                             "5A.-80.D3" = "D3 at -80C A",
                             "5B.-80.D3" = "D3 at -80C B",
                             "6A.RT.D8" = "D8 at 20C A",
                             "6B.RT.D8" = "D8 at 20C B",
                             "7A.4C.D8" = "D8 at 4C A",
                             "7B.4C.D8" = "D8 at 4C B",
                             "8A.-20.D8" = "D8 at -20C A",
                             "8B.-20.D8" = "D8 at -20C B",
                             "9A.-80.D8" = "D8 at -80C A",
                             "9B.-80.D8" = "D8 at -80C B",
                             "10A.RT.D11" = "D11 at 20C A",
                             "10B.RT.D11" = "D11 at 20C B",
                             "11A.4C.D11" = "D11 at 4C A",
                             "11B.4C.D11" = "D11 at 4C B",
                             "12A.-20.D11" = "D11 at -20C A",
                             "12B.-20.D11" = "D11 at -20C B",
                             "13A.-80.D11" = "D11 at -80C A",
                             "13B.-80.D11" = "D11 at -80C B"
    ),
    samplabel = recode(samplabel,
                       "1A" = "T1",
                       "1B" = "T2",
                       "TY.1_D0_NT" = "TY2",
                       "TY.19_D0_DNA_UV" = "TY19",
                       "TY.20_D6_DNA_UV" = "TY20",
                       "TY.21_D9_DNA_UV" = "TY21",
                       "TY.8_D6_AC" = "TY8",
                       "TY.9_D9_AC" = "TY9",
                       "TY.7_D0_AC" = "TY7",
                       "TY.16_D0_AC_UV" = "TY16",
                       "TY.17_D6_AC_UV" = "TY17",
                       "TY.18_D9_AC_UV" = "TY18",
                       "TY.2_D6_Dry" = "TY2",
                       "TY.3_D9_Dry" = "TY3",
                       "TY.4_D0_75C" = "TY4",
                       "TY.5_D6_75C" = "TY5",
                       "TY.6_D9_75C" = "TY6",
                       "TY.10_D0_UV" = "TY10",
                       "TY.11_D6_UV" = "TY11",
                       "TY.12_D9_UV" = "TY12",
                       "TY.13_D0_75C_UV" = "TY13",
                       "TY.14_D6_75C_UV" = "TY14",
                       "TY.15_D9_75C_UV" = "TY15",
                       "2A.RT.D3" = "T3",
                       "2B.RT.D3" = "T4",
                       "3A.4C.D3" = "T5",
                       "3B.4C.D3" = "T6",
                       "4A.-20.D3" = "T7",
                       "4B.-20.D3" = "T8",
                       "5A.-80.D3" = "T9",
                       "5B.-80.D3" = "T10",
                       "6A.RT.D8" = "T11",
                       "6B.RT.D8" = "T12",
                       "7A.4C.D8" = "T13",
                       "7B.4C.D8" = "T14",
                       "8A.-20.D8" = "T15",
                       "8B.-20.D8" = "T16",
                       "9A.-80.D8" = "T17",
                       "9B.-80.D8" = "T18",
                       "10A.RT.D11" = "T19",
                       "10B.RT.D11" = "T20",
                       "11A.4C.D11" = "T21",
                       "11B.4C.D11" = "T22",
                       "12A.-20.D11" = "T23",
                       "12B.-20.D11" = "T24",
                       "13A.-80.D11" = "T25",
                       "13B.-80.D11" = "T26"))
  return(cleaned)
}
create_samp_df = function(phy_object) {
  qPCR = read.csv("/Users/Tyler/Documents/Research/IntegrityProject/TY-qPCR.csv", stringsAsFactors = FALSE) %>%
    filter(Dilution.Number != "#N/A") %>%
    select(Sample.ID:Average.CT.Value) %>%
    slice(1:44) %>%
    rename('samplabel' = 'Sample.ID', 
           'seqs_total' = 'X16S.Adjusted') %>%
    mutate(seqs_total = as.numeric(seqs_total))
  
  allsamps = phy_object %>%
    get.samp(stats = TRUE) %>%
    add_clean_names() %>%
    mutate(
      label = if_else(experiment < 1.5, "one", "two"),
      sample = fct_relevel(sample,
                           "1A", "1B", "2A.RT.D3", "2B.RT.D3", "3A.4C.D3", "3B.4C.D3",
                           "4A.-20.D3", "4B.-20.D3", "5A.-80.D3", "5B.-80.D3", "6A.RT.D8",
                           "6B.RT.D8", "7A.4C.D8", "7B.4C.D8", "8A.-20.D8", "8B.-20.D8",
                           "9A.-80.D8", "9B.-80.D8", "10A.RT.D11", "10B.RT.D11", "11A.4C.D11",
                           "11B.4C.D11", "12A.-20.D11", "12B.-20.D11", "13A.-80.D11", "13B.-80.D11",
                           "TY.1_D0_NT", "TY.2_D6_Dry", "TY.3_D9_Dry", "TY.4_D0_75C", "TY.5_D6_75C",
                           "TY.6_D9_75C", "TY.10_D0_UV", "TY.11_D6_UV", "TY.12_D9_UV", "TY.13_D0_75C_UV",
                           "TY.14_D6_75C_UV", "TY.15_D9_75C_UV", "TY.7_D0_AC", "TY.8_D6_AC", "TY.9_D9_AC",
                           "TY.16_D0_AC_UV", "TY.17_D6_AC_UV", "TY.18_D9_AC_UV", "TY.19_D0_DNA_UV",
                           "TY.20_D6_DNA_UV", "TY.21_D9_DNA_UV"),
      sampname = fct_relevel(sampname,
                             "Control 1", "Control 2", "D3 at 20C A", "D3 at 20C B", "D3 at 4C A", "D3 at 4C B",
                             "D3 at -20C A", "D3 at -20C B", "D3 at -80C A", "D3 at -80C B", "D8 at 20C A",
                             "D8 at 20C B", "D8 at 4C A", "D8 at 4C B", "D8 at -20C A", "D8 at -20C B",
                             "D8 at -80C A", "D8 at -80C B", "D11 at 20C A", "D11 at 20C B", "D11 at 4C A",
                             "D11 at 4C B", "D11 at -20C A", "D11 at -20C B", "D11 at -80C A", "D11 at -80C B",
                             "Control 3", "D6 dried", "D9 dried", "D0 UV pre-extraction", "D6 UV pre-extraction", "D9 UV pre-extraction",
                             "D0 75C heat", "D6 75C heat", "D9 75C heat", "D0 75C + UV pre-ex", "D6 75C + UV pre-ex", "D9 75C + UV pre-ex", 
                             "D0 autoclave", "D6 autoclave", "D9 autoclave", "D0 autoclave + UV pre-ex", "D6 autoclave + UV pre-ex",
                             "D9 autoclave + UV pre-ex", "D0 UV post-extraction", "D6 UV post-extraction", "D9 UV post-extraction"),
      heat = if_else(experiment == 1, "no heat", heat),
      treatment = if_else(experiment == 1, "none", treatment),
      origin = case_when(
        is.na(experiment) ~ "Other",
        experiment == 1 ~ "Year 1 Samples",
        heat == "autoclave" ~ "Autoclave",
        uv == "UV DNA" ~ "UV DNA",
        TRUE ~ "Year 2 Samples")) %>%
    left_join(qPCR, by = "samplabel") 
  return(allsamps)
}
create_otu_df = function(phy_object, phy_samp_df) {
  physub = phy_object %>% 
    phy.collapse()  # combine OTUs of the same taxonomic classification
  otusub = physub %>% 
    get.otu.melt(sample_data = FALSE) %>%
    tax.plot(data = TRUE) %>%  
    left_join(phy_samp_df) %>%
    mutate(absseqs = pctseqs * seqs_total)
}

# add beta metric, using samp.comparator

# remove additional ranks
phy.tyler <- phy.tyler %>% 
  select(otu,Superkingdom,Phylum,Class,Order,Family,Genus,Species)
# modify sample_data
s.tyler <- phy.tyler %>% get.samp() %>%
  left_join(trace_tbl.tyler,by="sample") %>%
  mutate(temp=ifelse(sample %in% c("1A","1B"),"-80C",temp),
         temp=factor(temp,levels=c("-80C", "-20C", "4C", "room temp")),
         days=as.numeric(str_replace(time,"day ","")),
         time=fct_reorder(time,days),
         time=fct_relabel(time,~str_replace(.x,"day","Day")),
         treatment=factor(treatment,levels=c("none", "75C", "UV", "75C+UV", "autoclave", "autoclave+UV", "UV DNA")),
         sample_number=order(order(treatment,temp,time)),
         letter=str_extract(sample,"(?<=^[0-9]{1,2})[AB]"),
         sample2=paste(treatment,temp,storage,time,sample_number,sep="|"),
         sample2=fct_reordern(sample2,sample_number),
         sample.comparator=case_when(
           experiment==1 ~ "1A",
           experiment==2 ~ "TY.1_D0_NT"
         )) %>%
  group_by(experiment) %>%
  # make a new sample label
  mutate(ord=order(order(experiment,temp,treatment,time,letter)),
         lbl=paste0(experiment,LETTERS[ord]),
         lbl.comparator=lbl[match(sample.comparator,sample)]) %>%
  ungroup() 
sample_data(phy.tyler) <- s.tyler %>% set.samp()

# custom palette
pal <- list("Bacteroidota (phylum)"=Phylum %in% c("Bacteroidetes", "Bacteroidota") ~ shades("#51AB9B", variation = 0.25),
            "Lachnospiraceae (family)"=Family == "Lachnospiraceae" ~ shades("#EC9B96", variation = 0.25),
            "Oscillospiraceae (family)"=Family %in% c("Ruminococcaceae", "Oscillospiraceae") ~ shades("#9AAE73", variation = 0.25),
            "Eubacteriales (order)"=Order %in% c("Clostridiales", "Eubacteriales") ~ shades("#9C854E", variation = 0.25),
            "Actinomycetota (phylum)"=Phylum %in% c("Actinobacteria", "Actinomycetota") ~ shades("#A77097", variation = 0.25),
            # "Enterococcus (genus)"=Genus == "Enterococcus" ~ shades("#129246", variation = 0.15),
            # "Streptococcus (genus)"=Genus == "Streptococcus" ~ shades("#9FB846", variation = 0.15),
            # "Staphylococcus (genus)"=Genus == "Staphylococcus" ~ shades("#f1eb25", variation = 0.15),
            # "Lactobacillus (genus)"=Genus == "Lactobacillus" ~ shades("#3b51a3", variation = 0.15),
            "Pseudomonadota (phylum)"=Phylum %in% c("Proteobacteria", "Pseudomonadota") ~ shades("red", variation = 0.4),
            "Other Bacteria"=TRUE ~ shades("gray", variation = 0.25))


add_dist <- function(phy,method) {
  # phy=phy1;method="horn";sample0="1A"
  varname <- paste0("dist_",method)
  s <- get.samp(phy)
  samp1 <- s$sample.comparator
  samp2 <- s$sample
  comp <- unique(samp1)
  pw <- calc.pairwise(sample1=samp1,sample2=samp2,method=method,phy=phy)
  new.s <- s %>% mutate(!!varname:=pw)
  sample_data(phy) <- new.s %>% set.samp()
  # if (length(comp)!=1) {
  #   cli::cli_abort("YTError: should be one sample comparator!")
  # }
  cli::cli_alert("added var: {.pkg {varname}} (comparator={comp})")
  return(phy)
}

phy1 <- phy.tyler %>% 
  filter(experiment==1) %>%
  mutate_sample_data(baseline=sample==sample.comparator,
                     celsius=ifelse(temp=="room temp",20,as.numeric(str_replace(temp,"C",""))),
                     time.num=days,
                     time.rank=dense_rank(days),
                     temp.num=celsius,
                     temp.rank=dense_rank(temp.num),
                     # qpcr.totalseqs=coalesce(qpcr.totalseqs,100),
                     # for formatting purposes
                     timelabel=ifelse(time=="Day 0","","Storage Time"),
                     templabel="Storage Temperature",
                     time=fct_recode(time,"Day 0 (Baseline)"="Day 0"),
                     letter.rev=fct_rev(letter)) %>%
  add_dist("pct.bray") %>%
  add_dist("horn") %>%
  add_dist("mean.horn") %>%
  add_dist("unfold.horn")

phy2 <- phy.tyler %>%
  filter(experiment==2) %>%
  mutate(baseline=sample==sample.comparator,
         time.num=days,
         time.rank=dense_rank(time.num),
         uv=fct_relevel(uv,"no UV"),
         heat=fct_relevel(heat,"no heat"),
         heat.75C=heat=="75C",
         heat.autoclave=heat=="autoclave",
         uv.dna=uv=="UV DNA",
         # qpcr.totalseqs=coalesce(qpcr.totalseqs,100),
         # for formatting purposes
         timelabel="Storage Time",
         treatmentlabel="Pre-treatment",
         letter=factor(1)) %>%
  add_dist("pct.bray") %>%
  add_dist("horn") %>%
  add_dist("mean.horn") %>%
  add_dist("unfold.horn")

# Samp and OTU dataframes 
allsamps = create_samp_df(phy.tyler)  # get SAMP dataframe
otusub = create_otu_df(phy.tyler, allsamps)  # get OTU dataframe (with mapping/samp data)

firstsamps = allsamps %>%
  dplyr::filter(experiment == 1) %>%
  mutate(replicate = if_else(grepl("A", sample), 1, 2))  # for Figure 1 pairing
secondsamps = allsamps %>%
  dplyr::filter(experiment == 2)

firstotu = otusub %>% 
  dplyr::filter(experiment == 1) %>%
  mutate(replicate = if_else(grepl("A", sample), 1, 2)) 
secondotu = otusub %>%
  dplyr::filter(experiment == 2)

# FIGURE 1A: Overview Stacked Barplot -----------

otu1 <- phy1 %>% 
  # phy.collapse() %>% 
  get.otu.melt()
s1 <- phy1 %>% get.samp(stats=TRUE)

width <- 0.75
g1a <- ggplot() +
  # geom_taxonomy(data=otu1,aes(x=letter.rev, y=pctseqs, fill=otu, label=Species),width=width) +
  geom_taxonomy(data=otu1,aes(x=letter.rev, y=pctseqs, fill=otu),width=width) +
  # expand_limits(x=4) +
  # dashed line over comparator sample
  geom_col(data=filter(s1,baseline),aes(x=letter.rev,y=1),
           width=width,linewidth=0.75,linetype="longdash",color="blue",fill=NA) +
  # sample label text: 1A, ... 1Z
  geom_text(data=s1,aes(x=letter.rev,y=0,label=lbl),hjust=1.2) +
  # make space for sample label text
  expand_limits(y=-0.08) +
  facet_nested(templabel+temp ~ timelabel+time,
               strip=strip_nested( #overarching facets are blank background
                 background_x=list(element_blank(),element_blank(),element_rect(),element_rect(),element_rect(),element_rect()),
                 background_y=list(element_blank(),element_rect(),element_rect(),element_rect(),element_rect())
               )) +
  scale_fill_taxonomy(name="Bacterial Taxa", tax.palette = pal, data=otu1, fill=otu) +
  theme(aspect.ratio=0.6,
        legend.position=c(0.15,0.4),
        axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title = element_blank(), panel.background = element_blank(),
        # creates space with baseline samps
        panel.spacing.x = unit(c(40,5.5,5.5),"points")) +
  coord_flip()
g1a

# with beta diversity on top
littlewidth <- 0.1
g1a.alt <- g1a +
  geom_col(data=s1,aes(x=stage(letter.rev,after_stat=x+(width+littlewidth)/2),y=dist_horn),width=littlewidth,fill="blue") +
  geom_text(data=s1,aes(x=stage(letter.rev,after_stat=x+(width+littlewidth)/2),y=dist_horn,
                        label=str_glue("Horn={sprintf('%.3f',dist_horn)}")),hjust=0,vjust=0.2,size=3)
g1a.alt

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/fig1A_overview.pdf', width = 13, height = 7)
g1a.alt
dev.off()

# FIGURE 1B: PCOAs -----------------------

phy.others = read_rds("/Users/Tyler/Documents/Research/IntegrityProject/other.samps.rds")
others = get.samp(phy.others)

phy.tyler.first <- phy.tyler %>% subset_samples(experiment==1) %>% prune_unused_taxa()

taxa_names(phy.tyler.first) = as.character(refseq(phy.tyler.first)) %>% unname()  # ASK ABOUT THIS
taxa_names(phy.others) = as.character(refseq(phy.others)) %>% unname()

phy.together = merge_phyloseq(otu_table(phy.tyler.first), otu_table(phy.others),
                              tax_table(phy.tyler.first),tax_table(phy.others))
phy.together_test = merge_phyloseq(phy.tyler.first, phy.others)

# MEAN BRAY
bray_dist = calc.distance(phy.together, "bray", mean = TRUE)
# UNFOLD HORN
horn_dist = calc.distance(phy.together_test, "horn", unfold = TRUE)

bray_ord = ordinate(phy.together, distance = bray_dist, method = "NMDS")
bray_data = bray_ord$points %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(allsamps, by = "sample") %>%
  mutate(origin = case_when(
    is.na(experiment) ~ "Other",
    experiment == 1 ~ "Our Study"),
    temp = if_else(origin == "Other", "none", temp)) %>%
  mutate(origin = fct_relevel(origin, "Our Study", "Other"))
first_bray = ggplot(bray_data, aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = temp, shape = origin)) +
  scale_color_manual(
    values = c("-80C" = "steelblue", 
               "-20C" = "cadetblue", 
               "4C" = "cadetblue3", 
               "room temp" = "indianred",
               "none" = "gray3"),
    breaks = c("-80C", "-20C", "4C", "room temp")) +
  # theme(legend.position = "none") +
  labs(title = "Bray (mean)", color = "Storage Temperature", shape = "Origin")
first_bray

# pdf('/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/bray_mean_1.pdf', width = 7, height = 6)
# first_bray
# dev.off()

horn_ord = ordinate(phy.together, distance = horn_dist, method = "NMDS")
horn_data = horn_ord$points %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(allsamps, by = "sample") %>%
  mutate(origin = case_when(
    is.na(experiment) ~ "Other",
    experiment == 1 ~ "Our Study"),
    temp = if_else(origin == "Other", "none", temp)) %>%
  mutate(origin = fct_relevel(origin, "Our Study", "Other"))

first_horn = ggplot(horn_data, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = temp, shape = origin)) +
  scale_color_manual(
    values = c("-80C" = "steelblue", 
               "-20C" = "cadetblue", 
               "4C" = "cadetblue3", 
               "room temp" = "indianred",
               "none" = "gray3"),
    breaks = c("-80C", "-20C", "4C", "room temp")) +
  # theme(legend.position = "none") +
  labs(title = "Horn (unfolded)", color = "Storage Temperature", shape = "Origin")
first_horn

# pdf('/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/horn_unfolded_1.pdf', width = 7, height = 6)
# first_horn
# dev.off()


first_bray_no_legend = first_bray +
  guides(color = "none", shape = "none")
pcoas = ggarrange(first_bray_no_legend, first_horn, nrow = 1, widths = c(3, 4)) 
  # %>% annotate_figure(top = text_grob("PCOAs", size = 16))
pcoas

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/fig1B_pcoas.pdf', width = 10, height = 4)
pcoas
dev.off()

# FIGURE 1C: Compare -------------------

samples.compare <- c("1C","1E","1G","1U", "1W", "1Y")
phy1.compare <- phy1 %>% 
  mutate(compare=lbl %in% samples.compare,
         lbl.compare=str_glue("{lbl} (vs. {lbl.comparator})")) %>%
  filter(compare|baseline)

otu1base <- phy1.compare %>% 
  filter(baseline,prune_unused_taxa=FALSE) %>%
  get.otu.melt(filter.zero=FALSE) %>%
  transmute(otu,pctseqs0=pctseqs)
otu1compare <- phy1.compare %>% 
  filter(compare,prune_unused_taxa=FALSE) %>% 
  get.otu.melt(filter.zero=FALSE) %>%
  left_join(otu1base,by="otu") %>%
  filter(pctseqs>0|pctseqs0>0) %>%
  group_by(sample) %>%
  arrange(desc(pctseqs0),desc(pctseqs)) %>%
  mutate(col=row_number(),
         extra=pctseqs0==0 & pctseqs>0) %>%
  ungroup()
s1compare <- phy1.compare %>% 
  filter(compare,prune_unused_taxa=FALSE) %>% 
  get.samp(stats=TRUE) %>%
  mutate(label=str_glue("InvSimpson={pretty_number(InvSimpson,digits=3)}\nHorn={sprintf('%.3f',dist_horn)}"))

g1.asv <- ggplot() +
  geom_col(data=otu1compare,aes(x=col,y=pctseqs,fill=otu), show.legend = FALSE) +
  geom_step(data=otu1compare,aes(x=col,y=pctseqs0),direction="mid") +
  geom_text(data=s1compare,aes(x=Inf,y=Inf,label=label),hjust=1,vjust=1,color="blue") +
  geom_rect(data=s1compare,aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,linetype=baseline),
            fill=NA,color="blue",show.legend=FALSE) +
  scale_linetype_manual(values=c("TRUE"="longdash","FALSE"=NA)) +
  # geom_bracket(data=filter(otu1compare,extra),
  #              aes(x=col,y=ave(pctseqs,sample,FUN=max),
  #                  fontsize=3,label="unique\nASVs"),tip="square") + 
  scale_fill_taxonomy(name="Bacterial Taxa",data=otu1compare,fill=otu) +
  scale_y_continuous("Relative Abundance",trans=log_epsilon_trans(0.001)) +
  facet_wrap(~lbl.compare,nrow=2) +
  theme(aspect.ratio=1,
        axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title = element_blank(), panel.background = element_blank())

g1.asv

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/fig1C_compare.pdf', width = 10, height = 4)
g1.asv
dev.off()

# FIGURE 1 --------------

bottom = ggarrange(pcoas, g1.asv, widths = c(2, 1))
bottom

figure1 = ggarrange(g1a.alt, bottom, nrow = 2, heights = c(6, 4)) %>% 
  annotate_figure(top = text_grob("Example Figure 1", size = 20))
figure1

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/fig1.pdf', width = 15, height = 13)
figure1
dev.off()

# FIGURE 2A: Tree -----------------------

samples.compare <- c("1A","1B","1I","1Q","1Y")
phy1.tree <- phy1 %>% 
  filter(lbl %in% samples.compare)
gg <- phy.prepare.ggtree(phy1.tree,sortby=lbl,
                         xmin.tip = 0.6,
                         radius.range = c(1.05, 1.4))


tree = gg$ggtree +
  geom_tile(data=gg$otu,aes(x=x,y=y,fill=otu,alpha=pctseqs),color="gray") +
  geom_segment(data=gg$tax,aes(x=x,xend=gg$x.ring.min,y=y,yend=y),color="gray",linetype="dotted") +
  geom_text(data=gg$samp,aes(x=x,y=gg$angle_to_y(90),label=lbl)) +
  # taxa names
  geom_text(data=gg$tax,aes(x=gg$x.ring.max,y=y,label=Species,hjust=hjust,angle=angle),
            color="dark gray",size=2) +
  scale_alpha_continuous(trans=log_epsilon_trans()) +
  scale_fill_taxonomy(data=gg$otu,fill=otu) + 
  guides(alpha = "none") +
  labs(fill = "")
tree

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/fig2a_tree.pdf', width = 13, height = 13)
tree
dev.off()

# FIGURE 2B: LEFSE -----------------------

phy.lefse = phy.tyler %>% 
  phy.collapse(taxranks=c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
lefsesamps = firstsamps %>%
  mutate(timegroup = ifelse(time == "day 0", "early", ifelse(time == "day 3", "early", "late"))) %>%
  mutate(tempgroup = ifelse(temp == "-80C", "Ultra-low temperature storage", "Other storage")) 
sample_data(phy.lefse) = lefsesamps %>% set.samp()

phy.lefse.20 = phy.lefse %>%
  filter(temp %in% c("-80C", "-20C"))
lda_20 = phy.lefse.20 %>%
  lda.effect(class = "temp", subclass = "time")
lda_plot_20 = lda_20 %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "-80C", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Species")  # probably should ask about this
ldatemp_20 = ggplot(lda_plot_20, aes(x = taxon, y = lda, fill = direction)) +
  geom_col() + coord_flip() +
  labs() +
  theme(plot.title = element_text(size = 18, hjust = 0.6),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) +
  scale_fill_manual(values = c("-80C" = "steelblue", "-20C" = "cadetblue"),
                    labels = c("-20C" = "Freezer (-20C)", "-80C" = "Ultra-low (-80C)")) +
  labs(y = "LDA Index", x = "Significantly Different Taxa",
       fill = "Temperature Group",
       title = "-80C vs -20C")
ldatemp_20

phy.lefse.4 = phy.lefse %>%
  filter(temp %in% c("-80C", "4C"))
lda_4 = phy.lefse.4 %>%
  lda.effect(class = "temp", subclass = "time")
lda_plot_4 = lda_4 %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "-80C", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Species")  # probably should ask about this
ldatemp_4 = ggplot(lda_plot_4, aes(x = taxon, y = lda, fill = factor(direction, levels = c("4C", "-80C")))) +
  geom_col() + coord_flip() +
  labs() +
  theme(plot.title = element_text(size = 18, hjust = 0.6),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) +
  scale_fill_manual(values = c("-80C" = "steelblue", "4C" = "cadetblue3"),
                    labels = c("4C" = "Refrigerator (4C)", "-80C" = "Ultra-low (-80C)")) +
  labs(y = "LDA Index", x = "Significantly Different Taxa",
       fill = "Temperature Group",
       title = "-80C vs 4C")
ldatemp_4

phy.lefse.RT = phy.lefse %>%
  filter(temp %in% c("-80C", "room temp"))
lda_RT = phy.lefse.RT %>%
  lda.effect(class = "temp", subclass = "time")
lda_plot_RT = lda_RT %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "-80C", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Species")  # probably should ask about this
ldatemp_RT = ggplot(lda_plot_RT, aes(x = taxon, y = lda, fill = factor(direction, levels = c("room temp", "-80C")))) +
  geom_col() + coord_flip() +
  labs() +
  theme(plot.title = element_text(size = 18, hjust = 0.6),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) +
  scale_fill_manual(values = c("-80C" = "steelblue", "room temp" = "indianred"),
                    labels = c("room temp" = "Room Temperature", "-80C" = "Ultra-low (-80C)")) +
  labs(y = "LDA Index", x = "Significantly Different Taxa",
       fill = "Temperature Group",
       title = "-80C vs RT")
ldatemp_RT

ldatemp_20_comb = ldatemp_20 + guides(fill = "none") 
ldatemp_4_comb = ldatemp_4 + guides(fill = "none") + labs(x = "")
ldatemp_RT_comb = ldatemp_RT + guides(fill = "none") + labs(x = "")
 
ldas = ggarrange(ldatemp_20_comb, ldatemp_4_comb, ldatemp_RT_comb, nrow = 1)
ldas 

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/fig2b_ldas.pdf', width = 15, height = 9)
ldas
dev.off()

# LEFSE Tables --------------------
lda_plot_20 %>% dt
lda_plot_4 %>% dt
lda_plot_RT %>% dt

pretty_freezer_lda = lda_plot_20 %>%
  mutate(test = "Freezer Storage") %>%
  select(test, taxon, direction, lda, kw.pvalue)
pretty_freezer_lda %>% dt

pretty_refrigerator_lda = lda_plot_4 %>%
  mutate(test = "Refrigerator Storage") %>%
  select(test, taxon, direction, lda, kw.pvalue)
pretty_refrigerator_lda %>% dt

pretty_RT_lda = lda_plot_RT %>%
  mutate(test = "Room Temp Storage") %>%
  select(test, taxon, direction, lda, kw.pvalue)
pretty_RT_lda %>% dt

pretty_ldas = rbind(pretty_freezer_lda, pretty_refrigerator_lda, pretty_RT_lda) %>%
  arrange(as.character(taxon))
pretty_ldas %>% dt

write.csv(pretty_freezer_lda, '/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/LDA Tables/freezer_lda.csv')
write.csv(pretty_refrigerator_lda, '/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/LDA Tables/refrigerator_lda.csv')
write.csv(pretty_RT_lda, '/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/LDA Tables/roomtemp_lda.csv')
write.csv(pretty_ldas, '/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/LDA Tables/pretty_ldas.csv')

pretty_ldas = read_csv('/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/LDA Tables/pretty_ldas.csv')
pretty_lda %>% dt

pretty_freezer_lda_genus = lda_20 %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "-80C", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Genus") %>%
  mutate(test = "Freezer Storage") %>%
  select(test, taxon, direction, lda, kw.pvalue)
pretty_freezer_lda_genus %>% dt

pretty_refrigerator_lda_genus = lda_4 %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "-80C", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Genus") %>%
  mutate(test = "Refrigerator Storage") %>%
  select(test, taxon, direction, lda, kw.pvalue)
pretty_refrigerator_lda_genus %>% dt

pretty_roomtemp_lda_genus = lda_RT %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "-80C", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Genus") %>%
  mutate(test = "Room Temperature Storage") %>%
  select(test, taxon, direction, lda, kw.pvalue)
pretty_roomtemp_lda_genus %>% dt

pretty_ldas_genus = rbind(pretty_freezer_lda_genus, pretty_refrigerator_lda_genus, pretty_roomtemp_lda_genus) %>%
  arrange(as.character(taxon))
pretty_ldas_genus %>% dt

write.csv(pretty_freezer_lda_genus, '/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/LDA Tables/genus_freezer_lda.csv')
write.csv(pretty_refrigerator_lda_genus, '/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/LDA Tables/genus_refrigerator_lda.csv')
write.csv(pretty_roomtemp_lda_genus, '/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/LDA Tables/genus_roomtemp_lda.csv')
write.csv(pretty_ldas_genus, '/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/LDA Tables/genus_pretty_ldas.csv')


# FIGURE 2 ---------------------------

figure2 = ggarrange(tree, ldas, nrow = 2, heights = c(4, 3)) %>% 
  annotate_figure(top = text_grob("Figure 2", size = 25))
figure2

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/fig2.pdf', width = 15, height = 18) 
figure2
dev.off()

figure2_hori = ggarrange(tree, ldas, nrow = 1, widths = c(6, 6)) %>% 
  annotate_figure(top = text_grob("Figure 2", size = 25))

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/fig2_hori.pdf', width = 30, height = 14) 
figure2_hori
dev.off()

# FIGURE 3A: Second Overview ---------------

two_otusub = secondotu %>%
  mutate(seqs_total = if_else(treatment == "UV DNA", 3000, seqs_total),  # this is arbitrary
         absseqs = pctseqs * seqs_total)

second_overview = ggplot() +
  geom_taxonomy(data = two_otusub, aes(x = sampname, y = absseqs, fill = Species), # create columns
                show.legend = TRUE, width = 0.8) +
  scale_y_continuous(trans = log_epsilon_trans(10000)) +
  theme(panel.spacing.x = unit(2, "lines")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(), axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0.5, size = 15), axis.title.x = element_text(size = 24), # shows sample name on x-axis
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_text(size = 24),
        plot.title = element_text(hjust = 0.5, size = 35)) +
  geom_text(data = secondsamps, aes(x = sampname, y = 1.05, label = short_number(nseqs)),angle = -90) +
  labs(x = "Sample", y = "Microbial Composition", title = "", fill = "Taxa")
second_overview

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/fig3a_overview.pdf', width = 15, height = 7)
second_overview
dev.off()

# FIGURE 3B: Second PCOAs ----------------

qwer <- function(method,distance,...) {
  ord <- phy.ordinate(phy2,method=method,distance=distance,...)
  opts <- enexprs(...)
  opts_lbl <- paste(names(opts),opts,sep="=",collapse=", ")
  title <- str_glue("{method} {distance} {opts_lbl}")
  ggplot(ord$data,aes(x=axis1,y=axis2)) +
    geom_point(aes(color=treatment),size=1) +
    # geom_text(aes(label=treatment),size=2,vjust=1) +
    theme(aspect.ratio=1) +
    ggtitle(title)
}

second_bray = qwer("PCoA","pct.bray")
second_bray_edits = second_bray + 
  labs(title = "Bray", color = "Treatment") + 
  scale_color_manual(
    values = c("none" = "black", 
               "75C" = "darkcyan", 
               "UV" = "purple", 
               "75C+UV" = "darkolivegreen",
               "autoclave" = "red",
               "autoclave+UV" = "brown",
               "UV DNA" = "violet")) 
second_bray_edits

second_horn = qwer("PCoA","horn")
second_horn_edits = second_horn + 
  labs(title = "Horn", color = "Treatment") + 
  scale_color_manual(
    values = c("none" = "black", 
               "75C" = "darkcyan", 
               "UV" = "purple", 
               "75C+UV" = "darkolivegreen",
               "autoclave" = "red",
               "autoclave+UV" = "brown",
               "UV DNA" = "violet")) 
second_horn_edits

second_bray_no_legend = second_bray_edits +
  guides(color = "none")
second_pcoas = ggarrange(second_bray_no_legend, second_horn_edits, nrow = 1, widths = c(3, 4)) 
# %>% annotate_figure(top = text_grob("PCOAs", size = 16))
second_pcoas

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Spring Graphs/fig3B_second_pcoas.pdf', width = 10, height = 4)
second_pcoas
dev.off()

# FIGURE 3C: Uniques -------------------

compare.uniques <- function(sample1,sample2,phy=phy.tyler) {
  physub <- prune_samples(c(sample1,sample2),phy)
  ranks <- rank_names(physub)
  otusub <- physub %>% get.otu.melt() %>%
    group_by(otu) %>%
    mutate(unique=n()==1) %>%
    ungroup() %>% 
    mutate(otu=fct_reordern(otu,!!!syms(ranks)),
           col=as.numeric(otu),
           sample = factor(sample, levels = unique(c(!!sample1, !!sample2)))) 
  
  ggplot(otusub,aes(x=col,y=pctseqs)) + 
    geom_col(aes(alpha=unique,fill=otu)) +
    geom_point(data=filter(otusub,unique),aes(x=col,y=pctseqs),color="green",size=0.5) +
    scale_fill_taxonomy(data=otusub,fill=otu) +
    # scale_color_manual(values=c("TRUE"="red","FALSE"=NA)) +
    scale_alpha_manual(values=c("TRUE"=1,"FALSE"=0.4)) +
    scale_y_continuous(trans=log_epsilon_trans(0.001)) +
    facet_grid(sample ~ .)
}
compare.uniques.new = function(group1, group2, group1_label, group2_label, phy = phy.tyler) {
  all_targets = c(group1, group2)
  physub = prune_samples(all_targets, phy)
  ranks = rank_names(physub)
  
  otusub = physub %>% 
    get.otu.melt() %>%
    mutate(group_label = case_when(
      sample %in% group1 ~ group1_label,
      sample %in% group2 ~ group2_label)) %>%
    group_by(group_label, otu, across(all_of(ranks))) %>% 
    summarise(pctseqs = mean(pctseqs, na.rm = TRUE), .groups = "drop") %>%
    group_by(otu) %>%
    mutate(unique = n() == 1) %>% 
    ungroup() %>% 
    mutate(otu = fct_reordern(otu, !!!syms(ranks)),
           col = as.numeric(otu),
           group_label = factor(group_label, levels = c(group1_label, group2_label)))
  ggplot(otusub, aes(x = col, y = pctseqs)) + 
    geom_col(aes(alpha = unique, fill = otu)) +
    geom_point(data = filter(otusub, unique), aes(x = col, y = pctseqs), color = "green", size = 0.5) +
    scale_fill_taxonomy(data = otusub, fill = otu) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
    scale_y_continuous(trans = log_epsilon_trans(0.001)) +
    # Change facet to use the new group_label
    facet_grid(group_label ~ .)
}

dry_uniques = compare.uniques("TY.1_D0_NT", "TY.3_D9_Dry") +
  guides(alpha = "none", fill = "none")
heat_uniques = compare.uniques("TY.1_D0_NT", "TY.6_D9_75C") +
  guides(alpha = "none", fill = "none")
uv_uniques = compare.uniques("TY.1_D0_NT", "TY.12_D9_UV") +
  guides(alpha = "none", fill = "none")
heatuv_uniques = compare.uniques("TY.1_D0_NT", "TY.15_D9_75C_UV") +
  guides(alpha = "none", fill = "none")
autoclave_uniques = compare.uniques("TY.1_D0_NT", "TY.9_D9_AC") +
  guides(alpha = "none", fill = "none")
uvdna_uniques = compare.uniques("TY.1_D0_NT", "TY.21_D9_DNA_UV") +
  guides(alpha = "none", fill = "none")

test = ggarrange(dry_uniques, heat_uniques, uv_uniques,
                 heatuv_uniques, autoclave_uniques, uvdna_uniques)
test

# First Stacked Overview ---------------
first_overview = ggplot() +
  geom_taxonomy(data = firstotu, aes(x = sampname, y = absseqs, fill = Species), # create columns
                show.legend = TRUE, width = 0.8) +
  scale_y_continuous(trans = log_epsilon_trans(10000)) +
  theme(panel.spacing.x = unit(2, "lines")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(), axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0.5, size = 15), axis.title.x = element_text(size = 24), # shows sample name on x-axis
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_text(size = 24),
        plot.title = element_text(hjust = 0.5, size = 35)) +
  geom_text(data = firstsamps, aes(x = sampname, y = 1.05, label = short_number(nseqs)),angle = -90) +
  labs(x = "", y = "Microbial Composition", title = "", fill = "Taxa")
first_overview

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/overview_1.pdf', width = 13, height = 7)
first_overview
dev.off()

# First PCAs ---------------
phy.others = read_rds("/Users/Tyler/Documents/Research/IntegrityProject/other.samps.rds")
others = get.samp(phy.others)

phy.tyler.first <- phy.tyler %>% subset_samples(experiment==1) %>% prune_unused_taxa()

taxa_names(phy.tyler.first) = as.character(refseq(phy.tyler.first)) %>% unname()  # ASK ABOUT THIS
taxa_names(phy.others) = as.character(refseq(phy.others)) %>% unname()

phy.together = merge_phyloseq(otu_table(phy.tyler.first), otu_table(phy.others),
                              tax_table(phy.tyler.first),tax_table(phy.others))
phy.together_test = merge_phyloseq(phy.tyler.first, phy.others)

# MEAN BRAY
bray_dist = calc.distance(phy.together, "bray", mean = TRUE)
# UNFOLD HORN
horn_dist = calc.distance(phy.together_test, "horn", unfold = TRUE)

bray_ord = ordinate(phy.together, distance = bray_dist, method = "NMDS")
bray_data = bray_ord$points %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(allsamps, by = "sample") %>%
  mutate(origin = case_when(
    is.na(experiment) ~ "Other",
    experiment == 1 ~ "Our Study"),
    temp = if_else(origin == "Other", "none", temp)) %>%
  mutate(origin = fct_relevel(origin, "Our Study", "Other"))
first_bray = ggplot(bray_data, aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = temp, shape = origin)) +
  scale_color_manual(
    values = c("-80C" = "steelblue", 
               "-20C" = "cadetblue", 
               "4C" = "cadetblue3", 
               "room temp" = "indianred",
               "none" = "gray3"),
    breaks = c("-80C", "-20C", "4C", "room temp")) +
  # theme(legend.position = "none") +
  labs(title = "Bray (mean)", color = "Storage Temperature", shape = "Origin")
first_bray

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/bray_mean_1.pdf', width = 7, height = 6)
first_bray
dev.off()

horn_ord = ordinate(phy.together, distance = horn_dist, method = "NMDS")
horn_data = horn_ord$points %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(allsamps, by = "sample") %>%
  mutate(origin = case_when(
    is.na(experiment) ~ "Other",
    experiment == 1 ~ "Our Study"),
    temp = if_else(origin == "Other", "none", temp)) %>%
  mutate(origin = fct_relevel(origin, "Our Study", "Other"))

first_horn = ggplot(horn_data, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = temp, shape = origin)) +
  scale_color_manual(
    values = c("-80C" = "steelblue", 
               "-20C" = "cadetblue", 
               "4C" = "cadetblue3", 
               "room temp" = "indianred",
               "none" = "gray3"),
    breaks = c("-80C", "-20C", "4C", "room temp")) +
  # theme(legend.position = "none") +
  labs(title = "Horn (unfolded)", color = "Storage Temperature", shape = "Origin")
first_horn

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/horn_unfolded_1.pdf', width = 7, height = 6)
first_horn
dev.off()


first_bray_no_legend = first_bray +
  guides(color = "none", shape = "none")
pcoas = ggarrange(first_bray_no_legend, first_horn, nrow = 1, widths = c(3, 4)) %>% 
  annotate_figure(top = text_grob("PCOAs", size = 16))
pcoas

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/pcoas_1.pdf', width = 10, height = 4)
pcoas
dev.off()

# First Alpha Diversity Bar Graph ---------------

first_divsamps = firstsamps %>%
  mutate(time = fct_relevel(time, "day 0", "day 3", "day 8", "day 11"),
         sampname = fct_relevel(sampname, "Control 1", "Control 2", 
                                "D3 at -80C A", "D3 at -80C B", "D3 at -20C A", "D3 at -20C B", "D3 at 4C A", "D3 at 4C B", "D3 at 20C A", "D3 at 20C B",
                                "D8 at -80C A", "D8 at -80C B", "D8 at -20C A", "D8 at -20C B", "D8 at 4C A", "D8 at 4C B", "D8 at 20C A", "D8 at 20C B",
                                "D11 at -80C A", "D11 at -80C B", "D11 at -20C A", "D11 at -20C B", "D11 at 4C A", "D11 at 4C B", "D11 at 20C A", "D11 at 20C B"),
         temp = case_when(time == "day 0" ~ "control",
                          TRUE ~ temp)) 
simpson_plot = ggplot(first_divsamps, aes(x = sampname, y = InvSimpson, fill = fct_relevel(temp, "-80C", "-20C", "4C", "room temp"))) + 
  geom_col() + 
  facet_grid(~ time, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = c("-80C" = "steelblue", 
               "-20C" = "cadetblue", 
               "4C" = "cadetblue3", 
               "room temp" = "indianred")) +
  theme(legend.position = "right", 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_text(size = 15),
        axis.text.y = element_text(), axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.3, size = 20),
        legend.title = element_text(size = 13), legend.text = element_text(size = 12)) +
  labs(x = "", y = "Inverse Simpson Index", fill = "Storage Temperature") +
  ylim(0, 40)
simpson_plot

shannon_plot = ggplot(first_divsamps, aes(x = sampname, y = Shannon, fill = fct_relevel(temp, "-80C", "-20C", "4C", "room temp"))) + 
  geom_col() + 
  facet_grid(~ time, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = c("-80C" = "steelblue", 
               "-20C" = "cadetblue", 
               "4C" = "cadetblue3", 
               "room temp" = "indianred")) +
  theme(legend.position = "right", 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_text(size = 15),
        axis.text.y = element_text(), axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.3, size = 20),
        legend.title = element_text(size = 13), legend.text = element_text(size = 12)) +
  labs(x = "", y = "Shannon Index", fill = "Storage Temperature") +
  ylim(0, 9)
shannon_plot

alpha = ggarrange(simpson_plot, shannon_plot, nrow = 2, heights = c(1, 1))
alpha

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/alpha_1.pdf', width = 8, height = 6)
alpha
dev.off()














# First LEfSe ---------------------

phy.lefse = phy.tyler %>% 
  phy.collapse(taxranks=c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
lefsesamps = firstsamps %>%
  mutate(timegroup = ifelse(time == "day 0", "early", ifelse(time == "day 3", "early", "late"))) %>%
  mutate(tempgroup = ifelse(temp == "-80C", "Ultra-low temperature storage", "Other storage")) 
sample_data(phy.lefse) = lefsesamps %>% set.samp()

phy.lefse.20 = phy.lefse %>%
  filter(temp %in% c("-80C", "-20C"))
lda_20 = phy.lefse.20 %>%
  lda.effect(class = "temp", subclass = "time")
lda_plot_20 = lda_20 %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "-80C", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Species")  # probably should ask about this
ldatemp_20 = ggplot(lda_plot_20, aes(x = taxon, y = lda, fill = direction)) +
  geom_col() + coord_flip() +
  labs() +
  theme(plot.title = element_text(size = 18, hjust = 0.6),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) +
  scale_fill_manual(values = c("-80C" = "steelblue", "-20C" = "cadetblue"),
                    labels = c("-20C" = "Freezer (-20C)", "-80C" = "Ultra-low (-80C)")) +
  labs(y = "LDA Index", x = "Significantly Different Taxa",
       fill = "Temperature Group",
       title = "-80C vs -20C")
ldatemp_20

phy.lefse.4 = phy.lefse %>%
  filter(temp %in% c("-80C", "4C"))
lda_4 = phy.lefse.4 %>%
  lda.effect(class = "temp", subclass = "time")
lda_plot_4 = lda_4 %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "-80C", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Species")  # probably should ask about this
ldatemp_4 = ggplot(lda_plot_4, aes(x = taxon, y = lda, fill = factor(direction, levels = c("4C", "-80C")))) +
  geom_col() + coord_flip() +
  labs() +
  theme(plot.title = element_text(size = 18, hjust = 0.6),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) +
  scale_fill_manual(values = c("-80C" = "steelblue", "4C" = "cadetblue3"),
                    labels = c("4C" = "Refrigerator (4C)", "-80C" = "Ultra-low (-80C)")) +
  labs(y = "LDA Index", x = "Significantly Different Taxa",
       fill = "Temperature Group",
       title = "-80C vs 4C")
ldatemp_4

phy.lefse.RT = phy.lefse %>%
  filter(temp %in% c("-80C", "room temp"))
lda_RT = phy.lefse.RT %>%
  lda.effect(class = "temp", subclass = "time")
lda_plot_RT = lda_RT %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "-80C", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Species")  # probably should ask about this
ldatemp_RT = ggplot(lda_plot_RT, aes(x = taxon, y = lda, fill = factor(direction, levels = c("room temp", "-80C")))) +
  geom_col() + coord_flip() +
  labs() +
  theme(plot.title = element_text(size = 18, hjust = 0.6),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) +
  scale_fill_manual(values = c("-80C" = "steelblue", "room temp" = "indianred"),
                    labels = c("room temp" = "Room Temperature", "-80C" = "Ultra-low (-80C)")) +
  labs(y = "LDA Index", x = "Significantly Different Taxa",
       fill = "Temperature Group",
       title = "-80C vs RT")
ldatemp_RT

ldas = ggarrange(ldatemp_20, ldatemp_4, ldatemp_RT, nrow = 1)
ldas 

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/ldas_temp_1.pdf", width = 23, height = 12)
ldas
dev.off()

# Compare Uniques ----------------

compare.uniques <- function(sample1,sample2,phy=phy.tyler) {
  physub <- prune_samples(c(sample1,sample2),phy)
  ranks <- rank_names(physub)
  otusub <- physub %>% get.otu.melt() %>%
    group_by(otu) %>%
    mutate(unique=n()==1) %>%
    ungroup() %>% 
    mutate(otu=fct_reordern(otu,!!!syms(ranks)),
           col=as.numeric(otu),
           sample = factor(sample, levels = unique(c(!!sample1, !!sample2)))) 
  
  ggplot(otusub,aes(x=col,y=pctseqs)) + 
    geom_col(aes(alpha=unique,fill=otu)) +
    geom_point(data=filter(otusub,unique),aes(x=col,y=pctseqs),color="green",size=0.5) +
    scale_fill_taxonomy(data=otusub,fill=otu) +
    # scale_color_manual(values=c("TRUE"="red","FALSE"=NA)) +
    scale_alpha_manual(values=c("TRUE"=1,"FALSE"=0.4)) +
    scale_y_continuous(trans=log_epsilon_trans(0.001)) +
    facet_grid(sample ~ .)
}
compare.uniques.new = function(group1, group2, group1_label, group2_label, phy = phy.tyler) {
  all_targets = c(group1, group2)
  physub = prune_samples(all_targets, phy)
  ranks = rank_names(physub)
  
  otusub = physub %>% 
    get.otu.melt() %>%
    mutate(group_label = case_when(
      sample %in% group1 ~ group1_label,
      sample %in% group2 ~ group2_label)) %>%
    group_by(group_label, otu, across(all_of(ranks))) %>% 
    summarise(pctseqs = mean(pctseqs, na.rm = TRUE), .groups = "drop") %>%
    group_by(otu) %>%
    mutate(unique = n() == 1) %>% 
    ungroup() %>% 
    mutate(otu = fct_reordern(otu, !!!syms(ranks)),
           col = as.numeric(otu),
           group_label = factor(group_label, levels = c(group1_label, group2_label)))
  ggplot(otusub, aes(x = col, y = pctseqs)) + 
    geom_col(aes(alpha = unique, fill = otu)) +
    geom_point(data = filter(otusub, unique), aes(x = col, y = pctseqs), color = "green", size = 0.5) +
    scale_fill_taxonomy(data = otusub, fill = otu) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
    scale_y_continuous(trans = log_epsilon_trans(0.001)) +
    # Change facet to use the new group_label
    facet_grid(group_label ~ .)
}

standards = c("1A", "1B")
first80 = compare.uniques.new(standards, c("13A.-80.D11", "13B.-80.D11"), "Standards", "Ultra-low Freezer") +
  guides(alpha = "none", fill = "none")
first20 = compare.uniques.new(standards, c("12A.-20.D11", "12B.-20.D11"), "Standards", "Freezer") +
  guides(alpha = "none", fill = "none")
first4 = compare.uniques.new(standards, c("11A.4C.D11", "11B.4C.D11"), "Standards", "Refrigerator") +
  guides(alpha = "none", fill = "none")
firstRT = compare.uniques.new(standards, c("10A.RT.D11", "10B.RT.D11"), "Standards", "Room Temperature") +
  guides(alpha = "none", fill = "none")
firstuniques = ggarrange(first80, first20, first4, firstRT) %>%
  annotate_figure(top = text_grob("Compare Uniques", size = 20))
firstuniques

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/uniques_1.pdf", width = 15, height = 10)
firstuniques
dev.off()

first80_cherry = compare.uniques("1A","13A.-80.D11") +
  guides(alpha = "none", fill = "none")
first20_cherry = compare.uniques("1A","12A.-20.D11") +
  guides(alpha = "none", fill = "none")
first4_cherry = compare.uniques("1A","11A.4C.D11") +
  guides(alpha = "none", fill = "none")
firstRT_cherry = compare.uniques("1A","10A.RT.D11") +
  guides(alpha = "none", fill = "none")
firstuniques_cherry = ggarrange(first80_cherry, first20_cherry, first4_cherry, firstRT_cherry) %>%
  annotate_figure(top = text_grob("Compare Uniques (cherrypicked replicate)", size = 20))
firstuniques_cherry

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/uniques_cherry_1.pdf", width = 15, height = 10) 
firstuniques_cherry
dev.off()

# Depict ----------------------------

depict_ty <- function(sample1, sample2, phy=phy.tyler) {
  
  subsamps <- c(sample1, sample2)
  physub <- phy %>% 
    filter(sample %in% subsamps)
  
  # 1. Get the base data
  otu_raw <- physub %>% get.otu.melt()
  
  # 2. Create the explicit sorting order based ONLY on sample1
  # We group by OTU and grab the abundance specifically where sample == sample1.
  # If an OTU exists in sample2 but is 0 in sample1, this ensures it gets sorted as 0.
  otu_order <- otu_raw %>%
    group_by(otu) %>%
    summarize(sort_val = sum(pctseqs[sample == sample1], na.rm = TRUE)) %>%
    arrange(sort_val) %>% # Use desc(sort_val) if you want high-to-low
    pull(otu)
  
  # 3. Apply this order to the main dataframe
  otu <- otu_raw %>%
    mutate(first = sample == sample1,
           y = pctseqs,
           # Convert OTU to a factor with the specific levels we calculated above
           otu_ordered = factor(otu, levels = otu_order),
           # Convert to numeric for the x-axis mapping
           x = as.numeric(otu_ordered))
  
  # Prepare the layers
  # We extend the step line slightly to the left/right for aesthetics
  otu1 <- otu %>% 
    filter(first) %>%
    bind_rows(tibble(x = range(.$x, na.rm = TRUE) + c(-1, 1), y = c(0, 0)))
  
  otu2 <- otu %>% filter(!first)
  
  # Plot
  ggplot() +
    geom_col(data = otu2, aes(x = x, y = y, fill = otu)) +
    geom_step(data = otu1, aes(x = x, y = y), direction = "mid") +
    scale_fill_taxonomy(data = otu, fill = otu) +
    scale_y_continuous(trans = log_epsilon_trans(0.001)) +
    theme()
}

depict3("1A","1B")
depict3("1A","11B.4C.D11")

# Depict Full ------------

phy1 <- phy.tyler %>% filter(experiment==1) %>%
  mutate(baseline=sample=="1A",
         templabel="Temp",
         timelabel="Time")


otu1base <- phy1 %>% 
  filter(baseline,prune_unused_taxa=FALSE) %>%
  get.otu.melt(filter.zero=FALSE) %>%
  transmute(otu,pctseqs0=pctseqs)

otu1 <- phy1 %>% 
  get.otu.melt(filter.zero=FALSE) %>%
  left_join(otu1base,by="otu") %>%
  filter(pctseqs>0|pctseqs0>0) %>%
  group_by(sample) %>%
  arrange(desc(pctseqs0),desc(pctseqs)) %>%
  mutate(col=row_number(),
         extra=pctseqs0==0 & pctseqs>0) %>%
  ungroup()
s1 <- phy1 %>% get.samp()

g1.asv <- ggplot() +
  geom_col(data=otu1,aes(x=col,y=pctseqs,fill=otu)) +
  geom_step(data=otu1,aes(x=col,y=pctseqs0),direction="mid") +
  # geom_step(data=otu1,aes(x=col,y=pctseqs0),direction="mid") +
  geom_text(data=s1,aes(x=Inf,y=Inf,label=str_glue("{sample}\n{short_number(nseqs)} seqs")),hjust=1,vjust=1,color="blue") +
  # geom_text(data=s2,aes(x=Inf,y=Inf,label=str_glue("{sample}\n{short_number(nseqs)} seqs")),hjust=1,vjust=1,color="blue") +
  geom_rect(data=filter(otu1,baseline),
            aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf),
            fill=NA,color="blue",linetype="longdash") +
  geom_bracket(data=filter(otu1,extra),
               aes(x=col,y=ave(pctseqs,sample,FUN=max),
                   fontsize=3,label="unique\nASVs"),tip="square") + 
  scale_fill_taxonomy(data=otu1,fill=otu) +
  scale_y_continuous(trans=log_epsilon_trans(0.001)) +
  facet_nested(templabel+temp+letter~timelabel+time) + 
  labs(title = "Depict Overview")
g1.asv

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/depict_overview_1.pdf", width = 13, height = 18) 
g1.asv
dev.off()





