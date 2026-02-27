
# set up data -------------------------------------------------------------

library(yingtools2)
library(tidyverse)
library(phyloseq)
library(ggtree)
library(ggh4x)
rm(list=ls())

load("data/phy.tyler.RData")
load("data/phy.tyler.additional.RData")
phy.others = read_rds("data/other.samps.rds")


phy.tyler <- phy.tyler %>% 
  select(otu,Superkingdom,Phylum,Class,Order,Family,Genus,Species)
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
         lbl=paste0(experiment,LETTERS[ord])) %>%
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
  if (length(comp)!=1) {
    cli::cli_abort("YTError: should be one sample comparator!")
  }
  cli::cli_alert("added var: {.pkg {varname}} (comparator={comp})")
  return(phy)
}

phy1 <- phy.tyler %>% 
  filter(experiment==1) %>%
  phy.collapse.bins() %>% 
  mutate_sample_data(baseline=sample==sample.comparator,
                     #baseline=sample=="1A",
                     # for formatting purposes
                     timelabel=ifelse(time=="Day 0","","Storage Time"),
                     templabel="Storage Temperature",
                     time=fct_recode(time,"Day 0 (Baseline)"="Day 0"),
                     letter=fct_rev(letter),
                     one=factor(1)) %>%
  add_dist("pct.bray") %>%
  add_dist("horn") %>%
  add_dist("mean.horn") %>%
  add_dist("unfold.horn")
otu1 <- phy1 %>% get.otu.melt()
s1 <- phy1 %>% get.samp(stats=TRUE)



# fig 1A: experiment 1 stackplot -------------------------------------------------------------------



g1a <- ggplot() +
  geom_taxonomy(data=otu1,aes(x=letter, y=pctseqs, fill=otu, label=Species)) +
  # dashed line over comparator sample
  geom_col(data=filter(s1,baseline),aes(x=letter,y=1),
           width=0.95,linewidth=0.75,linetype="longdash",color="blue",fill=NA) +
  # sample label text: 1A, ... 1Z
  geom_text(data=s1,aes(x=letter,y=0,label=lbl),hjust=1) +
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





# fig 2A: experiment 2 stackplot ------------------------------------------


phy2 <- phy.tyler %>%
  filter(experiment==2) %>%
  phy.collapse.bins() %>% 
  mutate(baseline=sample=="TY.1_D0_NT",
         timelabel="Storage Time",
         treatmentlabel="Pre-treatment",
         letter=factor(1))

otu2 <- phy2 %>% get.otu.melt() %>%
  group_by(treatment,time) %>%
  mutate(sample.number=as.numeric(factor(sample))) %>%
  ungroup() 
s2 <- phy2 %>% get.samp()

g2a <- ggplot() + 
  geom_taxonomy(data=otu2,aes(x=letter,y=pctseqs,fill=otu,label=Species)) +
  geom_col(data=filter(s2,baseline),aes(x=letter,y=1),
           width=0.95,linewidth=0.75,linetype="longdash",color="blue",fill=NA) +
  geom_text(data=s2,aes(x=letter,y=0,label=lbl),hjust=1) +
  expand_limits(y=-0.08) +
  facet_nested(treatmentlabel+treatment~timelabel+time,
               strip=strip_nested( #overarching facets are blank background
                 background_x=list(element_blank(),element_rect(),element_rect(),element_rect()),
                 background_y=list(element_blank(),element_rect(),element_rect(),element_rect(),element_rect(),element_rect(),element_rect(),element_rect())
               )) +
  scale_fill_taxonomy(data=otu2,fill=otu,
                      # guide=guide_taxonomy(ncol=4,downwards=TRUE),
                      tax.palette = pal) +
  theme(aspect.ratio=0.3,
        legend.position="bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank()) +
  labs(fill="Bacterial Taxa") +
  coord_flip()
g2a







