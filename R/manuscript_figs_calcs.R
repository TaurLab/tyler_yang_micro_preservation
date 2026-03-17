
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


# generate phy1 and phy2 -----------------------------------------------------------


phy1 <- phy.tyler %>% 
  filter(experiment==1) %>%
  mutate_sample_data(baseline=sample==sample.comparator,
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
         timelabel="Storage Time",
         treatmentlabel="Pre-treatment",
         letter=factor(1)) %>%
  add_dist("pct.bray") %>%
  add_dist("horn") %>%
  add_dist("mean.horn") %>%
  add_dist("unfold.horn")


# fig 1A: experiment 1 stackplot -------------------------------------------------------------------

otu1 <- phy1 %>% 
  phy.collapse() %>% 
  get.otu.melt()
s1 <- phy1 %>% get.samp(stats=TRUE)

width <- 0.75
g1a <- ggplot() +
  geom_taxonomy(data=otu1,aes(x=letter.rev, y=pctseqs, fill=otu),width=width) +
  # expand_limits(x=4) +
  # dashed line over comparator sample
  geom_col(data=filter(s1,baseline),aes(x=letter.rev,y=1),
           width=width,linewidth=0.75,linetype="longdash",color="blue",fill=NA) +
  # sample label text: 1A, ... 1Z
  geom_text(data=s1,aes(x=letter.rev,y=0,label=lbl),hjust=1) +
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
                        label=str_glue("Horn={sprintf('%.3f',dist_horn)}")),hjust=0,size=3)
g1a.alt


# step compare ------------------------------------------------------------

samples.compare <- c("1A","1C","1L","1Z")
phy1.compare <- phy1 %>% 
  mutate(compare=lbl %in% samples.compare) %>%
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
  mutate(label=str_glue("InvSimpson={pretty_number(InvSimpson)}\nHorn={sprintf('%.3f',dist_horn)}"))


g1.asv <- ggplot() +
  geom_col(data=otu1compare,aes(x=col,y=pctseqs,fill=otu)) +
  geom_step(data=otu1compare,aes(x=col,y=pctseqs0),direction="mid") +
  geom_text(data=s1compare,aes(x=Inf,y=Inf,label=label),hjust=1,vjust=1,color="blue") +
  # geom_text(data=s2,aes(x=Inf,y=Inf,label=str_glue("{sample}\n{short_number(nseqs)} seqs")),hjust=1,vjust=1,color="blue") +
  geom_rect(data=s1compare,aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,linetype=baseline),
            fill=NA,color="blue") +
  scale_linetype_manual(values=c("TRUE"="longdash","FALSE"=NA)) +
  geom_bracket(data=filter(otu1compare,extra),
               aes(x=col,y=ave(pctseqs,sample,FUN=max),
                   fontsize=3,label="unique\nASVs"),tip="square") + 
  scale_fill_taxonomy(data=otu1compare,fill=otu) +
  scale_y_continuous(trans=log_epsilon_trans(0.001)) +
  facet_wrap(~lbl,nrow=1)

g1.asv

g1a.alt

# fig 2A: experiment 2 stackplot ------------------------------------------


otu2 <- phy2 %>% 
  # for speed
  phy.collapse() %>%
  get.otu.melt() %>%
  group_by(treatment,time) %>%
  mutate(sample.number=as.numeric(factor(sample))) %>%
  ungroup() 
s2 <- phy2 %>% get.samp()

width <- 0.75
g2a <- ggplot() +
  geom_taxonomy(data=otu2,aes(x=letter,y=pctseqs,fill=otu),width=width) +
  geom_col(data=filter(s2,baseline),aes(x=letter,y=1),
           width=width,linewidth=0.75,linetype="longdash",color="blue",fill=NA) +
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

# with beta diversity on top
littlewidth <- 0.1
g2a.alt <- g2a +
  geom_col(data=s2,aes(x=stage(letter,after_stat=x+(width+littlewidth)/2),y=dist_horn),width=littlewidth,fill="steelblue") +
  geom_text(data=s2,aes(x=stage(letter,after_stat=x+(width+littlewidth)/2),y=dist_horn,
                        label=str_glue("Horn={sprintf('%.3f',dist_horn)}")),hjust="inward",size=3)
g2a.alt



