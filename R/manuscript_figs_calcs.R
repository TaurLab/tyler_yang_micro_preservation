
# set up data -------------------------------------------------------------

library(yingtools2)
library(tidyverse)
library(phyloseq)
library(ggtree)
library(ggh4x)
library(vegan)
library(pairwiseAdonis)
library(ggpubr)
library(ggrepel)
library(decontam)
rm(list=ls())

load("data/phy.tyler.RData")
load("data/phy.tyler.additional.RData")
phy.others = read_rds("data/other.samps.rds")

# correct negative phylo tree lengths
tr <- phy_tree(phy.tyler)
tr$edge.length[tr$edge.length<0] <- 0
phy_tree(phy.tyler) <- tr
rm(tr)

# add beta metric, using samp.comparator

# remove additional ranks
phy.tyler <- phy.tyler %>% 
  select(otu,Superkingdom,Phylum,Class,Order,Family,Genus,Species)



# modify sample_data
s.tyler <- phy.tyler %>% get.samp() %>%
  left_join(trace_tbl.tyler,by="sample") %>%
  mutate(temp=ifelse(sample %in% c("1A","1B"),"n/a",temp),
         temp=factor(temp,levels=c("n/a","-80C", "-20C", "4C", "room temp")),
         days=as.numeric(str_replace(time,"day ","")),
         time=fct_reorder(time,days),
         time=fct_relabel(time,~str_replace(.x,"day","Day")),
         treatment=factor(treatment,levels=c("none", "75C", "UV", "75C+UV", "autoclave", "autoclave+UV", "UV DNA")),
         temp.lbl=fct_recode(temp,"italic('(n/a)')"="n/a","-80*degree*C"="-80C","-20*degree*C"="-20C",
                             "4*degree*C"="4C","'room temp'"="room temp"),
         time.lbl=fct_relabel(time,function(x) paste0("'",x,"'")),
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
         lbl.comparator=lbl[match(sample.comparator,sample)],
         temp.abbrev=case_match(temp,"room temp"~"RT",.default=temp),
         treatment.abbrev=str_replace(treatment,"autoclave","auto"),
         time.abbrev=paste0("d",days),
         lbl2=str_glue("{lbl} ({treatment.abbrev}|{time.abbrev},{temp.abbrev})"),
         lbl3=ifelse(experiment==1,
                     str_glue("{lbl} ({time.abbrev},{temp.abbrev})"),
                     str_glue("{lbl} ({treatment.abbrev}|{time.abbrev})"))) %>%
  ungroup() 
sample_data(phy.tyler) <- s.tyler %>% set.samp()

# custom palette
pal <- list("Bacteroidota (phylum)"=Phylum %in% c("Bacteroidetes", "Bacteroidota") ~ shades("#51AB9B", variation = 0.25),
            "Lachnospiraceae (family)"=Family == "Lachnospiraceae" ~ shades("#EC9B96", variation = 0.25),
            "Oscillospiraceae (family)"=Family %in% c("Ruminococcaceae", "Oscillospiraceae") ~ shades("#9AAE73", variation = 0.25),
            "Other Eubacteriales (order)"=Order %in% c("Clostridiales", "Eubacteriales") ~ shades("#9C854E", variation = 0.25),
            "Actinomycetota (phylum)"=Phylum %in% c("Actinobacteria", "Actinomycetota") ~ shades("#A77097", variation = 0.25),
            # "Enterococcus (genus)"=Genus == "Enterococcus" ~ shades("#129246", variation = 0.15),
            # "Streptococcus (genus)"=Genus == "Streptococcus" ~ shades("#9FB846", variation = 0.15),
            # "Staphylococcus (genus)"=Genus == "Staphylococcus" ~ shades("#f1eb25", variation = 0.15),
            # "Lactobacillus (genus)"=Genus == "Lactobacillus" ~ shades("#3b51a3", variation = 0.15),
            "Pseudomonadota (phylum)"=Phylum %in% c("Proteobacteria", "Pseudomonadota") ~ shades("red", variation = 0.4),
            "Other Bacteria"=TRUE ~ shades("gray", variation = 0.25))

add_dist <- function(phy,method,comparator=NULL,varname=NULL) {
  # phy=phy1;method="horn";sample0="1A"
  varname <- enquo(varname)
  if (rlang::quo_is_null(varname)) {
    varname <- paste0("dist_",method)  
  } else {
    varname <- as_label(varname)
  }
  s <- get.samp(phy)
  samp2 <- s$sample
  if (is.null(comparator)) {
    samp1 <- s$sample.comparator  
  } else {
    samp1 <- rep_along(samp2,comparator)
  }
  pw <- calc.pairwise(sample1=samp1,sample2=samp2,method=method,phy=phy)
  new.s <- s %>% mutate(!!varname:=pw)
  sample_data(phy) <- new.s %>% set.samp()
  cli::cli_alert("added var: {.pkg {varname}} (comparator={unique(samp1)})")
  return(phy)
}

# generate phy1 and phy2 -----------------------------------------------------------

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
                     # time=fct_recode(time,"Day 0 (Baseline)"="Day 0"),
                     letter.rev=fct_rev(letter),
                     lbl=fct_reordern(lbl,lbl),
                     lbl2=fct_reordern(lbl2,lbl),
                     lbl3=fct_reordern(lbl3,lbl)) %>%
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
         uv.sample=uv=="UV",
         # qpcr.totalseqs=coalesce(qpcr.totalseqs,100),
         # for formatting purposes
         timelabel="Storage Time",
         treatmentlabel="Treatment",
         letter=factor(1))

phy.pcrneg.control <- read_rds("data/phy.pcrneg.control.rds") %>%
  select_tax_table(otu,Superkingdom,Phylum,Class,Order,Family,Genus,Species) %>%
  mutate_sample_data(sample=str_replace(sample,"..pool1161",""),
                     lbl="PCRNeg",lbl2="PCRNeg",lbl3="PCRNeg",
                     treatment="PCRNeg")
phy2a <- phy.combine(phy2,phy.pcrneg.control) %>%
  mutate(lbl=fct_reordern(lbl,lbl),
         lbl2=fct_reordern(lbl2,lbl),
         lbl3=fct_reordern(lbl3,lbl),
         is.neg.control=sample=="PCRNeg4") %>%
  add_dist("horn",comparator="TY.1_D0_NT",varname=dist_2A) %>%
  add_dist("mean.horn",comparator="TY.1_D0_NT",varname=dist_2A_meanhorn) %>%
  add_dist("horn",comparator="PCRNeg4",varname=dist_PCRNeg)

# extract taxids from tax.blast.tyler and add to phy2a
tax.taxid <- tax.blast.tyler %>%
  group_by(otu) %>% arrange(evalue.rank,staxid) %>% slice(1) %>% ungroup() %>%
  transmute(otu,taxid=as.character(staxid))
tax.2a <- phy2a %>% get.tax() %>% left_join(tax.taxid,by="otu")
tax_table(phy2a) <- tax.2a %>% set.tax()


# fig 1a and 1b: experiment 1 stackplot ---------------------------------------------------------------

phy1a <- phy1 %>%
  mutate(x=as.numeric(lbl),
         baseline=lbl=="1A",
         dist=dist_mean.horn,
         top_label="Storage Temperature / Time")
otu1 <- phy1a %>% get.otu.melt()
s1a <- phy1a %>% get.samp(stats=TRUE)



panel.spacing.x <- s1a %>% distinct(temp,time) %>%
  arrange(temp,time) %>% 
  group_by(temp) %>%
  transmute(x=ifelse(row_number()==n(),3.5,1.5)) %>%
  ungroup() %>% slice(-n()) %>% pull(x) %>% unit("points")
facet <- facet_nested(. ~ temp.lbl+time.lbl,space="free_x",scales="free_x",
                      nest_line=element_line(),
                      resect=unit(3,"pt"),
                      solo_line=TRUE,
                      labeller=label_parsed)
dtext <- s1a %>% filter(sample=="1A") %>%
  mutate(label="baseline\nsample",
         x=1.5,y=0.02,
         xend=1,yend=0.003)
dtext.layer <- list(geom_text(data=dtext,aes(x=x,y=y,label=label),vjust=-0.02,lineheight=0.8),
                    geom_segment(data=dtext,aes(x=x,y=y,xend=xend,yend=yend)))
width <- 0.95

g1.dist <- ggplot(s1a) +
  geom_col(aes(x=lbl,y=dist),fill="steelblue",width=width) +
  geom_text(aes(x=lbl,y=dist,label=pretty_number(dist),color=baseline),
            vjust=0,size=3,show.legend = FALSE) +
  dtext.layer +
  scale_y_continuous(name="Distance\n(compared with baseline)")  +
  scale_color_manual(values=c("TRUE"="red","FALSE"="black")) +
  xlab("Sample") +
  facet +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        panel.spacing.x = panel.spacing.x)


# annotate("text",label="asdf",x=1,y=0.05,group="Day 0")

g1.tax <- ggplot(otu1) +
  geom_taxonomy(aes(x=lbl,y=pctseqs,fill=otu,label=Species),width=width,label.split=TRUE) +
  geom_col(data=filter(s1a,lbl=="1A"),aes(x=lbl,y=1,color="baseline sample (1A)"),
           linewidth=1,linetype="longdash",fill=NA,width=width) +
  scale_y_continuous("Relative abundance",
                     expand=FALSE,labels=scales::label_percent()) +
  scale_color_manual(name="Sample Comparison",values=c("baseline sample (1A)"="red")) +
  scale_fill_taxonomy(name="Bacterial Taxa",data=otu1,tax.palette=pal,fill=otu) +
  xlab("Sample") +
  facet +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        legend.key.size = unit(0.85,"lines"),
        panel.spacing.x = panel.spacing.x)

g1.tax
g1.dist

# g.fig1a <- gg.stack2(g1.dist,g1.tax,heights=c(2,3))
# 
# g1.tax
# g1.dist
# g.fig1a
# 
# pdf("plots/fig1a - exp1 stackplot.pdf",width=12,height=8)
# g.fig1a
# dev.off()
# shell.exec("plots/fig1a - exp1 stackplot.pdf")
# 
# g.fig1a %>% copy.to.clipboard.gg(width=12,height=8)


# get.dist <- function(method,phy=phy1a) {
#   phy <- phy %>% add_dist(method,varname=dist)
#   s <- get.samp(phy)
#   ggplot(s) +
#     geom_col(aes(x=lbl3,y=dist),fill="steelblue",width=width) +
#     geom_text(aes(x=lbl3,y=dist,label=pretty_number(dist),color=baseline),
#               vjust=0,size=3,show.legend = FALSE) +
#     scale_y_continuous(name=str_glue("{str_to_title(method)} Distance\n(compared with baseline)"))  +
#     scale_color_manual(values=c("TRUE"="red","FALSE"="black")) +
#     facet +
#     theme(panel.spacing.x = panel.spacing.x)
# }
# g1.horn <- get.dist("horn")
# g1.mean.horn <- get.dist("mean.horn")
# g1.unfold.horn <- get.dist("unfold.horn")
# g1.pct.bray <- get.dist("pct.bray")
# gg.stack2(g1.dist,
#           g1.horn,
#           g1.unfold.horn,
#           g1.mean.horn,
#           g1.pct.bray,
#           g1.tax,heights=c(2,2,2,2,2,3))

# g1.invsimpson <- ggplot(s1a) +
#   geom_col(aes(x=lbl3,y=InvSimpson),fill="steelblue",width=width) +
#   geom_text(aes(x=lbl3,y=InvSimpson,label=pretty_number(InvSimpson,digits=3),color=baseline),
#             vjust=0,size=3,show.legend = FALSE) +
#   scale_y_continuous(name="Inverse Simpson")  +
#   scale_color_manual(values=c("TRUE"="red","FALSE"="black")) +
#   facet +
#   theme(panel.spacing.x = panel.spacing.x)
# 
# g1.qpcr <- ggplot(s1a) +
#   geom_col(aes(x=lbl3,y=qpcr.totalseqs),fill="steelblue",width=width) +
#   geom_text(aes(x=lbl3,y=qpcr.totalseqs,label=short_number(qpcr.totalseqs),color=baseline),
#             vjust=0,size=3,show.legend = FALSE) +
#   scale_y_continuous(name="Total 16S qPCR",
#                      trans=log_epsilon_trans(10000),labels=pretty_power10)  +
#   scale_color_manual(values=c("TRUE"="red","FALSE"="black")) +
#   facet +
#   theme(panel.spacing.x = panel.spacing.x)
# g1.invsimpson
# g1.qpcr




# fig 1b: exp 1, pcoa and permanova ---------------------------------------------------------------

phy1b <- phy1 %>%
  mutate(storage=days!=0,
         time=ordered(time),
         temp=ordered(temp))
dist1 <- calc.distance(phy1b,"horn")
s1b <- get.samp(phy1b)

# valid
permanova <- adonis2(dist1 ~ temp*time, data=s1b)
permanova

# beta dispersion not signif, therefore differences are not from spread, can proceed
bd <- betadisper(dist1, s1b$temp)
permutest(bd)
# pairwise contrasts of time
# pairwise0 <- pairwise.adonis(dist1,factors=s1$time, p.adjust.m="BH")
pairwise <- pairwise.adonis2(dist1 ~ time, data=s1b) %>%
  keep(is.data.frame) %>%
  imap(~{
    .x %>% rownames_to_column("element") %>% 
      rename(p.value=`Pr(>F)`) %>% 
      mutate(pair=.y) %>%
      filter(!is.na(F)) %>% as_tibble()
  }) %>% list_rbind()

# generate pcoa plot with permanova results
ord1 <- phy.ordinate(phy1b,method="PCoA",distance=dist1)
g.fig1b.pcoa <- ggplot(arrange(ord1$data,lbl),aes(x=axis1,y=axis2,shape=temp)) +
  geom_point(aes(fill=time,color=time),size=4) + 
  # geom_text_repel(size=3,show.legend = FALSE) +
  # geom_text_repel(aes(label=lbl),size=3,vjust=1.4) +
  geom_text(aes(label=lbl),size=3,vjust=1.4) +
  # scale_fill_brewer(type="seq",palette=3) +
  scale_shape_manual(values=c("n/a"=21, "-80C"=22, "-20C"=23, "4C"=24, "room temp"=25)) +
  xlab("PC1") + ylab("PC2") +
  theme(aspect.ratio=1)
ptbl <- permanova %>% as.data.frame() %>% rownames_to_column("var") %>%
  filter(!is.na(`Pr(>F)`)) %>%
  mutate(var=str_replace_all(var,":"," %*% ")) %>%
  select("Predictor"=var,"R^2"=R2,"italic(P)*'-value'"=`Pr(>F)`) %>% 
  mutate(across(.cols=where(is.numeric),.fns=pretty_number))

target.pairs <- c("Day 0_vs_Day 3" = "Day 0 vs Day 3", "Day 3_vs_Day 8" = "Day 3 vs Day 8", "Day 11_vs_Day 8" = "Day 8 vs Day 11")
ptbl2 <- pairwise %>%
  filter(pair %in% names(target.pairs)) %>%
  mutate(pair=recode2(pair,target.pairs,as.factor=TRUE),
         adjusted.pval=p.adjust(p.value,"BH")) %>%
  arrange(pair) %>%
  select("Pair"=pair,"R^2"=R2,"italic(P)*'-value (adj)'"=adjusted.pval) %>%
  mutate(across(.cols=where(is.numeric),.fns=pretty_number))

tt <- ttheme(base_style="light",base_size=9,
             tbody.style = tbody_style(size=9,fill=NA,parse=TRUE),
             colnames.style = colnames_style(size=9,fill=NA,parse=TRUE))
gtbl1 <- ggtexttable(ptbl,rows=NULL,theme=tt)
gtbl2 <- ggtexttable(ptbl2,rows=NULL,theme=tt)


pos1 <- c(0.75,0.35)
pos2 <- c(0.75,0.15)
g.fig1b <- g.fig1b.pcoa + 
  patchwork::inset_element(gtbl1,pos1[1],pos1[2],pos1[1],pos1[2]) +
  patchwork::inset_element(gtbl2,pos2[1],pos2[2],pos2[1],pos2[2])
g.fig1b

pdf("plots/fig1b - exp1 pcoa.pdf",width=8,height=8)
g.fig1b
dev.off()
shell.exec("plots/fig1b - exp1 pcoa.pdf")

# fig 1c: step compare exp1 ------------------------------------------------------------

samples.compare <- c("1B","1E","1Q","1Z")
phy1.compare <- phy1 %>% 
  mutate(compare=lbl %in% samples.compare,
         lbl.compare=str_glue("{lbl3} (vs. {lbl.comparator})")) %>%
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
  # mutate(label=str_glue("InvSimpson={pretty_number(InvSimpson,digits=3)}\nHorn={sprintf('%.3f',dist_horn)}"))
  mutate(label=str_glue("InvSimpson={pretty_number(InvSimpson,digits=3)}\nmean Horn={sprintf('%.3f',dist_mean.horn)}"))


g1.asv <- ggplot() +
  geom_col(data=otu1compare,aes(x=col,y=pctseqs,fill=otu)) +
  geom_step(data=otu1compare,aes(x=col,y=pctseqs0),direction="mid") +
  geom_text(data=s1compare,aes(x=Inf,y=Inf,label=label),hjust=1,vjust=1,color="blue") +
  geom_rect(data=s1compare,aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,linetype=baseline),
            fill=NA,color="blue",show.legend=FALSE) +
  scale_linetype_manual(values=c("TRUE"="longdash","FALSE"=NA)) +
  geom_bracket(data=filter(otu1compare,extra),
               aes(x=col,y=ave(pctseqs,sample,FUN=max),
                   fontsize=3,label="unique\nASVs"),tip="square") + 
  scale_fill_taxonomy(name="Bacterial Taxa",data=otu1compare,fill=otu,tax.palette=pal) +
  scale_y_continuous("Relative Abundance",trans=log_epsilon_trans(0.001)) +
  facet_wrap(~lbl.compare,nrow=1) +
  theme(aspect.ratio=1,
        legend.key.size = unit(0.85,"lines"),
        axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title = element_blank(), panel.background = element_blank())

g1.asv

pdf("plots/fig1c - exp1 step compare.pdf",width=10,height=5)
g1.asv
dev.off()

shell.exec("plots/fig1c - exp1 step compare.pdf")



# fig 2A: stackplot ---------------------------------------------------------


s2a <- get.samp(phy2a)
s2a.long <- s2a %>% pivot_longer(cols=c(dist_2A,dist_PCRNeg)) %>%
  mutate(name=case_match(name,"dist_2A"~"2A","dist_PCRNeg"~"PCRNeg"))
#calculate asv in/out data
base <- phy2a %>% filter(sample=="TY.1_D0_NT") %>% get.otu.melt(sample_data=FALSE) %>% select(otu,taxid)
ctrl <- phy2a %>% filter(lbl=="PCRNeg") %>% get.otu.melt(sample_data=FALSE) %>% select(otu,taxid)
otu2a <- get.otu.melt(phy2a) %>%
  mutate(taxid.in.base=taxid %in% base$taxid,
         taxid.in.ctrl=taxid %in% ctrl$taxid,
         taxid.status=case_when(
           taxid.in.base & taxid.in.ctrl ~ "both",
           taxid.in.base & !taxid.in.ctrl ~ "in 2A",
           !taxid.in.base & taxid.in.ctrl ~ "in PCRNeg",
           TRUE ~ "neither"),
         Species=str_replace_all(Species,"\\[|\\]",""))
samp2a <- otu2a %>%
  group_by(!!!syms(sample_variables(phy2a))) %>%
  summarize(pct.taxid.base=mean(taxid.in.base),
            pct.taxid.ctrl=mean(taxid.in.ctrl),
            .groups="drop") %>%
  mutate(x=as.numeric(lbl3),
         xmin=x-0.475,xmax=x+0.475,
         qpcr.lbl=case_when(
           lbl3=="PCRNeg" ~ "(n/a)",
           is.na(qpcr.totalseqs) ~ "(undetectable)",
           TRUE ~ NA_character_))
samp2a.long <- samp2a %>% pivot_longer(cols=c(pct.taxid.base,pct.taxid.ctrl)) %>%
  mutate(name=case_match(name,"pct.taxid.base"~"2A","pct.taxid.ctrl"~"PCRNeg"))
treat.groups <- samp2a.long %>%
  group_by(treatment) %>%
  summarize(xmin=min(xmin),xmax=max(xmax),
            .groups="drop") %>%
  mutate(b.xmin=xmin+0.05,b.xmax=xmax-0.05,
         b.xmin2=b.xmin+0.05,b.xmax2=b.xmax-0.05)

panel.spacing.x <- s2a %>% distinct(treatment,time) %>%
  arrange(treatment,time) %>% 
  group_by(treatment) %>%
  transmute(x=ifelse(row_number()==n(),3.5,1.5)) %>%
  ungroup() %>% slice(-n()) %>% pull(x) %>% unit("points")

facet <- facet_nested(. ~ treatment+time,scales="free_x",space="free_x",
                      nest_line=element_line(),
                      resect=unit(3,"pt"),
                      solo_line=TRUE)

width <- 0.95
g.qpcr <- ggplot(samp2a) + 
  geom_col(aes(x=lbl3,y=qpcr.totalseqs),width=width) + 
  geom_text(aes(x=lbl3,y=0,label=qpcr.lbl),hjust=0,angle=90) +
  scale_y_continuous("Total abundance by 16S qPCR",trans=log_epsilon_trans(1000),
                     labels=pretty_power10) +
  facet +
  theme(panel.spacing.x=panel.spacing.x)

g.asv.compare <- ggplot(samp2a.long) +
  geom_col(aes(x=lbl3,y=value,fill=name),position="dodge") +
  facet +
  scale_fill_manual("Compared sample",values=c("PCRNeg"="blue","2A"="red")) +
  theme(panel.spacing.x=panel.spacing.x)
g.asv.compare

g.tax.compare <- ggplot(otu2a) +
  geom_taxonomy(aes(x=lbl3,y=pctseqs,label=Species,fill=otu),tax.palette=pal,label.split=TRUE,fontsize=3) +
  # dashed boxes
  geom_col(data=filter(s2a,lbl %in% c("2A","PCRNeg")),
            aes(x=lbl3,y=1,color=lbl),fill=NA,width=width,
           linetype="longdash",linewidth=0.75,show.legend = FALSE) +
  scale_color_manual("",values=c("PCRNeg"="blue","2A"="red")) +
  scale_x_discrete("Sample") +
  scale_y_continuous("Bacterial Composition\nRelative abundance",expand=FALSE) +
  scale_fill_taxonomy(name="Bacterial taxa",data=otu2a,fill=otu,tax.palette=pal) +
  facet + 
  theme(legend.key.size = unit(0.85,"lines"),
        axis.text.x=element_text(hjust=0,vjust=0.5,angle=90),
        axis.text.y=element_blank(),
        axis.ticks=element_blank())

g.fig2a <- gg.stack2(g.qpcr,g.asv.compare,g.tax.compare)
g.fig2a



# g.dist <- ggplot(s2a.long,aes(x=lbl3,y=value,color=name,group=name)) +
#   geom_point() + geom_line() +
#   scale_color_manual("Compared sample",values=c("PCRNeg"="blue","2A"="red")) +
#   theme(axis.text.x=element_text(angle=90))
# g.dist


# exp1 tree  --------------------------------------------------------

samples.compare <- c("1A","1B","1C","1L","1Z")
phy1.tree <- phy1 %>% 
  filter(lbl %in% samples.compare) 

gg <- phy.prepare.ggtree(phy1.tree,sortby=lbl,
                         xmin.tip = 0.6,
                         radius.range = c(1.05, 1.4))

gg$ggtree +
  geom_tile(data=gg$otu,aes(x=x,y=y,fill=otu,alpha=pctseqs),color="gray") +
  geom_segment(data=gg$tax,aes(x=x,xend=gg$x.ring.min,y=y,yend=y),color="gray",linetype="dotted") +
  geom_text(data=gg$samp,aes(x=x,y=gg$angle_to_y(90),label=lbl)) +
  # taxa names
  geom_text(data=gg$tax,aes(x=gg$x.ring.max,y=y,label=Species,hjust=hjust,angle=angle),
            color="dark gray",size=2) +
  scale_alpha_continuous(trans=log_epsilon_trans()) +
  scale_fill_taxonomy(data=gg$otu,fill=otu)

# exp1 tree diffs (circular) ---------------------------------------------------------

samples.compare <- c("1A","1B","1C","1L","1Z")

phy1.tree <- phy1 %>%
  mutate(compare=lbl %in% samples.compare) %>%
  filter(compare|baseline) %>%
  summarize_tax_table(in.baseline=any(baseline & numseqs>0),
                      baseline.pct=pctseqs[baseline])

otu1base <- phy1.tree %>%
  filter(baseline,prune_unused_taxa=FALSE) %>%
  get.otu.melt(filter.zero=FALSE) %>%
  transmute(otu,pctseqs0=pctseqs)
gg <- phy1.tree %>%
  filter(compare,prune_unused_taxa=FALSE) %>%
  phy.prepare.ggtree(sortby=lbl,
                     xmin.tip = 0.6,
                     radius.range = c(1.05, 1.4))
gg$otu <- gg$otu %>%
  left_join(otu1base,by="otu") %>%
  mutate(pct.diff=pctseqs-pctseqs0,
         pct.ratio=pctseqs/pctseqs0,
         log.pct.ratio=log(pct.ratio))

gt.base <- gg$ggtree +
  geom_segment(data=gg$tax,aes(x=x,xend=gg$x.ring.min,y=y,yend=y,
                               color=in.baseline,
                               linetype=in.baseline)) +
  geom_point(data=gg$tax,aes(x=x,y=y,color=in.baseline,
                             size=as.numeric(baseline.pct)),alpha=0.6) +
  geom_text(data=gg$tax,aes(x=gg$x.ring.max,y=y,label=Species,hjust=hjust,angle=angle,
                            color=in.baseline),size=2) +
  geom_text(data=gg$samp,aes(x=x,y=gg$angle_to_y(90),label=lbl)) +
  # taxa names
  scale_size_continuous(transform=log_epsilon_trans(0.01)) +
  scale_color_manual(values=c("TRUE"="pink","FALSE"="dark gray")) +
  scale_linetype_manual(values=c("TRUE"="solid","FALSE"="dotted"))

# regular plot
g1.tree <- gt.base + 
  geom_tile(data=filter(gg$otu,numseqs>0),
            aes(x=x,y=y,fill=otu,alpha=pctseqs),color="gray") +
  geom_text(data=gg$samp,aes(x=x,y=gg$angle_to_y(90),label=lbl)) +
  scale_alpha_continuous(trans=log_epsilon_trans()) +
  scale_fill_taxonomy(data=gg$otu,fill=otu)
g1.tree
# pct.diff
g1.tree.diff <- gt.base + 
  geom_tile(data=filter(gg$otu,numseqs>0),
            aes(x=x,y=y,fill=pct.diff),color="gray") +
  geom_text(data=gg$samp,aes(x=x,y=gg$angle_to_y(90),label=lbl)) +
  scale_fill_gradient2(transform=log_epsilon_trans(0.001),limits=c(-1,1))
g1.tree.diff

# pct.ratio
g1.tree.ratio <- gt.base + 
  geom_tile(data=filter(gg$otu,!is.nan(log.pct.ratio)),
            aes(x=x,y=y,fill=log.pct.ratio),color="gray",size=0.5) +
  geom_text(data=gg$samp,aes(x=x,y=gg$angle_to_y(90),label=lbl)) +
  scale_fill_gradient2(limits=c(-5,5))
g1.tree.ratio

pdf("plots/g1.tree.pdf",width=20,height=20)
g1.tree
g1.tree.diff
g1.tree.ratio
dev.off()

shell.exec("plots/g1.tree.pdf")

# exp1 tree diffs (rectangular) ---------------------------------------------------------

samples.compare <- c("1A","1B","1C","1L","1Z")

phy1.tree <- phy1 %>%
  mutate(compare=lbl %in% samples.compare) %>%
  filter(compare|baseline) %>%
  summarize_tax_table(in.baseline=any(baseline & numseqs>0),
                      baseline.pct=pctseqs[baseline])

otu1base <- phy1.tree %>%
  filter(baseline,prune_unused_taxa=FALSE) %>%
  get.otu.melt(filter.zero=FALSE) %>%
  transmute(otu,pctseqs0=pctseqs)
gg <- phy1.tree %>%
  filter(compare,prune_unused_taxa=FALSE) %>%
  phy.prepare.ggtree(layout="rectangular",
                     sortby=lbl,
                     xmin.tip = 0,
                     radius.range = c(1.05, 1.4))
gg$otu <- gg$otu %>%
  left_join(otu1base,by="otu") %>%
  mutate(pct.diff=pctseqs-pctseqs0,
         pct.ratio=pctseqs/pctseqs0,
         log.pct.ratio=log(pct.ratio))




gt.base <- gg$ggtree +
  expand_limits(x=1) +
  geom_segment(data=gg$tax,aes(x=x,xend=gg$x.ring.min,y=y,yend=y,
                               color=in.baseline,
                               linetype=in.baseline)) +
  geom_point(data=gg$tax,aes(x=x,y=y,color=in.baseline,
                             size=as.numeric(baseline.pct)),alpha=0.6) +
  geom_text(data=gg$tax,aes(x=gg$x.ring.max,y=y,label=Species,hjust=hjust,angle=angle,
                            color=in.baseline),size=2) +
  geom_text(data=gg$samp,aes(x=x,y=gg$angle_to_y(90),label=lbl)) +
  # taxa names
  scale_size_continuous(transform=log_epsilon_trans(0.01)) +
  scale_color_manual(values=c("TRUE"="pink","FALSE"="dark gray")) +
  scale_linetype_manual(values=c("TRUE"="solid","FALSE"="dotted"))

# regular plot
g1.tree <- gt.base + 
  geom_tile(data=filter(gg$otu,numseqs>0),
            aes(x=x,y=y,fill=otu,alpha=pctseqs),color="gray") +
  geom_text(data=gg$samp,aes(x=x,y=gg$angle_to_y(90),label=lbl)) +
  scale_alpha_continuous(trans=log_epsilon_trans()) +
  scale_fill_taxonomy(data=gg$otu,fill=otu)
g1.tree
# pct.diff
g1.tree.diff <- gt.base + 
  geom_tile(data=filter(gg$otu,numseqs>0),
            aes(x=x,y=y,fill=pct.diff),color="gray") +
  geom_text(data=gg$samp,aes(x=x,y=gg$angle_to_y(90),label=lbl)) +
  scale_fill_gradient2(transform=log_epsilon_trans(0.001),limits=c(-1,1))
g1.tree.diff

# pct.ratio
g1.tree.ratio <- gt.base + 
  geom_tile(data=filter(gg$otu,!is.nan(log.pct.ratio)),
            aes(x=x,y=y,fill=log.pct.ratio),color="gray",size=0.5) +
  geom_text(data=gg$samp,aes(x=x,y=gg$angle_to_y(90),label=lbl)) +
  scale_fill_gradient2(limits=c(-5,5))
g1.tree.ratio


pdf("plots/g1.tree.rect.pdf",width=10,height=20)
g1.tree
g1.tree.diff
g1.tree.ratio
dev.off()

shell.exec("plots/g1.tree.rect.pdf")

# exp 2 tree diffs -------------------------------------------------------------

samples.compare <- c("2A","2H","2L","2N","2S","2U")
# samples.compare <- c("2A","2H","2L","2N","2S")
phy2.tree <- phy2 %>%
  mutate(compare=lbl %in% samples.compare) %>%
  filter(compare|baseline) %>%
  summarize_tax_table(in.baseline=any(baseline & numseqs>0),
                      baseline.pct=pctseqs[baseline])

otu2base <- phy2.tree %>%
  filter(baseline,prune_unused_taxa=FALSE) %>%
  get.otu.melt(filter.zero=FALSE) %>%
  transmute(otu,pctseqs0=pctseqs)
gg <- phy2.tree %>%
  filter(compare,prune_unused_taxa=FALSE) %>%
  phy.prepare.ggtree(sortby=lbl,
                     xmin.tip = 0.6,
                     radius.range = c(1.05, 1.4))
gg$otu <- gg$otu %>%
  left_join(otu2base,by="otu") %>%
  mutate(pct.diff=pctseqs-pctseqs0,
         pct.ratio=pctseqs/pctseqs0,
         log.pct.ratio=log(pct.ratio))

gt.base <- gg$ggtree +
  geom_segment(data=gg$tax,aes(x=x,xend=gg$x.ring.min,y=y,yend=y,
                               color=in.baseline,
                               linetype=in.baseline)) +
  geom_point(data=gg$tax,aes(x=x,y=y,color=in.baseline,
                             size=as.numeric(baseline.pct)),alpha=0.6) +
  # taxa names
  geom_text(data=gg$tax,aes(x=gg$x.ring.max,y=y,label=Species,hjust=hjust,angle=angle,
                            color=in.baseline),size=2) +
  scale_size_continuous(transform=log_epsilon_trans(0.01)) +
  scale_color_manual(values=c("TRUE"="pink","FALSE"="dark gray")) +
  scale_linetype_manual(values=c("TRUE"="solid","FALSE"="dotted"))

# regular plot
g2.tree <- gt.base + 
  geom_tile(data=filter(gg$otu,numseqs>0),
            aes(x=x,y=y,fill=otu,alpha=pctseqs),color="gray") +
  geom_text(data=gg$samp,aes(x=x,y=gg$angle_to_y(90),label=lbl)) +
  scale_alpha_continuous(trans=log_epsilon_trans()) +
  scale_fill_taxonomy(data=gg$otu,fill=otu)
g2.tree
# pct.diff
g2.tree.diff <- gt.base + 
  geom_tile(data=filter(gg$otu,numseqs>0),
            aes(x=x,y=y,fill=pct.diff),color="gray") +
  geom_text(data=gg$samp,aes(x=x,y=gg$angle_to_y(90),label=lbl)) +
  scale_fill_gradient2(transform=log_epsilon_trans(0.001),limits=c(-1,1))
g2.tree.diff
# pct.ratio
g2.tree.ratio <- gt.base + 
  geom_tile(data=filter(gg$otu,!is.nan(log.pct.ratio)),
            aes(x=x,y=y,fill=log.pct.ratio),color="gray",size=0.5) +
  geom_text(data=gg$samp,aes(x=x,y=gg$angle_to_y(90),label=lbl)) +
  scale_fill_gradient2(limits=c(-5,5))
g2.tree.ratio





# exp1 tree diffs (old)---------------------------------------------------------


library(glue)

phy1.tree <- phy1 %>%
  phy.collapse() %>%
  mutate(compare=lbl %in% samples.compare) %>%
  filter(compare|baseline)
s1.tree <- phy1.tree %>%
  get.samp()
tr <- phy_tree(phy1.tree)

# ggtree object
gt <- ggtree(tr,layout="circular") %<+% get.tax(phy1.tree)
gd <- gt$data %>%
  filter(isTip) %>%
  mutate(otu=label,
         Species=str_replace_all(Species,"\\[|\\]",""),
         hjust=ifelse(is.between(angle,90,270),1,0),
         angle=ifelse(is.between(angle,90,270),angle+180,angle))

# otu
ydict <- gd %>% select(otu,y)
xlim <- max(gd$x)*c(1.3,1.5)
# samples
xdict <- tibble(lbl=samples.compare) %>% 
  left_join(s1.tree,by="lbl") %>%
  mutate(xring=1:n(),
         x=scales::rescale(xring,to=xlim),
         samplabel=str_glue("{lbl}: time={time}, temp={temp}"))
otu1base <- phy1.tree %>%
  filter(baseline,prune_unused_taxa=FALSE) %>%
  get.otu.melt(filter.zero=FALSE) %>%
  transmute(otu,pctseqs0=pctseqs) %>%
  inner_join(ydict,by="otu")

otu.subset <- phy1.tree %>%
  filter(compare,prune_unused_taxa=FALSE) %>%
  get.otu.melt(filter.zero=TRUE,sample_data=FALSE,tax_data=TRUE) %>%
  left_join(otu1base,by="otu") %>%
  inner_join_replace(ydict,by="otu") %>%
  inner_join_replace(xdict,by="sample") %>%
  mutate(in.baseline=pctseqs0>0,
         pct.diff=pctseqs-pctseqs0,
         pct.ratio=pctseqs/pctseqs0)

gg.tyler.tree <- gt +
  # geom_point2(data=gd,aes(subset=isTip,color=Phylum2)) +
  geom_point(data=gd,aes(color=otu)) +
  # geom_text(data=gd,aes(label=Species,angle=angle,hjust=hjust),size=2.5) +
  # hilight.clade(gd,Family,"Lachnospiraceae",fill.color="#EC9B96",alpha=0.1,xmax=1) +
  scale_color_taxonomy(data=gd,color=otu) +
  expand_limits(x=-0.5)

# gg.tyler.tree.data <- gg.tyler.tree +
#   geom_tile(data=otu.subset,aes(x=x,y=y,fill=otu,alpha=pctseqs),color="darkgray") +
#   scale_alpha_continuous(trans=log_epsilon_trans(0.001)) +
#   scale_fill_taxonomy(data=otu.subset,fill=otu) +
#   geom_text(data=xdict,aes(x=x,y=max(gd$y)/4,label=samplabel),size=3)

gg.tyler.tree.data <- gg.tyler.tree +
  geom_tile(data=otu.subset,aes(x=x,y=y,fill=pct.diff),color="darkgray") +
  geom_tile(data=filter(otu.subset,in.baseline),aes(x=x,y=y,color=in.baseline),fill=NA,color="red") +
  # scale_alpha_continuous(trans=log_epsilon_trans(0.001)) +
  scale_fill_gradient2(transform=log_epsilon_trans(0.001)) +
  geom_text(data=xdict,aes(x=x,y=max(gd$y)/4,label=samplabel),size=3)

gg.tyler.tree.data






# exp2 tree diffs? (old)--------------------------------------------------------


# g2a.alt

samples.compare <- c("2A","2H","2L","2N","2S","2U")

phy2.tree <- phy2 %>%
  phy.collapse() %>%
  mutate(compare=lbl %in% samples.compare) %>%
  filter(compare|baseline)
s2.tree <- phy2.tree %>%
  get.samp()
tr <- phy_tree(phy2.tree)

# ggtree object
gt <- ggtree(tr,layout="circular") %<+% get.tax(phy2.tree)
gd <- gt$data %>%
  filter(isTip) %>%
  mutate(otu=label,
         Species=str_replace_all(Species,"\\[|\\]",""),
         hjust=ifelse(is.between(angle,90,270),1,0),
         angle=ifelse(is.between(angle,90,270),angle+180,angle))

# otu
ydict <- gd %>% select(otu,y)
xlim <- max(gd$x)*c(1.1,1.5)
# samples
xdict <- tibble(lbl=samples.compare) %>% 
  left_join(s2.tree,by="lbl") %>%
  mutate(xring=1:n(),
         x=scales::rescale(xring,to=xlim),
         samplabel=str_glue("{lbl}"))
otu2base <- phy2.tree %>%
  filter(baseline,prune_unused_taxa=FALSE) %>%
  get.otu.melt(filter.zero=FALSE) %>%
  transmute(otu,pctseqs0=pctseqs) %>%
  inner_join(ydict,by="otu")

otu.subset <- phy2.tree %>%
  filter(compare,prune_unused_taxa=FALSE) %>%
  get.otu.melt(filter.zero=TRUE,sample_data=FALSE,tax_data=TRUE) %>%
  left_join(otu2base,by="otu") %>%
  inner_join_replace(ydict,by="otu") %>%
  inner_join_replace(xdict,by="sample") %>%
  mutate(in.baseline=pctseqs0>0,
         pct.diff=pctseqs-pctseqs0,
         pct.ratio=pctseqs/pctseqs0)
gg.tyler.tree <- gt +
  # geom_point2(data=gd,aes(subset=isTip,color=Phylum2)) +
  geom_point(data=gd,aes(color=otu)) +
  # geom_text(data=gd,aes(label=Species,angle=angle,hjust=hjust),size=2.5) +
  # hilight.clade(gd,Family,"Lachnospiraceae",fill.color="#EC9B96",alpha=0.1,xmax=1) +
  scale_color_taxonomy(data=gd,color=otu) +
  expand_limits(x=-0.5)
# gg.tyler.tree.data <- gg.tyler.tree +
#   geom_tile(data=otu.subset,aes(x=x,y=y,fill=otu,alpha=pctseqs),color="darkgray") +
#   scale_alpha_continuous(trans=log_epsilon_trans(0.001)) +
#   scale_fill_taxonomy(data=otu.subset,fill=otu) +
#   geom_text(data=xdict,aes(x=x,y=max(gd$y)/4,label=samplabel),size=3)

gg.tyler.tree.data2 <- gg.tyler.tree +
  geom_tile(data=otu.subset,aes(x=x,y=y,fill=pct.diff),color="darkgray") +
  geom_tile(data=filter(otu.subset,in.baseline),aes(x=x,y=y,color=in.baseline),fill=NA,color="red") +
  # scale_alpha_continuous(trans=log_epsilon_trans(0.001)) +
  scale_fill_gradient2(transform=log_epsilon_trans(0.001)) +
  geom_text(data=xdict,aes(x=x,y=max(gd$y)/4,label=samplabel),size=3)

gg.tyler.tree.data2









# testing exp 1 ----------------------------------------------------------

s1 <- phy1 %>% get.samp(stats=TRUE)
testit1 <- function(sampdata,var) {
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

test1 <- function(yvar,timevar=time,tempvar=temp,sampdata=s1) {
  yvar <- ensym(yvar)
  timevar <- ensym(timevar)
  tempvar <- ensym(tempvar)
  findp <- function(tbl,text) {
    pvals <- tbl %>% filter(term %ilike% text) %>% pull(p.value) 
    if (length(pvals)==0) {
      return(NA)
    } else {
      return(any(pvals<=0.05))
    }
  }
  formula <- rlang::inject(!!yvar ~ !!timevar + !!tempvar)
  model <- lm(formula,data=sampdata)
  regtable <- broom::tidy(model)
  tbl <- tibble(test=as_label(formula),
                model=list(model),
                table=list(regtable),
                temp=findp(regtable,"temp"),
                time=findp(regtable,"time"))
  tbl
}


testall1 <- function(yvar,sampdata=s1) {
  yvar <- ensym(yvar)
  bind_rows(test1(!!yvar,time,temp),
            test1(!!yvar,time.num,temp.num),
            test1(!!yvar,time.rank,temp.rank))
}

test1(InvSimpson) # no diff
test1(qpcr.totalseqs) # *time*
test1(dist_horn) # *time*
test1(dist_pct.bray) # *time*

testall1(InvSimpson) # no diff
testall1(qpcr.totalseqs) # time
testall1(dist_horn) # time
testall1(dist_pct.bray) # time and temp
# testing exp 1 (new) ----------------------------------------------------------

s1 <- phy1 %>% get.samp(stats=TRUE)
testit1 <- function(sampdata,var) {
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

model <- lm(InvSimpson ~ time + temp,data=s1)
tbl <- broom::tidy(model)
tbl


test1 <- function(yvar,timevar=time,tempvar=temp,sampdata=s1) {
  yvar <- ensym(yvar)
  timevar <- ensym(timevar)
  tempvar <- ensym(tempvar)
  findp <- function(tbl,text) {
    pvals <- tbl %>% filter(term %ilike% text) %>% pull(p.value) 
    if (length(pvals)==0) {
      return(NA)
    } else {
      return(any(pvals<=0.05))
    }
  }
  formula <- rlang::inject(!!yvar ~ !!timevar + !!tempvar)
  model <- lm(formula,data=sampdata)
  regtable <- broom::tidy(model)
  tbl <- tibble(test=as_label(formula),
                model=list(model),
                table=list(regtable),
                temp=findp(regtable,"temp"),
                time=findp(regtable,"time"))
  tbl
}


testall1 <- function(yvar,sampdata=s1) {
  yvar <- ensym(yvar)
  bind_rows(test1(!!yvar,time,temp),
            test1(!!yvar,time.num,temp.num),
            test1(!!yvar,time.rank,temp.rank))
}

test1(InvSimpson) # no diff
test1(qpcr.totalseqs) # *time*
test1(dist_horn) # *time*
test1(dist_pct.bray) # *time*

testall1(InvSimpson) # no diff
testall1(qpcr.totalseqs) # time
testall1(dist_horn) # time
testall1(dist_pct.bray) # time and temp


# testing exp 2 -----------------------------------------------------------


s2 <- phy2 %>% get.samp(stats=TRUE)
test2 <- function(yvar,timevar=time,uvvar=uv,heatvar=heat,sampdata=s2) {
  yvar <- ensym(yvar)
  timevar <- ensym(timevar)
  uvvar <- ensym(uvvar)
  heatvar <- ensym(heatvar)
  findp <- function(tbl,text) {
    pvals <- tbl %>% filter(term %ilike% text) %>% pull(p.value) 
    if (length(pvals)==0) {
      return(NA)
    } else {
      return(any(pvals<=0.05))
    }
  }
  formula <- rlang::inject(!!yvar ~ !!timevar + !!uvvar + !!heatvar)
  model <- lm(formula,data=sampdata)
  regtable <- broom::tidy(model)
  tbl <- tibble(test=as_label(formula),
                model=list(model),
                table=list(regtable),
                time=findp(regtable,"time"),
                heat.75c=findp(regtable,"heat75C"),
                heat.autoclave=findp(regtable,"autoclave"),
                uv.regular=findp(regtable,"uvUV$"),
                uv.dna=findp(regtable,"dna"))
  tbl
}

testall2 <- function(yvar,sampdata=s2) {
  yvar <- ensym(yvar)
  bind_rows(test2(!!yvar,time,uv,heat),
            test2(!!yvar,time.rank,uv,heat),
            test2(!!yvar,time,uv,heat.autoclave),
            test2(!!yvar,time.rank,uv,heat.autoclave))
}

test2(InvSimpson) # no diff
test2(qpcr.totalseqs) # *autoclave* and *uv.dna*
test2(dist_horn) # *autoclave* and *uv.dna*
test2(dist_pct.bray) # *time*, *autoclave* and *uv.dna*

# test all versions
testall2(InvSimpson) # no diff
testall2(qpcr.totalseqs) # time and temp
testall2(dist_horn) # no diff
testall2(dist_pct.bray) # time and temp







# permanova exp 2 ---------------------------------------------------------

dist2 <- calc.distance(phy2,"horn")
s2 <- get.samp(phy2)

adonis2(dist2 ~ time + uv + heat, data=s2)
# adonis2(dist2 ~ time + uv.sample + uv.dna + heat.75C + heat.autoclave, data=s2)

# best?
adonis2(dist2 ~ time + 
          uv.sample + uv.dna + 
          heat.75C + heat.autoclave, data=s2)

adonis2(dist2 ~ time + 
          time:uv.sample + time:uv.dna + 
          time:heat.75C + time:heat.autoclave, data=s2)




# (old) fig 1A: experiment 1 stackplot -------------------------------------------------------------------

phy1a <- phy1 %>%
  mutate_sample_data(temp=fct_recode(temp,"-80C"="n/a"))

otu1 <- phy1a %>% 
  get.otu.melt()
s1 <- phy1a %>% get.samp(stats=TRUE)

width <- 0.75
g1a <- ggplot() +
  # geom_taxonomy(data=otu1,aes(x=letter.rev, y=pctseqs, fill=otu, label=Species),width=width) +
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
                 background_x=list(element_blank(),element_blank(),element_rect(),
                                   element_rect(),element_rect(),element_rect()),
                 background_y=list(element_blank(),element_rect(),element_rect(),
                                   element_rect(),element_rect())
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
  geom_col(data=s1,aes(x=stage(letter.rev,after_stat=x+(width+littlewidth)/2),y=s1$dist_mean.horn),
           width=littlewidth,fill="blue") +
  geom_text(data=s1,aes(x=stage(letter.rev,after_stat=x+(width+littlewidth)/2),y=s1$dist_mean.horn,
                        label=str_glue("mean-Horn={sprintf('%.3f',s1$dist_mean.horn)}")
                        # label=str_glue("Horn={sprintf('%.3f',dist_horn)}\nUnfold-Horn={sprintf('%.3f',dist_unfold.horn)}")
                        # label=str_glue("Unfold-Horn={sprintf('%.3f',dist_unfold.horn)}")
  ),hjust=0,size=3)
g1a.alt

# pdf("plots/fig1a - exp1 stackplot.pdf",width=10,height=8)
# g1a.alt
# dev.off()
# 
# shell.exec("plots/fig1a - exp1 stackplot.pdf")




# (old) fig 1B: PCOA exp 1 ------------------------------------------------------------


### add phy.other data ###



phy.part1 <- phy.tyler %>% 
  filter(experiment==1) %>% 
  mutate_sample_data(origin="Our Study") %>%
  phy.use.refseq.as.taxanames()
phy.part2 <- phy.others %>% 
  select(otu,Superkingdom,Phylum,Class,Order,Family,Genus,Species) %>%
  mutate_sample_data(origin="Other") %>%
  phy.use.refseq.as.taxanames()
phy.together <- phy.combine(phy.part1,phy.part2)
rm(phy.part1,phy.part2)


asdf <- function(method,distance,...) {
  ord <- phy.ordinate(phy.together,method=method,distance=distance,...)
  opts <- enexprs(...)
  opts_lbl <- paste(names(opts),opts,sep="=",collapse=", ")
  title <- str_glue("{method} {distance} {opts_lbl}")
  ggplot(ord$data,aes(x=axis1,y=axis2,color=origin)) +
    geom_point() +
    theme(aspect.ratio=1) +
    ggtitle(title)
}

asdf("PCoA","horn") + asdf("PCoA","pct.bray")


### phy1 by itself ###
asdf2 <- function(method,distance,...) {
  ord <- phy.ordinate(phy1,method=method,distance=distance,...)
  ord$data <- ord$data %>% arrange(desc(time))
  opts <- enexprs(...)
  opts_lbl <- paste(names(opts),opts,sep="=",collapse=", ")
  title <- str_glue("{method} {distance} {opts_lbl}")
  ggplot(ord$data,aes(x=axis1,y=axis2)) +
    geom_point(aes(color=time),size=3) +
    geom_text(aes(label=time),size=2,vjust=1) +
    theme(aspect.ratio=1) +
    ggtitle(title)
}


asdf2("PCoA","horn") + asdf2("PCoA","pct.bray")




asdf3 <- function(method,distance,...) {
  
  
  ord <- phy.ordinate(phy.tyler,method=method,distance=distance,...)
  ord$data <- ord$data %>% 
    mutate(color=ifelse(experiment==1,as.character(time),"xanchor"),
           label=paste(time,temp,sep="|"),
           label=ifelse(experiment==1,label,"xanchor"))
  opts <- enexprs(...)
  opts_lbl <- paste(names(opts),opts,sep="=",collapse=", ")
  title <- str_glue("{method} {distance} {opts_lbl}")
  ggplot(ord$data,aes(x=axis1,y=axis2)) +
    geom_point(aes(color=color),size=4) +
    geom_text(aes(label=label),size=2,vjust=1) +
    theme(aspect.ratio=1) +
    ggtitle(title)
}
asdf3("PCoA","horn") +
  asdf3("PCoA","pct.bray")


### phy1a, with exp 2 baseline ###
phy1a <- phy.tyler %>% 
  filter(experiment==1 | sample==sample.comparator) %>%
  mutate_sample_data(exp2anchor=experiment==2,
                     label=paste(time,temp,sep="|"),
                     color=time)

asdf4 <- function(method,distance,...) {
  ord <- phy.ordinate(phy1a,method=method,distance=distance,...)
  opts <- enexprs(...)
  opts_lbl <- paste(names(opts),opts,sep="=",collapse=", ")
  title <- str_glue("{method} {distance} {opts_lbl}")
  ggplot(ord$data,aes(x=axis1,y=axis2)) +
    geom_point(aes(color=color),size=4) +
    geom_text(aes(label=label),size=2,vjust=1) +
    theme(aspect.ratio=1) +
    ggtitle(title)
}

asdf4("PCoA","horn") +
  asdf4("PCoA","pct.bray")





phy.part1 <- phy.tyler %>% 
  filter(experiment==1) %>% 
  mutate_sample_data(origin="Our Study") %>%
  phy.use.refseq.as.taxanames()
phy.part2 <- phy.others %>% 
  select(otu,Superkingdom,Phylum,Class,Order,Family,Genus,Species) %>%
  mutate_sample_data(origin="Other") %>%
  phy.use.refseq.as.taxanames()
phy.together <- phy.combine(phy.part1,phy.part2)
rm(phy.part1,phy.part2)


asdf <- function(method,distance,...) {
  ord <- phy.ordinate(phy.together,method=method,distance=distance,...)
  opts <- enexprs(...)
  opts_lbl <- paste(names(opts),opts,sep="=",collapse=", ")
  title <- str_glue("{method} {distance} {opts_lbl}")
  ggplot(ord$data,aes(x=axis1,y=axis2,color=origin)) +
    geom_point() +
    theme(aspect.ratio=1) +
    ggtitle(title)
}

asdf("PCoA","horn") + asdf("PCoA","pct.bray")


### phy1 by itself ###
asdf2 <- function(method,distance,...) {
  ord <- phy.ordinate(phy1,method=method,distance=distance,...)
  ord$data <- ord$data %>% arrange(desc(time))
  opts <- enexprs(...)
  opts_lbl <- paste(names(opts),opts,sep="=",collapse=", ")
  title <- str_glue("{method} {distance} {opts_lbl}")
  ggplot(ord$data,aes(x=axis1,y=axis2)) +
    geom_point(aes(color=time),size=3) +
    geom_text(aes(label=time),size=2,vjust=1) +
    theme(aspect.ratio=1) +
    ggtitle(title)
}


asdf2("PCoA","horn") + asdf2("PCoA","pct.bray")




asdf3 <- function(method,distance,...) {
  
  
  ord <- phy.ordinate(phy.tyler,method=method,distance=distance,...)
  ord$data <- ord$data %>% 
    mutate(color=ifelse(experiment==1,as.character(time),"xanchor"),
           label=paste(time,temp,sep="|"),
           label=ifelse(experiment==1,label,"xanchor"))
  opts <- enexprs(...)
  opts_lbl <- paste(names(opts),opts,sep="=",collapse=", ")
  title <- str_glue("{method} {distance} {opts_lbl}")
  ggplot(ord$data,aes(x=axis1,y=axis2)) +
    geom_point(aes(color=color),size=4) +
    geom_text(aes(label=label),size=2,vjust=1) +
    theme(aspect.ratio=1) +
    ggtitle(title)
}
asdf3("PCoA","horn") +
  asdf3("PCoA","pct.bray")


### phy1a, with exp 2 baseline ###
phy1a <- phy.tyler %>% 
  filter(experiment==1 | sample==sample.comparator) %>%
  mutate_sample_data(exp2anchor=experiment==2,
                     label=paste(time,temp,sep="|"),
                     color=time)

asdf4 <- function(method,distance,...) {
  ord <- phy.ordinate(phy1a,method=method,distance=distance,...)
  opts <- enexprs(...)
  opts_lbl <- paste(names(opts),opts,sep="=",collapse=", ")
  title <- str_glue("{method} {distance} {opts_lbl}")
  ggplot(ord$data,aes(x=axis1,y=axis2)) +
    geom_point(aes(color=color),size=4) +
    geom_text(aes(label=label),size=2,vjust=1) +
    theme(aspect.ratio=1) +
    ggtitle(title)
}

asdf4("PCoA","horn") +
  asdf4("PCoA","pct.bray")




# (old) fig 2A: experiment 2 stackplot ------------------------------------------

otu2 <- phy2 %>% 
  # for speed
  phy.collapse.bins() %>%
  get.otu.melt() %>%
  group_by(treatment,time) %>%
  mutate(sample.number=as.numeric(factor(sample))) %>%
  ungroup() 
s2 <- phy2 %>% get.samp(stats=TRUE) %>%
  mutate(qpcr.label=ifelse(is.na(qpcr.totalseqs),"undetectable",short_number(qpcr.totalseqs)))

width <- 0.75
g2a <- ggplot() +
  geom_taxonomy(data=otu2,aes(x=letter,y=pctseqs,fill=otu,label=Species),width=width) +
  # geom_taxonomy(data=otu2,aes(x=letter,y=pctseqs,fill=otu),width=width) +
  geom_col(data=filter(s2,baseline),aes(x=letter,y=1),
           width=width,linewidth=0.75,linetype="longdash",color="blue",fill=NA) +
  geom_text(data=s2,aes(x=letter,y=0,label=lbl),hjust=1) +
  expand_limits(y=-0.08) +
  facet_nested(treatmentlabel+treatment~timelabel+time,
               strip=strip_nested( #overarching facets are blank background
                 background_x=list(element_blank(),element_rect(),element_rect(),element_rect()),
                 background_y=list(element_blank(),element_rect(),element_rect(),element_rect(),
                                   element_rect(),element_rect(),element_rect(),element_rect())
               )) +
  scale_fill_taxonomy(data=otu2,fill=otu,
                      # guide=guide_taxonomy(ncol=4,downwards=TRUE),
                      tax.palette = pal) +
  theme(aspect.ratio=0.3,
        # legend.position="bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank()) +
  labs(fill="Bacterial Taxa") +
  coord_flip()
g2a

# with beta diversity on top
littlewidth <- 0.1
xtop <- (width+littlewidth)/2
g2a.alt <- g2a +
  geom_col(data=s2,aes(x=stage(letter,after_stat=x+xtop),y=dist_horn),width=littlewidth,fill="steelblue") +
  geom_text(data=s2,aes(x=stage(letter,after_stat=x+xtop),y=dist_horn,
                        label=str_glue("Horn={sprintf('%.3f',dist_horn)}")),hjust="inward",size=3)
# geom_text(data=s2,aes(x=stage(letter,after_stat=x-xtop),y=0,
#                       label=str_glue("qPCR={qpcr.label}")),hjust="inward",size=3) +
# geom_text(data=s2,aes(x=stage(letter,after_stat=x-xtop),y=1,
#                       label=str_glue("ntaxa={Observed}")),hjust="inward",size=3)

g2a.alt


pdf("plots/fig2a - exp2 stackplot.pdf",width=10,height=8)
g2a.alt
dev.off()

shell.exec("plots/fig2a - exp2 stackplot.pdf")
# (old) fig2A: errant microbes ---------------------------------------------------------



base <- phy2a %>% filter(sample=="TY.1_D0_NT") %>% get.otu.melt(sample_data=FALSE) %>% select(otu,taxid)
ctrl <- phy2a %>% filter(lbl=="PCRNeg") %>% get.otu.melt(sample_data=FALSE) %>% select(otu,taxid)
otu2a <- get.otu.melt(phy2a) %>%
  mutate(taxid.in.base=taxid %in% base$taxid,
         taxid.in.ctrl=taxid %in% ctrl$taxid,
         taxid.status=case_when(
           taxid.in.base & taxid.in.ctrl ~ "both",
           taxid.in.base & !taxid.in.ctrl ~ "in 2A",
           !taxid.in.base & taxid.in.ctrl ~ "in PCRNeg",
           TRUE ~ "neither"),
         Species=str_replace_all(Species,"\\[|\\]",""))
samp2a <- otu2a %>%
  group_by(!!!syms(sample_variables(phy2a))) %>%
  summarize(pct.taxid.base=mean(taxid.in.base),
            pct.taxid.ctrl=mean(taxid.in.ctrl),
            .groups="drop") %>%
  mutate(x=as.numeric(lbl3),
         xmin=x-0.475,xmax=x+0.475,
         qpcr.lbl=case_when(
           lbl3=="PCRNeg" ~ "(n/a)",
           is.na(qpcr.totalseqs) ~ "(undetectable)",
           TRUE ~ NA_character_))
samp2a.long <- samp2a %>% pivot_longer(cols=c(pct.taxid.base,pct.taxid.ctrl)) %>%
  mutate(name=case_match(name,"pct.taxid.base"~"2A","pct.taxid.ctrl"~"PCRNeg"))
treat.groups <- samp2a.long %>% 
  group_by(treatment) %>%
  summarize(xmin=min(xmin),xmax=max(xmax),
            .groups="drop") %>%
  mutate(b.xmin=xmin+0.05,b.xmax=xmax-0.05,
         b.xmin2=b.xmin+0.05,b.xmax2=b.xmax-0.05)

group.shading <- geom_rect(data=treat.groups,aes(xmin=b.xmin,xmax=b.xmax,ymin=-Inf,ymax=Inf),alpha=0.1)
group.brackets <- function(y) {
  geom_bracket(data=treat.groups,
               aes(x=b.xmin2,xend=b.xmax2,y=y,label=treatment),
               fontsize=3,tip="bare")
}

eps <- 1000
g.qpcr <- ggplot(samp2a) + 
  geom_col(aes(x=lbl3,y=qpcr.totalseqs)) + 
  geom_text(aes(x=lbl3,y=0,label=qpcr.lbl),hjust=0,angle=90) +
  group.shading + 
  group.brackets(-eps) +
  scale_y_continuous("Total abundance by 16S qPCR",trans=log_epsilon_trans(eps),
                     labels=pretty_power10)




phy2a <- phy2a %>%
  select_sample_data(-starts_with("dist")) %>%
  add_dist("mean.horn",comparator="TY.1_D0_NT",varname=dist_2A) %>%
  add_dist("mean.horn",comparator="PCRNeg4",varname=dist_PCRNeg)
s2a <- phy2a %>% get.samp()
s2a.long <- s2a %>% pivot_longer(cols=c(dist_2A,dist_PCRNeg)) %>%
  mutate(name=case_match(name,"dist_2A"~"2A","dist_PCRNeg"~"PCRNeg"))



g.asv.compare <- ggplot(samp2a.long) +
  geom_point(aes(x=lbl3,y=value,color=name,group=name)) + 
  geom_line(aes(x=lbl3,y=value,color=name,group=name)) +
  group.shading + 
  group.brackets(-0.1) +
  # # dashed boxes
  geom_rect(data=filter(samp2a,lbl %in% c("2A","PCRNeg")),
            aes(xmin=xmin,xmax=xmax,ymin=0,ymax=1,color=lbl),
            linetype="longdash",linewidth=0.75,fill=NA,show.legend = FALSE) +
  scale_color_manual("Compared sample",values=c("PCRNeg"="blue","2A"="red")) +
  scale_y_continuous("Percent ASV with similarity\nto compared sample",label=scales::label_percent()) +
  theme(axis.text.x=element_text(angle=90))

g.tax.compare <- ggplot(otu2a) +
  geom_taxonomy(aes(x=lbl3,y=pctseqs,label=Species,fill=otu),tax.palette=pal,label.split=TRUE,fontsize=3) +
  # dashed boxes
  geom_rect(data=filter(samp2a,lbl %in% c("2A","PCRNeg")),
            aes(xmin=xmin,xmax=xmax,ymin=0,ymax=1,color=lbl),
            linetype="longdash",linewidth=0.75,fill=NA,show.legend = FALSE) +
  group.shading +
  scale_color_manual("",values=c("PCRNeg"="blue","2A"="red")) +
  scale_x_discrete("Sample") +
  scale_y_continuous("Bacterial Composition\nRelative abundance",expand=FALSE) +
  scale_fill_taxonomy(name="Bacterial taxa",data=otu2a,fill=otu,tax.palette=pal) +
  theme(axis.text.x=element_text(hjust=0,vjust=0.5,angle=90),
        axis.text.y=element_blank(),
        axis.ticks=element_blank())




g.fig2a <- gg.stack2(g.qpcr,g.asv.compare,g.tax.compare)
g.fig2a



# g.dist <- ggplot(s2a.long,aes(x=lbl3,y=value,color=name,group=name)) +
#   geom_point() + geom_line() +
#   scale_color_manual("Compared sample",values=c("PCRNeg"="blue","2A"="red")) +
#   theme(axis.text.x=element_text(angle=90))
# g.dist



# (old) fig 1b: pcoa  -----------------------------------------------------------

library(ggrepel)

ord1 <- phy.ordinate(phy1,method="PCoA",distance="mean.horn")
g.fig1b <- ggplot(ord1$data,aes(x=axis1,y=axis2,fill=temp,color=temp,shape=temp,label=lbl3)) +
  geom_point(size=4) + 
  geom_text(aes(label=days),size=3,color="black") +
  geom_text_repel(size=3,show.legend = FALSE) +
  scale_shape_manual(values=c("n/a"=21, "-80C"=22, "-20C"=23, "4C"=24, "room temp"=25)) +
  theme(aspect.ratio=1)

g.fig1b

pdf("plots/fig1b - exp1 pcoa.pdf",width=7,height=7)
g.fig1b
dev.off()
shell.exec("plots/fig1b - exp1 pcoa.pdf")


g.fig1b %>% copy.to.clipboard.gg(width=7,height=7)



# (old) fig 2C? exp 2 step compare ------------------------------------------------------


# g2a.alt

samples.compare <- c("2A","2H","2L","2N","2S","2U")

phy2.compare <- phy2 %>%
  mutate_sample_data(compare=lbl %in% samples.compare,
                     lbl.compare=str_glue("{lbl} (vs. {lbl.comparator})")) %>%
  filter(compare|baseline)
otu2base <- phy2.compare %>% 
  filter(baseline,prune_unused_taxa=FALSE) %>%
  get.otu.melt(filter.zero=FALSE) %>%
  transmute(otu,pctseqs0=pctseqs)
otu2compare <- phy2.compare %>% 
  filter(compare,prune_unused_taxa=FALSE) %>% 
  get.otu.melt(filter.zero=FALSE) %>%
  left_join(otu2base,by="otu") %>%
  filter(pctseqs>0|pctseqs0>0) %>%
  group_by(sample) %>%
  arrange(desc(pctseqs0),desc(pctseqs)) %>%
  mutate(col=row_number(),
         extra=pctseqs0==0 & pctseqs>0) %>%
  ungroup()

s2compare <- phy2.compare %>% 
  filter(compare,prune_unused_taxa=FALSE) %>% 
  get.samp(stats=TRUE) %>%
  mutate(label=str_glue("qpcr.totalseqs={pretty_number(qpcr.totalseqs)}\nInvSimpson={pretty_number(InvSimpson,digits=3)}\nHorn={sprintf('%.3f',dist_horn)}"))


g2.asv <- ggplot() +
  geom_col(data=otu2compare,aes(x=col,y=pctseqs,fill=otu)) +
  geom_step(data=otu2compare,aes(x=col,y=pctseqs0),direction="mid") +
  geom_text(data=s2compare,aes(x=Inf,y=Inf,label=label),hjust=1,vjust=1,color="blue") +
  geom_rect(data=s2compare,aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,linetype=baseline),
            fill=NA,color="blue",show.legend=FALSE) +
  scale_linetype_manual(values=c("TRUE"="longdash","FALSE"=NA)) +
  geom_bracket(data=filter(otu2compare,extra),
               aes(x=col,y=ave(pctseqs,sample,FUN=max),
                   fontsize=3,label="unique\nASVs"),tip="square") + 
  scale_fill_taxonomy(name="Bacterial Taxa",data=otu2compare,fill=otu) +
  scale_y_continuous("Relative Abundance",trans=log_epsilon_trans(0.001)) +
  facet_wrap(~lbl.compare,ncol=1) +
  theme(#aspect.ratio=1,
    axis.text = element_blank(), axis.ticks = element_blank(), 
    axis.title = element_blank(), panel.background = element_blank())

g2.asv








# contam ------------------------------------------------------------------



xx <- isNotContaminant(phy2a,neg="is.neg.control",detailed=TRUE)
xx %>% glimpse()

isContaminant(phy2a,neg="is.neg.control")
# freq prev p.freq p.prev  p contaminant
isContaminant()
isNotContaminant()




base <- phy2a %>% filter(sample=="TY.1_D0_NT") %>% get.otu.melt(sample_data=FALSE) %>% select(otu,taxid)
ctrl <- phy2a %>% filter(lbl=="PCRNeg") %>% get.otu.melt(sample_data=FALSE) %>% select(otu,taxid)

phy2ax <- phy2a %>% 
  mutate_tax_table(not.contam=isNotContaminant(.,neg="is.neg.control"),
                   taxid.in.base=taxid %in% base$taxid,
                   taxid.in.ctrl=taxid %in% ctrl$taxid,
                   asv.in.base=otu %in% base$otu,
                   asv.in.ctrl=otu %in% ctrl$otu,
                   taxid.status=case_when(
                     taxid.in.base & taxid.in.ctrl ~ "both",
                     taxid.in.base & !taxid.in.ctrl ~ "in 2A",
                     !taxid.in.base & taxid.in.ctrl ~ "in PCRNeg",
                     TRUE ~ "neither"),
                   NULL) %>%
  summarize_sample_data(pct.taxid.base=mean(as.logical(taxid.in.base)),
                        pct.taxid.ctrl=mean(as.logical(taxid.in.ctrl)),
                        pct.notcontam=mean(as.logical(not.contam)),
                        pct.contam=mean(!as.logical(not.contam)),
                        pct.otu.base=mean(as.logical(asv.in.base)),
                        pct.otu.ctrl=mean(as.logical(asv.in.ctrl)),
                        pctseqs.base=sum(pctseqs[as.logical(taxid.in.base)]),
                        pctseqs.ctrl=sum(pctseqs[as.logical(taxid.in.ctrl)]),
                        pctseqs.contam=sum(pctseqs[!as.logical(not.contam)]),
                        pctseqs.noncontam=sum(pctseqs[as.logical(not.contam)]),
                        filter.zero=TRUE) %>%
  add_dist("horn",comparator="TY.1_D0_NT",varname=dist_2A) %>%
  add_dist("mean.horn",comparator="TY.1_D0_NT",varname=dist_2A_meanhorn) %>%
  add_dist("horn",comparator="PCRNeg4",varname=dist_PCRNeg)

otu2ax <- get.otu.melt(phy2ax)
samp2ax <- get.samp(phy2ax)
otu.nocontam <- phy2ax %>% filter(not.contam=="TRUE") %>%
  get.otu.melt(filter.zero=FALSE)


samp2ax.long <- samp2ax %>%
  pivot_longer(cols=c(qpcr.totalseqs,
                      dist_2A,dist_PCRNeg,
                      dist_2A_meanhorn,
                      pctseqs.base,
                      pctseqs.ctrl,
                      pctseqs.contam,
                      pctseqs.noncontam,
                      pct.contam,
                      pct.notcontam,
                      # pct.otu.base,
                      # pct.otu.ctrl,
                      pct.taxid.ctrl,pct.taxid.base))

ggplot(samp2ax.long,aes(x=lbl3,y=value,group=name,color=name)) + 
  geom_point() + geom_line() +
  theme(axis.text.x=element_text(hjust=0,vjust=0.5,angle=90),
        axis.ticks=element_blank()) +
  facet_wrap(~name,scales="free_y",nrow=4)


# tax regular
g.tax <- ggplot(otu2ax) +
  geom_taxonomy(aes(x=lbl3,y=pctseqs,
                    label=Species,fill=otu,
                    alpha=not.contam),
                tax.palette=pal,
                label.split=TRUE,fontsize=3) +
  scale_x_discrete("Sample") +
  scale_y_continuous("Bacterial Composition\nRelative abundance",expand=FALSE) +
  scale_fill_taxonomy(name="Bacterial taxa", data=otu2ax, fill=otu,tax.palette=pal) +
  theme(axis.text.x=element_text(hjust=0,vjust=0.5,angle=90),
        axis.text.y=element_blank(),
        axis.ticks=element_blank())
g.tax

# tax without contam
g.tax.nocontam <- ggplot(otu.nocontam) +
  geom_taxonomy(aes(x=lbl3,y=pctseqs,
                    label=Species,fill=otu),
                tax.palette=pal,
                label.split=TRUE,fontsize=3) +
  scale_x_discrete("Sample") +
  scale_y_continuous("Bacterial Composition\nRelative abundance",expand=FALSE) +
  scale_fill_taxonomy(name="Bacterial taxa", data=otu.nocontam, 
                      fill=otu,tax.palette=pal) +
  theme(axis.text.x=element_text(hjust=0,vjust=0.5,angle=90),
        axis.text.y=element_blank(),
        axis.ticks=element_blank())
g.tax.nocontam


g.qpcr <- ggplot(samp2ax,aes(x=lbl3,y=qpcr.totalseqs)) + 
  geom_col() +
  scale_y_continuous(trans=log_epsilon_trans(1000),
                     labels=pretty_power10)
g.qpcr



breakdown <- otu2ax %>% 
  group_by(sample,lbl3,not.contam) %>%
  summarize(pctseqs=sum(pctseqs),
            n.asv=n_distinct(otu),
            .groups="drop")
g.breakdown <- ggplot(breakdown) +
  geom_col(aes(x=lbl3,y=pctseqs,fill=not.contam))



gg.stack(g.breakdown,g.qpcr,g.tax,g.tax.nocontam)











# lefse -------------------------------------------------------------------

phy1.lefse = phy.tyler %>%
  phy.collapse(taxranks=c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  filter(experiment==1) %>%
  mutate(time.group=fct_recode(time,"Day 0-3"="Day 0","Day 0-3"="Day 3","Day 8-11"="Day 8","Day 8-11"="Day 11"),
         time.group2=fct_recode(time,"Day 0-3"="Day 0","Day 0-3"="Day 3"))
s1 <- phy1.lefse %>% get.samp()

lda <- lda.effect(phy1.lefse,class="time.group")
lda.plot(lda) + lda.clado(lda)

lda.b <- lda.effect(phy1.lefse,class="time.group2")
lda.plot(lda.b) + lda.clado(lda.b)




