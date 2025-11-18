
# load data ---------------------------------------------------------------

library(yingtools2)
library(tidyverse)
library(phyloseq)
library(ggtree)
library(glue)
rm(list=ls())

load("data/phy.tyler.RData")
phy.tyler <- phy.tyler %>% 
  select(otu,Superkingdom,Phylum,Class,Order,Family,Genus,Species) 

phy.tyler <- phy.tyler %>% 
  mutate(temp=ifelse(sample %in% c("1A","1B"),"-80C",temp),
         temp=factor(temp,levels=c("-80C", "-20C", "4C", "room temp")),
         time=factor(time,levels=c("day 0", "day 3", "day 6", "day 8", "day 9", "day 11")),
         treatment=factor(treatment,levels=c("none", "75C", "UV", "75C+UV", "autoclave", "autoclave+UV", "UV DNA")),
         sample_number=order(order(treatment,temp,time)),
         letter=str_extract(sample,"(?<=^[0-9]{1,2})[AB]"),
         sample2=paste(treatment,temp,storage,time,sample_number,sep="|"),
         sample2=fct_reordern(sample2,sample_number),
         sample.comparator=ifelse(experiment==1,"1A","TY.1_D0_NT"))

s <- phy.tyler %>% get.samp()

phy.species <- phy.tyler %>% 
  phy.collapse(taxranks=c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

# tax stack plot ----------------------------------------------------------

otu <- phy.species %>% get.otu.melt() %>%
  mutate(qpcr.numseqs=qpcr.totalseqs*pctseqs)


g.taxstack <- ggplot(data=otu,aes(x=sample2,y=pctseqs,fill=otu,label=Species)) +
  geom_taxonomy() +
  theme(axis.text.x=element_text(angle=90))
g.taxstack


g.16sqpcr <- ggplot(data=otu,aes(x=sample2,y=qpcr.numseqs,fill=otu,label=Species)) +
  geom_taxonomy() +
  theme(axis.text.x=element_text(angle=90)) + 
  scale_y_continuous(trans=log_epsilon_trans(10000))
g.16sqpcr


# tax stack faceted -------------------------------------------------------


otu1 <- phy.tyler %>% filter(experiment==1) %>%
  get.otu.melt()
g.exp1 <- ggplot(otu1,aes(x=letter,y=pctseqs,fill=otu,label=Species)) + 
  geom_taxonomy() +
  facet_grid(temp~time)
otu2 <- phy.tyler %>% filter(experiment==2) %>%
  get.otu.melt() %>%
  group_by(treatment,time) %>%
  mutate(sample.number=as.numeric(factor(sample))) %>%
  ungroup() 

g.exp2 <- ggplot(otu2,aes(x=factor(1),y=pctseqs,fill=otu,label=Species)) + 
  geom_taxonomy() +
  facet_grid(treatment~time)

g.exp1 / g.exp2


# show distances -------------------------------------------------------------------

show_dist <- function(method) {
  ss1 <- s %>% filter(experiment==1) %>%
    mutate(unfold.horn=map2_dbl(sample,sample.comparator,~{
      physub <- phy.tyler %>% filter(sample %in% c(.x,.y))
      dist <- calc.distance(physub,method)[1]
      dist
    }))
  g.dist1 <- ggplot(ss1,aes(x=letter,y=unfold.horn)) + 
    geom_col() +
    facet_grid(temp~time) +
    ggtitle(method)
  
  ss2 <- s %>% filter(experiment==2) %>%
    mutate(unfold.horn=map2_dbl(sample,sample.comparator,~{
      physub <- phy.tyler %>% filter(sample %in% c(.x,.y))
      dist <- calc.distance(physub,method)[1]
      dist
    }))
  
  g.dist2 <- ggplot(ss2,aes(x=factor(1),y=unfold.horn)) + 
    geom_col() +
    facet_grid(treatment~time)
  
  g.dist1 / g.dist2
  
}

show_dist("unfold.horn")
show_dist("horn")

show_dist("bray")

# ordination and hclust --------------------------------------------------------------


dist.unfold.horn <- calc.distance(phy.tyler,"unfold.horn")
hc <- hclust(dist.unfold.horn)
sample.order <- hc$labels[hc$order]
ord <- ordinate(phy.tyler,method="MDS",distance=dist.unfold.horn)
pca <- ord$vectors %>% as.data.frame() %>% 
  rownames_to_column("sample") %>%
  left_join(s,by="sample") 
g.pca <- ggplot(pca,aes(x=Axis.1,y=Axis.2)) + 
  geom_point(size=5,alpha=0.7) + 
  coord_fixed() + 
  theme(aspect.ratio=1)
g.pca


# phylo tree --------------------------------------------------------------

selected.samples <- c("1A","TY.12_D9_UV","TY.9_D9_AC")

phy.tree <- phy.species %>%
  filter(sample %in% selected.samples)

# phylo tree
tr <- phy_tree(phy.tree)
# ggtree object
gt <- ggtree(tr,layout="circular") %<+% get.tax(phy.tree)

gd <- gt$data %>%
  filter(isTip) %>%
  mutate(otu=label,
         Phylum2=fct_lump_n(fct_infreq(Phylum),n=4,other_level="Other Phyla"),
         Species=str_replace_all(Species,"\\[|\\]",""),
         hjust=ifelse(is.between(angle,90,270),1,0),
         angle=ifelse(is.between(angle,90,270),angle+180,angle))


ydict <- gd %>% select(otu,y)
xlim <- max(gd$x)*c(1.3,1.5)
xdict <- tibble(sample=selected.samples) %>% 
  left_join(s,by="sample") %>%
  mutate(xring=1:n(),
         x=scales::rescale(xring,to=xlim))

otu.subset <- phy.tree %>%
  get.otu.melt(filter.zero=FALSE) %>%
  inner_join(ydict,by="otu") %>%
  inner_join(xdict,by="sample") 

gg.tyler.tree <- gt + 
  # geom_point2(data=gd,aes(subset=isTip,color=Phylum2)) +
  geom_point(data=gd,aes(color=otu)) +
  # geom_text(data=gd,aes(label=Species,angle=angle,hjust=hjust),size=2.5) + 
  # hilight.clade(gd,Family,"Lachnospiraceae",fill.color="#EC9B96",alpha=0.1,xmax=1) +
  scale_color_taxonomy(data=gd,color=otu) +
  expand_limits(x=-0.5)

gg.tyler.tree

gg.tyler.tree.data <- gg.tyler.tree +
  geom_tile(data=otu.subset,aes(x=x,y=y,fill=otu,alpha=pctseqs),color="darkgray") +
  scale_alpha_continuous(trans=log_epsilon_trans(0.001)) +
  scale_fill_taxonomy(data=otu.subset,fill=otu) +
  geom_text(data=xdict,aes(x=x,y=max(gd$y)/4,label=sample2),size=3)

gg.tyler.tree.data

pdf("plots/tyler.tree.pdf",width=20,height=20)
gg.tyler.tree.data
dev.off()

shell.exec("plots/tyler.tree.pdf")


# compare taxa ------------------------------------------------------------

compare_samps <- function(samps) {
  physub <- phy.tyler %>%
    filter(sample %in% samps)
  otusub <- physub %>% get.otu.melt(filter.zero = FALSE) %>%
    mutate(sign=ifelse(sample==samps[1],1,-1),
           y=sign*pctseqs) %>%
    group_by(otu) %>%
    mutate(y1=max(y),y2=min(y)) %>%
    ungroup() %>% 
    mutate(x=fct_reordern(otu,y1,y2))
  ssub <- physub %>% get.samp()
  unfold.horn <- calc.distance(physub,"unfold.horn")[1]
  title <- str_glue("{ssub$sample[1]}={pretty_number(ssub$nseqs[1])}\n{ssub$sample[2]}={pretty_number(ssub$nseqs[2])}
unfold.horn={round(unfold.horn,3)}")
  g.compare <- ggplot(otusub) + 
    geom_col(aes(x=x,y=y,fill=otu),show.legend=FALSE) +
    expand_limits(y=c(-1,1)) +
    scale_y_continuous(trans=log_epsilon_trans(0.001)) +
    scale_fill_taxonomy(data=otusub,fill=otu) +
    ggtitle(title)
  g.compare
}

compare_samps2 <- function(samps) {
  physub <- phy.tyler %>%
    filter(sample %in% samps)
  otusub <- physub %>% get.otu.melt(filter.zero = FALSE) %>%
    mutate(sign=ifelse(sample==samps[1],1,-1),
           y=sign*pctseqs*qpcr.totalseqs) %>%
    group_by(otu) %>%
    mutate(y1=max(y),y2=min(y)) %>%
    ungroup() %>% 
    mutate(x=fct_reordern(otu,y1,y2))
  ssub <- physub %>% get.samp()
  unfold.horn <- calc.distance(physub,"unfold.horn")[1]
  title <- str_glue("{ssub$sample[1]}={pretty_number(ssub$nseqs[1])}\n{ssub$sample[2]}={pretty_number(ssub$nseqs[2])}
unfold.horn={round(unfold.horn,3)}")
  g.compare <- ggplot(otusub) + 
    geom_col(aes(x=x,y=y,fill=otu),show.legend=FALSE) +
    expand_limits(y=c(-1,1)) +
    scale_y_continuous(trans=log_epsilon_trans(10000)) +
    scale_fill_taxonomy(data=otusub,fill=otu) +
    ggtitle(title)
  g.compare
}

# exp 1 
samples <- c("1A","1B")
compare_samps(samples)
compare_samps2(samples)

# exp 1, -80 storage
samples <- c("1A","13A.-80.D11")
compare_samps(samples)
compare_samps2(samples)

# exp 1, room temp
samples <- c("1A","10A.RT.D11")
compare_samps(samples)
compare_samps2(samples)


# exp 1, room temp + day 11
samples <- c("1A","13A.-80.D11")
compare_samps(samples)
compare_samps2(samples)



# exp 2, none vs. 75C
samples <- c("TY.1_D0_NT","TY.4_D0_75C")
compare_samps(samples)
compare_samps2(samples)


samples <- c("TY.1_D0_NT","TY.21_D9_DNA_UV")
compare_samps(samples)
compare_samps2(samples)



samples <- c("TY.1_D0_NT","TY.15_D9_75C_UV")
compare_samps(samples)
compare_samps2(samples)



exp1_plots <- s %>% filter(experiment==1) %>%
  mutate(time=fct_drop(time),
         temp.letter=str_glue("{temp} ({letter})"),
         temp.letter=fct_reordern(temp.letter,temp,letter)) %>%
  select(experiment,temp.letter,time,temp,sample,sample.comparator) %>%
  complete(temp.letter,time) %>%
  mutate(g=map2(sample,sample.comparator,~{
    if (!is.na(.y)) {
      gg <- compare_samps(c(.x,.y))
    } else {
      gg <- patchwork::plot_spacer()
    }
    return(gg)
  }),
  g2=map2(sample,sample.comparator,~{
    if (!is.na(.y)) {
      gg <- compare_samps2(c(.x,.y))
    } else {
      gg <- patchwork::plot_spacer()
    }
    return(gg)
  })) %>%
  arrange(temp.letter,time)

pdf("plots/exp1_compare.pdf",width=12,height=20)
patchwork::wrap_plots(exp1_plots$g,nrow=8,ncol=4)
dev.off()

pdf("plots/exp1_qpcrcompare.pdf",width=12,height=20)
patchwork::wrap_plots(exp1_plots$g2,nrow=8,ncol=4)
dev.off()

shell.exec("plots/exp1_compare.pdf")
shell.exec("plots/exp1_qpcrcompare.pdf")


exp2_plots <- s %>% filter(experiment==2) %>%
  select(time,treatment,sample,sample.comparator) %>%
  mutate(g=map2(sample,sample.comparator,~{
    if (!is.na(.y)) {
      gg <- compare_samps(c(.x,.y))
    } else {
      gg <- patchwork::plot_spacer()
    }
    return(gg)
  }),
  g2=map2(sample,sample.comparator,~{
    if (!is.na(.y)) {
      gg <- compare_samps2(c(.x,.y))
    } else {
      gg <- patchwork::plot_spacer()
    }
    return(gg)
  })) %>%
  arrange(time,treatment)

# patchwork::wrap_plots(exp2_plots$g,nrow=7,ncol=3)

pdf("plots/exp2_compare.pdf",width=12,height=20)
patchwork::wrap_plots(exp2_plots$g,nrow=7,ncol=3)
dev.off()

pdf("plots/exp2_qpcrcompare.pdf",width=12,height=20)
patchwork::wrap_plots(exp2_plots$g2,nrow=7,ncol=3)
dev.off()

shell.exec("plots/exp2_compare.pdf")
shell.exec("plots/exp2_qpcrcompare.pdf")






exp1_qpcrplots <- s %>% filter(experiment==1) %>%
  mutate(time=fct_drop(time),
         temp.letter=str_glue("{temp} ({letter})"),
         temp.letter=fct_reordern(temp.letter,temp,letter)) %>%
  select(experiment,temp.letter,time,temp,sample,sample.comparator) %>%
  complete(temp.letter,time) %>%
  mutate(g=map2(sample,sample.comparator,~{
    if (!is.na(.y)) {
      gg <- compare_samps2(c(.x,.y))
    } else {
      gg <- patchwork::plot_spacer()
    }
    return(gg)
  })) %>%
  arrange(temp.letter,time)

# patchwork::wrap_plots(exp1_plots$g,nrow=8,ncol=4)

pdf("plots/exp1_qpcrcompare.pdf",width=12,height=20)
patchwork::wrap_plots(xx$g,nrow=8,ncol=4)
dev.off()
shell.exec("plots/exp1_qpcrcompare.pdf")


exp2_qpcrplots <- s %>% filter(experiment==2) %>%
  select(time,treatment,sample,sample.comparator) %>%
  mutate(g=map2(sample,sample.comparator,~{
    if (!is.na(.y)) {
      gg <- compare_samps2(c(.x,.y))
    } else {
      gg <- patchwork::plot_spacer()
    }
    return(gg)
  })) %>%
  arrange(time,treatment)

# patchwork::wrap_plots(exp2_plots$g,nrow=7,ncol=3)

pdf("plots/exp2_qpcrcompare.pdf",width=12,height=20)
patchwork::wrap_plots(exp2_qpcrplots$g,nrow=7,ncol=3)
dev.off()
shell.exec("plots/exp2_qpcrcompare.pdf")









# pw <- calc.distance(phy.tyler,"unfold.horn") %>% get.pairwise()
# 
# pw %>% arrange(desc(dist))
# 
# s <- get.samp(phy.tyler)



# 
# depict <- function(subsamps) {
#   physub <- phy %>% filter(sample %in% subsamps)
#   bray <- calc.distance(physub,"bray")[1] %>% round(3)
#   bray.pct <- calc.distance(physub,"pct.bray")[1] %>% round(3)
#   horn <- calc.distance(physub,"horn")[1] %>% round(3)
#   taxhorn <- calc.distance(physub,"mean.horn")[1] %>% round(3)
#   # taxhorn <- calc.mean.distance2(physub,"horn") %>% round(3)
#   unfoldhorn <- calc.distance(physub,"unfold.horn")[1] %>% round(3)
#   
#   
#   title <- str_glue("bray={bray}\npct.bray={bray.pct}\nhorn={horn},\ntaxhorn={taxhorn}\nunfoldhorn={unfoldhorn}")  
#   otusub <- get.otu.melt(physub,filter.zero = FALSE) %>%
#     mutate(sign=ifelse(sample==subsamps[1],1,-1),
#            y=sign*pctseqs) %>%
#     group_by(otu) %>%
#     mutate(y1=pctseqs[which(sample==subsamps[1])[1]],
#            y2=pctseqs[which(sample==subsamps[2])[1]]) %>%
#     ungroup() %>% 
#     arrange(y1,-y2) %>%
#     mutate(x=fct_inorder(otu))
#   
#   g.compare <- ggplot(otusub,aes(x=x,y=y,fill=otu)) + geom_col() +
#     scale_fill_taxonomy(data=otusub,fill=otu) +
#     scale_y_continuous(trans=log_epsilon_trans(0.001))
#   g.tax <- ggplot(otusub,aes(x=sample,y=pctseqs,fill=otu,label=Species)) + 
#     geom_taxonomy(width=0.5,show.ribbon = TRUE) +
#     ggtitle(title)
#   g.tax / g.compare
# }



g

g








