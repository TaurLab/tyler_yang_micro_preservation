
# load data ---------------------------------------------------------------

library(yingtools2)
library(tidyverse)
library(phyloseq)
library(ggtree)
library(glue)
rm(list=ls())

load("data/phy.tyler.RData")

phy.tyler <- phy.tyler %>% 
  mutate(temp=factor(temp,levels=c("-80C", "-20C", "4C", "room temp")),
         time=factor(time,levels=c("day 0", "day 3", "day 6", "day 8", "day 9", "day 11")),
         sample_number=order(order(treatment,temp,time)),
         sample2=paste(treatment,temp,storage,time,sample_number,sep="|"),
         sample2=fct_reordern(sample2,sample_number),
         group=xxxx)

s <- phy.tyler %>% get.samp()
phy.species <- phy.tyler %>% phy.collapse(taxranks=c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))



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


# 
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
  geom_point2(data=gd,aes(subset=isTip,color=otu)) +
  geom_text2(data=gd,aes(subset=isTip,label=Species,angle=angle,hjust=hjust),size=2.5) + 
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
