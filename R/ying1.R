library(phyloseq)
library(yingtools2)
library(tidyverse)

load("data/phy.tyler.RData")




phy <- phy.tyler %>% 
  select(Superkingdom,Phylum,Class,Order,Family,Genus,Species,otu) %>%
  mutate(group=case_when(
    uv=="no UV" & heat=="no heat" ~ "normal",
    TRUE ~ paste(uv,heat,sep="+")
  ))

s <- phy %>% get.samp(stats=TRUE)













# ordination and hclust --------------------------------------------------------------



dist.unfold.horn <- calc.distance(phy,"unfold.horn")
hc <- hclust(dist.unfold.horn)
sample.order <- hc$labels[hc$order]
ord <- ordinate(phy,method="MDS",distance=dist.unfold.horn)
pca <- ord$vectors %>% as.data.frame() %>% 
  rownames_to_column("sample") %>%
  left_join(s,by="sample") 
g.pca <- ggplot(pca,aes(x=Axis.1,y=Axis.2,color=group)) + geom_point(size=5,alpha=0.7) + coord_fixed() + theme(aspect.ratio=1)
g.pca


otu <- phy %>% get.otu.melt() %>%
  mutate(sample.hc=factor(sample,levels=sample.order),
         absolute.abundance=pctseqs * qpcr.totalseqs)

g.hc <- ggplot(otu,aes(x=sample.hc,y=pctseqs,fill=otu,label=Species)) + geom_taxonomy() + 
  theme(axis.text.x=element_text(angle=90))
g.hc

s %>% glimpse
g.tax <- ggplot(otu,aes(x=sample.hc,y=pctseqs,fill=otu,label=Species)) + geom_taxonomy() + 
  theme(axis.text.x=element_text(angle=90)) +
  facet_grid(group ~ .)
g.tax
g.tax2 <- ggplot(otu,aes(x=sample.hc,y=absolute.abundance,fill=otu,label=Species)) + geom_taxonomy() + 
  theme(axis.text.x=element_text(angle=90)) + 
  scale_y_continuous(trans=log_epsilon_trans(10000))
g.tax2

otu %>% filter(experiment==2) %>%
  ggplot(aes(x=factor(1),y=absolute.abundance,fill=otu)) + 
  geom_taxonomy() + 
  scale_y_continuous(trans=log_epsilon_trans(10000)) +
  facet_grid(treatment ~ time)





