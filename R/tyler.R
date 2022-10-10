
# load data ---------------------------------------------------------------

library(yingtools2)
library(tidyverse)
library(phyloseq)
library(ggfittext)

rm(list=ls())
load("data/phy.tyler.RData")

s <- get.samp(phy.tyler) %>%
  mutate(temp=factor(temp,levels=c("room temp", "4C", "-20C", "-80C")),
         time=factor(time,levels=c("day 0", "day 3", "day 6", "day 8", "day 9", "day 11"))) %>%
  arrange(treatment,temp,time,uv) %>%
  mutate(sample_number=row_number(),
         sample2=paste(treatment,temp,storage,time,sample_number,sep="|"),
         sample2=fct_inorder(sample2))
sample_data(phy.tyler) <- s %>% set.samp()

# overview stack ----------------------------------------------------------

phy.tyler.species <- phy.collapse(phy.tyler)

otu <- phy.tyler.species %>% get.otu.melt()
s <- get.samp(phy.tyler.species)
pal <- get.yt.palette2(otu)

g.stack <- ggplot() +
  geom_col(data=otu,aes(x=sample2,y=pctseqs,fill=Species),position="fill") +
  geom_bar_text(data=otu,
                aes(x=sample2,y=pctseqs,fill=Species,label=Species),
                position="fill",
                reflow=TRUE,lineheight=0.75,place="centre",angle=-90,
                color="black",contrast=FALSE,size=7,min.size=7) +
  geom_text(data=s,aes(x=sample2,y=1.05,label=short_number(nseqs)),angle=-90) +
  scale_fill_manual(values=pal) +
  theme(legend.position="none",axis.text.x=element_text(angle=-90,vjust=0.5,hjust=0.5)) +
  facet_grid(. ~ treatment,space="free_x",scales="free_x")

pdf("plots/taxstack.pdf",width=20,height=12)
g.stack
dev.off()
# shell.exec("plots\\taxstack.pdf")




# PCA ---------------------------------------------------------------------


library(vegan)
library(ggrepel)

phy.pca <- phy.tyler
s.pca <- get.samp(phy.pca) %>%
  mutate(uv.temp=paste(uv,temp,sep="|"),
         label=paste(treatment,time,sep="|"),
         ac.or.uv=uv=="UV"|temp=="autoclave")
pcadata <- ordinate(phy.pca, method = "NMDS", distance = "horn") %>%
  scores(display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(s.pca, by = "sample")

g.pca <- ggplot(pcadata, aes(x=NMDS1, y=NMDS2, label=label, color=treatment)) +
  geom_point(size=5,alpha=0.5) +
  geom_text_repel(color="black",max.overlaps=20) +
  theme(aspect.ratio = 1)
g.pca

pdf("plots/pca.pdf",width=12,height=12)
g.pca
dev.off()
# shell.exec("plots\\pca.pdf")


# diversity ---------------------------------------------------------------

phy.species <- phy.tyler %>% phy.collapse()

s <- get.samp(phy.tyler,stats=TRUE)
# s$InvSimpson

s %>% arrange(desc(InvSimpson))



# separate bar (large) ------------------------------------------------------------

phy.tyler.species <- phy.tyler %>% phy.collapse()

otu <- phy.tyler.species %>% get.otu.melt(filter.zero = FALSE) %>%
  mutate(Phylum2=fct_lump_n(Phylum,n=4))
pal <- get.yt.palette2(otu)


g.bar <- ggplot() +
  geom_col(data=otu,aes(x=Species,y=pctseqs,fill=Species)) +
  scale_fill_manual(values=pal) +
  scale_y_continuous(trans=log_epsilon_trans()) +
  theme(legend.position="none",axis.text.x=element_text(angle=-90,vjust=0.5,hjust=0.5),
        strip.text.y = element_text(angle=0)) +
  facet_grid(sample2 ~ Phylum2,space="free_x",scales="free_x")

pdf("plots/g.bar.pdf",width=50,height=30)
g.bar
dev.off()
shell.exec("plots\\g.bar.pdf")





