
library(yingtools2)
library(tidyverse)
library(phyloseq)
library(ggfittext)
library(ggtree)
library(glue)
rm(list=ls())
load("data/phy.tyler.RData")

# collapse down to species
phy.tyler.species <- phy.collapse(phy.tyler,taxranks=c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

# phylo tree
tr <- phy_tree(phy.tyler.species)
# ggtree object
gt <- ggtree(tr,layout="circular") %<+% get.tax(phy.tyler.species)

gd <- gt$data %>%
  filter(isTip) %>%
  mutate(Phylum2=fct_lump_n(fct_infreq(Phylum),n=4,other_level="Other Phyla"),
         Species=str_replace_all(Species,"\\[|\\]",""),
         hjust=ifelse(is.between(angle,90,270),1,0),
         angle=ifelse(is.between(angle,90,270),angle+180,angle))

gg.tyler.tree <- gt + 
  geom_point2(data=gd,aes(subset=isTip,color=Phylum2)) +
  geom_text2(data=gd,aes(subset=isTip,label=Species,angle=angle,hjust=hjust),size=2.5) + 
  # hilight.clade(gd,Family,"Lachnospiraceae",fill.color="#EC9B96",alpha=0.1,xmax=1) +
  expand_limits(x=-0.5)


phy.subset <- phy.tyler.species %>%
  filter(sample %in% c("1A",""))



s$sample

otu <- phy.tyler.species %>% get.otu.melt()
s <- phy.tyler.species %>% get.samp()



s


gg.tyler.tree



pdf("plots/tyler.tree.pdf",width=25,height=25)
gg.tyler.tree
dev.off()


shell.exec("plots/tyler.tree.pdf")

# phy.tyler
# 
# 
# s <- get.samp(phy.tyler) %>%
#   mutate(temp=factor(temp,levels=c("room temp", "4C", "-20C", "-80C")),
#          time=factor(time,levels=c("day 0", "day 3", "day 6", "day 8", "day 9", "day 11"))) %>%
#   arrange(treatment,temp,time,uv) %>%
#   mutate(sample_number=row_number(),
#          sample2=paste(treatment,temp,storage,time,sample_number,sep="|"),
#          sample2=fct_inorder(sample2))
# sample_data(phy.tyler) <- s %>% set.samp()
# 
# # overview stack ----------------------------------------------------------
# 
# phy.tyler.species <- phy.collapse(phy.tyler)
# 
# otu <- phy.tyler.species %>% get.otu.melt()
# s <- get.samp(phy.tyler.species)
# pal <- get.yt.palette2(otu)
# 
# g.stack <- ggplot() +
#   geom_col(data=otu,aes(x=sample2,y=pctseqs,fill=Species),position="fill") +
#   geom_bar_text(data=otu,
#                 aes(x=sample2,y=pctseqs,fill=Species,label=Species),
#                 position="fill",
#                 reflow=TRUE,lineheight=0.75,place="centre",angle=-90,
#                 color="black",contrast=FALSE,size=7,min.size=7) +
#   geom_text(data=s,aes(x=sample2,y=1.05,label=short_number(nseqs)),angle=-90) +
#   scale_fill_manual(values=pal) +
#   theme(legend.position="none",axis.text.x=element_text(angle=-90,vjust=0.5,hjust=0.5)) +
#   facet_grid(. ~ treatment,space="free_x",scales="free_x")
# 
# pdf("plots/taxstack.pdf",width=20,height=12)
# g.stack
# dev.off()
# # shell.exec("plots\\taxstack.pdf")
# 
# 
# 
