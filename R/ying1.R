
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

s <- phy.tyler %>% get.samp()o

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


# pairwise comparisons
# obligate anaerobe 
# testing

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




ggplot(otu2,aes(x=factor(1),y=pctseqs,fill=otu,label=Species)) + 
  geom_taxonomy() +
  facet_grid(treatment~time)

# bels and returns a list or data frame of character vectors. Each input column corresponds to one factor. 
# Thus there will be more than one with vars(cyl, am). Each output column gets displayed as one separate line 
# in the strip label. This function should inherit from the "labeller" S3 class for compatibility with labeller(). 
# You can use different labeling functions for different kind of labels, for example use label_parsed() for 
# formatting facet labels. label_value() is used by default, check it for more details and pointers to other options.

# t <- tibble(x=c(0,3,3,8,8,11,11),
#        dist=c(0.01,0.2,0.21,0.4,0.42,0.6,0.7))
# 
# 
# c(0.01,0.2,0.21,0.4,0.42,0.6,0.7) %>% rank()
# 
# x <- runif(40)
# y <- 2*x+runif(n=40,min=1,max=20)
# t <- tibble(x,y)
# t %>% ggplot(aes(x=x,y=y)) + geom_point()
# cor.test(x,y,method="spearman")
# permanova, manova


# tax stack 3 (more polished) -------------------------------------------------------------

phy1 <- phy.tyler %>% 
  filter(experiment==1)
otu1 <- phy1 %>% 
  get.otu.melt() %>%
  mutate(time=recode2(time,c("day 0"="day 0\n(baseline)")))

pal <- yt.palette3 %>% setNames(c("Bacteroidota (phylum)", "Lachnospiraceae (family)", "Oscillospiraceae (family)", 
                                  "Eubacteriales (order)", "Actinomycetota (phylum)", "Enterococcus (genus)", 
                                  "Streptococcus (genus)", "Staphylococcus (genus)", "Lactobacillus (genus)", 
                                  "Pseudomonadota (phylum)", "Other Bacteria"))

g1 <- ggplot(otu1) +
  geom_taxonomy(aes(x = letter, y = pctseqs, fill = otu, label = Species)) +
  facet_grid(temp ~ time) +
  scale_fill_taxonomy(data=otu1,fill=otu,
                      guide=guide_taxonomy(ncol=4,downwards=TRUE),
                      tax.palette = pal) +
  
  theme(aspect.ratio=2,
        legend.key.size = unit(0.85,"lines"), # smaller if necessary
        legend.byrow=TRUE,
        legend.title.position="right",
        legend.title=element_text(angle=-90),
        legend.text.position="bottom",
        legend.text=element_text(angle=-90,hjust=0,vjust=0.5),
        legend.key.spacing.x=unit(0,"lines"),
        legend.key.spacing.y=unit(0.5,"lines"),
        # legend.position = c(0.125,0.375),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank()) +
  labs(fill="Bacterial Taxa")
g1

g1t <- gtable::as.gtable(g1)
g1t$widths[8] <- grid::unit(35,"points")
plot(g1t)



phy2 <- phy.tyler %>%
  filter(experiment==2)
otu2 <- phy2 %>% get.otu.melt() %>%
  group_by(treatment,time) %>%
  mutate(sample.number=as.numeric(factor(sample))) %>%
  ungroup() 
g2
g2 <- ggplot(otu2,aes(x=factor(1),y=pctseqs,fill=otu,label=Species)) + 
  geom_taxonomy() +
  facet_grid(treatment~time) +
  scale_fill_taxonomy(data=otu2,fill=otu,
                      guide=guide_taxonomy(ncol=4,downwards=TRUE),
                      tax.palette = pal) +
  theme(aspect.ratio=2,
        legend.key.size = unit(0.85,"lines"), # smaller if necessary
        legend.byrow=TRUE,
        legend.title.position="right",
        legend.title=element_text(angle=-90),
        legend.text.position="bottom",
        legend.text=element_text(angle=-90,hjust=0,vjust=0.5),
        legend.key.spacing.x=unit(0,"lines"),
        legend.key.spacing.y=unit(0.5,"lines"),
        # legend.position = c(0.125,0.375),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank()) +
  labs(fill="Bacterial Taxa")
g2
g2t <- gtable::as.gtable(g2)

gridExtra::grid.arrange(g1t,g2,nrow=1)


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











# lefse -------------------------------------------------------------------


phy.add <- phy.tyler %>%
  add.abundance(totalseqs=TRUE,totalbacteroidetes=Phylum=="Bacteroidetes") %>%
  add.abundance(pctbacteroidetes=Phylum=="Bacteroidetes",denom=TRUE)

# store relative abundance
phy.pct <- phy.tyler %>%
  transform_sample_counts(function(x) x/sum(x))
phy.tyler

phy.pct


transform_sample_counts2 <- function(phy,fun) {
  fun <- as_mapper(fun)
  otu <- get.otu(phy,as.matrix=TRUE)
  samp <- get.samp(phy)
  
  # split otu_table into list of columns
  otu.list <- lapply(seq_len(ncol(otu)), function(i) otu[,i])
  samp.list <- lapply(seq_len(nrow(samp)), function(i) samp[i,])
  new.otu <- map2(otu.list,samp.list,function(x,s) {
    # add s variables to environment
    list2env(s, envir = environment())
    new.x <- eval_tidy(fun(x), env = environment())
    return(new.x)
  })
}

phy.qpcr <- phy.tyler %>%
  transform_sample_counts2(function(x) {
    x/sum(x)*qpcr.totalseqs
  })




phy.pct %>% get.otu.melt() %>% glimpse
phy.add %>% get.samp() %>% glimpse
phy.pct <- phy.tyler

phyloseq::
otu <- get.otu.melt(phy.tyler)

yingtools2:::sam


otu$numseqs

otu

otu[1:15,1:15]

phy.lefse <- phy.tyler %>%
  filter(experiment==2) %>%
  mutate(disrupt=ifelse(treatment %ilike% "UV|autoclave","disrupt","none"),
         time01=ifelse(time %in% c("day 6","day 9"),"after","before"))

phy.lefse %>% get.samp() %>% glimpse

lda <- lda.effect(phy.lefse,class="disrupt",subclass="time01")









# testing experiment 1 -----------------------------------------------------------------



add_dist1 <- function(sampdata,method,sample0="1A",phy=phy1) {
  varname <- paste0("dist_",method)
  dist <- phy.tyler %>%
    calc.distance(method=method)  
  pw <- dist %>% get.pairwise() %>%
    filter(sample1==sample0|sample2==sample0) %>%
    mutate(sample=coalesce(na_if(sample1,sample0),sample2)) %>%
    select(sample,!!varname:=dist)  
  sampdata %>% left_join(pw,by="sample")
}
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

phy1 <- phy.tyler %>% filter(experiment==1)
s1 <- phy1 %>% 
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

testit1(s1,InvSimpson)
testit1(s1,qpcr.totalseqs) # time
testit1(s1,dist_mean.bray) # temp
testit1(s1,dist_horn) # time
testit1(s1,dist_mean.horn) # temp+time
testit1(s1,dist_unfold.horn) # temp+time





# testing experiment 2 ----------------------------------------------------



add_dist2 <- function(sampdata,method,sample0="TY.1_D0_NT",phy=phy2) {
  varname <- paste0("dist_",method)
  dist <- phy.tyler %>%
    calc.distance(method=method)  
  pw <- dist %>% get.pairwise() %>%
    filter(sample1==sample0|sample2==sample0) %>%
    mutate(sample=coalesce(na_if(sample1,sample0),sample2)) %>%
    select(sample,!!varname:=dist)  
  sampdata %>% left_join(pw,by="sample")
}


phy2 <- phy.tyler %>% filter(experiment==2) 
s2 <- phy2 %>% get.samp(stats=TRUE) %>%
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
         qpcr.totalseqs=coalesce(qpcr.totalseqs,1))


testit2 <- function(sampdata,var) {
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
           uv.regular=findp(model,"uvUV"),
           uv.dna=findp(model,"dna"))
  return(tbl)
}
plotit2 <- function(sampdata,var,eps=0) {
  var <- enquo(var)
  ggplot(sampdata) + 
    geom_col(aes(x=time,y=!!var,fill=treatment)) +
    facet_grid(uv~heat) + 
    scale_y_continuous(trans=log_epsilon_trans(eps))
}



testit2(s2,InvSimpson)
testit2(s2,qpcr.totalseqs) # heat.autoclave
testit2(s2,dist_mean.bray) # uv.reg, uv.dna
testit2(s2,dist_horn) # time, heat.autoclave, uv.reg, uv.dna
testit2(s2,dist_mean.horn) # time, heat.autoclave, uv.reg, uv.dna
testit2(s2,dist_unfold.horn) # time, heat.autoclave, uv.reg, uv.dna

plotit2(s2,InvSimpson)
plotit2(s2,qpcr.totalseqs,eps=1000) 
plotit2(s2,dist_mean.bray) 
plotit2(s2,dist_horn)  
plotit2(s2,dist_mean.horn) 
plotit2(s2,dist_unfold.horn)  


# visualize comparisons ---------------------------------------------------


depict0 <- function(subsamps,phy=phy.tyler) {
  
  physub <- phy %>% filter(sample %in% subsamps)
  bray <- calc.distance(physub,"bray")[1] %>% round(3)
  bray.pct <- calc.distance(physub,"pct.bray")[1] %>% round(3)
  horn <- calc.distance(physub,"horn")[1] %>% round(3)
  taxhorn <- calc.distance(physub,"mean.horn")[1] %>% round(3)
  # taxhorn <- calc.mean.distance2(physub,"horn") %>% round(3)
  unfoldhorn <- calc.distance(physub,"unfold.horn")[1] %>% round(3)
  
  title <- str_glue("bray={bray}\npct.bray={bray.pct}\nhorn={horn},\ntaxhorn={taxhorn}\nunfoldhorn={unfoldhorn}")  
  otusub <- get.otu.melt(physub,filter.zero = FALSE) %>%
    mutate(sign=ifelse(sample==subsamps[1],1,-1),
           y=sign*pctseqs) %>%
    group_by(otu) %>%
    mutate(y1=pctseqs[which(sample==subsamps[1])[1]],
           y2=pctseqs[which(sample==subsamps[2])[1]]) %>%
    ungroup() %>% 
    arrange(y1,-y2) %>%
    mutate(x=fct_inorder(otu))
  
  g.compare <- ggplot(otusub,aes(x=x,y=y,fill=otu)) + geom_col() +
    scale_fill_taxonomy(data=otusub,fill=otu) +
    scale_y_continuous(trans=log_epsilon_trans(0.001))
  g.tax <- ggplot(otusub,aes(x=sample,y=pctseqs,fill=otu,label=Species)) + 
    geom_taxonomy(width=0.5,show.ribbon = TRUE) +
    ggtitle(title)
  g.tax / g.compare
}


depict <- function(subsamps,phy=phy.tyler) {
  
  physub <- phy %>% filter(sample %in% subsamps)
  # bray <- calc.distance(physub,"bray")[1] %>% round(3)
  # bray.pct <- calc.distance(physub,"pct.bray")[1] %>% round(3)
  # horn <- calc.distance(physub,"horn")[1] %>% round(3)
  # taxhorn <- calc.distance(physub,"mean.horn")[1] %>% round(3)
  # taxhorn <- calc.mean.distance2(physub,"horn") %>% round(3)
  # unfoldhorn <- calc.distance(physub,"unfold.horn")[1] %>% round(3)
  
  # title <- str_glue("bray={bray}\npct.bray={bray.pct}\nhorn={horn},\ntaxhorn={taxhorn}\nunfoldhorn={unfoldhorn}")  
  otusub <- get.otu.melt(physub,filter.zero = FALSE) %>%
    mutate(sign=ifelse(sample==subsamps[1],1,-1),
           y=sign*pctseqs) %>%
    group_by(otu) %>%
    mutate(y1=pctseqs[which(sample==subsamps[1])[1]],
           y2=pctseqs[which(sample==subsamps[2])[1]]) %>%
    ungroup() %>% 
    arrange(y1,-y2) %>%
    mutate(x=fct_inorder(otu))

  s.counts <- otusub %>% 
    mutate(sample=factor(sample,levels=subsamps)) %>%
    arrange(sample) %>%
    group_by(sample) %>%
    summarize(seqs=sum(numseqs),
              .groups="drop") %>%
    mutate(text=str_glue("{sample}={pretty_number(seqs)}"))
  title <- paste(s.counts$text,collapse="; ")

  g.compare <- ggplot(otusub,aes(x=x,y=y,fill=otu)) + 
    geom_col(show.legend=FALSE) +
    scale_fill_taxonomy(data=otusub,fill=otu) +
    scale_y_continuous(trans=log_epsilon_trans(0.001)) +
    ggtitle(title)
  # g.tax <- ggplot(otusub,aes(x=sample,y=pctseqs,fill=otu,label=Species)) + 
  #   geom_taxonomy(width=0.5,show.ribbon = TRUE) +
  #   ggtitle(title)
  # g.tax / g.compare
  g.compare
}

s <- phy.tyler %>% get.samp()
#samp samp

subsamps <- c("1A","1B")
depict(subsamps)

depict0(subsamps)




subsamps <- c("1A","TY.1_D0_NT")
depict(subsamps)

subsamps <- c("1A","TY.8_D6_AC")
depict(subsamps)


depict3 <- function(sample1,sample2,phy=phy.tyler) {
  # phy <- phy.tyler
  # subsamps <- c("1A","TY.1_D0_NT")
  subsamps <- c(sample1,sample2)
  physub <- phy %>% 
    filter(sample %in% subsamps)
  otu <- physub %>% 
    get.otu.melt() %>%
    mutate(first=sample==sample1,
           # y2=ifelse(first,pctseqs,-pctseqs),
           y=pctseqs,
           x=fct_reordern(otu,first,pctseqs),
           x=as.numeric(x))
  # g <- ggplot(otu) +
  #   geom_col(aes(x=x,y=y2,fill=otu)) +
  #   scale_fill_taxonomy(data=otu,fill=otu) +
  #   scale_y_continuous(trans=log_epsilon_trans(0.001))
  # g
  
  
  otu1 <- otu %>% filter(first) %>%
    bind_rows(tibble(x=range(.$x)+c(-1,1),y=c(0,0)))
  otu2 <- otu %>% filter(!first)
  
  
  ggplot() +
    geom_col(data=otu2,aes(x=x,y=y,fill=otu)) +
    geom_step(data=otu1,aes(x=x,y=y),direction="mid") +
    scale_fill_taxonomy(data=otu,fill=otu) +
    scale_y_continuous(trans=log_epsilon_trans(0.001)) +
    theme()
}

depict3("1A","1B")
depict3("1A","TY.1_D0_NT")





sample1 <- "1A"
phy1 <- phy.tyler %>% 
  filter(experiment==1)

otu0 <- phy1 %>%
  filter(sample==sample1) %>%
  get.otu.melt() %>%
  arrange(pctseqs) %>%
  transmute(otu,pctseqs0=pctseqs)
otu <- phy1 %>% 
  filter(letter=="B"|sample==sample1) %>%
  mutate(time=ifelse(sample==sample1,"day -1",as.character(time))) %>%
  get.otu.melt(filter.zero=FALSE) %>%
  left_join(otu0,by="otu") %>%
  mutate(pctseqs0=coalesce(pctseqs0,0)) %>%
  filter(pctseqs>0|pctseqs0>0) %>%
  group_by(sample) %>%
  arrange(-pctseqs0,pctseqs) %>%
  mutate(col=row_number()) %>%
  ungroup() %>%
  complete(nesting(sample,time,temp),col,fill=list(pctseqs0=0,pctseqs=0)) %>%
  group_by(sample) %>%
  arrange(-pctseqs0,pctseqs) %>%
  mutate(col=row_number(),
         newotu=pctseqs0==0 & pctseqs>0,
         min.x=min2(col[newotu]),
         max.x=max2(col[newotu]),
         max.y=max2(pctseqs[newotu])) %>%
  ungroup()


# (x=18, y=225, yend=350, l

ggplot(otu) +
  geom_col(aes(x=col,y=pctseqs,fill=otu)) +
  geom_step(aes(x=col,y=pctseqs0),direction="mid") +
  geom_bracket(aes(x=min.x,xend=max.x,y=max.y,label="unique ASVs"),fontsize=3,tip="square") +
  scale_fill_taxonomy(data=otu,fill=otu) +
  scale_y_continuous(trans=log_epsilon_trans(0.001)) +
  facet_grid(temp~time)




# geom_bracket(aes())


library(tidyverse)
mt <- mtcars %>%
  mutate(mpg.group=cut2(mpg, n.splits=3, lvls=c("Low MPG", "Med MPG", "High MPG")))
g <- ggplot(mt) +
  geom_point(aes(x=mpg, y=hp, color=mpg.group), size=3)
g
# whole-data brackets
g + geom_bracket(aes(x=mpg, y=350), label="MPG values") +
  geom_bracket(aes(x=35, y=hp), label="HP values", flip=TRUE)

geom_bracket()













































