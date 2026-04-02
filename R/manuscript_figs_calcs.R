
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


# fig 1A: experiment 1 stackplot -------------------------------------------------------------------

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



# fig 1B: PCOA exp 1 ------------------------------------------------------------


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

# testall1(InvSimpson) # no diff
# testall1(qpcr.totalseqs) # time
# testall1(dist_horn) # time
# testall1(dist_pct.bray) # time and temp

# fig 1c: step compare exp1 ------------------------------------------------------------

samples.compare <- c("1B","1C","1L","1Z")
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
  geom_col(data=otu1compare,aes(x=col,y=pctseqs,fill=otu)) +
  geom_step(data=otu1compare,aes(x=col,y=pctseqs0),direction="mid") +
  geom_text(data=s1compare,aes(x=Inf,y=Inf,label=label),hjust=1,vjust=1,color="blue") +
  geom_rect(data=s1compare,aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,linetype=baseline),
            fill=NA,color="blue",show.legend=FALSE) +
  scale_linetype_manual(values=c("TRUE"="longdash","FALSE"=NA)) +
  geom_bracket(data=filter(otu1compare,extra),
               aes(x=col,y=ave(pctseqs,sample,FUN=max),
                   fontsize=3,label="unique\nASVs"),tip="square") + 
  scale_fill_taxonomy(name="Bacterial Taxa",data=otu1compare,fill=otu) +
  scale_y_continuous("Relative Abundance",trans=log_epsilon_trans(0.001)) +
  facet_wrap(~lbl.compare,nrow=1) +
  theme(aspect.ratio=1,
        axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title = element_blank(), panel.background = element_blank())

g1.asv

# fig 2A: experiment 2 stackplot ------------------------------------------


otu2 <- phy2 %>% 
  # for speed
  # phy.collapse() %>%
  get.otu.melt() %>%
  group_by(treatment,time) %>%
  mutate(sample.number=as.numeric(factor(sample))) %>%
  ungroup() 
s2 <- phy2 %>% get.samp() %>%
  mutate(qpcr.label=ifelse(is.na(qpcr.totalseqs),"undetectable",short_number(qpcr.totalseqs)))


width <- 0.75
g2a <- ggplot() +
  # geom_taxonomy(data=otu2,aes(x=letter,y=pctseqs,fill=otu,label=Species),width=width) +
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
                        label=str_glue("Horn={sprintf('%.3f',dist_horn)}")),hjust="inward",size=3) +
  geom_text(data=s2,aes(x=stage(letter,after_stat=x-xtop),y=0,
                        label=str_glue("qPCR={qpcr.label}")),hjust="inward",size=3)
g2a.alt


# fig 2B: PCOA exp 2 ------------------------------------------------------------


qwer <- function(method,distance,...) {
  ord <- phy.ordinate(phy2,method=method,distance=distance,...)
  opts <- enexprs(...)
  opts_lbl <- paste(names(opts),opts,sep="=",collapse=", ")
  title <- str_glue("{method} {distance} {opts_lbl}")
  ggplot(ord$data,aes(x=axis1,y=axis2)) +
    geom_point(aes(color=treatment),size=3) +
    geom_text(aes(label=treatment),size=2,vjust=1) +
    theme(aspect.ratio=1) +
    ggtitle(title)
}

qwer("PCoA","pct.bray")
qwer("PCoA","horn")






# fig 2C? exp 2 step compare ------------------------------------------------------


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
# testall2(InvSimpson) # no diff
# testall2(qpcr.totalseqs) # time and temp
# testall2(dist_horn) # no diff
# testall2(dist_pct.bray) # time and temp





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



# exp1 tree diffs ---------------------------------------------------------

library(glue)
samples.compare <- c("1A","1B","1C","1L","1Z")
phy1.tree <- phy1 %>%
  # phy.collapse() %>%
  mutate(compare=lbl %in% samples.compare) %>%
  filter(compare|baseline)

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
  mutate(in.baseline=pctseqs0>0,
         pct.diff=pctseqs-pctseqs0,
         pct.ratio=pctseqs/pctseqs0)


# pct.diff
gg$ggtree +
  geom_tile(data=gg$otu,aes(x=x,y=y,fill=pct.diff,color=in.baseline),size=0.5) +
  geom_segment(data=gg$tax,aes(x=x,xend=gg$x.ring.min,y=y,yend=y),color="gray",linetype="dotted") +
  geom_text(data=gg$samp,aes(x=x,y=gg$angle_to_y(90),label=lbl)) +
  scale_color_manual(values=c("TRUE"="pink","FALSE"=NA)) +
  geom_text(data=gg$tax,aes(x=gg$x.ring.max,y=y,label=Species,hjust=hjust,angle=angle),
            color="dark gray",size=2) +
  scale_fill_gradient2(transform=log_epsilon_trans(0.001),limits=c(-1,1))



# pct.ratio
gg$ggtree +
  geom_tile(data=gg$otu,aes(x=x,y=y,fill=pct.ratio,color=in.baseline),size=0.5) +
  geom_segment(data=gg$tax,aes(x=x,xend=gg$x.ring.min,y=y,yend=y),color="gray",linetype="dotted") +
  geom_text(data=gg$samp,aes(x=x,y=gg$angle_to_y(90),label=lbl)) +
  scale_color_manual(values=c("TRUE"="pink","FALSE"=NA)) +
  geom_text(data=gg$tax,aes(x=gg$x.ring.max,y=y,label=Species,hjust=hjust,angle=angle),
            color="dark gray",size=2) +
  scale_fill_gradient2(transform="log")




# exp1 tree diffs ---------------------------------------------------------


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






# exp2 tree diffs? --------------------------------------------------------


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










# permanova ---------------------------------------------------------------

library(vegan)

dist <- calc.distance(phy1,"horn")
s <- get.samp(phy1)
adonis2(dist ~ temp + time, data=s)

