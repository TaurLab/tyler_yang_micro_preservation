
# PRIOR fig 5, exp2 pcoa permanova ----------------------------------------------

phy2b <- phy2
dist2 <- calc.distance(phy2b,"horn")
s2b <- get.samp(phy2b)
# perm2 <- do.permanova(phy2b,dist2,~time+heat+uv)
perm2 <- do.permanova(phy2b,dist2, ~time + heat.75C + uv.sample + uv.dna + heat.autoclave)

# perm2a <- do.permanova(phy2b,dist2, ~time + heat.75C + uv.sample + uv.dna + heat.autoclave)
# perm2b <- do.permanova(phy2b,dist2, ~time + heat + uv)
# perm2c <- do.permanova(phy2b,dist2, ~time + treatment)

perm2$tbl.formatted <- perm2$tbl.formatted %>%
  select(-`beta~'disper'~italic(P)`) %>%
  mutate(Predictor=fct_recode(Predictor,"'time'" = "'time'", 
                              "'75'*degree*'C'" = "'heat.75C'", 
                              "'UV'" = "'uv.sample'", 
                              "'UV DNA'" = "'uv.dna'", 
                              "'autoclave'" = "'heat.autoclave'"))

gtbl2 <- perm2$tbl.formatted %>% ggtexttable(rows=NULL,theme=tt)
size.scaling2 <- s2b$days %>% unique() %>% sort() %>% 
  scales::rescale(to=c(3,5))
g.fig8.pcoa <- perm2$ord$data %>%
  arrange(lbl) %>%
  ggplot(aes(x=axis1,y=axis2)) +
  geom_point(aes(color=treatment,
                 fill=treatment,
                 shape=treatment,
                 size=time),alpha=0.7) + 
  geom_text_repel(aes(label=lbl),size=3,vjust=1.4,max.overlaps = Inf) +
  scale_color_brewer(type="qual",palette=3) +
  scale_fill_brewer(type="qual",palette=3) +
  xlab(perm2$ord.axes[[1]]) + ylab(perm2$ord.axes[[2]]) +
  scale_shape_manual(values=c("none"=21, "75C"=22, "UV"=23, "75C+UV"=24,
                              "autoclave"=23, "autoclave+UV"=24, "UV DNA"=25)) +
  scale_size_manual(values=size.scaling2) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  theme(!!!theme.settings,
        aspect.ratio=1)

pos1 <- c(0.65,0.55)

g.fig8 <- g.fig8.pcoa +
  patchwork::inset_element(gtbl2,pos1[1],pos1[2],pos1[1],pos1[2])
g.fig8



# PRIOR Fig S5: exp2 pcoa permanova using bray ----------------------------------


phy2b <- phy2
dist2bray <- calc.distance(phy2b,"pct.bray")
s2b <- get.samp(phy2b)
# perm2 <- do.permanova(phy2b,dist2,~time+heat+uv)
perm2 <- do.permanova(phy2b,dist2bray, ~time + heat.75C + uv.sample + heat.autoclave + uv.dna)

# perm2a <- do.permanova(phy2b,dist2, ~time + heat.75C + uv.sample + uv.dna + heat.autoclave)
# perm2b <- do.permanova(phy2b,dist2, ~time + heat + uv)
# perm2c <- do.permanova(phy2b,dist2, ~time + treatment)

perm2$tbl.formatted <- perm2$tbl.formatted %>%
  select(-`beta~'disper'~italic(P)`) %>%
  mutate(Predictor=fct_recode(Predictor,"'time'" = "'time'", 
                              "'75'*degree*'C'" = "'heat.75C'", 
                              "'UV'" = "'uv.sample'", 
                              "'UV DNA'" = "'uv.dna'", 
                              "'autoclave'" = "'heat.autoclave'"))

gtbl2 <- perm2$tbl.formatted %>% ggtexttable(rows=NULL,theme=tt)

size.scaling2 <- s2b$days %>% unique() %>% sort() %>% 
  scales::rescale(to=c(3,5))

g.fig.pcoa.bray <- perm2$ord$data %>%
  arrange(lbl) %>%
  ggplot(aes(x=axis1,y=axis2)) +
  geom_point(aes(color=treatment,
                 fill=treatment,
                 shape=treatment,
                 size=time),alpha=0.7) + 
  # geom_text(aes(x=axis1,y=axis2,label=days),size=3,color="#616161",alpha=0.8) +
  geom_text_repel(aes(label=lbl),size=3,vjust=1.4,max.overlaps = Inf) +
  scale_color_brewer(type="qual",palette=3) +
  scale_fill_brewer(type="qual",palette=3) +
  xlab(perm2$ord.axes[[1]]) + ylab(perm2$ord.axes[[2]]) +
  # scale_shape_manual(values=c(1:7)) +
  scale_shape_manual(values=c("none"=21, "75C"=22, "UV"=23, "75C+UV"=24,
                              "autoclave"=23, "autoclave+UV"=24, "UV DNA"=25)) +
  scale_size_manual(values=size.scaling2) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  theme(!!!theme.settings,
        aspect.ratio=1)

g.fig.pcoa.bray

pos1 <- c(0.65,0.55)

g.fig.bray <- g.fig.pcoa.bray +
  patchwork::inset_element(gtbl2,pos1[1],pos1[2],pos1[1],pos1[2])
g.fig.bray


g.fig8 / g.fig.bray


# PRIOR fig S6: exp 2 alpha diversity ---------------------------------------------------------


s2 <- phy2 %>% get.samp(stats=TRUE)

itbl2 <- aov_contrast(InvSimpson ~ time + heat.75C + uv.sample + uv.dna + heat.autoclave, data = s2) %>% 
  mutate(term=fct_recode(term,"time" = "time", 
                         "75'*degree*'C" = "heat.75C", 
                         "UV" = "uv.sample", 
                         "UV DNA" = "uv.dna", 
                         "autoclave" = "heat.autoclave"))
itext2 <- anova_oneline(itbl2)



g2.invsimpson <- ggplot(s2,aes(x=lbl,y=InvSimpson)) +
  geom_col(fill="steelblue") +
  exp2.facet + 
  theme(!!!theme.settings,
        panel.spacing.x = exp2.panel.spacing.x) +
  ggplot2::labs(title="Inverse Simpson index",
                x="Sample", 
                y="Inverse Simpson index",
                caption=parse(text=itext2))

g2.invsimpson

# stbl2 <- aov_contrast(Shannon ~ time + heat.75C + uv.sample + uv.dna + heat.autoclave, data = s2) %>% 
#   mutate(term=fct_recode(term,"time" = "time", 
#                          "75'*degree*'C" = "heat.75C", 
#                          "UV" = "uv.sample", 
#                          "UV DNA" = "uv.dna", 
#                          "autoclave" = "heat.autoclave"))
# stext2 <- anova_oneline(stbl2)
# 
# g2.shannon <- ggplot(s2,aes(x=lbl,y=Shannon)) +
#   expand_limits(y=5.5) +
#   geom_col(fill="purple") +
#   exp2.facet + 
#   theme(!!!theme.settings,
#         panel.spacing.x = exp2.panel.spacing.x) +
#   ggplot2::labs(title="Shannon index",
#                 x="Sample", 
#                 y="Shannon index",
#                 caption=parse(text=stext2))
# 
# g2.alpha.diversity <- g2.invsimpson / g2.shannon
# g2.alpha.diversity






# PRIOR fig S7: exp2 qpcr ----------------------------------------------------------------


s2 <- get.samp(phy2) %>%
  mutate(qpcr.lbl=ifelse(is.na(qpcr.totalseqs),"(undetectable)",NA_character_),
         qpcr.totalseqs.impute=coalesce(qpcr.totalseqs,1000),
         log.qpcr.totalseqs.impute=log(qpcr.totalseqs.impute))

aov2.qpcr <- aov_contrast(log.qpcr.totalseqs.impute ~ time + heat.75C + uv.sample + uv.dna + heat.autoclave, data=s2) %>%
  mutate(term=fct_recode(term,"time" = "time", 
                         "75'*degree*'C" = "heat.75C", 
                         "UV" = "uv.sample", 
                         "UV DNA" = "uv.dna", 
                         "autoclave" = "heat.autoclave"))

qtext <- anova_oneline(aov2.qpcr)

g2.qpcr <- ggplot(s2) + 
  geom_col(aes(x=lbl,y=qpcr.totalseqs),fill="steelblue",width=width) + 
  geom_text(aes(x=lbl,y=0,label=qpcr.lbl),hjust=0,angle=90) +
  scale_x_discrete("Sample",expand=expansion(add=0.5)) +
  scale_y_continuous("Total abundance by 16S qPCR",trans=log_epsilon_trans(1000),
                     # expand=FALSE,
                     expand=expansion(mult=0.025),
                     labels=pretty_power10) +
  exp2.facet +
  theme(panel.spacing.x=exp2.panel.spacing.x) +
  labs(caption=parse(text=qtext))


g2.qpcr








# (OLD) fig 9: exp 2 pcoa and permanova -----------------------------------------

phy2b <- phy2
dist2 <- calc.distance(phy2b,"horn")
s2b <- get.samp(phy2b)
perm2 <- do.permanova(phy2b,dist2,~time+heat+uv)


perm2x <- do.permanova(phy2b,dist2,~time+heat*uv)
perm2x <- do.permanova(phy2b,dist2, ~time + uv.sample + uv.dna + heat.75C + heat.autoclave)

# 75C_vs_autoclave


perm2$heat.pair.formatted <- perm2$tbl$contrasts$heat %>%
  adjust.pairs(c("75C vs no heat" = "no heat_vs_75C", 
                 # "no heat_vs_autoclave" = "no heat_vs_autoclave", 
                 "autoclave vs 75C" = "75C_vs_autoclave"))

perm2$uv.pair.formatted <- perm2$tbl$contrasts$uv %>%
  adjust.pairs(c("no UV vs UV" = "UV_vs_no UV", 
                 # "UV DNA_vs_no UV" = "UV DNA_vs_no UV",
                 "UV vs UV DNA" = "UV_vs_UV DNA"))

perm2$tbl.formatted <- perm2$tbl.formatted %>%
  mutate(Predictor=str_replace(Predictor,"uv","UV"))

perm2$heat.pair.formatted <- perm2$heat.pair.formatted %>%
  mutate(Pair=str_replace(Pair,"75C","75'*degree*'C"))

gtbl2a <- perm2$tbl.formatted %>% ggtexttable(rows=NULL,theme=tt)
gtbl2b <- perm2$heat.pair.formatted %>% ggtexttable(rows=NULL,theme=tt)
gtbl2c <- perm2$uv.pair.formatted %>% ggtexttable(rows=NULL,theme=tt)

g.fig8.pcoa <- perm2$ord$data %>%
  arrange(lbl) %>%
  ggplot(aes(x=axis1,y=axis2)) +
  geom_point(aes(color=treatment),size=4,alpha=0.8) + 
  geom_text_repel(aes(label=lbl),size=3,vjust=1.4,max.overlaps = Inf) +
  scale_color_brewer(type="qual",palette=3) +
  xlab(perm2$ord.axes[[1]]) + ylab(perm2$ord.axes[[2]]) +
  theme(aspect.ratio=1)

pos1 <- c(0.75,0.55)
pos2 <- c(0.75,0.4)
pos3 <- c(0.75,0.25)

g.fig8 <- g.fig8.pcoa +
  patchwork::inset_element(gtbl2a,pos1[1],pos1[2],pos1[1],pos1[2]) +
  patchwork::inset_element(gtbl2b,pos2[1],pos2[2],pos2[1],pos2[2]) +
  patchwork::inset_element(gtbl2c,pos3[1],pos3[2],pos3[1],pos3[2])
g.fig8





# OLD fig 1: exp1 taxplot ---------------------------------------------------------------


otu1 <- phy1 %>% get.otu.melt()
s1 <- phy1 %>% get.samp()

g1.tax <- ggplot(otu1) +
  geom_taxonomy(aes(x=lbl,y=pctseqs,fill=otu,label=Species),width=width,label.split=TRUE) +
  geom_col(data=filter(s1,lbl=="1A"),aes(x=lbl,y=1,color="baseline sample (1A)"),
           linewidth=1,linetype="longdash",fill=NA,width=width) +
  scale_y_continuous("Relative abundance",
                     expand=FALSE,labels=scales::label_percent()) +
  scale_x_discrete(name="Sample",expand=expansion(add = 0.5)) + 
  scale_color_manual(name="Sample Comparison",values=c("baseline sample (1A)"="red")) +
  scale_fill_taxonomy(name="Bacterial Taxa",data=otu1,tax.palette=pal,fill=otu) +
  xlab("Sample") +
  exp1.facet +
  theme(legend.key.size = unit(0.85,"lines"),
        strip.clip="on",
        panel.spacing.x = exp1.panel.spacing.x)

g1.tax


# OLD fig 7: exp2 stackplot ---------------------------------------------------

s2 <- get.samp(phy2)
otu2 <- get.otu.melt(phy2)

g2.tax <- ggplot(otu2) +
  geom_taxonomy(aes(x=lbl,y=pctseqs,label=Species,fill=otu),
                tax.palette=pal,label.split=TRUE,fontsize=3,width=width) +
  geom_col(data=filter(s2,baseline),
           aes(x=lbl,y=1,color="baseline sample (2A)"),fill=NA,width=width,
           linetype="longdash",linewidth=0.75) +
  scale_x_discrete("Sample",expand=expansion(add=0.5)) +
  scale_y_continuous("Bacterial Composition\nRelative abundance",expand=FALSE) +
  scale_color_manual("Sample Comparison",values=c("baseline sample (2A)"="red")) +
  scale_fill_taxonomy(name="Bacterial taxa",data=otu2,fill=otu,tax.palette=pal) +
  exp2.facet +
  theme(legend.key.size = unit(0.85,"lines"),
        panel.spacing.x = exp2.panel.spacing.x,
        axis.text.y=element_blank(),
        axis.ticks=element_blank())
g2.tax


# OLD fig s3: distance 2 --------------------------------------------------------------

phy1a <- phy1 %>%
  mutate(x=as.numeric(lbl),
         baseline=lbl=="1A",
         dist=dist_mean.horn)

phy2b <- phy2 %>%
  mutate(x=as.numeric(lbl),
         baseline=lbl=="1A",
         dist=dist_2A)
dist2 <- calc.distance(phy2b,"horn")
s2b <- get.samp(phy2b) %>%
  mutate(dist=dist_2A)

dtext <- s2b %>% filter(sample=="2A") %>%
  mutate(label="baseline\nsample",
         x=1.5,y=0.02,
         xend=1,yend=0.003)
dtext.layer <- list(geom_text(data=dtext,aes(x=x,y=y,label=label),vjust=-0.02,lineheight=0.8),
                    geom_segment(data=dtext,aes(x=x,y=y,xend=xend,yend=yend)))
# 
g2.dist <- ggplot(s2b) +
  geom_col(aes(x=lbl,y=dist),fill="steelblue",width=width) +
  geom_text(aes(x=lbl,y=dist,label=pretty_number(dist),color=baseline),
            vjust=0,size=3,show.legend = FALSE) +
  scale_y_continuous(name="Distance\n(compared with baseline)")  +
  scale_color_manual(values=c("TRUE"="red","FALSE"="black")) +
  dtext.layer +
  xlab("Sample") +
  exp2.facet +
  theme(panel.spacing.x = exp2.panel.spacing.x)


g2.dist







# OLD fig 2: distance ---------------------------------------------------------

phy1a <- phy1 %>%
  mutate(x=as.numeric(lbl),
         # dist=dist_mean.horn,
         dist=dist_horn,
         baseline=lbl=="1A")

# otu1 <- phy1a %>% get.otu.melt()
s1a <- phy1a %>% get.samp(stats=TRUE)
dtext <- s1a %>% filter(sample=="1A") %>%
  mutate(label="baseline\nsample",
         x=1.5,y=0.02,
         xend=1,yend=0.003)
dtext.layer <- list(geom_text(data=dtext,aes(x=x,y=y,label=label),vjust=-0.02,lineheight=0.8),
                    geom_segment(data=dtext,aes(x=x,y=y,xend=xend,yend=yend)))

g1.dist <- ggplot(s1a) +
  geom_col(aes(x=lbl,y=dist),fill="steelblue",width=width) +
  geom_text(aes(x=lbl,y=dist,label=pretty_number(dist),color=baseline),
            vjust=0,size=3,show.legend = FALSE) +
  dtext.layer +
  scale_y_continuous(name="Distance\n(compared with baseline)")  +
  scale_color_manual(values=c("TRUE"="red","FALSE"="black")) +
  xlab("Sample") +
  exp1.facet +
  theme(panel.spacing.x = exp1.panel.spacing.x)

g1.dist


# OLD: fig 9: exp2 contam ------------------------------------------------------

phy2ax <- phy2a %>%
  summarize_sample_data(pct.notcontam=mean(as.logical(not.contam)),
                        pct.contam=mean(!as.logical(not.contam)),
                        n.asv.contam=sum(!as.logical(not.contam)),
                        pctseqs.contam=sum(pctseqs[!as.logical(not.contam)]),
                        pctseqs.noncontam=sum(pctseqs[as.logical(not.contam)]),
                        filter.zero=TRUE) %>%
  add_dist("horn",comparator="TY.1_D0_NT",varname=dist_2A) %>%
  add_dist("mean.horn",comparator="TY.1_D0_NT",varname=dist_2A_meanhorn) %>%
  add_dist("horn",comparator="PCRNeg4",varname=dist_PCRNeg)

otu2ax <- get.otu.melt(phy2ax)
s2ax <- get.samp(phy2ax)



exp2ax.panel.spacing.x <- s2ax %>% distinct(treatment,time) %>%
  arrange(treatment,time) %>% 
  group_by(treatment) %>%
  transmute(x=ifelse(row_number()==n(),3.5,1.5)) %>%
  ungroup() %>% slice(-n()) %>% pull(x) %>% unit("points")
exp2.facet <- facet_nested(. ~ treatment+time,scales="free_x",space="free_x",
                           nest_line=element_line(),
                           resect=unit(3,"pt"),
                           solo_line=TRUE)

g.contam <- ggplot(s2ax) + 
  geom_col(aes(x=lbl,y=pctseqs.contam),fill="dark gray") +
  scale_x_discrete("Sample",expand=expansion(add=0.5)) +
  scale_y_continuous("Proportion of Percent sequences contaminant",
                     expand=expansion(mult=0.025),
                     labels=scales::label_percent(1)) +
  exp2.facet +
  theme(panel.spacing.x=exp2ax.panel.spacing.x)

g.contam




# OLD ---------------------------------------------------------------------


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

g.qpcr
g.asv.compare
g.tax.compare

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











# testing exp 1 (old) ----------------------------------------------------------

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
# testing exp 1 (old2) ----------------------------------------------------------

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


# base <- phy2a %>% filter(sample=="TY.1_D0_NT") %>% get.otu.melt(sample_data=FALSE) %>% select(otu,taxid)
# ctrl <- phy2a %>% filter(lbl=="PCRNeg") %>% get.otu.melt(sample_data=FALSE) %>% select(otu,taxid)

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

















# tyler clarifying exp 2  -------------------------------------------------



phy2b <- phy2
dist2 <- calc.distance(phy2b,"horn")
s2b <- get.samp(phy2b)

# method 1: separate terms (heat.75 is signif which is not intuitive)
adonis2(dist2 ~ time + heat.75C + uv.sample + uv.dna + heat.autoclave, data=s2b)

# method 2: one big term (doesn't work, nothing is signif in pairwise)
adonis2(dist2 ~ time + treatment, data=s2b)

pairwise.adonis(dist2,s2b$treatment)

# method 3: two grouping terms (can switch to this, possibly)
adonis2(dist2 ~ time + heat + uv, data=s2b)
pairwise.adonis(dist2,s2b$heat)
pairwise.adonis(dist2,s2b$uv)

pp <- phy2 %>% filter(treatment %in% c("none","autoclave"))
ss <- pp %>% get.samp()
dd <- calc.distance(pp,"horn")
adonis2(dd ~ time + treatment, data=ss)

s2b
















