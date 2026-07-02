library(yingtools2)
library(tidyverse)
library(ytrecipes)
library(rlang)
library(yingtools2)
library(vegan)
library(pairwiseAdonis)
library(phyloseq)
d <- tribble(
  ~treatment,     ~uv.dna, ~uv.sample, ~uv,      ~heat,       ~time,   ~heat.75C, ~heat.autoclave,
  "none",         FALSE, FALSE,    "no UV",  "no heat",   "Day 0", FALSE,   FALSE,
  "none",         FALSE, FALSE,    "no UV",  "no heat",   "Day 6", FALSE,   FALSE,
  "none",         FALSE, FALSE,    "no UV",  "no heat",   "Day 9", FALSE,   FALSE,
  "75C",          FALSE, FALSE,    "no UV",  "75C",       "Day 0", TRUE,    FALSE,
  "75C",          FALSE, FALSE,    "no UV",  "75C",       "Day 6", TRUE,    FALSE,
  "75C",          FALSE, FALSE,    "no UV",  "75C",       "Day 9", TRUE,    FALSE,
  "UV",           FALSE, TRUE,     "UV",     "no heat",   "Day 0", FALSE,   FALSE,
  "UV",           FALSE, TRUE,     "UV",     "no heat",   "Day 6", FALSE,   FALSE,
  "UV",           FALSE, TRUE,     "UV",     "no heat",   "Day 9", FALSE,   FALSE,
  "75C+UV",       FALSE, TRUE,     "UV",     "75C",       "Day 0", TRUE,    FALSE,
  "75C+UV",       FALSE, TRUE,     "UV",     "75C",       "Day 6", TRUE,    FALSE,
  "75C+UV",       FALSE, TRUE,     "UV",     "75C",       "Day 9", TRUE,    FALSE,
  "autoclave",    FALSE, FALSE,    "no UV",  "autoclave", "Day 0", FALSE,   TRUE,
  "autoclave",    FALSE, FALSE,    "no UV",  "autoclave", "Day 6", FALSE,   TRUE,
  "autoclave",    FALSE, FALSE,    "no UV",  "autoclave", "Day 9", FALSE,   TRUE,
  "autoclave+UV", FALSE, TRUE,     "UV",     "autoclave", "Day 0", FALSE,   TRUE,
  "autoclave+UV", FALSE, TRUE,     "UV",     "autoclave", "Day 6", FALSE,   TRUE,
  "autoclave+UV", FALSE, TRUE,     "UV",     "autoclave", "Day 9", FALSE,   TRUE,
  "UV DNA",       TRUE,  FALSE,    "UV DNA", "no heat",   "Day 0", FALSE,   FALSE,
  "UV DNA",       TRUE,  FALSE,    "UV DNA", "no heat",   "Day 6", FALSE,   FALSE,
  "UV DNA",       TRUE,  FALSE,    "UV DNA", "no heat",   "Day 9", FALSE,   FALSE
)

b.uv.sample <- 0
b.heat.75C <- 0
b.uv.dna <- c(10,0,0)
b.heat.autoclave <- c(0,9,0)


b.uv.sample <- b.uv.sample %>% rep_len(3)
b.heat.75C <- b.heat.75C %>% rep_len(3)
b.uv.dna <- b.uv.dna %>% rep_len(3)
b.heat.autoclave <- b.heat.autoclave %>% rep_len(3)


dd <- d %>%
  mutate(x=1 + 
           b.uv.sample[1]*uv.sample + 
           b.uv.dna[1]*uv.dna + 
           b.heat.autoclave[1]*heat.autoclave + 
           b.heat.75C[1]*heat.75C,
         y=1 + 
           b.uv.sample[2]*uv.sample + 
           b.uv.dna[2]*uv.dna + 
           b.heat.autoclave[2]*heat.autoclave + 
           b.heat.75C[2]*heat.75C,
         z=1 + 
           b.uv.sample[3]*uv.sample + 
           b.uv.dna[3]*uv.dna + 
           b.heat.autoclave[3]*heat.autoclave + 
           b.heat.75C[3]*heat.75C)

ddd <- rep(list(dd),10) %>% list_rbind() %>% 
  mutate(x=x + rnorm(n=n(),mean=0,sd=0.1),
         y=y + rnorm(n=n(),mean=0,sd=0.1),
         z=z + rnorm(n=n(),mean=0,sd=0.1),
         sample=paste0("x",1:n()))

ggplot(ddd,aes(x=x,y=y,color=treatment)) + geom_point()

atbl <- aov_contrast(x ~ time + uv.dna + uv.sample + heat.autoclave + heat.75C,data=ddd)
atbl

##### phy #####

o <- ddd %>% select(sample,x,y,z) %>% 
  column_to_rownames("sample") %>% 
  t() %>% set.otu()
s <- ddd %>% select(sample,treatment,uv.dna,uv.sample,uv,heat,time,heat.75C,heat.autoclave) %>%
  set.samp()
pp <- phyloseq(o,s)
ppdist <- calc.distance(pp,"horn")
ord <- phy.ordinate(pp,method="PCoA",distance=ppdist)
ggplot(ord$data,aes(x=axis1,y=axis2,color=treatment)) +
  geom_point()

perm1 <- do.permanova(pp,ppdist,
                     ~ time + uv.dna + uv.sample + heat.autoclave + heat.75C,
                     permutations = 999)
perm1$tbl
perm1$tbl.formatted

perm2 <- do.permanova(pp,ppdist,
                     ~ time + treatment,
                     permutations = 999)
perm2$tbl
perm2$tbl.formatted
perm2$tbl$contrasts


perm3 <- do.permanova(pp,ppdist,
                     ~ time + uv + heat,
                     permutations = 999)
perm3$tbl.formatted
perm3$tbl



