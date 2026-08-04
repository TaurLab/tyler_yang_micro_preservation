

# 

library(yingtools2)
library(tidyverse)
library(ytpipeline)

rm(list=ls())

load("data/phy.tyler.RData")
load("data/phy.tyler.additional.RData")
# phy.others = read_rds("data/other.samps.rds")

t.orig <- phy.tyler %>% get.tax()
branks1 <- c("Superkingdom","Phylum","Class","Order","Family","Genus","Species")
branks2 <- c("Domain","Phylum","Class","Order","Family","Genus","Species")

# old tax, with taxid
t1 <- tax.blast.tyler %>%
  group_by(otu) %>%
  arrange(evalue,as.numeric(staxid)) %>%
  dplyr::slice(1) %>% 
  ungroup() %>%
  select(otu,!!!syms(branks1),taxid=staxid) %>%
  arrange(otu)

# new tax, with taxid
t2 <- t1 %>%
  select(otu,taxid) %>%
  mutate(full_taxonomy=get_full_taxonomy(taxid)) %>%
  get_ranks_from_full_taxonomy(taxranks=branks2,keep_taxonomy_var = FALSE) %>%
  mutate(across(.cols=branks2, .fns=~coalesce(.x,"xxxx"))) %>%
  arrange(otu) 

# compare and show changes
x1 <- t1 %>% pivot_longer(-c(otu,taxid))
x2 <- t2 %>% rename(Superkingdom=Domain) %>% pivot_longer(-c(otu,taxid))

x.compare <- inner_join(x1,x2,by=c("otu","taxid","name")) %>%
  group_by(name,value.x,value.y) %>%
  summarize(n=n(),
            taxid=list(unique(taxid)),
            .groups="drop") %>%
  mutate(changed=value.x!=value.y) %>%
  group_by(name,value.x) %>%
  mutate(total=sum(n),
         pct=n/total) %>%
  ungroup() %>%
  mutate(status=case_when(
    !changed ~ "unchanged",
    changed & (n==total) ~ "changed",
    changed & (n!=total) ~ "changed for some",
    TRUE ~ "???"
  )) %>%
  select(name,value.x,value.y,n,changed,total,pct,status,everything()) %>%
  mutate(name=factor(name,levels=branks1))

# show changes
x.changed1 <- x.compare %>% filter(changed) %>% 
  arrange(factor(name,levels=branks1),desc(n)) %>%
  mutate(pct=scales::label_percent(0.1)(pct))
x.changed1 %>% group_by(name) %>% dt()

save(t1,t2,x.compare,file="data/new_taxonomy.RData")







