

library(yingtools2)
library(tidyverse)
library(phyloseq)
rm(list=ls())

# copied from iski0021
load("data-raw/dada2_full_pipeline.RData")
# usethis::use_git_ignore(c("data-raw/*.RData","data-raw/*.rds"))


s <- get.samp(phy.dada2) %>%
  filter(!grepl("^(AT|ANG|PT|DT)",sample),
         !grepl("^(EB4|PCRNeg4|Zymo4)",sample)) %>%
  mutate(experiment=ifelse(grepl("pool1161",sample),2,1),
         time=as.numeric(str_extract(sample,middle.pattern("D","[0-9]+"))),
         time=coalesce(time,0),
         time=paste("day",time),
         temp=recode.grep(sample,c("RT|NT|Dry"="room temp","4C"="4C","-80"="-80C","-20"="-20C"),else.value="room temp"),
         storage=ifelse(experiment==1,"in cryovial","open"),
         uv=recode.grep(sample,c("DNA_UV"="UV DNA","UV"="UV"),else.value="no UV"),
         heat=case_when(
           experiment==1 ~ "75C",
           TRUE ~ recode.grep(sample,c("NT|Dry"="no heat","75C"="75C","AC"="autoclave"),else.value="no heat")
         ),
         treatment=paste2(na_if(heat,"no heat"),na_if(uv,"no UV"),sep="+"),
         treatment=ifelse(is.na(treatment),"none",treatment),
         treatment=factor(treatment,levels=c("none", "75C", "autoclave", "UV", "UV DNA", "75C+UV", "autoclave+UV")))

phy.tyler <- prune_samples(s$sample,phy.dada2)
sample_data(phy.tyler) <- s %>% set.samp()
sample_names(phy.tyler) <- sub("\\.\\.pool1161","",sample_names(phy.tyler))


qpcr <- read_csv("data-raw/qpcrdata.csv")
sample_data(phy.tyler) <- phy.tyler %>% get.samp() %>% left_join(qpcr,by="sample") %>% set.samp()


save(phy.tyler,file="data/phy.tyler.RData")


trace <- readRDS("data-raw/trace_list_tyler.rds")
tax.blast.tyler <- tax.blast.full
# setequal(tax.blast.tyler$otu,taxa_names(phy.tyler))
trace_tbl.tyler <- trace_tbl %>% mutate(sample=sub("\\.\\.pool1161","",sample)) %>%
  filter(sample %in% sample_names(phy.tyler))
# setequal(sample_names(phy.tyler),trace_tbl.tyler$sample)



save(tax.blast.tyler,trace_tbl.tyler,file="data/phy.tyler.additional.RData")



# trace -------------------------------------------------------------------



t <- trace_tbl %>%
  mutate(asv.ratio=nochim.asvs/seqtab.asvs,
         seq.ratio=nochim/seqtab)



ggplot(t,aes(x=sample,y=asv.ratio)) + geom_col() +
  theme(axis.text.x=element_text(angle=-90)) +
  expand_limits(y=c(0,1))

ggplot(t,aes(x=sample,y=seq.ratio)) + geom_col() +
  theme(axis.text.x=element_text(angle=-90)) +
  expand_limits(y=c(0,1))





# examine trace file ------------------------------------------------------

library(yingtools2)
library(tidyverse)
library(phyloseq)

# trace <- read_rds("C:/Users/Ying/Desktop/tyler/trace_list_tyler.rds")
# tr <- trace[c("seqtab","seqtab_nochim")]
# saveRDS(tr,"tyler_trace.rds")
setwd("C:/Users/Ying/Desktop/tyler")
rm(list=ls())
trace <- read_rds("tyler_trace.rds")
load("phy.tyler.RData")

s1 <- trace$seqtab
s2 <- trace$seqtab_nochim

asvs.kept <- colnames(s2)
asvs.all <- colnames(s1)
asvs.reject <- setdiff(asvs.all,asvs.kept)

otu <- s1 %>% as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols=-sample,names_to="asv",values_to="nseqs") %>%
  filter(nseqs>0) %>%
  mutate(kept=asv %in% colnames(s2)) %>%
  group_by(sample) %>%
  mutate(pctseqs=nseqs/sum(nseqs),
         sample=sub("..pool1161","",sample,fixed=TRUE)) %>%
  ungroup()

sum.samps <- function(otu) {
  samps <- otu %>% 
    group_by(sample) %>%
    summarize(totalseqs=sum(nseqs),
              n.asv=sum(nseqs>0),
              n.asv.singleton=sum(nseqs==1),
              pct.asv.singleton=n.asv.singleton/n.asv,
              pct.nseqs.singleton=n.asv.singleton/totalseqs,
              q50.pctseqs=median(pctseqs),
              q25.pctseqs=quantile(pctseqs,probs=0.25),
              NULL) %>%
    ungroup()
  samps
}

samps.all <- otu %>% sum.samps() %>% mutate(source="pre")
samps.nochim <- otu %>% filter(kept) %>% sum.samps() %>% mutate(source="post")
s.tyler <- get.samp(phy.tyler)

samps.long <- bind_rows(samps.all,samps.nochim) %>%
  pivot_longer(cols=-c(sample,source),names_to="metric") %>%
  inner_join(s.tyler,by="sample")

samps <- samps.long %>%
  pivot_wider(id_cols=c(sample,metric,!!!syms(names(s.tyler))),
              names_from=source,
              values_from=value) %>%
  arrange(treatment) %>% 
  mutate(sample2=paste(coalesce(treatment,sample),dense_rank(sample),sep="|"),
         sample2=fct_inorder(sample2),
         delta=(post-pre)/pre)

samps %>% filter(metric %like% "n.asv") %>% 
  group_by(metric) %>% arrange(treatment) %>% dt

ggplot(samps,aes(x=sample2,y=delta,label=scales::percent(delta,accuracy=1,big.mark=""),fill=treatment)) +
  geom_col() + 
  geom_text(vjust=0,size=3) +
  facet_grid(metric~.,scales="free_y") +
  theme(axis.text.x = element_text(angle=-90))








