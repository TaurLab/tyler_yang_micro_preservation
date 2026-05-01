

# load packages -----------------------------------------------------------


library(yingtools2)
library(tidyverse)
library(phyloseq)
library(rentrez)
library(XML)
rm(list=ls())




# pull all biosample ids from rentrez -------------------------------------

if (FALSE) {
  library(yingtools2)
  library(tidyverse)
  library(phyloseq)
  library(rentrez)
  library(XML)
  rm(list=ls())
  
  load("data/phy.tyler.RData")
  load("data/phy.tyler.additional.RData")
  phy.pcrneg.control <- read_rds("data/phy.pcrneg.control.rds")
  
  # all otus, from phy.tyler and phy.pcrneg.control
  otus <- c(taxa_names(phy.tyler),taxa_names(phy.pcrneg.control))
  # find all taxids from tax.blast.tyler
  tax.taxid <- tax.blast.tyler %>%
    filter(otu %in% otus) %>%
    group_by(otu) %>%
    arrange(evalue.rank,staxid) %>%
    slice(1) %>%
    ungroup() %>%
    distinct(Species,taxid=staxid) %>%
    mutate(taxid=as.character(taxid))
  tax.taxid
  # pulls sample IDs from entrez, given taxid
  get_sample_ids <- function(taxid) {
    tryCatch({
      link <- entrez_link(dbfrom="taxonomy", id=taxid, db="biosample")
      ids <- link$links$taxonomy_biosample
      return(ids)
    }, error=function(e) {
      return(character())
    })
  }

  # create df with taxid and list-col of sample_ids
  bact.attributes <- tax.taxid %>% 
    mutate(sample_ids=map(taxid,get_sample_ids,.progress=TRUE),
           n_samples=map_int(sample_ids,length))

  saveRDS(bact.attributes,file="data/bact.attributes.rds")  
  bact.attributes <- readRDS("data/bact.attributes.rds")

  bact.attributes %>% dt
  bact.attributes %>% arrange(desc(n_samples))
}

bact.attributes <- readRDS("data/bact.attributes.rds")
# head(bact.attributes)



# download xml data from some biosamples ----------------------------------


if (FALSE) {
  # create attr, an empty df to hold xml data
  library(yingtools2)
  library(tidyverse)
  library(phyloseq)
  library(rentrez)
  library(XML)
  rm(list=ls())
  
  bact.attributes <- readRDS("data/bact.attributes.rds")
  attr <- bact.attributes %>%
    unnest(sample_ids) %>%
    mutate(xml=NA_character_)
  saveRDS(attr,file="data/attributes.rds")  
}

# given subset of attr, fetch xml data and populate the xml column
fetch_samples <- function(subattr) {
  if (!all(is.na(subattr$xml))) {
    cli::cli_abort("YTError: not all xml empty")
  }
  new_subattr <- subattr %>%
    mutate(xml = map_chr(sample_ids,~{
      tryCatch({
        entrez_fetch(db="biosample",id=.x,rettype="xml")  
      },error=function(e) {
        cli::cli_warn("YTWarning: failed to download for sample_id={sample_id}")
        "error"
      })
    },.progress=TRUE))
}
# after using fetch_samples(subattr), 
# run this to update xml column in attr with new xml values.
update_attr <- function(attr,newsub) {
  nstart <- sum(!is.na(attr$xml))
  newattr <- attr %>% left_join_replace(newsub,by=c("Species","taxid","sample_ids"))
  nstop <- sum(!is.na(newattr$xml))
  cli::cli_alert_info("attr updated: {nstart} -> {nstop}")
  return(newattr)
}
# reports status of pull
status <- function() {
  sum <- attr %>% mutate(status=case_when(
    is.na(attr$xml) ~ "not done",
    !is.na(attr$xml) & attr$xml=="error" ~ "error",
    TRUE ~ "xml done"
  )) %>% group_by(Species,taxid) %>%
    summarize(n.samples=n(),
              n.done=sum(status=="xml done"),
              n.error=sum(status=="error"),
              .groups="drop") %>%
    mutate(pct.done=n.done/n.samples,
           complete=n.done==n.samples) %>%
    arrange(desc(n.samples))
  all <- sum %>%
    summarize(
      n.taxids=n(),
      max.done=max(n.done),
      total.samples.done=sum(n.done),
      total.samples=sum(n.samples),
      median.done=median(n.done),
      min.samples=min(n.samples),
      max.samples=max(n.samples),
      n.complete=sum(complete)
    )
  cli::cli_alert_info("n.taxids={all$n.taxids}, 
  total.samples.done={pretty_number(all$total.samples.done)}
  total.samples={pretty_number(all$total.samples)}
  n.taxids.complete={all$n.complete}
  max.samples={all$max.samples},
  max.samples.done={all$max.done}")
}

attr <- readRDS("data/attributes.rds")
status()

# **** run up to here to start pulling ****




# **** run this below to get more data and save ****
if (FALSE) {
  # find biosamples to pull. 
  sub <- attr %>% filter(is.na(xml)) %>%
    group_by(taxid) %>%
    slice(1:2) %>%
    ungroup()
  # fetch the biosample data
  newsub <- sub %>% fetch_samples()
  
  # update attr with data
  oldattr <- attr
  attr <- attr %>% update_attr(newsub)

  # save data
  saveRDS(attr,file="data/attributes.rds")  

}

# develop habitat extracting code -----------------------------------------


# attr <- readRDS("data/attributes.rds")

# # df bact: remove blanks and parse xml
# bact <- attr %>% filter(!is.na(xml)) %>%
#   mutate(xml=map(xml,xmlInternalTreeParse))
# 
# # all extract all attributes from xml data
# # cols: Species,taxid,sample_ids,n_samples,xml,attributes
# all_attribute_combos <- bact %>%
#   mutate(attributes=map(xml,~{
#   xpathApply(.x,"/BioSampleSet/BioSample/Attributes/Attribute",
#              function(x) {
#                setNames(xmlValue(x),xmlGetAttr(x,"attribute_name"))
#              }) %>% list_c()
# })) 
# 
# all_attribute_table <- all_attribute_combos %>%
#   select(-n_samples,-xml) %>%
#   mutate(attributes=map(attributes,~{
#     tibble(attribute_name=names(.x),value=unname(.x))
#   })) %>% unnest(attributes)



# atable: pull xml from attr, parse xml text into xml object, 
# pull *all* Attributes from each biosample, unnest
atable <- attr %>% 
  filter(!is.na(xml)) %>%
  mutate(attributes=map(xml,~{
    xml <- xmlInternalTreeParse(.x)
    attrs <- xpathApply(xml,"/BioSampleSet/BioSample/Attributes/Attribute",
                        function(x) {
                          setNames(xmlValue(x),xmlGetAttr(x,"attribute_name"))
                        }) %>% list_c()
    tibble(attribute_name=names(attrs),value=unname(attrs))
  },.progress="parsing xml and pulling attributes")) %>% 
  select(-n_samples,-xml) %>% 
  unnest(attributes,keep_empty = TRUE)

# tally of attribute_names
# (NA if biosample does not have attributes)
atable %>% count(attribute_name,sort=TRUE)
# samples with no attributes:
atable %>% filter(is.na(attribute_name))


# peruse attribute_name / value combos
# atable %>% 
#   filter(attribute_name %ilike% "source") %>%
#   count(attribute_name,value) %>%
#   group_by(attribute_name) %>%
#   arrange(desc(n)) %>%
#   slice(1:10) %>% dt()

# patterns for attribute_name, corresponding to habitat var
attr_patterns <- list(
  isolation_source="isolat.*source",
  env_local_scale="env_local_scale",
  env_broad_scale="env_broad_scale",
  env_medium="env_medium",
  source_material_id="source_material_id",
  sample_type="sample_type"
)
# patterns for value
patterns <- list(
  gut="stool|gut|fa?eces|fa?ecal|ca?ecum|ca?ecal|intestin|rectal|defecate|digestive tract",
  oral="sputum|nasopharynx|oral|saliva|tooth|dental|nasal",
  skin="skin",
  genitourinary="urine|urinary|urethra|vagina",
  other.body="ocular|blood|respiratory",
  animal="salmon|honeybee|ixodes",
  plant="endosphere|rhizosphere|maize|conifer",
  other.organism="mushroom",
  food="food|dairy|milk|cheese",
  equipment="equipment|hospital|microplastic",
  urban="groundwater|urban|plumbing|landfill",
  earth="soil|cave|air|seawater|river", 
  missing="missing|not provided|not applicable|^$"
)

# add search to atable
asearch <- atable %>%
  replace_grep_data(value,patterns,hits=value.hits,collapse.fn=~paste(unique(names(.x)),collapse="|")) %>%
  mutate(hit=map_chr(str_split(value.hits,"\\|"),first),
         hit=factor(hit,levels=names(patterns)),
         attr=recode.grep(attribute_name,attr_patterns,else.value="other",as.factor=TRUE),
         n.hits=str_count(value.hits,"\\|")+1,
         n.hits=coalesce(n.hits,0),
         is.missing=hit %in% "missing",
         is.habitat=attr!="other",
         usable=is.habitat & !is.missing)

# collapse sample level
asample <- asearch %>%
  group_by(Species,taxid,sample_ids) %>%
  arrange(attr) %>%
  summarize(s.hit=hit[usable][1],
            # same text met multiple patterns
            s.hit.detail=value.hits[usable][1],
            # different attrs in biosample showed different results
            s.hit.detail2=paste2(unique(hit[usable]),collapse="|"),
            n.usable=sum(usable),
            n.attrs=sum(!is.na(attribute_name)),
            .groups="drop")

# collapse to taxid level
ataxid <- asample %>%
  group_by(Species,taxid) %>%
  mutate(n.samps=n_distinct(sample_ids),
         n.samps.with.usable=sum(n.usable>0),
         n.samps.with.usable.hits=sum(!is.na(s.hit))) %>%
  group_by(n.samps,
           n.samps.with.usable,
           n.samps.with.usable.hits,
           .add=TRUE) %>%
  summarize(s.hit.breakdown=ifelse(n.samps.with.usable.hits[1],
                                   tab(s.hit[!is.na(s.hit)],as.char = TRUE,pct=FALSE),
                                   NA_character_),
            .groups="drop") %>%
  mutate(s.hit.simple=map_chr(str_split(s.hit.breakdown,"\n"),first),
         s.hit.simple=str_replace(s.hit.simple," \\(n=[0-9]+\\)",""),
         s.hit.gut.group=case_when(
           is.na(s.hit.breakdown) ~ NA_character_,
           s.hit.breakdown %ilike% "gut" ~ "gut",
           TRUE ~ "non-gut"
         ))

write_rds(ataxid,file="data/ataxid.rds")








