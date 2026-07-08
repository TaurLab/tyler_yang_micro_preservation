
# set up data (run) -------------------------------------------------------------

library(yingtools2)
library(tidyverse)
library(phyloseq)
library(ggtree)
library(ggh4x)
library(vegan)
library(pairwiseAdonis)
library(ggpubr)
library(ggrepel)
library(decontam)
library(emmeans) # for contrasts
rm(list=ls())

load("data/phy.tyler.RData")
load("data/phy.tyler.additional.RData")
phy.others = read_rds("data/other.samps.rds")

# correct negative phylo tree lengths
tr <- phy_tree(phy.tyler)
tr$edge.length[tr$edge.length<0] <- 0
phy_tree(phy.tyler) <- tr
rm(tr)

# add beta metric, using samp.comparator

# remove additional ranks
phy.tyler <- phy.tyler %>% 
  select(otu,Superkingdom,Phylum,Class,Order,Family,Genus,Species)

# modify sample_data
s.tyler <- phy.tyler %>% get.samp() %>%
  left_join(trace_tbl.tyler,by="sample") %>%
  mutate(temp=ifelse(sample %in% c("1A","1B"),"n/a",temp),
         temp=factor(temp,levels=c("n/a","-80C", "-20C", "4C", "room temp")),
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
         lbl.comparator=lbl[match(sample.comparator,sample)],
         temp.abbrev=case_match(temp,"room temp"~"RT",.default=temp),
         treatment.abbrev=str_replace(treatment,"autoclave","auto"),
         time.abbrev=paste0("d",days),
         lbl2=str_glue("{lbl} ({treatment.abbrev}|{time.abbrev},{temp.abbrev})"),
         lbl3=ifelse(experiment==1,
                     str_glue("{lbl} ({time.abbrev},{temp.abbrev})"),
                     str_glue("{lbl} ({treatment.abbrev}|{time.abbrev})"))) %>%
  ungroup() 
sample_data(phy.tyler) <- s.tyler %>% set.samp()

# custom palette
pal <- list("Bacteroidota (phylum)"=Phylum %in% c("Bacteroidetes", "Bacteroidota") ~ shades("#51AB9B", variation = 0.25),
            "Lachnospiraceae (family)"=Family == "Lachnospiraceae" ~ shades("#EC9B96", variation = 0.25),
            "Oscillospiraceae (family)"=Family %in% c("Ruminococcaceae", "Oscillospiraceae") ~ shades("#9AAE73", variation = 0.25),
            "Other Eubacteriales (order)"=Order %in% c("Clostridiales", "Eubacteriales") ~ shades("#9C854E", variation = 0.25),
            "Actinomycetota (phylum)"=Phylum %in% c("Actinobacteria", "Actinomycetota") ~ shades("#A77097", variation = 0.25),
            # "Enterococcus (genus)"=Genus == "Enterococcus" ~ shades("#129246", variation = 0.15),
            # "Streptococcus (genus)"=Genus == "Streptococcus" ~ shades("#9FB846", variation = 0.15),
            # "Staphylococcus (genus)"=Genus == "Staphylococcus" ~ shades("#f1eb25", variation = 0.15),
            # "Lactobacillus (genus)"=Genus == "Lactobacillus" ~ shades("#3b51a3", variation = 0.15),
            "Pseudomonadota (phylum)"=Phylum %in% c("Proteobacteria", "Pseudomonadota") ~ shades("red", variation = 0.4),
            "Other Bacteria"=TRUE ~ shades("gray", variation = 0.25))

add_dist <- function(phy,method,comparator=NULL,varname=NULL) {
  # phy=phy1;method="horn";sample0="1A"
  varname <- enquo(varname)
  if (rlang::quo_is_null(varname)) {
    varname <- paste0("dist_",method)  
  } else {
    varname <- as_label(varname)
  }
  s <- get.samp(phy)
  samp2 <- s$sample
  if (is.null(comparator)) {
    samp1 <- s$sample.comparator  
  } else {
    samp1 <- rep_along(samp2,comparator)
  }
  pw <- calc.pairwise(sample1=samp1,sample2=samp2,method=method,phy=phy)
  new.s <- s %>% mutate(!!varname:=pw)
  sample_data(phy) <- new.s %>% set.samp()
  cli::cli_alert("added var: {.pkg {varname}} (comparator={unique(samp1)})")
  return(phy)
}

# generate phy1 and phy2 (run) -----------------------------------------------------------

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
                     # time=fct_recode(time,"Day 0 (Baseline)"="Day 0"),
                     letter.rev=fct_rev(letter),
                     xtitle.lbl="'Storage Temperature / Storage Time'",
                     time.lbl=fct_relabel(time,function(x) paste0("'",x,"'")),
                     temp.lbl=fct_recode(temp,"italic('(n/a)')"="n/a",
                                         "-80*degree*C"="-80C",
                                         "-20*degree*C"="-20C",
                                         "4*degree*C"="4C",
                                         "'room temp'"="room temp"),
                     lbl=fct_reordern(lbl,lbl),
                     lbl2=fct_reordern(lbl2,lbl),
                     lbl3=fct_reordern(lbl3,lbl)) %>%
  add_dist("pct.bray") %>%
  add_dist("horn") %>%
  add_dist("mean.horn") %>%
  add_dist("unfold.horn")


phy2.contam <- phy.tyler %>%
  filter(experiment==2) %>%
  mutate(baseline=sample==sample.comparator,
         time.num=days,
         time.rank=dense_rank(time.num),
         uv=fct_relevel(uv,"no UV"),
         heat=fct_relevel(heat,"no heat"),
         heat.75C=heat=="75C",
         heat.autoclave=heat=="autoclave",
         uv.dna=uv=="UV DNA",
         uv.sample=uv=="UV",
         # qpcr.totalseqs=coalesce(qpcr.totalseqs,100),
         # for formatting purposes
         timelabel="Storage Time",
         treatmentlabel="Treatment",
         letter=factor(1),
         xtitle.lbl="'Sample Treatment / Open Storage Time'",
         time.lbl=fct_relabel(time,function(x) paste0("'",x,"'")),
         treatment.lbl=fct_recode(treatment,
                                  "italic(none)"="none",
                                  "75*degree*C"="75C",
                                  "'UV'"="UV",
                                  "75*degree*C+UV"="75C+UV",
                                  "'autoclave'"="autoclave",
                                  "'autoclave'+'UV'"="autoclave+UV",
                                  "'UV DNA'"="UV DNA"))

phy.pcrneg.control <- read_rds("data/phy.pcrneg.control.rds") %>%
  select_tax_table(otu,Superkingdom,Phylum,Class,Order,Family,Genus,Species) %>%
  mutate_sample_data(sample=str_replace(sample,"..pool1161",""),
                     lbl="PCRNeg",lbl2="PCRNeg",lbl3="PCRNeg",
                     treatment="PCRNeg")

phy2.with.pcrneg <- phy.combine(phy2.contam,phy.pcrneg.control) %>%
  mutate(lbl=fct_reordern(lbl,lbl),
         lbl2=fct_reordern(lbl2,lbl),
         lbl3=fct_reordern(lbl3,lbl),
         treatment=factor(treatment,levels=c("none", "75C", "UV", "75C+UV", "autoclave", "autoclave+UV", "UV DNA", "PCRNeg")),
         is.neg.control=sample=="PCRNeg4") %>%
  mutate_tax_table(not.contam=isNotContaminant(.,neg="is.neg.control"))

# extract taxids from tax.blast.tyler and add to phy2a
# tax.taxid <- tax.blast.tyler %>%
#   group_by(otu) %>% arrange(evalue.rank,staxid) %>% slice(1) %>% ungroup() %>%
#   transmute(otu,taxid=as.character(staxid))
# tax.2a <- phy2a %>% get.tax() %>% left_join(tax.taxid,by="otu")
# tax_table(phy2a) <- tax.2a %>% set.tax()

phy2 <- phy2.with.pcrneg %>% 
  filter(as.logical(not.contam),!is.neg.control) %>%
  select(-not.contam) %>%
  add_dist("horn",comparator="TY.1_D0_NT",varname=dist_2A) %>%
  add_dist("mean.horn",comparator="TY.1_D0_NT",varname=dist_2A_meanhorn) %>%
  add_dist("horn",comparator="PCRNeg4",varname=dist_PCRNeg)


# functions (run) ---------------------------------------------------------------

p.value.asterisk <- function(p.value) {
  case_when(is.na(p.value) ~ NA_character_,
            p.value>=0.05 ~ "ns",
            p.value>=0.01 ~ "*",
            p.value>=0.001 ~ "**",
            TRUE ~ "***")
  
}

do.permanova.old <- function(phy,dist,form,seed=1) {
  set.seed(seed)
  if (is.character(dist)) {
    .dist <- calc.distance(phy,dist)
  } else {
    .dist <- dist
  }
  s <- get.samp(phy)
  form <- as.formula(paste(".dist",deparse(form)))
  permanova.test <- adonis2(form, data=s, permutations=1e5)
  permanova.tbl <- permanova.test %>% broom::tidy() %>%
    filter(!is.na(statistic)) %>%
    mutate(signif=p.value<0.05,
           stars=p.value.asterisk(p.value))
  
  # beta dispersion test for anything significant
  beta.test <- permanova.tbl %>%
    # filter(signif) %>%
    pull(term) %>%
    setNames(.,.) %>%
    map(~{
      if (.x %notin% names(s)) {
        cli::cli_abort("YTError: term not found")
      }
      bd <- betadisper(.dist, s[[.x]])
      permutest <- permutest(bd)
      permutest
    })
  beta.tbl <- beta.test %>%
    imap(~{
      bd.pval <- .x[["tab"]] %>% 
        rownames_to_column("var") %>%
        filter(var=="Groups") %>% pull(`Pr(>F)`)
      tibble(term=.y,betadisp.p.value=bd.pval)
    }) %>% list_rbind()
  # pairwise tests for anything significant and multilevel
  pairwise.tests <- permanova.tbl %>%
    mutate(n.lvls=map_int(term,~{
      if (.x %in% names(s)) {
        return(n_distinct(s[[.x]]))
      } else {
        return(NA_integer_)
      }
    })) %>% 
    mutate(do.pairwise=signif & !is.na(n.lvls) & n.lvls>2) %>%
    filter(do.pairwise) %>% pull(term) %>%
    setNames(.,.) %>%
    map(~{
      form <- as.formula(paste0(".dist ~ ",.x))
      pairwise <- pairwise.adonis2(form, data=s)
      return(pairwise)
    }) 
  pairwise.tbl <- pairwise.tests %>%
    imap(~{
      .x %>% 
        keep(is.data.frame) %>%
        imap(~{
          .x %>% rownames_to_column("element") %>%
            rename(p.value=`Pr(>F)`) %>%
            mutate(pair=.y) %>%
            filter(!is.na(F)) %>% 
            as_tibble()
        }) %>% list_rbind()
    })
  tbl <- permanova.tbl %>% 
    left_join(beta.tbl,by="term")
  
  tbl.formatted <- tbl %>% 
    mutate(across(.cols=ends_with("p.value"),.fns=~scales::pvalue(.x)),
           term=str_replace_all(term,":","'%*%'")) %>%
    select("Predictor"=term,"R^2"=R2,"italic(P)*'-value'"=p.value,
           betadisp.p.value) %>%
    mutate(across(.cols=where(is.numeric),.fns=pretty_number),
           across(.cols=everything(),.fns=~paste0("'",.x,"'")))
  
  
  obj <- list(
    ord=ord,
    ord.axes=axes,
    tbl=tbl,
    tbl.formatted=tbl.formatted,
    pairwise=pairwise.tbl,
    permanova.test=permanova.test,
    beta.test=beta.test,
    pairwise.tests=pairwise.tests
  )
  return(obj)
}



do.permanova.old2 <- function(phy,dist,form,by="terms",seed=1,permutations=1e5) {
  set.seed(seed)
  if (is.character(dist)) {
    .dist <- calc.distance(phy,dist)
  } else {
    .dist <- dist
  }
  s <- get.samp(phy)
  form <- as.formula(paste(".dist",deparse(form)))
  permanova.test <- adonis2(form, data=s, permutations=permutations,by=by)
  permanova.tbl <- permanova.test %>% broom::tidy() %>%
    filter(!is.na(statistic)) %>%
    mutate(signif=p.value<0.05,
           stars=p.value.asterisk(p.value))
  
  permanova.tbl2 <- permanova.tbl %>%
    mutate(
      beta.dispersion.pvalue=map_dbl(term,function(term) {
        if (!(term %in% names(s))) {return(NA_real_)}
        bd <- betadisper(.dist, s[[term]])
        permutest <- permutest(bd)
        bd.pval <- permutest$tab %>%
          rownames_to_column("var") %>%
          filter(var=="Groups") %>% pull(`Pr(>F)`)
        return(bd.pval)
      }))
  
  tbl <- permanova.tbl2 %>%
    mutate(n.lvls=map_int(term,~{
      if (.x %in% names(s)) {
        return(n_distinct(s[[.x]]))
      } else {
        return(NA_integer_)
      }
    }),
    do.pairwise=signif & !is.na(n.lvls) & n.lvls>2,
    contrasts=map2(term,do.pairwise,function(term,do) {
      form <- as.formula(paste0(".dist ~ ",term))
      pairwise <- pairwise.adonis2(form, data=s)
      tbl <- pairwise %>% 
        keep(is.data.frame) %>%
        imap(~{
          .x %>% rownames_to_column("element") %>%
            rename(p.value=`Pr(>F)`) %>%
            mutate(pair=.y) %>%
            filter(!is.na(F)) %>% 
            as_tibble()
        }) %>% list_rbind()
      return(tbl)
    }))
  names(tbl$contrasts) <- tbl$term
  # perm$tbl$beta.dispersion.pvalue

  tbl.formatted <- tbl %>% 
    mutate(across(.cols=ends_with("p.value"),.fns=~scales::pvalue(.x)),
           # beta.dispersion.pvalue=str_replace_all(beta.dispersion.pvalue,"","")
           term=str_replace_all(term,":","'%*%'")) %>%
    select("Predictor"=term,"R^2"=R2,"italic(P)*'-value'"=p.value,
           "beta~'disper'~italic(P)"=beta.dispersion.pvalue) %>%
    mutate(across(.cols=where(is.numeric),.fns=pretty_number),
           across(.cols=everything(),.fns=~paste0("'",.x,"'")),
           across(.cols=everything(),.fns=~ifelse(.x=="'NA'","''",.x)))

  ord <- phy.ordinate(phy,method="PCoA",distance=dist)
  axes <- ord$obj$values$Rel_corr_eig %>% 
    imap(~{
      var <- paste0("PC",.y)
      pct <- scales::label_percent(accuracy=0.1)(.x)
      str_glue("{var} ({pct})")
    })
  
  list(ord=ord,
       ord.axes=axes,
       tbl=tbl,
       tbl.formatted=tbl.formatted)
}


do.permanova <- function(phy,dist,form,by="terms",seed=1,permutations=1e5) {
  set.seed(seed)
  if (is.character(dist)) {
    .dist <- calc.distance(phy,dist)
  } else {
    .dist <- dist
  }
  s <- get.samp(phy)
  form <- as.formula(paste(".dist",deparse(form)))
  permanova.test <- adonis2(form, data=s, permutations=permutations,by=by)
  permanova.tbl <- permanova.test %>% broom::tidy() %>%
    filter(!is.na(statistic)) %>%
    mutate(signif=p.value<0.05,
           stars=p.value.asterisk(p.value))
  allterms <- permanova.tbl$term
  
  permanova.tbl2 <- permanova.tbl %>%
    mutate(
      beta.dispersion.pvalue=map_dbl(term,function(term) {
        if (!(term %in% names(s))) {return(NA_real_)}
        bd <- betadisper(.dist, s[[term]])
        permutest <- permutest(bd)
        bd.pval <- permutest$tab %>%
          rownames_to_column("var") %>%
          filter(var=="Groups") %>% pull(`Pr(>F)`)
        return(bd.pval)
      }))
  
  tbl <- permanova.tbl2 %>%
    mutate(n.lvls=map_int(term,~{
      if (.x %in% names(s)) {
        return(n_distinct(s[[.x]]))
      } else {
        return(NA_integer_)
      }
    }),
    do.pairwise=signif & !is.na(n.lvls) & n.lvls>2,
    contrasts=map2(term,do.pairwise,function(term,do) {
      # first model is what gets pairwise
      other.terms <- setdiff(allterms,term) 
      model.terms <- c(term,other.terms) %>% paste(collapse=" + ")
      form <- as.formula(paste0(".dist ~ ",model.terms))
      pairwise <- pairwise.adonis2(form, data=s,nperm=permutations)
      tbl <- pairwise %>% 
        keep(is.data.frame) %>%
        imap(~{
          .x %>% rownames_to_column("element") %>%
            rename(p.value=`Pr(>F)`) %>%
            mutate(pair=.y) %>%
            # filter(!is.na(F)) %>% 
            filter(element==term) %>%
            as_tibble()
        }) %>% list_rbind()
      return(tbl)
    }))
  
  names(tbl$contrasts) <- tbl$term
  # perm$tbl$beta.dispersion.pvalue
  
  tbl.formatted <- tbl %>% 
    mutate(across(.cols=ends_with("p.value"),.fns=~scales::pvalue(.x)),
           # beta.dispersion.pvalue=str_replace_all(beta.dispersion.pvalue,"","")
           term=str_replace_all(term,":","'%*%'")) %>%
    select("Predictor"=term,"R^2"=R2,"italic(P)*'-value'"=p.value,
           "beta~'disper'~italic(P)"=beta.dispersion.pvalue) %>%
    mutate(across(.cols=where(is.numeric),.fns=pretty_number),
           across(.cols=everything(),.fns=~paste0("'",.x,"'")),
           across(.cols=everything(),.fns=~ifelse(.x=="'NA'","''",.x)))
  
  ord <- phy.ordinate(phy,method="PCoA",distance=dist)
  axes <- ord$obj$values$Rel_corr_eig %>% 
    imap(~{
      var <- paste0("PC",.y)
      pct <- scales::label_percent(accuracy=0.1)(.x)
      str_glue("{var} ({pct})")
    })
  
  list(ord=ord,
       ord.axes=axes,
       tbl=tbl,
       tbl.formatted=tbl.formatted)
}

adjust.pairs <- function(pairwise.table,pairs=NULL) {
  if (is.null(pairs)) {
    pairs <- pairwise.table$pair
  }
  if (!all(pairs %in% pairwise.table$pair)) {
    cli::cli_abort("YTError: not all pairs found")
  }
  total.pairs <- nrow(pairwise.table)
  npairs <- length(pairs)
  cli::cli_alert_info("Total {total.pairs} pairs, choosing {npairs}: {unname(pairs)}")
  
  ptbl <- pairwise.table %>% 
    filter(pair %in% pairs) %>%
    mutate(pair=factor(pair,levels=pairs),
           adjusted.pval=p.adjust(p.value,"BH"))
  
  if (!is.null(names(pairs))) {
    ptbl <- ptbl %>% mutate(pair=fct_recode(pair,!!!pairs))
  }
  ptbl <- ptbl %>%
    arrange(pair) %>%
    mutate(adjusted.pval=scales::pvalue(adjusted.pval)) %>%
    select("Pair"=pair,"R^2"=R2,"italic(P)*'-value (adj)'"=adjusted.pval) %>%
    mutate(across(.cols=where(is.numeric),.fns=pretty_number),
           across(.cols=everything(),.fns=~paste0("'",.x,"'")))  
  return(ptbl)
}

make.contrasts <- function(vec,pairs) {
  lvls <- levels(vec)
  all.values <- pairs %>% unlist()
  if (!all(all.values %in% lvls)) {
    cli::cli_abort("YTError: not all values found in vector")
  }
  con <- pairs %>%
    map(~{
      case_when(
        lvls==.x[1] ~ -1,
        lvls==.x[2] ~ 1,
        TRUE ~ 0
      )
    })
  return(con)
}

aov_contrast <- function(f,data,...) {
  contrasts <- list(...)
  # perform anova
  qmod <- aov(formula = f, data = data)
  yvar <- rlang::f_lhs(f) %>% rlang::as_label()
  tbl <- broom::tidy(qmod) %>%
    mutate(n.values=map_int(term,~{
      n_distinct(data[[.x]])
    })) %>%
    filter(term!="Residuals") %>%
    mutate(do.contrast=p.value<0.05 & n.values>2,
           contrast=map2(term,do.contrast,function(term,do.contrast) {
             if (!do.contrast) {
               return(tibble(contrasts=list(NULL),tukey=list(NULL)))
             }
             ff <- as.formula(paste(yvar,"~",term))
             qmod.single <- aov(formula=ff,data=data)
             eformula <- paste("~",term) %>% as.formula()
             emm <- emmeans(qmod.single,eformula)
             tukey <- pairs(emm,adjust="tukey")
             planned.contrasts <- contrasts[[term]]
             vec <- data[[term]]
             if (is.null(planned.contrasts)) {
               cli::cli_warn("YTWarning: no contrast specified for {term}, including all pairs")
               planned.contrasts <- combn(levels(vec),2,simplify = FALSE)
               names(planned.contrasts) <- planned.contrasts %>% map_chr(~paste(.x,collapse=" vs "))
             }
             code <- make.contrasts(vec,planned.contrasts)
             con <- contrast(emm,code) %>% as_tibble()
             tibble(contrasts=list(con),tukey=list(tukey))
           })) %>%
    unnest(contrast,keep_empty = TRUE)
  names(tbl$tukey) <- tbl$term
  names(tbl$contrasts) <- tbl$term
  return(tbl)
}


p.value.label <- function(p.value) {
  # p.value <- tbl$p.value
  p.formatted <- scales::label_pvalue()(p.value)
  p.operator <- p.formatted %>% str_replace("(^)(?!(<))","=")
  p.final <- str_glue("italic(P)*'{p.operator}'")
  p.final
}

anova_oneline <- function(tbl,use.stars=FALSE) {
  
  if (use.stars) {
    pconvert <- p.value.asterisk
  } else {
    pconvert <- p.value.label
  }
  t2 <- tbl %>% 
    select(term,p.value) %>%
    mutate(p.value=pconvert(p.value),
           # p.value=paste0("'",p.value,"'"),
           term=str_replace(term,":","'%*%'"),
           term=str_glue("'{term}: '"),
           text=paste0(term,"*",p.value))
  text <- paste(t2$text,collapse="*'; '*")
  text2 <- paste("'ANOVA  '*",text)
  text2
}

contrast_oneline <- function(tbl,heading="Contrasts",use.stars=FALSE) {
  if (use.stars) {
    pconvert <- p.value.asterisk
  } else {
    pconvert <- p.value.label
  }
  t2 <- tbl %>% 
    select(contrast,p.value) %>%
    mutate(p.value=pconvert(p.value),
           contrast=str_glue("'{contrast}: '"),
           text=paste0(contrast,"*",p.value))
  text <- paste(t2$text,collapse="*'; '*")
  text2 <- str_glue("'{heading}: '*{text}")
  text2
}

make.step <- function(samples,phy,dist,rows=2,aspect.ratio=1) {
  
  # samples = c("1B","1E","1Q","1Z");phy=phy1;dist=quo(dist_horn);rows=2
  dist <- enquo(dist)
  
  phy1.compare <- phy %>% 
    mutate(compare=lbl %in% samples,
           dist=!!dist,
           lbl.compare=str_glue("{lbl3} (vs. {lbl.comparator})")) %>%
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
    mutate(label=str_glue("Horn={sprintf('%.3f',dist)}"))
  
  g1.asv <- ggplot() +
    geom_col(data=otu1compare,aes(x=col,y=pctseqs,fill=otu)) +
    geom_step(data=otu1compare,aes(x=col,y=pctseqs0),direction="mid") +
    geom_text(data=s1compare,aes(x=Inf,y=Inf,label=label),hjust=1,vjust=1,color="blue") +
    geom_rect(data=s1compare,aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,linetype=baseline),
              fill=NA,color="blue",show.legend=FALSE) +
    scale_linetype_manual(values=c("TRUE"="longdash","FALSE"=NA)) +
    yingtools2::geom_bracket(data=filter(otu1compare,extra),
                             aes(x=col-0.5,xend=col+0.5,y=ave(pctseqs,sample,FUN=max),
                                 label="unique\nASVs"), fontsize=3,
                             tip="square") +
    scale_fill_taxonomy(name="Bacterial Taxa",data=otu1compare,fill=otu,tax.palette=pal) +
    scale_y_continuous("Relative Abundance",trans=log_epsilon_trans(0.001)) +
    facet_wrap(~lbl.compare,nrow=rows) +
    ggplot2::labs(x="Bacterial strain (ASV)") +
    theme(aspect.ratio=aspect.ratio,
          legend.key.size = unit(0.85,"lines"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          # axis.title = element_blank(), 
          panel.background = element_blank())
  
  g1.asv
}

make.step2 <- function(samples,phy,dist,aspect.ratio=1) {
  # samples = c("1B","1U","1W","1Y");phy=phy1;dist=quo(dist_horn);rows=2;aspect.ratio=1
  # samples <- samples.compare.all
  # dist=quo(dist_horn)
  # phy=phy1
  # rows=2
  # aspect.ratio=1
  dist <- enquo(dist)
  phy1.compare <- phy %>%
    mutate(compare=lbl %in% samples,
           dist=!!dist,
           lbl.compare=str_glue("{lbl3} (vs. {lbl.comparator})")) %>%
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
  maxcol <- max(otu1compare$col)
  
  s1compare <- phy1.compare %>%
    filter(compare,prune_unused_taxa=FALSE) %>%
    get.samp(stats=TRUE) %>%
    mutate(label=str_glue("Sample {lbl}\nHorn={sprintf('%.3f',dist)}"))
  # facetvars1 <- c("xtitle.lbl","temp.lbl","time.lbl")
  # form1 <- paste(facetvars1,collapse="+") %>% paste(".~",.) %>% as.formula()
  # exp1.facet1 <- facet_nested(letter~xtitle.lbl+temp.lbl+time.lbl+sample, space="free_x",scales="free_x",
  #                            nest_line=element_line(),
  #                            labeller=label_parsed,
  #                            resect=unit(3,"pt"),
  #                            strip = strip_nested(background_x = exp1.panel.striprect.x),
  #                            solo_line=TRUE)
  g1.asv <- ggplot() +
    geom_col(data=otu1compare,aes(x=col,y=pctseqs,fill=otu)) +
    geom_step(data=otu1compare,aes(x=col,y=pctseqs0),direction="mid") +
    geom_text(data=s1compare,aes(x=Inf,y=Inf,label=label),hjust=1,vjust=1,color="blue") +
    geom_rect(data=s1compare,aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,linetype=baseline),
              fill=NA,color="blue",show.legend=FALSE) +
    expand_limits(x=maxcol,y=1) +
    scale_linetype_manual(values=c("TRUE"="longdash","FALSE"=NA)) +
    yingtools2::geom_bracket(data=filter(otu1compare,extra),
                             aes(x=col-0.5,xend=col+0.5,y=ave(pctseqs,sample,FUN=max),
                                 label="unique\nASVs"), fontsize=3,
                             tip="square") +
    scale_fill_taxonomy(name="Bacterial Taxa",data=otu1compare,fill=otu,tax.palette=pal) +
    scale_y_continuous("Relative Abundance",trans=log_epsilon_trans(0.001)) +
    # exp1.facet +
    # facet_nested(facet_form,
    #              labeller=label_parsed,
    #              solo_line=TRUE,
    #              # nest_line=element_line()
    #              resect=unit(3,"pt")) +
    # facet_wrap(~lbl.compare,nrow=rows) +
    ggplot2::labs(x="Bacterial strain (ASV)") +
    theme(!!!theme.settings,
          legend.key.size = unit(0.85,"lines"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          # axis.title = element_blank(),
          panel.background = element_blank())
  # panel.spacing.x = exp1.panel.spacing.x
  g1.asv
}




# universal plotting elements (run) ---------------------------------------------

# exp1 nested facets
facetvars1 <- c("xtitle.lbl","temp.lbl","time.lbl")
form1 <- paste(facetvars1,collapse="+") %>% paste(".~",.) %>% as.formula()
# nested facet 1
exp1.panel.striprect.x <- seq_along(facetvars1) %>%
    map(~{
      vars <- facetvars1[1:.x]
      phy1 %>% get.samp() %>% distinct(!!!syms(vars)) %>% arrange(!!!syms(vars)) %>%
        mutate(label=!!sym(vars[.x]),
               level=.x)
    }) %>% bind_rows() %>%
    mutate(rect=map(level,~{
      if (.x==1) {return(element_rect(fill=NA))}
      else {return(element_rect())}
    })) %>% pull(rect)
# nested facet
exp1.facet <- facet_nested(form1, space="free_x",scales="free_x",
                           nest_line=element_line(),
                           labeller=label_parsed,
                           resect=unit(3,"pt"),
                           strip = strip_nested(background_x = exp1.panel.striprect.x),
                           solo_line=TRUE)
# blank version of facet
exp1.facet.blank <- facet_nested(form1, space="free_x",scales="free_x",
                                 labeller=label_parsed,
                                 resect=unit(3,"pt"),
                                 strip = strip_nested(background_x = element_blank(),
                                                      text_x = element_blank()),
                                 solo_line=TRUE)
# facet spacing 1
exp1.panel.spacing.x <- get.samp(phy1) %>% distinct(!!!syms(facetvars1)) %>% arrange(!!!syms(facetvars1)) %>%
  mutate(spacing=case_when(temp.lbl!=lead(temp.lbl) ~ 3.5, TRUE ~ 1.5)) %>%
  slice(-n()) %>% pull(spacing) %>% unit("points")



width <- 0.95



facetvars2 <- c("xtitle.lbl","treatment.lbl","time.lbl")
form2 <- paste(facetvars2,collapse="+") %>% paste(".~",.) %>% as.formula()
exp2.panel.spacing.x <- get.samp(phy2) %>% distinct(!!!syms(facetvars2)) %>% arrange(!!!syms(facetvars2)) %>%
  mutate(spacing=case_when(treatment.lbl!=lead(treatment.lbl) ~ 3.5, TRUE ~ 1.5)) %>%
  slice(-n()) %>% pull(spacing) %>% unit("points")
exp2.panel.striprect.x <- seq_along(facetvars2) %>%
  map(~{
    vars <- facetvars2[1:.x]
    get.samp(phy2) %>% distinct(!!!syms(vars)) %>% arrange(!!!syms(vars)) %>%
      mutate(label=!!sym(vars[.x]),
             level=.x)
  }) %>% bind_rows() %>%
  mutate(rect=map(level,~{
    if (.x==1) {return(element_rect(fill=NA))} 
    else {return(element_rect())}
  })) %>% pull(rect)

exp2.facet <- facet_nested(form2,scales="free_x",space="free_x",
                           nest_line=element_line(),
                           labeller=label_parsed,
                           resect=unit(3,"pt"),
                           strip = strip_nested(background_x = exp2.panel.striprect.x),
                           solo_line=TRUE)



exp2.facet.blank <- facet_nested(form2,scales="free_x",space="free_x",
                                 labeller=label_parsed,
                                 resect=unit(3,"pt"),
                                 strip = strip_nested(background_x = element_blank(),
                                                      text_x = element_blank()),
                                 solo_line=TRUE)

theme.settings <- list(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())

width <- 0.95

tt <- ttheme(base_style="light",base_size=9,
             tbody.style = tbody_style(size=9,fill=NA,parse=TRUE),
             colnames.style = colnames_style(size=9,fill=NA,parse=TRUE))

# ********* permanova, manual *********
# set.seed(1)
# permanova <- adonis2(dist1 ~ temp*time, data=s1b, permutations = 100000)
# permanova
# # beta dispersion not signif, therefore differences are not from spread, can proceed
# bd <- betadisper(dist1, s1b$time)
# permutest(bd)
# # pairwise contrasts of time
# # pairwise0 <- pairwise.adonis(dist1,factors=s1$time, p.adjust.m="BH")
# pairwise <- pairwise.adonis2(dist1 ~ time, data=s1b) %>%
#   keep(is.data.frame) %>%
#   imap(~{
#     .x %>% rownames_to_column("element") %>%
#       rename(p.value=`Pr(>F)`) %>%
#       mutate(pair=.y) %>%
#       filter(!is.na(F)) %>% as_tibble()
#   }) %>% list_rbind()
# pairwise
# ptbl <- permanova %>% broom::tidy() %>%
#   filter(!is.na(p.value)) %>%
#   mutate(term=str_replace_all(term,":"," %*% ")) %>%
#   select("Predictor"=term,"R^2"=R2,"italic(P)*'-value'"=p.value) %>%
#   mutate(across(.cols=where(is.numeric),.fns=pretty_number))
# target.pairs <- c("Day 0_vs_Day 3" = "Day 0 vs Day 3",
#                   "Day 3_vs_Day 8" = "Day 3 vs Day 8",
#                   "Day 11_vs_Day 8" = "Day 8 vs Day 11")
# ptbl2 <- pairwise %>%
#   filter(pair %in% names(target.pairs)) %>%
#   mutate(pair=recode2(pair,target.pairs,as.factor=TRUE),
#          adjusted.pval=p.adjust(p.value,"BH")) %>%
#   arrange(pair) %>%
#   select("Pair"=pair,"R^2"=R2,"italic(P)*'-value (adj)'"=adjusted.pval) %>%
#   mutate(across(.cols=where(is.numeric),.fns=pretty_number))
# ptbl2
# ***************************************


# fig 1: exp1 taxplot and distance compare --------------------------------

fig.1 <- local({

  phy1a <- phy1 %>%
    mutate(x=as.numeric(lbl),
           dist=dist_horn,
           baseline=lbl=="1A")
  otu1a <- phy1a %>% get.otu.melt()
  s1a <- phy1a %>% get.samp()
  
  g1a.tax <- ggplot(otu1a) +
    geom_taxonomy(aes(x=lbl,y=pctseqs,fill=otu,label=Species),width=width,label.split=TRUE) +
    geom_col(data=filter(s1a,lbl=="1A"),aes(x=lbl,y=1,color="baseline sample (1A)"),
             linewidth=1,linetype="longdash",fill=NA,width=width,
             show.legend = FALSE) +
    scale_y_continuous("Relative abundance",
                       expand=FALSE,labels=scales::label_percent()) +
    scale_x_discrete(name="Sample",expand=expansion(add = 0.5)) + 
    scale_color_manual(name="Sample Comparison",values=c("baseline sample (1A)"="red")) +
    scale_fill_taxonomy(name="Bacterial Taxa",data=otu1a,tax.palette=pal,fill=otu) +
    xlab("Sample") +
    exp1.facet.blank +
    theme(!!!theme.settings,
          axis.text.y=element_blank(),axis.ticks.y=element_blank(),
          legend.key.size = unit(0.85,"lines"),
          strip.clip="on",
          panel.spacing.x = exp1.panel.spacing.x)
  # g1a.tax
  
  dtext1 <- s1a %>% filter(sample=="1A") %>%
    mutate(label="baseline\nsample",
           x=1.5,y=0.2,
           xend=1,yend=0.003)
  dtext1.layer <- list(geom_text(data=dtext1,aes(x=x,y=y,label=label),
                                 vjust=-0.02,lineheight=0.8),
                       geom_segment(data=dtext1,aes(x=x,y=y,xend=xend,yend=yend)))
  
  g1a.dist <- ggplot(s1a) +
    expand_limits(y=1) +
    geom_col(aes(x=lbl,y=dist),fill="steelblue",width=width) +
    geom_text(aes(x=lbl,y=dist,label=pretty_number(dist),color=baseline),
              vjust=0,size=3,show.legend = FALSE) +
    dtext1.layer +
    scale_y_continuous(name="Distance\n(from baseline)",
                       expand=expansion(mult=c(0.03,0.05)))  +
    scale_color_manual(values=c("TRUE"="red","FALSE"="black")) +
    xlab("Sample") +
    exp1.facet +
    theme(!!!theme.settings,
          panel.spacing.x = exp1.panel.spacing.x)
  
  g1.tax.dist <- gg.stack2(g1a.dist,g1a.tax,heights=c(1,3))
  environment()
})

# fig.1$g1.tax.dist

# pdf("plots/fig 1 - exp1 taxdist.pdf",width=14,height=7)
# fig.1$g1.tax.dist
# dev.off()

# shell.exec("plots/fig 1 - exp1 taxdist.pdf")


  win.metafile("clipboard", width = 14, height = 7)  # width/height in inches
  print(g)
  dev.off()


g <- fig.1$g1.tax.dist


win.metafile("clipboard", width = 3, height = 5)
print(dev.cur())      # should show a device name/number, confirms metafile is active
# plot(1:10)             # watch closely for any warning/error here
print(fig.1$g1.tax.dist)
dev.off()


sessionInfo()

win.metafile("clipboard", width = 3, height = 5)
plot(1:10)
dev.off()

copy.to.clipboard.gg(fig.1$g1.tax.dist,width=14,height=7)


ggsave("plots/fig S1 - exp1 invsimpson.pdf",
       fig.s1$g.invsimpson, width=12,height=6)


# fig s1: invsimpson ------------------------------------------------------------


fig.s1 <- local({
  s1 <- phy1 %>% get.samp(stats=TRUE)
  itbl <- aov_contrast(InvSimpson ~ time * temp, data = s1)
  itext <- anova_oneline(itbl)
  
  g.invsimpson <- ggplot(s1,aes(x=lbl,y=InvSimpson)) +
    geom_col(fill="steelblue",width=width) +
    exp1.facet + 
    theme(!!!theme.settings,
          panel.spacing.x = exp1.panel.spacing.x) +
    ggplot2::labs(x="Sample", 
                  y="Inverse Simpson index",
                  caption=parse(text=itext))
  environment()
})


# fig.s1$g.invsimpson

# pdf("plots/fig S1 - exp1 invsimpson.pdf",width=12,height=6)
# fig.s1$g.invsimpson
# dev.off()

# shell.exec("plots/fig S1 - exp1 invsimpson.pdf")



# fig 2: exp 1, pcoa and permanova with horn ---------------------------------------------------------------


fig.2 <- local({
  phy1b <- phy1 %>%
    mutate(storage=days!=0,
           time=ordered(time),
           temp=ordered(temp),
           temp.lbl=ordered(temp.lbl),
           time.lbl=ordered(time.lbl))
  dist1 <- calc.distance(phy1b,"horn")
  s1b <- get.samp(phy1b)
  
  
  perm <- do.permanova(phy1b, dist=dist1, ~temp*time, by="terms")
  perm$time.pair.formatted <- perm$tbl$contrasts$time %>%
    adjust.pairs(c("Day 0 vs 3" = "Day 0_vs_Day 3", 
                   "Day 3 vs 8" = "Day 3_vs_Day 8", 
                   "Day 8 vs 11" = "Day 11_vs_Day 8"))
  
  gtbla <- perm$tbl.formatted %>% 
    select(-`beta~'disper'~italic(P)`) %>%
    ggtexttable(rows=NULL,theme=tt)
  gtblb <- perm$time.pair.formatted %>% 
    rename(`'Time pair'`=Pair) %>%
    ggtexttable(rows=NULL,theme=tt)
  
  size.scaling <- s1b$days %>% unique() %>% sort() %>% 
    scales::rescale(to=c(2,6))
  
  
  # perm$ord$data$temp.lbl

  values <- 21:25
  breaks <- levels(perm$ord$data$temp.lbl)
  labels <- parse(text=breaks)  
  
  g.fig2.pcoa <- perm$ord$data %>%
    arrange(lbl) %>%
    ggplot(aes(x=axis1,y=axis2)) +
    geom_point(aes(fill=temp.lbl,
                   color=temp.lbl,
                   shape=temp.lbl,
                   size=time)) + 
    geom_text(aes(label=lbl),size=3,vjust=1.4) +
    scale_size_manual("Storage Time",values=size.scaling) +
    scale_shape_manual("Storage Temp",values=values,breaks=breaks,labels=labels) +
    scale_color_ordinal("Storage Temp",breaks=breaks,labels=labels) +
    scale_fill_ordinal("Storage Temp",breaks=breaks,labels=labels) +
    xlab(perm$ord.axes[[1]]) + ylab(perm$ord.axes[[2]]) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    theme(!!!theme.settings,
          aspect.ratio=1)

  g.fig2.pcoa

  pos1 <- c(0.75,0.37)
  pos2 <- c(0.75,0.17)
  g.fig2 <- g.fig2.pcoa + 
    patchwork::inset_element(gtbla,pos1[1],pos1[2],pos1[1],pos1[2]) +
    patchwork::inset_element(gtblb,pos2[1],pos2[2],pos2[1],pos2[2])
  g.fig2
  environment()
})


# fig.2$g.fig2

# pdf("plots/fig 2 - exp1 pcoa permanova horn.pdf",width=9,height=9)
# fig.2$g.fig2
# dev.off()
# shell.exec("plots/fig 2 - exp1 pcoa permanova horn.pdf")



# fig s2, pcoa and permanova using bray curtis -------------------------------


fig.s2 <- local({
  
  phy1b <- phy1 %>%
    mutate(storage=days!=0,
           time=ordered(time),
           temp=ordered(temp),
           temp.lbl=ordered(temp.lbl),
           time.lbl=ordered(time.lbl))
  dist1.bray <- calc.distance(phy1b,"pct.bray")
  s1b <- get.samp(phy1b)
  
  perm.bray <- do.permanova(phy1b, dist=dist1.bray, ~temp*time)
  
  
  perm.bray$time.pair.formatted <- perm.bray$tbl$contrasts$time %>%
    adjust.pairs(c("Day 0 vs 3" = "Day 0_vs_Day 3", 
                   "Day 3 vs 8" = "Day 3_vs_Day 8", 
                   "Day 8 vs 11" = "Day 11_vs_Day 8"))
  gtbla.bray <- perm.bray$tbl.formatted %>% 
    select(-`beta~'disper'~italic(P)`) %>%
    ggtexttable(rows=NULL,theme=tt)
  gtblb.bray <- perm.bray$time.pair.formatted %>% 
    rename(`'Time pair'`=Pair) %>%
    ggtexttable(rows=NULL,theme=tt)
  
  size.scaling <- s1b$days %>% unique() %>% sort() %>% 
    scales::rescale(to=c(2,6))
  
  values <- 21:25
  breaks <- levels(perm.bray$ord$data$temp.lbl)
  labels <- parse(text=breaks) 
  
  g.fig.s2.pcoa.bray <- perm.bray$ord$data %>%
    arrange(lbl) %>%
    ggplot(aes(x=axis1,y=axis2)) +
    geom_point(aes(fill=temp.lbl,
                   color=temp.lbl,
                   shape=temp.lbl,
                   size=time)) + 
    geom_text(aes(label=lbl),size=3,vjust=1.4) +
    scale_size_manual("Storage Time",values=size.scaling) +
    # scale_shape_manual(values=c("n/a"=21, "-80C"=22, "-20C"=23, "4C"=24, "room temp"=25)) +
    scale_shape_manual("Storage Temp",values=values,breaks=breaks,labels=labels) +
    scale_color_ordinal("Storage Temp",breaks=breaks,labels=labels) +
    scale_fill_ordinal("Storage Temp",breaks=breaks,labels=labels) +
    xlab(perm.bray$ord.axes[[1]]) + ylab(perm.bray$ord.axes[[2]]) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    theme(!!!theme.settings,
          aspect.ratio=1)

  # g.fig.s2.pcoa.bray
  pos1 <- c(0.75,0.91)
  pos2 <- c(0.75,0.73)
  
  g.fig.s2.bray <- g.fig.s2.pcoa.bray +
    patchwork::inset_element(gtbla.bray,pos1[1],pos1[2],pos1[1],pos1[2]) +
    patchwork::inset_element(gtblb.bray,pos2[1],pos2[2],pos2[1],pos2[2])
  
  g.fig.s2.bray
  environment()
})

# fig.s2$g.fig.s2.bray
# pdf("plots/fig S2 - exp1 pcoa permanova bray.pdf",width=9,height=9)
# fig.s2$g.fig.s2.bray
# dev.off()
# shell.exec("plots/fig S2 - exp1 pcoa permanova bray.pdf")


# Fig s3, QPCR --------------------------------------------------------------------

fig.s3 <- local({
  s1b <- get.samp(phy1,stats=TRUE) %>%
    mutate(log.qpcr.totalseqs=log(qpcr.totalseqs))
  aov.qpcr <- aov_contrast(log.qpcr.totalseqs ~ time * temp, data = s1b,
                           time=list("Day 0 vs 3" = c("Day 0","Day 3"),
                                     "Day 3 vs 8" = c("Day 3","Day 8"),
                                     "Day 8 vs 11" = c("Day 8","Day 11")))
  qtext <- anova_oneline(aov.qpcr)
  qtext2 <- contrast_oneline(aov.qpcr$contrasts$time,"Time contrasts")
  qtext12 <- str_glue("atop({qtext},{qtext2})")
  
  
  g1.qpcr <- ggplot(s1b) +
    geom_col(aes(x=lbl,y=qpcr.totalseqs),fill="steelblue",width=width) +
    exp1.facet + 
    scale_y_continuous("Total bacterial load:\n16S qPCR (copies/mL)",
                       trans=log_epsilon_trans(1000),
                       limits=c(0,1e10),
                       # expand=expansion(mult=0.025),
                       labels=pretty_power10) +
    theme(!!!theme.settings,
          panel.spacing.x = exp1.panel.spacing.x) +
    ggplot2::labs(x="Sample", 
                  caption=parse(text=qtext12))
  environment()
})


# fig.s3$g1.qpcr

# pdf("plots/fig S3 - exp1 qPCR.pdf",width=12,height=6)
# fig.s3$g1.qpcr
# dev.off()
# shell.exec("plots/fig S3 - exp1 qPCR.pdf")


# fig 3: step compare exp1 -------------------------------------------------



fig.3 <- local({
  samples.compare <- c("1B","1U","1X","1Y")
  g1.asv <- make.step2(samples.compare,dist=dist_horn,phy=phy1) + 
    exp1.facet
  environment()

})


# fig.3$g1.asv

# pdf("plots/fig 3 - exp1 step selected.pdf",width=12,height=3.5)
# fig.3$g1.asv
# dev.off()

# shell.exec("plots/fig 3 - exp1 step selected.pdf")


# fig S4: step compare exp1: all ------------------------------------------------------------


fig.s4 <- local({
  samples.compare.all <- phy1 %>% get.samp() %>% pull(lbl)
  form1 <- letter ~ xtitle.lbl + temp.lbl + time.lbl
  # exp1.facet.all <- facet_nested(form1, 
  #                                # space="free_x",scales="free_x",
  #                                # nest_line=element_line(),
  #                                labeller=label_parsed,
  #                                resect=unit(3,"pt"),
  #                                strip = strip_nested(background_x = exp1.panel.striprect.x,
  #                                                     background_y = element_blank(),
  #                                                     text_y=element_blank()),
  #                                solo_line=TRUE)
  # g1.asv.all <- make.step2(samples.compare.all,dist=dist_horn,phy=phy1) +
  #   exp1.facet.all
  title.keep <- 1
  title.blanks <- 2:4
  blanks <- c(23:48)
  background_x <- rep(list(element_rect()),48)
  background_x[blanks] <- rep(list(element_blank()),length(blanks))
  background_x[title.blanks] <- rep(list(element_blank()),length(title.blanks))
  background_x[title.keep] <- rep(list(element_blank()),length(title.keep))
  text_x <- rep(list(element_text()),48)
  text_x[blanks] <- rep(list(element_blank()),length(blanks))
  text_x[title.blanks] <- rep(list(element_blank()),length(title.blanks))
  
  g1.asv.all <- make.step2(samples.compare.all,dist=dist_horn,phy=phy1) +
    facet_manual(. ~ xtitle.lbl + temp.lbl + time.lbl + letter, 
                 design=design,
                 labeller=label_parsed,
                 trim_blank = F,
                 strip=strip_nested(background_x = background_x,
                                    text_x = text_x))
  # g1.asv.all
  environment()
})


# fig.s4$g1.asv.all

# fig.s4$g1.asv.all




# ggsave("plots/fig S4 - exp1 step all.pdf",fig.s4$g1.asv.all,width=22,height=14)
# 
# shell.exec("plots/fig S4 - exp1 step all.pdf")

# shell.exec("plots/fig S4 - exp1 step all.pdf")


# table S1 exp1 lefse -------------------------------------------------------------------


tbl.s1 <- local({
  phy1.lefse <- phy1 %>%   
    mutate(time.group=fct_recode(time,
                                 "Day 0-8"="Day 0",
                                 "Day 0-8"="Day 3",
                                 "Day 0-8"="Day 8",
                                 "Day 11"="Day 11"))
  
  # s1 <- phy1.lefse %>% get.samp()
  
  set.seed(1)
  lda <- lda.effect(phy1.lefse,class="time.group",lda.cutoff=3) 
  
  # the table which can be printed
  lda.tbl <- lda %>%
    filter(pass) %>%
    arrange(direction) %>% 
    select(direction,rank,taxonomy,taxon,taxrank,lda,kw.pvalue) %>%
    mutate(across(.cols=where(is.numeric), .fns = ~sprintf("%.3f",.x)))
  # (just for summarizing)
  # info <- get_taxonomy_info(phy1.lefse,pct.cutoff = 0.99)
  # lda.info <- lda %>% left_join(select(info,taxonomy,grouper,mean.pct.of.parent),by="taxonomy")
  # lda.tbl.final <- lda.info %>% group_by(direction,grouper) %>% 
  #   slice(which.max(rank))
  # lda.tbl.final %>% group_by(direction) %>% arrange(taxonomy) %>% dt()

  environment()
})


# tbl.s1$lda %>% lda.plot()
tbl.s1$lda %>% lda.clado()

# tbl.s1$lda.tbl %>% group_by(direction) %>% dt(fontsize=10)
# tbl.s1$lda.tbl %>% copy.to.clipboard()


# fig 4: exp2 taxplot and distance compare -------------------------------

fig.4 <- local({
  phy2a <- phy2 %>%
    mutate(x=as.numeric(lbl),
           baseline=lbl=="1A",
           dist=dist_2A)
  s2a <- get.samp(phy2a)
  otu2a <- get.otu.melt(phy2a)
  
  g2a.tax <- ggplot(otu2a) +
    geom_taxonomy(aes(x=lbl,y=pctseqs,label=Species,fill=otu),
                  tax.palette=pal,label.split=TRUE,fontsize=3,width=width) +
    geom_col(data=filter(s2a,lbl=="2A"),aes(x=lbl,y=1,color="baseline sample (2A)"),
             linewidth=1,linetype="longdash",fill=NA,width=width,
             show.legend = FALSE) +
    scale_x_discrete("Sample",expand=expansion(add=0.5)) +
    scale_y_continuous("Bacterial Composition\nRelative abundance",expand=FALSE) +
    scale_color_manual("Sample Comparison",values=c("baseline sample (2A)"="red")) +
    scale_fill_taxonomy(name="Bacterial taxa",data=otu2a,fill=otu,tax.palette=pal) +
    exp2.facet.blank +
    theme(!!!theme.settings,
          legend.key.size = unit(0.85,"lines"),
          panel.spacing.x = exp2.panel.spacing.x,
          axis.text.y=element_blank(),
          axis.ticks=element_blank())
  # g2a.tax
  
  s2a %>% filter(sample=="TY.1_D0_NT") %>% pull(sample)
  dtext2 <- s2a %>% filter(lbl=="2A") %>%
    mutate(label="baseline\nsample",
           x=1.5,y=0.2,
           xend=1,yend=0.003)
  dtext2.layer <- list(geom_text(data=dtext2,aes(x=x,y=y,label=label),vjust=-0.02,lineheight=0.8),
                       geom_segment(data=dtext2,aes(x=x,y=y,xend=xend,yend=yend)))
  g2a.dist <- ggplot(s2a) +
    geom_col(aes(x=lbl,y=dist),fill="steelblue",width=width) +
    geom_text(aes(x=lbl,y=dist,label=pretty_number(dist),color=baseline),
              vjust=0,size=3,show.legend = FALSE) +
    dtext2.layer +
    scale_y_continuous(name="Distance\n(from baseline)",
                       expand=expansion(mult=c(0,0.05)))  +
    scale_color_manual(values=c("TRUE"="red","FALSE"="black")) +
    xlab("Sample") +
    exp2.facet +
    theme(!!!theme.settings,
          panel.spacing.x = exp2.panel.spacing.x)
  # g2a.dist
  
  g2.tax.dist <- gg.stack2(g2a.dist,g2a.tax,heights=c(1,3))
  environment()
})

# fig.4$g2.tax.dist

# pdf("plots/fig 4 - exp2 taxdist.pdf",width=12,height=7)
# fig.4$g2.tax.dist
# dev.off()

# shell.exec("plots/fig 4 - exp2 taxdist.pdf")

# fig 5, exp2 pcoa permanova horn ----------------------------------------------


fig.5 <- local({
  phy2b <- phy2
  dist2 <- calc.distance(phy2b,"horn")
  s2b <- get.samp(phy2b)

  
  # perm2 <- do.permanova(phy2b,dist2, ~time + heat + uv, by="margin")
  perm2 <- do.permanova(phy2b,dist2, ~time + heat.75C + uv.sample + heat.autoclave + uv.dna, by="margin")
  
  perm2$tbl.formatted <- perm2$tbl.formatted %>%
    select(-`beta~'disper'~italic(P)`) %>%
    mutate(Predictor=fct_recode(Predictor,
                                "'75'*degree*'C'" = "'heat.75C'", 
                                "'autoclave'" = "'heat.autoclave'", 
                                "'time'" = "'time'", 
                                "'UV DNA'" = "'uv.dna'", 
                                "'UV'" = "'uv.sample'"
    ))
  gtbl2a <- perm2$tbl.formatted %>% ggtexttable(rows=NULL,theme=tt)
  # gtbl2b <- perm2$heat.pair.formatted %>% ggtexttable(rows=NULL,theme=tt)
  # gtbl2c <- perm2$uv.pair.formatted %>% ggtexttable(rows=NULL,theme=tt)
  size.scaling2 <- s2b$days %>% unique() %>% sort() %>% 
    scales::rescale(to=c(3,5))
  
  
  values <- c(21,22,23,24,23,24,25)
  breaks <- levels(s2b$treatment.lbl)
  labels <- parse(text=breaks)
  
  g.fig8.pcoa <- perm2$ord$data %>%
    arrange(lbl) %>%
    ggplot(aes(x=axis1,y=axis2)) +
    geom_point(aes(color=treatment.lbl,
                   fill=treatment.lbl,
                   shape=treatment.lbl,
                   size=time),alpha=0.7) + 
    geom_text_repel(aes(label=lbl),size=3,vjust=1.4,max.overlaps = Inf) +
    scale_color_brewer("Treatment",type="qual",palette=3,breaks=breaks,labels=labels) +
    scale_fill_brewer("Treatment",type="qual",palette=3,breaks=breaks,labels=labels) +
    scale_shape_manual("Treatment",values=values,breaks=breaks,labels=labels) +
    xlab(perm2$ord.axes[[1]]) + ylab(perm2$ord.axes[[2]]) +
    scale_size_manual("Open Storage Time",values=size.scaling2) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    theme(!!!theme.settings,
          aspect.ratio=1)
  pos1 <- c(0.6,0.67)
  g.fig8 <- g.fig8.pcoa +
    patchwork::inset_element(gtbl2a,pos1[1],pos1[2],pos1[1],pos1[2])

  g.fig8

  
  environment()
})




# fig.5$g.fig8

# pdf("plots/fig 5 - exp2 pcoa permanova horn.pdf",width=9,height=9)
# fig.5$g.fig8
# dev.off()

# shell.exec("plots/fig 5 - exp2 pcoa permanova horn.pdf")



# Fig S5: exp2 pcoa permanova using bray ----------------------------------

fig.s5 <- local({
  
  phy2b <- phy2
  dist2 <- calc.distance(phy2b,"pct.bray")
  s2b <- get.samp(phy2b)
  
  
  # perm2 <- do.permanova(phy2b,dist2, ~time + heat + uv, by="margin")
  perm2 <- do.permanova(phy2b,dist2, ~time + heat.75C + uv.sample + heat.autoclave + uv.dna, by="margin")
  
  perm2$tbl.formatted <- perm2$tbl.formatted %>%
    select(-`beta~'disper'~italic(P)`) %>%
    mutate(Predictor=fct_recode(Predictor,
                                "'75'*degree*'C'" = "'heat.75C'", 
                                "'autoclave'" = "'heat.autoclave'", 
                                "'time'" = "'time'", 
                                "'UV DNA'" = "'uv.dna'", 
                                "'UV'" = "'uv.sample'"
    ))
  gtbl2a <- perm2$tbl.formatted %>% ggtexttable(rows=NULL,theme=tt)
  # gtbl2b <- perm2$heat.pair.formatted %>% ggtexttable(rows=NULL,theme=tt)
  # gtbl2c <- perm2$uv.pair.formatted %>% ggtexttable(rows=NULL,theme=tt)
  size.scaling2 <- s2b$days %>% unique() %>% sort() %>% 
    scales::rescale(to=c(3,5))
  
  
  values <- c(21,22,23,24,23,24,25)
  breaks <- levels(s2b$treatment.lbl)
  labels <- parse(text=breaks)
  
  g.fig8.pcoa <- perm2$ord$data %>%
    arrange(lbl) %>%
    ggplot(aes(x=axis1,y=axis2)) +
    geom_point(aes(color=treatment.lbl,
                   fill=treatment.lbl,
                   shape=treatment.lbl,
                   size=time),alpha=0.7) + 
    geom_text_repel(aes(label=lbl),size=3,vjust=1.4,max.overlaps = Inf) +
    scale_color_brewer("Treatment",type="qual",palette=3,breaks=breaks,labels=labels) +
    scale_fill_brewer("Treatment",type="qual",palette=3,breaks=breaks,labels=labels) +
    scale_shape_manual("Treatment",values=values,breaks=breaks,labels=labels) +
    xlab(perm2$ord.axes[[1]]) + ylab(perm2$ord.axes[[2]]) +
    scale_size_manual("Open Storage Time",values=size.scaling2) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    theme(!!!theme.settings,
          aspect.ratio=1)
  pos1 <- c(0.6,0.67)
  g.fig8 <- g.fig8.pcoa +
    patchwork::inset_element(gtbl2a,pos1[1],pos1[2],pos1[1],pos1[2])
  
  g.fig8

  environment()
})


# fig.s5$g.fig8


# pdf("plots/fig S5 - exp2 pcoa permanova bray.pdf",width=9,height=9)
# fig.s5$g.fig8
# dev.off()

# shell.exec("plots/fig S5 - exp2 pcoa permanova bray.pdf")


# fig S6: exp 2 invsimpson ---------------------------------------------------------

fig.s6 <- local({
  s2 <- phy2 %>% get.samp(stats=TRUE)
  
  aov2.invsimp <- aov_contrast(InvSimpson ~ time + heat.75C + uv.sample + heat.autoclave + uv.dna, data = s2) %>% 
    mutate(term=fct_recode(term,"75'*degree*'C" = "heat.75C", 
                           "autoclave" = "heat.autoclave", 
                           "time" = "time", 
                           "UV DNA" = "uv.dna", 
                           "UV" = "uv.sample"))
  itext2 <- anova_oneline(aov2.invsimp)
  
  g2.invsimpson <- ggplot(s2,aes(x=lbl,y=InvSimpson)) +
    geom_col(fill="steelblue",width=width) +
    exp2.facet + 
    theme(!!!theme.settings,
          panel.spacing.x = exp2.panel.spacing.x) +
    ggplot2::labs(title="Inverse Simpson index",
                  x="Sample", 
                  y="Inverse Simpson index",
                  caption=parse(text=itext2))
  environment()
})


# fig.s6$g2.invsimpson


# pdf("plots/fig S6 - exp2 invsimpson.pdf",width=12,height=6)
# fig.s6$g2.invsimpson
# dev.off()

# shell.exec("plots/fig S6 - exp2 invsimpson.pdf")




# fig S7: exp2 qpcr ----------------------------------------------------------------

fig.s7 <- local({
  s2 <- get.samp(phy2) %>%
    mutate(qpcr.lbl=ifelse(is.na(qpcr.totalseqs),"(undetectable)",NA_character_),
           qpcr.totalseqs.impute=coalesce(qpcr.totalseqs,1000),
           log.qpcr.totalseqs.impute=log(qpcr.totalseqs.impute))
  
  aov2.qpcr <- aov_contrast(log.qpcr.totalseqs.impute ~ time + heat.75C + uv.sample + heat.autoclave + uv.dna, data=s2) %>% 
    mutate(term=fct_recode(term,"75'*degree*'C" = "heat.75C", 
                           "autoclave" = "heat.autoclave", 
                           "time" = "time", 
                           "UV DNA" = "uv.dna", 
                           "UV" = "uv.sample"))
  
  qtext2a <- anova_oneline(aov2.qpcr)
  
  g2.qpcr <- ggplot(s2) + 
    geom_col(aes(x=lbl,y=qpcr.totalseqs),fill="steelblue",width=width) + 
    geom_text(aes(x=lbl,y=0,label=qpcr.lbl),hjust=0,angle=90) +
    exp2.facet +
    # scale_x_discrete("Sample",expand=expansion(add=0.5)) +
    scale_y_continuous("Total bacterial load:\n16S qPCR (copies/mL)",
                       trans=log_epsilon_trans(1000),
                       limits=c(0,1e10),
                       # expand=FALSE,
                       # expand=expansion(mult=0.025),
                       labels=pretty_power10) +
    theme(!!!theme.settings,
          panel.spacing.x=exp2.panel.spacing.x) +
    ggplot2::labs(x="Sample",
                  caption=parse(text=qtext2a))
  environment()
})




# fig.s7$g2.qpcr


# pdf("plots/fig S7 - exp2 qpcr.pdf",width=12,height=6)
# fig.s7$g2.qpcr
# dev.off()

# shell.exec("plots/fig S7 - exp2 qpcr.pdf")


# Fig s8, seq loss ------------------------------------------------------


fig.s8 <- local({
  tbl2.seqs <- phy2 %>% get.samp(stats=TRUE) %>%
    select(sample,lbl,lbl3,treatment,time,input,filtered,denoised1,denoised2,merged,seqtab,nochim,seqtab.asvs,nochim.asvs,nseqs,everything()) %>%
    arrange(lbl)
  
  tbl.vars <- c("input","filtered","seqtab","nochim","nseqs")
  
  
  tbl2.long <- tbl2.seqs %>%
    pivot_longer(cols=all_of(tbl.vars)) %>%
    mutate(label=ifelse(name=="nseqs",as.character(lbl3),NA_character_),
           name=factor(name,levels=tbl.vars),
           name=fct_recode(name,"demultiplexing" = "input", 
                           "filtering" = "filtered", 
                           "denoising" = "seqtab", 
                           "chimera removal" = "nochim",
                           "contaminant removal"="nseqs"))

  # scale_fill_brewer("Treatment",type="qual",palette=3,breaks=breaks,labels=labels) +
  
  breaks <- tbl2.seqs$treatment.lbl %>% levels()
  labels <- parse(text=breaks)

  g.seqloss <- 
    ggplot(tbl2.long,aes(x=name,y=value,color=treatment.lbl,label=label,group=lbl)) +
    geom_point(size=4) +
    expand_limits(y=c(0,1e7)) +
    geom_line(linewidth=1) +
    # geom_text(hjust=0) +
    geom_text_repel(hjust=0,nudge_y=0,color="black",xlim=c(5.25),
                    max.overlaps=Inf) +
    expand_limits(x=6.25) +
    xlab("Sample Sequence Processing Step") +
    scale_y_continuous("Number of sequences",trans=log_epsilon_trans(1e4),
                       labels=pretty_power10) + 
    scale_color_brewer("Treatment",type="qual",palette=3,breaks=breaks,labels=labels) +
    theme(#!!!theme.settings,
      panel.grid.minor=element_blank()
    )

  environment()  
})

# fig.s8$g.seqloss


# pdf("plots/fig S8 - exp2 seqloss.pdf",width=10,height=10)
# fig.s8$g.seqloss
# dev.off()

# shell.exec("plots/fig S8 - exp2 seqloss.pdf")


# fig 6 exp2 step selected --------------------------------------------------------------


fig.6 <- local({
  # day 9 samples
  samples <- c("2I", "2U", "2C", "2F", "2O")
  g.fig <- make.step2(samples,phy=phy2,dist=dist_2A) + exp2.facet
  environment()
})

# fig.6$g.fig

# pdf("plots/fig 6 - exp2 step selected.pdf",width=20,height=4)
# fig.6$g.fig
# dev.off()

# shell.exec("plots/fig 6 - exp2 step selected.pdf")


# fig S9 exp2 step all --------------------------------------------------------------

fig.s9 <- local({
  # samples.all <- phy2 %>% get.samp() %>% pull(lbl)
  # exp2.facet.fixed <- facet_nested(treatment.lbl ~ time.lbl,
  #                                  # scales="free_x",space="free_x",
  #                                  nest_line=element_line(),
  #                                  labeller=label_parsed,
  #                                  resect=unit(3,"pt"),
  #                                  # strip = strip_nested(background_x = exp2.panel.striprect.x),
  #                                  solo_line=TRUE)
  # g.fig <- make.step2(samples.all,phy=phy2,dist=dist_2A) +
  #   exp2.facet.fixed
  
  samples.all <- phy2 %>% get.samp() %>% pull(lbl)
  
  design <- "
    ABCDEF
    GHIJKL
    MNOPQR
    STU###
  "
  title.keep <- 1
  title.blanks <- 2:4
  
  background_x <- rep(list(element_rect()),32)
  background_x[title.blanks] <- rep(list(element_blank()),length(title.blanks))
  background_x[title.keep] <- rep(list(element_blank()),length(title.keep))
  
  text_x <- rep(list(element_text()),32)
  text_x[title.blanks] <- rep(list(element_blank()),length(title.blanks))

  g.fig <- make.step2(samples.all,phy=phy2,dist=dist_2A) +
    facet_manual( . ~ xtitle.lbl + treatment.lbl + time.lbl,
                  design=design,
                  labeller=label_parsed,
                  strip=strip_nested(background_x = background_x,
                                     text_x = text_x))
  environment()
})


# fig.s9$g.fig

# ggsave("plots/fig S9 - exp2 step all.pdf",fig.s9$g.fig,width=25,height=12)
# 
# # pdf("plots/fig S9 - exp2 step all.pdf",width=15,height=15)
# # fig.s9$g.fig
# # dev.off()
# 
# shell.exec("plots/fig S9 - exp2 step all.pdf")



# table s2: lefse exp2 --------------------------------------------------------------

# s2b$heat
# no heat 75C autoclave


tbl.s2 <- local({
  phy2.lefse <- phy2 %>% phy.collapse()
  
  lda.75c <- phy2.lefse %>% 
    filter(treatment %in% c("none","75C")) %>%
    mutate(treatment=as.character(treatment)) %>%
    lda.effect(class="treatment",n_boots=NULL,lda.cutoff = 3)
  # lda.uv <- phy2.lefse %>% 
  #   filter(treatment %in% c("none","UV")) %>%
  #   mutate(treatment=as.character(treatment)) %>%
  #   lda.effect(class="treatment",n_boots=NULL,lda.cutoff = 3)
  lda.autoclave <- phy2.lefse %>% 
    filter(treatment %in% c("none","autoclave")) %>%
    mutate(treatment=as.character(treatment)) %>%
    lda.effect(class="treatment",n_boots=NULL,lda.cutoff = 3)
  
  lda.uvdna <- phy2.lefse %>% 
    filter(treatment %in% c("none","UV DNA")) %>%
    mutate(treatment=as.character(treatment)) %>%
    lda.effect(class="treatment",n_boots=NULL,lda.cutoff = 3)
  
  
  
  make.lda.table <- function(lda) {
    lda %>%
      filter(pass,
             !str_detect(taxonomy,"xxxx$")) %>%
      mutate(taxonomy=str_replace_all(taxonomy,"(\\[|\\]xxxx\\|)","")) %>%
      arrange(direction,taxonomy) %>% 
      select(direction,taxonomy,taxrank,lda,kw.pvalue) %>%
      mutate(across(.cols=where(is.numeric), .fns = ~sprintf("%.3f",.x)))
  }
  
  tbl.autoclave <- lda.autoclave %>% make.lda.table() %>% arrange(desc(direction))
  tbl.uvdna <- lda.uvdna %>% make.lda.table() %>% arrange(desc(direction))
  
  environment()  
})

# tbl.s2$lda.75c %>% lda.plot(tax.label="taxonomy")
# tbl.s2$lda.autoclave %>% lda.plot(tax.label="taxonomy")
# tbl.s2$lda.uvdna %>% lda.plot(tax.label="taxonomy")
# tbl.s2$lda.75c %>% lda.clado()
# tbl.s2$lda.autoclave %>% lda.clado()
# tbl.s2$lda.uvdna %>% lda.clado()
# 
# 
# tbl.autoclave %>% copy.to.clipboard()
# tbl.uvdna %>% copy.to.clipboard()





# generate files ----------------------------------------------------------


if (FALSE) {
  # Fig 1-6, PDFs
  ggsave("plots/fig 1 - exp1 taxdist.pdf",
         fig.1$g1.tax.dist, width=14,height=7)
  ggsave("plots/fig 2 - exp1 pcoa permanova horn.pdf",
         fig.2$g.fig2, width=9,height=9)
  ggsave("plots/fig 3 - exp1 step selected.pdf",
         fig.3$g1.asv,width=12,height=3.5)
  ggsave("plots/fig 4 - exp2 taxdist.pdf",
         fig.4$g2.tax.dist, width=12,height=7)
  ggsave("plots/fig 5 - exp2 pcoa permanova horn.pdf",
         fig.5$g.fig8, width=9,height=9)
  ggsave("plots/fig 6 - exp2 step selected.pdf",
         fig.6$g.fig, width=20,height=4)
  
  # Fig 1-6 EPS version
  ggsave("plots/fig 1 - exp1 taxdist.eps",
         fig.1$g1.tax.dist, width=14,height=7, device = cairo_ps)
  ggsave("plots/fig 2 - exp1 pcoa permanova horn.eps",
         fig.2$g.fig2, width=9,height=9, device = cairo_ps)
  ggsave("plots/fig 3 - exp1 step selected.eps",
         fig.3$g1.asv,width=12,height=3.5, device = cairo_ps)
  ggsave("plots/fig 4 - exp2 taxdist.eps",
         fig.4$g2.tax.dist, width=12,height=7, device = cairo_ps)
  ggsave("plots/fig 5 - exp2 pcoa permanova horn.eps",
         fig.5$g.fig8, width=9,height=9, device = cairo_ps)
  ggsave("plots/fig 6 - exp2 step selected.eps",
         fig.6$g.fig, width=20,height=4, device = cairo_ps)
  
  

  # fig S1-S9, PDF
  ggsave("plots/fig S1 - exp1 invsimpson.pdf",
         fig.s1$g.invsimpson, width=12,height=6)
  ggsave("plots/fig S2 - exp1 pcoa permanova bray.pdf",
         fig.s2$g.fig.s2.bray, width=9,height=9)
  ggsave("plots/fig S3 - exp1 qPCR.pdf",
         fig.s3$g1.qpcr, width=12,height=6)
  ggsave("plots/fig S4 - exp1 step all.pdf",
         fig.s4$g1.asv.all, width=22,height=14)
  ggsave("plots/fig S5 - exp2 pcoa permanova bray.pdf",
         fig.s5$g.fig8, width=9,height=9)
  ggsave("plots/fig S6 - exp2 invsimpson.pdf",
         fig.s6$g2.invsimpson, width=12,height=6)
  ggsave("plots/fig S7 - exp2 qpcr.pdf",
         fig.s7$g2.qpcr, width=12,height=6)
  ggsave("plots/fig S8 - exp2 seqloss.pdf",
         fig.s8$g.seqloss, width=10,height=10)
  ggsave("plots/fig S9 - exp2 step all.pdf",
         fig.s9$g.fig,width=25,height=12)

}













