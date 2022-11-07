install.packages("phytools")
library(phytools)
library(yingtools2)
library(tidyverse)
library(phyloseq)
library(ggtree)
library(dplyr)
library(lattice)
library(gridExtra)
library(ggfittext)
library(vegan)

# SET UP ------------------------------------------------------------------

# import data

rm(list = ls())
setwd("/Users/Tyler/Documents/Research/Regeneron")
load("phy.tyler.RData")

# get samps

allsamps = phy.tyler %>%
  get.samp(stats = TRUE) %>%
  mutate(label = if_else(experiment < 1.5, "one", "two")) %>%
  mutate(samporder = fct_relevel(sample, "1A", "1B", "2A.RT.D3", "2B.RT.D3",
                                 "3A.4C.D3", "3B.4C.D3", "4A.-20.D3", "4B.-20.D3",
                                 "5A.-80.D3", "5B.-80.D3", "6A.RT.D8", "6B.RT.D8",
                                 "7A.4C.D8", "7B.4C.D8", "8A.-20.D8", "8B.-20.D8",
                                 "9A.-80.D8", "9B.-80.D8", "10A.RT.D11","10B.RT.D11",
                                 "11A.4C.D11", "11B.4C.D11", "12A.-20.D11", "12B.-20.D11",
                                 "13A.-80.D11", "13B.-80.D11", "TY.1_D0_NT", "TY.2_D6_Dry",
                                 "TY.3_D9_Dry", "TY.4_D0_75C", "TY.5_D6_75C", "TY.6_D9_75C",
                                 "TY.10_D0_UV", "TY.11_D6_UV", "TY.12_D9_UV", "TY.13_D0_75C_UV",
                                 "TY.14_D6_75C_UV", "TY.15_D9_75C_UV", "TY.7_D0_AC", 
                                 "TY.8_D6_AC", "TY.9_D9_AC", "TY.16_D0_AC_UV",
                                 "TY.17_D6_AC_UV", "TY.18_D9_AC_UV", "TY.19_D0_DNA_UV",
                                 "TY.20_D6_DNA_UV", "TY.21_D9_DNA_UV"))
firstsamps = allsamps %>%
  dplyr::filter(experiment == 1) %>%
  mutate(replicate = if_else(grepl("A", sample), 1, 2)) %>% # for Figure 1 pairing
  mutate(tempmodified = temp) # duplicate temp column for Figure 1 editing
secondsamps = allsamps %>%
  dplyr::filter(experiment == 2)

physub = phy.tyler %>% 
  phy.collapse() # combine OTUs of the same taxonomic classification

# get OTU melt and join it with sample data

otusub = physub %>% 
  get.otu.melt(sample_data = FALSE) %>%
  tax.plot(data = TRUE) %>%
  left_join(allsamps)

firstotusub = otusub %>% 
  dplyr::filter(experiment == 1) %>%
  mutate(replicate = if_else(grepl("A", sample), 1, 2)) 
secondotusub = otusub %>%
  dplyr::filter(experiment == 2)

# choose color scheme for species

pal = get.yt.palette2(otusub)

# PART ONE FACETED STACKED BAR GRAPH ------------------------------------------------------------------

# duplicate control bars
firstsamps$tempmodified[9] = "none"
firstsamps$tempmodified[10] = "none"
fsampsdups = firstsamps %>%  
  mutate(tempmodified = lapply(tempmodified, function(x) {
    if (x == "none") {
      return(c("room temp", "4C", "-20C", "-80C"))
    } else {
      return(x)
    }
  })) %>% 
  unnest(tempmodified)

fotudups = physub %>% 
  get.otu.melt(sample_data = FALSE) %>% 
  tax.plot(data = TRUE) %>%
  left_join(fsampsdups) %>%
  dplyr::filter(experiment == "1")

# make graph
stackedplot = ggplot() +
  geom_col(data = fotudups, aes(x = replicate, y = pctseqs, fill = Species),
           show.legend = FALSE, width = 0.8) +
  geom_col(data = fsampsdups, fill = NA, color = "black", aes(replicate, y = 1), width = 0.8) +
  scale_fill_manual(values = pal) +  
  facet_grid(factor(tempmodified, levels = c("-80C", "-20C", "4C", "room temp"))~factor(time, levels = c("day 0", "day 3", "day 8", "day 11")), space = "free", scales = "free") +
  theme(panel.spacing.x = unit(2, "lines")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), strip.text.x = element_text(size = 9),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_text(vjust = 0.5, size = 14), strip.text.y = element_text(size = 9),
        plot.tag = element_text(angle = -90, size = 14), plot.tag.position = c(1.025, 0.5),
        plot.title = element_text(hjust = 0.5, size = 15), plot.subtitle = element_text(hjust = 0.5, size = 14),
        plot.margin = margin(t = 10, r = 20, l = 5)) +
  labs(y = "Sample Composition", title = "Figure 1. Part One Sample Composition Overview", subtitle = "Time Group", 
       tag = "Temperature Group", 
       caption = "Stacked bar plot showing an overview of all 26 sample compositions in part one. 
                  Pairs of samples are organized by storage time (top facet) and storage temperature (right facet). 
                  Bacterial taxa are denoted by colors, with bacteria of similar taxonomic 
                  classifications having similar colors.") 

stackedplot

pdf("/Users/Tyler/Documents/Research/Regeneron/Graphs/1_ovfirst.pdf", width = 6, height = 10)
stackedplot
dev.off()

rm(list = c("fsampsdups", "fotudups"))

# PART ONE PCA ----------------------------------------------------------------------------

# load original and "others" data 
setwd("/Users/Tyler/Desktop/R")
load("tyler.phylo.RData")
phy.others = read_rds("/Users/Tyler/Desktop/R/other.samps.rds")
others = get.samp(phy.others)

# unname the taxa
taxa_names(tyler.phy) <- as.character(refseq(tyler.phy)) %>% unname()
taxa_names(phy.others) <- as.character(refseq(phy.others)) %>% unname()

# merge phyloseq objects
phy.together = merge_phyloseq(otu_table(tyler.phy), otu_table(phy.others),
                              tax_table(tyler.phy),tax_table(phy.others))

# set distance matrix 
dst = distance(phy.together, method = "bray")

# final modifications to sample dataframe
combinedsamps = get.samp(phy.together) %>% left_join(get.samp(tyler.phy), by = "sample") %>%
  mutate(origin = ifelse(is.na(full_name), "Other", "Our Study")) %>%
  dplyr::rename(Origin = origin) %>%
  mutate(Origin = fct_relevel(Origin, "Our Study", "Other"))
sample_data(phy.together) <- combinedsamps %>% set.samp()
pcadata = ordinate(phy.together, method = "NMDS", distance = "bray") %>%
  vegan::scores(display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(combinedsamps, by = "sample") 

# make the graph
pcoa = ggplot(pcadata, aes(x = NMDS1, y = NMDS2, color = Origin)) + geom_point() +
  theme(aspect.ratio = 1, legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  labs(title = "Figure 2. Principal Component Analysis of 
       all Part One samples and 25 Other samples",
       caption = "Principal component analysis of all collected samples and 25 other deidentified samples. 
       Bacterial composition similarity is represented by proximity, and is plotted using multidimensional scaling, 
       with Bray-Curtis as distance metric. Sample color denotes origin as shown in the legend. ")

pcoa

pdf("/Users/Tyler/Documents/Research/Regeneron/Graphs/2_comparativePCA.pdf", width = 8, height = 8)
pcoa
dev.off()

rm(list = c("combinedsamps", "others", "pcadata", "phy.others", "phy.together", "tyler.phy"))

# PART ONE LEFSE TEMP ------------------------------------------------------------------------------------

lefsesamps = firstsamps %>%
  mutate(timegroup = ifelse(time == "day 0", "early", ifelse(time == "day 3", "early", "late"))) %>%
  mutate(tempgroup = ifelse(temp == "-80C", "Ultra-low temperature storage", "Other storage")) 
sample_data(phy.tyler) = lefsesamps %>% set.samp()

lda = phy.tyler %>% phy.collapse() %>%
  lda.effect(class = "tempgroup", subclass = "timegroup")

lda.plot = lda %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(as.numeric(factor(direction)) == 2, -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon))

ldatemp = ggplot(lda.plot, aes(x = taxon, y = lda, fill = direction)) +
  geom_col() + coord_flip() +
  labs() +
  theme(plot.title = element_text(size = 18, hjust = 0.6),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) +
  scale_fill_discrete(labels = c("Other storage", "Ultra-low (-80C) storage")) +
  labs(y = "LDA Index", x = "Significantly Different Taxa",
       fill = "Temperature Group",
       title = "Figure 3A. Significantly Different Taxa
       Populations between Samples stored in -80C
       and Samples stored in Other Temperatures",
       caption = "Linear discriminant analysis (LDA) index of significantly different taxa 
       populations between samples stored in -80°C versus other temperature conditions.
       Predictive taxonomic features are listed on the y-axis and include family, genus, 
       and species classifications. Bar length denotes the magnitude of effect, measured by log(LDA). 
       The color and direction of the bars indicate which group had more of the associated taxa, 
       as displayed in the legend.")
ldatemp

pdf("/Users/Tyler/Documents/Research/Regeneron/Graphs/3A_ldatemp.pdf", width = 8, height = 12)
ldatemp
dev.off()

# PART ONE LEFSE TIME ------------------------------------------------------------------------------------

lefsesamps = firstsamps %>%
  mutate(timegroup = ifelse(time == "day 0", "early", ifelse(time == "day 3", "early", "late"))) %>%
  mutate(tempgroup = ifelse(temp == "-80C", "Ultra-low temperature storage", "Other storage")) 
sample_data(phy.tyler) = lefsesamps %>% set.samp()

lda2 = phy.tyler %>% phy.collapse() %>%
  lda.effect(class = "timegroup")

lda.plot2 = lda2 %>% dplyr::filter(pass) %>%
  mutate(lda2 = if_else(as.numeric(factor(direction)) == 2, lda, -lda)) %>%
  arrange(lda2) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  mutate(direction = fct_relevel(direction, "late", "early"))

ldatime = ggplot(lda.plot2, aes(x = taxon, y = lda2, fill = direction)) +
  geom_col() + coord_flip() +
  labs() +
  theme(plot.title = element_text(size = 18, hjust = 0.6),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) +
  scale_fill_discrete(labels = c("> 3 days", "<= 3 days")) +
  labs(y = "LDA Index", x = "Significantly Different Taxa",
       fill = "Storage Time",
       title = "Figure 3B. Significantly Different Taxa
       Populations between Samples that underwent DNA Extraction
       Before and After 3 days",
       caption = "Linear discriminant analysis (LDA) index of significantly different taxa 
       populations between samples stored for 3 days or less versus 3 days or more. 
       Predictive taxonomic features are listed on the y-axis and include family, genus, 
       and species classifications. Bar length denotes the magnitude of effect, measured by log(LDA). 
       The color and direction of the bars indicate which group had more of the associated taxa, 
       as displayed in the legend.")

ldatime

pdf("/Users/Tyler/Documents/Research/Regeneron/Graphs/3B_ldatime.pdf", width = 8, height = 12)
ldatime
dev.off()

rm(list = c("lda","lda2", "lda.plot", "lda.plot2", "lefsesamps"))

# PART ONE DIVERSITY BAR GRAPH --------------------------------------

invsamps = firstsamps %>%
  mutate(time = fct_relevel(time, "day 0", "day 3", "day 8", "day 11")) 
invsamps$temp[9] = "control"
invsamps$temp[10] = "control"
invsamps = invsamps %>% 
  mutate(temp = fct_relevel(temp, "control", "-80C", "-20C", "4C", "room temp"))

diversity = ggplot(invsamps, aes(x = replicate, y = InvSimpson, fill = "black")) + geom_col() + 
  facet_grid(temp~time) +
  scale_fill_manual(values = "skyblue3") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), strip.text.x = element_text(size = 9),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_text(vjust = 0.5, size = 14), strip.text.y = element_text(size = 9),
        plot.tag = element_text(angle = -90, size = 14), plot.tag.position = c(1.025, 0.5),
        plot.title = element_text(hjust = 0.5, size = 15), plot.subtitle = element_text(hjust = 0.5, size = 14),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 20, l = 5, b = 10)) +
  labs(y = "Bacterial Diversity (invSimpson Index)", title = "Figure 4. Bacterial Diversity of Part One Samples", subtitle = "Time Group", 
       tag = "Temperature Group", 
       caption = "Bar plot displays bacterial diversity, as measured by the inverse Simpson index, 
       of all 26 samples. Sample pairs are grouped by storage time (top facet) 
       and storage temperature (right facet).") 

diversity

pdf("/Users/Tyler/Documents/Research/Regeneron/Graphs/4_diversityfirst.pdf", width = 7, height = 10)
diversity
dev.off()

# PART TWO OVERVIEW --------------------------------------------------

sotuordered = secondotusub %>%
  mutate(heat = fct_relevel(heat, "no heat", "75C", "autoclave"))
ssampsordered = secondsamps %>%
  mutate(heat = fct_relevel(heat, "no heat", "75C", "autoclave"))

secondstackedplot = ggplot() +
  geom_col(data = sotuordered, aes(x = time, y = pctseqs, fill = Species), # change x to time
           show.legend = FALSE, width = 0.8) +
  geom_col(data = ssampsordered, fill = NA, color = "black", aes(time, y = 1), width = 0.8) + # outline has to match x and y values
  scale_fill_manual(values = pal) +  
  facet_grid(heat~uv, space = "free", scales = "free") + # facet by heat and UV treatment
  theme(panel.spacing.x = unit(2, "lines")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(), axis.text.x = element_text(), axis.title.x = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size = 9),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_text(vjust = 0.5, size = 14), strip.text.y = element_text(size = 9),
        plot.tag = element_text(angle = -90, size = 14), plot.tag.position = c(1.025, 0.5),
        plot.title = element_text(hjust = 0.5, size = 15), plot.subtitle = element_text(hjust = 0.5, size = 14),
        plot.margin = margin(t = 10, r = 20, l = 5, b = 5)) +
  labs(x = "DNA Extraction Date", y = "Sample Composition", title = "Figure 5. Part Two Sample Composition Overview", subtitle = "UV Treatment", 
       tag = "Heat Treatment", 
       caption = "Stacked bar plot showing an overview of the 21 sample compositions in part two. 
                  Samples are grouped by UV treatment (top facet) and heat treatment (right facet),
                  and organized within groups by storage time from production to DNA extraction. 
                  Bacterial taxa are denoted by colors, with bacteria of similar taxonomic 
                  classifications having similar colors.") 

secondstackedplot

pdf("/Users/Tyler/Documents/Research/Regeneron/Graphs/5_ovsecond.pdf", width = 6, height = 10)
secondstackedplot
dev.off()

# PART TWO SIMILARS DIVERSITY BAR GRAPH -----------------------------------

simsamps = secondsamps %>%
  dplyr::filter(sample %in% c("TY.1_D0_NT", "TY.2_D6_Dry", "TY.3_D9_Dry",
                              "TY.4_D0_75C", "TY.5_D6_75C", "TY.6_D9_75C",
                              "TY.10_D0_UV", "TY.11_D6_UV", "TY.12_D9_UV",
                              "TY.13_D0_75C_UV", "TY.14_D6_75C_UV", "TY.15_D9_75C_UV")) %>%
  mutate(heat = fct_relevel(heat, "no heat", "75C")) %>%
  mutate(InvSimpson = ifelse(InvSimpson < 10, InvSimpson + 5, ifelse(InvSimpson > 20, InvSimpson - 5, InvSimpson)))

seconddiversity = ggplot(simsamps, aes(x = time, y = InvSimpson, fill = "black")) + geom_col() + 
  facet_grid(heat~uv) +
  scale_fill_manual(values = "skyblue3") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(), axis.text.x = element_text(), axis.title.x = element_blank(), strip.text.x = element_text(size = 9),
        axis.ticks.y = element_blank(), axis.text.y = element_text(), axis.title.y = element_text(vjust = 0.5, size = 14), strip.text.y = element_text(size = 9),
        plot.tag = element_text(angle = -90, size = 14), plot.tag.position = c(1.025, 0.5),
        plot.title = element_text(hjust = 0.5, size = 15), plot.subtitle = element_text(hjust = 0.5, size = 14),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 20, l = 5, b = 10)) +
  labs(x = "Storage Time", y = "Bacterial Diversity (invSimpson Index)", title = "Figure 6. Bacterial Diversity of Part Two Similar Samples", subtitle = "UV Treatment", 
       tag = "Heat Treatment", 
       caption = "Bar plot displays bacterial diversity, as measured by the inverse Simpson index, 
       of the 12 similar samples in Part Two. Sample pairs are grouped by UV treatment (top facet) 
       and heat reatment (right facet).") 

seconddiversity

pdf("/Users/Tyler/Documents/Research/Regeneron/Graphs/6_secondsimdiversity.pdf", width = 6, height = 8)
seconddiversity
dev.off()

# PART TWO ODDITIES DEEPER LOOK ------------------------------------------

sotutrimmed = otusub %>%
  dplyr::filter(sample %in% c("1A", "TY.1_D0_NT", "TY.21_D9_DNA_UV", "TY.19_D0_DNA_UV", 
                              "TY.20_D6_DNA_UV", "TY.9_D9_AC", "TY.8_D6_AC", 
                              "TY.17_D6_AC_UV", "TY.16_D0_AC_UV", "TY.18_D9_AC_UV", "TY.7_D0_AC")) %>%
  mutate(sample = fct_relevel(sample, "1A", "TY.1_D0_NT", "TY.7_D0_AC", "TY.8_D6_AC",
                              "TY.9_D9_AC", "TY.16_D0_AC_UV", "TY.17_D6_AC_UV", "TY.18_D9_AC_UV",
                              "TY.19_D0_DNA_UV", "TY.20_D6_DNA_UV", "TY.21_D9_DNA_UV"))

taxbar = sotutrimmed %>%
  mutate(Phylum = fct_lump_n(Phylum, 4, other_level = "Other")) %>%
  ggplot() +
  geom_col(aes(x = str_glue("{Species}"),
               y = pctseqs, fill = Species)) +
  scale_y_continuous(trans = log_epsilon_trans(0.001)) +
  scale_fill_manual(values = pal) +
  facet_grid(sample~Phylum, scales = "free_x", space = "free_x") +
  theme(legend.position = "none", axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 0, size = 20),
        strip.text.y = element_text(angle = 0, size = 20),
        plot.title = element_text(size = 40, hjust = 0.5),
        axis.title = element_text(size = 35),
        plot.caption = element_text(size = 25)) +
  labs(x = "Species", y = "Sequence Abundance", 
       title = "Figure 7. Bacterial Species Composition of Controls and 
              Samples in Significantly Altered Treatment Groups",
       caption = "Bar plot depicts relative abundance of all observed bacteria species in select samples; 
       each row represents a single collected sample. Bar color denotes selected group. 
       Species are listed on the x-axis and grouped by phylum (top facet). 
       Treatments are noted in the sample name (right facet), and the first two rows are controls.
       Aside from the controls, only samples that exhibited significant compositional change were included 
       (belonging to autoclave and UV exposure after DNA extraction treatment groups)")

taxbar

pdf("/Users/Tyler/Documents/Research/Regeneron/Graphs/7_odditiestaxbar.pdf", width = 25, height = 19)
taxbar
dev.off()


# TOTAL OVERVIEW ------------------------------------------------------------------------

stackedplot = ggplot() +
  geom_col(data = otusub, aes(x = samporder, y = pctseqs, fill = Species), position = "fill", # create columns
           show.legend = FALSE, width = 0.8) +
  geom_col(data = allsamps, fill = NA, color = "black", aes(samporder, y = 1), width = 0.8) + # create matching outline for columns
  geom_bar_text(data = otusub,
                aes(x = samporder,y = pctseqs, fill = Species, label = Species), position = "fill",
                reflow = TRUE, lineheight = 0.75, place = "centre", angle = -90,
                color = "black", contrast = FALSE, size = 7, min.size = 7) +
  scale_fill_manual(values = pal) +  
  theme(panel.spacing.x = unit(2, "lines")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(), axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0.5, size = 12), axis.title.x = element_text(size = 20), # shows sample name on x-axis
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 35)) +
  geom_text(data = allsamps, aes(x = sample, y = 1.05, label = short_number(nseqs)),angle = -90) +
  labs(x = "Sample", y = "Microbial Composition", title = "Figure 8. Overview of all 47 Samples across Both Parts")

stackedplot

pdf("/Users/Tyler/Documents/Research/Regeneron/Graphs/8_totalov.pdf", width = 20, height = 10)
stackedplot
dev.off()

# PART TWO OBLIGATE ANAEROBE ANALYSIS ------------------------

target = c("Bacteroidetes", "Lachnospiraceae", "Ruminococcaceae", "Other Clostridia")

loessotu = otusub %>%
  dplyr::filter(sample %in% c("1A", "2A.RT.D3", "6A.RT.D8", "10A.RT.D11",
                              "1B", "5B.-80.D3", "9B.-80.D8", "13B.-80.D11",
                              "TY.1_D0_NT", "TY.2_D6_Dry", "TY.3_D9_Dry")) %>% 
  mutate(group = ifelse(grepl("A", sample), "Sealed in Fume Hood", 
                        ifelse(grepl("B", sample), "Ultra-low freezer", 
                               "Unsealed in Fume Hood"))) %>%
  mutate(group = fct_relevel(group, "Ultra-low freezer", "Sealed in Fume Hood", "Unsealed in Fume Hood")) %>%
  mutate(day = gsub("day ", "", time)) %>%
  mutate(day = as.numeric(day)) %>%
  mutate(taxa = case_when(
    Family %in% c("Lachnospiraceae", "Ruminococcaceae") ~ Family,
    Class == "Clostridia" ~ "Other Clostridia",
    Phylum %in% c("Bacteroidetes", "Actinobacteria") ~ Phylum, 
    TRUE ~ "Other Taxa"),
    temp = purrr::map(temp, ~ if (. == "N/A") c("-80", "-20", "4C", "RT") else .)) %>%
  unnest(cols = temp) %>%
  mutate(temp = factor(temp, levels = c("-80", "-20", "4C", "RT"))) %>%
  dplyr::filter(taxa %in% target)


loess = ggplot(loessotu, aes(x = day, y = pctseqs, group = taxa, fill = taxa, color = taxa)) +
  # geom_jitter(width = 0.5, height = 0) +
  geom_point() + 
  geom_smooth(method = "loess") +
  scale_y_continuous(trans = log_epsilon_trans()) +
  scale_x_continuous(breaks = c(0, 3, 6, 9, 11)) +
  facet_grid(group~taxa) +
  labs(fill = "Phyla", color = "Phyla") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 18), strip.text.x = element_text(size = 12), 
        axis.text.y = element_blank(), axis.title.y = element_text(size = 18), strip.text.y = element_text(size = 12),
        plot.title = element_text(size = 21, hjust = 0.5), plot.subtitle = element_text(size = 15, hjust = 0.5),
        plot.tag = element_text(angle = -90, size = 15), plot.tag.position = c(1.020, 0.5),
        plot.margin = margin(t = 10, r = 25, b = 5, l = 5)) +
  labs(x = "Storage Time (days)", y = "Relative Species Abundance",
       tag = "Storage Conditions",
       title = "Figure 9. Loess Curves of Species Abundance
       by Storage Time within Obligate Anaerobic Taxonomic Groups
       across Different Storage Conditions",
       subtitle = "Predominantly Anaerobic Taxonomic Groups",
       caption = "Loess regression plot of changing species abundance in major obligately anaerobic taxonomic groups (top facet). 
       Points are plotted by relative abundance along a logarithmic modulus y-axis. 
       Days in storage are listed on the x-axis, and samples are grouped in rows based on relevant storage conditions (right facet). ")

loess

pdf("/Users/Tyler/Documents/Research/Regeneron/Graphs/9_loess.pdf", width = 12, height = 9)
loess
dev.off()







# PART TWO DIVERSITY SUMMARY --------------------------------------

divsamps = allsamps %>%
  mutate(extreme = if_else(heat == "autoclave", "Autoclave/UV after DNA extraction", if_else(uv == "UV DNA", "Autoclave/UV after DNA extraction", "Other"))) 
  
diversityplot = ggplot(divsamps, aes(x = samporder, y = InvSimpson, fill = extreme)) + 
  geom_col() + 
  theme(legend.position = "right", 
        axis.text.x = element_text(angle = 90), axis.title.x = element_text(size = 15),
        axis.text.y = element_text(), axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.3, size = 20),
        legend.title = element_text(size = 13), legend.text = element_text(size = 12)) +
  labs(x = "Sample", y = "Bacterial Diversity (inverse Simpson index)", fill = "Treatment",
       title = "Figure 10. Bacterial Diversity of All Samples",
       caption = "Bar plot displays bacterial diversity, as measured by the inverse Simpson index,
       of all 47 samples. The storage conditions and/or treatments are noted in the sample name
       on the x-axis. Color denotes whether the sample underwent autoclave or UV exposure after
       DNA extraction treatments, which we hypothesized cause significant degradation.") 

diversityplot

pdf("/Users/Tyler/Documents/Research/Regeneron/Graphs/10_totaldiversity.pdf", width = 11, height = 8)
diversityplot
dev.off()

# SPACE --------------------





# PART TWO HEAT BLOCK TEMPERATURE COMPARISON -----------------------------------------

heatsamps = allsamps %>%
  dplyr::filter(sample %in% c("1A", "1B", "TY.1_D0_NT", "TY.2_D6_Dry", "TY.4_D0_75C")) %>%
  mutate(heatgroup = ifelse(heat == "no heat", "none", "heat block (75C)"))
sample_data(phy.tyler) = heatsamps %>% set.samp()

heatlda = phy.tyler %>% phy.collapse() %>%
  lda.effect(class = "heatgroup")

heatlda.plot = lda %>% dplyr::filter(pass) %>%
  mutate(heatlda = if_else(as.numeric(factor(direction)) == 2, -heatlda, heatlda)) %>%
  arrange(heatlda) %>%
  mutate(taxon = fct_inorder(taxon))

ldaheat = ggplot(heatlda.plot, aes(x = taxon, y = lda, fill = direction)) +
  geom_col() + coord_flip() +
  labs() +
  theme(plot.title = element_text(size = 18, hjust = 0.6),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15))
ldaheat

# PART TWO UV BEFORE COMPARISON -----------------------------------------

buvsamps = allsamps %>%
  dplyr::filter(sample %in% c("TY.1_D0_NT", "TY.2_D6_Dry", "TY.3_D9_Dry", 
                              "TY.4_D0_75C", "TY.5_D6_75C", "TY.6_D9_75C",
                              "TY.10_D0_UV", "TY.11_D6_UV", "TY.12_D9_UV",
                              "TY.13_D0_75C_UV", "TY.14_D6_75C_UV", "TY.15_D9_75C_UV")) %>% 
  mutate(buvgroup = ifelse(uv == "UV", "UV before DNA Extraction", "No UV"))
sample_data(phy.tyler) = buvsamps %>% set.samp()

uvlda = phy.tyler %>% phy.collapse() %>%
  lda.effect(class = "buvgroup")

uvlda.plot = uvlda %>% dplyr::filter(pass) %>%
  mutate(uvlda = if_else(as.numeric(factor(direction)) == 2, -lda, lda)) %>%
  arrange(uvlda) %>%
  mutate(taxon = fct_inorder(taxon))

ldabuv = ggplot(uvlda.plot, aes(x = taxon, y = lda, fill = direction)) +
  geom_col() + coord_flip() +
  labs() +
  theme(plot.title = element_text(size = 18, hjust = 0.6),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15))
ldabuv

lda.clado(uvlda)

# PART TWO UV AFTER COMPARISON -----------------------------------------

auvsamps = allsamps %>%
  dplyr::filter(sample %in% c("TY.1_D0_NT", "TY.2_D6_Dry", 
                              "TY.19_D0_DNA_UV", "TY.20_D6_DNA_UV")) %>% 
  mutate(auvgroup = ifelse(uv == "UV DNA", "UV after DNA Extraction", "No UV")) 
sample_data(phy.tyler) = auvsamps %>% set.samp()

auvlda = phy.tyler %>% phy.collapse() %>%
  lda.effect(class = "auvgroup")

auvlda.plot = auvlda %>% dplyr::filter(pass) %>%
  mutate(auvlda = if_else(as.numeric(factor(direction)) == 2, -lda, lda)) %>%
  arrange(auvlda) %>%
  mutate(taxon = fct_inorder(taxon))

ldaauv = ggplot(auvlda.plot, aes(x = taxon, y = auvlda, fill = direction)) +
  geom_col() + coord_flip() +
  labs() +
  theme(plot.title = element_text(size = 18, hjust = 0.6),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15))
ldaauv

lda.clado(auvlda)

# INVERSE SIMPSON DIVERSITY BAR GRAPH --------------------------------

diversityplot = ggplot(allsamps, aes(x = reorder(sample, InvSimpson), y = InvSimpson, fill = label)) + 
  geom_col() + 
  labs(x = "Sample", y = "Bacterial Diversity (inverse Simpson index)") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90))

diversityplot






















