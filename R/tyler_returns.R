# SET UP -------
# install.packages(c("openssl", "systemfonts", "gdtools", "ggiraph"))
# install.packages("remotes")
# remotes::install_github("ying14/yingtools2")
# install.packages("RcppArmadillo", type = "source")
# install.packages("ade4", type = "source")
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
install.packages('ggtree', repos = c('https://bioc.r-universe.dev', 'https://cloud.r-project.org'))
library(phyloseq)
library(yingtools2)
library(phytools)
library(tidyverse)
library(ggtree)
library(lattice)
library(grid)
library(gridExtra)
library(ggfittext)
library(ggpubr)
library(ggnewscale)
library(vegan)
library(readxl)
library(glue)

# IMPORT DATA ------
rm(list = ls())
setwd("/Users/Tyler/Documents/Research/IntegrityProject")  # change path name if working in a different directory
load("phy.tyler.RData")

# Functions that really just condense ugly code for dataframe cleaning
add_clean_names = function(phy_samp_df) {
  cleaned = phy_samp_df %>%
    mutate(sampname = sample,
           samplabel = sample) %>%
    # Manual renaming of special controls
    mutate(sampname = recode(sampname,
                             "1A" = "Control 1",
                             "1B" = "Control 2",
                             "TY.1_D0_NT" = "Control 3",
                             "TY.19_D0_DNA_UV" = "D0 UV post-extraction",
                             "TY.20_D6_DNA_UV" = "D6 UV post-extraction",
                             "TY.21_D9_DNA_UV" = "D9 UV post-extraction",
                             "TY.8_D6_AC" = "D6 autoclave",
                             "TY.9_D9_AC" = "D9 autoclave",
                             "TY.7_D0_AC" = "D0 autoclave",
                             "TY.16_D0_AC_UV" = "D0 autoclave + UV pre-ex",
                             "TY.17_D6_AC_UV" = "D6 autoclave + UV pre-ex",
                             "TY.18_D9_AC_UV" = "D9 autoclave + UV pre-ex",
                             "TY.2_D6_Dry" = "D6 dried",
                             "TY.3_D9_Dry" = "D9 dried",
                             "TY.4_D0_75C" = "D0 75C heat",
                             "TY.5_D6_75C" = "D6 75C heat",
                             "TY.6_D9_75C" = "D9 75C heat",
                             "TY.10_D0_UV" = "D0 UV pre-extraction",
                             "TY.11_D6_UV" = "D6 UV pre-extraction",
                             "TY.12_D9_UV" = "D9 UV pre-extraction",
                             "TY.13_D0_75C_UV" = "D0 75C + UV pre-ex",
                             "TY.14_D6_75C_UV" = "D6 75C + UV pre-ex",
                             "TY.15_D9_75C_UV" = "D9 75C + UV pre-ex",
                             "2A.RT.D3" = "D3 at 20C A",
                             "2B.RT.D3" = "D3 at 20C B",
                             "3A.4C.D3" = "D3 at 4C A",
                             "3B.4C.D3" = "D3 at 4C B",
                             "4A.-20.D3" = "D3 at -20C A",
                             "4B.-20.D3" = "D3 at -20C B",
                             "5A.-80.D3" = "D3 at -80C A",
                             "5B.-80.D3" = "D3 at -80C B",
                             "6A.RT.D8" = "D8 at 20C A",
                             "6B.RT.D8" = "D8 at 20C B",
                             "7A.4C.D8" = "D8 at 4C A",
                             "7B.4C.D8" = "D8 at 4C B",
                             "8A.-20.D8" = "D8 at -20C A",
                             "8B.-20.D8" = "D8 at -20C B",
                             "9A.-80.D8" = "D8 at -80C A",
                             "9B.-80.D8" = "D8 at -80C B",
                             "10A.RT.D11" = "D11 at 20C A",
                             "10B.RT.D11" = "D11 at 20C B",
                             "11A.4C.D11" = "D11 at 4C A",
                             "11B.4C.D11" = "D11 at 4C B",
                             "12A.-20.D11" = "D11 at -20C A",
                             "12B.-20.D11" = "D11 at -20C B",
                             "13A.-80.D11" = "D11 at -80C A",
                             "13B.-80.D11" = "D11 at -80C B"
    ),
    samplabel = recode(samplabel,
                       "1A" = "T1",
                       "1B" = "T2",
                       "TY.1_D0_NT" = "TY2",
                       "TY.19_D0_DNA_UV" = "TY19",
                       "TY.20_D6_DNA_UV" = "TY20",
                       "TY.21_D9_DNA_UV" = "TY21",
                       "TY.8_D6_AC" = "TY8",
                       "TY.9_D9_AC" = "TY9",
                       "TY.7_D0_AC" = "TY7",
                       "TY.16_D0_AC_UV" = "TY16",
                       "TY.17_D6_AC_UV" = "TY17",
                       "TY.18_D9_AC_UV" = "TY18",
                       "TY.2_D6_Dry" = "TY2",
                       "TY.3_D9_Dry" = "TY3",
                       "TY.4_D0_75C" = "TY4",
                       "TY.5_D6_75C" = "TY5",
                       "TY.6_D9_75C" = "TY6",
                       "TY.10_D0_UV" = "TY10",
                       "TY.11_D6_UV" = "TY11",
                       "TY.12_D9_UV" = "TY12",
                       "TY.13_D0_75C_UV" = "TY13",
                       "TY.14_D6_75C_UV" = "TY14",
                       "TY.15_D9_75C_UV" = "TY15",
                       "2A.RT.D3" = "T3",
                       "2B.RT.D3" = "T4",
                       "3A.4C.D3" = "T5",
                       "3B.4C.D3" = "T6",
                       "4A.-20.D3" = "T7",
                       "4B.-20.D3" = "T8",
                       "5A.-80.D3" = "T9",
                       "5B.-80.D3" = "T10",
                       "6A.RT.D8" = "T11",
                       "6B.RT.D8" = "T12",
                       "7A.4C.D8" = "T13",
                       "7B.4C.D8" = "T14",
                       "8A.-20.D8" = "T15",
                       "8B.-20.D8" = "T16",
                       "9A.-80.D8" = "T17",
                       "9B.-80.D8" = "T18",
                       "10A.RT.D11" = "T19",
                       "10B.RT.D11" = "T20",
                       "11A.4C.D11" = "T21",
                       "11B.4C.D11" = "T22",
                       "12A.-20.D11" = "T23",
                       "12B.-20.D11" = "T24",
                       "13A.-80.D11" = "T25",
                       "13B.-80.D11" = "T26"))
  return(cleaned)
}
create_samp_df = function(phy_object) {
  qPCR = read.csv("/Users/Tyler/Documents/Research/IntegrityProject/TY-qPCR.csv", stringsAsFactors = FALSE) %>%
    filter(Dilution.Number != "#N/A") %>%
    select(Sample.ID:Average.CT.Value) %>%
    slice(1:44) %>%
    rename('samplabel' = 'Sample.ID', 
           'seqs_total' = 'X16S.Adjusted') %>%
    mutate(seqs_total = as.numeric(seqs_total))
  
  allsamps = phy_object %>%
    get.samp(stats = TRUE) %>%
    add_clean_names() %>%
    mutate(
      label = if_else(experiment < 1.5, "one", "two"),
      sample = fct_relevel(sample,
                           "1A", "1B", "2A.RT.D3", "2B.RT.D3", "3A.4C.D3", "3B.4C.D3",
                           "4A.-20.D3", "4B.-20.D3", "5A.-80.D3", "5B.-80.D3", "6A.RT.D8",
                           "6B.RT.D8", "7A.4C.D8", "7B.4C.D8", "8A.-20.D8", "8B.-20.D8",
                           "9A.-80.D8", "9B.-80.D8", "10A.RT.D11", "10B.RT.D11", "11A.4C.D11",
                           "11B.4C.D11", "12A.-20.D11", "12B.-20.D11", "13A.-80.D11", "13B.-80.D11",
                           "TY.1_D0_NT", "TY.2_D6_Dry", "TY.3_D9_Dry", "TY.4_D0_75C", "TY.5_D6_75C",
                           "TY.6_D9_75C", "TY.10_D0_UV", "TY.11_D6_UV", "TY.12_D9_UV", "TY.13_D0_75C_UV",
                           "TY.14_D6_75C_UV", "TY.15_D9_75C_UV", "TY.7_D0_AC", "TY.8_D6_AC", "TY.9_D9_AC",
                           "TY.16_D0_AC_UV", "TY.17_D6_AC_UV", "TY.18_D9_AC_UV", "TY.19_D0_DNA_UV",
                           "TY.20_D6_DNA_UV", "TY.21_D9_DNA_UV"),
      sampname = fct_relevel(sampname,
                             "Control 1", "Control 2", "D3 at 20C A", "D3 at 20C B", "D3 at 4C A", "D3 at 4C B",
                             "D3 at -20C A", "D3 at -20C B", "D3 at -80C A", "D3 at -80C B", "D8 at 20C A",
                             "D8 at 20C B", "D8 at 4C A", "D8 at 4C B", "D8 at -20C A", "D8 at -20C B",
                             "D8 at -80C A", "D8 at -80C B", "D11 at 20C A", "D11 at 20C B", "D11 at 4C A",
                             "D11 at 4C B", "D11 at -20C A", "D11 at -20C B", "D11 at -80C A", "D11 at -80C B",
                             "Control 3", "D6 dried", "D9 dried", "D0 UV pre-extraction", "D6 UV pre-extraction", "D9 UV pre-extraction",
                             "D0 75C heat", "D6 75C heat", "D9 75C heat", "D0 75C + UV pre-ex", "D6 75C + UV pre-ex", "D9 75C + UV pre-ex", 
                             "D0 autoclave", "D6 autoclave", "D9 autoclave", "D0 autoclave + UV pre-ex", "D6 autoclave + UV pre-ex",
                             "D9 autoclave + UV pre-ex", "D0 UV post-extraction", "D6 UV post-extraction", "D9 UV post-extraction"),
      heat = if_else(experiment == 1, "no heat", heat),
      treatment = if_else(experiment == 1, "none", treatment),
      origin = case_when(
        is.na(experiment) ~ "Other",
        experiment == 1 ~ "Year 1 Samples",
        heat == "autoclave" ~ "Autoclave",
        uv == "UV DNA" ~ "UV DNA",
        TRUE ~ "Year 2 Samples")) %>%
    left_join(qPCR, by = "samplabel") 
  return(allsamps)
}
create_otu_df = function(phy_object, phy_samp_df) {
  physub = phy_object %>% 
    phy.collapse()  # combine OTUs of the same taxonomic classification
  otusub = physub %>% 
    get.otu.melt(sample_data = FALSE) %>%
    tax.plot(data = TRUE) %>%  
    left_join(phy_samp_df) %>%
    mutate(absseqs = pctseqs * seqs_total)
}

allsamps = create_samp_df(phy.tyler)  # get SAMP dataframe
otusub = create_otu_df(phy.tyler, allsamps)  # get OTU dataframe (with mapping/samp data)

firstsamps = allsamps %>%
  dplyr::filter(experiment == 1) %>%
  mutate(replicate = if_else(grepl("A", sample), 1, 2))  # for Figure 1 pairing
secondsamps = allsamps %>%
  dplyr::filter(experiment == 2)

firstotusub = otusub %>% 
  dplyr::filter(experiment == 1) %>%
  mutate(replicate = if_else(grepl("A", sample), 1, 2)) 
secondotusub = otusub %>%
  dplyr::filter(experiment == 2)

# Figure 1A Stacked Overview ----------------
first_overview = ggplot() +
  geom_taxonomy(data = firstotusub, aes(x = sampname, y = absseqs, fill = Species), # create columns
                show.legend = TRUE, width = 0.8) +
  scale_y_continuous(trans = log_epsilon_trans(10000)) +
  theme(panel.spacing.x = unit(2, "lines")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(), axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0.5, size = 15), axis.title.x = element_text(size = 24), # shows sample name on x-axis
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_text(size = 24),
        plot.title = element_text(hjust = 0.5, size = 35)) +
  geom_text(data = firstsamps, aes(x = sampname, y = 1.05, label = short_number(nseqs)),angle = -90) +
  labs(x = "", y = "Microbial Composition", title = "", fill = "Taxa")
first_overview

# Figure 1B PCAs ----------------

phy.others = read_rds("/Users/Tyler/Documents/Research/IntegrityProject/other.samps.rds")
others = get.samp(phy.others)

phy.tyler.first <- phy.tyler %>% subset_samples(experiment==1) %>% prune_unused_taxa()

taxa_names(phy.tyler.first) = as.character(refseq(phy.tyler.first)) %>% unname()  # ASK ABOUT THIS
taxa_names(phy.others) = as.character(refseq(phy.others)) %>% unname()

phy.together = merge_phyloseq(otu_table(phy.tyler.first), otu_table(phy.others),
                              tax_table(phy.tyler.first),tax_table(phy.others))

bray_dist = calc.distance(phy.together, "pct.bray")
horn_dist = calc.distance(phy.together, "horn", mean = TRUE)

bray_ord = ordinate(phy.together, distance = bray_dist, method = "NMDS")
bray_data = bray_ord$points %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(allsamps, by = "sample") %>%
  mutate(origin = case_when(
    is.na(experiment) ~ "Other",
    experiment == 1 ~ "Our Study")) %>%
  mutate(origin = fct_relevel(origin, "Our Study", "Other"))
bray_g1 = ggplot(bray_data, aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = origin)) +
  scale_color_manual(
    values = c("Our Study" = "steelblue", 
               "Other" = "gray3")) +
  # theme(legend.position = "none") +
  labs(title = "Bray (pct)", color = "Origin")
bray_g1

horn_ord = ordinate(phy.together, distance = horn_dist, method = "NMDS")
horn_data = horn_ord$points %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(allsamps, by = "sample") %>%
  mutate(origin = case_when(
    is.na(experiment) ~ "Other",
    experiment == 1 ~ "Our Study")) %>%
  mutate(origin = fct_relevel(origin, "Our Study", "Other"))

horn_g1 = ggplot(horn_data, aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = origin)) +
  scale_color_manual(
    values = c("Our Study" = "steelblue", 
               "Other" = "gray3")) +
  labs(title = "Horn (mean)", color = "Origin") 
horn_g1

pcoas = ggarrange(bray_g2, horn_g2, nrow = 1, widths = c(3, 4)) 
  # %>% annotate_figure(top = text_grob("PCOAs", size = 10))
pcoas

# Figure 1C Global Diversity Measures ---------------

first_divsamps = firstsamps %>%
  mutate(time = fct_relevel(time, "day 0", "day 3", "day 8", "day 11"),
         sampname = fct_relevel(sampname, "Control 1", "Control 2", 
                                "D3 at -80C A", "D3 at -80C B", "D3 at -20C A", "D3 at -20C B", "D3 at 4C A", "D3 at 4C B", "D3 at 20C A", "D3 at 20C B",
                                "D8 at -80C A", "D8 at -80C B", "D8 at -20C A", "D8 at -20C B", "D8 at 4C A", "D8 at 4C B", "D8 at 20C A", "D8 at 20C B",
                                "D11 at -80C A", "D11 at -80C B", "D11 at -20C A", "D11 at -20C B", "D11 at 4C A", "D11 at 4C B", "D11 at 20C A", "D11 at 20C B"),
         temp = case_when(time == "day 0" ~ "control",
                          TRUE ~ temp)) 
simpson_plot = ggplot(first_divsamps, aes(x = sampname, y = InvSimpson, fill = fct_relevel(temp, "-80C", "-20C", "4C", "room temp"))) + 
  geom_col() + 
  facet_grid(~ time, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = c("-80C" = "steelblue", 
               "-20C" = "cadetblue", 
               "4C" = "cadetblue3", 
               "room temp" = "indianred")) +
  theme(legend.position = "right", 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_text(size = 15),
        axis.text.y = element_text(), axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.3, size = 20),
        legend.title = element_text(size = 13), legend.text = element_text(size = 12)) +
  labs(x = "", y = "Inverse Simpson Index", fill = "Storage Temperature") +
  ylim(0, 40)
simpson_plot

shannon_plot = ggplot(first_divsamps, aes(x = sampname, y = Shannon, fill = fct_relevel(temp, "-80C", "-20C", "4C", "room temp"))) + 
  geom_col() + 
  facet_grid(~ time, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = c("-80C" = "steelblue", 
               "-20C" = "cadetblue", 
               "4C" = "cadetblue3", 
               "room temp" = "indianred")) +
  theme(legend.position = "right", 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_text(size = 15),
        axis.text.y = element_text(), axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.3, size = 20),
        legend.title = element_text(size = 13), legend.text = element_text(size = 12)) +
  labs(x = "", y = "Shannon Index", fill = "Storage Temperature") +
  ylim(0, 9)
shannon_plot



invsamps = firstsamps %>%
  mutate(time = fct_relevel(time, "day 0", "day 3", "day 8", "day 11")) 
invsamps$temp[9] = "control"
invsamps$temp[10] = "control"
invsamps = invsamps %>% 
  mutate(temp = fct_relevel(temp, "control", "-80C", "-20C", "4C", "room temp"))

diversity = ggplot(invsamps, aes(x = replicate, y = Shannon, fill = "black")) + geom_col() + 
  facet_grid(temp~time) +
  scale_fill_manual(values = "skyblue3") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), strip.text.x = element_text(size = 9),
        # axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_text(vjust = 0.5, size = 14), 
        strip.text.y = element_text(size = 9),
        plot.tag = element_text(angle = -90, size = 14), plot.tag.position = c(1.025, 0.5),
        plot.title = element_text(hjust = 0.5, size = 15), plot.subtitle = element_text(hjust = 0.5, size = 14),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 20, l = 5, b = 10)) +
  labs(y = "Bacterial Diversity (invSimpson Index)", subtitle = "Time Group", 
       tag = "Temperature Group") 

diversity

# Figure 1 Combining ------------------

pcoas = ggarrange(bray_g1, horn_g1, nrow = 2)
diversities = ggarrange(simpson_plot, shannon_plot, nrow = 2)
bottom = ggarrange(pcoas, diversities, nrow = 1, widths = c(2, 5))
figure1 = ggarrange(first_overview, bottom, nrow = 2, heights = c(8, 6)) %>% 
  annotate_figure(top = text_grob("Figure 1", size = 30))

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/figure1_legend.pdf", width = 20, height = 16)
figure1
dev.off()


bottom = ggarrange(pcoas, diversity, nrow = 1, widths = c(2, 1))
total = ggarrange(first_overview, bottom, nrow = 2, heights = c(8, 7)) %>% 
  annotate_figure(top = text_grob("Figure 1", size = 30))
total

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/figure1.pdf", width = 20, height = 16)
total
dev.off()

# Figure 2 Species Analysis ----------------------------

phy.lefse = phy.tyler %>% 
  phy.collapse(taxranks=c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
lefsesamps = firstsamps %>%
  mutate(timegroup = ifelse(time == "day 0", "early", ifelse(time == "day 3", "early", "late"))) %>%
  mutate(tempgroup = ifelse(temp == "-80C", "Ultra-low temperature storage", "Other storage")) 
sample_data(phy.lefse) = lefsesamps %>% set.samp()

lda = phy.lefse %>%
  lda.effect(class = "tempgroup", subclass = "time")

lda.plot = lda %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(as.numeric(factor(direction)) == 2, -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Species")  # probably should ask about this

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

# LEfSe with modified phyloseq

phy.lefse = phy.tyler %>% 
  phy.collapse(taxranks=c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
lefsesamps = firstsamps %>%
  mutate(timegroup = ifelse(time == "day 0", "early", ifelse(time == "day 3", "early", "late"))) %>%
  mutate(tempgroup = ifelse(temp == "-80C", "Ultra-low temperature storage", "Other storage")) 
sample_data(phy.lefse) = lefsesamps %>% set.samp()

phy.lefse.20 = phy.lefse %>%
  filter(temp %in% c("-80C", "-20C"))
lda_20 = phy.lefse.20 %>%
  lda.effect(class = "temp", subclass = "time")
lda_plot_20 = lda_20 %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "-80C", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Species")  # probably should ask about this
ldatemp_20 = ggplot(lda_plot_20, aes(x = taxon, y = lda, fill = direction)) +
  geom_col() + coord_flip() +
  labs() +
  theme(plot.title = element_text(size = 18, hjust = 0.6),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) +
  scale_fill_manual(values = c("-80C" = "steelblue", "-20C" = "cadetblue"),
                    labels = c("-20C" = "Freezer (-20C)", "-80C" = "Ultra-low (-80C)")) +
  labs(y = "LDA Index", x = "Significantly Different Taxa",
       fill = "Temperature Group",
       title = "-80C vs -20C")
ldatemp_20

phy.lefse.4 = phy.lefse %>%
  filter(temp %in% c("-80C", "4C"))
lda_4 = phy.lefse.4 %>%
  lda.effect(class = "temp", subclass = "time")
lda_plot_4 = lda_4 %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "-80C", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Species")  # probably should ask about this
ldatemp_4 = ggplot(lda_plot_4, aes(x = taxon, y = lda, fill = factor(direction, levels = c("4C", "-80C")))) +
  geom_col() + coord_flip() +
  labs() +
  theme(plot.title = element_text(size = 18, hjust = 0.6),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) +
  scale_fill_manual(values = c("-80C" = "steelblue", "4C" = "cadetblue3"),
                    labels = c("4C" = "Refrigerator (4C)", "-80C" = "Ultra-low (-80C)")) +
  labs(y = "LDA Index", x = "Significantly Different Taxa",
       fill = "Temperature Group",
       title = "-80C vs 4C")
ldatemp_4

phy.lefse.RT = phy.lefse %>%
  filter(temp %in% c("-80C", "room temp"))
lda_RT = phy.lefse.RT %>%
  lda.effect(class = "temp", subclass = "time")
lda_plot_RT = lda_RT %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "-80C", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Species")  # probably should ask about this
ldatemp_RT = ggplot(lda_plot_RT, aes(x = taxon, y = lda, fill = factor(direction, levels = c("room temp", "-80C")))) +
  geom_col() + coord_flip() +
  labs() +
  theme(plot.title = element_text(size = 18, hjust = 0.6),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) +
  scale_fill_manual(values = c("-80C" = "steelblue", "room temp" = "indianred"),
                    labels = c("room temp" = "Room Temperature", "-80C" = "Ultra-low (-80C)")) +
  labs(y = "LDA Index", x = "Significantly Different Taxa",
       fill = "Temperature Group",
       title = "-80C vs RT")
ldatemp_RT

ldas = ggarrange(ldatemp_20, ldatemp_4, ldatemp_RT, nrow = 1)
ldas 

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/ldas.pdf", width = 23, height = 12)
ldas
dev.off()

# lefse tables
lda_plot_20 %>% dt
lda_plot_4 %>% dt
lda_plot_RT %>% dt

pretty_freezer_lda = lda_plot_20 %>%
  mutate(test = "Freezer Storage") %>%
  select(test, taxon, direction, lda, kw.pvalue)
pretty_freezer_lda %>% dt

pretty_refrigerator_lda = lda_plot_4 %>%
  mutate(test = "Refrigerator Storage") %>%
  select(test, taxon, direction, lda, kw.pvalue)
pretty_refrigerator_lda %>% dt

pretty_RT_lda = lda_plot_RT %>%
  mutate(test = "Room Temp Storage") %>%
  select(test, taxon, direction, lda, kw.pvalue)
pretty_RT_lda %>% dt

pretty_ldas = rbind(pretty_freezer_lda, pretty_refrigerator_lda, pretty_RT_lda) %>%
  arrange(as.character(taxon))
pretty_ldas %>% dt

write.csv(pretty_freezer_lda, '/Users/Tyler/Documents/Research/IntegrityProject/tyler_yang_micro_preservation/lda/pretty ldas/freezer_lda.csv')
write.csv(pretty_refrigerator_lda, '/Users/Tyler/Documents/Research/IntegrityProject/tyler_yang_micro_preservation/lda/pretty ldas/refrigerator_lda.csv')
write.csv(pretty_RT_lda, '/Users/Tyler/Documents/Research/IntegrityProject/tyler_yang_micro_preservation/lda/pretty ldas/roomtemp_lda.csv')
write.csv(pretty_ldas, '/Users/Tyler/Documents/Research/IntegrityProject/tyler_yang_micro_preservation/lda/pretty ldas/pretty_ldas.csv')

# GENUS
lda_plot_20 = lda_20 %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "-80C", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Genus")  # probably should ask about this
pretty_freezer_lda = lda_plot_20 %>%
  mutate(test = "Freezer Storage") %>%
  select(test, taxon, direction, lda, kw.pvalue)
pretty_freezer_lda %>% dt

lda_plot_4 = lda_4 %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "-80C", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Genus") 
pretty_refrigerator_lda = lda_plot_4 %>%
  mutate(test = "Refrigerator Storage") %>%
  select(test, taxon, direction, lda, kw.pvalue)
pretty_refrigerator_lda %>% dt

lda_plot_RT = lda_RT %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "-80C", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Genus")
pretty_RT_lda = lda_plot_RT %>%
  mutate(test = "Room Temp Storage") %>%
  select(test, taxon, direction, lda, kw.pvalue)
pretty_RT_lda %>% dt

pretty_ldas = rbind(pretty_freezer_lda, pretty_refrigerator_lda, pretty_RT_lda) %>%
  arrange(as.character(taxon))
pretty_ldas %>% dt

write.csv(pretty_freezer_lda, '/Users/Tyler/Documents/Research/IntegrityProject/tyler_yang_micro_preservation/lda/pretty ldas/freezer_genus_lda.csv')
write.csv(pretty_refrigerator_lda, '/Users/Tyler/Documents/Research/IntegrityProject/tyler_yang_micro_preservation/lda/pretty ldas/refrigerator_genus_lda.csv')
write.csv(pretty_RT_lda, '/Users/Tyler/Documents/Research/IntegrityProject/tyler_yang_micro_preservation/lda/pretty ldas/roomtemp_genus_lda.csv')
write.csv(pretty_ldas, '/Users/Tyler/Documents/Research/IntegrityProject/tyler_yang_micro_preservation/lda/pretty ldas/pretty_genus_ldas.csv')


# tax comparison

subsamps_gold = list("1A", "1B")
physub_gold = phy.tyler %>%
  filter(sample %in% c("1A", "1B"))
otusub_gold = get.otu.melt(physub_gold, filter.zero = FALSE) %>%
  mutate(sign = ifelse(sample == subsamps_gold[1],1,-1),
         y = sign * pctseqs) %>%
  group_by(otu) %>%
  mutate(y1 = pctseqs[which(sample == subsamps_gold[1])[1]],
         y2 = pctseqs[which(sample == subsamps_gold[2])[1]]) %>%
  ungroup() %>%
  arrange(y1, -y2) %>%
  mutate(x = fct_inorder(otu))

compare_gold = ggplot(otusub_gold, aes(x = x, y = y, fill = otu)) +
  geom_col() + 
  scale_fill_taxonomy(data = otusub_gold, fill = otu) + 
  scale_y_continuous(trans = log_epsilon_trans(0.001)) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(title = "Comparing Gold Standards")
compare_gold

subsamps_RT = list("1A", "10B.RT.D11")
physub_RT = phy.tyler %>%
  filter(sample %in% c("1A", "10B.RT.D11"))
otusub_RT = get.otu.melt(physub_RT, filter.zero = FALSE) %>%
  mutate(sign = ifelse(sample == subsamps_RT[1],1,-1),
         y = sign * pctseqs) %>%
  group_by(otu) %>%
  mutate(y1 = pctseqs[which(sample == subsamps_RT[1])[1]],
         y2 = pctseqs[which(sample == subsamps_RT[2])[1]]) %>%
  ungroup() %>%
  arrange(y1, -y2) %>%
  mutate(x = fct_inorder(otu)) 

compare_RT = ggplot(otusub_RT, aes(x = x, y = y, fill = otu)) +
  geom_col() + 
  scale_fill_taxonomy(data = otusub_RT, fill = otu) + 
  scale_y_continuous(trans = log_epsilon_trans(0.001)) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(title = "Comparing Standard with Room Temp Day 11")
compare_RT

create_compare = function(sample1, sample2, title = TRUE) {
  subsamps = list(sample1, sample2)
  physub = phy.tyler %>%
    filter(sample %in% subsamps)
  otusub = get.otu.melt(physub, filter.zero = FALSE) %>%
    mutate(sign = ifelse(sample == subsamps[1],1,-1),
           y = sign * pctseqs) %>%
    group_by(otu) %>%
    mutate(y1 = pctseqs[which(sample == subsamps[1])[1]],
           y2 = pctseqs[which(sample == subsamps[2])[1]]) %>%
    ungroup() %>%
    arrange(y1, -y2) %>%
    mutate(x = fct_inorder(otu)) 
  if (title == TRUE) {
    this_title = paste0("Comparing ", sample1, " with ", sample2)
  } else {
    this_title = ""
  }
  
  compare = ggplot(otusub, aes(x = x, y = y, fill = otu)) +
    geom_col() + 
    scale_fill_taxonomy(data = otusub, fill = otu) + 
    scale_y_continuous(trans = log_epsilon_trans(0.001)) +
    theme(legend.position = "none") +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(title = this_title)
  return(compare)
}

compare_gold = create_compare("13A.-80.D11", "13B.-80.D11", title = FALSE)
compare_gold
compare_20 = create_compare("13A.-80.D11", "12B.-20.D11", title = FALSE)
compare_20
compare_4 = create_compare("13A.-80.D11", "11B.4C.D11", title = FALSE)
compare_4
compare_RT = create_compare("13A.-80.D11", "10B.RT.D11", title = FALSE)
compare_RT

compare_meta = ggarrange(
  compare_gold, compare_20, compare_4, compare_RT,
  ncol = 2,
  nrow = 2,                  # 4 plots → 2×2 grid
  labels = c("Gold", "-20°C", "4°C", "RT"),  # optional
  common.legend = TRUE,      # optional: share one legend
  legend = "none"           # position of common legend
)
compare_meta

figure2 = ggarrange(compare_meta, ldas, nrow = 2, heights = c(2, 3)) %>% 
  annotate_figure(top = text_grob("Figure 2", size = 30))
pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/figure2.pdf", width = 30, height = 20)
figure2
dev.off()

# Clado --------
clado_20 = lda.clado(lda_20) +
  scale_fill_manual(values = c("-80C" = "steelblue", "-20C" = "cadetblue"),
                    labels = c("-20C" = "Freezer (-20C)", "-80C" = "Ultra-low (-80C)")) +
  labs(fill = "Storage Preference") +
  guides(size = "none")
clado_4 = lda.clado(lda_4) +
  scale_fill_manual(values = c("-80C" = "steelblue", "4C" = "cadetblue3"),
                    labels = c("4C" = "Refrigerator (4C)", "-80C" = "Ultra-low (-80C)")) +
  labs(fill = "Storage Preference") +
  guides(size = "none")
clado_RT = lda.clado(lda_RT) +
  scale_fill_manual(values = c("-80C" = "steelblue", "room temp" = "indianred"),
                    labels = c("room temp" = "Room Temperature", "-80C" = "Ultra-low (-80C)")) +
  labs(fill = "Storage Preference") +
  guides(size = "none")

lda_all = phy.lefse %>%
  lda.effect(class = "temp", subclass = "time")
clado_all = lda.clado(lda_all) + # class = "temp", subclass = "time
  scale_fill_manual(
    values = c("-80C" = "steelblue",
               "-20C" = "cadetblue",
               "4C" = "cadetblue3",
               "room temp" = "indianred"))

clados = ggarrange(clado_20, clado_4, clado_RT, nrow = 1)
clados

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/clados.pdf", width = 21, height = 7)
clados
dev.off()

 pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/clado_20.pdf", width = 10, height = 10)
clado_20
dev.off()

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/clado_4.pdf", width = 10, height = 10)
clado_4
dev.off()

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/clado_RT.pdf", width = 10, height = 10)
clado_RT
dev.off()

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/clado_all.pdf", width = 10, height = 10)
clado_all
dev.off()

# Clado Combined --------------

lda_diff_RT = lda_RT %>% 
  filter(info == "PASS")
lda_diff_20 = lda_20 %>% 
  filter(info == "PASS")
lda_diff_4 = lda_4 %>%
  filter(info == "PASS")
lda_diff_RT %>% dt
lda_diff_20 %>% dt
lda_diff_4 %>% dt

lda_alldiffs = rbind(lda_diff_RT, lda_diff_20, lda_diff_4) %>%
  mutate(
    direction = factor(direction,
                       levels = c("-20C", "4C", "room temp", "-80C")
    )
  ) %>%
  group_by(taxon) %>%
  mutate(multiple = n() > 1) %>%
  slice_min(order_by = direction, n = 1, with_ties = FALSE) 

taxa_to_remove = lda_alldiffs$taxon  # or unique(df1$taxon)

lda_all_filtered = lda_all %>%
  mutate(multiple = FALSE) %>%
  filter(!taxon %in% taxa_to_remove) %>%
  rbind(lda_alldiffs) %>%
  mutate(direction = factor(direction, levels = c("-80C", "-20C", "4C", "room temp")))

lda_all_filtered %>% dt


clado_combined = lda.clado(lda_all_filtered) + # class = "temp", subclass = "time
  scale_fill_manual(
    values = c("-80C" = "steelblue",
               "-20C" = "cadetblue",
               "4C" = "cadetblue3",
               "room temp" = "indianred")) +
  guides(size = "none") +
  labs(fill = "Storage Temp Preference")
clado_combined 

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/clado_combined.pdf", width = 10, height = 7)
clado_combined
dev.off()

lda_enrichment = lda_alldiffs %>%
  ungroup() %>%
  filter(direction != "-80C")
taxa_to_remove_enrich = lda_enrichment$taxon
lda_all_enrich = lda_all %>%
  mutate(multiple = FALSE) %>%
  filter(!taxon %in% taxa_to_remove_enrich) %>%
  rbind(lda_enrichment) %>%
  mutate(direction = factor(direction, levels = c("-80C", "-20C", "4C", "room temp")))

lda_all_enrich %>% dt


clado_enrich = lda.clado(lda_all_enrich) + # class = "temp", subclass = "time
  scale_fill_manual(
    values = c("-80C" = "steelblue",
               "-20C" = "steelblue",
               "4C" = "cadetblue3",
               "room temp" = "indianred")) +
  guides(size = "none") +
  labs(fill = "Storage Temp Enrichment")
clado_enrich 

lda_diff_RT_decr = lda_RT %>% 
  filter(info == "PASS") %>%
  filter(direction == "-80C") %>%
  select(-direction) %>%
  mutate(direction = "room temp")
lda_diff_20_decr = lda_20 %>% 
  filter(info == "PASS") %>%
  filter(direction == "-80C") %>%
  select(-direction) %>%
  mutate(direction = "-20C")
lda_diff_4_decr = lda_4 %>%
  filter(info == "PASS") %>%
  filter(direction == "-80C") %>%
  select(-direction) %>%
  mutate(direction = "4C")

lda_diff_decr_all = rbind(lda_diff_RT_decr, lda_diff_20_decr, lda_diff_4_decr) %>%
  mutate(
    direction = factor(direction, levels = c("-20C", "4C", "room temp", "-80C"))) %>%
  group_by(taxon) %>%
  mutate(multiple = n() > 1) %>%
  slice_min(order_by = direction, n = 1, with_ties = FALSE) 
taxa_to_remove_decr = lda_diff_decr_all$taxon

lda_all_decr = lda_all %>%
  mutate(multiple = FALSE) %>%
  filter(!taxon %in% taxa_to_remove_decr) %>%
  rbind(lda_diff_decr_all) %>%
  mutate(direction = factor(direction, levels = c("-80C", "-20C", "4C", "room temp")))

lda_all_decr %>% dt
clado_decr = lda.clado(lda_all_decr) + # class = "temp", subclass = "time
  scale_fill_manual(
    values = c("-80C" = "steelblue",
               "-20C" = "steelblue",
               "4C" = "cadetblue3",
               "room temp" = "indianred")) +
  guides(size = "none") +
  labs(fill = "Storage Temp Aversion")
clado_decr 
clado_enrich

clado_enrich_decr = ggarrange(clado_enrich, clado_decr, nrow = 1)
clado_enrich_decr

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/clado_enrich.pdf", width = 10, height = 7)
clado_enrich
dev.off()

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/clado_decr.pdf", width = 10, height = 7)
clado_decr
dev.off()

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/clado_preference.pdf", width = 20, height = 7)
clado_enrich_decr
dev.off()
 # Figure 3 Part Two ----------------
two_otusub = secondotusub %>%
  mutate(seqs_total = if_else(treatment == "UV DNA", 3000, seqs_total),  # this is arbitrary
         absseqs = pctseqs * seqs_total)

second_overview = ggplot() +
  geom_taxonomy(data = two_otusub, aes(x = sampname, y = absseqs, fill = Species), # create columns
                show.legend = FALSE, width = 0.8) +
  scale_y_continuous(trans = log_epsilon_trans(10000)) +
  theme(panel.spacing.x = unit(2, "lines")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(), axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0.5, size = 15), axis.title.x = element_text(size = 24), # shows sample name on x-axis
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_text(size = 24),
        plot.title = element_text(hjust = 0.5, size = 35)) +
  geom_text(data = secondsamps, aes(x = sampname, y = 1.05, label = short_number(nseqs)),angle = -90) +
  labs(x = "Sample", y = "Microbial Composition", title = "", fill = "Taxa")
second_overview

phy.others = read_rds("/Users/Tyler/Documents/Research/IntegrityProject/other.samps.rds")
others = get.samp(phy.others)

phy.tyler.second <- phy.tyler %>% subset_samples(experiment==2) %>% prune_unused_taxa()

taxa_names(phy.tyler.second) = as.character(refseq(phy.tyler.second)) %>% unname()  # ASK ABOUT THIS
taxa_names(phy.others) = as.character(refseq(phy.others)) %>% unname()

phy.together.second = merge_phyloseq(otu_table(phy.tyler.second), otu_table(phy.others),
                              tax_table(phy.tyler.second),tax_table(phy.others))

bray_dist = calc.distance(phy.together.second, "pct.bray")
horn_dist = calc.distance(phy.together.second, "horn", mean = TRUE)

bray_ord2 = ordinate(phy.together.second, distance = bray_dist, method = "NMDS")
bray_data2 = bray_ord2$points %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(allsamps, by = "sample") %>%
  mutate(origin = case_when(
    is.na(experiment) ~ "Other Samples",
    treatment == "autoclave" ~ "Autoclave", 
    treatment == "autoclave+UV" ~ "Autoclave + UV",
    treatment == "UV DNA" ~ "UV DNA", 
    TRUE ~ "Other Treatments")) %>%
  mutate(origin = fct_relevel(origin, "UV DNA", "Autoclave + UV", "Autoclave", "Other Treatments", "Other Samples"))
bray_g2 = ggplot(bray_data2, aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = origin)) +
  scale_color_manual(
    values = c("UV DNA" = "purple3",
               "Autoclave + UV" = "pink2",
               "Autoclave" = "indianred2",
               "Other Treatments" = "steelblue",
               "Other Samples" = "gray3")) +
  theme(legend.position = "none") +
  labs(title = "Bray (pct)", color = "Treatment")
bray_g2

horn_ord2 = ordinate(phy.together.second, distance = horn_dist, method = "NMDS")
horn_data2 = horn_ord2$points %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(allsamps, by = "sample") %>%
  mutate(origin = case_when(
    is.na(experiment) ~ "Other Samples",
    treatment == "autoclave" ~ "Autoclave", 
    treatment == "autoclave+UV" ~ "Autoclave + UV",
    treatment == "UV DNA" ~ "UV DNA", 
    TRUE ~ "Other Treatments")) %>%
  mutate(origin = fct_relevel(origin, "UV DNA", "Autoclave + UV", "Autoclave", "Other Treatments", "Other Samples"))

horn_g2 = ggplot(horn_data2, aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = origin)) +
  scale_color_manual(
    values = c("UV DNA" = "purple3",
               "Autoclave + UV" = "pink2",
               "Autoclave" = "indianred2",
               "Other Treatments" = "steelblue",
               "Other Samples" = "gray3")) +
  theme(legend.position = "right") +
  labs(title = "Horn (mean)", color = "Origin") 
horn_g2

second_pcoas = ggarrange(bray_g2, horn_g2, nrow = 1, widths = c(3, 4)) 
# %>% annotate_figure(top = text_grob("PCOAs", size = 10))
second_pcoas

simsamps = secondsamps %>%
  mutate(heat = fct_relevel(heat, "no heat", "75C"))

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
  labs(x = "Storage Time", y = "Bacterial Diversity (invSimpson Index)", subtitle = "UV Treatment",
       tag = "Heat Treatment")
seconddiversity

bottom = ggarrange(second_pcoas, seconddiversity, nrow = 1, widths = c(2, 1))
total_two = ggarrange(second_overview, bottom, nrow = 2, heights = c(8, 4)) %>% 
  annotate_figure(top = text_grob("Figure 3", size = 30))

total_two

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/figure3.pdf", width = 20, height = 14)
total_two
dev.off()

# Figure 4 Proteobacteria Big Figure ----------------------------

select_samps = list("1A", "TY.1_D0_NT", "TY.21_D9_DNA_UV", "TY.19_D0_DNA_UV", 
                    "TY.20_D6_DNA_UV", "TY.9_D9_AC", "TY.8_D6_AC", 
                    "TY.17_D6_AC_UV", "TY.16_D0_AC_UV", "TY.18_D9_AC_UV", "TY.7_D0_AC")
sotutrimmed = otusub %>%
  dplyr::filter(sample %in% select_samps) %>%
  mutate(sample = fct_relevel(sample, "1A", "TY.1_D0_NT", "TY.7_D0_AC", "TY.8_D6_AC",
                              "TY.9_D9_AC", "TY.16_D0_AC_UV", "TY.17_D6_AC_UV", "TY.18_D9_AC_UV",
                              "TY.19_D0_DNA_UV", "TY.20_D6_DNA_UV", "TY.21_D9_DNA_UV"))

taxbar = sotutrimmed %>%
  mutate(Phylum = fct_lump_n(Phylum, 4, other_level = "Other")) %>%
  ggplot() +
  geom_col(aes(x = str_glue("{Species}"),
               y = pctseqs, fill = Species)) +
  scale_y_continuous(trans = log_epsilon_trans(0.001)) +
  scale_fill_taxonomy(data=otusub,fill=Species, guide=guide_taxonomy(ncol=1)) +
  facet_grid(sample~Phylum, scales = "free_x", space = "free_x") +
  theme(legend.position = "none", axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 0, size = 20),
        strip.text.y = element_text(angle = 0, size = 20),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title = element_text(size = 40, hjust = 0.5),
        axis.title = element_text(size = 35),
        plot.caption = element_text(size = 25)) +
  labs(x = "Species", y = "Sequence Abundance")

taxbar

pie_samps = list("TY.1_D0_NT", "TY.2_D6_Dry", "TY.3_D9_Dry",
                 "TY.7_D0_AC", "TY.8_D6_AC", "TY.9_D9_AC",
                 "TY.16_D0_AC_UV", "TY.17_D6_AC_UV", "TY.18_D9_AC_UV",
                 "TY.19_D0_DNA_UV", "TY.20_D6_DNA_UV", "TY.21_D9_DNA_UV")

pie_df = otusub %>%
  filter(sample %in% select_samps) %>%
  mutate(sample = fct_relevel(sample, "1A", "TY.1_D0_NT", "TY.7_D0_AC", "TY.8_D6_AC",
                              "TY.9_D9_AC", "TY.16_D0_AC_UV", "TY.17_D6_AC_UV", "TY.18_D9_AC_UV",
                              "TY.19_D0_DNA_UV", "TY.20_D6_DNA_UV", "TY.21_D9_DNA_UV"),
         pie_size = case_when(
           sample %in% c("1A", "TY.1_D0_NT") ~ 20,
           sample %in% c("TY.7_D0_AC") ~ 12,
           sample == "TY.8_D6_AC" ~ 15,
           sample == "TY.9_D9_AC" ~ 11,
           sample == "TY.16_D0_AC_UV" ~ 15,
           sample == "TY.17_D6_AC_UV" ~ 17,
           sample == "TY.18_D9_AC_UV" ~ 11,
           sample == "TY.19_D0_DNA_UV" ~ 8,
           sample == "TY.20_D6_DNA_UV" ~ 6,
           sample == "TY.21_D9_DNA_UV" ~ 4
         ))

pie_df = otusub %>%
  filter(sample %in% select_samps) %>%
  mutate(sample = fct_relevel(sample, "1A", "TY.1_D0_NT", "TY.7_D0_AC", "TY.8_D6_AC",
                              "TY.9_D9_AC", "TY.16_D0_AC_UV", "TY.17_D6_AC_UV", "TY.18_D9_AC_UV",
                              "TY.19_D0_DNA_UV", "TY.20_D6_DNA_UV", "TY.21_D9_DNA_UV"),
         pie_size = case_when(
           sample %in% c("1A", "TY.1_D0_NT") ~ 20,
           sample %in% c("TY.7_D0_AC") ~ 17,
           sample == "TY.8_D6_AC" ~ 19,
           sample == "TY.9_D9_AC" ~ 16,
           sample == "TY.16_D0_AC_UV" ~ 19,
           sample == "TY.17_D6_AC_UV" ~ 18,
           sample == "TY.18_D9_AC_UV" ~ 16,
           sample == "TY.19_D0_DNA_UV" ~ 15,
           sample == "TY.20_D6_DNA_UV" ~ 14,
           sample == "TY.21_D9_DNA_UV" ~ 13
         ))

red_pie = ggplot(pie_df, aes(x=pie_size / 2, y=pctseqs, fill=Species, width = pie_size)) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +
  facet_wrap(~ sample, ncol = 1) +
  scale_fill_taxonomy(data=otusub,fill=Species, guide=guide_taxonomy(ncol=1)) +
  theme_void() +
  theme(legend.position = "none", strip.text.x = element_blank()) +
  labs(fill = "Taxa")
red_pie

figure4 = ggarrange(taxbar, red_pie, nrow = 1, widths = c(5, 1))
pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/figure4.pdf", width = 20, height = 14)
figure4
dev.off()

# Figure 5 Obligate Anaerobes ----------------
target = c("Bacteroidetes", "Lachnospiraceae", "Oscillospiraceae", "Other Clostridia")

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
    Family %in% c("Lachnospiraceae", "Oscillospiraceae") ~ Family,
    Class == "Clostridia" ~ "Other Clostridia",
    Phylum %in% c("Bacteroidetes", "Actinobacteria") ~ Phylum, 
    TRUE ~ "Other Taxa")) %>%
  dplyr::filter(taxa %in% target)


figure5 = ggplot(loessotu, aes(x = day, y = pctseqs, group = taxa, fill = taxa, color = taxa)) +
  geom_jitter(width = 0.5, height = 0) +
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
       title = "Figure 5",
       subtitle = "Predominantly Anaerobic Taxonomic Groups")

figure5

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/figure5.pdf", width = 7, height = 9)
figure5
dev.off()

# Figure 6 Comparing Heat Block and Pre-UV ----------------------
phy.lefse.two = phy.tyler %>% 
  phy.collapse(taxranks=c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
sample_data(phy.lefse.two) = secondsamps %>% set.samp()

phy.lefse.heat = phy.lefse.two %>%
  filter(treatment %in% c("none", "UV", "75C", "75C+UV"))
lda_heat = phy.lefse.heat %>%
  lda.effect(class = "uv", subclass = "time")
lda_plot_heat = lda_heat %>% dplyr::filter(pass) %>%
  mutate(lda = if_else(direction == "UV", -lda, lda)) %>%
  arrange(lda) %>%
  mutate(taxon = fct_inorder(taxon)) %>%
  filter(taxrank == "Species")  # probably should ask about this
ldatemp_uv = ggplot(lda_plot_heat, aes(x = taxon, y = lda, fill = direction)) +
  geom_col() + 
  labs() +
  ylim(0, 4) + 
  theme(plot.title = element_text(size = 18, hjust = 0.6),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("UV" = "steelblue", "no UV" = "cadetblue"),
                    labels = c("UV" = "UV exposure pre-extraction", "no UV" = "No UV")) +
  labs(y = "LDA Index", x = "Significantly Different Taxa",
       fill = "Temperature Group",
       title = "UV vs No UV")
ldatemp_uv
pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/ldauv.pdf", width = 7, height = 9)
ldatemp_uv
dev.off()
