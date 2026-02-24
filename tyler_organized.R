# LIBRARIES ------------------------
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
# DATA IMPORT ---------------
rm(list = ls())
setwd("/Users/Tyler/Documents/Research/IntegrityProject/tyler_yang_micro_preservation")  # change path name if working in a different directory

load("data/phy.tyler.RData")
load("data/phy.tyler.additional.RData")

# Reset the samp
phy.tyler = phy.tyler %>% 
  select(otu,Superkingdom,Phylum,Class,Order,Family,Genus,Species)

s.tyler = phy.tyler %>% get.samp() %>%
  left_join(trace_tbl.tyler,by="sample") %>%
  mutate(temp=ifelse(sample %in% c("1A","1B"),"-80C",temp),
         temp=factor(temp,levels=c("-80C", "-20C", "4C", "room temp")),
         time=factor(time,levels=c("day 0", "day 3", "day 6", "day 8", "day 9", "day 11")),
         treatment=factor(treatment,levels=c("none", "75C", "UV", "75C+UV", "autoclave", "autoclave+UV", "UV DNA")),
         sample_number=order(order(treatment,temp,time)),
         letter=str_extract(sample,"(?<=^[0-9]{1,2})[AB]"),
         sample2=paste(treatment,temp,storage,time,sample_number,sep="|"),
         sample2=fct_reordern(sample2,sample_number),
         sample.comparator=ifelse(experiment==1,"1A","TY.1_D0_NT"))
sample_data(phy.tyler) <- s.tyler %>% set.samp()

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


# Samp and OTU dataframes 
allsamps = create_samp_df(phy.tyler)  # get SAMP dataframe
otusub = create_otu_df(phy.tyler, allsamps)  # get OTU dataframe (with mapping/samp data)

firstsamps = allsamps %>%
  dplyr::filter(experiment == 1) %>%
  mutate(replicate = if_else(grepl("A", sample), 1, 2))  # for Figure 1 pairing
secondsamps = allsamps %>%
  dplyr::filter(experiment == 2)

firstotu = otusub %>% 
  dplyr::filter(experiment == 1) %>%
  mutate(replicate = if_else(grepl("A", sample), 1, 2)) 
secondotu = otusub %>%
  dplyr::filter(experiment == 2)
# First Stacked Overview ---------------
first_overview = ggplot() +
  geom_taxonomy(data = firstotu, aes(x = sampname, y = absseqs, fill = Species), # create columns
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

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/overview_1.pdf', width = 13, height = 7)
first_overview
dev.off()

# First PCAs ---------------
phy.others = read_rds("/Users/Tyler/Documents/Research/IntegrityProject/other.samps.rds")
others = get.samp(phy.others)

phy.tyler.first <- phy.tyler %>% subset_samples(experiment==1) %>% prune_unused_taxa()

taxa_names(phy.tyler.first) = as.character(refseq(phy.tyler.first)) %>% unname()  # ASK ABOUT THIS
taxa_names(phy.others) = as.character(refseq(phy.others)) %>% unname()

phy.together = merge_phyloseq(otu_table(phy.tyler.first), otu_table(phy.others),
                              tax_table(phy.tyler.first),tax_table(phy.others))
phy.together_test = merge_phyloseq(phy.tyler.first, phy.others)

# MEAN BRAY
bray_dist = calc.distance(phy.together, "bray", mean = TRUE)
# UNFOLD HORN
horn_dist = calc.distance(phy.together_test, "horn", unfold = TRUE)

bray_ord = ordinate(phy.together, distance = bray_dist, method = "NMDS")
bray_data = bray_ord$points %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(allsamps, by = "sample") %>%
  mutate(origin = case_when(
    is.na(experiment) ~ "Other",
    experiment == 1 ~ "Our Study"),
    temp = if_else(origin == "Other", "none", temp)) %>%
  mutate(origin = fct_relevel(origin, "Our Study", "Other"))
first_bray = ggplot(bray_data, aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = temp, shape = origin)) +
  scale_color_manual(
    values = c("-80C" = "steelblue", 
               "-20C" = "cadetblue", 
               "4C" = "cadetblue3", 
               "room temp" = "indianred",
               "none" = "gray3"),
    breaks = c("-80C", "-20C", "4C", "room temp")) +
  # theme(legend.position = "none") +
  labs(title = "Bray (mean)", color = "Storage Temperature", shape = "Origin")
first_bray

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/bray_mean_1.pdf', width = 7, height = 6)
first_bray
dev.off()

horn_ord = ordinate(phy.together, distance = horn_dist, method = "NMDS")
horn_data = horn_ord$points %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(allsamps, by = "sample") %>%
  mutate(origin = case_when(
    is.na(experiment) ~ "Other",
    experiment == 1 ~ "Our Study"),
    temp = if_else(origin == "Other", "none", temp)) %>%
  mutate(origin = fct_relevel(origin, "Our Study", "Other"))

first_horn = ggplot(horn_data, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = temp, shape = origin)) +
  scale_color_manual(
    values = c("-80C" = "steelblue", 
               "-20C" = "cadetblue", 
               "4C" = "cadetblue3", 
               "room temp" = "indianred",
               "none" = "gray3"),
    breaks = c("-80C", "-20C", "4C", "room temp")) +
  # theme(legend.position = "none") +
  labs(title = "Horn (unfolded)", color = "Storage Temperature", shape = "Origin")
first_horn

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/horn_unfolded_1.pdf', width = 7, height = 6)
first_horn
dev.off()


first_bray_no_legend = first_bray +
  guides(color = "none", shape = "none")
pcoas = ggarrange(first_bray_no_legend, first_horn, nrow = 1, widths = c(3, 4)) %>% 
  annotate_figure(top = text_grob("PCOAs", size = 16))
pcoas

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/pcoas_1.pdf', width = 10, height = 4)
pcoas
dev.off()

# First Alpha Diversity Bar Graph ---------------

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

alpha = ggarrange(simpson_plot, shannon_plot, nrow = 2, heights = c(1, 1))
alpha

pdf('/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/alpha_1.pdf', width = 8, height = 6)
alpha
dev.off()














# First LEfSe ---------------------

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

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Test Graphs/ldas_temp_1.pdf", width = 23, height = 12)
ldas
dev.off()

# Compare Uniques ----------------

compare.uniques <- function(sample1,sample2,phy=phy.tyler) {
  physub <- prune_samples(c(sample1,sample2),phy)
  ranks <- rank_names(physub)
  otusub <- physub %>% get.otu.melt() %>%
    group_by(otu) %>%
    mutate(unique=n()==1) %>%
    ungroup() %>% 
    mutate(otu=fct_reordern(otu,!!!syms(ranks)),
           col=as.numeric(otu),
           sample = factor(sample, levels = unique(c(!!sample1, !!sample2)))) 
  
  ggplot(otusub,aes(x=col,y=pctseqs)) + 
    geom_col(aes(alpha=unique,fill=otu)) +
    geom_point(data=filter(otusub,unique),aes(x=col,y=pctseqs),color="green",size=0.5) +
    scale_fill_taxonomy(data=otusub,fill=otu) +
    # scale_color_manual(values=c("TRUE"="red","FALSE"=NA)) +
    scale_alpha_manual(values=c("TRUE"=1,"FALSE"=0.4)) +
    scale_y_continuous(trans=log_epsilon_trans(0.001)) +
    facet_grid(sample ~ .)
}
compare.uniques.new = function(group1, group2, group1_label, group2_label, phy = phy.tyler) {
  all_targets = c(group1, group2)
  physub = prune_samples(all_targets, phy)
  ranks = rank_names(physub)
  
  otusub = physub %>% 
    get.otu.melt() %>%
    mutate(group_label = case_when(
      sample %in% group1 ~ group1_label,
      sample %in% group2 ~ group2_label)) %>%
    group_by(group_label, otu, across(all_of(ranks))) %>% 
    summarise(pctseqs = mean(pctseqs, na.rm = TRUE), .groups = "drop") %>%
    group_by(otu) %>%
    mutate(unique = n() == 1) %>% 
    ungroup() %>% 
    mutate(otu = fct_reordern(otu, !!!syms(ranks)),
           col = as.numeric(otu),
           group_label = factor(group_label, levels = c(group1_label, group2_label)))
  ggplot(otusub, aes(x = col, y = pctseqs)) + 
    geom_col(aes(alpha = unique, fill = otu)) +
    geom_point(data = filter(otusub, unique), aes(x = col, y = pctseqs), color = "green", size = 0.5) +
    scale_fill_taxonomy(data = otusub, fill = otu) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
    scale_y_continuous(trans = log_epsilon_trans(0.001)) +
    # Change facet to use the new group_label
    facet_grid(group_label ~ .)
}

standards = c("1A", "1B")
first80 = compare.uniques.new(standards, c("13A.-80.D11", "13B.-80.D11"), "Standards", "Ultra-low Freezer") +
  guides(alpha = "none", fill = "none")
first20 = compare.uniques.new(standards, c("12A.-20.D11", "12B.-20.D11"), "Standards", "Freezer") +
  guides(alpha = "none", fill = "none")
first4 = compare.uniques.new(standards, c("11A.4C.D11", "11B.4C.D11"), "Standards", "Refrigerator") +
  guides(alpha = "none", fill = "none")
firstRT = compare.uniques.new(standards, c("10A.RT.D11", "10B.RT.D11"), "Standards", "Room Temperature") +
  guides(alpha = "none", fill = "none")
firstuniques = ggarrange(first80, first20, first4, firstRT) %>%
  annotate_figure(top = text_grob("Compare Uniques", size = 20))
firstuniques

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/uniques_1.pdf", width = 15, height = 10)
firstuniques
dev.off()

first80_cherry = compare.uniques("1A","13A.-80.D11") +
  guides(alpha = "none", fill = "none")
first20_cherry = compare.uniques("1A","12A.-20.D11") +
  guides(alpha = "none", fill = "none")
first4_cherry = compare.uniques("1A","11A.4C.D11") +
  guides(alpha = "none", fill = "none")
firstRT_cherry = compare.uniques("1A","10A.RT.D11") +
  guides(alpha = "none", fill = "none")
firstuniques_cherry = ggarrange(first80_cherry, first20_cherry, first4_cherry, firstRT_cherry) %>%
  annotate_figure(top = text_grob("Compare Uniques (cherrypicked replicate)", size = 20))
firstuniques_cherry

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/uniques_cherry_1.pdf", width = 15, height = 10) 
firstuniques_cherry
dev.off()

# Depict ----------------------------

depict_ty <- function(sample1, sample2, phy=phy.tyler) {
  
  subsamps <- c(sample1, sample2)
  physub <- phy %>% 
    filter(sample %in% subsamps)
  
  # 1. Get the base data
  otu_raw <- physub %>% get.otu.melt()
  
  # 2. Create the explicit sorting order based ONLY on sample1
  # We group by OTU and grab the abundance specifically where sample == sample1.
  # If an OTU exists in sample2 but is 0 in sample1, this ensures it gets sorted as 0.
  otu_order <- otu_raw %>%
    group_by(otu) %>%
    summarize(sort_val = sum(pctseqs[sample == sample1], na.rm = TRUE)) %>%
    arrange(sort_val) %>% # Use desc(sort_val) if you want high-to-low
    pull(otu)
  
  # 3. Apply this order to the main dataframe
  otu <- otu_raw %>%
    mutate(first = sample == sample1,
           y = pctseqs,
           # Convert OTU to a factor with the specific levels we calculated above
           otu_ordered = factor(otu, levels = otu_order),
           # Convert to numeric for the x-axis mapping
           x = as.numeric(otu_ordered))
  
  # Prepare the layers
  # We extend the step line slightly to the left/right for aesthetics
  otu1 <- otu %>% 
    filter(first) %>%
    bind_rows(tibble(x = range(.$x, na.rm = TRUE) + c(-1, 1), y = c(0, 0)))
  
  otu2 <- otu %>% filter(!first)
  
  # Plot
  ggplot() +
    geom_col(data = otu2, aes(x = x, y = y, fill = otu)) +
    geom_step(data = otu1, aes(x = x, y = y), direction = "mid") +
    scale_fill_taxonomy(data = otu, fill = otu) +
    scale_y_continuous(trans = log_epsilon_trans(0.001)) +
    theme()
}

depict3("1A","1B")
depict3("1A","11B.4C.D11")

# Depict Full ------------

phy1 <- phy.tyler %>% filter(experiment==1) %>%
  mutate(baseline=sample=="1A",
         templabel="Temp",
         timelabel="Time")


otu1base <- phy1 %>% 
  filter(baseline,prune_unused_taxa=FALSE) %>%
  get.otu.melt(filter.zero=FALSE) %>%
  transmute(otu,pctseqs0=pctseqs)

otu1 <- phy1 %>% 
  get.otu.melt(filter.zero=FALSE) %>%
  left_join(otu1base,by="otu") %>%
  filter(pctseqs>0|pctseqs0>0) %>%
  group_by(sample) %>%
  arrange(desc(pctseqs0),desc(pctseqs)) %>%
  mutate(col=row_number(),
         extra=pctseqs0==0 & pctseqs>0) %>%
  ungroup()
s1 <- phy1 %>% get.samp()

g1.asv <- ggplot() +
  geom_col(data=otu1,aes(x=col,y=pctseqs,fill=otu)) +
  geom_step(data=otu1,aes(x=col,y=pctseqs0),direction="mid") +
  # geom_step(data=otu1,aes(x=col,y=pctseqs0),direction="mid") +
  geom_text(data=s1,aes(x=Inf,y=Inf,label=str_glue("{sample}\n{short_number(nseqs)} seqs")),hjust=1,vjust=1,color="blue") +
  # geom_text(data=s2,aes(x=Inf,y=Inf,label=str_glue("{sample}\n{short_number(nseqs)} seqs")),hjust=1,vjust=1,color="blue") +
  geom_rect(data=filter(otu1,baseline),
            aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf),
            fill=NA,color="blue",linetype="longdash") +
  geom_bracket(data=filter(otu1,extra),
               aes(x=col,y=ave(pctseqs,sample,FUN=max),
                   fontsize=3,label="unique\nASVs"),tip="square") + 
  scale_fill_taxonomy(data=otu1,fill=otu) +
  scale_y_continuous(trans=log_epsilon_trans(0.001)) +
  facet_nested(templabel+temp+letter~timelabel+time) + 
  labs(title = "Depict Overview")
g1.asv

pdf("/Users/Tyler/Documents/Research/IntegrityProject/Organized Graphs/depict_overview_1.pdf", width = 13, height = 18) 
g1.asv
dev.off()





