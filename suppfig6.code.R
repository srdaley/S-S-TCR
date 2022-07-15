# Graphs and summaries in Supp Fig 6 are outputs of this script. 
# In RStudio, install and/or load "tidyverse" 
if (!require("tidyverse"))
  install.packages("tidyverse")
library(tidyverse)
# Set working directory to the folder that contains data
data_folder <- "" # **Instruction: Make "data_folder" the path to the folder containing the file called "fig3.clonotype.table.txt"
 data_folder <- "/Users/daleys6/OneDrive - Queensland University of Technology/Manuscripts/200801_dsb/NatComms/post_acceptance/code/" 
# The preceding line is the path to the folder on the author's computer
setwd (data_folder) # set the working directory to "data_folder"

# Load reference TCR dataset from WT Yae62beta-tg mice (PMID: 29726044)
pp <- read.delim(file = "suppfig6.clonotype.table.reference.txt", header = T)

# Count the number of times each TCR clonotype was detected in each lineage
# in order to determine a “Distribution in reference TCRa catalogs” for each TCR clonotype 
# LOI = Lineage Of Interest
LOI <- c("iel", "reg", "c4", "c8")
pp$lin <- factor("iel", levels = LOI)
pp$lin[pp$subset == "t_r" | pp$subset == "s_r"] <- "reg"
pp$lin[pp$subset == "t_4" | pp$subset == "s_4"] <- "c4"
pp$lin[pp$subset == "t_8" | pp$subset == "s_8"] <- "c8"

ppr <- pp %>% 
  group_by(v, cdr3aa, lin) %>%
  summarise(
    clones = n()
  )

# Determine the number of clones in each catalog
ppc <- ppr %>%
  group_by(lin) %>%
  summarise(
    catalog = n() 
  )


ppw <- ppr %>%
  pivot_wider(
    names_from = lin,
    values_from = clones,
    values_fill = 0
  ) 

ppw$iel_freq <- ppw$iel / ppc$catalog[ppc$lin == "iel"]
ppw$reg_freq <- ppw$reg / ppc$catalog[ppc$lin == "reg"]
ppw$c4_freq <- ppw$c4 / ppc$catalog[ppc$lin == "c4"]
ppw$c8_freq <- ppw$c8 / ppc$catalog[ppc$lin == "c8"]


ppw$total_freq <- ppw$iel_freq + ppw$reg_freq + ppw$c4_freq + ppw$c8_freq
ppw$iel_norm <- 100 * ppw$iel_freq / ppw$total_freq
ppw$reg_norm <- 100 * ppw$reg_freq / ppw$total_freq
ppw$c4_norm <- 100 * ppw$c4_freq / ppw$total_freq
ppw$c8_norm <- 100 * ppw$c8_freq / ppw$total_freq


ppl <- ppw %>%
  select(v, cdr3aa, iel_norm, reg_norm, c4_norm, c8_norm) %>%
  pivot_longer(
    names_to = "lineage",
    values_to = "rel_freq",
    cols = iel_norm:c8_norm
  )

# ppl contains a “Distribution in reference TCRa catalogs” for each TCR clonotype 

# Load the test datasets from Zap70-wt and Zap70-mutant Yae62beta-tg mice 
dt <- read.delim(file = "suppfig6.clonotype.table.test.txt", header = T)

# Plot overlaps at the TCR catalog level 
# Change the order of levels within BOI according to lineage being analyzed
# Set a colour scheme for plots
pal <- c("c8_norm" = "black", "c4_norm" = "royalblue1", 
         "reg_norm" = "green", "iel_norm" = "firebrick1")
BOI <- c("iel_norm", "c4_norm","reg_norm",  "c8_norm")
BOI_plot <- rev(BOI)
SGOI <- "wt"
SSOI <- "t_1"

# Get a list of clonotypes present in SGOI and SSOI 
# and attach info about distribution in published dataset
dt1 <- dt %>%
  filter(genotype == SGOI, subset == SSOI) %>% 
  distinct(v, cdr3aa) %>%
  merge(ppl)
dt1$lineage <- factor(dt1$lineage, levels = BOI)

################################################################################
# establish order for plotting along x-axis                                    #  
################################################################################

dt2 <- dt1 %>%
  arrange(v, cdr3aa, desc(rel_freq)) %>%
  group_by(v, cdr3aa) %>%
  summarise(
    lineage = first(lineage))

dt3 <- merge(dt1,dt2)
dt3$lineage <- factor(dt3$lineage, levels = BOI)
dt4 <- dt3 %>% filter(lineage == BOI[1]) %>% arrange(desc(rel_freq))
dt5 <- dt3 %>% filter(!(lineage == BOI[1])) %>% arrange(lineage,rel_freq)
dt6 <- bind_rows(dt4,dt5)
dt6$x <- row.names(dt6)
dt6$x <- as.numeric(dt6$x)

# Attach dt6$x to dt1
dt7 <- merge(dt1, dt6 %>% select(v, cdr3aa, x)) 
dt7$lineage <- factor(dt7$lineage, levels = BOI_plot)

# Plot it
pdf(file = paste("suppfig6",SGOI, SSOI,"pdf", sep = "."), height = 3, width = 4)
p <- ggplot(dt7, aes(x = x, y = rel_freq, color = lineage, fill = lineage))+
  geom_col(aes(fill = lineage)) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(title = paste (SGOI, SSOI, sep = "_"),
       x = "Individual TCR",
       Y = "Relative frequency (%)") +
  scale_x_continuous(
    limits = c(0, length(BOI) + nrow(dt7)/length(BOI)),
    breaks = c(0, nrow(dt7)/(2 * length(BOI)), nrow(dt7)/length(BOI))) + 
  annotate("text",x = nrow(dt7)/(5 * length(BOI)), y = 20, label = paste(nrow(dt7)/length(BOI)), colour = "white")
p
dev.off()
# This plot has labels on the axes and a legend.


pdf(file = paste("suppfig6",SGOI, SSOI,"clean.pdf", sep = "."), height = 1, width = 1)
p <- ggplot(dt7, aes(x = x, y = rel_freq, color = lineage, fill = lineage))+
  geom_col(aes(fill = lineage)) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  scale_x_continuous(
    limits = c(0, length(BOI) + nrow(dt7)/length(BOI)),
    breaks = c(0, nrow(dt7)/(2 * length(BOI)), nrow(dt7)/length(BOI))) + 
  annotate("text",x = nrow(dt7)/(4 * length(BOI)), y = 20, 
           label = paste(nrow(dt7)/length(BOI)), colour = "white", size = 3) +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_line(size = .2),
    axis.text = element_blank(),
    legend.position = 'none',
    panel.grid = element_blank()
  ) 
p
dev.off()
# This plot lacks labels on the axes and a legend.
# To make the figure for the paper, the axes and legend were annotated manually using Adobe Illustrator.    

# For statistical analysis, calculate area represented by each "lineage".
# Get a list of clonotypes present in each specific sample
# and attach info about distribution in published TCR dataset from WT Yae62beta-tg mice (PMID: 29726044)
dt$sample_iden <- paste(dt$mouse, dt$subset, sep = "_")
ESOI <- unique(dt$sample_iden)
invdl_c8 <- numeric(length = length(ESOI))
invdl_c4 <- numeric(length = length(ESOI))
invdl_iel <- numeric(length = length(ESOI))
invdl_reg <- numeric(length = length(ESOI))
overlapping_clones <- numeric(length = length(ESOI))
ot <- tibble(ESOI, invdl_iel, invdl_c4, invdl_c8, invdl_reg, overlapping_clones)

for(sn in 1:length(ESOI)){
  
  ef1 <- dt %>%
    filter(sample_iden == ESOI[sn]) %>% 
    distinct(v, cdr3aa) %>%
    merge(ppl)
  
  # establish order for plotting along x-axis called "conv_bias"
  ef2 <- ef1 %>%
    pivot_wider(
      names_from = lineage,
      values_from = rel_freq
    )
  ef2$conv_bias <- ef2$c4_norm / (ef2$c4_norm + ef2$c8_norm + ef2$iel_norm + ef2$reg_norm)
  ef2 <- arrange(ef2, conv_bias)
  ef2$x <- rownames(ef2)
  ef3 <- select(ef2, v, cdr3aa, x)
  ef3$x <- as.integer(ef3$x)
  # attach "conv_bias" column to ef1
  ef4 <- merge(ef1, ef3) %>%
    arrange(x)
  ef4$lineage <- factor(ef4$lineage, levels = BOI)
  ef5 <- ef4 %>%
    pivot_wider(names_from = lineage,
                values_from = rel_freq)
  ot$invdl_iel[sn] <-  sum(ef5$iel_norm) / nrow(ef5)
  ot$invdl_c4[sn] <-  sum(ef5$c4_norm) / nrow(ef5)
  ot$invdl_c8[sn] <- sum(ef5$c8_norm) / nrow(ef5)
  ot$invdl_reg[sn] <- sum(ef5$reg_norm) / nrow(ef5)
  ot$overlapping_clones[sn] <- nrow(ef1)/4
  }

# attach info
info <- dt %>% 
  group_by(sample_iden, genotype, subset) %>%
  summarise(
    clones = n()
  )
out <- merge(ot, info, by.x = "ESOI", by.y = "sample_iden")

out$subset <- factor(out$subset, levels = c("t_1", "t_4", "t_8", "t_r", "g8a", "s_4", "s_8", "s_r"))
out$genotype <- factor(out$genotype, levels = c("wt", "zap"))
out <- arrange(out, subset, genotype)
out1 <- out[, c(1,7,8,2:6,9)]
write.table(out1, file = "suppfig6.sourcedata.txt", sep = "\t", row.names = F, quote = F)
