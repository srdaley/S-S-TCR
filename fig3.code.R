# Graphs in Fig 3 are outputs of this script. 
# In RStudio, install and/or load "tidyverse" 
if (!require("tidyverse"))
  install.packages("tidyverse")
library(tidyverse)
# Set working directory to the folder that contains data
data_folder <- "" # **Instruction: Make "data_folder" the path to the folder containing the file called "fig3.clonotype.table.txt"
 data_folder <- "/Users/daleys6/OneDrive - Queensland University of Technology/Manuscripts/200801_dsb/NatComms/post_acceptance/code/" 
# The preceding line is the path to the folder on the author's computer
setwd (data_folder) # set the working directory to "data_folder"

# Load the data called "fig3.clonotype.table.txt"
df <- read.delim(file = "fig3.clonotype.table.txt", header = T)

# Set factors and levels of factors
COI <- c("a", "b") # COI = Chains Of Interest
GOI <- c("wt", "dko", "zap") # GOI = Genotypes Of Interest
SOI <- c("t_p", "g8a", "s_4", "s_8") # SOI = Subsets Of Interest

# Exclude non-functional clones and clones encountered only once or twice 
min_count <- 3
df$nf <- str_count(df$cdr3aa, "[^[:upper:]]")
dfy <- df %>% 
  filter(count >= min_count, nf == 0) 
# Some samples were from the same T cell subset of the same Zap70-mutant mouse,
# either because multiple aliquots of T cells were sampled independently (thymic samples)
# or because the same cDNA sample was subjected to PCR amplification more than once (2 SI CD8aa IEL samples). 
# Each result must be a unique combo of mouse, genotype, subset, chain
dfx <- dfy %>% 
  group_by(mouse, genotype, subset, chain, v, cdr3nt, cdr3aa) %>% 
  summarise(instances = n())

###################################################################
# Check cysteine usage within 2 positions of CDR3 apex
# As Trbv1 sequences can have a germline-encoded Cys at CDR3 position 2, 
# which is within 2 positions of the apex of CDR3 sequences < 8 amino acids long, 
# exclude Trbv1 sequences with a CDR3 length < 8 amino acids. 
VOI <- c("TRBV1", "TRBV1*01")
dfx$cdr3aa <- as.character(dfx$cdr3aa)
dfx$cdr3_trimmed <- substr(dfx$cdr3aa, 2, nchar(dfx$cdr3aa)-1)
dfx$len <- as.integer(nchar(dfx$cdr3_trimmed))
dfs <- filter(dfx, !(v %in% VOI & len < 8))
dfs["n2"] <- as.character("-")
dfs["n1"] <- as.character("-")
dfs["top"] <- as.character("-")
dfs["c1"] <- as.character("-")
dfs["c2"] <- as.character("-")
dfs["top"] <- ifelse(dfs$len>0,(substr(dfs$cdr3_trimmed, floor((nchar(dfs$cdr3_trimmed)/2)+1),floor((nchar(dfs$cdr3_trimmed)/2)+1))),"-")
dfs["c1"] <- ifelse(dfs$len>2,(substr(dfs$cdr3_trimmed, floor((nchar(dfs$cdr3_trimmed)/2)+2),floor((nchar(dfs$cdr3_trimmed)/2)+2))),"-")
dfs["n1"] <- ifelse(dfs$len>1,(substr(dfs$cdr3_trimmed, floor((nchar(dfs$cdr3_trimmed)/2)-0),floor((nchar(dfs$cdr3_trimmed)/2)-0))),"-")
dfs["c2"] <- ifelse(dfs$len>4,(substr(dfs$cdr3_trimmed, floor((nchar(dfs$cdr3_trimmed)/2)+3),floor((nchar(dfs$cdr3_trimmed)/2)+3))),"-")
dfs["n2"] <- ifelse(dfs$len>3,(substr(dfs$cdr3_trimmed, floor((nchar(dfs$cdr3_trimmed)/2)-1),floor((nchar(dfs$cdr3_trimmed)/2)-1))),"-")
dfs$cencys <- dfs$n2 == "C" | dfs$n1 == "C" | dfs$top == "C" |dfs$c1 == "C" | dfs$c2 == "C"
###################################################################
dfo <- dfs %>%
  group_by(genotype, mouse, chain, subset, cencys) %>%
  summarise(clones = n()) %>%
  pivot_wider(names_from = cencys, values_from = clones, values_fill = 0) %>%
  mutate(
    sequences = `FALSE` + `TRUE`,
    freq = `TRUE` / sequences
  )

dfo$freq_corr <- dfo$freq + (dfo$`TRUE`==0)/(dfo$sequences)  
dfo$log_freq <- log10(dfo$freq_corr)
dfo$lod <- dfo$`TRUE` == 0

dfo <- rename(dfo, Cys_F = 'FALSE')
dfo <- rename(dfo, Cys_T = 'TRUE')

d <- dfo %>% 
  summarise(
    Cys_F = sum(Cys_F),
    Cys_T = sum(Cys_T)
  ) %>%
  mutate(
    sequences = Cys_F + Cys_T,
    freq = Cys_T / sequences
  )

d$freq_corr <- d$freq + (d$Cys_T==0)/(d$sequences)  
d$log_freq <- log10(d$freq_corr)
d$lod <- d$Cys_T == 0


# plot it for figure 1

pal <- c("t_p" = "gray50", "g8a" = "firebrick1",
         "s_4" = "royalblue1", "s_8" = "black")

# plot Cysteine index results by subset
ht = .9
unique(d$genotype) # Result of this line: "wt"  "dko" "zap" 
# SGOI = Special Genotype Of Interest 
SGOI <- "wt" # change SGOI to "wt",  "dko", or "zap" to change the genotype being plotted
dfp <- d %>% filter(genotype == SGOI)
dfp$subset <- factor(dfp$subset, levels = SOI)
dfp$chain <- factor(dfp$chain, levels = COI)
dfp$xanc <- match(dfp$subset, SOI)
shi <- c(-.2, .2)
dfp$x <- dfp$xanc + shi[dfp$chain]

# get group means
dfm <- dfp %>% 
  group_by(subset, chain) %>%
  summarise(
    log_freq = mean(log_freq)
  )

dfm$subset <- factor(dfm$subset, levels = SOI)
dfm$chain <- factor(dfm$chain, levels = COI)
dfm$xanc <- match(dfm$subset, SOI)
dfm$x <- dfm$xanc + shi[dfm$chain]


pdf(file = paste("fig3", SGOI, "pdf", sep = "."),
    height = 3, width = 4)
g <- ggplot(dfp, aes(x=x, y= log_freq, color = subset, shape = chain)) +
  geom_jitter(height = 0, width = .1) +
  scale_color_manual(values = pal) +
  scale_shape_manual(values = c(21,4)) +
  theme_bw() +
  scale_y_continuous(limits = c(-4,-.5),
                     breaks = seq(-4,-1,1)) +
  scale_x_continuous(breaks = seq(1, length(SOI)),
                     labels = SOI) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(title = SGOI) +
  geom_segment(data = dfm, aes(x = x - .1, xend = x + .1, y = log_freq, yend = log_freq))
g
dev.off()
# This plot has labels on the axes and a legend.
###############################################################################

pdf(file = paste("fig3", SGOI,"clean.pdf", sep = "."), height = ht, width =  ht,
    useDingbats = FALSE)
g <- ggplot(dfp, aes(x=x, y= log_freq, color = subset, shape = chain)) +
  geom_jitter(width = .1, size = .8, stroke = .3, height = 0) +
  scale_color_manual(values = pal) +
  scale_shape_manual(values = c(21,4)) +
  theme_bw() +
  scale_y_continuous(limits = c(-4,-.1),
                     breaks = seq(-4,-1,1)) +
  scale_x_continuous(breaks = seq(1, length(SOI)),
                     labels = SOI) +
  theme(panel.grid.major.x  = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none',
        axis.ticks.y = element_line(size = .2),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(size = .2)) +
  geom_segment(data = dfm, aes(x = x - .2, xend = x + .2, y = log_freq, yend = log_freq), 
               alpha = .6, size = .5)
g
dev.off()  
# This plot lacks labels on the axes and a legend.
# To make the figure for the paper, the axes and legend were annotated manually using Adobe Illustrator.    

# Write file for export to Prism for statistical analyses.
dp <- d 
dp$genotype <- factor(dp$genotype, levels = GOI)
dp$subset <- factor(dp$subset, levels = SOI)
dp <- arrange(dp, genotype, chain, subset)

#write.table(dp, file = paste(pre, doa, "prism.txt", sep = "."), sep = "\t", row.names = F, quote = F)  

# Write file for export to the Source Data file to accompany the online version of the paper
ed <- dp %>%
  select(genotype,  subset, mouse, chain, Cys_T, Cys_F, freq)
NOI <- c("Pre-sel'n", "CD8aa IEL", "CD4 Tconv", "CD8 Tconv")
ed$T_cell_subset <- NOI[match(ed$subset, SOI)] 
ed$subset <- NULL
ed$Cys_index <- ed$freq * 100
ed$freq <- NULL
colnames(ed) <- c( "Genotype", "MouseID", "TCR_chain", "Cys_present", 
                   "Cys_absent", "T_cell_subset",  "Cys_index")
ed1 <- ed[,c(1,2,6,3:5,7)]
write.table(ed1, file = "fig3.sourcedata.txt", sep = "\t", row.names = F, quote = F)


