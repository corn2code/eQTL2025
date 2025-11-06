###############################################
setwd("/Users/vladimir/Downloads")
#library(readxl)
library(data.table)
#library(writexl)
#library(ggvenn)
library(tidyverse)
library(ggrastr)
library(patchwork)
library(ggpubr)
library(ggpmisc)
library(svglite)
library(cowplot) 
library(multcompView) # anova_res

theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
             plot.title = element_text(hjust = 0.5), plot.subtitle=element_text(hjust=0.5))


# load data
maize <- fread("totalcombined_snps_summary_NS_maize_693_V2.csv", data.table = F)
sorghum.rna <- fread("totalcombined_snps_summary_NS_sorghum_2_V2.csv", data.table = F)

### remove peaks with less than three snps, we have plotted the boxplot and checked how the look.
maize2 <- maize[which(maize$num_SNPs >= 3),]
sorghum2.rna <- sorghum.rna[which(sorghum.rna$num_SNPs >= 3),]

trans.maize<- maize2[maize2$cis_trans == "trans",]
trans.sorghum <- sorghum2.rna[sorghum2.rna$cis_trans == "trans",]

length(unique(trans.maize$traits))
length(unique(trans.sorghum$traits))

mean(maize2$pLength)
mean(sorghum2.rna$pLength)

# Create a horizontal boxplot without the x-axis
box1 <- ggplot(maize2, aes(x = "", y = pLength/1000)) +  # Use an empty x-axis
  geom_boxplot(fill = "lightblue", color = "black") +
  ylim(0,60000) +
  labs(title = "",
       x = "maize (WGS)",
       y = "Peak length (Kb)") +
  theme(axis.title.x = element_blank(),  # Remove x-axis title
        axis.text.x = element_blank(),  # Remove x-axis title
        axis.text.y = element_blank(),   # Remove x-axis text
        axis.ticks.y = element_blank()) + # Remove x-axis ticks
  coord_flip()  # Flip the plot to make it horizontal

# Create a horizontal boxplot without the x-axis
box3 <- ggplot(sorghum2.rna, aes(x = "", y = pLength/1000)) +  # Use an empty x-axis
  geom_boxplot(fill = "lightblue", color = "black") +
  ylim(0,60000) +
  labs(title = "",
       x = "sorghum (RNA-Seq)",
       y = "Peak length (Kb)") + # Remove x-axis ticks
  coord_flip()  # Flip the plot to make it horizontal


summary(maize2$pLength/1000)
summary(sorghum2.rna$pLength/1000)

box1 / box3

#### number of peaks 
#### leave only genes with 1 to 10 indivuals peaks

genesPeak <- data.frame(table(maize2$traits))
pm <- ggplot(genesPeak, aes(x = Freq)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  labs(x = "eQTLs per gene", y = "Count", title = "maize (WGS)")

genesPeakS <- data.frame(table(sorghum2$traits))
ps <- ggplot(genesPeakS, aes(x = Freq)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  labs(x = "eQTLs per gene", y = "Count", title = "sorghum (WGS)")


genesPeakS.rna <- data.frame(table(sorghum2.rna$traits))
ps.rna <- ggplot(genesPeakS.rna, aes(x = Freq)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  labs(x = "eQTLs per gene", y = "Count", title = "sorghum (RNA-Seq)")

pm + ps.rna


toremove <- genesPeak$Var1[which(genesPeak$Freq>20)] # we are not filtering
maize3 <- maize2[!maize2$traits %in% toremove,]
nrow(maize3)

23.9+9.8+44.4

######
eQTL.M <- round(length(unique(maize2$traits))/24585, 2)*100
eQTL.M
#78 % of genes have a eQTL in maize

# Count occurrences of each gene in cis and trans
cis_counts <- table(maize2$traits[maize2$cis_trans == "cis"])
trans_counts <- table(maize2$traits[maize2$cis_trans == "trans"])



# Count genes in each category
only_cis <- sum(cis_counts > 0 & !(names(cis_counts) %in% names(trans_counts)))
Trans_cis_trans <- sum(names(cis_counts) %in% names(trans_counts))
only_trans <- sum(trans_counts > 0 & !(names(trans_counts) %in% names(cis_counts)))

# Print results
cat("Cis:", only_cis, "\n")
cat("Cis + Trans:", Trans_cis_trans, "\n")
cat("Trans:", only_trans, "\n")


nrow(sorghum2)
eQTL.S <- round(length(unique(sorghum2$traits))/21430, 2)*100
eQTL.S
#26 % of genes have a eQTL in sorghum

# Count occurrences of each gene in cis and trans
cis_countsS <- table(sorghum2$traits[sorghum2$cis_trans == "cis"])
trans_countsS <- table(sorghum2$traits[sorghum2$cis_trans == "trans"])

# Count genes in each category
only_cisS <- sum(cis_countsS > 0 & !(names(cis_countsS) %in% names(trans_countsS)))
Trans_cis_transS <- sum(names(cis_countsS) %in% names(trans_countsS))
only_transS <- sum(trans_countsS > 0 & !(names(trans_countsS) %in% names(cis_countsS)))

# Print results
cat("Cis:", only_cisS, "\n")
cat("Cis + Trans:", Trans_cis_transS, "\n")
cat("Trans:", only_transS, "\n")

nrow(sorghum2)
eQTL.S.rna <- round(length(unique(sorghum2.rna$traits))/21430, 2)*100
eQTL.S.rna
#51 % of genes have a eQTL in sorghum RNA

#Count occurrences of each gene in cis and trans
cis_countsS.rna <- table(sorghum2.rna$traits[sorghum2.rna$cis_trans == "cis"])
trans_countsS.rna <- table(sorghum2.rna$traits[sorghum2.rna$cis_trans == "trans"])

# Count genes in each category
only_cisS.rna <- sum(cis_countsS.rna > 0 & !(names(cis_countsS.rna) %in% names(trans_countsS.rna)))
Trans_cis_transS.rna <- sum(names(cis_countsS.rna) %in% names(trans_countsS.rna))
only_transS.rna <- sum(trans_countsS.rna > 0 & !(names(trans_countsS.rna) %in% names(cis_countsS.rna)))


# Genes only in cis (not in trans)
only_cisS.rna2 <- names(cis_countsS.rna)[cis_countsS.rna > 0 & !(names(cis_countsS.rna) %in% names(trans_countsS.rna))]

# Genes in Cis + Trans cis and trans
Trans_cis_transS.rna2 <- intersect(names(cis_countsS.rna), names(trans_countsS.rna))

# Genes only in trans (not in cis)
only_transS.rna2 <- names(trans_countsS.rna)[trans_countsS.rna > 0 & !(names(trans_countsS.rna) %in% names(cis_countsS.rna))]


# Print results
cat("Cis:", only_cisS.rna, "\n")
cat("Cis + Trans:", Trans_cis_transS.rna, "\n")
cat("Trans:", only_transS.rna, "\n")


###########################
# Create a data frame for ggplot
df_plot <- data.frame(
  category = c("None","Cis", "Cis + Trans", "Trans"),
  specie = c("maize","maize","maize", "maize"),
  count = c((24585-sum(c(only_cis, Trans_cis_trans, only_trans))), only_cis, Trans_cis_trans, only_trans))

df_plot

# Calculate percentages
df_plot <- df_plot %>%
  group_by(category) %>%
  mutate(percentage = round(count / sum(df_plot$count) * 1, 2))


df_plotS.rna <- data.frame(
  category = c("None","Cis", "Cis + Trans", "Trans"),
  specie = c("sorghum","sorghum","sorghum","sorghum"),
  count = c((21430-sum(c(only_cisS.rna, Trans_cis_transS.rna, only_transS.rna))),only_cisS.rna, Trans_cis_transS.rna, only_transS.rna)
)


# Calculate percentages
df_plotS.rna <- df_plotS.rna %>%
  group_by(category) %>%
  mutate(percentage = round(count / sum(df_plotS.rna$count) * 1, 2))
df_plotS.rna


# Plot
# Define the order of categories
df_plot$category <- factor(df_plot$category, levels = c("None", "Cis + Trans", "Trans", "Cis"))

M <- ggplot(df_plot, aes(x = "Genes", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(count, " (", percentage, "%)")), 
            position = position_stack(vjust = 0.5), size = 5) + # Labels inside the bars
  scale_fill_manual(values = c("grey", "skyblue", "orange", "red")) + # Custom colors
  labs(title = "Gene Classification: Cis vs Trans", y = "Number of Genes", x = "") +
  theme(
    axis.text.x = element_blank(), # Hide x-axis labels
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), # Hide y-axis labels
    axis.ticks.y = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 14)
  )

# Define the order of categories
df_plotS.rna$category <- factor(df_plotS.rna$category, levels = c("None", "Cis + Trans", "Trans", "Cis"))


S.rna <- ggplot(df_plotS.rna, aes(x = "Genes", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(count, " (", percentage, "%)")), 
            position = position_stack(vjust = 0.5), size = 5) + # Labels inside the bars
  #annotate("text", x = 1, y = total_genesS.rna + 2, label = paste0("Total:", total_genesS.rna, " (", eQTL.S.rna, "%)"), size = 6, fontface = "bold") +
  scale_fill_manual(values = c("grey","skyblue", "orange", "red")) + # Custom colors
  labs(title = "Gene Classification: Cis vs Trans", y = "Number of Genes", x = "") +
  theme(
    axis.text.x = element_blank(), # Hide x-axis labels
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), # Hide y-axis labels
    axis.ticks.y = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 14)
  )


M + S.rna

# Combine the data frames into one
df_combined <- bind_rows(
  df_plot %>% mutate(Plot = "Maize"),
  df_plotS.rna %>% mutate(Plot = "Sorghum")
)


combined_plot <- ggplot(df_combined, aes(x = "Genes", y = percentage, fill = category)) +
  geom_bar(stat = "identity", width = 1.1) +
  geom_text(aes(label = count ), 
            position = position_stack(vjust = 0.5), size = 5) +
  facet_wrap(~ Plot, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = c("gray", "orange", "red", "skyblue")) +
  labs(y = "Proportion of genes", x = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1.0)) + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.85),
    strip.background = element_blank(),
    strip.placement = "outside"
  )

combined_plot


# heritability 

h2.maize <- fread("maize_693_genemodels_variance.csv", data.table = F)
h2.sorghum <- fread("sorghum_811_genemodels_variance.csv", data.table = F)

head(h2.maize$Heritability)
head(h2.sorghum)

mean(h2.maize$Heritability)
mean(h2.sorghum$Heritability)


h2.maize$eQTL <- ifelse(h2.maize$Trait %in% maize0$traits, "yes", "no")
h2.sorghum$eQTL <- ifelse(h2.sorghum$Trait %in% sorghum0$traits, "yes", "no")
table(h2.maize$eQTL)
table(h2.sorghum$eQTL)

# Add species label
h2.maize$Species <- "Maize"
h2.sorghum$Species <- "Sorghum"


# Combine datasets
h2.all <- bind_rows(h2.maize, h2.sorghum)
h2.all$linetype = ifelse(h2.all$eQTL == "yes", "solid", "dotted")

# Custom colors
custom_colors <- c("Maize" = "blue", 
                   "Sorghum" = "green4")

# Plot
h.2_eQTLs <- ggplot(h2.all, aes(x = Heritability, color = Species, linetype = linetype)) +
  geom_density(size = 1) +
  scale_color_manual(values = custom_colors, name = "eQTL Type") +
  scale_linetype_identity(name = "Species", labels = c("Maize" = "solid", "Sorghum" = "dotted")) +
  labs(x = "Heritability (h2)",
       y = "Density") +
  theme(legend.position = "inside",
    legend.position.inside = c(0.85, 0.85),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.box = "vertical"
  )

h.2_eQTLs

#relation of heritability with the number of trans-eQTLS


cis_counts.df <- data.frame(cis_counts)
cis_countsS.rna.df <- data.frame(cis_countsS.rna)

trans_counts.df <- data.frame(trans_counts)
trans_countsS.rna.df <- data.frame(trans_countsS.rna)


hist(trans_counts.df$Freq)
hist(trans_countsS.rna.df$Freq)

trans_countsS.rna.df2 <- merge(trans_countsS.rna.df, h2.sorghum, by=1)
trans_counts.df2 <- merge(trans_counts.df, h2.maize, by=1)


cor.test(trans_counts.df2$Freq, trans_counts.df2$Heritability)
cor.test(trans_countsS.rna.df2$Freq, trans_countsS.rna.df2$Heritability)

# Add species column
trans_counts.df2$Species <- "Maize"
trans_countsS.rna.df2$Species <- "Sorghum"

# Combine them into a single dataframe for plotting
trans_counts.all <- bind_rows(trans_counts.df2, trans_countsS.rna.df2)

head(trans_counts.all)

freq_jitter <- ggplot(trans_counts.all, aes(x = Species, y = Freq, color = Species)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_color_manual(values = custom_colors) +
  labs(x = "Species",
       y = "Number of eQTL per gene",
       color = "Species")

freq_jitter


h2vTeqtls <- ggplot(trans_counts.all, aes(x = Freq, y = Heritability, fill = Species)) +
  rasterise(geom_jitter(
    shape = 21,
    color = "black",
    alpha = 0.1,
    size = 2.5,
    width = 0.2,
    height = 0.1
  ), dpi = 600, dev = "ragg_png")  +
  geom_smooth(
    aes(color = Species),
    method = "lm",
    se = FALSE
  ) +
  scale_fill_manual(values = c("Maize" = "blue", "Sorghum" = "green4")) +
  scale_color_manual(values = c("Maize" = "blue", "Sorghum" = "green4")) +
  labs(
    x = "Number of trans-eQTLs",
    y = "Heritability"
  ) +
  theme(legend.position = "inside",
    legend.position.inside = c(.99, 0.01),       # bottom right inside
    legend.justification = c("right", "bottom") )
h2vTeqtls


fancyplots <- plot_grid(h.2_eQTLs, combined_plot, h2vTeqtls, rel_widths = c(1, 1, 1),
                 nrow = 1, labels = c('A', 'B', 'C'),hjust = 0, vjust = 1, label_size = 18)
fancyplots

ggsave("Fig1_v2.svg", width = 12, height = 4) # to save text as editable
