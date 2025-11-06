#Figure #2
setwd("/Users/vladimir/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/Previous_downloads")

library(GenomicRanges)
library(ggplot2)
library(ggnewscale)
library(data.table)
library(tidyverse)
library(ggrastr)
library(fuzzyjoin)
library(data.table)
library(patchwork)
library(ggpubr)
library(ggpmisc)
library(svglite)
library(cowplot) 
library(dplyr)
library(scales)
library(stringr)
library(forcats)  # Needed for fct_reorder


theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
             plot.title = element_text(hjust = 0.5), plot.subtitle=element_text(hjust=0.5))

#maize regions
maize <- fread("totalcombined_snps_summary_NS_maize_693_V2_2.csv", data.table = F)
maize0 <- maize[which(maize$num_SNPs >= 3),]

length(unique(maize0$traits))

maize.trans <- maize0[which(maize0$cis_trans == "trans"),]
head(maize.trans[,c(2,3,5)])

maize.trans.freqMarkers <- data.frame(table(maize.trans$top_SNP))


#-----------------------
# 1. Create GRanges for maize genome
#-----------------------
maize_genome <- data.frame(
  chr = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
  start = rep(1, 10),
  end = c(308452471, 243675191, 238017767, 250330460, 226353449, 
          181357234, 185808916, 182411202, 163004744, 152435371)
)

genome_gr <- GRanges(
  seqnames = maize_genome$chr,
  ranges = IRanges(start = maize_genome$start, end = maize_genome$end)
)
seqlengths(genome_gr) <- width(genome_gr)

#-----------------------
# 2. Create sliding windows: 1Mb size, 100kb step
#-----------------------
# Create a vector of all chr names and their lengths
all_chr <- maize_genome$chr
chr_lengths <- setNames(maize_genome$end, maize_genome$chr)

window_size <- 100000
step_size <- 100

make_windows <- function(chr, chr_len) {
  starts <- seq(1, chr_len - window_size, by = step_size)
  ends <- starts + window_size - 1
  GRanges(seqnames = chr, ranges = IRanges(start = starts, end = ends),
          seqlengths = chr_lengths)
}

# Apply to all chromosomes
window_list <- mapply(make_windows,
                      chr = all_chr,
                      chr_len = chr_lengths,
                      SIMPLIFY = FALSE)

# Flatten to one GRanges, safe merge
windows <- do.call(c, unname(window_list))

#-----------------------
# 3. Convert eQTL peaks to GRanges
#-----------------------
eqtl_gr <- GRanges(
  seqnames = maize.trans$CHROM,
  ranges = IRanges(start = maize.trans$POS, end = maize.trans$POS),
  seqlengths = chr_lengths
)

#-----------------------
# 4. Count overlaps between windows and eQTL peaks
#-----------------------
counts <- countOverlaps(windows, eqtl_gr)
mcols(windows)$eqtl_count <- counts

#-----------------------
# 5. Convert to data.frame and plot
#-----------------------
hotspot_df <- as.data.frame(windows)
hotspot_df$start_Mb <- hotspot_df$start / 1e6
hotspot_df$seqnames <- factor(hotspot_df$seqnames, levels = maize_genome$chr)

head(hotspot_df)

hotspot_plot <- hotspot_df

nCHR <- length(unique(hotspot_plot$seqnames))
hotspot_plot$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(hotspot_plot$seqnames))){
  nbp[i] <- max(hotspot_plot[hotspot_plot$seqnames == i,]$start)
  hotspot_plot[hotspot_plot$seqnames == i,"BPcum"] <- hotspot_plot[hotspot_plot$seqnames == i,"start"] + s
  s <- s + nbp[i]}

axis.set <- hotspot_plot %>% 
  group_by(seqnames) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

table(hotspot_plot$eqtl_count)

hotspot_plot.maize <- ggplot(hotspot_plot %>% filter(eqtl_count >= 3)) +
  ggrastr::rasterize(geom_segment(aes(x = BPcum, xend = BPcum,
                   y = 3,     yend = eqtl_count,
                   colour = factor(seqnames, levels = 1:10)),
               linewidth = 1.2), dpi = 900, dev = "ragg_png" ) +
  scale_color_manual(values = c("#000000","#E69F00","#56B4E9",
                                "#009E73","#F0E442","#0072B2",
                                "#D55E00","#CC79A7","#999999","#A6CEE3")) +
  scale_x_continuous(breaks = axis.set$center, labels = axis.set$seqnames) +
  ylim(3,100) +
  labs(y = "eQTL per 100 kb", x = "Chromosome") +
  theme(legend.position = "none", axis.title.x = element_blank())

#hotspot_plot.maize

#sorghum regions
sorghum <- fread("totalcombined_snps_summary_NS_sorghum_2_V2.csv", data.table = F)
sorghum0 <- sorghum[which(sorghum$num_SNPs >= 3),]

length(unique(sorghum0$traits))

sorghum.trans <- sorghum0[which(sorghum0$cis_trans == "trans"),]
head(sorghum.trans[,c(2,3,5)])

sorghum.trans.freqMarkers <- data.frame(table(sorghum.trans$top_SNP))

#-----------------------
# 1. Create GRanges for maize genome
#-----------------------
sorghum_genome <- data.frame(
  chr = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
  start = rep(1, 10),
  end = c(85112863, 79114963, 80873341, 71215609, 77058072, 
          62713908, 68911884, 65779274, 63277606, 62870657)
)


genome_gr <- GRanges(
  seqnames = sorghum_genome$chr,
  ranges = IRanges(start = sorghum_genome$start, end = sorghum_genome$end)
)
seqlengths(genome_gr) <- width(genome_gr)

#-----------------------
# 2. Create sliding windows: 1Mb size, 100kb step
#-----------------------
# Create a vector of all chr names and their lengths
all_chr <- sorghum_genome$chr
chr_lengths <- setNames(sorghum_genome$end, sorghum_genome$chr)

window_size <- 100000
step_size <- 100

make_windows <- function(chr, chr_len) {
  starts <- seq(1, chr_len - window_size, by = step_size)
  ends <- starts + window_size - 1
  GRanges(seqnames = chr, ranges = IRanges(start = starts, end = ends),
          seqlengths = chr_lengths)
}

# Apply to all chromosomes
window_list <- mapply(make_windows,
                      chr = all_chr,
                      chr_len = chr_lengths,
                      SIMPLIFY = FALSE)

# Flatten to one GRanges, safe merge
windows <- do.call(c, unname(window_list))

#-----------------------
# 3. Convert eQTL peaks to GRanges
#-----------------------
eqtl_gr <- GRanges(
  seqnames = sorghum.trans$CHROM,
  ranges = IRanges(start = sorghum.trans$POS, end = sorghum.trans$POS),
  seqlengths = chr_lengths
)

#-----------------------
# 4. Count overlaps between windows and eQTL peaks
#-----------------------
counts <- countOverlaps(windows, eqtl_gr)
mcols(windows)$eqtl_count <- counts

#-----------------------
# 5. Convert to data.frame and plot
#-----------------------
hotspot_dfS <- as.data.frame(windows)
hotspot_dfS$start_Mb <- hotspot_dfS$start / 1e6
hotspot_dfS$seqnames <- factor(hotspot_dfS$seqnames, levels = sorghum_genome$chr)


nCHR <- length(unique(hotspot_dfS$seqnames))
hotspot_dfS$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(hotspot_dfS$seqnames))){
  nbp[i] <- max(hotspot_dfS[hotspot_dfS$seqnames == i,]$start)
  hotspot_dfS[hotspot_dfS$seqnames == i,"BPcum"] <- hotspot_dfS[hotspot_dfS$seqnames == i,"start"] + s
  s <- s + nbp[i]}

axis.set <- hotspot_dfS %>% 
  group_by(seqnames) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

table(hotspot_dfS$eqtl_count)

hotspot_plot.sorghum <- ggplot(hotspot_dfS %>% filter(eqtl_count >= 3)) +
  ggrastr::rasterize(geom_segment(aes(x = BPcum, xend = BPcum,
                             y = 3,     yend = eqtl_count,
                             colour = factor(seqnames, levels = 1:10)),
                         linewidth = 1.2), dpi = 900, dev = "ragg_png" ) +
  scale_color_manual(values = c("#000000","#E69F00","#56B4E9",
                                "#009E73","#F0E442","#0072B2",
                                "#D55E00","#CC79A7","#999999","#A6CEE3")) +
  scale_x_continuous(breaks = axis.set$center, labels = axis.set$seqnames) +
  labs(y = "eQTL per 100 kb", x = "Chromosome") +
  #geom_hline(yintercept = 30, col="red", lwd=1, lty=2) +
  theme(legend.position = "none")

#hotspot_plot.maize / hotspot_plot.sorghum

#define interesting regions
hotspot_df_selected <- hotspot_dfS[which(hotspot_dfS$eqtl_count > 30),]

# Convert to GRanges object
gr.hot <- GRanges(seqnames = hotspot_df_selected$seqnames,
                  ranges = IRanges(start = hotspot_df_selected$start, end = hotspot_df_selected$end),
                  eqtl_count = hotspot_df_selected$eqtl_count)

# Reduce ranges allowing a 1kb gap (i.e., merge regions ≤1kb apart)
reduced <- GenomicRanges::reduce(gr.hot, min.gapwidth = 10000)
reduced
reduced_df <- as.data.frame(reduced)
reduced_df$region <- c(LETTERS, paste0("A", LETTERS))[1:nrow(reduced_df)]
reduced_df

#str(reduced_df)
#str(hotspot_dfS)

reduced_df$eqtl_count <- NA
reduced_df$BPcum <- NA

for (i in 1:nrow(reduced_df)) {
  chr=reduced_df[i,1]
  start=reduced_df[i,2]
  BPcum <- hotspot_dfS$BPcum[which(hotspot_dfS$seqnames == chr & hotspot_dfS$start == start)]
  eqtl_count <- hotspot_dfS$eqtl_count[which(hotspot_dfS$seqnames == chr & hotspot_dfS$start == start)]
  reduced_df$BPcum[which(reduced_df$seqnames == chr & reduced_df$start == start)] <- BPcum 
  reduced_df$eqtl_count[which(reduced_df$seqnames == chr & reduced_df$start == start)] <- eqtl_count 
}

#zoom-in three major regions
# Create a vector of all chr names and their lengths
reduced_chr <- reduced_df$seqnames
reduced_chr_lengths <- setNames(reduced_df$end, reduced_df$seqnames)

# Define the window and step size
window_size2 <- 100
step_size2 <- 10

# Define a function that creates windows *within a region*
make_windows_region <- function(chr, start_pos, end_pos) {
  region_len <- end_pos - start_pos + 1
  if (region_len < window_size2) return(NULL)  # skip if too small
  starts <- seq(start_pos, end_pos - window_size2 + 1, by = step_size2)
  ends <- starts + window_size2 - 1
  GRanges(seqnames = chr, ranges = IRanges(start = starts, end = ends))
}

# Apply to each row of reduced_df
window_list2 <- mapply(make_windows_region,
                       chr = reduced_df$seqnames,
                       start_pos = reduced_df$start,
                       end_pos = reduced_df$end,
                       SIMPLIFY = FALSE)

# Combine into a single GRanges object
all_windows <- do.call(c, window_list2)

#-----------------------
# 4. Count overlaps between windows and eQTL peaks
#-----------------------
counts_regions <- countOverlaps(all_windows, eqtl_gr)
mcols(all_windows)$eqtl_count <- counts_regions

hotspot_df_regions <- as.data.frame(all_windows)
hotspot_df_regions$start_Mb <- hotspot_df_regions$start / 1e6
hotspot_df_regions$seqnames <- factor(hotspot_df_regions$seqnames, levels = sorghum_genome$chr)
hotspot_df_regions

hotspot_df_regions <- hotspot_df_regions[order(-hotspot_df_regions$eqtl_count),]

head(hotspot_df_regions)

#############
hotspot_df_regions2 <- hotspot_df_regions[which(hotspot_df_regions$eqtl_count > 10),]

# Convert to GRanges object
gr.hot_regions <- GRanges(seqnames = hotspot_df_regions2$seqnames,
                  ranges = IRanges(start = hotspot_df_regions2$start, end = hotspot_df_regions2$end),
                  eqtl_count = hotspot_df_regions2$eqtl_count)

# Reduce ranges allowing a 10000 gap (i.e., merge regions ≤1000bp apart)
reduced_regions <- GenomicRanges::reduce(gr.hot_regions, min.gapwidth = 10000)
reduced_regions


# Find overlaps between reduced merged regions and original ones
hits_regions <- findOverlaps(reduced_regions, gr.hot_regions)
hits_regions

# Summarize: for each merged region, calculate start, end, and mean eqtl_count
merged_df_regions <- as.data.frame(reduced_regions)[, c("seqnames", "start", "end")]

colnames(merged_df_regions)[1] <- "chr"
merged_df_regions$width <- round((merged_df_regions$end-merged_df_regions$start),2)
merged_df_regions <- merged_df_regions[order(merged_df_regions$chr),]
colnames(merged_df_regions)[1] <- "chr"


# View result
head(merged_df_regions)

# this is to add group back (labels/letters)
merged_df_regions_letters <- merge(merged_df_regions, reduced_df, by = 1,suffixes= c("_small", "_big") )
merged_df_regions_letters2 <- merged_df_regions_letters[with(merged_df_regions_letters, start_small >= start_big & end_small <= end_big), ]
merged_df_regions_letters2

merged_df_regions_letters2 <- merged_df_regions_letters2 %>%
  select(-start_big, -end_big, -width_big, -strand, -eqtl_count)


# add how many eQTLs are controlled by such regions, max is the number of eqtl is controlled by a single marker
#but that does not include all eQTL per region
regions <- merged_df_regions_letters2 %>%
  dplyr::rename(
    CHROM = chr,
    start = start_small,
    end = end_small
  ) %>%
  mutate(CHROM = as.character(CHROM))

# Prepare gene data
genes <- sorghum.trans %>%
  mutate(CHROM = as.character(CHROM))

# Add interval columns for POS (treating it as a 1-bp "interval")
genes_for_join <- genes %>%
  mutate(start = POS, end = POS)

head(genes_for_join)

# Perform interval join (POS within region_start–region_end on same chr)
joined_data <- genome_inner_join(
  genes_for_join, regions,
  by = c("CHROM", "start", "end")
)

controlledGenes <- data.frame(table(joined_data$region))
merged_df_regions_letters2 <- merge(merged_df_regions_letters2, controlledGenes, by.x= "region", by.y = "Var1")
colnames(merged_df_regions_letters2)[which(colnames(merged_df_regions_letters2) == "Freq")] <- "controlledGenes"
merged_df_regions_letters2

merged_df_regions_letters3 <- merged_df_regions_letters2[merged_df_regions_letters2$controlledGenes >= 30,]
merged_df_regions_letters3$regionRAW <- merged_df_regions_letters3$region
nrow(merged_df_regions_letters3)

merged_df_regions_letters3 <- merged_df_regions_letters3[order(merged_df_regions_letters3$controlledGenes, decreasing = T),]

merged_df_regions_letters3$region <- c(LETTERS, paste0("A", LETTERS))[1:nrow(merged_df_regions_letters3)]

fwrite(merged_df_regions_letters3[,-8], "hotspots_sorghumRNA.csv")

hotspot_plot.sorghum.regions <- hotspot_plot.sorghum +
  # Add region letters
  geom_text(
    data = merged_df_regions_letters3,
    aes(x = BPcum, y = controlledGenes, label = region),
    color = "black",
    size = 6,
    fontface = "bold"
  )


#hotspot_plot.sorghum.regions


joined_data2 <- joined_data %>%
  dplyr::select(region, CHROM.x, POS, top_SNP, top_Pvalue, pStart, pStop, pLength, num_SNPs, gene, chr,  start.x, end.x, cis_trans) %>%
  dplyr::rename(region.old = region, chr_eQTL = CHROM.x, pos = POS, gene_start = start.x, gene_end = end.x, chr_gene = chr)
joined_data2

table(joined_data2$region.old)
# remove region menor of 30 = K, M and O

joined_data3 <- joined_data2[!joined_data2$region %in% c("K","M","O"),]
table(joined_data3$region.old)

#change names of regions
joined_data3$region <- NA
#i="F"
for (i in merged_df_regions_letters3$regionRAW) {
  letter <- merged_df_regions_letters3$region[which(merged_df_regions_letters3$regionRAW == i)]
  joined_data3$region[which(joined_data3$region.old == i)] <- letter
}
table(joined_data3$region)

fwrite(joined_data3, "hotspots_sorghumRNA_geneInregions.csv")

#############################
################
#plot 

hotspot_plot.maize / hotspot_plot.sorghum.regions


### zoom in genes

merged_df_regions_letters2

group="A"

region.old <- merged_df_regions_letters3[which(merged_df_regions_letters3$region == group),"regionRAW"]

reduced_df4plot <- merged_df_regions_letters[merged_df_regions_letters$region == region.old,]
reduced_df4plot <- reduced_df4plot[1,]
merged_df_regions_letters2_4plot <- merged_df_regions_letters3[which(merged_df_regions_letters3$region == group),]
merged_df_regions_letters2_4plot


reduced_df4plot2 <- sorghum.trans %>%
  filter(CHROM == reduced_df4plot$chr,
         POS>=reduced_df4plot$start_big, POS<= reduced_df4plot$end_big)

merged_df_regions_letters2_4plot
#this one is to plot the zoom in window
merged_df_regions_letters2_4plot2 <- sorghum.trans %>%
  filter(CHROM == merged_df_regions_letters2_4plot$chr[1], 
         POS>=merged_df_regions_letters2_4plot$start_small[1]-1000, POS<= merged_df_regions_letters2_4plot$end_small[1]+1000)

#to to barplots later on
merged_df_regions_letters2_4plot2.1 <- sorghum.trans %>%
  filter(CHROM == merged_df_regions_letters2_4plot$chr[1], 
         POS>=merged_df_regions_letters2_4plot$start_small[1], POS<= merged_df_regions_letters2_4plot$end_small[1])

#to to barplots later on
merged_df_regions_letters2_4plot2.2 <- sorghum.trans %>%
  filter(CHROM == merged_df_regions_letters2_4plot$chr[2], 
         POS>=merged_df_regions_letters2_4plot$start_small[2], POS<= merged_df_regions_letters2_4plot$end_small[2])

#to to barplots later on
merged_df_regions_letters2_4plot2.3 <- sorghum.trans %>%
  filter(CHROM == merged_df_regions_letters2_4plot$chr[3], 
         POS>=merged_df_regions_letters2_4plot$start_small[3], POS<= merged_df_regions_letters2_4plot$end_small[3])


merged_df_regions_letters2_4plot3 <- merged_df_regions_letters2_4plot2 %>%
  group_by(top_SNP)  %>% 
  summarise(genes_afected = sum(!is.na(gene)), .groups = 'drop')

merged_df_regions_letters2_4plot3.1 <- merged_df_regions_letters2_4plot2.1 %>%
  group_by(top_SNP)  %>% 
  summarise(genes_afected = sum(!is.na(gene)), .groups = 'drop')

merged_df_regions_letters2_4plot3.2 <- merged_df_regions_letters2_4plot2.2 %>%
  group_by(top_SNP)  %>% 
  summarise(genes_afected = sum(!is.na(gene)), .groups = 'drop')

merged_df_regions_letters2_4plot3.3 <- merged_df_regions_letters2_4plot2.3 %>%
  group_by(top_SNP)  %>% 
  summarise(genes_afected = sum(!is.na(gene)), .groups = 'drop')

merged_df_regions_letters2_4plot3.1
merged_df_regions_letters2_4plot3.2
merged_df_regions_letters2_4plot3.3

# get windows in the region to plot them
reduced_df4plot2.2 <- hotspot_dfS %>%
  filter(seqnames == reduced_df4plot$chr,
         start>=reduced_df4plot$start_big, end<= reduced_df4plot$end_big)

#order them
reduced_df4plot2.2 <- reduced_df4plot2.2[order(reduced_df4plot2.2$start),]
reduced_df4plot2.2$level <- 1:nrow(reduced_df4plot2.2)


#### get exon information
## coding versus no coding in maize?
# regions <- fread("Sbicolor_730_v5.1.gene.gff3.gz")
# head(regions)
# regions$ID <- sub("ID=([^.]+\\.[^.]+\\.[^\\.]+)\\..*", "\\1", regions$V9)
# 
# regions2 <- regions[,c(10,3,1,4,5,7)]
# head(regions2)
# 
# #filter transcripts used in the quantification
# transcripts <- fread("Sbicolor_730_v5.1.transcript_primaryTranscriptOnly.fa", header = F)
# head(transcripts)
# 
# # Filter rows that start with ">"
# transcript_headers <- transcripts[grep("^>", transcripts$V1), ]
# 
# # Extract only the transcript name (between '>' and first space)
# transcript_names <- sub("^>(\\S+).*", "\\1", transcript_headers$V1)
# 
# # Convert to a df
# filtered_transcripts <- data.frame(Transcript_ID = transcript_names)
# 
# regions3 <- regions2[regions2$ID %in% filtered_transcripts$Transcript_ID,]
# colnames(regions3) <- c("gene","type","chr","start","end","string")
# 
# regions3$gene <- gsub("\\.[^.]*$","", regions3$gene)
# regions3$chr <- gsub("Chr_0","", regions3$chr)
# regions3$chr <- as.numeric(as.character(regions3$chr))
# regions3 <- na.omit(regions3)
# 
# 
# fwrite(regions3, "Sbicolor_730_v5.1_primary_exons.csv")

# get the information from genes
Sbicolor_730_v5.1_primary_exons <- fread("Sbicolor_730_v5.1_primary_exons.csv", data.table = F)
#change names
Sbicolor_730_v5.1_primary_exons$type[which(Sbicolor_730_v5.1_primary_exons$type == "five_prime_UTR")] <- "5'UTR"
Sbicolor_730_v5.1_primary_exons$type[which(Sbicolor_730_v5.1_primary_exons$type == "three_prime_UTR")] <- "3'UTR"

head(Sbicolor_730_v5.1_primary_exons)
Sbicolor_730_v5.1_primary_exons_4plot2 <- Sbicolor_730_v5.1_primary_exons %>%
  filter(chr == reduced_df4plot$chr,
         start>=reduced_df4plot$start_big, end<= reduced_df4plot$end_big, type %in% c("5'UTR", "CDS", "3'UTR"))



####
#Calculate LD
#singleMarker <- merged_df_regions_letters2_4plot2$top_SNP[which(merged_df_regions_letters2_4plot2$num_SNPs ==
#                                                  max(merged_df_regions_letters2_4plot2$num_SNPs))]

merged_df_regions_letters2_4plot3.1
merged_df_regions_letters2_4plot3.2
merged_df_regions_letters2_4plot3.3 #select the one with more snps
singleMarker <- "Chr03_63696436" #Chr03_63696487 Chr03_63696436

#SNP.all <- unique(reduced_df4plot2$top_SNP)
# ### calculate LD
#SNP.all <- data.frame(c(maize0.e340$V1, e340_eSummary_filt.csv$snps))
#write.table(SNP.all,paste0("SNP_maize",group,".txt"),row.names=FALSE,sep="\t", quote = FALSE, col.names = F)
#system(paste0("wc -l ",paste0("SNP_maize",group,".txt")))


plinkPath <- "plink"
#name <- paste0("SNP_maize",group,".txt")
#name
#plink_command <- paste0(plinkPath," --bfile vla_karla_811geno_filter2_maf_het --r2 --ld-window-r2 0 --ld-window-kb 1000 --ld-window 9999999  --extract ", name," --out ", name,"LD")
plink_command <- paste0(plinkPath," --bfile vla_karla_811geno_filter2_maf_het --r2 --ld-window-r2 0 --ld-window-kb 10000000 --ld-window 99999999  --ld-snp ", singleMarker , " --out ", singleMarker)

system(plink_command)


#load plink results
#LD.true <- fread(paste0(name,"LD.ld"),data.table = F)
LD.true <- fread(paste0(singleMarker,".ld"),data.table = F)

### this will tell us the LD between the reference marker (the one with more hits in group one) against the other markers in the other region
LD.true[LD.true$SNP_B %in% merged_df_regions_letters2_4plot3.1$top_SNP,]
merged_df_regions_letters2_4plot

LD.true_for_plot1 <- LD.true[which(LD.true$BP_B >= merged_df_regions_letters2_4plot$start_small[1]-5000 &
                                     LD.true$BP_B <= merged_df_regions_letters2_4plot$end_small[1]+5000),]

LD.true_for_plot2 <- LD.true[which(LD.true$BP_B >= merged_df_regions_letters2_4plot$start_small[3]-5000 &
                                     LD.true$BP_B <= merged_df_regions_letters2_4plot$end_small[3]+5000),]

# file: plots/ggplot_zoom_from_POS.R
# Minimal, safe zoom centered on a focal marker POS (in base pairs)
# why: coord_cartesian zooms the view without dropping data used by geoms/stats

# ---- adjustable parameters ----

focus_pos_bp   <- 63696436        # <- set your focal marker POS (in bp)
half_window_kb <- 25             # <- set half window size (in kb) around POS
half_window_bp <- half_window_kb * 1000

pw1 <- ggplot(reduced_df4plot2, aes(x = POS/1000000, y = nrow(reduced_df4plot2.2)/2)) +
  # Background rectangles
  # geom_rect(data = reduced_df4plot2.2,
  #           inherit.aes = FALSE,
  #           aes(xmin = start/1000000, xmax = end/1000000,
  #               ymin = level - 0.4, ymax = level + 0.4),
  #           fill = "lightblue", alpha = 0.5) +
  # Vertical lines
  geom_vline(xintercept = merged_df_regions_letters2_4plot$start_small/1000000, 
             color = "red", linewidth = 0.5) +
  geom_vline(xintercept = merged_df_regions_letters2_4plot$end_small/1000000, 
             color = "red", linewidth = 0.5) +
  # Gene structure: Exons (UTRs and CDS as boxes)
  geom_rect(data = Sbicolor_730_v5.1_primary_exons_4plot2,
            inherit.aes = FALSE,
            aes(xmin = start/1000000, xmax = end/1000000,
                ymin = 0, ymax = nrow(reduced_df4plot2.2)*.05,
                fill = type),
            color = "black") +  # black border around exons
  # Gene structure: Introns (lines between exons)
  geom_segment(data = Sbicolor_730_v5.1_primary_exons_4plot2 %>%
                 dplyr::group_by(gene) %>%
                 dplyr::arrange(start) %>%
                 dplyr::mutate(next_start = dplyr::lead(start),
                               next_end = dplyr::lead(end)) %>%
                 dplyr::filter(!is.na(next_start)),
               inherit.aes = FALSE,
               aes(x = end/1000000,
                   xend = next_start/1000000,
                   y = (nrow(reduced_df4plot2.2)*.05)/2,
                   yend = (nrow(reduced_df4plot2.2)*.05)/2),
               color = "black",
               linewidth = 1) +
  # Jitter points ON TOP
  geom_jitter(shape = 6, height = nrow(reduced_df4plot2.2)/4, width = 0, size = 5) +
  # Axis and labels (kept; view will be zoomed by coord_cartesian below)
  scale_x_continuous(
    limits = c(min(reduced_df4plot$start_big/1000000), max(reduced_df4plot$end_big/1000000)),
    expand = c(0, 0),
    breaks = round(seq(from = min(reduced_df4plot$start_big/1000000), 
                       to = max(reduced_df4plot$end_big/1000000), 
                       by = 0.01), 2)
  ) +
  scale_fill_manual(values = c("5'UTR" = "skyblue", 
                               "3'UTR" = "grey", 
                               "CDS" = "blue")) +
  # ---- the zoom line (adjust focus_pos_bp / half_window_kb above) ----
coord_cartesian(xlim = c((focus_pos_bp - half_window_bp)/1e6,
                         (focus_pos_bp + half_window_bp)/1e6)) +
  labs(x = "Genomic Position (Mb)", y = "", fill = "Feature") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

pw1

# Optional: if you prefer auto-pretty ticks at any zoom level, replace the scale_x_continuous() above with:
# scale_x_continuous(expand = c(0, 0), breaks = scales::pretty_breaks(n = 8))


# Calculate shared y-axis limits
ymax <- max(
  max(merged_df_regions_letters2_4plot3.1$genes_afected, na.rm = TRUE),
  max(merged_df_regions_letters2_4plot3.2$genes_afected, na.rm = TRUE),
  max(merged_df_regions_letters2_4plot3.3$genes_afected, na.rm = TRUE)
)

# Modify the data to create new labels
# Prepare the data: format labels and reorder top_SNP
plot_df <- merged_df_regions_letters2_4plot3.1 %>%
  mutate(
    chr = str_extract(top_SNP, "\\d+"),
    pos = str_extract(top_SNP, "\\d+$"),
    formatted_label = paste0(as.integer(chr), ":", comma(as.numeric(pos))),
    top_SNP = fct_reorder(top_SNP, genes_afected, .desc = TRUE)  # reorder factor
  )

# Create the plot
# Filter top 5 based on genes_afected
top5_df <- plot_df %>%
  arrange(desc(genes_afected)) %>%
  dplyr::slice(1:5)

# Plot only the top 5
genes_afected.1 <- ggplot(data = top5_df, 
                          aes(x = fct_reorder(top_SNP, genes_afected), 
                              y = genes_afected)) +
  geom_bar(stat = "identity", color = "black", fill = "red") +
  geom_text(aes(y = 5, label = formatted_label), 
            angle = 0, hjust = -0.1, size = 7) +
  ylab("Genes with the strongest association") +
  xlab("Marker") +
  ylim(0, max(top5_df$genes_afected) * 1.1) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip()

#pw1 / (pw2 + genes_afected.1)
# A.2 <- plot_grid(pw1,middle,bottom, ncol = 1, hjust = 0, vjust = 1, label_size = 18,  labels = c('A','',''))
# A.2 <- plot_grid(pw2, genes_afected.1, rel_widths = c(1,.6),
#                     nrow = 1, labels = c('D','E'),hjust = 0, vjust = 1, label_size = 18)
# A <- plot_grid(pw1, A.2, nrow = 2, labels = c('C',''),hjust = 0, vjust = 1, label_size = 18)
# A
A <- plot_grid(pw1, genes_afected.1, nrow = 2, rel_heights = c(1,1))
A



######################

group="B"

region.old <- merged_df_regions_letters3[which(merged_df_regions_letters3$region == group),"regionRAW"]

reduced_df4plot <- merged_df_regions_letters[merged_df_regions_letters$region == region.old,]
reduced_df4plot <- reduced_df4plot[1,]
merged_df_regions_letters2_4plot <- merged_df_regions_letters3[which(merged_df_regions_letters3$region == group),]
merged_df_regions_letters2_4plot


reduced_df4plot2 <- sorghum.trans %>%
  filter(CHROM == reduced_df4plot$chr,
         POS>=reduced_df4plot$start_big, POS<= reduced_df4plot$end_big)

merged_df_regions_letters2_4plot
#this one is to plot the zoom in window
merged_df_regions_letters2_4plot2 <- sorghum.trans %>%
  filter(CHROM == merged_df_regions_letters2_4plot$chr[1], 
         POS>=merged_df_regions_letters2_4plot$start_small[1]-1000, POS<= merged_df_regions_letters2_4plot$end_small[1]+1000)

#to to barplots later on
merged_df_regions_letters2_4plot2.1 <- sorghum.trans %>%
  filter(CHROM == merged_df_regions_letters2_4plot$chr[1], 
         POS>=merged_df_regions_letters2_4plot$start_small[1], POS<= merged_df_regions_letters2_4plot$end_small[1])

#to to barplots later on
merged_df_regions_letters2_4plot2.2 <- sorghum.trans %>%
  filter(CHROM == merged_df_regions_letters2_4plot$chr[2], 
         POS>=merged_df_regions_letters2_4plot$start_small[2], POS<= merged_df_regions_letters2_4plot$end_small[2])

#to to barplots later on
merged_df_regions_letters2_4plot2.3 <- sorghum.trans %>%
  filter(CHROM == merged_df_regions_letters2_4plot$chr[3], 
         POS>=merged_df_regions_letters2_4plot$start_small[3], POS<= merged_df_regions_letters2_4plot$end_small[3])


merged_df_regions_letters2_4plot3 <- merged_df_regions_letters2_4plot2 %>%
  group_by(top_SNP)  %>% 
  summarise(genes_afected = sum(!is.na(gene)), .groups = 'drop')

merged_df_regions_letters2_4plot3.1 <- merged_df_regions_letters2_4plot2.1 %>%
  group_by(top_SNP)  %>% 
  summarise(genes_afected = sum(!is.na(gene)), .groups = 'drop')

merged_df_regions_letters2_4plot3.2 <- merged_df_regions_letters2_4plot2.2 %>%
  group_by(top_SNP)  %>% 
  summarise(genes_afected = sum(!is.na(gene)), .groups = 'drop')

merged_df_regions_letters2_4plot3.3 <- merged_df_regions_letters2_4plot2.3 %>%
  group_by(top_SNP)  %>% 
  summarise(genes_afected = sum(!is.na(gene)), .groups = 'drop')

merged_df_regions_letters2_4plot3.1
merged_df_regions_letters2_4plot3.2
merged_df_regions_letters2_4plot3.3

# get windows in the region to plot them
reduced_df4plot2.2 <- hotspot_dfS %>%
  filter(seqnames == reduced_df4plot$chr,
         start>=reduced_df4plot$start_big, end<= reduced_df4plot$end_big)

#order them
reduced_df4plot2.2 <- reduced_df4plot2.2[order(reduced_df4plot2.2$start),]
reduced_df4plot2.2$level <- 1:nrow(reduced_df4plot2.2)


# get the information from genes
Sbicolor_730_v5.1_primary_exons <- fread("Sbicolor_730_v5.1_primary_exons.csv", data.table = F)
#change names
Sbicolor_730_v5.1_primary_exons$type[which(Sbicolor_730_v5.1_primary_exons$type == "five_prime_UTR")] <- "5'UTR"
Sbicolor_730_v5.1_primary_exons$type[which(Sbicolor_730_v5.1_primary_exons$type == "three_prime_UTR")] <- "3'UTR"

head(Sbicolor_730_v5.1_primary_exons)
Sbicolor_730_v5.1_primary_exons_4plot2 <- Sbicolor_730_v5.1_primary_exons %>%
  filter(chr == reduced_df4plot$chr,
         start>=reduced_df4plot$start_big, end<= reduced_df4plot$end_big, type %in% c("5'UTR", "CDS", "3'UTR"))



####
#Calculate LD
#singleMarker <- merged_df_regions_letters2_4plot2$top_SNP[which(merged_df_regions_letters2_4plot2$num_SNPs ==
#                                                  max(merged_df_regions_letters2_4plot2$num_SNPs))]

merged_df_regions_letters2_4plot3.1
merged_df_regions_letters2_4plot3.2
merged_df_regions_letters2_4plot3.3 #select the one with more snps
singleMarker <- "Chr10_53899820" #Chr03_63696487 Chr03_63696436

plinkPath <- "plink"
#name <- paste0("SNP_maize",group,".txt")
#name
#plink_command <- paste0(plinkPath," --bfile vla_karla_811geno_filter2_maf_het --r2 --ld-window-r2 0 --ld-window-kb 1000 --ld-window 9999999  --extract ", name," --out ", name,"LD")
plink_command <- paste0(plinkPath," --bfile vla_karla_811geno_filter2_maf_het --r2 --ld-window-r2 0 --ld-window-kb 100000 --ld-window 999999  --ld-snp ", singleMarker , " --out ", singleMarker)

system(plink_command)


#load plink results
#LD.true <- fread(paste0(name,"LD.ld"),data.table = F)
LD.true <- fread(paste0(singleMarker,".ld"),data.table = F)

### this will tell us the LD between the reference marker (the one with more hits in group one) against the other markers in the other region
LD.true[LD.true$SNP_B %in% merged_df_regions_letters2_4plot3.1$top_SNP,]
merged_df_regions_letters2_4plot

LD.true_for_plot1 <- LD.true[which(LD.true$BP_B >= merged_df_regions_letters2_4plot$start_small[1]-5000 &
                                     LD.true$BP_B <= merged_df_regions_letters2_4plot$end_small[1]+5000),]

LD.true_for_plot2 <- LD.true[which(LD.true$BP_B >= merged_df_regions_letters2_4plot$start_small[3]-5000 &
                                     LD.true$BP_B <= merged_df_regions_letters2_4plot$end_small[3]+5000),]

max(reduced_df4plot2$POS)-min(reduced_df4plot2$POS)



focus_pos_bp   <- 53899820        # <- set your focal marker POS (in bp)
half_window_kb <- 25             # <- set half window size (in kb) around POS
half_window_bp <- half_window_kb * 1000

pw1 <- ggplot(reduced_df4plot2, aes(x = POS/1000000, y = nrow(reduced_df4plot2.2)/2)) +
  # Background rectangles
  # geom_rect(data = reduced_df4plot2.2,
  #           inherit.aes = FALSE,
  #           aes(xmin = start/1000000, xmax = end/1000000,
  #               ymin = level - 0.4, ymax = level + 0.4),
  #           fill = "lightblue", alpha = 0.5) +
  # Vertical lines
  geom_vline(xintercept = merged_df_regions_letters2_4plot$start_small/1000000, 
             color = "red", linewidth = 0.5) +
  geom_vline(xintercept = merged_df_regions_letters2_4plot$end_small/1000000, 
             color = "red", linewidth = 0.5) +
  # Gene structure: Exons (UTRs and CDS as boxes)
  geom_rect(data = Sbicolor_730_v5.1_primary_exons_4plot2,
            inherit.aes = FALSE,
            aes(xmin = start/1000000, xmax = end/1000000,
                ymin = 0, ymax = nrow(reduced_df4plot2.2)*.05,
                fill = type),
            color = "black") +  # black border around exons
  
  # Gene structure: Introns (lines between exons)
  geom_segment(data = Sbicolor_730_v5.1_primary_exons_4plot2 %>%
                 group_by(gene) %>%
                 arrange(start) %>%
                 mutate(next_start = lead(start),
                        next_end = lead(end)) %>%
                 filter(!is.na(next_start)),
               inherit.aes = FALSE,
               aes(x = end/1000000,
                   xend = next_start/1000000,
                   y = (nrow(reduced_df4plot2.2)*.05)/2,
                   yend = (nrow(reduced_df4plot2.2)*.05)/2),
               color = "black",
               linewidth = 1) +
  # Jitter points ON TOP
  geom_jitter(shape = 6, height = nrow(reduced_df4plot2.2)/4, width = 0, size = 5) +
  # Axis and labels
  scale_x_continuous(
    limits = c(min(reduced_df4plot$start_big/1000000), max(reduced_df4plot$end_big/1000000)),
    expand = c(0, 0),
    breaks = round(seq(from = min(reduced_df4plot$start_big/1000000), 
                       to = max(reduced_df4plot$end_big/1000000), 
                       by = 0.01), 2)
  ) +
  scale_fill_manual(values = c("5'UTR" = "skyblue", 
                               "3'UTR" = "grey", 
                               "CDS" = "blue")) +
  # ---- the zoom line (adjust focus_pos_bp / half_window_kb above) ----
coord_cartesian(xlim = c((focus_pos_bp - half_window_bp)/1e6,
                         (focus_pos_bp + half_window_bp)/1e6)) +
  labs(x = "Genomic Position (Mb)", y = "", fill = "Feature") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
pw1

### Zoom-in plot with gene models
# pw2 <- ggplot(LD.true_for_plot1[which(!LD.true_for_plot1$R2 == 1),], aes(x = BP_B/1000000, y = R2)) +
#   rasterize(geom_point(shape = 19, size = 5, color="grey"), dpi = 900, dev = "ragg_png") +
#   scale_x_continuous(
#     expand = c(0, 0),
#     breaks = seq(from = (merged_df_regions_letters2_4plot$start_small[1]-6000)/1000000, 
#                  to = (merged_df_regions_letters2_4plot$end_small[1]+6000)/1000000, 
#                  by = 0.001)
#   ) +
#   ylim(0,1) +
#   # Vertical lines
#   geom_vline(xintercept = merged_df_regions_letters2_4plot$start_small/1000000, 
#              color = "red", linewidth = 0.5) +
#   geom_vline(xintercept = merged_df_regions_letters2_4plot$end_small/1000000, 
#              color = "red", linewidth = 0.5) +
#   
#   # Gene structure: Exons (UTRs and CDS as boxes)
#   geom_rect(data = Sbicolor_730_v5.1_primary_exons_4plot2,
#             inherit.aes = FALSE,
#             aes(xmin = start/1000000, xmax = end/1000000,
#                 ymin = 0, ymax = .1,
#                 fill = type),
#             color = "black") +  # black border around exons
#   
#   # Gene structure: Introns (lines between exons)
#   geom_segment(data = Sbicolor_730_v5.1_primary_exons_4plot2 %>%
#                  group_by(gene) %>%
#                  arrange(start) %>%
#                  mutate(next_start = lead(start),
#                         next_end = lead(end)) %>%
#                  filter(!is.na(next_start)),
#                inherit.aes = FALSE,
#                aes(x = end/1000000,
#                    xend = next_start/1000000,
#                    y = .05,
#                    yend = 0.05),
#                color = "black",
#                linewidth = 4) +
#   # Jitter points ON TOP
#   geom_jitter(data = merged_df_regions_letters2_4plot2, aes(x= POS/1e6, y=.5) , shape = 6, height = .1, width = 0, size = 5) +
#   scale_fill_manual(values = c("5'UTR" = "skyblue", 
#                                "3'UTR" = "grey", 
#                                "CDS" = "blue")) +
#   
#   labs(x = "Genomic Position (Mb)", y = "LD", fill = "Feature") +
#   # ggnewscale::new_scale_fill() +  # This allows adding a *new* fill scale
#   # geom_rect(data = LD.true_for_plot,
#   #           aes(xmin = BP_B - 0.0001, xmax = BP_B + 0.0001,
#   #               ymin = nrow(reduced_df4plot2.2)*.90,
#   #               ymax = nrow(reduced_df4plot2.2)*.95,
#   #               fill = R2),
#   #           inherit.aes = FALSE) +
#   # scale_fill_gradient(low = "blue", high = "red", name = expression(R^2)) +
#   theme(legend.position = "inside", axis.text.x = element_text(angle = 15, hjust = 1), 
#         legend.position.inside = c(.3, .99),       # bottom right inside
#         legend.justification = c("right", "top")) +
#   coord_cartesian(xlim = c((merged_df_regions_letters2_4plot$start_small[1]-6000)/1e6,
#                            (merged_df_regions_letters2_4plot$end_small[1]+6000)/1e6), clip="on")
# 
# pw2


#pw1/pw2

# Calculate shared y-axis limits
ymax <- max(
  max(merged_df_regions_letters2_4plot3.1$genes_afected, na.rm = TRUE),
  max(merged_df_regions_letters2_4plot3.2$genes_afected, na.rm = TRUE),
  max(merged_df_regions_letters2_4plot3.3$genes_afected, na.rm = TRUE)
)

# Modify the data to create new labels
# Prepare the data: format labels and reorder top_SNP
plot_df <- merged_df_regions_letters2_4plot3.1 %>%
  mutate(
    chr = str_extract(top_SNP, "\\d+"),
    pos = str_extract(top_SNP, "\\d+$"),
    formatted_label = paste0(as.integer(chr), ":", comma(as.numeric(pos))),
    top_SNP = fct_reorder(top_SNP, genes_afected, .desc = TRUE)  # reorder factor
  )


# Create the plot
# Filter top 5 based on genes_afected
top5_df <- plot_df %>%
  arrange(desc(genes_afected)) %>%
  dplyr::slice(1:5)

# Plot only the top 5
genes_afected.1 <- ggplot(data = top5_df, 
                          aes(x = fct_reorder(top_SNP, genes_afected), 
                              y = genes_afected)) +
  geom_bar(stat = "identity", color = "black", fill = "red") +
  geom_text(aes(y = 5, label = formatted_label), 
            angle = 0, hjust = -0.1, size = 7) +
  ylab("Genes with the strongest association") +
  xlab("Marker") +
  ylim(0, max(top5_df$genes_afected) * 1.1) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip()


# pw1 / (pw2 + genes_afected.1)
# 
# B.2 <- plot_grid(pw2, genes_afected.1, rel_widths = c(1,.6),
#                  nrow = 1, labels = c('G','H'),hjust = 0, vjust = 1, label_size = 18)
# B <- plot_grid(pw1, B.2, nrow = 2, labels = c('F',''),hjust = 0, vjust = 1, label_size = 18)
# B
B <- plot_grid(pw1, genes_afected.1, nrow = 2, rel_heights = c(1,1))
B


######################

group="C"

region.old <- merged_df_regions_letters3[which(merged_df_regions_letters3$region == group),"regionRAW"]

reduced_df4plot <- merged_df_regions_letters[merged_df_regions_letters$region == region.old,]
reduced_df4plot <- reduced_df4plot[1,]
merged_df_regions_letters2_4plot <- merged_df_regions_letters3[which(merged_df_regions_letters3$region == group),]
merged_df_regions_letters2_4plot


reduced_df4plot2 <- sorghum.trans %>%
  filter(CHROM == reduced_df4plot$chr,
         POS>=reduced_df4plot$start_big, POS<= reduced_df4plot$end_big)

merged_df_regions_letters2_4plot
#this one is to plot the zoom in window
merged_df_regions_letters2_4plot2 <- sorghum.trans %>%
  filter(CHROM == merged_df_regions_letters2_4plot$chr[1], 
         POS>=merged_df_regions_letters2_4plot$start_small[1]-1000, POS<= merged_df_regions_letters2_4plot$end_small[1]+1000)

#to to barplots later on
merged_df_regions_letters2_4plot2.1 <- sorghum.trans %>%
  filter(CHROM == merged_df_regions_letters2_4plot$chr[1], 
         POS>=merged_df_regions_letters2_4plot$start_small[1], POS<= merged_df_regions_letters2_4plot$end_small[1])

#to to barplots later on
merged_df_regions_letters2_4plot2.2 <- sorghum.trans %>%
  filter(CHROM == merged_df_regions_letters2_4plot$chr[2], 
         POS>=merged_df_regions_letters2_4plot$start_small[2], POS<= merged_df_regions_letters2_4plot$end_small[2])

#to to barplots later on
merged_df_regions_letters2_4plot2.3 <- sorghum.trans %>%
  filter(CHROM == merged_df_regions_letters2_4plot$chr[3], 
         POS>=merged_df_regions_letters2_4plot$start_small[3], POS<= merged_df_regions_letters2_4plot$end_small[3])


merged_df_regions_letters2_4plot3 <- merged_df_regions_letters2_4plot2 %>%
  group_by(top_SNP)  %>% 
  summarise(genes_afected = sum(!is.na(gene)), .groups = 'drop')

merged_df_regions_letters2_4plot3.1 <- merged_df_regions_letters2_4plot2.1 %>%
  group_by(top_SNP)  %>% 
  summarise(genes_afected = sum(!is.na(gene)), .groups = 'drop')

merged_df_regions_letters2_4plot3.2 <- merged_df_regions_letters2_4plot2.2 %>%
  group_by(top_SNP)  %>% 
  summarise(genes_afected = sum(!is.na(gene)), .groups = 'drop')

merged_df_regions_letters2_4plot3.3 <- merged_df_regions_letters2_4plot2.3 %>%
  group_by(top_SNP)  %>% 
  summarise(genes_afected = sum(!is.na(gene)), .groups = 'drop')

merged_df_regions_letters2_4plot3.1
merged_df_regions_letters2_4plot3.2
merged_df_regions_letters2_4plot3.3

# get windows in the region to plot them
reduced_df4plot2.2 <- hotspot_dfS %>%
  filter(seqnames == reduced_df4plot$chr,
         start>=reduced_df4plot$start_big, end<= reduced_df4plot$end_big)

#order them
reduced_df4plot2.2 <- reduced_df4plot2.2[order(reduced_df4plot2.2$start),]
reduced_df4plot2.2$level <- 1:nrow(reduced_df4plot2.2)


# get the information from genes
Sbicolor_730_v5.1_primary_exons <- fread("Sbicolor_730_v5.1_primary_exons.csv", data.table = F)
#change names
Sbicolor_730_v5.1_primary_exons$type[which(Sbicolor_730_v5.1_primary_exons$type == "five_prime_UTR")] <- "5'UTR"
Sbicolor_730_v5.1_primary_exons$type[which(Sbicolor_730_v5.1_primary_exons$type == "three_prime_UTR")] <- "3'UTR"

head(Sbicolor_730_v5.1_primary_exons)
Sbicolor_730_v5.1_primary_exons_4plot2 <- Sbicolor_730_v5.1_primary_exons %>%
  filter(chr == reduced_df4plot$chr,
         start>=reduced_df4plot$start_big, end<= reduced_df4plot$end_big, type %in% c("5'UTR", "CDS", "3'UTR"))



####
#Calculate LD
#singleMarker <- merged_df_regions_letters2_4plot2$top_SNP[which(merged_df_regions_letters2_4plot2$num_SNPs ==
#                                                  max(merged_df_regions_letters2_4plot2$num_SNPs))]

merged_df_regions_letters2_4plot3.1
merged_df_regions_letters2_4plot3.2
merged_df_regions_letters2_4plot3.3 #select the one with more snps
singleMarker <- "Chr01_68265218" #Chr03_63696487 Chr03_63696436

plinkPath <- "plink"
#name <- paste0("SNP_maize",group,".txt")
#name
#plink_command <- paste0(plinkPath," --bfile vla_karla_811geno_filter2_maf_het --r2 --ld-window-r2 0 --ld-window-kb 1000 --ld-window 9999999  --extract ", name," --out ", name,"LD")
plink_command <- paste0(plinkPath," --bfile vla_karla_811geno_filter2_maf_het --r2 --ld-window-r2 0 --ld-window-kb 100000 --ld-window 999999  --ld-snp ", singleMarker , " --out ", singleMarker)

system(plink_command)


#load plink results
#LD.true <- fread(paste0(name,"LD.ld"),data.table = F)
LD.true <- fread(paste0(singleMarker,".ld"),data.table = F)

### this will tell us the LD between the reference marker (the one with more hits in group one) against the other markers in the other region
LD.true[LD.true$SNP_B %in% merged_df_regions_letters2_4plot3.1$top_SNP,]
merged_df_regions_letters2_4plot

LD.true_for_plot1 <- LD.true[which(LD.true$BP_B >= merged_df_regions_letters2_4plot$start_small[1]-5000 &
                                     LD.true$BP_B <= merged_df_regions_letters2_4plot$end_small[1]+5000),]

LD.true_for_plot2 <- LD.true[which(LD.true$BP_B >= merged_df_regions_letters2_4plot$start_small[3]-5000 &
                                     LD.true$BP_B <= merged_df_regions_letters2_4plot$end_small[3]+5000),]

max(reduced_df4plot2$POS)-min(reduced_df4plot2$POS)

focus_pos_bp   <- 68265218        # <- set your focal marker POS (in bp)
half_window_kb <- 25             # <- set half window size (in kb) around POS
half_window_bp <- half_window_kb * 1000

pw1 <- ggplot(reduced_df4plot2, aes(x = POS/1000000, y = nrow(reduced_df4plot2.2)/2)) +
  # Background rectangles
  # geom_rect(data = reduced_df4plot2.2,
  #           inherit.aes = FALSE,
  #           aes(xmin = start/1000000, xmax = end/1000000,
  #               ymin = level - 0.4, ymax = level + 0.4),
  #           fill = "lightblue", alpha = 0.5) +
  # Vertical lines
  geom_vline(xintercept = merged_df_regions_letters2_4plot$start_small/1000000, 
             color = "red", linewidth = 0.5) +
  geom_vline(xintercept = merged_df_regions_letters2_4plot$end_small/1000000, 
             color = "red", linewidth = 0.5) +
  # Gene structure: Exons (UTRs and CDS as boxes)
  geom_rect(data = Sbicolor_730_v5.1_primary_exons_4plot2,
            inherit.aes = FALSE,
            aes(xmin = start/1000000, xmax = end/1000000,
                ymin = 0, ymax = nrow(reduced_df4plot2.2)*.05,
                fill = type),
            color = "black") +  # black border around exons
  
  # Gene structure: Introns (lines between exons)
  geom_segment(data = Sbicolor_730_v5.1_primary_exons_4plot2 %>%
                 group_by(gene) %>%
                 arrange(start) %>%
                 mutate(next_start = lead(start),
                        next_end = lead(end)) %>%
                 filter(!is.na(next_start)),
               inherit.aes = FALSE,
               aes(x = end/1000000,
                   xend = next_start/1000000,
                   y = (nrow(reduced_df4plot2.2)*.05)/2,
                   yend = (nrow(reduced_df4plot2.2)*.05)/2),
               color = "black",
               linewidth = 1) +
  # Jitter points ON TOP
  geom_jitter(shape = 6, height = nrow(reduced_df4plot2.2)/4, width = 0, size = 5) +
  # Axis and labels
  scale_x_continuous(
    limits = c(min(reduced_df4plot$start_big/1000000), max(reduced_df4plot$end_big/1000000)),
    expand = c(0, 0),
    breaks = round(seq(from = min(reduced_df4plot$start_big/1000000), 
                       to = max(reduced_df4plot$end_big/1000000), 
                       by = 0.01), 2)
  )  +
  scale_fill_manual(values = c("5'UTR" = "skyblue", 
                               "3'UTR" = "grey", 
                               "CDS" = "blue")) +
  # ---- the zoom line (adjust focus_pos_bp / half_window_kb above) ----
coord_cartesian(xlim = c((focus_pos_bp - half_window_bp)/1e6,
                         (focus_pos_bp + half_window_bp)/1e6)) +
  labs(x = "Genomic Position (Mb)", y = "", fill = "Feature") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
pw1

### Zoom-in plot with gene models
# pw2 <- ggplot(LD.true_for_plot1[which(!LD.true_for_plot1$R2 == 1),], aes(x = BP_B/1000000, y = R2)) +
#   rasterize(geom_point(shape = 19, size = 5, color="grey"), dpi = 900, dev = "ragg_png") +
#   scale_x_continuous(
#     expand = c(0, 0),
#     breaks = seq(from = (merged_df_regions_letters2_4plot$start_small[1]-6000)/1000000, 
#                  to = (merged_df_regions_letters2_4plot$end_small[1]+6000)/1000000, 
#                  by = 0.001)
#   ) +
#   ylim(0,1) +
#   # Vertical lines
#   geom_vline(xintercept = merged_df_regions_letters2_4plot$start_small/1000000, 
#              color = "red", linewidth = 0.5) +
#   geom_vline(xintercept = merged_df_regions_letters2_4plot$end_small/1000000, 
#              color = "red", linewidth = 0.5) +
#   
#   # Gene structure: Exons (UTRs and CDS as boxes)
#   geom_rect(data = Sbicolor_730_v5.1_primary_exons_4plot2,
#             inherit.aes = FALSE,
#             aes(xmin = start/1000000, xmax = end/1000000,
#                 ymin = 0, ymax = .1,
#                 fill = type),
#             color = "black") +  # black border around exons
#   
#   # Gene structure: Introns (lines between exons)
#   geom_segment(data = Sbicolor_730_v5.1_primary_exons_4plot2 %>%
#                  group_by(gene) %>%
#                  arrange(start) %>%
#                  mutate(next_start = lead(start),
#                         next_end = lead(end)) %>%
#                  filter(!is.na(next_start)),
#                inherit.aes = FALSE,
#                aes(x = end/1000000,
#                    xend = next_start/1000000,
#                    y = .05,
#                    yend = 0.05),
#                color = "black",
#                linewidth = 4) +
#   # Jitter points ON TOP
#   geom_jitter(data = merged_df_regions_letters2_4plot2, aes(x= POS/1e6, y=.5) , shape = 6, height = .1, width = 0, size = 5) +
#   scale_fill_manual(values = c("5'UTR" = "skyblue", 
#                                "3'UTR" = "grey", 
#                                "CDS" = "blue")) +
#   
#   labs(x = "Genomic Position (Mb)", y = "LD", fill = "Feature") +
#   # ggnewscale::new_scale_fill() +  # This allows adding a *new* fill scale
#   # geom_rect(data = LD.true_for_plot,
#   #           aes(xmin = BP_B - 0.0001, xmax = BP_B + 0.0001,
#   #               ymin = nrow(reduced_df4plot2.2)*.90,
#   #               ymax = nrow(reduced_df4plot2.2)*.95,
#   #               fill = R2),
#   #           inherit.aes = FALSE) +
#   # scale_fill_gradient(low = "blue", high = "red", name = expression(R^2)) +
#   theme(legend.position = "inside", axis.text.x = element_text(angle = 15, hjust = 1), 
#         legend.position.inside = c(.3, .99),       # bottom right inside
#         legend.justification = c("right", "top")) +
#   coord_cartesian(xlim = c((merged_df_regions_letters2_4plot$start_small[1]-6000)/1e6,
#                            (merged_df_regions_letters2_4plot$end_small[1]+6000)/1e6), clip="on")
# 
# pw2



# Calculate shared y-axis limits
ymax <- max(
  max(merged_df_regions_letters2_4plot3.1$genes_afected, na.rm = TRUE),
  max(merged_df_regions_letters2_4plot3.2$genes_afected, na.rm = TRUE),
  max(merged_df_regions_letters2_4plot3.3$genes_afected, na.rm = TRUE)
)

# Modify the data to create new labels
# Prepare the data: format labels and reorder top_SNP
plot_df <- merged_df_regions_letters2_4plot3.1 %>%
  mutate(
    chr = str_extract(top_SNP, "\\d+"),
    pos = str_extract(top_SNP, "\\d+$"),
    formatted_label = paste0(as.integer(chr), ":", comma(as.numeric(pos))),
    top_SNP = fct_reorder(top_SNP, genes_afected, .desc = TRUE)  # reorder factor
  )

# Plot with rotated labels at y = 1
# Filter top 5 based on genes_afected
top5_df <- plot_df %>%
  arrange(desc(genes_afected)) %>%
  dplyr::slice(1:5)

# Plot only the top 5
genes_afected.1 <- ggplot(data = top5_df, 
                          aes(x = fct_reorder(top_SNP, genes_afected), 
                              y = genes_afected)) +
  geom_bar(stat = "identity", color = "black", fill = "red") +
  geom_text(aes(y = 5, label = formatted_label), 
            angle = 0, hjust = -0.1, size = 7) +
  ylab("Genes with the strongest association") +
  xlab("Marker") +
  ylim(0, max(top5_df$genes_afected) * 1.1) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip()

genes_afected.1


C <- plot_grid(pw1, genes_afected.1, nrow = 2, rel_heights = c(1,1))
C

#### Upsetplot

hotspots_sorghumRNA_geneInregions <- fread("hotspots_sorghumRNA_geneInregions.csv", data.table = F)

table(hotspots_sorghumRNA_geneInregions$region)

head(hotspots_sorghumRNA_geneInregions.working)


# Get all distinct gene-region pairs
binary_matrix <- hotspots_sorghumRNA_geneInregions %>%
  distinct(gene, region) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = region, values_from = value, values_fill = 0)

region_cols <- setdiff(names(binary_matrix), "gene")


# Columns that must be all zeros
cols_zero <- c("E","D","H","I","J","M","G","K","L","F")

# Base solution: rows where A, B, C are 1 AND all listed columns are 0
ABC <- binary_matrix %>%
  filter(A == 1, B == 1, C == 1) %>%
  filter(if_all(all_of(cols_zero), ~ .x == 0))

fwrite(ABC,"ABC_131genes.csv")
fwrite(binary_matrix,"All_131genes.csv")

intersections <- binary_matrix %>%
  rowwise() %>%
  mutate(combination = paste(sort(region_cols[which(c_across(all_of(region_cols)) == 1)]), collapse = "+")) %>%
  ungroup() %>%
  count(combination, sort = TRUE)


top_n <- 20  # or fewer if you want

upset <- ggplot(intersections %>% slice_max(n, n = top_n),
       aes(x = reorder(combination, -n), y = n)) +
  geom_col(fill = "green4", color = "black") +
  geom_text(aes(label = n), vjust = -0.3, size = 6) +
  labs(x = "Region Combination", y = "Number of Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



##################
top <- plot_grid(hotspot_plot.maize, hotspot_plot.sorghum.regions, nrow = 2, labels = c('A','B'),
                 align = "v", axis = "l"  ,hjust = 0, vjust = 1, label_size = 18)
bottom <- plot_grid(C, A, B, nrow = 1, labels = c('C','D','E'),hjust = 0, vjust = 1,  label_size = 18)

fancyplots <- plot_grid(top,bottom, ncol = 1, hjust = 0, rel_heights = c(1, 1), labels = c('',''), vjust = 1,  label_size = 18)
fancyplots

################

ggsave("Fig2_v2.svg", width = 12, height = 10) # to save text as editable
upset

ggsave(plot = upset,"Figsuup_upset.svg", width = 12, height = 5) # to save text as editable
