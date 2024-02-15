# script to get annotation summary statistics from a few files
library(ggplot2)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)
gff = args[1]
introns = args[2]
genestats = args[3]
peps = args[4]
species = args[5]

#setwd("H:/KrabbLab/Genome_annotations/Mmel/test")
#gff = "Caur.longest_isoform.simple.gff"
#introns = "introns.simple.txt"
#genestats = "Caur.genestats.txt"
#peps = "Caur.longest_isoform.protein.lens.txt"
#species="Caur"

data <- read.table(gff)
data$V10 <- data$V5 - data$V4
exons <- subset(data, data$V3=="exon")
mean_exon_len <- mean(exons$V10)
median_exon_len <- median(exons$V10)
stddev_exon_len <- sd(exons$V10)

in_data <- read.table(introns)
in_data$V10 <- in_data$V5 - in_data$V4
mean_intron_len <- mean(in_data$V10)
median_intron_len <- median(in_data$V10)
stddev_intron_len <- sd(in_data$V10)

gs_data <- read.table(genestats, header = TRUE)
mean_gene_len <- mean(gs_data$transcript_length)
median_gene_len <- median(gs_data$transcript_length)
stddev_gene_len <- sd(gs_data$transcript_length)

mean_exons <- mean(gs_data$exons)
median_exons <- median(gs_data$exons)
stddev_exons <- sd(gs_data$exons)

ngenes_with_one_exon <- nrow(subset(gs_data, exons==1))

mean_introns <- mean(gs_data$introns)
median_introns <- median(gs_data$introns)
stddev_introns <- sd(gs_data$introns)

pep_data <- read.table(peps)
mean_pep_len <- mean(pep_data[,ncol(pep_data)])
median_pep_len <- median(pep_data[,ncol(pep_data)])
stddev_pep_len <- sd(pep_data[,ncol(pep_data)])

zz <- file(paste0(species, ".EVM.summaryStats.txt"), open = "wt")
sink(zz, type = "message")
message(paste0("number of genes: ", nrow(gs_data)))
message("")
message(paste0("mean gene length: ", round(mean_gene_len, digits=2)))
message(paste0("median gene length: ", median_gene_len))
message(paste0("std dev gene length: ", round(stddev_gene_len, digits=2)))
message("")
message(paste0("mean number of exons: ", round(mean_exons, digits=2)))
message(paste0("median number of exons: ", median_exons))
message(paste0("std dev number of exons: ", round(stddev_exons, digits=2)))
message(paste0("number of genes with one exon: ", ngenes_with_one_exon))
message("")
message(paste0("mean exon length: ", round(mean_exon_len, digits=2)))
message(paste0("median exon length: ", median_exon_len))
message(paste0("std dev exon length: ", round(stddev_exon_len, digits=2)))
message("")
message(paste0("mean number of introns: ", round(mean_introns, digits=2)))
message(paste0("median number of introns: ", median_introns))
message(paste0("std dev number of introns: ", round(stddev_introns, digits=2)))
message("")
message(paste0("mean intron length: ", round(mean_intron_len, digits=2)))
message(paste0("median intron length: ", median_intron_len))
message(paste0("std dev intron length: ", round(stddev_intron_len, digits=2)))
message("")
message(paste0("mean peptide length: ", round(mean_pep_len, digits=2)))
message(paste0("median peptide length: ", median_pep_len))
message(paste0("std dev peptide length: ", round(stddev_pep_len, digits=2)))

closeAllConnections()

gene_len_p <- ggplot(data=gs_data, aes(x=transcript_length/1000)) +
  xlab("gene length (kb)") +
  ylab("genes") +
  geom_histogram(position="identity", color="white",
                 bins=50, boundary=0)+ 
  theme_minimal()

gene_len_p_crop <- gene_len_p + xlim(c(0,20))


pep_len_p <- ggplot(data=pep_data, aes(x=pep_data[,ncol(pep_data)])) +
  xlab("protein length (aa)") +
  ylab("proteins") +
  geom_histogram(position="identity", color="white",
                 bins=50, boundary=0)+ 
  theme_minimal()

pep_len_p_crop <- pep_len_p + xlim(c(0,3000))


intron_len_p <- ggplot(data=in_data, aes(x=V10/1000)) +
  xlab("intron length (kb)") +
  ylab("introns") +
  geom_histogram(position="identity", color="white",
                 bins=50, boundary=0)+ 
  theme_minimal()

intron_len_p_crop <- ggplot(data=in_data, aes(x=V10)) +
  xlab("intron length (bp)") +
  ylab("introns") +
  geom_histogram(position="identity", color="white",
                 bins=50, boundary=0)+ 
  theme_minimal() +
  xlim(c(0,1500))

exon_len_p <- ggplot(data=exons, aes(x=V10/1000)) +
  xlab("exon length (kb)") +
  ylab("exons") +
  geom_histogram(position="identity", color="white",
                 bins=50, boundary=0)+ 
  theme_minimal()

exon_len_p_crop <- ggplot(data=exons, aes(x=V10)) +
  xlab("exon length (bp)") +
  ylab("exons") +
  xlim(c(0,750)) +
  geom_histogram(position="identity", color="white",
                 bins=50, boundary=0)+ 
  theme_minimal()

p_grid <- plot_grid(gene_len_p, gene_len_p_crop,
          exon_len_p, exon_len_p_crop,
          intron_len_p, intron_len_p_crop,
          pep_len_p, pep_len_p_crop,
          ncol=2)

ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)
ggsave(p_grid, width = 11, height = 8, device = "pdf", units = "in",
       filename = paste0(species, ".EVM.summaryStats.pdf"))


