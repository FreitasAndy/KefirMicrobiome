# Kefir Microbiome
# By Anderson Freitas

library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(readxl)
library(dplyr)
library(magrittr)
library(microbiome)
library(microeco)
library(ggpubr)
library(corrplot)
theme_set(theme_bw())


Vol1 <- read_excel("Vol_Anderson_duplicata.xlsx", 
                   sheet = "Vol1")
Vol2 <- read_excel("Vol_Anderson_duplicata.xlsx", 
                   sheet = "Vol2")
Vol1.s <- Vol1[,c(9,8,4)]
Vol2.s <- Vol2[,c(9,8,4)]
vol3 <- rbind(Vol1.s, Vol2.s)
write.csv(vol3, file = "./vol_ok.csv")


# DADA2
path <- "./16S"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnRs[1:6])
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(300,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[4]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
taxa <- assignTaxonomy(seqtab.nochim, "./tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "./tax/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

samdf <- read_excel("./Map_T.xlsx", sheet = "16S")
samdf <- samdf[order(samdf$SampleID, decreasing = F),]
samdf <- as.data.frame(samdf)
row.names(samdf) <- c("S13", "S14", "S15", "S16", "S17", "S18")
otu <- as.data.frame(t(seqtab.nochim))
ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE), 
               sample_data(samdf), 
               tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

plot_richness(ps, x="Time", measures=c("Shannon", "Observed"), color="Time")
# Diversity decreases over time
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Time", title="Beta Diversity")
# Samples are completelly different over time
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:100]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Time", fill="Genus") + facet_wrap(~Time, scales="free_x")

# alpha diversity
summarize_phyloseq(ps)
inputR <- rarefy_even_depth(ps, sample.size = 16983)
richness <- estimate_richness(inputR)
alpha = cbind(richness, meta(inputR))
colMeans(alpha[,c(1,8)])
ababaxe <- reshape2::melt(alpha[,c(1,8,10,11)])
ababax  <- reframe(alpha, by = Time)
ababaxe %>%
  ggplot(aes(x=Time, y=value, group = variable, color = variable)) +
  geom_point() +
  geom_smooth(method=lm)+
  labs(y = "Diversity Values", x = "Maturation Period")
dev.print(tiff, "./Alpha_Diversity_Bacteria.tiff", compression = "lzw", res=600, height=5, width=8, units="in")

# Beta diversity
input.clr = microbiome::transform(ps, "clr")
df        = as(sample_data(input.clr), "data.frame")
ds        = phyloseq::distance(input.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Time, data = df, permutations = 999)
permanova
input_ord = ordinate(input.clr, "PCoA" , "euclidean") 
p4 = plot_ordination(input.clr, input_ord, color = "Time")
p1.bac = p4 + geom_point(size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = -5, y = -2.5, hjust = 0.2 , 
           label = bquote('PERMANOVA:'~R^2~'= 0.66  |  p = 0.067'), size = 3)+
  theme_bw()
p1.bac
dev.print(tiff, "./Beta_Diversity_Bacteria.tiff", compression = "lzw", res=600, height=6, width=8, units="in")

# Genera distribution
dataset <- file2meco::phyloseq2meco(ps)
dataset$tax_table %<>% tidy_taxonomy
dataset$cal_abund()
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 50)
aaaa <- t1$plot_bar(others_color = "grey70", facet = "Time", xtext_keep = FALSE, legend_text_italic = T) +
  theme(legend.position = "none", axis.title.y = element_text(size = 12))
bbbb <- t1$plot_donut(label = F, legend_text_italic = T) +
  theme(legend.position = "none") 
ggpubr::ggarrange(bbbb, aaaa, nrow = 2, legend = NULL, common.legend = F)
dev.print(tiff, "./Phyla_abundance_Bacteria.tiff", compression = "lzw", res=600, height=6, width=9, units="in")

#Differential abundance
# Need it?

# FAPROTAX
t1 <- trans_func$new(dataset)
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = TRUE)
genes = cbind(as.data.frame(t1$res_spe_func_perc),sample_data(ps))
colnames(genes)
allgenes <- genes[-c(3,4,5)] %>% 
  group_by(Time) %>%
  summarise_all("mean")
genes2 <- reshape2::melt(allgenes)
genes2 %>%
  ggplot(aes(x=Time, y=value, group = variable, color = variable)) +
  geom_line(size = 3) +
  labs(y = "Function Abundance (%)", x = "Maturation Period")
dev.print(tiff, "./Gene_abundance_Bacteria.tiff", compression = "lzw", res=600, height=6, width=9, units="in")


###### Heatmap
ps
psf

agg.ps <- aggregate_rare(ps, level = "Genus", prevalence = 1/100, detection = 1/100)
agg.psf <- aggregate_rare(psf, level = "Genus", prevalence = 1/100, detection = 1/100)
otu.ps <- otu_table(agg.ps)
otu.psf <- otu_table(agg.psf)
Vol1.samples <- Vol1[,c(1,8,4)]
Vol2.samples <- Vol2[,c(1,8,4)]
vol3 <- rbind(Vol1.samples, Vol2.samples)
colnames(vol3) <- c("SampleID", "Compound", "Area")
vol4.samples <- dcast(data = vol3,formula = SampleID~Compound,fun.aggregate = sum,value.var = "Area")
vol5.f <- cbind(vol4.samples, t(otu.psf))
write.csv(vol5.f, file = "./vol.csv")
write.csv(otu.ps, file = "./otu_bac.csv")
df <- vol_6
correlation_df<-cor(df[,2:15], df[,16:32], method="spearman", use="pairwise.complete.obs")
corrplot(correlation_df,
         method="color")
heatmap(correlation_df)
pheatmap::pheatmap(correlation_df, cluster_rows = T, cluster_cols = T, show_rownames = TRUE, show_colnames = TRUE)



c_df <- Hmisc::rcorr(cor(correlation_df), type='spearman')
corrplot(corr=c_df$r[2:15, 16:32], p.mat=c_df$P[2:15, 16:32], sig.level=0.05, 
         method='color', diag=FALSE, addCoef.col=1, type='upper', insig='blank',
         number.cex=.8)
