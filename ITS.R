# Kefir Microbiome
# By Anderson Freitas

library(dada2); packageVersion("dada2")

path <- "./ITS"
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
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(300,280),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])
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
taxa <- assignTaxonomy(seqtab.nochim, "./tax/sh_general_release_dynamic_04.04.2024.fasta", multithread=TRUE)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

samdf <- read_excel("./Map_T.xlsx", sheet = "ITS")
#samdf <- samdf[order(samdf$SampleID, decreasing = F),]
samdf <- as.data.frame(samdf)
row.names(samdf) <- c("S81",  "S141", "S79",  "S139", "S77",  "S137")
otu <- as.data.frame(t(seqtab.nochim))
psf <- phyloseq(otu_table(otu, taxa_are_rows=TRUE), 
               sample_data(samdf), 
               tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(psf))
names(dna) <- taxa_names(psf)
psf <- merge_phyloseq(psf, dna)
taxa_names(psf) <- paste0("ASV", seq(ntaxa(psf)))
psf

plot_richness(psf, x="Time", measures=c("Shannon", "Observed"), color="Time")
# Diversity decreases over time
psf.prop <- transform_sample_counts(psf, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(psf.prop, method="NMDS", distance="bray")
plot_ordination(psf.prop, ord.nmds.bray, color="Time", title="Beta Diversity")
# Samples are completelly different over time
top20 <- names(sort(taxa_sums(psf), decreasing=TRUE))[1:100]
psf.top20 <- transform_sample_counts(psf, function(OTU) OTU/sum(OTU))
psf.top20 <- prune_taxa(top20, psf.top20)
plot_bar(psf.top20, x="Time", fill="Genus") + facet_wrap(~Time, scales="free_x")

# alpha diversity
summarize_phyloseq(psf)
inputR <- rarefy_even_depth(psf, sample.size = 6550)
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
dev.print(tiff, "./Alpha_Diversity_Fungi.tiff", compression = "lzw", res=600, height=5, width=8, units="in")

# Beta diversity
input.clr = microbiome::transform(psf, "clr")
df        = as(sample_data(input.clr), "data.frame")
ds        = phyloseq::distance(input.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Time, data = df, permutations = 999)
permanova
input_ord = ordinate(input.clr, "PCoA" , "euclidean") 
p4 = plot_ordination(input.clr, input_ord, color = "Time")
p1.bac = p4 + geom_point(size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = 3, y = 5, hjust = 0.2 , 
           label = bquote('PERMANOVA:'~R^2~'= 0.40  |  p = 0.73'), size = 3)+
  theme_bw()
p1.bac
dev.print(tiff, "./Beta_Diversity_Fungi.tiff", compression = "lzw", res=600, height=6, width=8, units="in")

# Genera distribution
dataset <- file2meco::phyloseq2meco(psf)
dataset$tax_table %<>% tidy_taxonomy
dataset$cal_abund()
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 50)
aaaa <- t1$plot_bar(others_color = "grey70", facet = "Time", xtext_keep = FALSE, legend_text_italic = T) +
  theme(legend.position = "none", axis.title.y = element_text(size = 12))
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 50, groupmean = "Time")
bbbb <- t1$plot_donut(label = F, legend_text_italic = T) +
  theme(legend.position = "none") 
ggpubr::ggarrange(bbbb, aaaa, nrow = 2, legend = NULL, common.legend = F)
dev.print(tiff, "./Phyla_abundance_Fungi.tiff", compression = "lzw", res=600, height=6, width=9, units="in")

#Differential abundance
# Need it?

 # FAPROTAX
t1 <- trans_func$new(dataset)
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = TRUE)
genes = cbind(as.data.frame(t1$res_spe_func_perc),sample_data(psf))
colnames(genes)
allgenes <- genes[-c(2,12)] %>% 
  group_by(Time) %>%
  summarise_all("mean")
genes2 <- reshape2::melt(allgenes)
genes2 %>%
  ggplot(aes(x=Time, y=value, group = variable, color = variable)) +
  geom_line(size = 3) +
  geom_jitter() +
  labs(y = "Function Abundance (%)", x = "Maturation Period")
dev.print(tiff, "./Gene_abundance_Fungi.tiff", compression = "lzw", res=600, height=6, width=9, units="in")


#Volateis
colnames(vol3) <- c("Day", "Compound", "Area")
vol3 %>%
  ggplot(aes(x=Day, y=Area, group = Compound, color = Compound)) +
  geom_line(size = 3) +
  geom_jitter() +
  labs(y = "Function Abundance (%)", x = "Maturation Period")
library(maditr)
vol4 <- dcast(data = vol3,formula = Day~Compound,fun.aggregate = NULL,value.var = "Area")
vol5 <- melt(vol4, id.vars = "Day")
vol5 %>%
  ggplot(aes(x=Day, y=value, group = variable, color = variable)) +
  geom_line(size = 2) +
  theme_bw()+
  labs(y = "Area Percentage (%)", x = "Maturation Period")+
  scale_color_manual(values = c("#72acba",
                               "#9440dc",
                               "#87d84d",
                               "#543c99",
                               "#cdb651",
                               "#cb53a9",
                               "#89d8a7",
                               "#ad3c50",
                               "#578342",
                               "#938cd2",
                               "#ce5f33",
                               "#4b3450",
                               "#cd9d99",
                               "#564b2f"))



