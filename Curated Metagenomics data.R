#BiocManager::install("curatedMetagenomicData")
library(curatedMetagenomicData)
library(tidyverse)
library(phyloseq)
?combined_metadata
combined_metadata[1, ]
View(combined_metadata)
View(combined_metadata$disease)
table(combined_metadata$disease == 'healthy')
table(combined_metadata$country)
table(combined_metadata$body_site)
#read depth of samples across all studies
dsranking <- combined_metadata %>% 
  as.data.frame() %>% 
  group_by(dataset_name) %>%
  summarize(mediandepth = median(number_reads) / 1e6) %>%
  mutate(dsorder = rank(mediandepth)) %>%
  arrange(dsorder)
dsranking
combined_metadata %>%
  as.data.frame() %>%
  mutate(ds = factor(combined_metadata$dataset_name, levels = dsranking$dataset_name)) %>%
  ggplot(aes(ds, number_reads / 1e6, fill = body_site)) +
  geom_boxplot() +
  theme(axis.title.x = element_text(angle=45, hjust = 1)) +
  labs(x="Dataset", y="Read Depth (millions)")

#?curatedMetagenomicData
#Accessing datasets
loman.eset = LomanNJ_2013.metaphlan_bugs_list.stool()
yes
loman <- curatedMetagenomicData("LomanNJ_2013.metaphlan_bugs_list.stool", dryrun = FALSE)
loman
loman.eset <- loman[[1]]
loman.eset

oral <- c("BritoIL_2016.metaphlan_bugs_list.oralcavity",
          "Castro-NallarE_2015.metaphlan_bugs_list.oralcavity")
esl <- curatedMetagenomicData(oral, dryrun = FALSE)
esl
esl[[1]]
esl[[2]]
curatedMetagenomicData("*metaphlan_bugs_list.stool*", dryrun = TRUE)
#merging datasets
eset <- mergeData(esl)
eset
#using expressionset object
experimentData(loman.eset)
head(pData(loman.eset))
exprs(loman.eset)[1:6, 1:5] #first 6 rows and 5 columns
#Estimating Absolute Raw Count Data
loman.counts = sweep(exprs(loman.eset), 2, loman.eset$number_reads / 100, "*")
loman.counts = round(loman.counts)
loman.counts[1:6, 1:5]
loman.eset2 = curatedMetagenomicData("LomanNJ_2013.metaphlan_bugs_list.stool",
                                     counts = TRUE, dryrun = FALSE)[[1]]
all.equal(exprs(loman.eset2), loman.counts)
#E.coli prevalence
grep("coli", rownames(loman.eset), value = TRUE)
x = exprs(loman.eset)[grep("s__Escherichia_coli$", rownames(loman.eset)), ]
summary(x)
hist(x, xlab = "Relative Abundance", main = "Prevalence of E. Coli", 
     breaks = "FD")
#Taxonomy-Aware Analysis using phyloseq
loman.pseq = ExpressionSet2phyloseq(loman.eset)
loman.pseq
loman.tree <- ExpressionSet2phyloseq(loman.eset, phylogenetictree = TRUE)
wt = UniFrac(loman.tree, weighted = TRUE, normalized = FALSE, 
             parallel = FALSE, fast = TRUE)
plot(hclust(wt), main = "Weighted UniFrac distances")

#Components of a phyloseq object
otu_table(loman.pseq)[1:6, 1:5]
sample_data(loman.pseq)[1:6, 1:5]
head(tax_table(loman.pseq))

#Subsetting/Pruning
rank_names(loman.pseq)
subset_taxa(loman.pseq, !is.na(Species))
subset_taxa(loman.pseq, is.na(Class) & !is.na(Phylum))
loman.bd = subset_taxa(loman.pseq, Phylum == "Bacteroidetes")
head(taxa_names(loman.bd))

#Advanced Pruning
keepotu = genefilter_sample(loman.pseq, filterfun_sample(topp(0.05)), A = 5)
summary(keepotu)
subset_taxa(loman.pseq, keepotu)

#Taxonomy Heatmap
loman.filt = subset_taxa(loman.pseq, keepotu & !is.na(Strain))
plot_heatmap(loman.filt, method = "PCoA", distance = "bray")

#Taxonomy Histogram
loman.sp = subset_taxa(loman.pseq, !is.na(Species) & is.na(Strain))
par(mar = c(20, 4, 0, 0) + 0.15)#increase margin size on the bottom
barplot(sort(taxa_sums(loman.sp), TRUE)[1:20] / nsamples(loman.sp),
        ylab = "Total counts", las = 2)
alphas = c("Shannon", "Simpson", "InvSimpson")
plot_richness(loman.sp, "stool_texture", measures = alphas)
pairs(estimate_richness(loman.sp, measures = alphas))#pairs tests as in qiime2

#Beta Diversity / Dissimilarity Clustering
mydist = distance(loman.sp, method = "bray")
myhclust = hclust(mydist)
plot(myhclust, main = "Bray-Curtis Dissimilarity", 
     method = "ward.D", xlab = "Samples", sub = "")

#Ordination Analysis
ordinated_taxa = ordinate(loman.sp, method = "PCoA", distance = "bray")
plot_ordination(loman.sp, ordinated_taxa, color = "stool_texture",
                title = "Bray-Curtis Principal Coordinates Analysis")
plot_scree(ordinated_taxa, title = "Screeplot")

#Browsing ExperimentHub
library(ExperimentHub)
eh = ExperimentHub()
myquery = query(eh, "curatedMetagenomicData")
myquery
View(mcols(myquery))
write.csv(mcols(myquery), file = "cureatedMetagenomicData_allrecords.csv")

#Advanced searching of ExperimentHub
head(myquery$tags)
grep("(?=HMP)(?=.*metaphlan)", myquery$title, perl = TRUE, val = TRUE)
grep("(?=.*stool)(?=.*metaphlan)", myquery$title, perl = TRUE, val = TRUE)
#a list of expression set objects containing all of the Loman products
(lomannames = grep("LomanNJ_2013", myquery$title, perl = TRUE, val = TRUE))
(idx = grep("LomanNJ_2013", myquery$title, perl = TRUE))
loman.list = lapply(idx, function(i){
  return(myquery[[i]])
})
names(loman.list) = lomannames
loman.list