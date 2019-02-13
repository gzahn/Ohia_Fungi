# ------------------------------------------------------------------#
# Ohia Fungal Endophyte Analyses - Depends on "Process_Raw_Reads.R" #
# Author: Geoffrey Zahn                                             #
# ------------------------------------------------------------------#


# Packages ####

# source("https://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")


# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library(DESeq2)
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(gridExtra)
library(ggpubr)
library(dada2); packageVersion("dada2")
library(stringr)
library(purrr)
library(corncob)
library(deseq2)
library(vegan)
library(dplyr)
library(reshape2)
library(ade4)
library(MASS)
library(ggbiplot)

# Load data ####
ps = readRDS(file = "./Output/clean_phyloseq_object.RDS")
meta=read.csv("metadata.csv")

# Remove empty ESVs and singletons and doublets ####
ps <- prune_taxa(colSums(otu_table(ps)) > 2, ps)


# Remove non-fungi ####
ps <- subset_taxa(ps, Kingdom == "k__Fungi")

# Remove NA values ####
ps = subset_samples(ps, Taxon != "#N/A")
# Remove any newly-empty samples ####
ps <- prune_samples(sample_sums(ps) != 0, ps)

sample_names(ps)

# Sanity check
which(sample_sums(ps) == 0)
plot(taxa_sums(ps))
plot(sample_sums(ps))


# Merge Samples by . . .    ####
psm = merge_samples(ps, "Taxon")

# repair merged values
sample_data(psm)$Taxon <- levels(sample_data(ps)$Taxon)
sample_data(psm)$Abaxial_Surface = plyr::mapvalues(sample_data(psm)$Abaxial_Surface, from = c(1,2), to=levels(sample_data(ps)$Abaxial_Surface) )
sample_data(ps)$Collection_Site


head(sample_data(psm))
 # Maybe?
# relative abundance transformation ####
psm_ra <- transform_sample_counts(psm, function(x) x / sum(x) )
ps_ra <- transform_sample_counts(ps, function(x) x / sum(x) )


# merge again by site
# Merge Samples by . . .    ####
psm_site = merge_samples(ps, "Collection_Site")

# repair merged values
sample_data(psm_site)$Collection_Site <- levels(sample_data(ps)$Collection_Site)
# sample_data(psm_site)$Abaxial_Surface = plyr::mapvalues(sample_data(psm_site)$Abaxial_Surface, from = c(1,2), to=levels(sample_data(ps)$Abaxial_Surface) )

head(sample_data(psm_site))

# relative abundance transformation ####
psm_ra_site <- transform_sample_counts(psm_site, function(x) x / sum(x) )


# taxonomy table
tax = as.data.frame(ps_ra@tax_table)



# OTU Richness per tree
df = as.data.frame(otu_table(ps))
df_pa = decostand(df,method = "pa")
tree_richness = rowSums(df_pa)

# write.csv(data.frame(TreeID = names(tree_richness), Richness = tree_richness), "./Output/Tree_Richness.csv", row.names = FALSE)  
  
# heatmap of fungal orders for each sample
ps_order = tax_glom(ps, "Order")
ps_order = merge_samples(ps_order, "Taxon")
main_orders = (which(colSums(ps_order@otu_table) > 1000))

df_order = as.data.frame(decostand(ps_order@otu_table,method = "total"))
orderlabs = ps_order@tax_table@.Data[,"Order"]
orderlabs = map(strsplit(orderlabs, "o__"),2)
orderlabs = map(strsplit(as.character(orderlabs), "_ord_"),1)
orderlabs = orderlabs[main_orders]
df_order = df_order[,main_orders]

 png("./Output/Heatmap_of_Order_by_Taxon.png",width = 8008, height = 8008, res = 300)
heatmap(as.matrix(df_order), labCol = orderlabs, col = gray.colors(50), Colv = NA, margins=c(10,16),
        cexRow = 2,cexCol = 2)
 dev.off()

 
 
 
# heatmap of fungal genera for each sample
ps_genus = tax_glom(ps, "Genus")
df_genus = as.data.frame(otu_table(ps_genus))
genuslabs = ps_genus@tax_table@.Data[,"Genus"]
genuslabs = map(strsplit(as.character(genuslabs), "g__"),2)

 png("./Output/Heatmap_of_Genus_by_Tree.png",width = 8008, height = 8008, res = 300)
 heatmap(as.matrix(df_genus), labCol = genuslabs, col = gray.colors(50), Colv = NA, margins=c(10,8),
         cexRow = 2)
 dev.off()



# Bar plot of relative abundance, split by taxon ####

plot_bar(psm_ra, fill = "Phylum",x="Taxon") +
  geom_bar(stat = "identity") + coord_flip()  + #facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Ohia Taxon",y="Relative Abundance")
# ggsave("./Output/BarPlot_Fungal_Phylum_by_Tree.png", height = 8, width = 12, dpi=300)

# Same, but merged by site
plot_bar(psm_ra_site, fill = "Phylum",x="Collection_Site") +
  geom_bar(stat = "identity") + coord_flip()  + #facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Collection Site",y="Relative Abundance")
ggsave("./Output/BarPlot_Fungal_Phylum_by_Site.png", dpi = 300)

plot_bar(psm_ra_site, fill = "Class",x="Collection_Site") +
  geom_bar(stat = "identity") + coord_flip()  + #facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Collection Site",y="Relative Abundance")
ggsave("./Output/BarPlot_Fungal_Class_by_Site.png", dpi = 300, height = 8, width = 12)



plot_bar(psm_ra, fill = "Phylum",x="Taxon") +
  geom_bar(stat = "identity") + coord_flip()  + facet_wrap(~(sample_data(psm_ra)$Abaxial_Surface)) +
  labs(x="Site",y="Relative Abundance") + lims(y=c(0,1)) + theme_bw()
ggsave("./Output/BarPlot_Fungal_Phylum_by_Tree_and_Surface.png", height = 8, width = 12,dpi=300)

names(sample_data(ps))

# Summarise taxa and plot ####
source("./summarize_taxa_Joey711.R")
order_summary = summarize_taxa(ps,"Order","Collection_Site")
# Remove NAs
order_summary = order_summary[order_summary$Order != "NA",]
# reorder
setorder(order_summary, -meanRA)
# plot
ggplot(order_summary, aes(x= meanRA, y=Order, facet=Collection_Site)) +
  geom_point() + geom_errorbarh(aes(xmin=meanRA-sdRA,xmax=meanRA+sdRA)) + labs(x="Relative Abundance") +
  facet_grid(~Collection_Site) +
  theme_bw() + theme(panel.spacing = unit(2, "lines"))
ggsave("./Output/Summarized_taxa_Order_by_Site.png", height = 8, width = 10, dpi = 300)

order_summary = summarize_taxa(ps,"Order","Taxon")
# Remove NAs
order_summary = order_summary[order_summary$Order != "NA",]
# reorder
setorder(order_summary, -meanRA)
# plot
ggplot(order_summary, aes(x= meanRA, y=Order, facet=Taxon)) +
  geom_point() + geom_errorbarh(aes(xmin=meanRA-sdRA,xmax=meanRA+sdRA)) + labs(x="Relative Abundance") +
  facet_grid(~Taxon) +
  theme_bw() + theme(panel.spacing = unit(0.5, "lines"))
ggsave("./Output/Summarized_taxa_Order_by_Taxon.png", height = 8, width = 10, dpi = 300)

order_summary = summarize_taxa(ps,"Order","Elevation")
# Remove NAs
order_summary = order_summary[order_summary$Order != "NA",]
# reorder
setorder(order_summary, -meanRA)
# plot
ggplot(order_summary, aes(x= meanRA, y=Order, facet=Elevation)) +
  geom_point() + geom_errorbarh(aes(xmin=meanRA-sdRA,xmax=meanRA+sdRA)) + labs(x="Relative Abundance") +
  facet_grid(~Elevation) +
  theme_bw() + theme(panel.spacing = unit(0.5, "lines"))
ggsave("./Output/Summarized_taxa_Order_by_Elevation.png", height = 8, width = 10, dpi = 300)




# Taxonomy table ... frequency of each taxa (number of samples each taxon is found in) ####
taxtable = as.data.frame(table(tax_table(ps)))
write.csv(taxtable[order(taxtable$Freq, decreasing = TRUE),],
          "./Output/Most_Widespread_Taxa.csv", row.names = FALSE, quote = FALSE)


# Ordination ####
# Phyloseq method trying several ordination methods
NMDS = ordinate(ps_ra, method = "NMDS")
DCA = ordinate(ps_ra, method = "DCA")
CCA = ordinate(ps_ra, method = "CCA")
RDA = ordinate(ps_ra, method = "RDA")
MDS = ordinate(ps_ra, method = "MDS")
PCoA = ordinate(ps_ra, method = "PCoA")

# Plot them  ####
nmdsplot=plot_ordination(ps_ra, NMDS, color = "Taxon",shape = "Abaxial_Surface") + theme_minimal() + labs(color="Taxon",shape="Leaf Type") + ggtitle("NMDS")
dcaplot=plot_ordination(ps_ra, DCA, color = "Taxon",shape = "Abaxial_Surface") + theme_minimal() + labs(color="Taxon",shape="Leaf Type") + ggtitle("DCA") + theme(legend.position = "none")
ccaplot=plot_ordination(ps_ra, CCA, color = "Taxon",shape = "Abaxial_Surface") + theme_minimal() + labs(color="Taxon",shape="Leaf Type") + ggtitle("CCA")+ theme(legend.position = "none")
rdaplot=plot_ordination(ps_ra, RDA, color = "Taxon",shape = "Abaxial_Surface") + theme_minimal() + labs(color="Taxon",shape="Leaf Type") + ggtitle("RDA")+ theme(legend.position = "none")
mdsplot=plot_ordination(ps_ra, MDS, color = "Taxon",shape = "Abaxial_Surface") + theme_minimal() + labs(color="Taxon",shape="Leaf Type") + ggtitle("MDS")+ theme(legend.position = "none")
pcoaplot=plot_ordination(ps_ra, PCoA, color = "Taxon",shape = "Abaxial_Surface") + theme_minimal() + labs(color="Taxon",shape="Leaf Type") + ggtitle("PCoA")+ theme(legend.position = "none")

 png(filename = "./Output/Ordination_Plots_taxon-and-leaftype.png",height = 1000,width = 1000,res = 100)
 grid.arrange(nmdsplot,dcaplot,ccaplot,rdaplot,mdsplot,pcoaplot)
 dev.off()

 ggsave(nmdsplot + stat_ellipse(alpha = 0.5), filename = "./Output/NMDS_taxon-and-leaftype_with_ellipses.png", dpi = 300)

#permanova ####
adonis.1 = adonis(otu_table(ps_ra) ~ sample_data(ps_ra)$Taxon +
                          sample_data(ps_ra)$Abaxial_Surface + 
                          sample_data(ps_ra)$Collection_Site)
adonis.2 = adonis(otu_table(ps_ra) ~ sample_data(ps_ra)$Taxon *
                    sample_data(ps_ra)$Abaxial_Surface * 
                    sample_data(ps_ra)$Collection_Site)

sink("./Output/Permanova_tables.txt")
print("additive model", quote = FALSE)
(adonis.1)
print("", quote = FALSE)
print("", quote = FALSE)
print("", quote = FALSE)
print("multiplicative model", quote = FALSE)
adonis.2
sink(NULL)

# Post-hoc test...where do differences manifest among taxa and sites?

# CCA plot and principle component biplot
CCA = cca(otu_table(ps_ra) ~ sample_data(ps_ra)$Taxon + sample_data(ps_ra)$Collection_Site)
plot(CCA, correlation = TRUE)
text(CCA, display="sites")

dim(df)
dim(dat)
dat = df[,colSums(df) > 0]
princomp = prcomp((dat),scale. = TRUE)
pc = as.data.frame(princomp$x[,1:2])

outliers = c(row.names(pc[which(pc$PC1 < 1),]),row.names(pc[which(pc$PC2 < .5),]),row.names(pc[which(pc$PC2 > 1),]))
outliers = unique(outliers)

sample_data(ps_ra)$outliers <- row.names(sample_data(ps_ra)) %in% outliers
meta$outliers <- row.names(sample_data(ps_ra)) %in% outliers
adonis.3 = adonis(otu_table(ps_ra) ~ sample_data(ps_ra)$outliers)

samdat = as.data.frame(sample_data(ps_ra))
mod2 = glm(outliers ~ Abaxial_Surface + Taxon + ALTITUDE + Collection_Site, family=binomial(link='logit'), data=meta)
mod1 = aov(outliers ~ Abaxial_Surface, data=meta)
# plot(mod1)

sink("./Output/Outlier_Predictors.txt")
summary(mod2)
sink(NULL)


sample_data(ps_ra)
sink(file = "./Output/Stat-Tests.txt", append = TRUE)
print("PermANOVA -- Additive")
adonis.1
print("PermANOVA -- Interactive")
adonis.2
print("PermANOVA -- Outliers")
adonis.3
sink(NULL)

# Find outlier trees
plot(princomp$rotation)
princomp$rotation[,2] > .4


eigenvals(princomp)


AIC(adonis.1,adonis.2)
# Mantel Test ####

spatial.dist = dist(cbind(ps_ra@sam_data$LONG, ps_ra@sam_data$LAT))
comm.dist =   vegdist(as.matrix(ps_ra@otu_table))

mantel.test = mantel.rtest(spatial.dist, comm.dist, nrepet = 9999)

 sink(file = "./Output/Stat-Tests.txt", append = TRUE)
 print("MANTEL TEST")
 mantel.test
 sink(NULL)

 png("./Output/Mantel_Plot.png")
 plot(mantel.test)
 dev.off()

# Rarefaction analyses ####
rarecurve(ps@otu_table,label = FALSE)

# Diversity analyses ####

shannon = diversity(ps@otu_table)

 write.csv(shannon, "./Output/Shannon_Diversity_by_Tree.csv", quote = FALSE)

div_aov = aov(shannon ~ ps@sam_data$Taxon * ps@sam_data$ALTITUDE)
 sink(file = "./Output/Stat-Tests.txt", append = TRUE)
 print("")
 print("Shannon Diversity ANOVA - Taxon*Altitude")
 summary(div_aov)
 sink(NULL)

ggplot(mapping = aes(x=ps@sam_data$Taxon,y=shannon, fill = ps@sam_data$Taxon)) +
  geom_boxplot() + theme_minimal() + labs(x="Taxon", fill="Taxon",y="Shannon Diversity") +
  theme(axis.text.x = element_text(angle=90))
 ggsave("./Output/Shannon_Diversity_by_Taxon.png", dpi=300)


# Model comparisons for Shannon diversity as dependent variable ####

# subset metadata to remaining sites
meta = meta[meta$IDENT %in% sample_data(ps_ra)$IDENT,]
meta$IDENT = factor(meta$IDENT)

mod = data.frame(Shannon = shannon, Altitude = meta$ALTITUDE,Taxon=meta$Taxon,Site=meta$Collection_Site,
           Elevation=meta$Elevation,Surface=meta$Abaxial_Surface)
names(mod)

# sink(file = "./Output/Stat-Tests.txt", append = TRUE)
print("Model Selection for Shannon Diversity")
m1 = aov(Shannon ~ Altitude+Taxon*Site+Elevation+Surface, data = mod)
summary(m1)
m2 = aov(Shannon ~ Surface+Taxon+Site, data = mod)
summary(m2)
print("Stepwise AIC")
stepAIC(m1,scope = ~Altitude+Taxon*Site+Elevation+Surface)
print("Best model for Diversity looks to be ~Altitude+Taxon+Site")
# sink(NULL)

# Plot of Shannon vs altitude+taxon ####
ggplot(mod, aes(x=Altitude,y=Shannon, color=Taxon)) +
  geom_point() + stat_smooth(se=FALSE, method = "lm") +
  theme_bw() + labs(x="Elevation (m)",y="Shannon Diversity")
 ggsave("Output/Diversity_vs_Altitude_w_Taxon.png", dpi=300)


 
 
 
 
 
 # Taxonomy tables and qualitative analyses ####
taxtab = sort(table(tax_table(ps)), decreasing = TRUE)
head(taxtab, 30)
