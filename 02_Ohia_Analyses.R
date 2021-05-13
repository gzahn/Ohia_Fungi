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
library(deseq2)
library(vegan)
library(dplyr)
library(reshape2)
library(ade4)
library(MASS)
library(ggbiplot)
library(fitdistrplus)
library(FSA)
library(indicspecies)


# Note: the plyr package must be installed, but should not be loaded

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

# relative abundance
ps_ra = transform_sample_counts(ps, function(x) x/sum(x))

# Host specificity ####
host.specif = (specnumber(otu_table(ps),sample_data(ps)$Taxon,2))
png("./Output/single-host_taxa_abundance.png")
plot((taxa_sums(ps)[host.specif == 1]),main = "Abundance of taxa found in only a single host", ylab = "Absolute Abundance")
dev.off()

abundant.1host.taxa = which(taxa_sums(ps)[host.specif == 1] > 1000)
host.specific.taxa = as.data.frame(tax_table(ps)[abundant.1host.taxa]) # names of taxa found on only one host with more than 1000 abs. abundance
write.csv(host.specific.taxa, "./Output/host-sepcific_taxa.csv",quote = FALSE,row.names = FALSE)


# indicspecies analysis ... are results similar?

asv <- ps@otu_table %>% as.data.frame()
taxa <- ps@sam_data$Taxon %>% unique %>% as.character()

for(i in taxa){
  indic <- indicators(asv,cluster=ps@sam_data$Taxon,group=i, enableFixed = TRUE,max.order = 1)
  indic <- indic$group.vec
  assign(paste0("indic.",i),value = indic,envir = .GlobalEnv)
}



hstaxa = specnumber(otu_table(ps), sample_data(ps)$Taxon,2)
hstaxa = which(hstaxa == 1)
hsotus = otu_table(ps)[,hstaxa]
hsotus[hsotus > 0] <- 1

# Plot host specificity against ALTITUDE
# add column to meta that shows number of host-specific taxa found in each sample

meta2 = meta
meta2 = (meta2[meta2$IDENT %in% names(hsotus@.Data[,2]),])
meta2 = arrange(meta2,IDENT)


hs.count = rowSums(hsotus)
meta2$No.HS_Taxa = hs.count
meta2$Proportion.HS = meta2$No.HS_Taxa / vegan::specnumber(otu_table(ps))
meta2$Proportion.HS.total = meta2$Proportion.HS / rowSums(otu_table(ps))


ggplot(meta2, aes(x=ALTITUDE,y=Proportion.HS)) +
  geom_point(aes(color=Taxon)) + stat_smooth( se=FALSE) +
  theme_bw() + labs(x="Elevation (m)",y="Proportion of host-specific taxa") +
  theme(legend.text = element_text(face="italic"))
ggsave("./Output/Host-Specificity_elevation.png", dpi=300, height = 10, width = 12)

# one for each collection site
sites <- c("Aiea Ridge","Konahuanui","Kuliouou","Mt. Kaala")
meta2$Collection_Site <- str_replace(meta2$Collection_Site,pattern = "`",replacement = "")

for(i in sites){
p=ggplot(meta2[meta2$Collection_Site == i,], aes(x=ALTITUDE,y=Proportion.HS)) +
  geom_point(aes(color=Taxon)) + stat_smooth( se=FALSE) +
  theme_bw() + labs(x="Elevation (m)",y="Proportion of host-specific taxa") +
  theme(legend.text = element_text(face="italic")) + ggtitle(i)
assign(i,p, envir = .GlobalEnv)
}
multiplots = arrangeGrob(`Aiea Ridge`,Konahuanui,Kuliouou,`Mt. Kaala`)
ggsave(file="./Output/Host-Specificity_Elevation_By_Site.png",multiplots, dpi=300, width = 10, height = 6)

ggplot(meta2, aes(x=ALTITUDE,y=Proportion.HS,color=Taxon)) +
  geom_point(aes(color=Taxon)) + stat_smooth( se=FALSE, method = "lm") +
  theme_bw() + labs(x="Elevation (m)",y="Proportion of host-specific taxa") +
  theme(legend.text = element_text(face="italic"))
ggsave("./Output/Host-Specificity_elevation2.png", dpi=300, height = 10, width = 12)

# anova ... does tree taxon affect number of host-specific fungi found?
hsotus = hsotus[,taxa_sums(hsotus) > 1]

sample_sums(hsotus)
# meta2 = meta2[meta$IDENT %in% sample_data(ps)$IDENT,]

# get number of reps in each taxon
tree.reps = meta2 %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(N=n())
tree.reps=tree.reps[1:8,]

meta2$reps = as.numeric(as.character(plyr::mapvalues(meta2$Taxon, from = tree.reps$Taxon, to = tree.reps$N)))

# Check distribution
plotdist(sqrt(sample_sums(hsotus)), histo = TRUE, demp = TRUE)
descdist(sqrt(sample_sums(hsotus)), discrete = FALSE, boot = 1000)
summary(fitdist(sqrt(sample_sums(hsotus)), "norm"))
descdist(sqrt(meta2$Proportion.HS),discrete = FALSE,boot = 1000)


# Make ANOVA model with sq-rt transformed data 
mod = aov(sqrt(sample_sums(hsotus)) ~ sample_data(ps)$ALTITUDE * sample_data(ps)$Taxon * meta2$reps)
summary(mod) # no...almost?

as.character(sample_data(ps)$IDENT)
as.character(meta2$IDENT)

meta2 = arrange(meta2,IDENT)
mod2 <- aov(meta2$Proportion.HS ~ sample_data(ps)$ALTITUDE * sample_data(ps)$Taxon * meta2$reps)
summary(mod2)


norm.fit = fitdistr(sample_sums(hsotus), "normal")
plot(sample_sums(hsotus))



kruskal.test(sample_sums(hsotus),meta2$reps)
plot(sample_sums(hsotus) ~ sample_data(ps)$ALTITUDE)

sink("./Output/host-specificity_anova_table.txt")
noquote(print("Prevalence of host-specific fungi"))
noquote(print(""))
noquote(print(""))
noquote(print("Model includes Ohia taxon, Altitude, and number of replicate trees."))
noquote(print(""))
summary(mod)
closeAllConnections()

# Change taxa names to remove rank tags
tax_table(ps)[,"Kingdom"] <- str_remove(tax_table(ps)[,"Kingdom"],"k__")
tax_table(ps)[,"Phylum"] <- str_remove(tax_table(ps)[,"Phylum"],"p__")
tax_table(ps)[,"Order"] <- str_remove(tax_table(ps)[,"Order"],"o__")
tax_table(ps)[,"Class"] <- str_remove(tax_table(ps)[,"Class"],"c__")
tax_table(ps)[,"Family"] <- str_remove(tax_table(ps)[,"Family"],"f__")
tax_table(ps)[,"Genus"] <- str_remove(tax_table(ps)[,"Genus"],"g__")
tax_table(ps)[,"Species"] <- str_remove(tax_table(ps)[,"Species"],"s__")

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


# re-order for updated plots of tree taxa
reordered_taxa <- psm_ra@sam_data$Taxon[c(4,3,8,5,2,6,7,1)]
psm_ra@sam_data$Taxon
psm_ra@sam_data$TaxonCode <- c("A","F","G","I","L","B","R","T")

sample_data(psm_ra)[c(4,3,8,5,2,6,7,1),]


# psm_ra@sam_data$Taxon <- factor(psm_ra@sam_data$Taxon)
# psm_ra@sam_data$Taxon <- factor(psm_ra@sam_data$Taxon, levels = reordered_taxa)
# c("I","G","T","L","F","B","R","A")
# psm_ra@sam_data$Taxon

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

write.csv(data.frame(TreeID = names(tree_richness), Richness = tree_richness), "./Output/Tree_Richness.csv", row.names = FALSE)



# heatmap of fungal orders for each sample
ps_order = tax_glom(ps, "Order")
ps_order = merge_samples(ps_order, "Taxon")
main_orders = (which(colSums(ps_order@otu_table) > 1000))

df_order = as.data.frame(decostand(ps_order@otu_table,method = "total"))
orderlabs = ps_order@tax_table@.Data[,"Order"]


# orderlabs = map(strsplit(orderlabs, "o__"),2)
# orderlabs = map(strsplit(as.character(orderlabs), "_ord_"),1)
orderlabs = orderlabs[main_orders]
orderlabs = str_remove(orderlabs,"_ord_Incertae_sedis")
df_order = df_order[,main_orders]

png("./Output/Heatmap_of_Order_by_Taxon.png",width = 8008, height = 8008, res = 300)
heatmap(as.matrix(df_order), labCol = orderlabs, col = gray.colors(50), Colv = NA, margins=c(10,16),
        cexRow = 2,cexCol = 2)
dev.off()

png("./Output/Heatmap_of_Order_by_Tree.png",width = 8008, height = 8008, res = 300)
heatmap(as.matrix(df_order), labCol = orderlabs, col = gray.colors(50), Colv = NA, margins=c(10,8),
        cexRow = 2)
dev.off()


 
# heatmap of fungal genera for each sample
ps_genus = tax_glom(ps, "Genus")
df_genus = as.data.frame(otu_table(ps_genus))
genuslabs = ps_genus@tax_table@.Data[,"Genus"]
# genuslabs = map(strsplit(as.character(genuslabs), "g__"),2)

png("./Output/Heatmap_of_Genus_by_Tree.png",width = 8008, height = 8008, res = 300)
heatmap(as.matrix(df_genus), labCol = genuslabs, col = gray.colors(50), Colv = NA, margins=c(10,8),
         cexRow = 2)
dev.off()



# Bar plot of relative abundance, split by taxon ####

plot_bar(psm_ra, fill = "Phylum",x="TaxonCode") +
  geom_bar(stat = "identity") + coord_flip()  + #facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Ohia taxon",y="Relative abundance")+ theme_bw() + 
  scale_x_discrete(limits = c("A","R","B","F","L","T","G","I")) + theme(axis.text = element_text(size=12,face="bold"),
                                                                        axis.title = element_text(size=16,face="bold"),
                                                                        legend.title = element_text(size=12,face="bold"),
                                                                        legend.text = element_text(size=10))
ggsave("./Output/BarPlot_Fungal_Phylum_by_Tree.png", height = 8, width = 12, dpi=300)


plot_bar(psm_ra, fill = "Class",x="TaxonCode") +
  geom_bar(stat = "identity") + coord_flip()  + #facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Ohia taxon",y="Relative abundance")+ theme_bw() + 
  scale_x_discrete(limits = c("A","R","B","F","L","T","G","I")) + theme(axis.text = element_text(size=12,face="bold"),
                                                                        axis.title = element_text(size=16,face="bold"),
                                                                        legend.title = element_text(size=12,face="bold"),
                                                                        legend.text = element_text(size=10))
ggsave("./Output/Barplot_Fungal_Class_by_Tree.png",height = 8,width = 12,dpi=300)



sample_data(psm_ra_site)$SiteCode <- c("Aiea","KHN","WWNKOO","MK")
# Same, but merged by site
plot_bar(psm_ra_site, fill = "Phylum",x="SiteCode") +
  geom_bar(stat = "identity") + 
  scale_x_discrete(limit = c("WWNKOO","KHN","Aiea","MK")) +
  coord_flip()  + #facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Collection site",y="Relative abundance")+ theme_bw() + theme(axis.text = element_text(size=12,face="bold"),
                                                                       axis.title = element_text(size=16,face="bold"),
                                                                       legend.title = element_text(size=12,face="bold"),
                                                                       legend.text = element_text(size=10))

ggsave("./Output/BarPlot_Fungal_Phylum_by_Site.png", dpi = 300,height = 8,width = 8)



plot_bar(psm_ra_site, fill = "Class",x="SiteCode") +
  geom_bar(stat = "identity") + 
  scale_x_discrete(limit = c("WWNKOO","KHN","Aiea","MK")) +
  coord_flip()  + #facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Collection site",y="Relative abundance")+ theme_bw() + theme(axis.text = element_text(size=12,face="bold"),
                                                                       axis.title = element_text(size=16,face="bold"),
                                                                       legend.title = element_text(size=12,face="bold"),
                                                                       legend.text = element_text(size=10))

ggsave("./Output/BarPlot_Fungal_Class_by_Site.png", dpi = 300, height = 8, width = 8)


plot_bar(psm_ra, fill = "Phylum",x="Taxon") +
  geom_bar(stat = "identity") + 
  coord_flip()  + facet_wrap(~(sample_data(psm_ra)$Abaxial_Surface)) +
  labs(x="Ohia Taxon",y="Relative Abundance") + lims(y=c(0,1)) + theme_bw()
ggsave("./Output/BarPlot_Fungal_Phylum_by_Tree_and_Surface.png", height = 8, width = 12,dpi=300)

names(sample_data(ps))

# Summarise taxa and plot ####
source("./summarize_taxa_Joey711.R")
order_summary = summarize_taxa(ps,"Order","Collection_Site")
# Remove NAs
order_summary = order_summary[order_summary$Order != "NA",]
# reorder
setorder(order_summary, -meanRA)

# rename factor levels and reorder for plotting
SiteCode = plyr::mapvalues(order_summary$Collection_Site,from=levels(order_summary$Collection_Site),to=c("Aiea Ridge", "Konahuanui", "WWNKOO", "Mt. Ka`ala"))
levels(SiteCode) <- c("Mt. Ka`ala","Aiea Ridge", "Konahuanui", "WWNKOO")
order_summary$SiteCode <- SiteCode

# plot
ggplot(order_summary, aes(x= meanRA, y=Order, facet=SiteCode)) +
  geom_point() + geom_errorbarh(aes(xmin=meanRA-sdRA,xmax=meanRA+sdRA)) + labs(x="Relative abundance") +
  facet_grid(~SiteCode) +
  theme_bw() + theme(panel.spacing = unit(2, "lines"),
                     axis.title = element_text(size=16,face="bold"))

ggsave("./Output/Summarized_taxa_Order_by_Site.png", height = 8, width = 12, dpi = 300)

order_summary = summarize_taxa(ps,"Order","Taxon")
# Remove NAs
order_summary = order_summary[order_summary$Order != "NA",]
# reorder
setorder(order_summary, -meanRA)

order_summary$TaxaCode <- plyr::mapvalues(order_summary$Taxon,from=levels(order_summary$Taxon),to=c("A","F","G","I","L","B","R","T"))
levels(order_summary$TaxaCode) <- c("I","G","T","F","L","B","R","A")

# plot
ggplot(order_summary, aes(x= meanRA, y=Order, facet=TaxaCode)) +
  geom_point() + geom_errorbarh(aes(xmin=meanRA-sdRA,xmax=meanRA+sdRA)) + labs(x="Relative abundance") +
  facet_grid(~TaxaCode) +
  theme_bw() + theme(panel.spacing = unit(0.5, "lines"),
                     axis.title = element_text(size=16,face="bold"))
ggsave("./Output/Summarized_taxa_Order_by_Taxon.png", height = 8, width = 12, dpi = 300)

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
sample_data(ps_ra)$TaxaCode <- plyr::mapvalues(sample_data(ps_ra)$Taxon,from=levels(sample_data(ps_ra)$Taxon),to=c("A","F","G","I","L","B","R","T"))
levels(sample_data(ps_ra)$TaxaCode)<- c("I","G","T","F","L","B","R","A")

NMDS = ordinate(ps_ra, method = "NMDS")
DCA = ordinate(ps_ra, method = "DCA")
CCA = ordinate(ps_ra, method = "CCA")
RDA = ordinate(ps_ra, method = "RDA")
MDS = ordinate(ps_ra, method = "MDS")
PCoA = ordinate(ps_ra, method = "PCoA")

# Plot them  ####
nmdsplot=plot_ordination(ps_ra, NMDS, color = "TaxaCode",shape = "Abaxial_Surface") + theme_minimal() + labs(color="Taxon",shape="Leaf Type") + ggtitle("NMDS")
dcaplot=plot_ordination(ps_ra, DCA, color = "TaxaCode",shape = "Abaxial_Surface") + theme_minimal() + labs(color="TaxaCode",shape="Leaf Type") + ggtitle("DCA") + theme(legend.position = "none")
ccaplot=plot_ordination(ps_ra, CCA, color = "TaxaCode",shape = "Abaxial_Surface") + theme_minimal() + labs(color="TaxaCode",shape="Leaf Type") + ggtitle("CCA")+ theme(legend.position = "none")
rdaplot=plot_ordination(ps_ra, RDA, color = "TaxaCode",shape = "Abaxial_Surface") + theme_minimal() + labs(color="TaxaCode",shape="Leaf Type") + ggtitle("RDA")+ theme(legend.position = "none")
mdsplot=plot_ordination(ps_ra, MDS, color = "TaxaCode",shape = "Abaxial_Surface") + theme_minimal() + labs(color="TaxaCode",shape="Leaf Type") + ggtitle("MDS")+ theme(legend.position = "none")
pcoaplot=plot_ordination(ps_ra, PCoA, color = "TaxaCode",shape = "Abaxial_Surface") + theme_minimal() + labs(color="TaxaCode",shape="Leaf Type") + ggtitle("PCoA")+ theme(legend.position = "none")

 png(filename = "./Output/Ordination_Plots_taxon-and-leaftype_renamed-taxa.png",height = 1000,width = 1000,res = 100)
 grid.arrange(nmdsplot,dcaplot,ccaplot,rdaplot,mdsplot,pcoaplot)
 dev.off()

 sample_data(ps_ra)$Collection_Site <- plyr::mapvalues(sample_data(ps_ra)$Collection_Site, from = "Wiliwilinui", to = "Kuliouou")
 
 nmdsplot=plot_ordination(ps_ra, NMDS, color = "Collection_Site",shape = "Abaxial_Surface") + theme_minimal() + labs(color="Site",shape="Surface") + ggtitle("NMDS") + theme(legend.position = "none")
 dcaplot=plot_ordination(ps_ra, DCA, color = "Collection_Site",shape = "Abaxial_Surface") + theme_minimal() + labs(color="Site",shape="Surface") + ggtitle("DCA")
 ccaplot=plot_ordination(ps_ra, CCA, color = "Collection_Site",shape = "Abaxial_Surface") + theme_minimal() + labs(color="Site",shape="Surface") + ggtitle("CCA")+ theme(legend.position = "none")
 rdaplot=plot_ordination(ps_ra, RDA, color = "Collection_Site",shape = "Abaxial_Surface") + theme_minimal() + labs(color="Site",shape="Surface") + ggtitle("RDA")+ theme(legend.position = "none")
 mdsplot=plot_ordination(ps_ra, MDS, color = "Collection_Site",shape = "Abaxial_Surface") + theme_minimal() + labs(color="Site",shape="Surface") + ggtitle("MDS")+ theme(legend.position = "none")
 pcoaplot=plot_ordination(ps_ra, PCoA, color = "Collection_Site",shape = "Abaxial_Surface") + theme_minimal() + labs(color="Site",shape="Surface") + ggtitle("PCoA")+ theme(legend.position = "none")
 
 png(filename = "./Output/Ordination_Plots.png",height = 3000,width = 3000,res = 300)
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
closeAllConnections()

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
meta2$outliers <- row.names(sample_data(ps_ra)) %in% outliers
adonis.3 = adonis(otu_table(ps_ra) ~ sample_data(ps_ra)$outliers)

sink(file = "./Output/PermANOVA_PCA_outlier_community_differences.txt")
adonis.3
closeAllConnections()


samdat = as.data.frame(sample_data(ps_ra))
mod2 = glm(outliers ~ Abaxial_Surface + Taxon + ALTITUDE + Collection_Site, family=binomial(link='logit'), data=meta2)
mod1 = aov(outliers ~ Abaxial_Surface, data=meta2)
# plot(mod1)

sink("./Output/Outlier_Predictors.txt")
summary(mod2)
closeAllConnections()


sample_data(ps_ra)
sink(file = "./Output/Stat-Tests.txt", append = TRUE)
print("PermANOVA -- Additive")
adonis.1
print("PermANOVA -- Interactive")
adonis.2
print("PermANOVA -- Outliers")
adonis.3
closeAllConnections()

# Find outlier trees
plot(princomp$rotation)
princomp$rotation[,2] > .4


eigenvals(princomp)


# Mantel Test ####

spatial.dist = dist(cbind(ps_ra@sam_data$LONG, ps_ra@sam_data$LAT))
comm.dist =   vegdist(as.matrix(ps_ra@otu_table))

mantel.test = mantel.rtest(spatial.dist, comm.dist, nrepet = 9999)

 sink(file = "./Output/Stat-Tests.txt", append = TRUE)
 print("MANTEL TEST")
 mantel.test
 closeAllConnections()

 png("./Output/Mantel_Plot.png")
 plot(mantel.test)
 dev.off()

# Mantel Tests for each site
sites = levels(ps_ra@sam_data$Collection_Site)

for(i in sites){
ps.temp = subset_samples(ps_ra, Collection_Site == i)
  spatial.dist = dist(cbind(ps.temp@sam_data$LONG, ps.temp@sam_data$LAT))
  comm.dist =   vegdist(as.matrix(ps.temp@otu_table))
  
  mantel.test = mantel.rtest(spatial.dist, comm.dist, nrepet = 9999)
  
  png(filename = paste0("./Output/mantel_plot_",i,".png"))
  plot(mantel.test)
  dev.off()
  
  sink("./Output/mantel_tests_each_site.txt", append = TRUE)
  print(i)
  print(mantel.test)
  print("")
  sink(NULL)
  }

 
 
# Rarefaction analyses ####
# rarecurve(ps@otu_table,label = FALSE)

# Diversity analyses ####

shannon = diversity(ps_ra@otu_table)

 write.csv(shannon, "./Output/Shannon_Diversity_by_Tree.csv", quote = FALSE)

div_aov = aov(shannon ~ ps_ra@sam_data$Taxon * ps_ra@sam_data$ALTITUDE + ps_ra@sam_data$outliers + meta2$reps)
 sink(file = "./Output/Stat-Tests.txt", append = TRUE)
 print("")
 print("Shannon Diversity ANOVA - Taxon*Altitude")
 summary(div_aov)
 sink(NULL)

shannon.df = cbind(shannon,sample_data(ps_ra)) 
 
shannon.df$TaxaCode <- plyr::mapvalues(shannon.df$Taxon,from=levels(shannon.df$Taxon),to=c("A","F","G","I","L","B","R","T"))
# levels(shannon.df$TaxaCode) <- c("I","G","T","F","L","B","R","A")
 
row.names(sample_data(ps_ra)[sample_data(ps_ra)$Taxon == "M.p. var. fakerug",])
mean(shannon[row.names(sample_data(ps_ra)[sample_data(ps_ra)$Taxon == "M.p. var. fakerug",])])


ggplot(shannon.df, mapping = aes(x=TaxaCode,y=shannon, fill = TaxaCode)) +
  geom_boxplot() + theme_minimal() + labs(x="Taxon", fill="Taxon",y="Shannon diversity") +
  scale_x_discrete(limits = c("I","G","T","L","F","B","R","A")) +
  theme(axis.title = element_text(size=16,face = "bold"),
        axis.text = element_text(size=12,face="bold"),
        legend.title = element_text(size=12,face="bold"))
ggsave("./Output/Shannon_Diversity_by_Taxon.png", dpi=300)

 

# Model comparisons for Shannon diversity as dependent variable ####

# subset metadata to remaining sites
meta = meta[meta$IDENT %in% sample_data(ps_ra)$IDENT,]
meta$IDENT = factor(meta$IDENT)

mod = data.frame(Shannon = shannon, Altitude = meta$ALTITUDE,Taxon=meta$Taxon,Site=newsites,
           Elevation=meta$Elevation,Surface=meta$Abaxial_Surface)
names(mod)

sink(file = "./Output/Stat-Tests.txt", append = TRUE)
print("Model Selection for Shannon Diversity")
m1 = aov(Shannon ~ Altitude+Taxon*Site+Elevation+Surface, data = mod)
summary(m1)
m2 = aov(Shannon ~ Surface+Taxon+Site, data = mod)
summary(m2)
m3 = aov(Shannon ~ Surface*Site*Altitude, data = mod)
summary(m3)
print("Stepwise AIC")
stepAIC(m1,scope = ~Altitude+Taxon*Site+Elevation+Surface)
sink(NULL)
?stepAIC
# Plot of Shannon vs altitude+taxon ####
ggplot(mod, aes(x=Altitude,y=Shannon, color=Taxon)) +
  geom_point() + stat_smooth(se=FALSE, method = "lm") +
  theme_bw() + labs(x="Elevation (m)",y="Shannon Diversity")
ggsave("Output/Diversity_vs_Altitude_w_Taxon.png", dpi=300)


mod$Site
# One for each site 
ggplot(mod, aes(x=Altitude,y=Shannon, color=Taxon)) +
  geom_point() + stat_smooth(se=FALSE, method = "lm") +
  theme_bw() + labs(x="Elevation (m)",y="Shannon Diversity")



# Richness
names(meta2)
meta2$Taxon
meta2$Richness <- tree_richness
meta2$Collection_Site <- newsites
meta2$Collection_Site <- factor(meta2$Collection_Site)


richness.taxon <- meta2 %>% dplyr::group_by(as.character(Taxon)) %>%
  dplyr::summarise(Mea.Richness = mean(Richness), StDev = sd(Richness)/sqrt(length(Richness)))

richness.taxon$TaxaCode <- c("A","F","G","I","L","B","R","T")


ggplot(richness.taxon, aes(x=TaxaCode,y=Mea.Richness)) +
  geom_bar(stat="identity") + geom_errorbar(aes(ymin=Mea.Richness-StDev,ymax=Mea.Richness+StDev)) + theme_bw() +
  theme(axis.title = element_text(size=16,face="bold"),axis.text = element_text(size = 12,face="bold")) + 
  labs(x="Taxon",y="Endophyte richness") +
  scale_x_discrete(limits = c("I","G","T","L","F","B","R","A"))
ggsave("./Output/Endophyte_Richness_by_Taxa-relabeled.png",dpi=300)


rich.mod1 = aov(Richness ~ Taxon, data = meta2)
summary(rich.mod1)

TukeyHSD(rich.mod1)

k1 = kruskal.test(Richness ~ Taxon, data = meta2)
dunnTest(Richness ~ Taxon, data = meta2, method = "bh")

k2 = kruskal.test(Richness ~ Collection_Site, data = meta2)
dunnTest(Richness ~ Collection_Site, data = meta2, method = "bh")
k1
