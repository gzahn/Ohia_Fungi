# deleted from 02_Ohia_Analyses.R
# Starting at line 64

# count number of taxa that are host-specific in each sample

otu.data <- as.data.frame(as.matrix(otu_table(ps)))
otu.data[otu.data>0] <- 1

# build data frame showing taxa that are only found in one tree species
filling <- data.frame(`M. prostrata`=1:ntaxa(ps),`M.p. var. fakerug`=1:ntaxa(ps),
                      `M.p. var. glaberrima`=1:ntaxa(ps), `M.p. var. icana`=1:ntaxa(ps),
                      `M.p. var. microglab`=1:ntaxa(ps), `M.p. var. orb`=1:ntaxa(ps),
                      `M. rugosa`=1:ntaxa(ps),`M. tremuloides`=1:ntaxa(ps))

j=1
for(i in levels(sample_data(ps)$Taxon)){
  taxon.otu.table=otu.data[sample_data(ps)$Taxon==i,]
  taxon.otu.sums=colSums(taxon.otu.table)
  taxon.otu.sums[taxon.otu.sums>0] <- 1
  filling[,j] <- taxon.otu.sums
  j=j+1
}

row.names(filling) <- taxa_names(ps)
filling = filling[rowSums(filling)==1,]
specific = subset_taxa(ps,taxa_names(ps)%in%row.names(filling))
specific <- merge_samples(ps,"Taxon")
sample_names(specific)
filling=as.data.frame(t(filling))
row.names(filling) <- sample_names(specific)


specific <- phyloseq(otu_table(filling,taxa_are_rows=FALSE),tax_table(specific),sample_data(specific))
