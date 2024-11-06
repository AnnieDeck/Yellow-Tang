##clear all objects from the work space and start with a clean environment 
rm(list=ls())

##load libraries 
library("ggplot2") 
library("phyloseq") 
library("tidyverse")
library("vegan")
library("pals")
library("treemap")
library("treemapify") 
library("ggfittext")
library("devtools")
library("pairwiseAdonis")
library("eulerr")
library("microbiome")
library("DivNet")
library("tidyverse")
library("dplyr")

##set wd
setwd("~/Desktop/Chapter3_data/Rstudio_highread")

##start to bring things into R!

##read in ASV(feature) table
asv_table_ch3 = read.csv("feature-table.csv", sep=",", row.names=1)
asv_table_ch3 = as.matrix(asv_table_ch3)

##read in taxonomy 
Taxonomy_table_ch3 = read.csv("taxtable_clean_v2.csv", sep=",", row.names=1)
Taxonomy_table_ch3 = as.matrix(Taxonomy_table_ch3)

                  
##read in metadata
##fill=TRUE fills character fields with empty strings instead of NAs
##header=TRUE tells R that the first row of the file contains the variable names 
Metadata_table_ch3 = read.csv("Ch3Metadata_subset.csv", header=TRUE, fill=TRUE)

##takes first column and makes it the row name
rownames(Metadata_table_ch3) <- Metadata_table_ch3$SampleID
Metadata_table_ch3$SampleID


##import as phyloseq objects
OTU_phy_object = otu_table(asv_table_ch3, taxa_are_rows = TRUE)
TAX_phy_object = tax_table(Taxonomy_table_ch3)
META_phy_object = sample_data(Metadata_table_ch3)

##read in tree
phy_tree = read_tree("tree.nwk")

##check metadata
sample_names(META_phy_object)
colnames(META_phy_object)

#make sure sample names match *exactly*
sample_names(META_phy_object)
sample_names(OTU_phy_object)

##create the phyloseq object
physeq_ch3 = phyloseq(OTU_phy_object, TAX_phy_object, phy_tree)

physeq_ch3_meta = merge_phyloseq(physeq_ch3, META_phy_object)
ntaxa(physeq_ch3_meta)
nsamples(physeq_ch3_meta)


##Create table, which shows the number of features for each tax rank
table(tax_table(physeq_ch3_meta)[, "Domain"], exclude = NULL)
table(tax_table(physeq_ch3_meta)[, "Phylum"], exclude = NULL)
table(tax_table(physeq_ch3_meta)[, "Genus"], exclude = NULL)

            ########removing ASVs unassigned at domain (5 ASVs) ##########
physeq_ch3_cleaned_domain <- subset_taxa(physeq_ch3_meta, !is.na(Domain) & !Domain %in% c("", "Unassigned")) 

table(tax_table(physeq_ch3_cleaned_domain)[, "Domain"], exclude = NULL)


            ########removing ASVs unassigned at phylum##########

##Before: 22 phyla, 1,657 ASVs
physeq_phyla <- tax_glom(physeq_ch3_cleaned_domain, taxrank= "ASV_ID.1")
ntaxa(physeq_phyla)

##Features (ASVs) with ambiguous phylum annotation are removed 
physeq_ch3_cleaned <- subset_taxa(physeq_ch3_cleaned_domain, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#After: still 22 phlya, filtered to 1,447 ASVs (210 ASVs removed)
physeq_phyla_cleaned <- tax_glom(physeq_ch3_cleaned, taxrank= "ASV_ID.1")
ntaxa(physeq_phyla_cleaned)

## Create table, which shows the number of features for domain
#16 archaea, 1431 bacteria
table(tax_table(physeq_ch3_cleaned)[, "Domain"], exclude = NULL)

##Read depth
mean(phyloseq::sample_sums(physeq_ch3_cleaned))
sum(phyloseq::sample_sums(physeq_ch3_cleaned))

##look at tax table to make sure all is well- looks good
#write.csv(tax_table(physeq_ch3_cleaned), "taxtable_filtered.csv")

             ##########cleaning labels of taxonomy table############

##looking at unassigned ASVs at the genus level- 522 (a lot)
##taxglom throws out ASVs unassigned at the genus level. This means I might be losing some fairly abundant ASVs just because we don't know exactly what they are.
Unassigned=subset_taxa(rel_abundances, Genus=="")
ntaxa(Unassigned)


test=tax_glom(rel_abundances, taxrank="Genus", NArm=FALSE)
ntaxa(test)

v2 <-test %>% psmelt() %>% dplyr::group_by(SampleID, Timepoint, SampleType, Genus) %>% dplyr::summarise(sum = sum(Abundance))
write.csv(v2, "v2.csv")

##step one re-naming all ASVs unassigned at the genus level

##pulling out tax table
Taxonomy_table_ch3 <- data.frame(tax_table(physeq_ch3_cleaned))

##change columns to characters
for (i in 1:8){ Taxonomy_table_ch3[,i] <- as.character(Taxonomy_table_ch3[,i])}

Taxonomy_table_ch3[is.na(Taxonomy_table_ch3)] <- ""


##re-naming 
for (i in 1:nrow(Taxonomy_table_ch3)){
  if (Taxonomy_table_ch3[i,2] == ""){
    kingdom <- paste("Kingdom_", Taxonomy_table_ch3[i,1], sep = "")
    Taxonomy_table_ch3[i, 2:7] <- kingdom
  } else if (Taxonomy_table_ch3[i,3] == ""){
    phylum <- paste("Phylum_", Taxonomy_table_ch3[i,2], sep = "")
    Taxonomy_table_ch3[i, 3:7] <- phylum
  } else if (Taxonomy_table_ch3[i,4] == ""){
    class <- paste("Class_", Taxonomy_table_ch3[i,3], sep = "")
    Taxonomy_table_ch3[i, 4:7] <- class
  } else if (Taxonomy_table_ch3[i,5] == ""){
    order <- paste("Order_", Taxonomy_table_ch3[i,4], sep = "")
    Taxonomy_table_ch3[i, 5:7] <- order
  } else if (Taxonomy_table_ch3[i,6] == ""){
    family <- paste("Family_", Taxonomy_table_ch3[i,5], sep = "")
    Taxonomy_table_ch3[i, 6:7] <- family
  } else if (Taxonomy_table_ch3[i,7] == ""){
    Taxonomy_table_ch3$Species[i] <- paste("Genus",Taxonomy_table_ch3$Genus[i], sep = "_")
  }
}


##exporting it to look at in excel- looks good
write.csv(Taxonomy_table_ch3, "taxtable_clean.csv")

##step 2 trying to get rid of "uncultured" now
for (i in 1:nrow(Taxonomy_table_ch3)){
  if (Taxonomy_table_ch3[i,2] == "p__uncultured"){
    kingdom <- paste("Kingdom_", Taxonomy_table_ch3[i,1], sep = "")
    Taxonomy_table_ch3[i, 2:7] <- kingdom
  } else if (Taxonomy_table_ch3[i,3] == "c__uncultured"){
    phylum <- paste("Phylum_", Taxonomy_table_ch3[i,2], sep = "")
    Taxonomy_table_ch3[i, 3:7] <- phylum
  } else if (Taxonomy_table_ch3[i,4] == "o__uncultured"){
    class <- paste("Class_", Taxonomy_table_ch3[i,3], sep = "")
    Taxonomy_table_ch3[i, 4:7] <- class
  } else if (Taxonomy_table_ch3[i,5] == "f__uncultured"){
    order <- paste("Order_", Taxonomy_table_ch3[i,4], sep = "")
    Taxonomy_table_ch3[i, 5:7] <- order
  } else if (Taxonomy_table_ch3[i,6] == "g__uncultured"){
    family <- paste("Family_", Taxonomy_table_ch3[i,5], sep = "")
    Taxonomy_table_ch3[i, 6:7] <- family
  } else if (Taxonomy_table_ch3[i,7] == "s__uncultured_bacterium"){
    Taxonomy_table_ch3$Species[i] <- paste("Genus",Taxonomy_table_ch3$Genus[i], sep = "_")
  }
}

##export to look at- looks good
write.csv(Taxonomy_table_ch3, "taxtable_clean_v2.csv")


              ####Physeq w/o chloroplast & mitochondria reads####

chloro= subset_taxa(physeq_ch3_meta, Order== "o__Chloroplast")
ntaxa(chloro)
#write.csv(otu_table(chloro),"Chloro_otu_table.csv" )


mito= subset_taxa(physeq_ch3_meta, Family== "f__Mitochondria")
ntaxa(mito)

##remove
NoMito <- subset_taxa(physeq_ch3_meta, !Family %in% c("f__Mitochondria"))
NoChloroNoMito <- subset_taxa(NoMito, !Order %in% c("o__Chloroplast"))

##10 taxa removed
ntaxa(NoChloroNoMito)
ntaxa(physeq_ch3_meta)

table(tax_table(NoChloroNoMito)[, "Order"], exclude = NULL)
table(tax_table(NoChloroNoMito)[, "Family"], exclude = NULL)

#rr=table(tax_table(NoChloroNoMito)[, "Phylum"], exclude = NULL)

physeq_larvae= subset_samples(NoChloroNoMito, SampleType== "Larvae")
nsamples(physeq_larvae)

##transform sample counts to relative abundance
rel_abundances = transform_sample_counts(NoChloroNoMito, function(OTU) OTU/sum(OTU))
ntaxa(rel_abundances)


##using psmelt to turn my phyloseq object into a data frame
mdf_2 = psmelt(rel_abundances)
names(mdf_2)
                         
          

##Identifying common number of ASVs between food and larvae

##timepoint 1- 102 shared ASVs
physeq_timepoint1= subset_samples(NoChloroNoMito, SampleType== "Larvae" | SampleType=="Food")
nsamples(physeq_timepoint1)

physeq_timepoint1_x= subset_samples(physeq_timepoint1, Timepoint2== "One")
nsamples(physeq_timepoint1_x)

phylo2 = filter_taxa(physeq_timepoint1_x, function(x) sum(x >= 1) == (2), TRUE)
otu_table(phylo2)

#timepoint 2- 125
physeq_timepoint2= subset_samples(physeq_timepoint1, Timepoint2== "Two")
nsamples(physeq_timepoint8)

phylo2 = filter_taxa(physeq_timepoint2, function(x) sum(x >= 1) == (2), TRUE)
otu_table(phylo2)

#timepoint 3- 167
physeq_timepoint3= subset_samples(physeq_timepoint1, Timepoint2== "Three")
nsamples(physeq_timepoint8)

phylo2 = filter_taxa(physeq_timepoint3, function(x) sum(x >= 1) == (2), TRUE)
otu_table(phylo2)

#timepoint 4- 121
physeq_timepoint4= subset_samples(physeq_timepoint1, Timepoint2== "Four")
nsamples(physeq_timepoint8)

phylo2 = filter_taxa(physeq_timepoint4, function(x) sum(x >= 1) == (2), TRUE)
otu_table(phylo2)

#timepoint 5- 158
physeq_timepoint5= subset_samples(physeq_timepoint1, Timepoint2== "Five")
nsamples(physeq_timepoint8)

phylo2 = filter_taxa(physeq_timepoint5, function(x) sum(x >= 1) == (2), TRUE)
otu_table(phylo2)

#timepoint 6- 37
physeq_timepoint6= subset_samples(physeq_timepoint1, Timepoint2== "Six")
nsamples(physeq_timepoint8)

phylo2 = filter_taxa(physeq_timepoint6, function(x) sum(x >= 1) == (2), TRUE)
otu_table(phylo2)


#timepoint 7- 36
physeq_timepoint7= subset_samples(physeq_timepoint1, Timepoint2== "Seven")
nsamples(physeq_timepoint8)

phylo2 = filter_taxa(physeq_timepoint7, function(x) sum(x >= 1) == (2), TRUE)
otu_table(phylo2)

##timepoint 8- 47 shared ASVs
physeq_timepoint8= subset_samples(physeq_timepoint1, Timepoint2== "Eight")
nsamples(physeq_timepoint8)

phylo2 = filter_taxa(physeq_timepoint8, function(x) sum(x >= 1) == (2), TRUE)
otu_table(phylo2)


                           #####making figures#####

##tree map data prep
test=tax_glom(rel_abundances, taxrank="Order")

v2 <-test %>% psmelt() %>% dplyr::group_by(SampleID, Timepoint, SampleType, Order) %>% dplyr::summarise(sum = sum(Abundance))

V3<-v2%>% dplyr::group_by(Timepoint, SampleType, Order) %>% dplyr::summarise(mean = mean(sum))

##genera that are greater than 5% of relative abundance 
V4= V3 %>% dplyr::filter(mean >0.05)



##re-naming and creating "other" category" (based on 3% threshold)
V5<-V4 %>% mutate(Order = case_when(Order== "Class_c__Bacteroidia" ~ "Unidentified Bacteroidia (Class)",
                                    Order== "Class_c__Gammaproteobacteria" ~ "Unidentified Gammaproteobacteria (Class)",
                                    Order=="o__Alteromonadales" ~ "Alteromonadales", 
                                    Order== "o__Bacillales" ~ "Bacillales",
                                    Order== "o__Bacteroidales" ~ "Bacteroidales",
                                    Order== "o__Burkholderiales" ~ "Burkholderiales",
                                    Order== "c__Cytophagales" ~ "Cytophagales", 
                                    Order == "o__Flavobacteriales" ~ "Flavobacteriales",
                                    Order=="o__Gammaproteobacteria_Incertae_Sedis" ~"Gammaproteobacteria Incertae Sedis",
                                    Order== "o__Oceanospirillales" ~ "Oceanospirillales", 
                                    Order == "o__Fusobacteriales" ~ "Fusobacteriales",
                                    Order == "o__Pseudomonadales" ~ "Pseudomonadales",
                                    Order == "o__Rhodobacterales" ~ "Rhodobacterales", 
                                    Order=="o__Sphingobacteriales" ~ "Sphingobacteriales", 
                                    Order== "Phylum_p__Firmicutes" ~ "Unidentified Firmicutes (Phylum)",
                                    Order == "Phylum_p__Proteobacteria" ~ "Unidentified Proteobacteria (Phylum)",
                                    Order == "o__Vibrionales"~ "Vibrionales", 
                                    TRUE~"zOther"))


##plotting tree map- genera that make up greater than 5% of the relative abundance

##set colors
all_taxa_color = c("#fa9fb5", "#c51b8a", "#fdbb84", "#ef3b2c", "#ffffcc", "#feb24c", "#cc4c02", "#67000d", "#c7e9b4", "palegreen2","yellow4", "#238b45", "#c6dbef", "#99d8c9","#02818a", "#253494", "#1d91c0","#54278f","#9e9ac8", "#bdbdbd")



ggplot(V5, aes(area= mean, fill= Order, subgroup=Order)) +
  geom_treemap(color = "white") +
  geom_treemap_subgroup_border(colour = "white", size = 2) +
  facet_grid(SampleType ~ Timepoint) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 12),
        strip.background = element_blank(), strip.text = element_text(color = "black", face = "bold", vjust = 0.5, hjust = 0, size = 12),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 12),
        panel.border = element_blank(),
        axis.title = element_text(face = "bold", size = 12)) +
  scale_y_continuous(expand = c(0,0)) + labs(x = "Sequence proportion by Genus", y = "") + scale_fill_manual(values = all_taxa_color)

all_taxa_color_protalve = c("#fa9fb5", "#c51b8a","#f0f0f0", "#67000d", "#ef3b2c", "#ffffcc", "#feb24c","#cc4c02", "#c7e9b4", "#238b45", "#c6dbef", "#99d8c9","#253494","#1d91c0", "#9e9ac8", "#6A3D9A", "#54278f", "#bdbdbd", "#252525", "#f0f0f0")


##facet by time point- phylum
fig1= ggplot(mdf_2, aes(x=SampleID, y= Abundance, fill= Phylum)) + facet_grid(~Timepoint, scales="free_x") + geom_col() + theme(axis.text.x = element_text(angle = 90)) + theme(legend.position="none")+ scale_fill_manual(values=as.vector(alphabet(26)))
fig1

fig1_better= ggplot(mdf_2, aes(x=SampleID, y= Abundance, fill= Genus)) + facet_grid(~Timepoint, scales="free_x") + geom_col() +theme(axis.text.x = element_blank(), axis.title.y=element_text(size= 20), axis.text.y = element_text(size= 18), strip.text = element_text(size=15, face= "bold"), legend.title=element_text(size=17), strip.background = element_rect(colour="black", fill="white"), panel.background = element_rect(fill = "white", colour = "black")) + ylab("Relative Abundance (%)") + xlab("") + scale_y_continuous(labels = scales::percent)
fig1_better

##facet by sample type- phlyum 
fig2= ggplot(mdf_2, aes(x=SampleID, y= Abundance, fill= Phylum)) + facet_grid(~SampleType, scales="free_x") + geom_col() + theme(axis.text.x = element_text(angle = 90))+ scale_fill_manual(values=as.vector(alphabet(26)))
fig2

##facet by sample type- class 
fig3= ggplot(mdf_2, aes(x=SampleID, y= Abundance, fill= Class)) + facet_grid(~SampleType, scales="free_x") + geom_col() + theme(axis.text.x = element_text(angle = 90))
fig3

##Subset for top 20 genus
physeq_genus <- tax_glom(physeq_nowater, taxrank= "Genus")
ntaxa(physeq_genus)
N <- 20
top20_genus <- names(sort(taxa_sums(physeq_genus), decreasing = TRUE))[1:N]
top20_genus

# Genus relative abundance
genus_rel_abundance <- transform_sample_counts(physeq_genus, function(x) x / sum(x) )

# Subset object to top N taxa
GP.genus.prop.top <- prune_taxa(top20_genus, genus_rel_abundance)

genus_df = psmelt(GP.genus.prop.top)



##facet by sample type- top 20 genera
fig4= ggplot(genus_df, aes(x=SampleID, y= Abundance, fill= Genus)) + facet_grid(~SampleType, scales="free_x") + geom_col() + theme(axis.text.x = element_text(angle = 90))+ scale_fill_manual(values=as.vector(polychrome(36)))
fig4


fig1_better_genus= ggplot(genus_df, aes(x=SampleID, y= Abundance, fill= Genus)) + facet_grid(~Timepoint, scales="free_x") + geom_col() +theme(axis.text.x = element_blank(), axis.title.y=element_text(size= 20), axis.text.y = element_text(size= 18), strip.text = element_text(size=15, face= "bold"), legend.title=element_text(size=17), strip.background = element_rect(colour="black", fill="white"), panel.background = element_rect(fill = "white", colour = "black")) + ylab("Relative Abundance (%)") + xlab("") + scale_fill_manual(values = all_taxa_color)+ scale_y_continuous(labels = scales::percent)
fig1_better_genus


physeq_nowater= subset_samples(NoChloroNoMito, SampleType== "Larvae" | SampleType=="Food")
nsamples(physeq_nowater)


fig1_better_genus= ggplot(genus_df, aes(x=SampleID, y= Abundance, fill= Genus)) + facet_grid(~Timepoint, scales="free_x") + geom_col() +theme(axis.text.x = element_blank(), axis.title.y=element_text(size= 20), axis.text.y = element_text(size= 18), strip.text = element_text(size=15, face= "bold"), legend.title=element_text(size=17), strip.background = element_rect(colour="black", fill="white"), panel.background = element_rect(fill = "white", colour = "black")) + ylab("Relative Abundance (%)") + xlab("") + scale_fill_manual(values=as.vector(alphabet(26)))+ scale_y_continuous(labels = scales::percent)
fig1_better_genus





                           ######beta diversity########
##NMDS bray-curtis
##Not using this one-just for comparison to Aitchison
GP.ord <- ordinate(NoChloroNoMito, "NMDS", "bray")
GP.ord

NMDS_bray = plot_ordination(NoChloroNoMito, GP.ord, type="samples", color="Timepoint2", shape="SampleType") + annotate("text", x=-1.2, y=1.1, label= "Stress= 0.151") + theme_bw() + geom_point(size=2.5) +
  stat_ellipse(aes(group = Timepoint), linetype = 2)


NMDS_bray

NMDS_taxa_color = c("#fa9fb5", "#c51b8a", "#ef3b2c", "#feb24c", "yellow4","#02818a", "#1d91c0", "#253494")

NMDS_taxa_color2 = c("#253494", "yellow4", "#feb24c","#fa9fb5", "#1d91c0","#02818a","#ef3b2c","#c51b8a")

df2$Timepoint2 <- factor(df2$Timepoint2, levels = c("One","Two","Three", "Four", "Five", "Six", "Seven", "Eight"))

##NMDS with CLR and euclidean data transformation
CLR_All_Phyto=microbiome::transform(NoChloroNoMito, transform="clr")
metadata <- as(sample_data(CLR_All_Phyto), "data.frame")
dist.uf <- phyloseq::distance(CLR_All_Phyto, method = "euclidean")

GP.ord.euclidean <- ordinate(CLR_All_Phyto, "NMDS", "euclidean")
GP.ord.euclidean


NMDS_allsamples = plot_ordination(CLR_All_Phyto, GP.ord.euclidean, type="samples", color="Timepoint2", shape="SampleType") + annotate("text", x=-78, y=83.5, label= "Stress= 0.151",fontface =2, size= 5.5) + guides(color = guide_legend(title = "Days Post Hatch"))+ geom_point(size=3.5) + theme_bw()+ theme(text = element_text(size = 16)) + stat_ellipse(aes(group = Timepoint2), linetype = 2) + scale_color_manual(values = NMDS_taxa_color2) 
NMDS_allsamples

class(metadata$Timepoint)
metadata$Timepoint2 <- factor(metadata$Timepoint2, levels = c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight"))

levels(metadata$Timepoint2)




###NMDS with only larval samples
CLR_All_Phyto=microbiome::transform(physeq_larvae, transform="clr")
metadata <- as(sample_data(CLR_All_Phyto), "data.frame")
dist.uf <- phyloseq::distance(CLR_All_Phyto, method = "euclidean")

GP.ord.euclidean <- ordinate(CLR_All_Phyto, "NMDS", "euclidean")
GP.ord.euclidean


NMDS_euclidean = plot_ordination(CLR_All_Phyto, GP.ord.euclidean, type="samples", color="Timepoint2") + annotate("text", x=-32, y=70.5, label= "Stress= 0.105", fontface =2) + geom_point(size=2.8) + guides(color = guide_legend(title = "Timepoint"))+ theme_bw() + scale_color_manual(values = NMDS_taxa_color)  
NMDS_euclidean


library(ggforce)
NMDS_larvae = plot_ordination(CLR_All_Phyto, GP.ord.euclidean, type="samples", color="Timepoint2") + annotate("text", x=-68, y=70.5, label= "Stress= 0.105") + geom_point(size=2.8) + theme_bw() + scale_color_manual(values = NMDS_taxa_color) 
NMDS_larvae

                       #######Alpha diversity########
##not working with this data set for some reason
##Error in pick_base(W) : Yikes! No taxa observed in all samples! Pick which taxon is to be the base

sample_data(physeq_ch3_cleaned)$Timepoint
sample_data(physeq_ch3_cleaned)$SampleType


divnet_pointestimate <- physeq_ch3_cleaned %>%  tax_glom("Phylum")%>% 
  divnet(X = "SampleType", ncores = 1)
divnet_pointestimate
plot(divnet_pointestimate)
testDiversity(divnet_pointestimate, "shannon")


dv_ch3<- divnet_pointestimate


ests <- sapply(dv_ch3[["shannon"]], function(x) x$estimate)
lci <- sapply(dv_ch3[["shannon"]], function(x) x$interval[1])
uci <- sapply(dv_ch3[["shannon"]], function(x) x$interval[2])

df_ch3 <- data.frame("SampleID" = names(ests),
                 "h0" = ests, lci, uci,
                 dv_ch3$X) 
df_ch3




df2_ch3=left_join(df_ch3, Metadata_table_ch3, by="SampleID")

plot(df2_ch3)


ff=ggplot(df2_ch3, aes(x = SampleID, xend = SampleID, color=SampleType)) +
  geom_point(aes(x = SampleID, y = ests)) +
  geom_segment(aes(y = lci, yend = uci)) +
  ylab(paste("Shannon estimate")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ff





##SampleType_Genus- not working yet

divnet_pointestimate <- physeq_ch3_cleaned %>%  tax_glom("Genus")%>% 
  divnet(X = "SampleType", ncores = 1)




##I want to look at species richness with breakaway
ba=breakaway(physeq_ch3_cleaned)
ba
plot(ba, physeq_ch3_cleaned, color = "SampleType")

By_Sample_Type<- betta(summary(ba)$estimate,
                               summary(ba)$error,
                               make_design_matrix(physeq_ch3_cleaned, "SampleType"))
By_Sample_Type_bet=as.data.frame(By_Sample_Type$table)
By_Sample_Type_bet
write.csv(By_Sample_Type_bet, "By_Sample_Type_Tucker.csv")



                      ######summary statistics######

##subsetting dataframe to only retain larval samples (3 replicates/timepoint)
physeq_larvae= subset_samples(NoChloroNoMito, SampleType== "Larvae" | SampleType=="Larvae")

subset8= subset_samples(NoChloroNoMito, Timepoint=="8")

##Adonis- timepoint significant 
CLR_All_Phyto=microbiome::transform(physeq_larvae, transform="clr")
metadata <- as(sample_data(CLR_All_Phyto), "data.frame")
dist.uf <- phyloseq::distance(CLR_All_Phyto, method = "euclidean")
ps.adonis <- adonis2(dist.uf ~ Timepoint, data = metadata, perm=999, by="margin")
ps.adonis



##Pairwise comparison between time points- seems to not be working when I run with only larvae samples
CLR_All_Phyto=microbiome::transform(subset8, transform="clr")
metadata <- as(sample_data(CLR_All_Phyto), "data.frame")
dist.uf <- phyloseq::distance(CLR_All_Phyto, method = "euclidean")
ps.adonis <- pairwise.adonis2(dist.uf ~ SampleType, data = metadata, perm=999, by="margin")
ps.adonis

###test SampleType
CLR_All_Phyto=microbiome::transform(NoChloroNoMito, transform="clr")
metadata <- as(sample_data(CLR_All_Phyto), "data.frame")
dist.uf <- phyloseq::distance(CLR_All_Phyto, method = "euclidean")
ps.adonis <- adonis2(dist.uf ~ SampleType, data = metadata, perm=999, by="margin")
ps.adonis


###SampleType was not statistically significant but it just looks at all groups against each other, can use pairwise.adonis2 to make pairwise comparisons and see that water vs larvae is significant 
CLR_All_Phyto=microbiome::transform(NoChloroNoMito, transform="clr")
metadata <- as(sample_data(CLR_All_Phyto), "data.frame")
dist.uf <- phyloseq::distance(CLR_All_Phyto, method = "euclidean")
ps.adonis <- pairwise.adonis2(dist.uf ~ SampleType, data = metadata, perm=999, by="margin")
ps.adonis



###test interaction between timepoint and sampletype
CLR_All_Phyto=microbiome::transform(physeq_ch3_cleaned, transform="clr")
metadata <- as(sample_data(CLR_All_Phyto), "data.frame")
dist.uf <- phyloseq::distance(CLR_All_Phyto, method = "euclidean")
ps.adonis <- adonis2(dist.uf ~ Timepoint+SampleType, data = metadata, perm=999, by="margin")
ps.adonis


##Sarah code 

mod <- betadisper(dist.uf , metadata$Year)
mod
anova(mod)
boxplot(mod)
plot(mod)

pathotype.disp.TukeyHSD <- TukeyHSD(mod)
pathotype.disp.TukeyHSD

anosim(dist.uf, metadata$Cluster_DeSeq_Spatial, permutations = 1000)


pathotype.disp.TukeyHSD <- TukeyHSD(mod)
pathotype.disp.TukeyHSD

###beta disperser
#null hypothesis is that the replicates are more similar to those within the same extraction group than those in a different extraction group 
mod <- betadisper(dist.uf, metadata$Timepoint)
mod
anova(mod)
boxplot(mod)
plot(mod)
##great no significant dispersion, you can use your results from the adonis plots

##how about for SampleType
mod <- betadisper(dist.uf, metadata$SampleType)
mod
anova(mod)
boxplot(mod)
plot(mod)
##great no significant dispersion, you can use your results from the adonis plots

##Tukey test
##So this is looking at the last betadisper run. If you had significance you could look at which groups are significant 
pathotype.disp.TukeyHSD <- TukeyHSD(mod)
pathotype.disp.TukeyHSD




                   ########summary stats with dplr##########  

##Relative abundance summaries by time point and sample type
##First step is to tax glom
test=tax_glom(rel_abundances, taxrank="Phylum")

v2 <-test %>% psmelt() %>% dplyr::group_by(SampleID, Timepoint, SampleType,Phylum) %>% dplyr::summarise(sum = sum(Abundance))

V3<-v2%>% dplyr::group_by(SampleType, Phylum) %>% dplyr::summarise(mean = mean(sum)*100, sd=sd(sum)*100)

##Relative abundance summaries by sample type
v2_overall <-rel_abundances %>% psmelt() %>% dplyr::group_by(SampleID, SampleType, Phylum) %>% dplyr::summarise(sum = sum(Abundance))

V3_overall<-v2%>% dplyr::group_by(SampleType, Phylum) %>% dplyr::summarise(mean = mean(sum)*100, sd=sd(sum)*100)


##exporting OTU table from rel_abundances into a CSV
#write.csv(as.data.frame(otu_table(rel_abundances)), "otu_rel2All_phyto.csv")



##chloroplast 
chloro=subset_taxa(physeq_ch3_cleaned, Order=="o__Chloroplast")
tax_table(chloro)
plot_bar(chloro, fill= "ASV_ID.1")


nochloro=subset_taxa(rel_abundances, Order!="o__Chloroplast")
tax_table(nochloro)

                 #####exploring ASVs within genera#######

##looking at vibrio
vibrio=subset_taxa(NoChloroNoMito, Genus== "g__Vibrio")
tax_table(vibrio)


##using psmelt to turn my phyloseq object into a data frame
mdf_2_vibrio = psmelt(vibrio)

vibrio_ASVs<- ggplot(mdf_2_vibrio, aes(x=SampleID, y= Abundance, fill= ASV_ID.1)) + facet_grid(~SampleType, scales="free_x") + geom_col() + theme(axis.text.x = element_text(angle = 90))+ scale_fill_manual(values=as.vector(alphabet(26))) +theme_bw() 
vibrio_ASVs

plot_tree(vibrio)

#looking at uncultured 
uncultured=subset_taxa(NoChloroNoMito, Genus== "g__uncultured")
ntaxa(uncultured)

plot_bar(uncultured, fill="ASV_ID.1", title="Uncultured ASVs")+facet_wrap(~SampleType)+theme(legend.position = "none")

##looking at proteobacteria
proteo=subset_taxa(NoChloroNoMito, Phylum== "p__Proteobacteria")
plot_bar(proteo, fill="Genus", title="Genus of Proteo")+facet_wrap(~SampleType)
psmelt(proteo)

# Calculate relative abundances
proteo_rel_abund <- transform_sample_counts(proteo, function(x) x / sum(x))

# Calculate total abundance of Proteobacteria
total_proteo_abund <- sum(taxa_sums(proteo_rel_abund))

# Subset to genus Vibrio
vibrio <- subset_taxa(proteo_rel_abund, Genus == "g__Vibrio")

# Calculate total abundance of Vibrio
total_vibrio_abund <- sum(taxa_sums(vibrio))

# Calculate percentage of Vibrio within Proteobacteria
vibrio_perc_within_proteo <- (total_vibrio_abund / total_proteo_abund) * 100

# Print the percentage
print(vibrio_perc_within_proteo)

                    #####DeSeq- not using this right now######
#Figure 4b
########
library(DESeq2)
glom=tax_glom(physeq_final, "Genus")
ntaxa(glom)
dex = phyloseq::phyloseq_to_deseq2(glom~ SampleType)
tax_table(glom)



diagdds_pol = DESeq2::DESeq(dex, test="Wald", fitType="parametric")
res = results(diagdds_pol)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab
test=as.data.frame(sigtab)
test$ASV<- rownames(test)
write.csv(test, "test.csv")

taxa_for_deseq=as.data.frame(tax_table(glom))
taxa_for_deseq$ASV<- rownames(taxa_for_deseq)
rownames(taxa_for_deseq) <- NULL
print(taxa_for_deseq)
write.csv(taxa_for_deseq, "taxa_for_deseq_glom.csv")


Deseq_taxa=left_join(test, taxa_for_deseq, by="ASV")
write.csv(Deseq_taxa, "Deseq_taxa.csv")



test$DeSeq_Spatial_Cluster=sigtab$Domain
test$DeSeq_Spatial_Cluster="DeSeq_Spatial_Variable"
sigtab2=as.data.frame(sigtab2)
paper=left_join(taxa_for_deseq,sigtab2,by="ASV")

write.csv(paper,"taxonmy_all_Deseqdata.csv")
left_join()


res <- results(diagdds_pol, contrast=c("SampleType", "Larvae","Food"))
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]

sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_final_final)[rownames(sigtab), ], "matrix"))
head(sigtab)
write.csv(sigtab, "Genus_DESeq_All_Phyto_Clustering_NEWCLR_Transition_Offshore.csv")


sigtab2 = cbind(as(sigtab, "data.frame"), as(as.data.frame(tax_table(physeq_final_final))[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab


                    ####Analyzing core community####


##Seeing how many ASVs are shared across all three sample types at each time point
table(meta(NoChloroNoMito)$Timepoint, useNA = "always")
          
pseq.rel <- microbiome::transform(NoChloroNoMito, "compositional")

Timepoint <- unique(as.character(meta(pseq.rel)$Timepoint))
print(Timepoint)

list_core <- c()


for (n in Timepoint){ # for each variable n in Sample Type
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, Timepoint == n) # Choose sample from Sample Type by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in at least 90% samples 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

##these are the core ASVs
print(list_core)
test_bubble_1<-as.data.frame(list_core[1])
test_bubble_2<-as.data.frame(list_core[2])
test_bubble_3<-as.data.frame(list_core[3])
test_bubble_4<-as.data.frame(list_core[4])
test_bubble_5<-as.data.frame(list_core[5])
test_bubble_6<-as.data.frame(list_core[6])
test_bubble_7<-as.data.frame(list_core[7])
test_bubble_8<-as.data.frame(list_core[8])

##combining above lists
##error "differing number of rows" 
?cbind
cbind(test_bubble_1, test_bubble_2)

##plotting ven diagram 
mycols <- c(nonCRC="#d6e2e9", CRC="#cbf3f0", H="#fcf5c7") 
plot(venn(list_core),
     fills = mycols)

?`eulerr-package`

proteo=subset_taxa(rel_abundances, Phylum== "p__Proteobacteria")


### lets look at the prevalence of ASVs over samples 
physeq_final_final_binary <- transform_sample_counts(physeq_ch3_cleaned, function(abund) 1*(abund>0))
binary_melt=psmelt(physeq_final_final_binary)
names(binary_melt)

##how many samples each ASV occurs in- lots of samples only occurring in 2 samples AKA not a ton of ASVs found across many of our samples.
tax_binary_table =as.data.frame((tax_table(physeq_final_final_binary)))
phy_binary_table =as.data.frame((otu_table(physeq_final_final_binary)))
sums=rowSums(phy_binary_table)
sums
hist(sums)

##Subsetting for each sample type- I believe this is for the ven diagram 
Food=subset_samples(physeq_final_final_binary, SampleType=="Food")
Food_binary_table =as.data.frame((otu_table(Food)))
Food_binary_table$sum=rowSums(Food_binary_table)
hist(Food_binary_table$sum)
Food_binary_table$Food[Food_binary_table$sum > 0] = 1
Food_binary_table$Food[is.na(Food_binary_table$Food)] <- 0
Food_binary_table$Food
Food_binary_table2=Food_binary_table%>%select(Food)


Larvae=subset_samples(physeq_final_final_binary, SampleType=="Larvae")
Larvae_binary_table =as.data.frame((otu_table(Larvae)))
Larvae_binary_table$sum=rowSums(Larvae_binary_table)
Larvae_binary_table$Larvae[Larvae_binary_table$sum > 0] = 1
Larvae_binary_table$Larvae[is.na(Larvae_binary_table$Larvae)] <- 0
Larvae_binary_table$Larvae
Larvae_binary_table2=Larvae_binary_table%>%select(Larvae)




Water=subset_samples(physeq_final_final_binary, SampleType=="Water")
Water_binary_table =as.data.frame((otu_table(Water)))
Water_binary_table$sum=rowSums(Water_binary_table)
Water_binary_table$Water[Water_binary_table$sum > 0] = 1
Water_binary_table$Water[is.na(Water_binary_table$Water)] <- 0
Water_binary_table$Water
Water_binary_table2=Water_binary_table%>%select(Water)



Water_binary_table3=filter(Water_binary_table2,Water==1)
write.csv(Water_binary_table3, "Water_binary_table3.csv")

Food_binary_table3=filter(Food_binary_table2,Food==1)
write.csv(Food_binary_table3, "Food__binary_table3.csv")


Larvae_binary_table3=filter(Larvae_binary_table2,Larvae==1)
write.csv(Larvae_binary_table3, "Larvae_binary_table3.csv")



water_names=as.vector(rownames(Water_binary_table3))
larvae_names=as.vector(rownames(Larvae_binary_table3))
food_names=as.vector(rownames(Food_binary_table3))


install.packages("VennDiagram")
library(VennDiagram)
venn.diagram(
  x = list(water_names, larvae_names, food_names),
  category.names = c("water" , "larvae" , "food"),
  filename = 'venn_diagramm.png',
  output=TRUE
)


