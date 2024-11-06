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
setwd("~/Desktop/Chapter3_data/Ch3_highread_all")

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
Metadata_table_ch3 = read.csv("Ch3Metadata.csv", header=TRUE, fill=TRUE)

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


##Create table, which shows the number of features for each phyla 
table(tax_table(physeq_ch3_meta)[, "Phylum"], exclude = NULL)
table(tax_table(physeq_ch3_meta)[, "Genus"], exclude = NULL)

## Create table, which shows the number of features for domain
table(tax_table(physeq_ch3_meta)[, "Domain"], exclude = NULL)

##Features with ambiguous phylum annotation are removed (16 removed)
##23 phyla and 1743 ASVs
physeq_phyla <- tax_glom(physeq_ch3_meta, taxrank= "ASV_ID_2")
ntaxa(physeq_phyla)

##filtering phyla
physeq_ch3_cleaned <- subset_taxa(physeq_ch3_meta, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#still 23 phyla and now 1,523 ASVs 
physeq_phyla <- tax_glom(physeq_ch3_cleaned, taxrank= "ASV_ID_2")
ntaxa(physeq_phyla)

## Create table, which shows the number of features for domain
#16 archaea, 1507 bacteria
table(tax_table(physeq_ch3_cleaned)[, "Domain"], exclude = NULL)


##Read depth
mean(phyloseq::sample_sums(physeq_ch3_cleaned))
sum(phyloseq::sample_sums(physeq_ch3_cleaned))


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

##12 taxa removed
ntaxa(NoChloroNoMito)
ntaxa(physeq_ch3_meta)

table(tax_table(NoChloroNoMito)[, "Domain"], exclude = NULL)



table(tax_table(NoChloroNoMito)[, "Order"], exclude = NULL)
table(tax_table(NoChloroNoMito)[, "Family"], exclude = NULL)

mean(phyloseq::sample_sums(NoChloroNoMito))
sum(phyloseq::sample_sums(NoChloroNoMito))
list(phyloseq::sample_sums(NoChloroNoMito))
sd(phyloseq::sample_sums(NoChloroNoMito))





##trying to determine number of ASVs per sample 
physeq_ASVs <- tax_glom(NoChloroNoMito, taxrank= "ASV_ID_2")
ntaxa(physeq_ASVs)

sort(phyloseq::sample_sums(physeq_ASVs))

sort(rowSums(asv_table_ch3 != 0))
#rr=table(tax_table(NoChloroNoMito)[, "Phylum"], exclude = NULL)

##physeq with just water and just larvae
physeq_water= subset_samples(NoChloroNoMito, SampleType== "Water")

physeq_no_larvae= subset_samples(NoChloroNoMito, SampleType== "Food"| SampleType== "FoodB" |SampleType== "Water")
nsamples(physeq_no_larvae)


##transform sample counts to relative abundance
rel_abundances = transform_sample_counts(physeq_no_larvae, function(OTU) OTU/sum(OTU))
ntaxa(rel_abundances)
names(rel_abundances)



##using psmelt to turn my phyloseq object into a data frame
mdf_2 = psmelt(physeq_no_larvae)
names(mdf_2)


##trying alpha diversity 
##code from Sarah on DivNet--- YT002 at the ASV level this time
divnet_test <- mdf_2_matrix %>%
  divnet(X= "SampleID", ncores = 1)
divnet_test


dv<- divnet_test

testDiversity(divnet_test, "shannon")

ests <- sapply(dv[["shannon"]], function(x) x$estimate)
lci <- sapply(dv[["shannon"]], function(x) x$interval[1])
uci <- sapply(dv[["shannon"]], function(x) x$interval[2])

df <- data.frame("SampleID" = names(ests),
                 "h0" = ests, lci, uci,
                 dv$X) 
df


df2=left_join(df, Metadata_table_YT, by="SampleID")



##vibrio plots
vibrio=subset_taxa(rel_abundances, Genus== "g__Vibrio")

ntaxa(vibrio)
mdf_2_vibrio = psmelt(vibrio)


##trying to get only ASVs that make up the top 1%

v2 <-vibrio %>% psmelt()

v2 <-vibrio %>% psmelt() %>% dplyr::group_by(SampleID, SampleType, ASV_ID_2) dplyr::summarise(sum = sum(Abundance))


##ASVs >1% of relative abundance 
V4= v2 %>% dplyr::filter(mean >0.015)

write.csv(mdf_2_vibrio, "mdf_2_vibrio")

NMDS_taxa_color = c("#fa9fb5", "#c51b8a", "#ef3b2c", "#feb24c", "yellow4","#02818a", "#1d91c0", "#253494")

##ordering sample type
mdf_2_vibrio_cs$SampleType <- factor(mdf_2_vibrio_cs$SampleType, levels = c("Larvae","Food","FoodB", "Water"))

##re-naming 
mdf_2_vibrio_cs$SampleType <- factor(mdf_2_vibrio_cs$SampleType, levels = c("larvae", "Food", "FoodB", "water"),
                     labels = c("Larvae", "Tank Food", "Food Stock", "Water"))

##Plotting 
vibrio_ASVs <- ggplot(mdf_2_vibrio_cs, aes(x=Timepoint, y= Abundance, fill= ASV_ID_2)) + facet_wrap(~ SampleType, ncol=8)+ geom_col() + theme(axis.text.x = element_text(size = 20))+theme(axis.title.x = element_text(size= 20), axis.text.x = element_text(angle = 45,size= 12, vjust = 1, hjust=1, face="bold"), axis.title.y=element_text(size= 15), axis.text.y = element_text(size= 18), strip.text = element_text(size=15, face= "bold"), legend.title=element_text(size=15), strip.background = element_rect(colour="black", fill="white"), panel.background = element_rect(fill = "white", colour = "black")) + scale_y_continuous(labels = scales::percent)+ ylab("Proportion of Vibrio")+ xlab("Days Post Hatch")+scale_fill_manual(values = NMDS_taxa_color) +scale_x_continuous(breaks = 1:8, labels = c("2", "5", "8", "13", "20", "36", "55", "76"))
vibrio_ASVs + guides(fill=guide_legend(title="Vibrio ASV"))

mdf_2_vibrio_cs$Timepoint

vibrio_ASVs_all <- ggplot(V4, aes(x=SampleID, y= mean, fill= ASV_ID_2)) + facet_wrap(~SampleType) + geom_col() +theme_bw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y=element_text(size= 20), axis.text.y = element_text(size= 18),plot.subtitle = element_text(size = 22, face="bold"), strip.text = element_text(size=15, face= "bold"), legend.title=element_text(size=15), strip.background = element_rect(colour="black", fill="white"), panel.background = element_rect(fill = "white", colour = "black"), panel.grid.minor.y = element_line()) + scale_y_continuous(labels = scales::percent)+ ylab("Relative Abundance (%)") + scale_fill_manual(values = NMDS_taxa_color) 
vibrio_ASVs_all + guides(fill=guide_legend(title="Vibrio ASV"))


write.csv(V4, "vibrio.csv")

##tree map 

test=tax_glom(rel_abundances, taxrank="Genus", NArm=FALSE)
ntaxa(test)


v2 <-test %>% psmelt() %>% dplyr::group_by(SampleID, Timepoint, SampleType, Genus) %>% dplyr::summarise(sum = sum(Abundance))

V3<-v2%>% dplyr::group_by(Timepoint,SampleType, Genus) %>% dplyr::summarise(mean = mean(sum))


##genera that are greater than 5% of relative abundance 
V4= V3 %>% dplyr::filter(mean >0.04)



##re-naming and creating "other" category" (based on 3% threshold)
V5<-V4 %>% mutate(Genus = case_when(Genus== "Class_c__Bacteroidia" ~ "Unidentified Bacteroidia (Class)",
                                    Genus=="g__Arcobacter" ~ "Arcobacter", 
                                    Genus== "Class_c__Gammaproteobacteria" ~ "Unidentified Gammaproteobacteria (Class)",
                                    Genus== "g__Alteromonas" ~ "Alteromonas",
                                    Genus == "Phylum_p__Proteobacteria" ~ "Unidentified Proteobacteria (Phylum)",
                                    Genus == "Family_f__Alteromonadaceae" ~ "Unidentified Alteromonadaceae (Family)",
                                    Genus== "Family_f__Bacillaceae" ~ "Unidentified Bacillaceae (Family)",
                                    Genus== "Family_f__Flavobacteriaceae" ~ "Unidentified Flavobacteriaceae (Family)",
                                    Genus== "Family_f__Rhodobacteraceae" ~ "Unidentified Rhodobacteraceae (Family)",
                                    Genus == "g__Colwellia" ~ "Colwellia",
                                    Genus == "g__Flavobacterium"~"Flavobacterium",
                                    Genus== "g__Fusobacterium" ~ "Fusobacterium",
                                    Genus== "g__Methylotenera" ~ "Methylotenera",
                                    Genus== "g__MWH-CFBk5" ~ "MWH-CFBk5",
                                    Genus== "g__NS2b_marine_group" ~ "NS2b Marine Group",
                                    Genus== "g__Oleiphilus" ~ "Oleiphilus",
                                    Genus== "g__Polaribacter" ~ "Polaribacter",
                                    Genus== "g__Pseudoalteromonas"~ "Pseudoalteromonas",
                                    Genus== "g__Pseudomonas"~ "Pseudomonas",
                                    Genus== "g__Vibrio" ~ "Vibrio",
                                    
                                    TRUE~"zOther"))



##re-naming and creating "other" category" (larvae phylum)
V5<-V4 %>% mutate(Phylum = case_when(Phylum== "p__Actinobacteriota" ~ "Actinobacteriota",
                                    Phylum=="p__Bacteroidota" ~ "Bacteroidota", 
                                    Phylum== "p__Bdellovibrionota" ~ "Bdellovibrionota",
                                    Phylum== "p__Campilobacterota" ~ "Campilobacterota",
                                    Phylum== "p__Chloroflexi" ~ "Chloroflexi", 
                                    Phylum== "p__Crenarchaeota" ~ "Crenarchaeota",
                                    Phylum== "p__Desulfobacterota" ~ "Desulfobacterota",
                                    Phylum== "p__Firmicutes" ~ "Firmicutes",
                                    Phylum== "p__Fusobacteriota" ~ "Fusobacteriota",
                                    Phylum== "p__Halobacterota" ~ "Halobacterota",
                                    Phylum== "p__Myxococcota"~ "Myxococcota", 
                                    Phylum== "p__Nitrospirota" ~ "Nitrospirota",
                                    Phylum== "p__Planctomycetota" ~ "Planctomycetota",
                                    Phylum== "p__Proteobacteria" ~ "Proteobacteria",
                                    TRUE~"zOther"))

##plotting tree map- genera that make up greater than 5% of the relative abundance

##set colors
all_taxa_color = c("#fa9fb5", "#c51b8a", "#67000d", "#ef3b2c", "#ffffcc", "#feb24c", "#cc4c02", "#fdbb84", "#c7e9b4", "palegreen2","yellow4", "#238b45", "#c6dbef", "#99d8c9","#02818a", "#1d91c0", "#253494", "#6A3D9A", "#54278f","#9e9ac8", "#bdbdbd",  "green1")

all_taxa_color_new= c("#fa9fb5", "#c51b8a", "#67000d", "#ef3b2c", "#ffffcc", "#feb24c", "#cc4c02","#fdbb84", "yellow4", "palegreen2","#c7e9b4", "#238b45","#99d8c9", "#c6dbef","#02818a", "#253494","#9e9ac8", "#bdbdbd")


taxa_phylum= c("#fa9fb5", "#c51b8a","#ffffcc","#feb24c", "#cc4c02","#fdbb84", "yellow4","#c7e9b4", "#238b45","#99d8c9","#c6dbef","#02818a", "#253494","#9e9ac8","#bdbdbd")

ggplot(V5, aes(area= mean, fill= Genus, subgroup=Genus)) +
  geom_treemap(color = "white") +
  geom_treemap_subgroup_border(colour = "white", size = 2) +
  facet_grid(SampleType ~ Timepoint) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 18),
        strip.background = element_blank(), strip.text = element_text(color = "black", face = "bold", vjust = 0.5, hjust = 0, size = 18),
        legend.position = "right",
        legend.title=element_text(size=20),
        legend.text = element_text(color = "black", size = 18),
        panel.border = element_blank(),
        axis.title = element_text(face = "bold", size = 18)) +
  scale_y_continuous(expand = c(0,0)) + labs(x = "", y = "")+ scale_fill_manual(values = all_taxa_color) 


                     #####adonis statistics######
###test SampleType- significant 
CLR_All_Phyto=microbiome::transform(NoChloroNoMito, transform="clr")
metadata <- as(sample_data(CLR_All_Phyto), "data.frame")
dist.uf <- phyloseq::distance(CLR_All_Phyto, method = "euclidean")
ps.adonis <- adonis2(dist.uf ~ SampleType, data = metadata, perm=999, by="margin")
ps.adonis

mod <- betadisper(dist.uf, metadata$SampleType)
mod
anova(mod)


###Food and FoodB significantly different 
CLR_All_Phyto=microbiome::transform(NoChloroNoMito, transform="clr")
metadata <- as(sample_data(CLR_All_Phyto), "data.frame")
dist.uf <- phyloseq::distance(CLR_All_Phyto, method = "euclidean")
ps.adonis <- pairwise.adonis2(dist.uf ~ SampleType, data = metadata, perm=999, by="margin")
ps.adonis


##plotting vibrio ASVs

##subsetting the phyloseq object to only include Vibrio
vibrio=subset_taxa(rel_abundances, Genus== "g__Vibrio")

##21 ASVs
ntaxa(vibrio)

##turn into a data frame to plot
mdf_2_vibrio = psmelt(vibrio)

##Plotting 
vibrio_ASVs<- ggplot(mdf_2_vibrio, aes(x=SampleID, y= Abundance, fill= ASV_ID_2)) + facet_grid(~SampleType, scales="free_x") + geom_col() + theme(axis.text.x = element_text(angle = 90))+theme_bw() + scale_y_continuous(labels = scales::percent)+ ylab("Relative Abundance (%)")+ scale_fill_manual(values = all_taxa_color)
vibrio_ASVs


##trying to get only ASVs that make up the top X%

v2 <-vibrio %>% psmelt() %>% dplyr::group_by(SampleID, Timepoint, SampleType, ASV_ID_2) %>% dplyr::summarise(sum = sum(Abundance))

V3<-v2%>% dplyr::group_by(SampleID, Timepoint, SampleType, ASV_ID_2) %>% dplyr::summarise(mean = mean(sum))

##ASVs >5% of relative abundance 
V4= V3 %>% dplyr::filter(mean >0.01)


##Plotting 
vibrio_ASVs<- ggplot(V4, aes(x=SampleID, y= mean, fill= ASV_ID_2)) + facet_grid(~SampleType, scales="free_x") + geom_col() + theme(axis.text.x = element_text(angle = 90))+ theme(axis.text.x = element_text(angle = 90))+theme(axis.text.x = element_blank(), axis.title.y=element_text(size= 20), axis.text.y = element_text(size= 18), strip.text = element_text(size=15, face= "bold"), legend.title=element_text(size=15), strip.background = element_rect(colour="black", fill="white"), panel.background = element_rect(fill = "white", colour = "black")) + scale_y_continuous(labels = scales::percent)+ ylab("Relative Abundance (%)")+ scale_fill_manual(values = Vibrio_taxa_color)
vibrio_ASVs+ guides(fill=guide_legend(title="Vibrio ASV"))

write.csv(V4, Vibrio_data)


all_taxa_color = c("#f0f0f0","#fa9fb5", "#c51b8a", "#67000d", "#ef3b2c", "#ffffcc", "#feb24c", "#cc4c02", "#fdbb84", "#c7e9b4", "palegreen2","yellow4", "#238b45", "#c6dbef", "#99d8c9","#02818a", "#1d91c0", "#253494", "#9e9ac8", "brown", "#54278f", "#bdbdbd", "#CC79A7", "#252525")
Vibrio_taxa_color = c("#ffffcc", "#fa9fb5", "#c51b8a", "#ef3b2c", "#feb24c", "yellow4","#02818a", "#1d91c0", "#253494")



                   #######NMDS with food, foodB, water######

##creating phyloseq object with just food, foodB 
physeq_larvae= subset_samples(NoChloroNoMito, SampleType== "Larvae" | SampleType=="Larvae")

NMDS_taxa_color = c("#fa9fb5", "#c51b8a", "#ef3b2c", "#feb24c", "yellow4","#02818a", "#1d91c0", "#253494")

##NMDS with CLR and euclidean data transformation
CLR_All_Phyto=microbiome::transform(physeq_larvae, transform="clr")
metadata <- as(sample_data(CLR_All_Phyto), "data.frame")
dist.uf <- phyloseq::distance(CLR_All_Phyto, method = "euclidean")

GP.ord.euclidean <- ordinate(CLR_All_Phyto, "NMDS", "euclidean")
GP.ord.euclidean


NMDS_nolarvae = plot_ordination(CLR_All_Phyto, GP.ord.euclidean, color="Timepoint2") + annotate("text", x=-25, y=70.5, label= "Stress= 0.117", size=5.5, fontface = "bold") + geom_point(size=3.5) + theme_bw() + scale_color_manual(values = NMDS_taxa_color) + guides(color = guide_legend(title = "Timepoint"))+ theme(text = element_text(size = 16), plot.subtitle = element_text(size = 22, face="bold"))+ labs(subtitle = "(A)")
NMDS_nolarvae



test=tax_glom(rel_abundances, taxrank="Genus")

v2 <-test %>% psmelt() %>% dplyr::group_by(SampleID, Timepoint, SampleType, Genus) %>% dplyr::summarise(sum = sum(Abundance))

V3<-v2%>% dplyr::group_by(Timepoint,SampleType, Genus) %>% dplyr::summarise(mean = mean(sum))


##Sigclust with all samples 
library(sigclust2)

mdf_2_matrix = as.matrix(mdf_2)
str(mdf_2_matrix)

otutable_shc_larvae<-as.matrix(as.data.frame(t(otu_table(CLR_All_Phyto))))


shc_result_larvae <- shc(otutable_shc_larvae, metric="euclidean", linkage="ward.D2", n_min=9)

plot(shc_result_larvae)


##sigclust with food A, water, larvae
physeq_food= subset_samples(NoChloroNoMito, SampleType== "Food" | SampleType=="Larvae")
nsamples(physeq_food)



#####core 
####Analyzing core community####

##just timepoints 1-3
physeq_early= subset_samples(NoChloroNoMito, Period== "Early")
nsamples(physeq_early)
                             
##Seeing how many ASVs are shared for early and late larvae stage 

table(meta(physeq_early)$SampleType, useNA = "always")

pseq.rel <- microbiome::transform(physeq_early, "compositional")

early <- unique(as.character(meta(pseq.rel)$SampleType))
print(early)

list_core <- c()


for (n in early){ # for each variable n in Sample Type
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, early == n) # Choose sample from Sample Type by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in at least 90% samples 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

##plotting ven diagram 
mycols <- c(nonCRC="#d6e2e9", CRC="#cbf3f0", H="#fcf5c7") 
plot(venn(list_core),
     fills = mycols)


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




