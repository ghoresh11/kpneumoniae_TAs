library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggsignif)

setwd("/Users/gh11/klebsiella_TAs/antitoxins_search/")

antitoxins = read.table("summary_per_antitoxin.csv", sep = ",",
                        header = T, stringsAsFactors = F) 

##
ids_to_remove = antitoxins$ID[which(antitoxins$Domains == " Phage_pRha")]
antitoxins = antitoxins[-which(antitoxins$Domains == " Phage_pRha"),]

## for text in paper
length(which(antitoxins$Num_new_hits>0)) ## num found unpaired
length(which(antitoxins$Num_new_hits>0 &
               antitoxins$Num_new_hits< 259/10 )) 
length(which(antitoxins$Num_new_hits>0 &
               antitoxins$Num_new_hits > (259/10 * 8) )) 
sum(antitoxins$Num_new_hits)


## distribution of lone antitoxins per antitoxin cluster
antitoxins.tmp = antitoxins[which(antitoxins$Num_new_hits > 0),]
ggplot(antitoxins.tmp, aes(x = Num_new_hits)) + geom_histogram(bins = 70, fill = "#145298", alpha = 0.8)+
  theme_classic(base_size = 16) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0.01,0)) +
  xlab("Number of orphan copies") + ylab("Antitoxin groups") 

## correlation between length and likelihood of finding the antitoxin?
ggplot(antitoxins, aes(x = Average_Length, y = Num_new_hits)) + geom_point(size = 3, shape = 1) +
  theme_bw(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(expand = c(0,0))

## connection between domain and finding it elsewhere
antitoxins$Domains = factor(antitoxins$Domains, antitoxins$Domains[order(antitoxins$Num_new_hits, decreasing = T)])
ggplot(antitoxins, aes(x = Domains, y = Num_new_hits)) + geom_boxplot()+
  theme_bw(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
## to me this suggests TA systems that were missed due to the cutoffs,
## the only way to validate this is to check if the upstream/downstream genes have toxins
## or to look at the discarded hits


####### summary with more info
antitoxins_full = read.table("full_summary_per_antitoxin.csv", sep = ",",
                             header = T, stringsAsFactors = T)
antitoxins_full = antitoxins_full[match(as.character(antitoxins$ID),
                                        as.character(antitoxins_full$ID)),][,-c(1,2)]
antitoxins = data.frame(antitoxins, antitoxins_full[,c(1,2)])
antitoxins.m = melt(antitoxins, id.vars = c("ID", "Domains","Average_Length","Num_copies","Num_new_hits"))
antitoxins.m$value = as.numeric(antitoxins.m$value)
# antitoxins.m = antitoxins.m[-which(antitoxins.m$Num_new_hits < 30),]
antitoxins.m = antitoxins.m[which(antitoxins.m$ID %in% c("12P","58P","204P","156P","124P","73P","236P","136P","122P")),]
o = order(antitoxins.m$Num_new_hits, decreasing = T)
antitoxins.m$ID = factor(antitoxins.m$ID, 
                         unique(antitoxins.m$ID[o]))

labels = paste (antitoxins.m$Domains[o][seq(1,494/2,2)]," (", unique(antitoxins.m$ID[o]), ")",sep="" )


ggplot(antitoxins.m, aes(x = ID, y = value, fill =  variable)) + geom_bar(stat = "identity")+
  theme_bw(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#2e2e2e","#b0b0b0")) + xlab("Antitoxin Pfam Domain") + ylab("Unpaired antitoxins") +
  scale_x_discrete(labels = labels) +  guides(fill=FALSE)


### Strains analysis with

strains = read.table("full_summary_per_strain.csv", sep = ",",
                     header = T, stringsAsFactors = T, comment.char = "")
metadata = read.table("../metadata.csv", sep = ",", header = T,
                      stringsAsFactors = T, comment.char = "")
metadata = metadata[-which(metadata$Phylogroup == "kpiv"),]
metadata = metadata[match( strains$Name, metadata$Strain),]
strains = data.frame(strains, Phylogroup = metadata$Phylogroup)

my_comparisons = list( c("kpi", "kpii"),c("kpii","kpiii"),  c("kpi","kpiii"))

ggplot(strains, aes(x = Phylogroup, y = Num_without_unfit)) + geom_boxplot() +
  theme_bw(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_signif(comparisons = my_comparisons,step_increase=0.1,
              map_signif_level=T) 

res = pairwise.wilcox.test(strains$Num_new, strains$Phylogroup,
                           paired = F,p.adjust.method = "fdr")$p.value

median(strains$Num_without_unfit[strains$Phylogroup == "kpi"])
median(strains$Num_without_unfit[strains$Phylogroup == "kpii"])
median(strains$Num_without_unfit[strains$Phylogroup == "kpiii"])

res = pairwise.wilcox.test(strains$Num_without_unfit, strains$Phylogroup,
                           paired = F,p.adjust.method = "fdr")$p.value
## I see differences in the number of extra antitoxins between the three groups. The only way to validate if they are
## indeed lone antitoxins which are for the protection of the toxins is to look at the region

my_comparisons = list(c("kpii","kpiii"), c("kpi", "kpii"))

strains.m = melt(strains, id.vars = c("Name", "Num_New","Phylogroup"))
strains.m = strains.m[which(strains.m$variable == "Num_without_unfit"),]
ggplot(strains.m, aes(x = Phylogroup, y = value)) + geom_violin(fill = "#b0b0b0" ) + 
  scale_x_discrete(labels=c("K. pneumoniae","K. quasipneumoniae","K. variicola")) +  
  theme_bw(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ylab("Unpaired antitoxins")+ 
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + xlab("")+  guides(fill=FALSE)+
  geom_signif(comparisons = my_comparisons,step_increase=0.1,
              map_signif_level=T)  +  
  geom_boxplot(width=0.05,fill="white") 


res = pairwise.wilcox.test(strains$Num_without_unfit, strains$Phylogroup,
                           paired = F)$p.value


## stratify by category
cats_toxins = read.table("../toxin_categories_with_conservation.csv", sep = ",",
                  header = T, comment.char = "", stringsAsFactors = F)
antitoxins_pairs = read.table("../complete_categories_new.csv", sep = ",", header =T,
                              comment.char ="", stringsAsFactors = F)
cats = data.frame(ID = antitoxins$ID,
                  cat = rep("", dim(antitoxins)[1]),
                  count = antitoxins$Num_new_hits,
                  stringsAsFactors = F)
for (i in 1:dim(antitoxins_pairs)[1]){
  at1 = antitoxins_pairs$Antitoxin1[i]
  at2 = antitoxins_pairs$Antitoxin2[i]
  t = antitoxins_pairs$Toxin_Cluster[i]
  cat = cats_toxins$Category[which(cats_toxins$Toxin_Cluster == t)]
  if (at1 != "-"){
    cats$cat[which(cats$ID == at1)] = cat
  }
  if (at2 != "-"){
    cats$cat[which(cats$ID == at2)] = cat
  }
}
cats$cat = factor(cats$cat , c("ubiq","lin","spo_amr_vir_plasmid","spo_no_assoc","rare"))

ggplot(cats, aes( x = cat, y = count)) + geom_violin(fill = "#b0b0b0" ) +  
  theme_bw(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ylab("Unpaired antitoxins")+ 
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + xlab("")+  guides(fill=FALSE)+ 
  geom_boxplot(width=0.05,fill="white") +
  scale_x_discrete(labels=c("Ubiquitous","Species associated","Sporadic, with associations","Sporadic, no associations", "Rare")) + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) 
res = pairwise.wilcox.test( cats$count,cats$cat,
                           paired = F,p.adjust.method = "fdr")$p.value

