library(ggplot2)
library(RColorBrewer)
library(ggsignif)
library(plyr)

PLOT = F

orig_dir = "/Users/gh11/klebsiella_TAs/GROUP/"
setwd(orig_dir)

lineage_colors = c("#E69F00","#56B4E9","#009E73")

#### Figure 1: Distribution of toxins in Klebsiella spp. ####


toxins = read.csv("hits.csv",row.names = 1, comment.char = "")
antitoxins = read.csv("partners.csv",row.names = 1, comment.char = "")
complete = read.csv("complete.csv",row.names = 1, comment.char = "")


## A. num toxins per strain ###

total_per_strain = rowSums(toxins)

if (PLOT){
  ggplot() + aes(total_per_strain) + geom_histogram(bins=20, fill="#d3d3d3",alpha=0.8) +  theme_bw(base_size = 16) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0.3)) +
    xlab("TA Systems") + ylab("Strains")
}

median(total_per_strain)
mean(total_per_strain)
min(total_per_strain)
max(total_per_strain)

## B. Tree with three different categories ##
### TODO in ITOL ###

## C. PCA plot of the toxins

## remove constant columns:
toxins.pca = toxins[,sapply(toxins, function(v) var(v, na.rm=TRUE)!=0)]
## calculate PCA
toxins.pca = prcomp(toxins.pca,center = TRUE, scale. = TRUE)
summary(toxins.pca)
toxins.pca = data.frame(toxins.pca$x)
## read and match metadata file
metadata = read.csv("/Users/gh11/klebsiella_TAs/metadata.csv", row.names = 1,header = T, comment.char = "")
metadata = metadata[match(rownames(toxins.pca),rownames(metadata)),]

### define colors
toxins.pca = cbind(toxins.pca, Phylogroup =   metadata$Phylogroup)


if (PLOT) {
  ggplot(toxins.pca, aes(x=PC1, y=PC2, color=Phylogroup))  +
    geom_point(size=3, shape=1, stroke =1.5) +  theme_bw(base_size = 16) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    scale_colour_manual(values = lineage_colors, name = "Species",
                        labels = c("K. pneumoniae","K. quasipneumoniae","K. variicola")) +
   # geom_text(aes(label = PC4)) +
    theme(legend.position="none")  
}

### D. Violin plot of number of systems per phylogroup
toxins_per_phylogroup = data.frame(count = total_per_strain,phylogroup = metadata$Phylogroup)

### significance lines
my_comparisons <- list( c("kpii", "kpiii"), c("kpi", "kpiii"))

if (PLOT) {

  ggplot(toxins_per_phylogroup, aes(x=phylogroup,y=count,fill=phylogroup)) + geom_violin()+  theme_bw(base_size = 16) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  +
    scale_fill_manual(values = lineage_colors) + theme(legend.position="none") +
    ylab("TA Systems per Strain") + xlab("Species") + 
    scale_x_discrete(labels=c("K. pneumoniae","K. quasipneumoniae","K. variicola")) +  
    geom_boxplot(width=0.05,fill="white") +
    geom_signif(comparisons = my_comparisons,step_increase=0.1,
                map_signif_level=T) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

### for writing the p-values into the paper
res = pairwise.wilcox.test(toxins_per_phylogroup$count,toxins_per_phylogroup$phylogroup,
                           paired = F,p.adjust.method = "fdr")$p.value

## stratify by the two common countries: USA and Vietnam, and compare frequencies, color by "type" of toxin
metadata = metadata[match(rownames(toxins),rownames(metadata)),]

cols = as.character(unlist(read.table("/Users/gh11/klebsiella_TAs/classification_toxins.csv",sep=",",
                  row.names = 1, stringsAsFactors = F)[,1]))
shapes = as.character(unlist(read.table("/Users/gh11/klebsiella_TAs/amr_and_vir.csv",sep=",",
                                      row.names = 1, stringsAsFactors = F)[,1]))

## colors of lin, rare, spo and ubiq taken from the tree
fill_cols = c("#e7298a","#d3d3d3","#d95f02","#7570b3")

countries = c("vietnam","usa","singapore","indonesia","australia","laos")
par(mfrow = c(6,6), mar=c(2,2,2,2))
for (country1 in countries) {
  for (country2 in countries){
    country_1 = colSums(country_1) / dim(country_1)[1]
    country_2 = colSums(country_2) / dim(country_2)[1]
    
    p = ggplot(mydf , aes(x=  country_1, y=country_2, color = cols,shape = shapes)) + geom_point(size = 3) +
      scale_color_manual(values = fill_cols,guide = F) + xlab(country1) + ylab(country2) + 
     theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    ggsave(filename = paste("country_correlations/",country1,"_",country2,".pdf"),device = "pdf",plot = p)
    
  }
}

## check if there are some STs which are more common in a country
for (country in countries) {
  STs = count(as.character(unlist(metadata$ST[which(metadata$Country == country)])))
  STs =  STs[with(STs, order(freq, decreasing = T)), ]
  STs$x = factor( STs$x , levels = STs$x)
  p = ggplot(STs, aes(x=x, y=freq) ) +
    geom_bar(width = 1, stat = "identity", fill = "#d9d9d9", color = "#737373") +
    theme_bw(base_size = 12) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_continuous(expand = c(0,0)) +
    xlab("ST") + ylab("Count") + ggtitle(country)
  ggsave(filename = paste("country_correlations/",country,"_STs.pdf",sep=""),device = "pdf",plot = p,
         height = 3.5, width = 6)
}


####### ubiq toxins : see if the number of informative SNPs correlates with unique structure
snps_summary = read.table("/Users/gh11/klebsiella_TAs/gene_trees/informative_snps_summary.csv",sep=",",
                          header = T, row.names = 1)
ggplot(snps_summary,aes(x = snp_per_base, y = num_structures, color = class,shape = class)) + 
  geom_point(size = 3,stroke=1.2) +
  scale_color_manual(values = fill_cols) + scale_shape_manual(values=c(0,1,2,5)) +
  xlab("Informative SNPs") + ylab("Unique operon structures") + 
  theme_bw(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = "lm",se=F) 


## E. Curve of how many genomes we samples and how many more toxins/antitoxins/complete operons we see
num_strains = dim(toxins)[1]
curve = data.frame(num_genomes=character(0),type = character(0),average = numeric(0),
                   min = numeric(0), max = numeric(0), stringsAsFactors = F)

# for (n in 1:num_strains){
#   print("Number of Strains:")
#   print(n)
#   values_toxins = c()
#   values_antitoxins = c()
#   values_complete = c()
#   for (i in 1:100){
#     rows = sample(seq(1,num_strains),n)
#     curr_toxins = toxins[rows,]
#     curr_antitoxins = antitoxins[rows,]
#     curr_complete = complete[rows,]
#     unique_toxins = c()
#     unique_antitoxins = c()
#     unique_complete = c()
#     for (r in 1:n){
#       unique_toxins = c(unique_toxins,which(curr_toxins[r,]>0))
#       unique_antitoxins = c(unique_antitoxins,which(curr_antitoxins[r,]>0))
#       unique_complete = c(unique_complete,which(curr_complete[r,]>0))
#     }
#     unique_toxins = length(unique(unique_toxins))
#     unique_antitoxins = length(unique(unique_antitoxins))
#     unique_complete = length(unique(unique_complete))
#     values_toxins =c(values_toxins,unique_toxins)
#     values_antitoxins =c(values_antitoxins,unique_antitoxins)
#     values_complete =c(values_complete,unique_complete)
#   }
#   curve = rbind(curve, c(n,"toxins",mean(values_toxins), sd(values_toxins)), stringsAsFactors = F)
#   curve = rbind(curve,c(n,"antitoxins",mean(values_antitoxins),sd(values_antitoxins)),stringsAsFactors = F)
#   curve = rbind(curve,c(n,"complete",mean(values_complete),sd(values_antitoxins)),stringsAsFactors = F)
# }
# 
# colnames(curve) = c("genomes","type","avg","sd")
# curve$avg = as.numeric(curve$avg)
# curve$genomes = as.numeric(curve$genomes)
# curve$sd = as.numeric(curve$sd)
# curve = cbind(curve, curve$avg - curve$sd) 
# curve = cbind(curve, curve$avg + curve$sd) 
# colnames(curve) = c("genomes","type","avg","sd","min","max")
# 
# 
# # curve=prev_curve
# ## first color: #0072B2
# ## second color: #D55E00
# curve = curve[-which(curve$type == "complete"),]
# 
# if (PLOT) {
#   ggplot(data = curve, aes(x=genomes,y=avg,ymin=min,ymax=max,fill=type,stroke = type)) + geom_line() +
#     geom_ribbon(alpha=0.5) +  theme_bw(base_size = 14) + 
#     theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
#     xlab("Klebsiella genomes") + ylab("Unique gene count") +
#     theme(legend.position="none") + scale_fill_manual(values=c("#4d4d50", "#4d4d50"))
# }
# #write.table(curve,file="accumilation_curve_sd.csv",quote=F,sep=",",row.names = F,col.names = T)
# curve = read.table("accumilation_curve_sd.csv", header = T,sep=",")

### Rare toxins -> do they distribute differently
rare = toxins
#metadata = read.csv("/Users/gh11/Klebsiella/metadata.csv", row.names = 1,header = T, comment.char = "")
#metadata = metadata[match(rownames(rare),rownames(metadata)),]
rare = colSums(toxins)
rare = rare[which(rare <= 8)]
rare = rare[-which(names(rare) == "X43H.HicA_toxin")] ## remove lineage specific rare
write.table(rowSums(toxins[names(rare)]),file="rare.txt",sep=",",quote = F,col.names = F)

### maybe randomly select one strain from each ST





