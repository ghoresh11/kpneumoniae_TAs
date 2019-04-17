library(gridExtra)
library(grid)
library(ggplot2)

setwd("/Users/gh11/klebsiella_TAs/")

categories = read.table("toxins_categories.csv", header = T, sep = ",") 

toxins = read.table("GROUP/hits.csv", comment.char = "", header = T, sep = ",")
metadata = read.table("metadata.csv", comment.char =  "", header = T, sep = ",")

## remove Kp4
metadata = metadata[-which(metadata$Phylogroup == "kpiv"),]

toxins = toxins[match(metadata$Strain,toxins$Strain),]

kp1s = as.character(unlist(metadata$Strain[(which(metadata$Phylogroup == "kpi"))]))
kp2s = as.character(unlist(metadata$Strain[(which(metadata$Phylogroup == "kpii"))]))
kp3s = as.character(unlist(metadata$Strain[(which(metadata$Phylogroup == "kpiii"))]))

porps = data.frame(toxin = character(0),
                   domain = character(0),
                   class = character(0),
                   count = numeric(0),
                   phylogroup = character(0),
                   porp = numeric(0))


## count the porportion of each phylogroup having a toxin
for (i in 2:dim(toxins)[2]) {
  
  ## retrive class and ID for each toxin
  toxin_toks = strsplit(x = colnames(toxins[i]), split = ".", fixed = T)[[1]]
  id = strsplit(x = toxin_toks[1], split = "X", fixed = T)[[1]][2]
  domain = toxin_toks[2]
  toxin_class = as.character(unlist(categories$Category[which(categories$Toxin_Cluster == id)]))
  
  ## get the porpotion of isolates having this toxin for each phylogroup
  strains = toxins$Strain[which(toxins[,i] > 0)]
  kp1_porp = length(which(strains %in% kp1s)) / length(kp1s)
  kp2_porp = length(which(strains %in% kp2s)) / length(kp2s)
  kp3_porp = length(which(strains %in% kp3s)) / length(kp3s)
  new_rows = data.frame(toxin = rep(id, 3),
                        domain = rep(domain, 3),
                        class = rep(toxin_class, 3),
                        count = rep(length(strains), 3),
                        phylogroup = c("kpi", "kpii", "kpiii"),
                        porp = c(kp1_porp, kp2_porp, kp3_porp))
  porps = rbind(porps, new_rows)
}

## remove rares from this plot
porps = porps[-which(porps$class == "rare"),]


## plot
plot_heatmap <- function(df, o) {
  df$toxin = factor(df$toxin, o)  
  x_labels = c()
  for (i in o) {
    x_labels = c(x_labels, as.character(unique(df$domain[which(df$toxin == i)])))
  }
  df$phylogroup = factor(df$phylogroup, c("kpi", "kpii","kpiii"))
  p = ggplot(df, aes(x = toxin, y = phylogroup, fill = porp)) + geom_tile(colour = "white") + 
    scale_fill_gradient(low = "#fbf7f5", high = "steelblue", limits = c(0,1)) +
    scale_y_discrete(labels = c( "K. pneumoniae sensu stricto\n(n=222)",
                                 "K. quasipneumoniae\n(n=19)",
                                 "K. variicola\n(n=18)"), expand = c(0,0)) + 
    theme_bw(base_size = 16) + theme(panel.border = element_blank(), 
                                     panel.grid.major = element_blank(), 
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("") + xlab("") + 
    scale_x_discrete(labels = x_labels) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) + coord_flip()
  return(p)
}  


ubiq = porps[which(porps$class == "ubiq"),]
lin = porps[which(porps$class == "lin"),]
spo_assoc = porps[which(porps$class == "spo_amr_vir_plasmid"),]
spo_no_assoc = porps[which(porps$class == "spo_no_assoc"),]

o1 = as.character(unlist(unique(ubiq$toxin[order(-ubiq$count)])))
o2 = c(              "87H",
                     "51H",
                     "33H",
                     "25H",
                     "64H",
                     "37H",
                     "42H",
                     "2H",
                     "9H",
                     "7H",
                     "21H",
                     "26H",
                     "13H",
                     "16H",
                     "17H")
o3 = as.character(unlist(unique(spo_assoc$toxin[order(-spo_assoc$count)])))
o4 = as.character(unlist(unique(spo_no_assoc$toxin[order(-spo_no_assoc$count)])))
plot_heatmap(porps, rev(c(o1,rev(o2),o3,o4)))


domains = as.character(unique(porps$domain))
count = rep(0,length(domains))
names(count) = domains




#### for figure 3C -> make specific examples
### need to make a DF with structure, lineage, freq
complete_clusters = read.table("GROUP/complete.csv", header = T,
                               comment.char = "", sep = ",")

IDs = c("22H", "11H","10H","51H","7H","51H","18H", "42H")
df = data.frame(struct = character(0),
                phylogroup = character(0),
                porp = numeric(0))
labels = read.table("figures/id_to_label.csv",sep=",",header = T)

for (ID in IDs) {
  name = as.character(labels$Figure_name[which(labels$ID == ID)])
  df = rbind(df, data.frame(struct = rep(name,3), porps[which(porps$toxin == ID),c(5,6)]))
  for (i in 1:length(colnames(complete_clusters))){
    f = colnames(complete_clusters)[i]
    toks = strsplit(f, split = "_", fixed = T)[[1]]
    
    if (ID %in% toks || paste("X",ID,sep="") %in% toks){
      strains = as.character(unlist(complete_clusters$Strain[which(complete_clusters[,i] > 0)]))
      kp1_porp = length(which(strains %in% kp1s)) / length(kp1s)
      kp2_porp = length(which(strains %in% kp2s)) / length(kp2s)
      kp3_porp = length(which(strains %in% kp3s)) / length(kp3s)
      
      new_rows = data.frame(struct = rep(f, 3),
                            phylogroup = c("kpi", "kpii", "kpiii"),
                            porp = c(kp1_porp, kp2_porp, kp3_porp))
      df = rbind(df, new_rows)
    }
    
  }
}



ggplot(df, aes(x = struct, y = phylogroup, fill = porp)) + geom_tile(colour = "white") + 
  scale_fill_gradient(low = "#fbf7f5", high = "steelblue", limits = c(0,1)) +
  scale_y_discrete(labels = c( "K. pneumoniae sensu stricto\n(n=222)",
                               "K. quasipneumoniae\n(n=19)",
                               "K. variicola\n(n=18)"), expand = c(0,0)) + 
  theme_bw(base_size = 16) + theme(panel.border = element_blank(), 
                                   panel.grid.major = element_blank(), 
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("") + xlab("") + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) 






