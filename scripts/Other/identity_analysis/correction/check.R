library(RColorBrewer)
library(ggplot2)
library(gridExtra)

setwd("/Users/gh11/klebsiella_TAs/identities/")


all_files = list.files(".", pattern = "hits_*", full.names = F)

all = data.frame(Var1 = character(0),Freq = numeric(0),identity = character(0), stringsAsFactors = F)

for (curr in all_files) {
  curr_identity = strsplit(strsplit(curr, "_", fixed = T)[[1]][2], split = ".", fixed = T)[[1]][1]
  curr_groups = read.table(curr, header = F, sep = ",", comment.char = "", stringsAsFactors = F, quote = "", row.names = 1)[1,]
  counts_per_domain = data.frame(table(sapply(sapply(curr_groups, strsplit, split = "H-", fixed = T),tail, 1)))
  counts_per_domain = cbind(counts_per_domain, identity = rep(curr_identity, dim(counts_per_domain)[1]))
  all = rbind(all, counts_per_domain)
}

counts_per_domain = counts_per_domain[-which(counts_per_domain$Var1 == "Phage_pRha"),]


## get the klebsiella species
md = read.table("../metadata.csv", sep = ",", comment.char = "", quote = "",
                stringsAsFactors = F, header = T)
md = md[-which(md$Phylogroup == "kpiv"),]
plot_domain <- function(name){
  all_counts = data.frame(identity = character(0), count = numeric(0), values = numeric(0), stringsAsFactors = F)
  for (curr in all_files){
    curr_identity = strsplit(strsplit(curr, "_", fixed = T)[[1]][2], split = ".", fixed = T)[[1]][1]
    curr_hits = read.table(curr, header = T, sep = ",", comment.char = "", stringsAsFactors = F, quote = "", row.names = 1)
    curr_hits = curr_hits[match(md$Strain, rownames(curr_hits)),]
    curr_groups = read.table(curr, header = F, sep = ",", comment.char = "", stringsAsFactors = F, quote = "", row.names = 1)[1,]
    curr_groups = sapply(sapply(curr_groups, strsplit, split = "H-", fixed = T),tail, 1)
    curr_hits = curr_hits[,grepl(x = as.character(curr_groups),pattern = name)]
    curr_counts = data.frame(identity = rep(curr_identity,3),
                             counts = 1:3,
                             values = rep(0,3), stringsAsFactors = F)
    
    curr_hits = data.frame(curr_hits)
    for (i in 1:dim(curr_hits)[2]) {
      cnt = 0
      if (sum(curr_hits[,i][which(md$Phylogroup == "kpi")]) > 0) { cnt = cnt + 1}
      if (sum(curr_hits[,i][which(md$Phylogroup == "kpii")]) > 0) {cnt = cnt + 1}
      if (sum(curr_hits[,i][md$Phylogroup == "kpiii"]) > 0)  {cnt = cnt + 1}
      curr_counts$values[curr_counts$counts == cnt] =  curr_counts$values[curr_counts$counts == cnt] + 1
    }
    # curr_counts$values = curr_counts$values / sum(curr_counts$values)
    all_counts = rbind(all_counts, curr_counts)
  }
  all_counts$counts = factor(all_counts$counts, rev(1:3))
  p = ggplot(all_counts, aes(x = identity, y = values, fill = counts, group = counts)) + geom_line() + 
    geom_point(shape = 21, color = "black", size = 3)+
    scale_fill_manual(values = rev(brewer.pal(3,"Blues")), 
                      labels = rev(c("contaning 1 species", "containing 2 species", "containing 3 species")), name = "")+ 
    theme_classic() + xlab("blastp identity threshold") + guides(fill=guide_legend(nrow=3, byrow = F)) +
    ylab("Number of toxin clusters")  + theme(legend.position = "bottom") +
    scale_y_continuous(expand = c(0,0,0.1,0), limits = c(0, max(all_counts$values) + 1)) + ggtitle(name) 
  return(p)
}

relE = plot_domain("RelE")
plot_domain("ParE")
plot_domain("Gp49")
plot_domain("YdaT")
fic = plot_domain("Fic")
pemk = plot_domain("PemK")
plot_domain("HigB")
plot_domain("CcdB")
grid.arrange(relE, pemk, fic, nrow = 1)

curr_hits = read.table("hits_75.csv", header = T, sep = ",", comment.char = "", stringsAsFactors = F, quote = "", row.names = 1)
curr_hits = curr_hits[match(md$Strain, rownames(curr_hits)),]
num_kp1 = length(which(md$Phylogroup == "kpi"))
num_kp2 = length(which(md$Phylogroup == "kpii"))
num_kp3 = length(which(md$Phylogroup == "kpiii"))
for (i in 1:dim(curr_hits)[2]) {
  cnt = 0
  if (sum(curr_hits[,i][which(md$Phylogroup == "kpi")]) / num_kp1 >= 0.70) { cnt = cnt + 1}
  if (sum(curr_hits[,i][which(md$Phylogroup == "kpii")]) / num_kp2 >= 0.70) {cnt = cnt + 1}
  if (sum(curr_hits[,i][md$Phylogroup == "kpiii"])/num_kp3 >= 0.70)  {cnt = cnt + 1}
  if (cnt == 3) {
    print(colnames(curr_hits)[i])
  }
}
