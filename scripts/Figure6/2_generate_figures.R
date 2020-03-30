library(ggplot2)
library(RColorBrewer)
library(VennDiagram)
library(dplyr)
library(reshape2)
library(ggsignif)

setwd("/Users/gh11/klebsiella_TAs/sporadic/")

summary = read.table("on_contig_summary_new.csv", sep = ",", header = T,stringsAsFactors = F)
## remove the phage inputs
phage = c("116H", "144H", "55H", "79H", "140H", "46H", "75H", "80H")
summary = summary[-which(summary$ID %in% phage),]

num_domains = length(unique(summary$Domain))

domains = read.table("../figures/id_to_label.csv", sep = ",", header = T)
summary = cbind(summary, Domain = 
                as.character(unlist(domains$Figure_name[match( summary$ID, domains$ID)])))

draw_barplot <- function(df, to_remove, pallete){
  ## if adding specials -> change their name
  if (length(specials)>0){
    for (i in 1:length(specials)){
      df$Domain[length(df$Domain)-i+1] = paste("(",as.character(unlist(df$Domain[length(df$Domain)-i+1])),"*)",sep="")
    }
  }
  to_remove = to_remove - 1
  df = df[,-to_remove]
  o = order(df$Count,decreasing = T)
  or = as.character(unlist(df$ID[o]))
  ls = as.character(unlist(df$Domain[o]))
  

  
  no = as.numeric(unlist(df[,2])) - as.numeric(unlist(df[,3]))
  for (i in 1:length(no)){
    if (no[i]<0){
      no[i] = 0
    }
  }
  
  df = cbind(df,no)
  df = df[,-2]
  df = melt(df, id.vars = c("Domain","ID"))
  
  df$ID = factor(df$ID, or)
  ggplot(df, aes(x=ID, y = value, fill = variable)) + 
    geom_bar(stat="identity", position = "stack") + scale_x_discrete(labels = ls) +
    theme_bw(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = brewer.pal(4,"Greys")[c(4,2)]) +
    xlab("") + ylab("Number of copies") + scale_y_continuous(expand = c(0,0))+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE)
  
  
  
}

draw_barplot_2<- function(df,name,col) {
  col = brewer.pal(4,"Greys")[4]
  colnames(df)[colnames(df) == name] = "var"
  o = order(df$var,decreasing = T)
  or = as.character(unlist(df$ID[o]))
  ls = as.character(unlist(df$Domain[o]))
  df$ID = factor(df$ID, or)
#  df$Plasmid[df$Plasmid>0] = 1
  ggplot(df, aes(x=ID, y = var)) + 
    geom_bar(stat="identity", fill=col) + scale_x_discrete(labels = ls) +
    theme_bw(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Unique gene count") + xlab("")+ scale_y_continuous(expand = c(0,0))+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE)
}
d = data.frame(count(summary,vars=Domain))
letters = c("a","b","c","d","e")
for (i in 1:dim(d)[1]){
  if (d$n[i]>1){
    j=1
    change = which(summary$Domain == d$vars[i])
    for (c in change){
      summary$Domain[c] = paste(summary$Domain[c], " (",letters[j],")",sep="" )
      j = j+1
    }
  }
}


## draw a boxplot to see if having more genes increases chances of association to Inc types
draw_boxplot <- function(df, name, col){
  colnames(df)[colnames(df) == name] = "var"
  df$Plasmid[df$Plasmid>0] = 1
  df$Plasmid = as.character(df$Plasmid)
  my_comparisons = list(c("0","1"))
  ggplot(df, aes(x = Plasmid, y = var)) + geom_violin() +  theme_bw(base_size = 16) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  +
    scale_fill_manual(values = lineage_colors) + theme(legend.position="none") +  
    geom_boxplot(width=0.05,fill="white") +
    geom_signif(comparisons = my_comparisons,step_increase=0.1,
                map_signif_level=F) 

}

#specials = which(summary$AMR==0 & summary$Vir==0 & summary$Plasmid==0)
specials=c()

amr = summary[c(which(summary$AMR>0),specials),]
draw_barplot(amr,c(4,6,7,8,9,10,11),"Oranges")
draw_barplot_2(amr, "AMR", brewer.pal(9,"Oranges")[9])

vir = summary[c(which(summary$Vir>0),specials),] 
draw_barplot(vir, c(4,5,6,8,9,10,11), "Blues")
draw_barplot_2(vir,"Vir", brewer.pal(9,"Blues")[9])

plasmid = summary[c(which(summary$Plasmid>0),specials),] 
draw_barplot(plasmid, c(4,5,6,7,8,10,11),"Greens")
draw_barplot_2(plasmid,"Plasmid", brewer.pal(9,"Greens")[9])

## plasmid Inc is always tops 2 genes





### check the contig lengths
lengths = read.csv("contig_lengths.csv",header = F, sep = ",")
## 
lengths = lengths[which(lengths$V1 %in% summary$ID),]

lengths$V1 = factor(lengths$V1,summary[with(summary, order(summary$Avg_Length, decreasing = F)), ]$ID)
labs = summary$Domain[match(levels(lengths$V1),summary$ID)]


ggplot(lengths, aes(x=V1, y=V2,color = V3)) + geom_boxplot(width=0.8)+ geom_point(size=3, alpha = 0.7)  + theme_bw(base_size = 16) + 
  theme(legend.position="none") +
  ylab("Contig length (log10(bp))") + xlab("") + scale_color_continuous(low ="#a4a4a4", high="#e04078") +
  scale_x_discrete(labels = labs) + theme(axis.text.x = element_text(angle = 20, hjust = 1)) + scale_y_log10()

my_comparisons = list(c("0","1"))
lengths$V3 = as.character(lengths$V3)
ggplot(lengths, aes(x=V3,y=V2))  + geom_violin(fill="#d3d3d3") + scale_y_log10()  + theme_bw(base_size = 16) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  + 
  theme(legend.position="none") +
  ylab("Contig lengths (log10(bp))") + xlab("") +  
  geom_boxplot(width=0.1,fill="white") + 
  scale_x_discrete(labels=c("Different Contig","Same Contig")) + 
  geom_signif(comparisons = my_comparisons,step_increase=0.1,
              map_signif_level=T)


