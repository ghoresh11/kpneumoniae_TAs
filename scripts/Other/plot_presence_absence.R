# library(phytools)
# library(RColorBrewer)
# 
# setwd("/Users/gh11/klebsiella_TAs/")
# 
# tree = read.tree(file = "tree/RAxML_bipartitions.tree.out", comment.char = "")
# tree = reroot(tree,node.number = which(tree$tip.label == "5235_5#1"))

## spo = 3.5 width
## toxins = 2 width (height also slightly lower)
## lin = 2.25
## ubiq = 1.3
## choose from "toxins", "ubiq","lineage","spo"
type = "spo"

if (type == "toxins") {
  ### to make the toxin tree
  X = read.table("GROUP/hits.csv", comment.char = "", stringsAsFactors = F,
                 header = T, row.names = 1, sep = ",")
} else {
  ## to make the structure trees
  X = read.table("GROUP/complete.csv", comment.char = "", stringsAsFactors = F,
                 header = T, row.names = 1, sep = ",")
}


## match X to the categories of the toxins
toxin_cats = read.table("/Users/gh11/Submissions/klebsiella_TAs/Supplementary/Tables/S2_toxin_cats.csv",
                        sep = ",", stringsAsFactors = F, comment.char = "", header = T)
toxin_cats$Category = factor(toxin_cats$Category, c("ubiq","lin","spo_amr_vir_plasmid","spo_no_assoc","rare"))
strains = read.table("/Users/gh11/Submissions/klebsiella_TAs/Supplementary/Tables/S1_strains.csv",
                     sep = ",", comment.char = "", stringsAsFactors = F, header = T)


X<-X[tree$tip.label,]
plotTree(tree, plot = F)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
## change here the size of the figure if it gets truncated
plotTree(tree,lwd=1,ylim=c(0,obj$y.lim[2]*1.2),xlim=c(0,obj$x.lim[2]*3.5),
         ftype="off")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
h<-max(obj$xx)
fsize<-0.3

lineage_colors = c("#E69F00","#56B4E9","#009E73")
strains = strains[match(tree$tip.label, strains$Strain),]
strains$Phylogroup[which(strains$Phylogroup == "kpi")] = lineage_colors[1]
strains$Phylogroup[which(strains$Phylogroup == "kpii")] = lineage_colors[2]
strains$Phylogroup[which(strains$Phylogroup == "kpiii")] = lineage_colors[3]
for(i in 1:Ntip(tree)){
  text(h,obj$yy[i],strains$Accession[i],cex=fsize,pos=4,font=3,offset=0.1, col = strains$Phylogroup[i])
}

s<-max(fsize*strwidth(tree$tip.label))
start.x<-1*h+s


### set the colors for the different categories
colors = brewer.pal(n = 10, "Paired")
cols_ubiq<-setNames(c("white",colors[2]),0:1)
cols_lineage<-setNames(c("white",colors[4]),0:1)
cols_sporadic_assoc<-setNames(c("white",colors[6]),0:1)
cols_sporadic_no_assoc<-setNames(c("white",colors[5]),0:1)
cols_rare<-setNames(c("white",colors[8]),0:1)


## get the column names
columns = c()
for (i in 1:ncol(X)){
  curr_toxin = substring(text = strsplit(x = colnames(X)[i], split = ".", fixed = T)[[1]][1], 2)
  columns = c(columns, curr_toxin)
}



## For toxins:
if (type == "toxins") {
  prev = "ubiq"
  toxin_cats = toxin_cats[order(toxin_cats$Category, toxin_cats$ID),]
  
  for (t in 1:nrow(toxin_cats)){
    curr_toxin = toxin_cats$ID[t]
    cat = toxin_cats$Category[t]
    i = which(columns == curr_toxin)
    if (cat == "ubiq") { cols = cols_ubiq}
    else if (cat == "lin") {cols = cols_lineage}
    else if (cat == "spo_amr_vir_plasmid") { cols = cols_sporadic_assoc}
    else if (cat == "spo_no_assoc") { cols = cols_sporadic_no_assoc}
    else {cols = cols_rare}
    
    if (cat != prev) {
      start.x<-start.x+asp
      prev = cat
    }
    
    text(start.x,max(obj$yy)+1,paste(curr_toxin),pos=4,srt=90,
         cex=0.3,offset=0)
    for(j in 1:nrow(X)){
      
      xy<-c(start.x,obj$yy[j])
      y<-c(xy[2]-0.5,xy[2]+0.5,xy[2]+0.5,xy[2]-0.5)
      asp<-(par()$usr[2]-par()$usr[1])/(par()$usr[4]-par()$usr[3])*
        par()$pin[2]/par()$pin[1]+ 0.01
      x<-c(xy[1]-0.5*asp,xy[1]-0.5*asp,xy[1]+0.5*asp,xy[1]+0.5*asp)
      if (X[j,i]>1) { X[j,i] = 1 }
      polygon(x,y,col=cols[as.character(X[j,i])], lwd = 0.05)
    }
    start.x<-start.x+asp
  }
} else {
  ## For structures:
  all_cols = c(brewer.pal(n = 8,"Dark2"), brewer.pal(n = 8,"Dark2"),
                      brewer.pal(n = 8,"Dark2"))
  if (type == "spo") {
    curr_hits = toxin_cats$ID[which(toxin_cats$Category == "spo_amr_vir_plasmid" | 
                                      toxin_cats$Category == "spo_no_assoc")]
  } else {
    curr_hits = toxin_cats$ID[which(toxin_cats$Category == type)]
  }
  cnt = 0
  for (h in curr_hits){
    cnt = cnt + 1
    cols = setNames(c("white",all_cols[cnt]),0:1)
    indexes = which(grepl(x = columns, pattern = paste("_",h,"_",sep=""), fixed = T))
    indexes = c(indexes, which(unlist(lapply(columns, startsWith, h))))
    indexes = c(indexes, which(unlist(lapply(columns, endsWith, paste("_",h,sep="")))))
    
    for (i in indexes) {
      start.x<-start.x+asp
      text(start.x,max(obj$yy)+1,columns[i],pos=4,srt=90,
           cex=0.7,offset=0)
      for(j in 1:nrow(X)){
        
        xy<-c(start.x,obj$yy[j])
        y<-c(xy[2]-0.5,xy[2]+0.5,xy[2]+0.5,xy[2]-0.5)
        asp<-((par()$usr[2]-par()$usr[1])/(par()$usr[4]-par()$usr[3])*
          par()$pin[2]/par()$pin[1]) + 0.1
        x<-c(xy[1]-0.5*asp,xy[1]-0.5*asp,xy[1]+0.5*asp,xy[1]+0.5*asp)
        if (X[j,i]>1) { X[j,i] = 1 }
        polygon(x,y,col=cols[as.character(X[j,i])], lwd = 0.05)
        
      }
    }
    start.x<-start.x+asp
  }
}
