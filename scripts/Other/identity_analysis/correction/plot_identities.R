library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(pheatmap)
library(grid)

setwd("/Users/gh11/klebsiella_TAs/identities/")
RES = "res3"


analyse_one_identity <- function(curr_file){
  curr_ident = strsplit(strsplit(curr_file, split = ".", fixed = T)[[1]][1], split = "_", fixed = T)[[1]][2]
  curr_hits = read.table(curr_file, header = T, row.names = 1, stringsAsFactors = F, comment.char = "", sep = ",")
  curr_hits = curr_hits[,-which(grepl(x = colnames(curr_hits), pattern = "Phage", fixed = T))]
  num_toxins = dim(curr_hits)[2]
  phylogroups = md$Phylogroup[match(rownames(curr_hits), md$Strain)]
  pvals_species_associated = c()
  freqs = colSums(curr_hits)/dim(curr_hits)[1]
  for (i in 1:dim(curr_hits)[2]){
    t=table(curr_hits[,i],phylogroups)
    if (dim(t)[1]<2 || dim(t)[2]<2 || colSums(curr_hits)[i] <= 4){ # no varibality in one of the groups, very rare
      pvals_species_associated = c(pvals_species_associated, 1)
      next
    }
    res = fisher.test(t,workspace=2e+9)$p.value
    tnorm <- t(t)/colSums(t) 
    # if (i == 3){
    # setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
    # pheatmap(tnorm,col = brewer.pal(9, "Blues"),cluster_cols = F,cluster_rows = F)
    # setHook("grid.newpage", NULL, "replace")
    # }
    pvals_species_associated = c(pvals_species_associated, res)
  }
  pvals_species_associated = p.adjust(pvals_species_associated, method = "fdr")
  species =colnames(curr_hits)[which(pvals_species_associated < 0.01)]
  ubiq = colnames(curr_hits)[which(freqs>0.8 & !(names(freqs) %in%  species))]
  spo = colnames(curr_hits)[ which(colSums(curr_hits) > 26 & !(colnames(curr_hits) %in% c(species, ubiq))) ]
  singletons = length(which(colSums(curr_hits) == 1))
  res = data.frame(identity = rep(curr_ident, 5),
                   type = c("ubiq", "lin", "spo", "rare", "singletons"),
                   count = c(length(ubiq), length(species), length(spo), 
                             dim(curr_hits)[2] - sum(length(ubiq), length(species), length(spo)), singletons), stringsAsFactors = F)
  res2 = data.frame(identity = rep(curr_ident, length(ubiq)),
                    ubiq = ubiq,
                    domain = as.character(sapply(sapply(ubiq, strsplit, split = ".", fixed = T), tail,1)),
                    count = colSums(curr_hits)[which(colnames(curr_hits) %in% ubiq)], stringsAsFactors = F)
  res3 = data.frame(identity = rep(curr_ident, length(species)),
                    species = species,
                    domain = as.character(sapply(sapply(species, strsplit, split = ".", fixed = T), tail,1)),
                    count = colSums(curr_hits)[which(colnames(curr_hits) %in% species)], stringsAsFactors = F)
 all = data.frame(identity = rep(curr_ident, dim(curr_hits)[2]),
                  domain_ID = colnames(curr_hits),
                  domain = as.character(sapply(sapply(colnames(curr_hits), strsplit, split = ".", fixed = T), tail,1)),
                  count = colSums(curr_hits), stringsAsFactors = F)

  if (RES == "res") {return(res)}
  if (RES == "res2") {return(res2)}
  if (RES == "res3") {return(res3)}
 if (RES == "all") {return(all)}
}


## B = how many ubiq/spo/lin/rare using different cut-offs
identity_files = list.files(path = '.', pattern = "hits*")
md = read.table("../metadata.csv", sep = ",", header = T, comment.char = "", stringsAsFactors = F)
dataframes = lapply(FUN = analyse_one_identity, identity_files)
if (RES == "res") {
  results = do.call(rbind, dataframes)
}
if (RES == "res2") {
  ubiq = do.call(rbind, dataframes)
}
if (RES == "res3") {
  species = do.call(rbind, dataframes)
}
if (RES == "all") {
  all = do.call(rbind,dataframes)
}

### for res 
if (RES == "res") {
  results$type = factor(results$type, c("ubiq", "lin","spo", "rare","singletons"))
  A = ggplot(results, aes(x = type, y = count, fill = identity)) + geom_bar(stat = "identity", position = "dodge", color = "black")+
    scale_fill_brewer(palette = "Purples", name = "blastp identity") + theme_bw(base_size = 14) + scale_y_continuous(expand = c(0,0, 0.05, 0)) +
    ylab("Number of toxin clusters") + xlab("") + scale_x_discrete(labels = c("Ubiquitous", "Species associated", "Sporadic","Rare","Singletons")) +
    theme(legend.position = "bottom")
}

if (RES == "res2") {
  ## for res2
  for (d in unique(ubiq$domain)) {
    for (i in unique(ubiq$identity)) {
      if (length(which(ubiq$domain == d & ubiq$identity == i)) == 0) {
        ubiq = rbind(ubiq, data.frame(identity = i,
                                      ubiq = "X",
                                      domain = d,
                                      count = 0,
                                      stringsAsFactors = F))
      }
    }
  }
  B =  ggplot(ubiq, aes(x = domain, y = count, fill = identity)) + geom_bar(stat = "identity", position = "dodge", color = "black") +
    scale_fill_brewer(palette = "Purples", name = "blastp identity ") + theme_classic(base_size = 14) + xlab("Toxin profile") + ylab("Count") +
    theme(legend.position = "bottom") + scale_y_continuous(expand = c(0,0,0.05,0))
}

if (RES == "res3") {
  ## for res2
  ## break it into the domains
  species$domain[which(species$domain %in% c("AntA","N"))] = "BroN"
  species$domain[which(species$domain %in% c("like_toxin"))] = "ParE_1"
  species$domain[which(species$domain %in% c("ParE_toxin"))] = "ParE_2"
  
  species$domain[which(species$domain == 'BroN' & species$count < 5)] = "BroN_1"
  species$domain[which(species$domain == 'BroN' & species$count > 5)] = "BroN_2"
  
  species$domain[which(species$domain == 'Fic' & species$count < 6)] = "Fic_1"
  species$domain[which(species$domain == 'Fic' & species$count %in% c(236,19))] = "Fic(ubiq)"
  species$domain[which(species$domain == 'Gp49' & species$count == 40)] = "Gp49_1"
  species$domain[which(species$domain == 'Gp49' & species$count != 40)] = "Gp49_2"
  species$domain[which(species$domain == 'RelE' & species$count %in% c(17,26))] = "RelE_1"
  species$domain[which(species$domain == 'RelE' & species$count == 26)] = "RelE_2"
  species = species[-which(species$domain == "RES"),]##sporadic because associated with virulence genes
  for (d in unique(species$domain)) {
    for (i in unique(species$identity)) {
      if (length(which(species$domain == d & species$identity == i)) == 0) {
        species = rbind(species, data.frame(identity = i,
                                            species = "X",
                                            domain = d,
                                            count = 0,
                                            stringsAsFactors = F))
      }
    }
  }
  
    C = ggplot(species, aes(x = domain, y = count, fill = identity)) + geom_bar(stat = "identity", position = "dodge", color = "black") +
      scale_fill_brewer(palette = "Purples", name = "blastp identity ") + theme_bw(base_size = 14) + xlab("Toxin profile") + ylab("Count") +
      theme(legend.position = "bottom") + scale_y_continuous(expand = c(0,0,0.05,0))
  
}



count_per_lineage <- function(curr_file){
  curr_ident = strsplit(strsplit(curr_file, split = ".", fixed = T)[[1]][1], split = "_", fixed = T)[[1]][2]
  curr_hits = read.table(curr_file, header = T, row.names = 1, stringsAsFactors = F, comment.char = "", sep = ",")
  curr_hits = curr_hits[,-which(grepl(x = colnames(curr_hits), pattern = "Phage", fixed = T))]
  num_toxins = dim(curr_hits)[2]
  phylogroups = md$Phylogroup[match(rownames(curr_hits), md$Strain)]
  df = data.frame(ID=character(0),identity=numeric(0),lin=character(0),counts=numeric(0))
  for (i in 1:dim(curr_hits)[2]){
    counts = c()
    for (phy in c("kpi","kpii","kpiii")){
      counts = c(counts, sum(curr_hits[,i][which(phylogroups == phy)]))
    }
    df = rbind(df, data.frame(ID = rep(colnames(curr_hits)[i],3),
                              identity = rep(curr_ident, 3),
                              lin = c("kpi","kpii","kpiii"),
                              counts = counts, stringsAsFactors = F))
  }
  return(df)
}

lineage_colors = c("#E69F00","#56B4E9","#009E73")

species_name = paste(species$species, species$identity, sep = "")
ubiq_names = paste(ubiq$ubiq, ubiq$identity, sep = "")

fic_ubiq = ubiq_names[which(grepl(x = ubiq$ubiq,pattern = "Fic", fixed =T))]
fic_ubiq = c(fic_ubiq,species_name[which(species$domain == "Fic(ubiq)")])

## worked hard to get the IDs of missing RelE and PemK in thresholds 85 and 95
rele = species_name[which(species$domain == "RelE_1")]
rele = c(rele, "X33H.RelE85","X34H.RelE95")

pemk = species_name[which(species$domain == "PemK_toxin")]
pemk = c(pemk, "X119H.PemK_toxin95","X57H.PemK_toxin85")

dataframes2 = lapply(FUN = count_per_lineage, identity_files)
per_phylo = do.call(rbind, dataframes2)
per_phylo = cbind(per_phylo, name = paste(per_phylo$ID, per_phylo$identity, sep = ""))


fic = per_phylo[which(per_phylo$name %in% fic_ubiq),]
fic$name = factor(fic$name, unique(fic$name[order(fic$identity)]))
D = ggplot(fic, aes(x = name, y = counts, fill = lin)) + geom_bar(stat = "identity", color = "black") + ggtitle("D") +labs(subtitle="Fic (ubiq)")+
  scale_fill_manual(values = lineage_colors, guide = F) + theme_classic(base_size = 14) +
  scale_x_discrete(labels = c("35","45","55","65","75","85","95(a)","95(b)")) + scale_y_continuous(expand = c(0,0,0.05,0)) +
  xlab("blastp identity") + ylab("Isolates")


rele = per_phylo[which(per_phylo$name %in% rele),]
rele$name = factor(rele$name, unique(rele$name[order(rele$identity)]))
E = ggplot(rele, aes(x = name, y = counts, fill = lin)) + geom_bar(stat = "identity", color = "black") + ggtitle("E")+labs(subtitle="RelE_1")+
  scale_fill_manual(values = lineage_colors, guide = F) + theme_classic(base_size = 14)  + scale_y_continuous(expand = c(0,0,0.05,0)) +
  scale_x_discrete(labels = c("35","45","55","65","75","85(a)","85(b)","95(a)","95(b)"))+
  xlab("blastp identity") + ylab("Isolates")


pemk = per_phylo[which(per_phylo$name %in% pemk),]
pemk$name = factor(pemk$name, unique(pemk$name[order(pemk$identity)]))
F_plot = ggplot(pemk, aes(x = name, y = counts, fill = lin)) + geom_bar(stat = "identity", color = "black") + ggtitle("F") +labs(subtitle="PemK")+
  scale_fill_manual(values = lineage_colors, guide = F) + theme_classic(base_size = 14)  + scale_y_continuous(expand = c(0,0,0.05,0)) +
  scale_x_discrete(labels = c("35","45","55","65","75","85(a)","85(b)","95(a)","95(b)"))+
  xlab("blastp identity") + ylab("Isolates")


## add the examples from the other analysis
A_no_legend =  A + theme_classic(base_size = 14)+ theme(legend.position = "None") + ggtitle("A") + ylab("Isolates")
B_no_legend =  B + theme_classic(base_size = 14)+ theme(legend.position = "None") + ggtitle("B") + ylab("Isolates")
C_no_legend =  C + theme_classic(base_size = 14)+ theme(legend.position = "None") + ggtitle("C") + ylab("Isolates")

lay = rbind(c(1,1,1),
            c(2,2,2),
            c(3,3,3),
            c(4,5,6))
grid.arrange(A_no_legend,B_no_legend,C_no_legend, D, E, layout_matrix = lay)
## save w:1000, heigh = ~900


