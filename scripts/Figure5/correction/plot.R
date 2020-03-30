library(ggplot2)
require(dplyr)

setwd("/Users/gh11/Submissions/klebsiella_TAs/revision_work1/orphan_toxins/")

files = list.files("/Users/gh11/klebsiella_TAs/GROUP/unfit_clusters/", full.names = T)

with_hit = files[which(grepl(pattern = "(", x = files, fixed = T))]
files = with_hit

contig_lengths = read.table("contigs.csv", sep = ",", comment.char = "", stringsAsFactors = F, quote = "",
                            header = T)

hit_length = data.frame(antitoxin_group = character(0), is_also_hit = numeric(0), mean_length = numeric(0), sd_length = numeric(0))
for (f in files){
  if (!grepl(x = f, pattern = "UNFIT",fixed = T, ignore.case = F)){ next }
  example = read.table(f, header = T, stringsAsFactors = F, comment.char = "%", quote = "", sep = ",")
  most_common_reason = (tail(names(sort(table(example$Reason1))), 1))
  is_also_hit = 0
  if (grepl(pattern = "(", x = f, fixed = T)) {
    is_also_hit = 1 }
  group = strsplit(x = tail(strsplit(x = f, split = "/", fixed = T)[[1]], n=1), split = "UNFIT", fixed = T)[[1]][1]
  
  if (most_common_reason == "hit length") {
    lengths = example$Hit_Length[which(example$Reason1 == "hit length")]
    curr =data.frame(antitoxin_group = group, is_also_hit = is_also_hit, mean_length = mean(lengths), sd_length = sd(lengths), stringsAsFactors = F)
    hit_length = rbind(hit_length, curr)
  }
}

count = 0
for (f in files) {
  if (!grepl(x = f, pattern = "(",fixed = T, ignore.case = F)){ next }
  example = read.table(f, header = T, stringsAsFactors = F, comment.char = "%", quote = "", sep = ",")
  most_common_reason = (tail(names(sort(table(example$Reason1))), 1))
  curr_lengths = contig_lengths$Length[match(example$Contig, contig_lengths$Contig)]
  dists = curr_lengths - example$Hit_Stop
  upstream = c(example$Hit_Start[which(example$Strand == "+")], dists[which(example$Strand == "-")])
  downstream = c(example$Hit_Start[which(example$Strand == "-")], dists[which(example$Strand == "+")])
  if (most_common_reason != "hit length" ){
    if (mean(downstream) < 500 || mean(upstream) < 500) {count = count + 1}
  }
}




import_files <- lapply(with_hit,read.table,stringsAsFactors =FALSE,header = T, comment.char = "%",sep = ",", quote = "")
df <- Reduce(rbind, import_files)
#df = df[grepl(pattern = "UNFIT", x = df$Unfit_Cluster, fixed = F),]
df$Unfit_Cluster[df$Unfit_Cluster == "100UNFIT(104H)"]=100
df$Unfit_Cluster[df$Unfit_Cluster == "14UNFIT(142H)"]=14
df$Unfit_Cluster[df$Unfit_Cluster == "29UNFIT(20H)"]=29
df$Unfit_Cluster[df$Unfit_Cluster == "47UNFIT(45H)"]=47

## fix reasons
reasons = paste(df$Reason1, "and",df$Reason2)
unique(reasons)
reasons[reasons %in% c("hit length and ", "hit length and NA")] = "Toxin length"
reasons[reasons %in% c("Downstream length and Upstream length", "Upstream length and Downstream length")] = "Adjacent CDSs' lengths"
reasons[reasons %in% c("No adjacent upstream ORF and No adjacent downstream ORF","No adjacent downstream ORF and No adjacent upstream ORF")] = "No adjacent CDSs"
reasons[reasons %in% c("Downstream length and No adjacent upstream ORF",  "No adjacent upstream ORF and Downstream length")] =
  "No upstream CDS and downstream CDS length"
reasons[reasons %in% c("Upstream length and No adjacent downstream ORF","No adjacent downstream ORF and Upstream length") ] = 
  "No downstream CDS and upstream CDS length"

count = data.frame(table(df$Unfit_Cluster,reasons))
count2 = data.frame( group = unique(count$Var1),
                     most_common_reason = rep("", length(unique(count$Var1))), stringsAsFactors = F)
for (t in count2$group) {
  curr = count[which(count$Var1 == t),]
  curr = curr[order(curr$Freq, decreasing = T),]
  count2$most_common_reason[which(count2$group == t)] = as.character(curr$reasons[1])
}
count3 = count2[count2$most_common_reason !="Toxin length",]


count2 = count %>% group_by(Var1, .drop = F) %>% summarise(Value = max(Freq))

count %>% 
  group_by(Var1, reasons) %>% 
  slice(which.max(EstMax))

## order 
for_order = data.frame(table(df$Unfit_Cluster))
for_order = for_order[order(for_order$Freq, decreasing = T),]
for_order = cbind(for_order, domain = paste(df$Hit_Clusters[match(for_order$Var1,df$Unfit_Cluster)],
                                            "(",df$Domain[match(for_order$Var1,df$Unfit_Cluster)] ,")", sep = ""))
count$Var1 = factor(count$Var1 , for_order$Var1)




ggplot(count, aes(y = Freq, x = Var1, fill = reasons)) + geom_bar(stat = "identity", color = "black") +
  scale_x_discrete(labels = for_order$domain) + scale_fill_brewer(palette = "Set3") +
  theme_classic(base_size = 16)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "bottom")+guides(fill=guide_legend(nrow=3,byrow=TRUE))


	
