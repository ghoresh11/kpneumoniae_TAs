library(ggplot2)
library(vegan)
##################### FUNCTIONS  ##################### 

count_known_vs_novel <- function(df){
  toxins = as.character(unlist(unique(df$Profiles)))
  new_toxins = c()
  for (t in toxins){
    toks = strsplit(t,"/",fixed=T)[[1]]
    new_toxins = c(new_toxins, toks)
  }
  toxins = unique(new_toxins)
  
  out = data.frame(domains = rep(toxins,2),
                   fill = c(rep("known",length(toxins)),
                            rep("novel",length(toxins))),
                   count = rep(0, length(toxins)*2),
                   total = rep(0, length(toxins)*2))
  
  for (i in 1:dim(df)[1]){
    curr_domains = strsplit(as.character(unlist(df$Profiles[i])),
                            "/",fixed=T)[[1]]
    for (d in curr_domains){
      rows = which(out$domains == d)
      if (df$Novel[i] == 0) {
        out$count[rows[1]] =  out$count[rows[1]] + 1
      } else {
        out$count[rows[2]] =  out$count[rows[2]] + 1
      }
      out$total[rows[1]] = out$total[rows[1]] + 1
      out$total[rows[2]] = out$total[rows[1]]
    }
  }
  out = out[order(out$total, decreasing = T),]
  out$domains = factor(out$domains, unique(out$domains))
  return (out)
}


##################### MAIN  ##################### 

setwd("/Users/gh11/klebsiella_TAs/novel_antitoxins/")

######## known/novel per domain, stratified by upstream / downstream 

antitoxins_summary_partial = read.table("novel_summary_new.csv", sep=",",
                                        header = T)

upstream = antitoxins_summary_partial[which(antitoxins_summary_partial$Upstream > 0),]
upstream = count_known_vs_novel(upstream)

downstream = antitoxins_summary_partial[which(antitoxins_summary_partial$Downstream > 0),]
downstream = count_known_vs_novel(downstream)

## add domains with 0 in that orientation
add_zeros <- function(df, df2) {
  to_add = unique(as.character(unlist(df2$domains[which(!df2$domains %in% df$domains)])))
  for (d in to_add){
    add = data.frame(domains = rep(d,2), 
                     fill = c("known","novel"),
                     count = c(0,0),
                     total = c(0,0))
    df = rbind(df, add)
  }
  return(df)
}

upstream = add_zeros(upstream, downstream)
downstream = add_zeros(downstream, upstream)

cols = brewer.pal(3,"Set2")[-1]
cols = c("#d7eef4","#00aad4")

for_order = rbind(upstream, downstream)
or = rev(unique(as.character(unlist(for_order$domains[order(for_order$total)]))))

upstream$domains = factor(upstream$domains, or)
downstream$domains = factor(downstream$domains, or)

ggplot(upstream, aes(x = domains, y = count, fill=fill)) +
  geom_bar(stat='identity', color = "#a9a9a9") + theme_bw(base_size = 14) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(expand = c(0,0)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Antitoxin candidates") + xlab("") + scale_fill_manual(values = rev(cols)) + coord_flip()


ggplot(downstream, aes(x = domains, y = count, fill=fill)) +
  geom_bar(stat='identity',  color = "#a9a9a9") + theme_bw(base_size = 14) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,20)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Antitoxin candidates") + xlab("") + scale_fill_manual(values = rev(cols)) + coord_flip() 


### correlation between num hit clusters and num partner clusters
hits_vs_partners = read.table("num_clusters_per_domain.csv", sep = ",",
                              header = T, stringsAsFactors = F)

ggplot(hits_vs_partners, aes(x = HitClusters, y=PartnerClusters, label = Domain)) + 
  geom_point(colour = "black", alpha = 0.6, size=5, shape = 19) +
  # geom_text(vjust = 0, nudge_y = 0.5) +
  scale_color_gradientn(colors =c(cols[1],"#e9e9e9",cols[2]),limits=c(0,1))+ theme_bw(base_size = 16) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method="glm",se=FALSE, color="#d3d3d3",lty=3) + 
  xlab("Toxin clusters") + ylab("Antitoxin clusters")

l = lm(hits_vs_partners$HitClusters ~ hits_vs_partners$PartnerClusters)


### how does novel vs known fit with the 5 categories

novel_known_cats = read.table("known_novel_per_toxin.csv", sep = ",", header = T)
novel_known_cats$Class = factor(novel_known_cats$Class, c("ubiq","lin","spo_amr_vir_plasmid","spo_no_assoc","rare"))

ggplot(novel_known_cats, aes(x = Class, y = Num, fill = Type))  + geom_boxplot()  + theme_bw(base_size = 16) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  + 
  ylab("Number of antitoxins") + xlab("") +  
  scale_x_discrete(labels=c("Ubiquitous","Species associated","Sporadic, with associations","Sporadic, no associations", "Rare")) + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + scale_fill_manual(values = rev(cols))


# downstream = downstream[match(upstream$domains, downstream$domains),]
# print("More upstream")
# length(which((upstream$total - downstream$total)[seq(1,dim(downstream)[1],2)] > 0))
# print("More downstream")
# length(which((upstream$total - downstream$total)[seq(1,dim(downstream)[1],2)] < 0))
# length(which((upstream$total - downstream$total)[seq(1,dim(downstream)[1],2)] == 0))
# 

## how many in percent are known upstream vs downstream
sum(downstream$count[which(downstream$fill == "known")]) / (sum(downstream$total)/2)
sum(upstream$count[which(upstream$fill == "known")]) / (sum(upstream$total)/2)


### VEGAN
df = data.frame(type = character(0),
                richness = numeric(0),
                genomes = numeric(0),
                sd = numeric(0))
files = c("/Users/gh11/klebsiella_TAs/GROUP/hits.csv",
          "/Users/gh11/klebsiella_TAs/GROUP/partners.csv")

for (f in files) {
  mydata <- read.table(f, header = T, row.names = 1, comment.char = "", quote = "", sep = ",")
  sp <- specaccum(mydata, "random", permutations=100)
  df = rbind(df, data.frame(type = rep(f, length(sp$sites)),
                            richness = sp$richness,
                            genomes = sp$sites,
                            sd = sp$sd))
}
df = cbind(df, min= df$richness-df$sd, max = df$richness+df$sd )
ggplot(df, aes(x = genomes, y = richness, color = type)) + geom_line(size = 2,alpha = 0.6) +
  theme_bw(base_size = 16) + xlab("Genomes") +
  ylab("Genes") + scale_color_manual(values = c("black","black"), guide = F) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=.2,alpha = 0.2,
                position=position_dodge(0.05))
