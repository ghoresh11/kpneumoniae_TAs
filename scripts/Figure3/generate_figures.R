library(ggplot2)
library(RColorBrewer)
library(ggsignif)

setwd("/Users/gh11/klebsiella_TAs/structures/")
summary = read.table("structures_summary_final.csv", sep = ",",
                     header = T, stringsAsFactors = F, comment.char = "")



### sort and refactor DF ####

summary$vals = as.numeric(summary$vals)
summary$num_structs = as.numeric(summary$num_structs)
summary$total = as.numeric(summary$num_structs)

summary$class = factor(summary$class, c("ubiq","lin","spo_amr_vir_plasmid","spo_no_assoc","rare"))
summary = summary[
  order( summary$class, summary$num_structs ),
  ]
summary$ID = factor(summary$ID, unique(summary$ID))



#### boxplot #######

num_structs = summary[match(as.character(unlist(unique(summary$ID))),summary$ID),]
num_structs$class = factor(num_structs$class, c("ubiq","lin","spo_amr_vir_plasmid","spo_no_assoc","rare"))

res = pairwise.wilcox.test(num_structs$num_structs,num_structs$class,
                           paired = F,p.adjust.method = "fdr")$p.value

my_comparisons <- rev(list( c("rare", "spo_no_assoc"),
                        c("spo_amr_vir_plasmid", "rare"),
                        c("lin","rare"),
                        c("ubiq","rare"),
                        c("spo_amr_vir_plasmid","spo_no_assoc"),
                        c("lin", "spo_no_assoc"), 
                        c("ubiq","spo_no_assoc"),
                        c("lin", "spo_amr_vir_plasmid")))

ggplot(data = num_structs, aes(x = class, y = num_structs)) + geom_violin(fill="#d3d3d3")  + theme_bw(base_size = 16) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  + 
  theme(legend.position="none") +
  ylab("Operon structures") + xlab("") +  
  geom_boxplot(width=0.03,fill="white") + 
  scale_x_discrete(labels=c("Ubiquitous","Species associated","Sporadic, with associations","Sporadic, no associations", "Rare")) + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) 
  # geom_signif(comparisons = my_comparisons,step_increase=0.1,
  #             map_signif_level=T) + theme(axis.text.x = element_text(angle = 20, hjust = 1))







###### point plot #########
summary = summary[-which(summary$class == "rare"),]
summary$color[which(summary$color == "#d3d3d3")] = "#a9a9a9"
labels = as.character(unlist(summary$ID[match(as.character(unique(summary$ID)), summary$ID)]))
colors = summary$label_color[match(as.character(unique(summary$ID)), summary$ID)]

  
ggplot(data = summary, aes(x = ID, y = vals, color = color, shape = class)) + geom_point(size = 3.5,  stroke=3) + 
  scale_x_discrete(labels = labels) +
 scale_color_manual(values = c("#a9a9a9", "#f7dfb4" , "#c6e6f8" , "#B0DFD0",  "#e69e01" ,"#54b3e8", "#009c73")) +
  theme_bw(base_size = 16) + 
  theme(panel.border = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  + 
  theme(legend.position="none") + 
 theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = colors)) + theme(panel.grid.minor.y=element_blank(),
                                                                              panel.grid.major.y=element_blank()) +
  xlab("") + ylab("Porportion of copies")





structs_summary = summary
#### run generate figures in sporadic for next section!!!
sporadic_summary = summary
amr_structs = structs_summary$total[match(amr$ID,structs_summary$ID)]
vir_structs = structs_summary$total[match(vir$ID,structs_summary$ID)]
plasmid_structs = structs_summary$total[match(plasmid$ID,structs_summary$ID)]
boxplot(amr_structs,vir_structs,plasmid_structs)
