library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)


setwd("/Users/gh11/Submissions/klebsiella_TAs/revision_work1/antitoxins_functions/")

summary = read.table("summary_all.csv", sep = ",", header = T, 
                     comment.char = "", stringsAsFactors = F, quote = "")
all_antitoxins = read.table("/Users/gh11/klebsiella_TAs/novel_antitoxins/novel_summary_new.csv", sep = ",",
                            header = T, stringsAsFactors = F, comment.char = "", quote = "")
all_antitoxins = all_antitoxins[-which(all_antitoxins$Profiles == "Phage_pRha"),]
all_antitoxins = all_antitoxins[-which(all_antitoxins$Novel == 0),]
all_antitoxins = cbind(all_antitoxins, summary[match(all_antitoxins$ID, summary$name),])


all_antitoxins$name[!is.na(all_antitoxins$name)] = "assigned"
all_antitoxins$name[is.na(all_antitoxins$name)] = "unassigned"

## an additional 11 DO map back to known antitoxins when using interpro-scan -> change the numbers in the rest of the figures
length(which(is.na(all_antitoxins$classification)))

all_cats = data.frame(table(all_antitoxins$classification))
all_cats = all_cats[order(all_cats$Freq, decreasing = T),]

all_antitoxins$classification[which(is.na(all_antitoxins$classification))] = "None"


#all_antitoxins$classification = factor(all_antitoxins$classification, c(as.character(unlist(all_cats$Var1)), "None"))

## factor so the order is correct and colours match as I like
ggplot(all_antitoxins, aes(x = name, fill = classification)) + geom_bar(color = "black") +
  theme_classic(base_size = 14) + xlab("Interpro-scan result") + ylab("Novel antitoxins") +
  scale_fill_manual(values = c(brewer.pal(4, "Set2"), brewer.pal(10,"Paired"), "#d3d3d3"), name = "Function") + scale_y_continuous(expand = c(0,0)) +
theme(legend.position="right") + guides(fill =  guide_legend(ncol=2,byrow=T))


## add to supplementary tables S6 and S8
s6 = read.table("s6.csv", sep = ",", skip = 1, comment.char = "", 
                stringsAsFactors = F, header = T)
s6$Known.Novel[which(s6$Known.Novel == "Novel")] = "No"
s6$Known.Novel[which(s6$Known.Novel == "Known")] = "Yes"
colnames(s6)[which(colnames(s6) == "Known.Novel")] = "In TADB"

classification = all_antitoxins$classification[match(s6$Name,all_antitoxins$ID)]
classification[which(is.na(classification))]="In TADB"
s6 = cbind(s6[,1:8], interpro = classification, s6[,9:11])
table(s6$Pfam.Profile[which(s6$`In TADB` == "No" & s6$interpro == "antitoxin")])
write.table(x = s6, file = "new_s6.csv", sep = ",", quote = F, row.names = F, col.names = T)



s8 = read.table("s8.csv", sep = ",", skip = 1, comment.char = "", 
                stringsAsFactors = F, header = T)
s6 = s6[match(s8$ID, s6$Name),]
s8$Novel.Known[which(s6$`In TADB` == "Yes" | s6$interpro == "antitoxin")] = "Known" 
write.table(x = s8, file = "new_s8.csv", sep = ",", quote = F, row.names = F, col.names = T)


for_text = read.table("new_s6.csv", sep = ",", header = T, comment.char = "", stringsAsFactors = F)
length(which(!for_text$Interpro %in% c("antitoxin", "In TADB"))               )
length(which(for_text$Interpro == "None"))
