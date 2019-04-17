library(ggplot2)

base_dir = "/Users/gh11/klebsiella_TAs/identities/"
setwd(base_dir)

idents = seq(35,95,10)

counts = data.frame(stringsAsFactors = F)

types = c("hits", "partners", "complete")

for (i in 1:length(idents)) {
  
  
  ident = idents[i]
  
  ## paste filename -> read table -> get dim -> extract second column
  nums = sapply ( lapply( lapply( lapply(types, paste,"_",ident, ".csv",sep=""),
                                  read.table, sep = ",", header = T, comment.char = "", row.names = 1),
                          dim), "[" , 2)
  counts = rbind(counts, data.frame(rep(ident,3), types, nums) )
}

setwd("/Users/gh11/klebsiella_TAs/GROUP/")
nums = sapply ( lapply( lapply( lapply(types, paste,".csv",sep=""),
                                read.table, sep = ",", header = T, comment.char = "", row.names = 1),
                        dim), "[" , 2)
names(counts) = c("identity", "type", "count")


singletons = clusters[which(clusters$size == 1),]
sng_df = data.frame(name = idents, num_singles = rep(0,length(idents)))
for (i in 1:length(idents)){
  sng_df[i,2]  = length(which(singletons$identity == idents[i]))
}
ggplot(data = sng_df, aes(x=name, y=num_singles)) + geom_line() + geom_line(size = 1.5) +
  xlab("Minimum BLAST identity") + ylab("Number of singletons") + 
  theme_bw(base_size = 16) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_vline(xintercept=75, color = "purple", lty = 3, size = 1)

ggplot(data = counts, aes(x = identity, y = count, color = type)) + geom_line(size = 1.5) +
  xlab("Minimum BLAST identity") + ylab("Number of clusters") + 
  scale_color_brewer("", palette = "Set2", labels = c("Complete", "Toxins", "Antitoxins")) + 
  theme_bw(base_size = 16) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  +
  geom_vline(xintercept=75, color = "purple", lty = 3, size = 1)





setwd(base_dir)
clusters = read.table("r_output.csv", sep = ",", header = T, stringsAsFactors = F)  


get_attributes<- function (id, p,identity, sizes){
  if (identity == 105){
    return (sizes)
  }
  curr = clusters[which( clusters$identity == identity & clusters$parent == p),]
  #  num_clusters = c(num_clusters, dim(curr)[1])
  for (i in 1:dim(curr)[1]){
    # print(identity)
    # print(curr$id[i])
    
    sizes = rbind(sizes, c(id, name = curr$id[i], identity = identity, size = curr$size[i]))
    sizes = get_attributes(id, curr$id[i], identity+10, sizes)
  }
  return(sizes)
}

parents = clusters$id[which(clusters$identity == 35)]
all_sizes = sizes[1,]
for (p in parents){
  sizes = data.frame(id = p, name=p, identity = 35, 
                     size = clusters$size[which(clusters$id == p & clusters$identity == 35)], 
                     stringsAsFactors = F)
  sizes = get_attributes(p, p, 45, sizes)
  sizes$size = as.numeric(sizes$size)
  sizes$identity = as.numeric(sizes$identity)
  plt = ggplot(sizes, aes(x=identity, y = size)) + geom_point()
  ggsave(filename = paste(p,".png", sep=""),plt)
  all_sizes = rbind(all_sizes, sizes)
}
all_sizes = all_sizes[-1,]

df = data.frame(id = "something", identity=35, num_clusters=0, stringsAsFactors = F)
for (p in parents){
  for (i in seq(35,95,10)){
    num = length(which(all_sizes$id == p & all_sizes$identity == i))
    df = rbind(df, c(id = p, identity = i, num_clusters = num))
  }
}
df = df[-1,]
df$identity = as.numeric(df$identity)
df$num_clusters = as.numeric(df$num_clusters)
ggplot(df, aes(x= identity, y=num_clusters, color=id)) + 
  geom_line() +
 geom_vline(xintercept=75, lty=2) + scale_y_continuous(limits = c(0,32))



df = df[-which(df$id %in% c("31H", "36H","27H")),]
ggplot(df, aes(x= identity, y=num_clusters, color=id)) + 
  geom_smooth(method='lm', se=F) +
  geom_vline(xintercept=75, lty=2) + scale_y_continuous(limits = c(0,10))
