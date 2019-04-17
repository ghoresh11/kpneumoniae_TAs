### script to read the results and metadata and run association tests between the TA systems and the metadata
library(RColorBrewer)
library(pheatmap)
library(grid)
library(ggplot2)
library(reshape2)


##########################  FUNCTIONS ################################
plot_association <- function(vec1,vec2,name1,name2,p_value,plots) {
  #vec1[which(vec1>0)] = 1
  t=table(vec1,vec2) # create contingency
  png(paste(plots,name1,"_",name2,".png",sep=""))
  tnorm <- t(t)/colSums(t) 
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  pheatmap(tnorm,col = brewer.pal(9, "Blues"),cluster_cols = F,cluster_rows = F)
  setHook("grid.newpage", NULL, "replace")
  grid.text(name1, y=-0.07, gp=gpar(fontsize=16))
  grid.text(name2, x=-0.07, rot=90, gp=gpar(fontsize=16))
  dev.off()
  
  if (name2=="ST" || name2 == "KpI_Lineage") {
    indexes = which(tnorm[,2] > 0.75) # for now, only interested in STs that are extremely common
    write.table(t(rownames(tnorm)[indexes]),paste(plots,name1,"_",name2,".txt",sep=""),col.names = F, row.names = F, quote = F,sep=",")
  }
}

run_assoication_test<- function(vec1,vec2, name, alternative){
  t=t(table(vec1,vec2)) # create contingency
  if (dim(t)[1]<2 || dim(t)[2]<2){ # no varibality in one of the groups
    return(1)
  }
  
  print(dim(t))
  if (dim(t)[1] == 12){
    print(t)
    return (1)
  }
  f=fisher.test(t,workspace=2e+9,alternative = alternative) # run fisher's exact test on the contingency table 
  
  ### for fisher's exact test, for metadata + TAs, we use "g" because I want to see if having a TA increases likelihood 
  ## of property p, for TAs I need a two sided test,
  return(f$p.value) 
}

run_complete_association <- function(df1, df2, type, alternative, symmetric ){
  
  plots = paste("association_plots_",type,"/",sep="")
  dir.create(file.path(plots), showWarnings = F)
  
  print("Running association tests and calculating p-values...")
  p_vals = matrix(0,nrow = dim(df1)[2], ncol = dim(df2)[2])
  
  rownames(p_vals) = colnames(df1)
  colnames(p_vals) = colnames(df2)
  
  for (i in 1:dim(df1)[2]){
    for (j in 1:dim(df2)[2]){
      if (symmetric && i >= j) { ## when running symmetric test, use half the p-values
        next
      }
      name = paste(colnames(df1)[i], colnames(df2)[j], sep=":")
      p_vals[i,j] = run_assoication_test(df1[,i], df2[,j], name, alternative)
      if (symmetric) {
        p_vals[j,i] = p_vals[i,j]
      }
    }
  }
  
  ## correct all the p_values
  print("Correcting for multiple testing...")
  
  
  ## save uncorrected p-values to file for future reference
  write.table(p_vals,paste("uncorrected_pvals_",type,".csv",sep=""), sep = ",",
              quote = F)
  
  ## to not rerun everything again - load the table here
  p_vals = read.table(paste("uncorrected_pvals_",type,".csv",sep=""), sep = ",",
                      header = T, row.names = 1)
  
  
  
  ## extract the relevant p-values from the p-value matrix
  if (symmetric) {
    p_vals_vector = p_vals[lower.tri(p_vals, diag = F)]
  } else {
    # indexes = diag(nrow = dim(p_vals)[1],ncol = dim(p_vals)[2])
    #  p_vals_vector = p_vals[which(indexes == 0),]
    p_vals_vector = as.numeric(unlist(p_vals))
  }
  
  ## correct p-values
  p_vals_vector_corrected = p.adjust(p_vals_vector,method="fdr")
  
  ## get number of significant associations
  significants = length(which(p_vals_vector_corrected<0.01))
  
  ## generate DF for the results
  results = data.frame(matrix(NA, nrow = significants, ncol = 4),stringsAsFactors = F)
  colnames(results) = c("Variable_1","Variable_2","FDR", "Copies")
  print("Writing output file...")
  
  
  indexes = data.frame(matrix(nrow = dim(p_vals)[1], 
                              ncol = dim(p_vals)[2], 
                              1:length(p_vals_vector)))
  index_out=1
  for (k in 1:length(p_vals_vector)){
    if (p_vals_vector_corrected[k] < 0.01){
      index_in = which(indexes == k,arr.ind=TRUE)
      if (symmetric) {
        index_in = index_in[1,]
      }
      i = index_in[[1]]
      j = index_in[[2]]

      results[index_out,] = c(rownames(p_vals)[i], 
                              colnames(p_vals)[j],
                              p_vals_vector_corrected[k],
                              length(which(df1[,i]>0)))
      index_out = index_out + 1
      plot_association(df1[,i],df2[,j],colnames(df1)[i],colnames(df2)[j],p_vals_vector_corrected[k],plots)
    }
  }
  write.table(results,paste("associations_",type,".csv",sep=""),
              quote = F, col.names = T,row.names = F,sep = ",")
}


########################### MAIN ############################

suffix = "complete"
orig_dir = "/Users/gh11/klebsiella_TAs/GROUP/"
outdir = paste(orig_dir,suffix,"_clusters/",sep="")
setwd(outdir)

ta_file = paste(orig_dir,suffix,".csv",sep="")
metadata_file = "/Users/gh11/klebsiella_TAs/metadata.csv"



print("Reading files...")
tas = read.csv(ta_file,row.names = 1)
metadata = read.csv(metadata_file, row.names = 1)
metadata = metadata[match(rownames(tas),rownames(metadata)),] # sort metadata to match matrix

## remove columns that are too complex from analysis -> run this analysis at another time?
metadata = metadata[,-which(colnames(metadata) == "ST")]
metadata = metadata[,-which(colnames(metadata) == "KpI_Lineage")]
metadata = metadata[,-which(colnames(metadata) == "Year")]
# metadata = metadata[,-which(colnames(metadata) == "Country")]
# metadata = metadata[,-which(colnames(metadata) == "Host")]
# metadata = metadata[,-which(colnames(metadata) == "Infection_status")]
metadata = metadata[,-which(colnames(metadata) == "Specimen_type")]

## 1.find associations between metadata and results
## gets stuck at i = 2, because of the STs and lineages... (too many variables)
run_complete_association(tas,metadata,"metadata","g",F)


## 2. find associations between TA systems to see if some always co-occur or are mutually exclusive
### NOT SURE IF THIS WOULD WORK NOW WITH NEW METHOD OF RETRIEVING INDEXES, need to change "indexes" df
# run_complete_association(tas,tas,"TAs", "two.sided",T)




# ###### 3 Look for linear regression with time ###
# years = as.numeric(as.character(metadata$Year))
# num_tests = dim(tas)[2]
# for (t in 1:dim(tas)[2]){
#   vec1 = unlist(tas[,t])
#   fit = lm(vec1 ~ years)
#   res = summary(fit)$r.squared
#   pVal = anova(fit)$'Pr(>F)'[1]
#   if (pVal <= (0.05/num_tests)) {
#     print(colnames(tas)[t])
#     print("Significant!")
#     plot(years, vec1, main = colnames(tas)[t], ylab="Count", xlab="Years")
#   }

