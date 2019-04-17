
setwd("/Users/gh11/klebsiella_TAs/GROUP/complete_clusters/")

h = "10H"

metadata = read.table("/Users/gh11/klebsiella_TAs/metadata.csv", header = T, 
                      comment.char = "", sep="," )
kpi = as.character(unlist(metadata$Strain[which(metadata$Phylogroup == "kpi")]))
kpii = as.character(unlist(metadata$Strain[which(metadata$Phylogroup == "kpii")]))
kpiii = as.character(unlist(metadata$Strain[which(metadata$Phylogroup == "kpiii")]))


files = list.files(".",pattern = "*.txt")

tot_kp1 = 0
tot_kp2 = 0
tot_kp3 = 0

for (f in files){
  if (grepl(h,f,fixed = T)) {
    print(f)
    strains = as.character(unlist(read.table(f, sep = ",", comment.char = "", header = T)[,1]))
    num_copies = length(strains)
    num_kpis = length(intersect(kpi,strains)) 
    num_kpiis = length(intersect(kpii,strains)) 
    num_kpiiis = length(intersect(kpiii,strains)) 
    tot_kp1 = tot_kp1 + num_kpis
    tot_kp2 = tot_kp2 + num_kpiis
    tot_kp3= tot_kp3 + num_kpiiis
  }
}


for (f in files){
  if (grepl(h,f,fixed = T)) {
    print(f)
    strains = as.character(unlist(read.table(f, sep = ",", comment.char = "", header = T)[,1]))
    num_copies = length(strains)
    num_kpis = length(intersect(kpi,strains)) 
    num_kpiis = length(intersect(kpii,strains)) 
    num_kpiiis = length(intersect(kpiii,strains)) 
    # print(num_copies)
    # print(num_kpis / tot_kp1)
    # print(num_kpiis / tot_kp2)
    # print(num_kpiiis / tot_kp3)

    print((num_kpiiis + num_kpiis + num_kpis) / 259)
  }
}

