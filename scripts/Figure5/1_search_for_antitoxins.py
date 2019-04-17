import os
import subprocess

## given all the antitoxin sequences
## get all of them -> create a blast DB that is all the antitoxins
## using the same BLAST cutoffs used in SLING
## blast the genomes against the antitoxins DB and see if I find any


## step 1 -> aggregate all antitoxin sequences
antitoxins_dir = "/lustre/scratch118/infgen/team216/gh11/klebsiella/newGroup_GROUP/partners_clusters/"
#antitoxins_dir = "/Users/gh11/klebsiella_TAs/GROUP/partners_clusters/"
antitoxin_files = os.listdir(antitoxins_dir)

antitoxins = open("all_antitoxins.fasta", "w")

for f in antitoxin_files:
	if not f.endswith(".txt"):
		continue
	name = f.replace(".txt","")
	line_num = 1
	with open(os.path.join(antitoxins_dir,f)) as f_open:
		for line in f_open:
			if line.startswith("#") or line.startswith("Strain"):
				continue
			toks = line.strip().split(",")
			antitoxins.write(">" + name + "_" + str(line_num) + "\n" + toks[-3] + "\n")
			line_num += 1

antitoxins.close()

## run cd-hit to remove identical sequences in DB
subprocess.call(["cd-hit","-i","all_antitoxins.fasta",\
		"-o","all_antitoxins_clustered.fasta","-c","0.9","-d","0"])

## build a BLAST database for this
subprocess.call(["makeblastdb","-in","all_antitoxins_clustered.fasta","-dbtype","prot"])


## run BLAST against each genome
proteins_dir = "/lustre/scratch118/infgen/team216/gh11/klebsiella/fixed_PREPARE/"
#proteins_dir = "/Users/gh11/klebsiella_TAs/functional_validation/antitoxins_search/lab_strains/proteins_PREPARE/"
proteins = os.listdir(proteins_dir)
for f in proteins:
	if not f.endswith(".fasta"):
		continue
	name = f.replace(".fasta","")
	subprocess.call(["blastp", "-db","all_antitoxins_clustered.fasta", "-query", \
	os.path.join(proteins_dir, f), "-out", name + "_blast.tab","-outfmt",\
	"6 qseqid sseqid pident length slen qlen","-evalue","0.01",\
	 "-num_threads", "16"])




