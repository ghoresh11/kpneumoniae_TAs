import os

## Step 1: retrieve the information on the antitoxins and strains from the GROUP results

group_dir = "/lustre/scratch118/infgen/team216/gh11/klebsiella/newGroup_GROUP/partners_clusters/" 
blast_dir = "/lustre/scratch118/infgen/team216/gh11/klebsiella/antitoxins_search/"

group_dir = "/Users/gh11/klebsiella_TAs/GROUP/partners_clusters/" 
blast_dir = "/Users/gh11/klebsiella_TAs/antitoxins_search/klebs_strains/"

antitoxins = {}
strains = {}

antitoxin_clusters = os.listdir(group_dir)

## generate a summary for all the partner clusters and all the strains to compare to the BLAST results
for file in antitoxin_clusters:
	if not file.endswith(".txt"):
		continue
	cluster = file.replace(".txt","")
	antitoxins[cluster] = {"hits": set(), "avg_length":0, "domains":[], "strains":[], "new_strains":[]}

	with open(os.path.join(group_dir, file)) as f_open:
		for line in f_open:
			if "Domains:" in line:
				line = line.strip().split(":")[-1]
				antitoxins[cluster]["domains"] = line.split(",")
				continue
			if "Partner Length" in line:
				line = line.strip().split(":")[-1]
				antitoxins[cluster]["avg_length"] = float(line)
				continue
			if line.startswith("Strain") or line.startswith("##"):
				continue
			toks = line.strip().split(",")
			antitoxins[cluster]["hits"].add(toks[1])
			antitoxins[cluster]["strains"].append(toks[0])

			if toks[0] not in strains:
				strains[toks[0]] = {"new_antitoxins":[]}

			if cluster not in strains[toks[0]]:
				strains[toks[0]][cluster] = {"contigs": [], "lengths": [], "starts":[], "stops":[], "strand":[]}

			strains[toks[0]][cluster]["contigs"].append(toks[5].split("|")[2])
			strains[toks[0]][cluster]["starts"].append(int(float(toks[9])))
			strains[toks[0]][cluster]["stops"].append(int(float(toks[10])))
			strains[toks[0]][cluster]["strand"].append(toks[6])

blast_results = os.listdir(blast_dir)


def add_cluser(antitoxins, cluster, strain, s, contig, start, stop, strand):


	antitoxins[cluster]["new_strains"].append(s)
	strain["new_antitoxins"].append(cluster+"(" + contig + ":" +strand + ":" + str(start) + ":" + str(stop) + ")")

	if cluster not in strain:
		strain[cluster] = {"contigs": [], "lengths": [], "starts":[], "stops":[], "strand":[]}
	
	strain[cluster]["contigs"].append(contig)
	strain[cluster]["starts"].append(start)
	strain[cluster]["stops"].append(stop)
	strain[cluster]["strand"].append(strand)
	return

for file in blast_results:
	if not file.endswith(".tab"):
		continue
	s = file.split(".")[0]

	strain = strains[s]

	with open(os.path.join(blast_dir,file)) as f_open:
		for line in f_open:
			toks = line.strip().split()
			cluster = toks[1].split("_")[0]
			identity = float(toks[2])
			alingment_length = int(toks[3])

			## calculate the protein length:
			toks2 = toks[0].split("|")
			strand = toks2[-3].replace("Strand:","")
			start = int(toks2[-2].replace("Start:",""))
			stop = int(toks2[-1].replace("Stop:",""))
			contig = toks2[5]
			
			length = (stop - start) / 3

			if length < 50 or length > 300 or identity < 75: # keeping with SLING identity threshold, allow for bigger antitoxins
				continue

			if alingment_length < 50:
				continue

			if cluster not in strain:
				## found a new cluster
				add_cluser(antitoxins, cluster, strain, s, contig, start, stop, strand)

			else:
				if contig not in strain[cluster]["contigs"]:
					add_cluser(antitoxins, cluster, strain, s, contig, start, stop, strand)
				else: ## same contig, check start stop and strand
					new = True
					indices = [i for i, x in enumerate(strain[cluster]["contigs"]) if x == contig]
					for i in indices:
						curr_start = strain[cluster]["starts"][i]
						curr_stop = strain[cluster]["stops"][i]
						curr_strand = strain[cluster]["strand"][i]
						if curr_strand == strand and (curr_start <= start + 100 and curr_start >= start - 100) \
						or (curr_stop <= stop + 100 and curr_stop >= stop - 100):
							# print("Stops: %d vs %d" % (curr_stop, stop))
							# print("Starts: %d vs %d" % (curr_start, start))
							# print("Strands: %s vs %s" % (curr_strand, strand))
							new = False
							continue ## it's the same ORF
					if new:
						add_cluser(antitoxins, cluster, strain, s, contig, start, stop, strand)
					

out = open("summary_per_antitoxin.csv","w")
out.write("ID, Domains, Num_copies, Average_Length, Num_new_hits\n")
for c in antitoxins:
	out.write(",".join(map(str,[c, ";".join(antitoxins[c]["domains"]), len(antitoxins[c]["strains"]),\
	 antitoxins[c]["avg_length"], len(antitoxins[c]["new_strains"]) ])) + "\n")
out.close()


out = open("summary_per_strain.csv","w")
out.write("Name, Num_new_hits, IDs\n" )
for s in strains:
	out.write(",".join(map(str, [s, len(strains[s]["new_antitoxins"]), ";".join(strains[s]["new_antitoxins"])])) + "\n")

out.close()
