import os

## 1. generate a dictionary of strains + their unfits
## for each unfit I need start, stop, strand and reason

unfit_dir = "/Users/gh11/klebsiella_TAs/GROUP/unfit_clusters/"
unfit_files = os.listdir(unfit_dir)

strains = {}

for f in unfit_files:
	if not f.endswith(".txt"):
		continue
	with open(os.path.join(unfit_dir, f)) as f_open:
		for line in f_open:
			if line.startswith("Strain") or line.startswith("##"):
				continue
			toks = line.strip().split(",")

			strain = toks[0]

			if strain not in strains:
				strains[strain] = { "unfits" : [], "contigs": [], "starts" : [], "stops": [], "strands":[], "reasons":[]}

			strains[strain]["unfits"].append(toks[1])
			strains[strain]["contigs"].append(toks[5].split("|")[-1])
			strains[strain]["starts"].append(int(float(toks[8])))
			strains[strain]["stops"].append(int(float(toks[9])))
			strains[strain]["strands"].append(toks[6])
			strains[strain]["reasons"].append([toks[11],toks[12]])


## 2. Go over all the extra hits of antitoxins in these strains, and see if the antitoxin
## was found by one of the discarded toxins

antitoxins_out = {}
strains_out = {}


with open("summary_per_strain.csv") as f_open:
	for line in f_open:
		if line.startswith("Name"):
			continue
		toks = line.strip().split(",")
		strain = strains[toks[0]]

		strains_out[toks[0]] = {"Num_with_unfit":0, "Num_without_unfit":0} 

		partners = toks[2].split(";")
		for p in partners:
			toks2 = p.split("(")
			ID = toks2[0]

			if ID not in antitoxins_out:
				antitoxins_out[ID] = {"Num_with_unfit":0, "Num_without_unfit":0, "hit length":0, "Upstream length":0,\
				 "No adjacent upstream ORF":0,"Downstream length":0, "No adjacent downstream ORF":0} 


			toks2 = toks2[1].split(":")
			contig = toks2[0]
			strand = toks2[1]
			start = int(toks2[2])
			stop = int(toks2[3].replace(")",""))
			flag = False
			for i in range(0, len(strain["contigs"])):
				if contig == strain["contigs"][i] and strand == strain["strands"][i] and \
				strain["starts"][i] - 5000 >= start  and start <= strain["starts"][i] + 5000:
					flag = True
					break
			if flag:
				antitoxins_out[ID]["Num_with_unfit"] += 1
				strains_out[toks[0]]["Num_with_unfit"] += 1
				for r in strain["reasons"][i]:
					if r == "":
						continue
					antitoxins_out[ID][r] += 1
			else:
				antitoxins_out[ID]["Num_without_unfit"] += 1
				strains_out[toks[0]]["Num_without_unfit"] += 1




## write output files to load into R
with open("full_summary_per_antitoxin.csv","w") as out:
	out.write("ID, Num_New, Num_with_unfit, Num_without_unfit, hit length, Upstream length, No adjacent upstream ORF, Downstream length, No adjacent downstream ORF\n")
	for ID in antitoxins_out:
		a = antitoxins_out[ID]
		out.write(",".join(map(str,[ID, a["Num_with_unfit"] + a["Num_without_unfit"],  a["Num_with_unfit"],\
			a["Num_without_unfit"], a["hit length"], a["Upstream length"],
			a["No adjacent upstream ORF"], a["Downstream length"], a["No adjacent downstream ORF"]] )) + "\n")


with open("full_summary_per_strain.csv", "w") as out:
	out.write("Name, Num_New, Num_with_unfit, Num_without_unfit\n")
	for ID in strains_out:
		s = strains_out[ID]
		out.write(",".join(map(str, [ID, s["Num_with_unfit"] + s["Num_without_unfit"],\
			s["Num_with_unfit"], s["Num_without_unfit"]])) + "\n")


