import math
import os
import numpy as np
## read all the association results into a single easy to read dataframe
## and whilt doing this classify each system into what ever group it should be in


#### FUNCTIONS ####


def get_lables(metadata_file):
	with open(metadata_file) as f:
		for line in f:
			if line.startswith("FIELD_LABELS") or line.startswith("name"):
				toks = line.strip().split(",")
				toks = toks[1:]
				res = []
				for t in toks:
					res.append(t.split(".")[0].replace("-",""))
				return res

def get_counts(counts_file):
	d = {}
	with open(counts_file) as f:
		for line in f:
			if line.startswith("ubiq"):
				continue
			toks = line.strip().split(",")
			ID = toks[0].split(".")[0].replace("X","")
			d[ID] = {}
			d[ID]["count"] = int(toks[1])
			d[ID]["vir"] = []
			d[ID]["amr"] = []
			d[ID]["plasmid"] = []
			d[ID]["contig_lengths"] = []
			d[ID]["contig_col"] = []
	return d

def summarise_association_file(association_file, d, vir, amr, plasmid):
	with open(association_file) as f:
		for line in f:
			if line.startswith("Variable"):
				continue
			toks = line.strip().split(",")
			ID = toks[0].split(".")[0].replace("X","")	
			var = toks[1].split(".")[0]

			orig_var = var 

			if var in vir:
				var = "vir"
			elif var in amr:
				var = "amr"
			elif var in plasmid:
				var = "plasmid"

			fdr = toks[2]
			
			if var not in d[ID]:
				d[ID][var] = []
			d[ID][var].append(orig_var)

def categorise_associations(d):
	for item in d:
		if "Phylogroup" in d[item]:
			d[item]["cat"] = "lin"
		elif len(d[item]["vir"])>0 or len(d[item]["amr"])>0 or len(d[item]["plasmid"])>0:
			d[item]["cat"] = "spo_amr_vir_plasmid"
		elif d[item]["count"] >= 207: # more than 80% of isolates is ubiq
			d[item]["cat"] = "ubiq"
		elif d[item]["count"] >= 25: ## in less than 10% of the isolates
			d[item]["cat"] = "spo_no_assoc"
		else:
			d[item]["cat"] = "rare"


## for each strain and each contig, check if an AMR/Vir/Plasmid is on the same contig as the TA system
## name = hits / partners/ complete
def check_contigs(d, name):

	results_dir = os.path.join("/Users/gh11/klebsiella_TAs/GROUP/", name + "_clusters")

	for item in d:
		d[item]["contigs"] = []
		d[item]["strains"] = []
		with open(os.path.join(results_dir, item + ".txt")) as f:
			for line in f:
				if line.startswith("%") or line.startswith("#") or line.startswith("Strain"):
					continue
				toks = line.strip().split(",")
				d[item]["strains"].append(toks[0])
				d[item]["contigs"].append(toks[5])

		
	for item in d:
		d[item]["on_contig"] = {}
		i = -1
		for strain in d[item]["strains"]:
			i += 1
			with open(os.path.join("/Users/gh11/klebsiella_TAs/plasmid_v_chromosome/summary", strain + "_summary.csv")) as f:
				for line in f:
					toks = line.strip().split(",")
					if toks[0] == d[item]["contigs"][i]:
						
						## toks[8] = AMR_score
						## toks[9] = Vir_score
						## toks[10] = plasmid_Score
						### correcting to only include scores higher than 200
						if (float(toks[8]) < 200):
							toks[5] = "-"
						if (float(toks[9]) < 200):
							toks[6] = "-"
						if (float(toks[10])< 200):
							toks[7] = "-"

						on = ";".join([toks[5],toks[6],toks[7]])



						if on not in d[item]["on_contig"]:
							d[item]["on_contig"][on] = 0
						d[item]["on_contig"][on] += 1
						d[item]["contig_lengths"].append(int(toks[1]))
						if toks[5] == "-" and toks[6]=="-" and toks[7]=="-":
							d[item]["contig_col"].append("0")
						else:
							d[item]["contig_col"].append("1")


		string = ""
		for key in d[item]["on_contig"]:
			string += "(" + key + ":" + str(d[item]["on_contig"][key]) + ") "
		d[item]["on_contig"] = string[:-1]
	return



def write_to_output(d, name):

	out = open(name + "_categories_new.csv","w")
	out.write("Toxin_Cluster, Antitoxin1, Antitoxin2,Count, Category, on_contig, vir, amr, plasmid, avg_contig_length, sd_contig_length\n")
	

	for item in d:
		avg_contig_length = str(np.mean(d[item]["contig_lengths"]))
		sd_contig_length = str(np.std(d[item]["contig_lengths"]))
		end = "," + avg_contig_length + "," + sd_contig_length + "\n"

		if name == "toxins":
			out.write(",".join(map(str,[item , "-", "-", d[item]["count"], d[item]["cat"],d[item]["on_contig"],
				 "/".join(d[item]["vir"]),"/".join(d[item]["amr"]),"/".join(d[item]["plasmid"])])) + end)
		elif name == "antitoxins":
			out.write(",".join(map(str,["-", item, "-", d[item]["count"], d[item]["cat"],d[item]["on_contig"],
				 "/".join(d[item]["vir"]),"/".join(d[item]["amr"]),"/".join(d[item]["plasmid"])])) + end)
		else:
			items = item.split("_")
			if len(items) == 3:
				out.write(",".join(map(str,[items[1], items[0], items[2], d[item]["count"], d[item]["cat"],d[item]["on_contig"],
				 "/".join(d[item]["vir"]),"/".join(d[item]["amr"]),"/".join(d[item]["plasmid"]) ])) + end)
			elif "H" in items[0]:
				out.write(",".join(map(str,[items[0], "-", items[1], d[item]["count"], d[item]["cat"],d[item]["on_contig"],
					"/".join(d[item]["vir"]),"/".join(d[item]["amr"]),"/".join(d[item]["plasmid"])])) + end)
			else:
				out.write(",".join(map(str,[items[1], "-", items[0], d[item]["count"], d[item]["cat"],d[item]["on_contig"],
					"/".join(d[item]["vir"]),"/".join(d[item]["amr"]),"/".join(d[item]["plasmid"])])) + end)
	out.close()
	return


def move_itol_files(d, name):

	itol_path = "/Users/gh11/klebsiella_TAs/GROUP/ITOL/" + name + "_clusters"
	itol_files = os.listdir(itol_path)

	if not os.path.exists(os.path.join(itol_path, "ubiq")):
		os.makedirs(os.path.join(itol_path, "ubiq"))
		os.makedirs(os.path.join(itol_path, "lin"))
		os.makedirs(os.path.join(itol_path, "rare"))
		os.makedirs(os.path.join(itol_path, "spo_no_assoc"))
		os.makedirs(os.path.join(itol_path, "spo_amr_vir_plasmid"))

	for item in d:
		for f in itol_files:
			if f.startswith(item):
				os.rename(os.path.join(itol_path,f), os.path.join(itol_path, d[item]["cat"], f))
				break
	return



#### MAIN ####
vir = get_lables("metadata/vir_for_gal.txt")
amr = get_lables("metadata/AMR_cleaned_forGal.txt")
plasmid = get_lables("metadata/plasmid.csv")

toxins = get_counts("GROUP/hits_counts.csv")
summarise_association_file("GROUP/hits_clusters/associations_metadata.csv", toxins, vir, amr, plasmid)
categorise_associations(toxins)
check_contigs(toxins,"hits")


antitoxins = get_counts("GROUP/partners_counts.csv")
summarise_association_file("GROUP/partners_clusters/associations_metadata.csv", antitoxins, vir, amr, plasmid)
categorise_associations(antitoxins)
check_contigs(antitoxins,"partners")

complete = get_counts("GROUP/complete_counts.csv")
summarise_association_file("GROUP/complete_clusters/associations_metadata.csv", complete, vir, amr, plasmid)
categorise_associations(complete)
check_contigs(complete,"complete")

move_itol_files(complete,"hits")
move_itol_files(antitoxins, "partners")
move_itol_files(complete,"complete")


write_to_output(toxins, "toxins")
write_to_output(antitoxins, "antitoxins")
write_to_output(complete, "complete")

quit()
## write output of all the contig lengths
out = open("contig_lengths.csv","w")
for item in toxins:
	for i in range(0,len(toxins[item]["contig_lengths"])):
		out.write(item + "," + str(toxins[item]["contig_lengths"][i]) + "," +  toxins[item]["contig_col"][i] + "\n")
out.close()

