from collections import Counter


def convert_lists(l):
	d = {}
	if len(l) == 0:
		return d
	for item in l:
		if item == "":
			continue
		item = item.split("_")[0]
		d[item] = 0
	return d

out = open("on_contig_summary_new.csv","w")
out.write("ID, Count, AMR, Num, Vir, Num , Plasmid, Num, Avg_Length, Sd_Length\n")

with open("../toxins_categories_new.csv") as f:
	for line in f:
		if line.startswith("Toxin"):
			continue
		toks = line.strip().split(",")
		cat = toks[4]
		if not "spo" in cat:
			continue



		associated_vir_genes = convert_lists(toks[6].split("/"))
		associated_amr_genes = convert_lists(toks[7].split("/"))
		associated_plasmid_genes = convert_lists(toks[8].split("/"))

		
		on_contig = toks[5].split()
		for v in on_contig:
			v = v.replace("(","")
			v = v.replace(")","")
			v = v.split(";")
			v = v[0:2] + v[2].split(":")

			for i in range(0,3):
				v[i] = v[i].split("_")[0]
				if v[i] != "-":
					v[i] = v[i].split("-")[0]
				v[i] = v[i].split("/")[0]
			
			for gene in associated_amr_genes:
				if gene.startswith(v[0]):	
					associated_amr_genes[gene] += int(v[3])
			
			for gene in associated_vir_genes:
				if gene.startswith(v[1]):	
					associated_vir_genes[gene] += int(v[3])
		
			for gene in associated_plasmid_genes:
				if gene.startswith(v[2]):	
					associated_plasmid_genes[gene] += int(v[3])


		out.write(toks[0] + "," + toks[3] + ",")
		total = float(toks[3])


		out.write( str(len(associated_amr_genes)) + "," + str(sum(associated_amr_genes.values())) + ",")		
		out.write(str(len(associated_vir_genes)) + "," + str(sum(associated_vir_genes.values())) + ",")
		out.write(str(len(associated_plasmid_genes)) + "," + str(sum(associated_plasmid_genes.values())) + "," + toks[-2] + "," + toks[-1] + "\n")

out.close()	





