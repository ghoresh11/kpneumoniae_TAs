import os


out = open("novel_summary_new.csv","w")
out.write("ID,Hit_clusters,Profiles,Copies,Upstream,Downstream,Novel\n")


partners = {}

hit_clusters = {}

for file in os.listdir("/Users/gh11/klebsiella_TAs/GROUP/partners_clusters/"):
	if not file.endswith(".txt"):
		continue
	with open("/Users/gh11/klebsiella_TAs/GROUP/partners_clusters/" + file) as f:
		partner_id = file.split("/")[-1]
		partner_id = partner_id.replace(".txt","")
		hits = set()
		for line in f:
			if line.startswith("##  Num Copies:"):
				num_copies = line.strip().split()[-1]
				continue
			if line.startswith("##  Domains:"):
				domains = line.strip().split()[-1]
				domains = domains.replace(",","/")
				continue
			if line.startswith("##  Order:"):
				order = line.strip().split()
				downstream = order[-1]
				upstream = order[-3]
				continue
			if line.startswith("##") or line.startswith("Strain"):
				continue
			toks = line.strip().split(",")

			hits.add(toks[1])
			if toks[1] not in hit_clusters:
				hit_clusters[toks[1]] = {"size" : 0, "domains" : set(), "partners" : set()}
			
			hit_clusters[toks[1]]["size"] += 1
			hit_clusters[toks[1]]["domains"].add(toks[2])
			hit_clusters[toks[1]]["partners"].add(partner_id)

		orientation = "upstream"
		if int(downstream) > 0 and int(upstream) > 0:
			orientation = "both"
		elif int(downstream) > 0:
			orientation = "downstream"
		partners[partner_id] = {"num_copies": num_copies, "hits": ";".join(list(hits)), "domains": domains, 
		"downstream": downstream, "upstream": upstream, "orientation": orientation}



for file in os.listdir("."):
	if not file.endswith("_result"):
		continue
	partner_id = file.replace(".fasta_result","")
	novel = 1
	with open(file) as f:
		for line in f:
			toks = line.strip().split()
			if float(toks[2]) > 75:  ## increased to 75% identity because of antitoxins only sharing a domain
				novel = 0
				break
	d = partners[partner_id]
	d["novel"] =  "Novel" if novel == 1 else "Known"
	out.write(",".join([partner_id,d["hits"],d["domains"],d["num_copies"],d["upstream"],d["downstream"],str(novel)]) + "\n")

out.close()

cyto_out = open("novel_summary_cyto_new.csv","w")
cyto_out.write("toxin_node,size,domains,orientation,color,antitoxin_node,size,domains,orientation,color\n")

for h in hit_clusters:
	for p in hit_clusters[h]["partners"]:
		cyto_out.write(",".join( map(str, [h, 
			hit_clusters[h]["size"], 
			"/".join(list(hit_clusters[h]["domains"])) + "  (" + h + ")", 
			"toxin",
			"toxin",
			p,
			partners[p]["num_copies"],
			partners[p]["num_copies"],
			partners[p]["orientation"],
			partners[p]["novel"]])) + "\n")

cyto_out.close()
