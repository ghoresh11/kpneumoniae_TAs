

domain_partners = {}
domain_hits = {}
domain_known = {}

with open("novel_summary_cyto_new.csv") as f:
	for line in f:
		if line.startswith("toxin"):
			continue
		toks = line.strip().split(",")
		domains = toks[2].split()[0]
		domains = domains.split("/")

		hit_cluster = toks[0]
		partner_cluser = toks[5]

		for d in domains:
			if d not in domain_hits:
				domain_hits[d] = set()
				domain_partners[d] = set()
				domain_known[d] = 0
			domain_hits[d].add(hit_cluster)
			domain_partners[d].add(partner_cluser)
			if toks[-1] == "Known":
				domain_known[d] +=  1


with open("num_clusters_per_domain.csv","w") as out:
	out.write("Domain,HitClusters,PartnerClusters,Known\n")
	for d in domain_partners:
		out.write(",".join(map(str,[d, len(domain_hits[d]), len(domain_partners[d]), float(domain_known[d]) / len(domain_partners[d]) ])) + "\n")