import os


for file in os.listdir("/Users/gh11/klebsiella_TAs/GROUP/partners_clusters/"):
	hit_id = file.split("/")[-1]
	hit_id = hit_id.replace(".txt","")

	out = open(hit_id + ".fasta","w")
	
	counter = 0
	if not file.endswith(".txt"):
		continue 
	with open("/Users/gh11/klebsiella_TAs/GROUP/partners_clusters/" + file) as f:
		for line in f:
			if line.startswith("#") or line.startswith("Strain"):
				continue
			toks = line.strip().split(",")
			out.write(">" + str(counter) + "|" + toks[1] + "\n" + toks[14] + "\n")
			counter += 1
	out.close()
