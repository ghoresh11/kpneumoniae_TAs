import os

d = {}

with open("/Users/gh11/klebsiella_TAs/complete_categories.csv") as f:
	for line in f:
		toks = line.strip().split(",")
		if toks[0] == "Toxin_Cluster":
			continue
		if toks[0] not in d:
			d[toks[0]] = {}
			d[toks[0]]["vals"] = []
			d[toks[0]]["operons"] = []
			d[toks[0]]["structs"] = 0
			d[toks[0]]["total"] = 0
		d[toks[0]]["vals"].append(int(toks[3]))
		d[toks[0]]["structs"] += 1
		d[toks[0]]["total"] += int(toks[3])

		d[toks[0]]["operons"].append(toks[1].replace("-","X") + "_" + toks[0] + "_" + toks[2].replace("-","X"))


with open("/Users/gh11/klebsiella_TAs/toxins_categories.csv") as f:
	for line in f:
		toks = line.strip().split(",")
		if toks[0] == "Toxin_Cluster":
			continue
		d[toks[0]]["class"] = toks[4]


files = os.listdir("/Users/gh11/klebsiella_TAs/GROUP/hits_clusters")
for f in files:
	if not f.endswith(".txt"):
		continue
	key = f.replace(".txt","")
	with open("/Users/gh11/klebsiella_TAs/GROUP/hits_clusters/" + f) as f_open:
		for line in f_open:
			if "Domains" in line:
				toks = line.strip().split()
				toks[-1]=toks[-1].replace(",","/")
				d[key]["domain"] = toks[-1]
				break

out = open("structures_summary.csv","w")
out.write("ID,domain,total,num_structs,color,class,operons,count,vals\n")
for toxin in d:
	prefix = ",".join(map(str,[toxin, d[toxin]["domain"], d[toxin]["total"], d[toxin]["structs"], "#d3d3d3", d[toxin]["class"]]))
	for i in range(0,len(d[toxin]["vals"])):
		v = d[toxin]["vals"][i]
		operon = d[toxin]["operons"][i]
		value = float(v) / float(d[toxin]["total"])
		out.write(prefix + "," + operon + ","+ str(v) + "," + str(value) + "\n")


out.close()
