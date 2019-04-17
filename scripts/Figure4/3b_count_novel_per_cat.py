import os

## go over the toxin categories and count how many novel and how many known for each of the five categories

## 1. classify each antitoxin
antitoxin_clusters = {}

blast_files = os.listdir("/Users/gh11/klebsiella_TAs/novel_antitoxins")
known = 0
new = 0

for f in blast_files:
	if not f.endswith("_result"):
		continue
	curr_antitoxin = f.split(".")[0]
	antitoxin_clusters[curr_antitoxin] = "new"
	new += 1 
	with open(f) as f_open:
		for line in f_open:
			toks = line.strip().split()
			if float(toks[2]) > 75:
				antitoxin_clusters[curr_antitoxin] = "known"
				new -= 1
				known += 1
				break


operon_files = os.listdir("/Users/gh11/klebsiella_TAs/GROUP/complete_clusters/")
toxins = {}
for f in operon_files:
	if not f.endswith(".txt"):
		continue
	f = f.replace(".txt","")
	toks = f.split("_")
	if len(toks) == 3:
		toxin = toks[1]
		antitoxins = [toks[0], toks[2]]
	elif "H" in toks[0]:
		toxin = toks[0]
		antitoxins = [toks[1]]
	else:
		toxin = toks[1]
		antitoxins = [toks[0]]

	if toxin not in toxins:
		toxins[toxin] = {"new":0, "known":0}

	for a in antitoxins:
		toxins[toxin][antitoxin_clusters[a]] += 1



with open("/Users/gh11/klebsiella_TAs/toxins_categories.csv") as f:
	for line in f:
		if line.startswith("Toxin"):
			continue
		toks = line.strip().split(",")
		toxin = toks[0]
		toxins[toxin]["class"] = toks[4]

out = open("known_novel_per_toxin.csv","w")
out.write("Toxin_ID,Class,Known,Novel\n")
for t in toxins:
	out.write(",".join(map(str, [t, toxins[t]["class"],  toxins[t]["known"],  toxins[t]["new"]])) + "\n")
out.close()




