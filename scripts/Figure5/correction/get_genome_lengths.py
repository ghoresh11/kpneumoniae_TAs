import os

print("Genome,Contig,Length")
files = os.listdir(".")

for f in files:
	if not f.endswith(".gff"):
		continue
	genome_name = f.replace(".gff","")
	with open(f) as f_open:
		for line in f_open:
			if line.startswith("##sequence-region"):
				toks = line.strip().split()
				print(genome_name + "," + toks[1] + "," + toks[-1])
				continue
			if not line.startswith("#"):
				break