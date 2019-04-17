
import os

files = os.listdir(".")

counts = {}

for f in files:
	if not f.endswith(".txt"):
		continue
	f = f.replace(".txt","")
	toks = f.split("_")
	for t in toks:
		if "H" in t:
			identifier = t
			if t not in counts:
				counts[t] = set()
				break
	for t in toks:
		if "P" in t:
			counts[identifier].add(t)

out = "c("
for c in counts:
	out += str(len(counts[c])) + ","

out = out[:-1]
out = out + ")"

print(out)