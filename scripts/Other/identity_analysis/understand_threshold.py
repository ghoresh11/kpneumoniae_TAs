import os


class Cluster:

	def __init__(self, name, identity):
		self.name = name
		self.identity = str(identity)
		self.domains = {}
		self.phylogroups = {"kpi":0, "kpii":0, "kpiii":0}
		self.parent = "None"
		self.rows = []
		self.children = set()

	def __str__(self):
		return ("\n".join(map(str, [self.name, self.identity, self.domains, self.phylogroups, self.parent, self.children ])))

	def get_domains(self):
		res=""
		for d in self.domains:
			res += d+":" + str(self.domains[d]) + "/"
		return res[:-1]

	def get_phylogroups(self):
		res = ""
		for k in self.phylogroups:
			res += k + ":" + str(self.phylogroups[k]) + "/"
		return res[:-1]


def get_ident_clusters(ident):
	hits_dir = "/Users/gh11/klebsiella_TAs/identities/" + str(ident) + "_GROUP/hits_clusters"
	hit_files = os.listdir(hits_dir)

	clusters = []

	for hf in hit_files:
		if not hf.endswith(".txt"):
			continue
		name = hf.replace(".txt","")
		curr_cluster = Cluster(name, ident)
		with open(hits_dir + "/" + hf) as f:
			for line in f:
				if line.startswith("#") or line.startswith("Strain"):
					continue
				toks = line.strip().split(",")
				curr_cluster.phylogroups[phylogroups[toks[0]]] += 1
				curr_cluster.rows.append(",".join(toks[3:]))
				if toks[3] not in curr_cluster.domains:
					curr_cluster.domains[toks[3]] = 0
				curr_cluster.domains[toks[3]] += 1
		clusters.append(curr_cluster)
	return clusters


def connect_clusters(parents, children):
	for p1 in parents:
		for c1 in children:
			for r in p1.rows:
				if r in c1.rows:
					c1.parent = p1.name
					p1.children.add(c1.name)
	return


def create_network_output(clusters):
	nodes = open("nodes.csv", "w")
	edges = open("edges.csv", "w")
	r_output = open("r_output.csv","w")

	nodes.write("name,domains,phylogroups\n")
	edges.write("source,target\n")
	r_output.write("id, parent, identity, num_clusters, size\n")

	for c in clusters:
		for node in c:
			nodes.write(",".join([node.name+ "_"+node.identity, node.get_domains(), node.get_phylogroups()]) + "\n")
			parent_identity = str(int(node.identity) - 10)
			
			if node.parent != "None":
				edges.write(",".join([node.parent + "_" + parent_identity, node.name + "_" + node.identity]) + "\n")
			
			r_output.write(",".join(map(str, [node.name, node.parent, node.identity, len(node.children), sum(node.phylogroups.values()) ] )) + "\n")

			child_identity = str(int(node.identity) + 10)
			for c in node.children:
				edges.write(",".join([node.name + "_" + node.identity, c + "_" + child_identity]) + "\n")

	nodes.close()
	edges.close()
	r_output.close()
	return


### get phylogroup numbers
phylogroups={}
with open("/Users/gh11/klebsiella_TAs/metadata.csv") as f:
	for line in f:
		toks = line.strip().split(",")
		phylogroups[toks[0]] = toks[1]

clusters = []
for i in range(35,105,10):
	clusters.append(get_ident_clusters(i))

for i in range(0,len(clusters) - 1):
	connect_clusters(clusters[i],clusters[i+1])

create_network_output(clusters)
