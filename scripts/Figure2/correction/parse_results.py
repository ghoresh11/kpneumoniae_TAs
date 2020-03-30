

## 1. constuct dictionary of
## interpro_scan -> gene name
## GO term -> gene names
## in the end make a table of gene_name -> GO terms (biological process, namespaces...) + interpro_scan (family, maybe domain, Homologous_superfamily)

class Entry:

    def __init__(self, name):
        self.name = name
        self.classification = "other general"
        self.pfam_signature = set()
        self.biological_process = set()
        self.molecular_function = set()
        self.cellular_component = set()
        self.family = set()
        self.domain = set()
        self.superfamily = set()
        return

    def __str__(self):
        return ",".join([self.name, self.classification, "/".join(list(self.pfam_signature)),
        "/".join(list(self.biological_process)), "/".join(list(self.molecular_function)),
        "/".join(list(self.cellular_component)), "/".join(list(self.superfamily)),
        "/".join(list(self.family)), "/".join(list(self.domain))]
        )


interpro_scan_ids = {}
go_terms = {}
genes = {}



with open("antitoxins.gff") as f:
    for line in f:
        if line.startswith("#"):
            continue
        if line.startswith(">"):
            break
        toks = line.strip().split("\t")
        entry_name = toks[0].split("_")[0]

        if entry_name not in genes:
            genes[entry_name] = Entry(entry_name)
        line = line.replace(",",";")
        if "InterPro" in line:
            ipr = line.split("IPR")[-1]
            ipr = ipr.split("\"")[0]
            ipr = "IPR" + ipr
            if ipr not in interpro_scan_ids:
                interpro_scan_ids[ipr] = set()
            interpro_scan_ids[ipr].add(entry_name)
            ## get the interpro id and save in interpro scan ids
        if "signature_desc" in line: ## pfam
            ## get pfam desciptor and add to entry
            pfam = line.split("signature_desc=")[-1]
            pfam = pfam.split(";")[0]
            genes[entry_name].pfam_signature.add(pfam)
        if "Ontology_term" in line:
            gos = line.strip().split("Ontology_term=")[-1]
            gos = gos.split(";")[0]
            gos = gos.replace("\"","")
            gos = gos.split(";")
            for go in gos:
                if go not in go_terms:
                    go_terms[go] = set()
                go_terms[go].add(entry_name)



## step 2: go over the interpro scan file to extract the family etc.
with open("050819_entry.list") as f:
    for line in f:
        toks = line.strip().split("\t")
        if toks[1] not in ["Homologous_superfamily","Family","Domain"]:
            continue
        if toks[0] not in interpro_scan_ids:
            continue
        for gene in interpro_scan_ids[toks[0]]:
            if toks[1] == "Homologous_superfamily":
                genes[gene].superfamily.add(toks[2])
            elif toks[1] == "Family":
                genes[gene].family.add(toks[2])
            else:
                genes[gene].domain.add(toks[2])


## step 3: go over the GO term file and extact the namespaces
flag = False
with open("/lustre/scratch118/infgen/pathogen/pathpipe/go_ontology/gene_ontology_20180201.obo") as f:
    for line in f:
        if line.startswith("id:"):
            curr_term = line.strip().split()[-1]
            if curr_term not in go_terms:
                flag = False
            else:
                flag = True
            continue
        if flag and line.startswith("name:"):
            curr_name = line.strip().replace("name: ","")
            continue
        if flag and line.startswith("namespace"):
            curr_namespace = line.strip().replace("namespace: ","")
            if curr_namespace == "biological_process":
                for gene in go_terms[curr_term]:
                    genes[gene].biological_process.add(curr_name)
            elif curr_namespace == "molecular_function":
                for gene in go_terms[curr_term]:
                    genes[gene].molecular_function.add(curr_name)
            else:
                for gene in go_terms[curr_term]:
                    genes[gene].cellular_component.add(curr_name)
            flag = False

defined_processes = ["transmembrane transport", "DNA recombination", "DNA replication initiation","DNA replication","methylation", "pathogenesis","cell adhesion"]
## classify genes
for gene in genes:
    gene_line = str(genes[gene]).lower()
    gene = genes[gene]
    if "phage" in gene_line or "viral" in gene_line or "virus" in gene_line or "capsid" in gene_line:
        gene.classification = "phage"
    elif "plasmid" in gene_line:
        gene.classification = "plasmid"
    elif "conjuga" in gene_line:
        gene.classification = "conjugation"
    elif "mobilisation" in gene_line or "mobilization" in gene_line:
        gene.classification = "mobilisation"
    elif "toxin-antitoxin" in gene_line or "toxin antitoxin" in gene_line:
        gene.classification = "toxin-antitoxin system"
    elif "abc transporter" in gene_line:
        gene.classification = "ABC transporter"
    elif ("secretion system" in gene_line or "secretion-system" in line):
        gene_line = gene_line.replace("-"," ")
        system_type = gene_line.split(" secretion system")[0]
        system_type = system_type.split("type ")[-1]
        gene.classification = "type " + system_type + " secretion system"
    elif gene.biological_process ==  "cell adhesion" and "fimbria" in gene_line:
        gene.biological = "fimbrial, cell adhesion"
    elif "transposase" in str(gene.pfam_signature).lower():
        gene.classification = "transposase"
    elif "resolvase" in str(gene.pfam_signature).lower():
        gene.classification = "resolvase"
    elif "membrane" in gene_line:
        gene.classification = "membrane"
    elif "crispr" in gene_line:
        gene.classification = "CRISPR"
    elif "toxin" in gene_line:
        gene.classification = "toxin"
    elif "pilus" in gene_line:
        gene.classification = "pilus"
    elif len(gene.biological_process & set(defined_processes)) > 0 :
        gene.classification = list(gene.biological_process & set(defined_processes))[0]
    elif len(gene.biological_process) > 0:
        gene.classification = "/".join(list(gene.biological_process))
    elif len(gene.molecular_function) > 0:
        gene.classification = "/".join(list(gene.molecular_function))
    elif "duf" in gene_line or "domain of unknown function" in gene_line:
        gene.classification = "DUF"


with open("summary_all.csv","w") as out:
    out.write("name,classification,pfam_signature,biological_process,molecular_function,cellular_component,superfamily,family,domain\n")
    for gene in genes:
        out.write(str(genes[gene]) + "\n")
