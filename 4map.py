import json

with open("./toronto_results/analysis_iv/Z_analysis/z_analysis_74.json") as f:
    data = json.load(f)
with open("./toronto_results/analysis_iii/iJO1366.json", "r") as f:
    core = json.load(f)
with open("./toronto_results/analysis_iii/iML1515.json", "r") as f:
    genome = json.load(f)

l = list(range(0,len(list(core["reactions"]))))
name_core = []
rxns_core = []
stoi_core = []
for i in l:
    name_core.append(core["reactions"][i]["name"])
    rxns_core.append(core["reactions"][i]["id"])
    stoi_core.append(core["reactions"][i]["metabolites"])

l = list(range(0,len(list(genome["reactions"]))))
name_genome = []
rxns_genome = []
stoi_genome = []
for i in l:
    name_genome.append(genome["reactions"][i]["name"])
    rxns_genome.append(genome["reactions"][i]["id"])
    stoi_genome.append(genome["reactions"][i]["metabolites"])


l = (list(range(0,len(list(data)))))
open_react = {}
for i in l:
    open_react[list(data)[i]] = data[list(data)[i]]['Z']

l = list(range(0,len(rxns_genome)))
for i in l:
    if list(rxns_genome)[i] in list(open_react):
        continue
    else:
        open_react[list(rxns_genome)[i]] = 1

l = list(range(0,len(list(open_react))))
id = []
value = []
for i in l:
    if stoi_genome[rxns_genome.index(list(open_react)[i])] in list(stoi_core):
        id.append(rxns_core[stoi_core.index(stoi_genome[rxns_genome.index(list(open_react)[i])])])
        value.append(open_react[list(open_react)[i]])
    else:
        id.append(list(open_react)[i])
        value.append(open_react[list(open_react)[i]])

# replace biomass name
id[id.index("BIOMASS_Ec_iML1515_core_75p37M")] = "BIOMASS_Ec_iJO1366_core_53p95M"

    # pair all PFK together
    # if "PFK" not in list(id):
    #     id.append("PFK")
    #     value.append(0)
    # if "PFK_2" in list(id):
    #     value[id.index("PFK")] = value[id.index("PFK")] + value[id.index("PFK_2")]
    #     value[id.index("PFK_2")] = 0
    # if "PFK_3" in list(id):
    #     value[id.index("PFK")] = value[id.index("PFK")] + value[id.index("PFK_3")]
    #     value[id.index("PFK_3")] = 0
    # make the actual dictionary

Flux = {}
l = list(range(0,len(list(id))))
for i in l:
    Flux[id[i]] = value[i]

with open("./toronto_results/analysis_iv/Z_analysis/z_analysis_74_figure.json", "w") as f:
    json.dump(Flux,f)