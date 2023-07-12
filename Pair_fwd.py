import json

with open("./toronto_results/analysis_iv/Z_analysis/z_analysis_74.json", "r") as f:
    data = json.load(f)

l = (list(range(0,len(list(data)))))
new_dict = {}
for i in l:
    if list(data)[i].endswith('_FWD'):
        if list(data)[i].startswith('EX_'):
            continue
        else:
            rxns = list(data)[i].replace("_FWD","")
            new_dict[rxns]=data[list(data)[i]]
            if data[list(data)[i]]['Z']=='1.0' and data[list(data)[i+1]]['Z']=='1.0':
                new_dict[rxns]['Z'] = '3'
            elif data[list(data)[i]]['Z']=='1.0':
                new_dict[rxns]['Z'] = '1'
            elif data[list(data)[i+1]]['Z']=='1.0':
                new_dict[rxns]['Z'] = '-2'
            else:
                new_dict[rxns]['Z'] = '0'
    elif list(data)[i].endswith('_REV'):
        continue
    elif list(data)[i] == "var_B":
        continue
    elif list(data)[i].startswith('EX_'):
        continue
    else:
        new_dict[list(data)[i]]=data[list(data)[i]]

with open("./toronto_results/analysis_iv/Z_analysis/z_analysis_74.json", "w") as f:
    json.dump(new_dict,f)