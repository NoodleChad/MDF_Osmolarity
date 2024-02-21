import math
from helper import ensure_folder_existence
from optimization import perform_variable_minimization, perform_variable_maximization
from model_class import Model
from toronto_set_up import cobra_model, optmdfpathway_base_milp_with_osmolarity, optmdfpathway_base_variables, biomass_reaction_id

MIN_MMDF = 0.01
Z_analysis = "max_growth,phi_val,rxns,Z\n"
print(f"==Performing analysis IV: Concentration variability analysis (CVA) @ OptMF>={MIN_MMDF} kJ/mol AND max µ at different phi===")

print(">Ensure existence of main analysis folder")
ensure_folder_existence("./toronto_results/analysis_iv/Z_analysis")

optmdfpathway_base_variables["ATPM"].bounds(6.86,1000)
optmdfpathway_base_variables["var_B"].bounds(MIN_MMDF, 100)
theta = -0.01
growth_data = [
    (0.87699, MIN_MMDF),
    (0.87075, MIN_MMDF),
    (0.84777, MIN_MMDF),
    (0.80446, MIN_MMDF),
    (0.71739, MIN_MMDF)
]

print(">Set preselected and partially rounded growth rates for this analysis")
phis_and_growth_rates = [] 

# Here we calculate the minimal phi associated to these growth rates
for data_input in growth_data:
    print(" Load data...")
    current_maxgrowth = data_input[0]
    print(f" Setting minimal µ to {current_maxgrowth}")
    optmdfpathway_base_variables[biomass_reaction_id].bounds(current_maxgrowth, 1e2)
    optmdfpathway_base_variables["var_B"].bounds(data_input[1], 100)
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(0, 100)
    print(" Perform OptMDFpathway...")
    optmdfpathway_result = perform_variable_minimization(optmdfpathway_base_milp_with_osmolarity, "osmolarity_sum_var")
    if optmdfpathway_result["status"] == "Optimal":
        sum_met = optmdfpathway_result['values']['osmolarity_sum_var']
        if current_maxgrowth == 0.87699:
            current_phi = (sum_met + theta)*1.2
            phen = (current_phi,current_maxgrowth)
            phis_and_growth_rates.append(phen)
        else:
            current_phi = (sum_met + theta)
            phen = (current_phi,current_maxgrowth)
            phis_and_growth_rates.append(phen)

# Here we find the location of each binairy variables and find what continuous variables are associated to them
to_opt_val = []
to_opt_rxns = []
l = (list(range(0,len(list((optmdfpathway_result['values']))))))
for i in l:
    reaction_id = list((optmdfpathway_result['values']))[i]
    if reaction_id.startswith("z_var_"):
        to_opt_val.append(reaction_id)
        to_opt_rxns.append(reaction_id.replace('z_var_',''))

# Here we will test wether the variable is thermodynimically feasible by seing if each invididual binary variable
# can be 1 (will be considered in the MDF problem if it is). Subsequently, to prevent mathematical instability leading
# to false postive, we also test that the continuous variable is also feasible.
for i in range(len(phis_and_growth_rates)):
    current_phi = phis_and_growth_rates[i][0]
    current_growth_rate = phis_and_growth_rates[i][1]
    print(" ~ mu set to ", current_growth_rate, "~")
    print(f" Setting maximal phi to {current_phi}")
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(current_phi, current_phi-theta)
    print(f">Set µ>={current_growth_rate} 1/h")
    optmdfpathway_base_variables[biomass_reaction_id].bounds(current_growth_rate, 100)
    for j in range(len(to_opt_val)):
        max_val = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, to_opt_val[j])
        val = max_val['values'][to_opt_val[j]]
        if max_val['status']=='Optimal' and max_val['values'][to_opt_val[j]] == 1.0:
            flux = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, to_opt_rxns[j])
            if flux['status']=='Optimal' and flux['values'][to_opt_rxns[j]] > 1e-10:
                val = '1.0'
            else:
                val = '0'
        mets = to_opt_val[j].replace("z_var_","")
        Z_analysis += f"{current_growth_rate},{current_phi},{mets},{val}\n"

# using this files we can tell what reactions are feasible or not and make the figure in escher by running convertCSVtoJSON.py,
# Pair_fwd.py and 4map.py
print(" Store data")
with open("./Results/Figure_3.csv", "w") as f:
    f.write(Z_analysis)