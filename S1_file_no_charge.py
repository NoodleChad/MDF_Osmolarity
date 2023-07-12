import math
from helper import ensure_folder_existence
from optimization import perform_variable_minimization, perform_variable_maximization
from toronto_set_up import cobra_model, optmdfpathway_base_milp_with_osmolarity, optmdfpathway_base_variables, biomass_reaction_id

MIN_MMDF = 1e-6
CVA_analysis = "max_growth,k_val,mets,min,max\n"
print(f"==Performing analysis IV: Concentration variability analysis (CVA) @ OptMF>={MIN_MMDF} kJ/mol AND max µ at different K===")

print(">Ensure existence of main analysis folder")
ensure_folder_existence("./toronto_results/analysis_iv")

print(f">Set OptMDF>={MIN_MMDF} kJ/mol")
optmdfpathway_base_variables["var_B"].bounds(MIN_MMDF, 1e12)

# Find the minimal value of K that allows the growth rate
growth_data = [
    (0.902, MIN_MMDF),
    (0.9008278, MIN_MMDF),
    (0.8966, MIN_MMDF),
    (0.87272, MIN_MMDF),
    (0.7432, MIN_MMDF)
]

print(">Set preselected and partially rounded Ks and growth rates from analysis II for this analysis")
ks_and_growth_rates = [] #first row should not be limiting here

for data_input in growth_data:
    print(" Load data...")
    current_maxgrowth = data_input[0]
    print(f" Setting minimal µ to {current_maxgrowth}")
    optmdfpathway_base_variables[biomass_reaction_id].bounds(current_maxgrowth, 1e12)
    optmdfpathway_base_variables["var_B"].bounds(data_input[1], 1e12)
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(0, 0.8)
    print(" Perform OptMDFpathway...")
    optmdfpathway_result = perform_variable_minimization(optmdfpathway_base_milp_with_osmolarity, "osmolarity_sum_var")
    if optmdfpathway_result["status"] == "Optimal":
        current_k = optmdfpathway_result['values']['osmolarity_sum_var']
        phen = (current_k,current_maxgrowth)
        ks_and_growth_rates.append(phen)

print(">Set loop variables")

print(">Perform CVA analysis loop")
to_opt_val = []
l = (list(range(0,len(list((optmdfpathway_result['values']))))))
for i in l:
    reaction_id = list((optmdfpathway_result['values']))[i]
    if reaction_id.startswith("x_") and reaction_id.endswith("_c"):
        to_opt_val.append(reaction_id)
    elif reaction_id == 'k_conc':
        to_opt_val.append(reaction_id)

for i in range(len(ks_and_growth_rates)):
    current_k = ks_and_growth_rates[i][0]
    current_growth_rate = ks_and_growth_rates[i][1]
    print(" ~ mu set to ", current_growth_rate, "~")
    print(f" Setting maximal K to {current_k}")
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(0, current_k)
    print(f">Set µ>={current_growth_rate} 1/h")
    optmdfpathway_base_variables[biomass_reaction_id].bounds(current_growth_rate, 1e12)
    for j in range(len(to_opt_val)):
        min_val = perform_variable_minimization(optmdfpathway_base_milp_with_osmolarity, to_opt_val[j])
        max_val = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, to_opt_val[j])
        mets = to_opt_val[j].replace("x_","").replace("_c","")
        if to_opt_val == 'k_conc':
            CVA_analysis += f"{current_growth_rate},{current_k},{mets},{min_val['values'][to_opt_val[j]]},{max_val['values'][to_opt_val[j]]}\n"
        else:
            CVA_analysis += f"{current_growth_rate},{current_k},{mets},{math.exp(min_val['values'][to_opt_val[j]])},{math.exp(max_val['values'][to_opt_val[j]])}\n"

print(" Store data")
with open("./toronto_results/analysis_iv/CVA_new_dG0_nomal_max_charge.csv", "w") as f:
    f.write(CVA_analysis)