import math
from helper import ensure_folder_existence
from optimization import perform_variable_minimization, perform_variable_maximization
from toronto_set_up import cobra_model, optmdfpathway_base_milp_with_osmolarity, optmdfpathway_base_variables, biomass_reaction_id

MIN_MMDF = 0.01
CVA_analysis = "max_growth,phi_val,mets,min,max\n"
print(f"==Performing analysis IV: Concentration variability analysis (CVA) @ OptMF>={MIN_MMDF} kJ/mol AND max µ at different phi===")

print(">Ensure existence of main analysis folder")
ensure_folder_existence("./toronto_results/analysis_iv")

print(f">Set OptMDF>={MIN_MMDF} kJ/mol")
optmdfpathway_base_variables["var_B"].bounds(MIN_MMDF, 1000)
optmdfpathway_base_variables["ATPM"].bounds(6.86,1000)
theta = -0.01

# Find the minimal value of phi that allows the growth rate
growth_data = [
    (0.87699, MIN_MMDF),
    (0.87075, MIN_MMDF),
    (0.84777, MIN_MMDF),
    (0.80446, MIN_MMDF),
    (0.71739, MIN_MMDF)
]

print(">Set preselected and partially rounded phis and growth rates from analysis II for this analysis")
phis_and_growth_rates = [] #first row should not be limiting here

for data_input in growth_data:
    print(" Load data...")
    current_maxgrowth = data_input[0]
    print(f" Setting minimal µ to {current_maxgrowth}")
    optmdfpathway_base_variables[biomass_reaction_id].bounds(current_maxgrowth, 1000)
    optmdfpathway_base_variables["var_B"].bounds(data_input[1], 1000)
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(0, 1000)
    print(" Perform OptMDFpathway...")
    optmdfpathway_result = perform_variable_minimization(optmdfpathway_base_milp_with_osmolarity, "osmolarity_sum_var")
    if optmdfpathway_result["status"] == "Optimal":
        sum_met = optmdfpathway_result['values']['osmolarity_sum_var']
        current_phi = sum_met + theta
        phen = (current_phi,current_maxgrowth)
        phis_and_growth_rates.append(phen)

print(">Set loop variables")

print(">Perform CVA analysis loop")
to_opt_val = []
l = (list(range(0,len(list((optmdfpathway_result['values']))))))
for i in l:
    reaction_id = list((optmdfpathway_result['values']))[i]
    if reaction_id.startswith("x_") and reaction_id.endswith("_c"):
        to_opt_val.append(reaction_id)
    elif reaction_id == 'phi_conc':
        to_opt_val.append(reaction_id)

for i in range(len(phis_and_growth_rates)):
    current_phi = phis_and_growth_rates[i][0]
    current_growth_rate = phis_and_growth_rates[i][1]
    print(" ~ mu set to ", current_growth_rate, "~")
    print(f" Setting maximal phi to {current_phi}")
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(current_phi, current_phi - theta)
    print(f">Set µ>={current_growth_rate} 1/h")
    optmdfpathway_base_variables[biomass_reaction_id].bounds(current_growth_rate, 1000)
    for j in range(len(to_opt_val)):
        min_val = perform_variable_minimization(optmdfpathway_base_milp_with_osmolarity, to_opt_val[j])
        max_val = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, to_opt_val[j])
        mets = to_opt_val[j].replace("x_","").replace("_c","")
        if to_opt_val == 'phi_conc':
            CVA_analysis += f"{current_growth_rate},{current_phi},{mets},{min_val['values'][to_opt_val[j]]},{max_val['values'][to_opt_val[j]]}\n"
        else:
            CVA_analysis += f"{current_growth_rate},{current_phi},{mets},{math.exp(min_val['values'][to_opt_val[j]])},{math.exp(max_val['values'][to_opt_val[j]])}\n"        

print(" Store data")
with open("./Results/Figure_2e.csv", "w") as f:
    f.write(mets_analysis)