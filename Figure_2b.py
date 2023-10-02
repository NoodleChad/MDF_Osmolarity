import math
from helper import ensure_folder_existence, json_write, save_xy_point_plot
from optimization import perform_variable_maximization, perform_variable_minimization
from optmdfpathway import get_thermodynamic_bottlenecks
from toronto_set_up import cobra_model, optmdfpathway_base_milp_with_osmolarity, optmdfpathway_base_variables, biomass_reaction_id, used_max_growth

min_mmdf = 1e-6
mets_analysis = "max_growth,k_val,mets,min,max\n"
optmdfpathway_base_variables["var_B"].bounds(min_mmdf, 1e6)
optmdfpathway_base_variables[biomass_reaction_id].bounds(0, 1)

current_k = 0.003
while current_k >= 0.002:
    print(" ~ Set K:", current_k, "~")
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(0, current_k)
    optmdfpathway_base_variables[biomass_reaction_id].bounds(0, 1)
    growth_optimization_result = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, biomass_reaction_id)
    max_growth = round(-1*growth_optimization_result['objective_value'],4)
    print(" ~ Set mu max:", max_growth)
    optmdfpathway_base_variables[biomass_reaction_id].bounds(max_growth, 1e2)
    print("Begin metabolite analysis")
    mets = "x_accoa_c"
    max  = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, mets)
    min  = perform_variable_minimization(optmdfpathway_base_milp_with_osmolarity, mets)
    mets_analysis += f"{max_growth},{current_k},{mets},{math.exp(min['values'][mets])},{math.exp(max['values'][mets])}\n"
    mets = "x_ac_c"
    max  = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, mets)
    min  = perform_variable_minimization(optmdfpathway_base_milp_with_osmolarity, mets)
    mets_analysis += f"{max_growth},{current_k},{mets},{math.exp(min['values'][mets])},{math.exp(max['values'][mets])}\n"
    mets = "x_atp_c"
    max  = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, mets)
    min  = perform_variable_minimization(optmdfpathway_base_milp_with_osmolarity, mets)
    mets_analysis += f"{max_growth},{current_k},{mets},{math.exp(min['values'][mets])},{math.exp(max['values'][mets])}\n"
    mets = "x_nad_c"
    max  = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, mets)
    min  = perform_variable_minimization(optmdfpathway_base_milp_with_osmolarity, mets)
    mets_analysis += f"{max_growth},{current_k},{mets},{math.exp(min['values'][mets])},{math.exp(max['values'][mets])}\n"
    current_k -= 0.0001
print(" Store data")
with open("./Results/Figure_1b.csv", "w") as f:
    f.write(mets_analysis)
