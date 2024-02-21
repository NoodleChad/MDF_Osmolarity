import math
from optimization import perform_variable_maximization, perform_variable_minimization
from toronto_set_up import optmdfpathway_base_milp_with_osmolarity, optmdfpathway_base_variables, biomass_reaction_id

min_mmdf = 0.01
mets_analysis = "max_growth,phi_val,mets,min,max\n"
optmdfpathway_base_variables["var_B"].bounds(min_mmdf, 100)
optmdfpathway_base_variables[biomass_reaction_id].bounds(0, 100)
optmdfpathway_base_variables["ATPM"].bounds(6.86,100)

current_phi = 0.17
theta = -0.01
while current_phi >= 0.140:
    print(" ~ Set phi:", current_phi, "~")
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(current_phi, current_phi-theta)
    optmdfpathway_base_variables[biomass_reaction_id].bounds(0, 1)
    growth_optimization_result = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, biomass_reaction_id)
    max_growth = -1*growth_optimization_result['objective_value']
    print(" ~ Set mu max:", max_growth)
    optmdfpathway_base_variables[biomass_reaction_id].bounds(max_growth, 1e2)
    print("Begin metabolite analysis")
    mets = "x_accoa_c"
    max  = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, mets)
    min  = perform_variable_minimization(optmdfpathway_base_milp_with_osmolarity, mets)
    mets_analysis += f"{max_growth},{current_phi},{mets},{math.exp(min['values'][mets])},{math.exp(max['values'][mets])}\n"
    mets = "x_ac_c"
    max  = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, mets)
    min  = perform_variable_minimization(optmdfpathway_base_milp_with_osmolarity, mets)
    mets_analysis += f"{max_growth},{current_phi},{mets},{math.exp(min['values'][mets])},{math.exp(max['values'][mets])}\n"
    mets = "x_atp_c"
    max  = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, mets)
    min  = perform_variable_minimization(optmdfpathway_base_milp_with_osmolarity, mets)
    mets_analysis += f"{max_growth},{current_phi},{mets},{math.exp(min['values'][mets])},{math.exp(max['values'][mets])}\n"
    mets = "x_nad_c"
    max  = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, mets)
    min  = perform_variable_minimization(optmdfpathway_base_milp_with_osmolarity, mets)
    mets_analysis += f"{max_growth},{current_phi},{mets},{math.exp(min['values'][mets])},{math.exp(max['values'][mets])}\n"
    current_phi -= 0.0005
print(" Store data")
with open("./Results/Figure_2d.csv", "w") as f:
    f.write(mets_analysis)