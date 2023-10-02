from optimization import perform_variable_minimization, perform_variable_maximization
from toronto_set_up import optmdfpathway_base_milp_with_osmolarity, optmdfpathway_base_variables, biomass_reaction_id

analysis_num_k = "max_growth,k_val,num_open,MDF\n"
analysis_data = "rxns,z_val\n"
MIN_MMDF = 1e-6
optmdfpathway_base_variables["var_B"].bounds(MIN_MMDF, 1e6)
initial_k = 0.03
current_k = initial_k
while current_k >= 0.001:
    print(f"K is {current_k}")
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(0, current_k)
    optmdfpathway_base_variables[biomass_reaction_id].bounds(0, 1e6)
    growth_optimization_result = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, biomass_reaction_id)
    print(growth_optimization_result["values"][biomass_reaction_id])
    if growth_optimization_result["status"] == "Optimal":
        print("Feasible!")
        current_maxgrowth = growth_optimization_result["values"][biomass_reaction_id]
        optmdfpathway_result = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, "z_sum_var")
        num_open = optmdfpathway_result["values"]["z_sum_var"]
        analysis_num_k += f"{current_maxgrowth},{current_k},{num_open},{MIN_MMDF}\n"
        

        current_k = current_k - 0.00005
    else:
        print("Infeasible!")
        break
with open("./Results/Figure_2a.csv", "w") as f:
            f.write(analysis_num_k)
    
