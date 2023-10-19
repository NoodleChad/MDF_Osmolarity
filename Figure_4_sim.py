from optimization import perform_variable_maximization
from model_class import Model
from toronto_set_up import optmdfpathway_base_milp_with_osmolarity, optmdfpathway_base_variables, biomass_reaction_id

analysis_num_k = "max_growth,k_val,num_open,MDF\n"
MIN_MMDF = 5.25
optmdfpathway_base_variables["var_B"].bounds(MIN_MMDF, 1e6)
optmdfpathway_base_variables["ATPM"].bounds(80,1000)
initial_k = 0.873
current_k = initial_k
while current_k >= 0.01:
    print(f"K is {current_k}")
    optmdfpathway_base_variables["potassium_pump_rate"].bounds(current_k,1000)
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(0, current_k)
    optmdfpathway_base_variables[biomass_reaction_id].bounds(0, 1e6)
    growth_optimization_result = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, biomass_reaction_id)
    print(growth_optimization_result["values"][biomass_reaction_id])
    if growth_optimization_result["status"] == "Optimal":
        print("Feasible!")
        current_maxgrowth = growth_optimization_result["values"][biomass_reaction_id]
        num_open = 1
        analysis_num_k += f"{current_maxgrowth},{current_k},{num_open},{MIN_MMDF}\n"
        current_k = current_k - 0.005
    else:
        print("Infeasible!")
        break
with open("./Experimental/KCl_simulated.csv", "w") as f:
            f.write(analysis_num_k)

analysis_num_k = "max_growth,k_val,num_open,MDF\n"
MIN_MMDF = 5.25
optmdfpathway_base_variables["var_B"].bounds(MIN_MMDF, 1e6)
optmdfpathway_base_variables["ATPM"].bounds(62.5,1000)
initial_k = 0.873
current_k = initial_k
while current_k >= 0.01:
    print(f"K is {current_k}")
    optmdfpathway_base_variables["potassium_pump_rate"].bounds(current_k,1000)
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(0, current_k)
    optmdfpathway_base_variables[biomass_reaction_id].bounds(0, 1e6)
    growth_optimization_result = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, biomass_reaction_id)
    print(growth_optimization_result["values"][biomass_reaction_id])
    if growth_optimization_result["status"] == "Optimal":
        print("Feasible!")
        current_maxgrowth = growth_optimization_result["values"][biomass_reaction_id]
        num_open = 1
        analysis_num_k += f"{current_maxgrowth},{current_k},{num_open},{MIN_MMDF}\n"
        current_k = current_k - 0.005
    else:
        print("Infeasible!")
        break
with open("../Experimental/NaCl_simulated.csv", "w") as f:
            f.write(analysis_num_k)
