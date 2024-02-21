from optimization import perform_variable_maximization
from model_class import Model
from toronto_set_up import optmdfpathway_base_milp_with_osmolarity, optmdfpathway_base_variables, biomass_reaction_id
import math

analysis_num_phi = "max_growth,phi_val,num_open,MDF\n"
MIN_MMDF = 0.01
optmdfpathway_base_variables["var_B"].bounds(MIN_MMDF, 1000)
optmdfpathway_base_variables["ATPM"].bounds(6.86,1000)
initial_phi = 0.17
current_phi = initial_phi
theta = -0.01
while current_phi >= 0.14560:
    print(f"phi is {current_phi}")
    optmdfpathway_base_variables["Phi_high_osm"].bounds(current_phi,1000)
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(current_phi, current_phi-theta)
    optmdfpathway_base_variables[biomass_reaction_id].bounds(0, 1000)
    growth_optimization_result = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, biomass_reaction_id)
    print(growth_optimization_result["values"][biomass_reaction_id])
    if growth_optimization_result["status"] == "Optimal":
        print("Feasible!")
        current_maxgrowth = growth_optimization_result["values"][biomass_reaction_id]
        optmdfpathway_base_variables[biomass_reaction_id].bounds(math.floor(growth_optimization_result["values"][biomass_reaction_id] * 10**4) / 10**4, 1000)
        optmdfpathway_result = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, "z_sum_var")
        num_open = optmdfpathway_result["values"]["z_sum_var"]
        analysis_num_phi += f"{current_maxgrowth},{current_phi},{num_open},{MIN_MMDF}\n"
        current_phi = round(current_phi - 0.00001,7)
    else:
        print("Infeasible!")
        current_phi = round(current_phi - 0.00001,7)
with open("./Results/Figure_2b.csv", "w") as f:
            f.write(analysis_num_k)
