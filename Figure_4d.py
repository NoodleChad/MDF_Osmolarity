from optimization import perform_variable_maximization
from model_class import Model
from toronto_set_up import optmdfpathway_base_milp_with_osmolarity, optmdfpathway_base_variables, biomass_reaction_id


## Set-up
analysis_num_phi = "max_growth,phi_val,MDF\n"
MIN_MMDF = 0.01
optmdfpathway_base_variables["var_B"].bounds(MIN_MMDF, 1000)
optmdfpathway_base_variables["ATPM"].bounds(6.86,1000)

## WT
initial_phi = 0.1457
current_phi = initial_phi
theta = -0.01
while current_phi <= 0.14690:
    print(f"phi is {current_phi}")
    optmdfpathway_base_variables["Phi_high_osm"].bounds(current_phi,1000)
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(current_phi, current_phi-theta)
    optmdfpathway_base_variables[biomass_reaction_id].bounds(0, 1000)
    growth_optimization_result = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, biomass_reaction_id)
    print(growth_optimization_result["values"][biomass_reaction_id])
    if growth_optimization_result["status"] == "Optimal":
        print("Feasible!")
        current_maxgrowth = growth_optimization_result["values"][biomass_reaction_id]
        analysis_num_phi += f"{current_maxgrowth},{current_phi},{MIN_MMDF}\n"
        current_phi = round(current_phi + 0.000005,8)
    else:
        print("Infeasible!")
        current_phi = round(current_phi + 0.000005,8)
with open("./Results/Figure_4d_WT.csv", "w") as f:
    f.write(Z_analysis)
    
## Delta EDD
analysis_num_phi = "max_growth,phi_val,MDF\n"
optmdfpathway_base_variables["EDD"].bounds(0,0)
initial_phi = 0.1457
current_phi = initial_phi
theta = -0.01
while current_phi <= 0.14690:
    print(f"phi is {current_phi}")
    optmdfpathway_base_variables["Phi_high_osm"].bounds(current_phi,1000)
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(current_phi, current_phi-theta)
    optmdfpathway_base_variables[biomass_reaction_id].bounds(0, 1000)
    growth_optimization_result = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, biomass_reaction_id)
    print(growth_optimization_result["values"][biomass_reaction_id])
    if growth_optimization_result["status"] == "Optimal":
        print("Feasible!")
        current_maxgrowth = growth_optimization_result["values"][biomass_reaction_id]
        analysis_num_phi += f"{current_maxgrowth},{current_phi},{MIN_MMDF}\n"
        current_phi = round(current_phi + 0.000005,8)
    else:
        print("Infeasible!")
        current_phi = round(current_phi + 0.000005,8)
with open("./Results/Figure_4d_EDD.csv", "w") as f:
    f.write(Z_analysis)


## Delta PGL
analysis_num_phi = "max_growth,phi_val,MDF\n"
optmdfpathway_base_variables["EDD"].bounds(0,1000)
optmdfpathway_base_variables["PGL"].bounds(0,0)
initial_phi = 0.1457
current_phi = initial_phi
theta = -0.01
while current_phi <= 0.14690:
    print(f"phi is {current_phi}")
    optmdfpathway_base_variables["Phi_high_osm"].bounds(current_phi,1000)
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(current_phi, current_phi-theta)
    optmdfpathway_base_variables[biomass_reaction_id].bounds(0, 1000)
    growth_optimization_result = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, biomass_reaction_id)
    print(growth_optimization_result["values"][biomass_reaction_id])
    if growth_optimization_result["status"] == "Optimal":
        print("Feasible!")
        current_maxgrowth = growth_optimization_result["values"][biomass_reaction_id]
        analysis_num_phi += f"{current_maxgrowth},{current_phi},{MIN_MMDF}\n"
        current_phi = round(current_phi + 0.000005,8)
    else:
        print("Infeasible!")
        current_phi = round(current_phi + 0.000005,8)
with open("./Results/Figure_4d_PGL.csv", "w") as f:
    f.write(Z_analysis)

## Delta GND
analysis_num_phi = "max_growth,phi_val,MDF\n"
optmdfpathway_base_variables["EDD"].bounds(0,1000)
optmdfpathway_base_variables["PGL"].bounds(0,1000)
optmdfpathway_base_variables["GND"].bounds(0,0)
initial_phi = 0.1457
current_phi = initial_phi
theta = -0.01
while current_phi <= 0.14690:
    print(f"phi is {current_phi}")
    optmdfpathway_base_variables["Phi_high_osm"].bounds(current_phi,1000)
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(current_phi, current_phi-theta)
    optmdfpathway_base_variables[biomass_reaction_id].bounds(0, 1000)
    growth_optimization_result = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, biomass_reaction_id)
    print(growth_optimization_result["values"][biomass_reaction_id])
    if growth_optimization_result["status"] == "Optimal":
        print("Feasible!")
        current_maxgrowth = growth_optimization_result["values"][biomass_reaction_id]
        analysis_num_phi += f"{current_maxgrowth},{current_phi},{MIN_MMDF}\n"
        current_phi = round(current_phi + 0.000005,8)
    else:
        print("Infeasible!")
        current_phi = round(current_phi + 0.000005,8)
with open("./Results/Figure_4d_GND.csv", "w") as f:
    f.write(Z_analysis)