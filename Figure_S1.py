from optimization import perform_variable_maximization
from model_class import Model
from toronto_set_up import optmdfpathway_base_milp_with_osmolarity, optmdfpathway_base_variables, biomass_reaction_id
import math

analysis_num_phi = "phi_val,mets,ln_met,exp_met\n"
MIN_MMDF = 0.01
optmdfpathway_base_variables["var_B"].bounds(MIN_MMDF, 1000)
optmdfpathway_base_variables["ATPM"].bounds(6.86,1000)
optmdfpathway_base_variables["osmolarity_sum_var"].bounds(0, 1000)
optmdfpathway_base_variables[biomass_reaction_id].bounds(0, 1000)
theta = -0.01
# Get the metabolite ID
optimization_result = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, biomass_reaction_id)
ln_met = []
exp_met = []
l = (list(range(0,len(list((optimization_result['values']))))))
for i in l:
    reaction_id = list((optimization_result['values']))[i]
    if reaction_id.startswith("x_") and reaction_id.endswith("_c"):
        ln_met.append(reaction_id)
        exp_met.append(reaction_id.replace('x_','expx_',1))

l = (list(range(0,len(ln_met))))
phi_list = [0.146,0.165,0.5,2]
for current_phi in phi_list:
    optmdfpathway_base_variables["Phi_high_osm"].bounds(current_phi,1000)
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(current_phi, current_phi - theta)
    optmdfpathway_base_variables[biomass_reaction_id].bounds(0, 1000)
    growth_optimization_result = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, biomass_reaction_id)
    if growth_optimization_result["status"] == "Optimal":
        current_maxgrowth = growth_optimization_result["values"][biomass_reaction_id]
        val_ln = []
        val_exp = []
        for i in l:
            if math.exp(growth_optimization_result['values'][ln_met[i]]) != 1:
                 val_ln = math.exp(growth_optimization_result['values'][ln_met[i]])
                 try:
                      val_exp = growth_optimization_result['values'][exp_met[i]]
                 except:
                      val_exp = 0
                 analysis_num_phi += f"{current_phi},{ln_met[i]},{val_ln},{val_exp}\n"
    else:
        break
with open("./Results/FigureS1.csv", "w") as f:
            f.write(analysis_num_k)