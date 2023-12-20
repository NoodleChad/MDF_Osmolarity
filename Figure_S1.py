from optimization import perform_variable_maximization
from model_class import Model
from toronto_set_up import optmdfpathway_base_milp_with_osmolarity, optmdfpathway_base_variables, biomass_reaction_id
import math

analysis_num_k = "k_val,mets,ln_met,exp_met\n"
MIN_MMDF = 1e-6
optmdfpathway_base_variables["var_B"].bounds(MIN_MMDF, 1000)
optmdfpathway_base_variables["ATPM"].bounds(0,1000)
optmdfpathway_base_variables["osmolarity_sum_var"].bounds(0, 0.1)
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
K_list = [0.0026,0.02,0.2,2]
for current_k in K_list:
    optmdfpathway_base_variables["potassium_pump_rate"].bounds(current_k,1000)
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(0, current_k)
    optmdfpathway_base_variables[biomass_reaction_id].bounds(0, 10)
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
                 analysis_num_k += f"{current_k},{ln_met[i]},{val_ln},{val_exp}\n"
        current_k = current_k - 0.005
    else:
        break
with open("./Results/FigureS1.csv", "w") as f:
            f.write(analysis_num_k)