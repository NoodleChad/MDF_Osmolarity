from optimization import perform_variable_minimization
from optmdfpathway import get_thermodynamic_bottlenecks
from toronto_set_up import cobra_model, optmdfpathway_base_milp_with_osmolarity, optmdfpathway_base_variables, biomass_reaction_id

MIN_MMDF = 0
optmdfpathway_base_variables["var_B"].bounds(MIN_MMDF, 1e12)
growth_data = [
    (0.902, MIN_MMDF),
    (0.9008278, MIN_MMDF),
    (0.8966, MIN_MMDF),
    (0.87272, MIN_MMDF),
    (0.7432, MIN_MMDF)
]

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
        print(f"The minimal total concentration is {current_k}")
        print(f" *Perform bottleneck analysis")
        print(" Perform bottleneck analysis using OptMDFpathway base calculation")
        optmdfpathway_base_variables["osmolarity_sum_var"].bounds(0, current_k)
        bi, br = get_thermodynamic_bottlenecks(
        cobra_model=cobra_model,
        optmdfpathway_base_problem=optmdfpathway_base_milp_with_osmolarity,
        optmdfpathway_result=optmdfpathway_result,
        )

        with open(f"./Results/Bottleneck/bottleneck_analysis_{current_k}_newdGO.txt", "w", encoding="utf-8") as f:
            f.write(br)
