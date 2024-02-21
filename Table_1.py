from optimization import perform_variable_minimization, perform_variable_maximization
from optmdfpathway import get_thermodynamic_bottlenecks
from toronto_set_up import cobra_model, optmdfpathway_base_milp_with_osmolarity, optmdfpathway_base_variables, biomass_reaction_id

MIN_MMDF = 0.01
optmdfpathway_base_variables["ATPM"].bounds(6.86,1000)
optmdfpathway_base_variables["var_B"].bounds(MIN_MMDF, 100)
theta = -0.01
growth_data = [
    (0.87699, MIN_MMDF),
    (0.87075, MIN_MMDF),
    (0.84777, MIN_MMDF),
    (0.80446, MIN_MMDF),
    (0.71739, MIN_MMDF)
]

for data_input in growth_data:
    print(" Load data...")
    current_maxgrowth = data_input[0]
    print(f" Setting minimal µ to {current_maxgrowth}")
    optmdfpathway_base_variables[biomass_reaction_id].bounds(current_maxgrowth, 100)
    optmdfpathway_base_variables["var_B"].bounds(data_input[1], 100)
    optmdfpathway_base_variables["osmolarity_sum_var"].bounds(0, 100)
    print(" Perform OptMDFpathway...")
    optmdfpathway_result = perform_variable_minimization(optmdfpathway_base_milp_with_osmolarity, "osmolarity_sum_var")
    if optmdfpathway_result["status"] == "Optimal":
        sum_met = optmdfpathway_result['values']['osmolarity_sum_var']
        current_phi = sum_met + theta
        print(f"The minimal total concentration is {current_phi}")
        print(f" *Perform bottleneck analysis")
        print(" Perform bottleneck analysis using OptMDFpathway base calculation")
        optmdfpathway_base_variables["osmolarity_sum_var"].bounds(current_phi, current_phi - theta)
        optmdfpathway_result_for_thermo = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, "var_B")
        bi, br = get_thermodynamic_bottlenecks(
        cobra_model=cobra_model,
        optmdfpathway_base_problem=optmdfpathway_base_milp_with_osmolarity,
        optmdfpathway_result=optmdfpathway_result_for_thermo,
        )

        with open(f"./Results/Bottleneck/bottleneck_analysis_{round(current_phi,5)}.txt", "w", encoding="utf-8") as f:
            f.write(br)