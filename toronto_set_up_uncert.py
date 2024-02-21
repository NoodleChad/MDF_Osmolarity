import cobra
from fba import get_fba_base_problem, perform_fba_flux_maximization
from helper import ensure_folder_existence, json_load, json_write, pickle_load, pickle_write, save_xy_point_plot
from optmdfpathway_uncert import STANDARD_R, STANDARD_T, add_osmolarity_constraints, get_optmdfpathway_base_problem, perform_optmdfpathway_mdf_maximization, get_thermodynamic_bottlenecks
from optimization import perform_optimization, perform_stepwise_variable_optimization, perform_variable_maximization, perform_variable_minimization, solve_current_problem

MIN_MMDF = 0.1

print("==Setting up iML model with osmolarity constraints==")
print(">Read SBML file")
cobra_model = cobra.io.read_sbml_model("resources/iML1515_irreversible_cleaned.xml")

print(">Lower minimal ATPM demand to 0 mmol/(gDW*h)")
cobra_model.reactions.get_by_id("ATPM").lower_bound = 0

print(">Perform test FBA...")
biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"
print(f" Selected biomass reaction: {biomass_reaction_id}")
fba_base_problem = get_fba_base_problem(cobra_model=cobra_model, extra_constraints=[])
fba_result = perform_fba_flux_maximization(base_problem=fba_base_problem, reaction_id=biomass_reaction_id)
precise_max_growth = fba_result["values"][biomass_reaction_id]
used_max_growth = float(str(precise_max_growth)[:5])

print(f" Precise max growth is {precise_max_growth}")
print(f" Used rounded and floored max growth is {used_max_growth}")

print(">Load precomputed eQuilibrator dG0 values")
dG0_values = json_load("./resources/dG0_iML1515_irreversible_cleaned.json")

print(">Delete precomputed multi-compartmental dG0 values")
dG0_value_ids = list(dG0_values.keys())
for dG0_value_id in dG0_value_ids:
    if dG0_values[dG0_value_id]["num_compartments"] > 1:
        del(dG0_values[dG0_value_id])

print(">Load predefined concentration range values")
concentration_values = {
    "DEFAULT": {
        "min": 1e-6,
        "max": 0.1,
    },
    "h2o_c": {
        "min": 1.0,
        "max": 1.0,
    },
    "h2o_p": {
        "min": 1.0,
        "max": 1.0,
    },
    "h2o_e": {
        "min": 1.0,
        "max": 1.0,
    },
    "h_c": {
        "min": 1.0,
        "max": 1.0,
    },
    "h_p": {
        "min": 1.0,
        "max": 1.0,
    },
    "h_e": {
        "min": 1.0,
        "max": 1.0,
    }
}

####
print(">Perform bottleneck test...")
optmdfpathway_base_problem_bottlenecktest = get_optmdfpathway_base_problem(
    cobra_model=cobra_model,
    dG0_values=dG0_values,
    metabolite_concentration_values=concentration_values,
    R=STANDARD_R,
    T= STANDARD_T,
    extra_constraints=[],
    add_optmdf_bottleneck_analysis=True,
)
optmdfpathway_base_variables = optmdfpathway_base_problem_bottlenecktest.variablesDict()
optmdfpathway_base_variables[biomass_reaction_id].bounds(
    used_max_growth,
    1e12
)
optmdfpathway_base_variables["var_B"].bounds(
    0.1,
    1e12
)

optmdfpathway_result = perform_variable_minimization(
    optmdfpathway_base_problem_bottlenecktest,
    "zb_sum_var"
)
print("ZB SUM ANALYSIS")
print("STATUS:", optmdfpathway_result["status"])
print("zb sum:", optmdfpathway_result["values"]["zb_sum_var"], "reactions")
print("var_B:", optmdfpathway_result["values"]["var_B"], "kJ/mol")


# These reactions are removed because they make growth infeasible in typical growth condition.
# We assume the incertainty on these values is to high to make them reliable and therefore, should not be considered in the solution.
del(dG0_values["ATPPRT"])
del(dG0_values["SHCHD2"])
del(dG0_values["MECDPS"])
del(dG0_values["DHPPDA2"])
del(dG0_values["MCTP1App"])
del(dG0_values["KDOCT2"])
del(dG0_values["IG3PS"])
del(dG0_values["AIRC3_REV"])
del(dG0_values["AIRC3_FWD"])
#del(dG0_values["MALCOAMT"])

print(">Get base OptMDFpathway MILP (without osmolarity constraints)...")
optmdfpathway_base_problem = get_optmdfpathway_base_problem(
    cobra_model=cobra_model,
    dG0_values=dG0_values,
    metabolite_concentration_values=concentration_values,
    R=STANDARD_R,
    T= STANDARD_T,
    extra_constraints=[]
)

print(">Add osmolarity constraints to OptMDFpathway base problem with a very high K")
optmdfpathway_base_milp_with_osmolarity = add_osmolarity_constraints(
    optmdfpathway_base_problem=optmdfpathway_base_problem,
    max_total_metabolite_concentration=100.0,
    metabolite_weights={},
    max_relative_exp_linearization_error=0.05,
)

print(">Load model variables dictionary")
optmdfpathway_base_variables = optmdfpathway_base_milp_with_osmolarity.variablesDict()

print(f">Test OptMDFpathway run by temporarily setting µ>{used_max_growth}...")
optmdfpathway_base_variables[biomass_reaction_id].bounds(used_max_growth, 1e12)
optmdfpathway_test_result = perform_variable_maximization(optmdfpathway_base_milp_with_osmolarity, "var_B")
print(" Status:", optmdfpathway_test_result["status"])
print(" OptMDF:", optmdfpathway_test_result["values"]["var_B"])

print(">Ensure existence of main analyses folder")
ensure_folder_existence("./toronto_results")