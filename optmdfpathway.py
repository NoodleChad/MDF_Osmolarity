#!/usr/bin/env python3
"""[summary]"""

## IMPORTS ##
# External
import cobra
import pickle
import copy
import concurrent.futures
import pulp
import ray
from math import exp, log
from typing import Any, List, Dict
# Internal
from fba import get_fba_base_problem
from linear_approximation import add_linear_function_approximation
from optimization import perform_variable_maximization, perform_variable_minimization, solve_current_problem


## CONSTANTS ##
M = 10_000
STANDARD_R = 8.314e-3  # kJ�K{1�mol{1 (standard value is in J�K{1�mol{1)
STANDARD_T = 298.15  # K


## PUBLIC FUNCTIONS ##
def add_osmolarity_constraints(optmdfpathway_base_problem: pulp.LpProblem,
                               max_total_metabolite_concentration: float,
                               metabolite_weights: Dict[str, float],
                               max_relative_exp_linearization_error: float
                               ):
    """[summary]

    Args:
        optmdfpathway_base_problem (pulp.LpProblem): [description]
        max_total_metabolite_concentration (float): [description]
        metabolite_weights (Dict[str, float]): [description]
        max_relative_exp_linearization_error (float): [description]

    Returns:
        [type]: [description]
    """
    optmdfpathway_base_variables: Dict[str, pulp.LpVariable] = optmdfpathway_base_problem.variablesDict()
    sum_expression = pulp.LpAffineExpression(name="concentration_sum")
    for variable in optmdfpathway_base_problem.variables():
        variable_name = variable.name
        if not variable_name.startswith("x_"):
            continue

        metabolite_id = variable_name[2:]
        ln_min_concentration = optmdfpathway_base_variables[variable_name].getLb()
        ln_max_concentration = optmdfpathway_base_variables[variable_name].getUb()

        if (ln_min_concentration == 0.0) and (ln_max_concentration == 0.0):
            continue

        # Define exp-approximating variable and add it to model
        optmdfpathway_base_problem, y_variable, num_approximations = add_linear_function_approximation(
            base_problem=optmdfpathway_base_problem,
            x_variable=variable,
            y_variable_name=f"expx_{metabolite_id}",
            function_to_approximate=lambda x: exp(x),
            function_derivative=lambda x: exp(x),
            min_x=ln_min_concentration,
            max_x=ln_max_concentration,
            max_relative_error=max_relative_exp_linearization_error,
            is_minimum_for_y=True,
        )

        # Determine the metabolite's osmolarity concentration weight
        if metabolite_id in metabolite_weights.keys():
            weight = metabolite_weights[metabolite_id]
        else:
            weight = 1.0

        # Add sum expression
        sum_expression += weight * y_variable

    # Sum constraint
    osmolarity_sum_var = pulp.LpVariable(
        name="osmolarity_sum_var",
        lowBound=0,
        upBound=max_total_metabolite_concentration,
    )

    osmolarity_sum_constraint = osmolarity_sum_var - sum_expression == 0
    optmdfpathway_base_problem += osmolarity_sum_constraint

    # Add the ion pump constraint
    psi = 1020 # Slope have to be -5 and 80 is -.3
    high_osm = pulp.LpAffineExpression(name="high_osm")
    Osmotic_demand = pulp.LpVariable(
        name="Osmotic_ATP_demand",
        lowBound=0,
        upBound=1000,
    )
    Phi_high_osm = pulp.LpVariable(
        name="Phi_high_osm",
        lowBound=0,
        upBound=1000,
    )
    Phi_threshold = pulp.LpVariable(
        name="no_cost_until",
        lowBound=0,
        upBound=psi*1000
    )
    float_h_osm = pulp.LpVariable(
        name="high_osm_float",
        lowBound=0,
        upBound=1000,
    )
    high_osm = -psi*Phi_high_osm + Osmotic_demand + 1*Phi_threshold  + -1*float_h_osm == 0

    optmdfpathway_base_problem += high_osm
    # define osmotic demand metabolites
    K_free_export = pulp.LpVariable(
        name="Free_transport_out",
        lowBound=0,
        upBound=1000,
    )
    optmdfpathway_base_problem.constraints['atp_c'] += -1*Osmotic_demand
    optmdfpathway_base_problem.constraints['h2o_c'] += -1*Osmotic_demand
    optmdfpathway_base_problem.constraints['adp_c'] += 1*Osmotic_demand
    optmdfpathway_base_problem.constraints['h_c'] += 1*Osmotic_demand
    optmdfpathway_base_problem.constraints['pi_c'] += 1*Osmotic_demand    
    return optmdfpathway_base_problem

def get_optmdfpathway_base_problem(cobra_model: cobra.Model, dG0_values: Dict[str, Dict[str, float]],
                                   metabolite_concentration_values: Dict[str, Dict[str, float]],
                                   R: float, T: float,
                                   extra_constraints: List[Dict[str, float]],
                                   sub_network_ids: List[str]=[],
                                   add_optmdf_bottleneck_analysis: bool=False) -> pulp.LpProblem:
    """[summary]

    Args:
        cobra_model (cobra.Model): [description]
        dG0_values (Dict[str, Dict[str, float]]): [description]
        metabolite_concentration_values (Dict[str, Dict[str, float]]): [description]
        R (float): [description]
        T (float): [description]
        extra_constraints (List[Dict[str, float]]): [description]

    Returns:
        pulp.LpProblem: [description]
    """
    # Concentrations in M
    # T in K
    # R in kJ�K{1�mol{1 (standard value is in J�K{1�mol{1)
    # Get FBA base
    base_problem = get_fba_base_problem(
        cobra_model=cobra_model,
        extra_constraints=extra_constraints,
    )
    fba_base_variables = base_problem.variablesDict()

    # Min/max concentrations
    concentration_vars: Dict[str, pulp.LpVariable] = {}
    for metabolite in cobra_model.metabolites:
        if metabolite.id in metabolite_concentration_values.keys():
            id_key = metabolite.id
        else:
            id_key = "DEFAULT"

        min_concentration = log(metabolite_concentration_values[id_key]["min"])
        max_concentration = log(metabolite_concentration_values[id_key]["max"])

        concentration_vars[metabolite.id] = pulp.LpVariable(
            name="x_"+metabolite.id,
            lowBound=min_concentration,
            upBound=max_concentration,
            cat=pulp.LpContinuous,
        )

    # df variables (f_r)
    var_B = pulp.LpVariable(
        name="var_B",
        lowBound=-float("inf"),
        upBound=float("inf"),
        cat=pulp.LpContinuous,
    )

    has_subnetwork = len(sub_network_ids) > 0
    if has_subnetwork:
        var_B2 = pulp.LpVariable(
            name="var_B2",
            lowBound=-float("inf"),
            upBound=float("inf"),
            cat=pulp.LpContinuous,
        )

    zb_sum: pulp.LpAffineExpression = 0.0
    z_sum: pulp.LpAffineExpression = 0.0
    for reaction in cobra_model.reactions:
        if reaction.id not in dG0_values.keys():
            continue

        current_f_variable = pulp.LpVariable(
            name=f"f_var_{reaction.id}",
            lowBound=-float("inf"),
            upBound=float("inf"),
            cat=pulp.LpContinuous,
        )

        # f_r = -(�G'0_r + RT * S^(t)_.,r * x)
        # <=> -�G'0_r - - f_r RT * S^(t)_.,r * x = 0
        dG0_value = dG0_values[reaction.id]["dG0"]
        dG0_randomfixed_variable = pulp.LpVariable(
            name=f"dG0_{reaction.id}",
            cat=pulp.LpContinuous,
        )
        dG0_randomfixed_variable.setInitialValue(dG0_value)
        dG0_randomfixed_variable.fixValue()

        # e.g.,
        #     RT * ln(([S]*[T])/([A]�*[B]))
        # <=> RT*ln([S]) + RT*ln([T]) - 2*RT*ln([A]) - RT*ln([B])
        reaction_f_expression = -current_f_variable
        for metabolite in reaction.metabolites:
            stoichiometry = reaction.metabolites[metabolite]
            reaction_f_expression -= stoichiometry * R * T * concentration_vars[metabolite.id]
        reaction_f_expression -= dG0_randomfixed_variable
        reaction_f_constraint: pulp.LpConstraint = reaction_f_expression == 0

        base_problem += reaction_f_constraint

        # OptMDFpathway binary variable introduction
        current_z_variable = pulp.LpVariable(
            name=f"z_var_{reaction.id}",
            cat=pulp.LpBinary,
        )
        z_sum += current_z_variable
        if add_optmdf_bottleneck_analysis:
            current_zb_variable = pulp.LpVariable(
                name=f"zb_var_{reaction.id}",
                cat=pulp.LpBinary,
            )

        # Find flux variable
        current_flux_variable = fba_base_variables[reaction.id]

        # z_r = 0 -> v_r = 0
        # In Big M form: v_r <= M*z_r <=> v_r - M*z_r <= 0
        z_zero_constraint: pulp.LpConstraint = current_flux_variable - M*current_z_variable <= 0
        # z_r = 1 -> f_r >= B == f_r - B >= 0
        # In Big M form: f_r + (1-z_i)*M >= B
        z_one_constraint: pulp.LpConstraint = 0.0
        if (not add_optmdf_bottleneck_analysis):
            z_one_constraint = var_B <= current_f_variable + (1-current_z_variable)*M
        else:
            z_one_constraint = var_B <= current_f_variable + (1-current_z_variable)*M + current_zb_variable*5_000
            zb_sum += current_zb_variable

        base_problem += z_zero_constraint
        base_problem += z_one_constraint

        if has_subnetwork:
            if reaction.id in sub_network_ids:
                z_one_constraint_2: pulp.LpConstraint = var_B2 <= current_f_variable + (1-current_z_variable)*M
                base_problem += z_one_constraint_2

    if add_optmdf_bottleneck_analysis:
        zb_sum_var = pulp.LpVariable(
            name=f"zb_sum_var",
            lowBound=-float("inf"),
            upBound=float("inf"),
            cat=pulp.LpContinuous,
        )
        base_problem += zb_sum_var == zb_sum
    z_sum_var = pulp.LpVariable(
        name="z_sum_var",
        lowBound=-float("inf"),
        upBound=float("inf"),
        cat=pulp.LpContinuous
        )
    
    base_problem += z_sum_var == z_sum



    return base_problem

def get_thermodynamic_bottlenecks(cobra_model: cobra.Model,
                                  optmdfpathway_base_problem: pulp.LpProblem,
                                  optmdfpathway_result: Dict[str, Dict[str, float]],
                                  **kwargs) -> List[str]:
    """[summary]

    Args:
        cobra_model (cobra.Model): [description]
        optmdfpathway_base_problem (pulp.LpProblem): [description]
        original_max_mdf (float): [description]

    Returns:
        List[str]: [description]
    """
    epsilon = 1e-3
    optmdfpathway_base_variables: Dict[str, pulp.LpVariable] = optmdfpathway_base_problem.variablesDict()

    original_max_mdf = optmdfpathway_result["values"]["var_B"]
    optmdfpathway_base_variables["var_B"].bounds(original_max_mdf, 10e9)

    checked_ids = []
    optmdfs = {}
    thermodynamic_bottleneck_ids: List[str] = []
    for reaction in cobra_model.reactions:
        if abs(optmdfpathway_result["values"][reaction.id]) < 1e-6:
            continue

        f_var_name = f"f_var_{reaction.id}"

        if f_var_name not in optmdfpathway_result["values"].keys():
            continue

        if abs(optmdfpathway_result["values"][f_var_name] - original_max_mdf) > 0.5:
            continue

        # Test result
        try:
            current_dG0_variable = optmdfpathway_base_variables[f"dG0_{reaction.id}"]
        except KeyError:
            continue

        original_dG0 = current_dG0_variable.lowBound
        current_dG0_variable.bounds(original_dG0-500, original_dG0-500)
        bottleneck_problem_result = perform_variable_maximization(optmdfpathway_base_problem, "var_B", optmdfpathway_base_variables, **kwargs)
        current_dG0_variable.bounds(original_dG0, original_dG0)

        if bottleneck_problem_result["status"] != "Optimal":
            # print(bottleneck_problem_result["status"])
            # print("ERROR in OptMDFpathway with ", reaction.id)
            continue
        checked_ids.append(reaction.id)
        new_optmdf = bottleneck_problem_result["values"]["var_B"]
        if new_optmdf > (original_max_mdf + epsilon):
            optmdfs[reaction.id] = new_optmdf
            thermodynamic_bottleneck_ids.append(reaction.id)

    reactions: List[cobra.Reaction] = [cobra_model.reactions.get_by_id(x) for x in thermodynamic_bottleneck_ids]
    metabolite_data: Dict[str, List[str]] = {}
    bottleneck_report = f"===BOTTLENECK REPORT===\n"
    bottleneck_report += f">Number of bottleneck reactions: {len(reactions)}\n"
    bottleneck_report += f">Original OptMDF: {original_max_mdf} kJ/mol\n"
    counter = 0
    for reaction in reactions:
        current_dG0_variable = optmdfpathway_base_variables[f"dG0_{reaction.id}"]
        bottleneck_report += f">Bottleneck reaction no. {counter+1}:\n"
        bottleneck_report += f" *ID: {reaction.id}\n"
        bottleneck_report += f" *�G'�: {round(current_dG0_variable.lowBound, 4)} kJ/mol\n"
        bottleneck_report += f" *Reached OptMDF without this bottleneck: {round(optmdfs[reaction.id], 4)} kJ/mol\n"
        bottleneck_report += f" *Reaction string: {reaction.reaction}\n"
        metabolite_ids = [x.id for x in reaction.metabolites]

        for metabolite_id in metabolite_ids:
            if metabolite_id not in metabolite_data.keys():
                metabolite_data[metabolite_id] = []
            metabolite_data[metabolite_id].append(reaction.id)
        counter += 1

    counter = 1
    bottleneck_report += f">Reaction-connecting metabolites:\n"
    for metabolite_id in metabolite_data.keys():
        reaction_ids = metabolite_data[metabolite_id]
        if len(reaction_ids) > 1:
            bottleneck_report += f" *Metabolite no. {counter}: {metabolite_id}\n"
            bottleneck_report += f"  in reactions: "
            bottleneck_report += ", ".join(reaction_ids)
            bottleneck_report += "\n"
            counter += 1

    reaction_clusters = []
    for reaction_ids in metabolite_data.values():
        cluster_expanded = False
        cluster_index = 0
        for reaction_cluster in reaction_clusters:
            if len(set(reaction_cluster).intersection(set(reaction_ids))) > 0:
                reaction_cluster_expanded = list(set(reaction_cluster + reaction_ids))
                reaction_clusters[cluster_index] = reaction_cluster_expanded
                cluster_expanded = True
            cluster_index += 1
        if not cluster_expanded:
            reaction_clusters.append(copy.deepcopy(reaction_ids))

        new_clusters = []
        for reaction_cluster_1 in reaction_clusters:
            new_cluster = copy.deepcopy(reaction_cluster_1)
            for reaction_cluster_2 in reaction_clusters:
                if len(set(reaction_cluster_1).intersection(set(reaction_cluster_2))) > 1:
                    new_cluster += set(new_cluster + reaction_cluster_2)
            new_clusters.append(sorted(list(set(new_cluster))))
        reaction_clusters = copy.deepcopy(new_clusters)

    reaction_clusters = list(set([tuple(x) for x in reaction_clusters]))

    bottleneck_report += ">Metabolite-connected reaction clusters:\n"
    counter = 1
    for cluster in reaction_clusters:
        bottleneck_report += f"*Cluster no. {counter}:\n"
        bottleneck_report += " " + ", ".join(cluster) + "\n"
        counter += 1

    #######
    if thermodynamic_bottleneck_ids == []:
        print("Searching cluster...")
        full_cluster = copy.deepcopy(checked_ids)

        original_dG0s = {}
        print(full_cluster)
        for reaction_id in full_cluster:
            current_dG0_variable = optmdfpathway_base_variables[f"dG0_{reaction_id}"]
            original_dG0s[reaction_id] = current_dG0_variable.lowBound
            current_dG0_variable.bounds(original_dG0s[reaction_id]-500, original_dG0s[reaction_id]-500)
            print(reaction_id)

        # Check all
        bottleneck_problem_result = perform_variable_maximization(optmdfpathway_base_problem, "var_B", optmdfpathway_base_variables, **kwargs)
        new_optmdf = bottleneck_problem_result["values"]["var_B"]
        print(new_optmdf)
        if not (new_optmdf > (original_max_mdf + epsilon)):
            print(original_max_mdf)
            print(new_optmdf)
            print("!ERROR 1!")

        # Check each reaction
        final_cluster = []
        for reaction_id in full_cluster:
            current_dG0_variable.bounds(original_dG0s[reaction_id], original_dG0s[reaction_id])
            bottleneck_problem_result = perform_variable_maximization(optmdfpathway_base_problem, "var_B", optmdfpathway_base_variables, **kwargs)
            new_optmdf = bottleneck_problem_result["values"]["var_B"]
            if bottleneck_problem_result["status"] != "Optimal":
                print("!ERROR 2!")
                break
            if not (new_optmdf > (original_max_mdf + epsilon)):
                final_cluster.append(reaction_id)
            current_dG0_variable.bounds(original_dG0s[reaction_id]-500, original_dG0s[reaction_id]-500)

        thermodynamic_bottleneck_ids = final_cluster
    #######

    #######
    optmdfpathway_base_variables["var_B"].bounds(-10e9, 10e9)

    return thermodynamic_bottleneck_ids, bottleneck_report


def perform_optmdfpathway_mdf_maximization(optmdfpathway_base_problem: pulp.LpProblem, **kwargs) -> Dict[str, Any]:
    results = perform_variable_maximization(
        base_problem=optmdfpathway_base_problem,
        variable_id="var_B",
        **kwargs,
    )
    return results


def add_differential_reactions_constraints(optmdfpathway_base_problem: pulp.LpProblem,
                                           in_vivo_solution: Dict[str, int]) -> pulp.LpProblem:
    optmdfpathway_base_variables: Dict[str, pulp.LpVariable] = optmdfpathway_base_problem.variablesDict()

    error_sum = 0.0
    for reaction_id in in_vivo_solution.keys():
        reaction_error_var = pulp.LpVariable(
            name=f"z_error_{reaction_id}",
            cat=pulp.LpBinary,
        )
        try:
            current_z_variable = optmdfpathway_base_variables[f"z_var_{reaction_id}"]
        except KeyError:
            print(">", reaction_id)
            continue
        current_flux_variable = optmdfpathway_base_variables[f"f_var_{reaction_id}"]
        # if in_vivo_solution[reaction_id] == 0:
        #     optmdfpathway_base_problem += reaction_error_var == z_var
        # else:
        #     optmdfpathway_base_problem += reaction_error_var == 1.0 - z_var
        # error_sum += reaction_error_var
        ####
        current_z2_variable = pulp.LpVariable(
            name=f"z2_var_{reaction_id}",
            cat=pulp.LpBinary,
        )

        current_error_variable = pulp.LpVariable(
            name=f"g_var_{reaction_id}",
            lowBound=-float("inf"),
            upBound=float("inf"),
            cat=pulp.LpContinuous,
        )
        z2_zero_constraint: pulp.LpConstraint = current_z2_variable - M*current_flux_variable <= 0
        if in_vivo_solution[reaction_id] == 0:
            error_z_z2_constraint: pulp.LpConstraint = current_error_variable - current_z_variable + current_z2_variable * M <= M
            error_z2_constraint_1: pulp.LpConstraint = -current_error_variable + current_z_variable + current_z2_variable * M <= M
        else:
            error_z_z2_constraint: pulp.LpConstraint = current_error_variable - (1.0-current_z_variable) + current_z2_variable * M <= M
            error_z2_constraint_1: pulp.LpConstraint = -current_error_variable + (1.0-current_z_variable) + current_z2_variable * M <= M
        error_z2_constraint_2: pulp.LpConstraint = current_error_variable - current_z2_variable * M <= 0
        error_z2_constraint_3: pulp.LpConstraint = -current_error_variable - current_z2_variable * M <= 0

        optmdfpathway_base_problem += z2_zero_constraint
        optmdfpathway_base_problem += error_z_z2_constraint
        optmdfpathway_base_problem += error_z2_constraint_1
        optmdfpathway_base_problem += error_z2_constraint_2
        optmdfpathway_base_problem += error_z2_constraint_3

        error_sum += current_error_variable
        ####

    error_sum_var = pulp.LpVariable(
        name=f"reaction_error_sum",
        cat=pulp.LpContinuous,
    )
    optmdfpathway_base_problem += error_sum_var == error_sum

    return optmdfpathway_base_problem


def get_z_variable_status(optmdfpathway_solution: Dict[str, Any], necessary_id_part: str="") -> Dict[str, int]:
    z_variable_status: Dict[str, int] = {}
    for var_id in optmdfpathway_solution["values"].keys():
        if var_id.startswith("z_var_"):
            if not necessary_id_part in var_id:
                continue
            z_float = optmdfpathway_solution["values"][var_id]
            if z_float < 1e-5:
                z_int = 0
            elif z_float > 0.99:
                z_int = 1
            else:
                print(f"ERROR in get_z_variable_status():"
                      f"Value {z_float} for {var_id} is no reasonable binary value")
                raise Exception
            z_variable_status[var_id] = z_int

    return z_variable_status


def add_differential_concentration_constraints(cobra_model: cobra.Model,
                                               optmdfpathway_base_problem: pulp.LpProblem,
                                               in_vivo_solution: Dict[str, float],
                                               **kwargs) -> pulp.LpProblem:
    optmdfpathway_base_variables: Dict[str, pulp.LpVariable] = optmdfpathway_base_problem.variablesDict()

    solve_current_problem(optmdfpathway_base_problem)

    error_sum = 0.0
    for metabolite in cobra_model.metabolites:
        if len(metabolite.reactions) == 0:
            continue

        if metabolite.id not in in_vivo_solution.keys():
            continue

        x_var_id = f"x_{metabolite.id}"
        try:
            optmdfpathway_base_variables[x_var_id]
        except KeyError:
            continue

        error_variable = pulp.LpVariable(
            name=f"error_{metabolite.id}",
            cat=pulp.LpContinuous,
        )
        optmdfpathway_base_problem += -error_variable <= log(in_vivo_solution[metabolite.id]) - optmdfpathway_base_variables[x_var_id]
        optmdfpathway_base_problem += log(in_vivo_solution[metabolite.id]) - optmdfpathway_base_variables[x_var_id] <= error_variable

        error_sum += 1.0 * error_variable

    error_sum_var = pulp.LpVariable(
        name=f"error_sum",
        cat=pulp.LpContinuous,
    )
    optmdfpathway_base_problem += error_sum_var == error_sum

    return optmdfpathway_base_problem


def perform_concentration_variability_analysis(cobra_model: cobra.Model,
                                               optmdfpathway_base_problem: pulp.LpProblem,
                                               min_mdf: float,
                                               selected_metabolites: List[str]=[],
                                               **kwargs) -> Dict[str, Dict[str, float]]:
    """[summary]

    Args:
        cobra_model (cobra.Model): [description]
        optmdfpathway_base_problem (pulp.LpProblem): [description]
        min_mdf (float): [description]

    Returns:
        Dict[str, Dict[str, float]]: [description]
    """
    optmdfpathway_base_variables: Dict[str, pulp.LpVariable] = optmdfpathway_base_problem.variablesDict()
    var_B = optmdfpathway_base_variables["var_B"]
    optmdfpathway_base_problem += var_B >= min_mdf

    solve_current_problem(optmdfpathway_base_problem)

    concentration_bounds: Dict[str, Dict[str, float]] = {}
    for metabolite in cobra_model.metabolites:
        if len(metabolite.reactions) == 0:
            continue
        if selected_metabolites != []:
            if metabolite.id not in selected_metabolites:
                continue

        x_var_id = f"x_{metabolite.id}"
        try:
            optmdfpathway_base_variables[x_var_id]
        except KeyError:
            print("KEYERROR")
            continue

        concentration_bounds[metabolite.id] = {}

        # Minimization
        min_result = perform_variable_minimization(
            base_problem=optmdfpathway_base_problem,
            variable_id=x_var_id,
            **kwargs
        )
        if min_result["status"] == "Optimal":
            min_value = min_result["values"][x_var_id]
        else:
            min_value = float("NaN")
        concentration_bounds[metabolite.id]["min"] = exp(min_value)

        # Maximization
        max_result = perform_variable_maximization(
            base_problem=optmdfpathway_base_problem,
            variable_id=x_var_id,
            **kwargs
        )
        if max_result["status"] == "Optimal":
            max_value = max_result["values"][x_var_id]
        else:
            max_value = float("NaN")
        concentration_bounds[metabolite.id]["max"] = exp(max_value)

    return concentration_bounds


def _mt_perform_concentration_variability_analysis_multi(metabolite: cobra.Metabolite, base_problem_dict: Dict[Any, Any]):
    variables_dict, optmdfpathway_base_problem = pulp.LpProblem.from_dict(base_problem_dict)

    x_var_id = f"x_{metabolite.id}"
    try:
        variables_dict[x_var_id]
    except KeyError:
        return (metabolite.id, "NO_VARIABLE", "NO_VARIABLE")

    # Minimization
    min_result = perform_variable_minimization(
        base_problem=optmdfpathway_base_problem,
        variable_id=x_var_id
    )
    if min_result["status"] == "Optimal":
        min_value = exp(min_result["values"][x_var_id])
    else:
        min_value = float("NaN")

    # Maximization
    max_result = perform_variable_maximization(
        base_problem=optmdfpathway_base_problem,
        variable_id=x_var_id
    )
    if max_result["status"] == "Optimal":
        max_value = exp(max_result["values"][x_var_id])
    else:
        max_value = float("NaN")

    range = abs(min_value - max_value)

    return (metabolite.id, min_value, max_value, range)


@ray.remote(num_cpus=2, max_retries=1)
def _ray_perform_concentration_variability_analysis_multi(metabolites: List[cobra.Metabolite], base_problem_dict: Dict[Any, Any]):
    with concurrent.futures.ThreadPoolExecutor(max_workers=None) as executor:
        futures = [executor.submit(_mt_perform_concentration_variability_analysis_multi, metabolite, base_problem_dict) for metabolite in metabolites]
        results = [future.result() for future in futures]

    return results


def perform_concentration_variability_analysis_multi(cobra_model: cobra.Model,
                                                     optmdfpathway_base_problem: pulp.LpProblem,
                                                     min_mdf: float) -> Dict[str, Dict[str, float]]:
    """[summary]

    Args:
        cobra_model (cobra.Model): [description]
        optmdfpathway_base_problem (pulp.LpProblem): [description]
        min_mdf (float): [description]

    Returns:
        Dict[str, Dict[str, float]]: [description]
    """
    optmdfpathway_base_variables: Dict[str, pulp.LpVariable] = optmdfpathway_base_problem.variablesDict()
    optmdfpathway_base_variables["var_B"].bounds(min_mdf, 10e9)

    solve_current_problem(optmdfpathway_base_problem)
    base_problem_dict = optmdfpathway_base_problem.to_dict()

    batch_size = 10
    loop_values = []
    start_value = 0
    while start_value < len(cobra_model.metabolites):
        loop_values.append(start_value)
        start_value += batch_size
    futures = [_ray_perform_concentration_variability_analysis_multi.remote(cobra_model.metabolites[loop_value:loop_value+batch_size], base_problem_dict)
               for loop_value in loop_values]
    results = ray.get(futures)

    out_results = []
    for result in results:
        out_results += result

    optmdfpathway_base_variables["var_B"].bounds(-10e9, 10e9)

    return out_results