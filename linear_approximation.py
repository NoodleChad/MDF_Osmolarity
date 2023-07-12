#!/usr/bin/env python3
"""[summary]"""

## IMPORTS ##
# External
import numpy
import pulp
from typing import Any, Callable, List, Tuple


## PUBLIC FUNCTIONS ##
def add_linear_function_approximation(base_problem: pulp.LpProblem,
                                      x_variable: pulp.LpVariable,
                                      y_variable_name: str,
                                      function_to_approximate: Callable[[float], float],
                                      function_derivative: Callable[[float], float],
                                      min_x: float,
                                      max_x: float,
                                      max_relative_error: float,
                                      is_minimum_for_y: bool):
    """[summary]

    Args:
        base_problem (pulp.LpProblem): [description]
        x_variable (pulp.LpVariable): [description]
        y_variable_name (str): [description]
        function_to_approximate (Callable[[float], float]): [description]
        function_derivative (Callable[[float], float]): [description]
        min_x (float): [description]
        max_x (float): [description]
        max_relative_error (float): [description]
        is_minimum_for_y (bool): [description]

    Returns:
        [type]: [description]
    """
    current_num_sections = 2
    is_above_error = True
    while is_above_error:
        section_xs = numpy.linspace(min_x, max_x, current_num_sections)
        step_size = section_xs[1] - section_xs[0]
        is_above_error = False
        linear_approximations: List[Tuple[float, float]] = []
        for section_x in section_xs:
            function_y = function_to_approximate(section_x)
            derivative_y = function_derivative(section_x)
            m_value = derivative_y
            b_value = function_y - derivative_y * section_x

            left_x = section_x - step_size/2
            right_x = section_x + step_size/2
            for test_x in (left_x, right_x):
                try:
                    function_y = function_to_approximate(test_x)
                except ValueError:
                    continue
                approximation_y = m_value * test_x + b_value
                error_absolute = function_y - approximation_y
                error_relative = abs(error_absolute / function_y)
                if error_relative > max_relative_error:
                    is_above_error = True
                    break

            if is_above_error:
                current_num_sections += 1
                break

            linear_approximations.append((m_value, b_value))

    y_variable = pulp.LpVariable(
        name=y_variable_name,
        cat=pulp.LpContinuous,
    )
    for linear_approximation in linear_approximations:
        m_value = linear_approximation[0]
        b_value = linear_approximation[1]
        if is_minimum_for_y:
            base_problem += y_variable >= m_value * x_variable + b_value
        else:
            base_problem += y_variable <= m_value * x_variable + b_value

    return base_problem, y_variable, len(linear_approximations)
