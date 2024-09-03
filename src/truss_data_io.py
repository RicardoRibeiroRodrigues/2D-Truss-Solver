# -*- coding: utf-8 -*-
"""
Module to handle the input and output of the truss data.
"""
from dataclasses import dataclass
import numpy as np
import xlrd


@dataclass
class SolverData:
    """
    Class to transfer data between the input file and the solver.
    """

    nodes_number: int
    nodes: np.ndarray
    elements_number: int
    incidence: np.ndarray
    loads_number: int
    loads: np.ndarray
    restraints_number: int
    restraints: np.ndarray


def import_data(input_name) -> SolverData:
    """
    Import the data from the input file and return a SolverData object.
    """
    file = xlrd.open_workbook(input_name)

    node_number, nodes = read_node_data(file)
    n_elements, incidence_matrix = read_element_data(file)
    n_loads, load_vector = read_load_data(file, node_number)
    num_restraints, restraints = read_restraint_data(file)

    output = SolverData(
        node_number,
        nodes,
        n_elements,
        incidence_matrix,
        n_loads,
        load_vector,
        num_restraints,
        restraints,
    )
    return output


def read_node_data(file):
    """
    Function to read the node data from the input file.
    """
    nos = file.sheet_by_name("Nos")
    node_number = int(nos.cell(1, 3).value)
    nodes = np.zeros((2, node_number))

    for c in range(node_number):
        nodes[0, c] = nos.cell(c + 1, 0).value
        nodes[1, c] = nos.cell(c + 1, 1).value
    return node_number, nodes


def read_element_data(file):
    """
    Function to read the element data from the input file.
    """
    incidence_sheet = file.sheet_by_name("Incidencia")
    n_elements = int(incidence_sheet.cell(1, 5).value)
    incidence_matrix = np.zeros((n_elements, 4))

    for c in range(n_elements):
        incidence_matrix[c, 0] = int(incidence_sheet.cell(c + 1, 0).value)
        incidence_matrix[c, 1] = int(incidence_sheet.cell(c + 1, 1).value)
        incidence_matrix[c, 2] = incidence_sheet.cell(c + 1, 2).value
        incidence_matrix[c, 3] = incidence_sheet.cell(c + 1, 3).value
    return n_elements, incidence_matrix


def read_load_data(file, node_number):
    """
    Function to read the load data from the input file.
    """
    load_sheet = file.sheet_by_name("Carregamento")

    n_loads = int(load_sheet.cell(1, 4).value)
    load_vector = np.zeros((node_number * 2, 1))

    for c in range(n_loads):
        no = load_sheet.cell(c + 1, 0).value
        xouy = load_sheet.cell(c + 1, 1).value
        gdl = int(no * 2 - (2 - xouy))
        load_vector[gdl - 1, 0] = load_sheet.cell(c + 1, 2).value
    return n_loads, load_vector


def read_restraint_data(file):
    """
    Function to read the restraint data from the input file.
    """
    restraints_sheet = file.sheet_by_name("Restricao")
    num_restraints = int(restraints_sheet.cell(1, 3).value)
    restrains = np.zeros((num_restraints, 1))

    for c in range(num_restraints):
        no = restraints_sheet.cell(c + 1, 0).value
        xouy = restraints_sheet.cell(c + 1, 1).value
        gdl = no * 2 - (2 - xouy)
        restrains[c, 0] = gdl - 1
    return num_restraints, restrains


def generate_output(
    output_name,
    reaction_forces,
    displacements,
    deformations,
    internal_forces,
    internal_tensions,
) -> None:
    """
    Write the output data to a text file.
    """
    output_name = output_name + ".txt"
    with open("saida.txt", "w+", encoding="utf-8") as f:
        f.write("Reacoes de apoio [N]\n")
        f.write(str(reaction_forces))
        f.write("\n\nDeslocamentos [m]\n")
        f.write(str(displacements))
        f.write("\n\nDeformacoes []\n")
        f.write(str(deformations))
        f.write("\n\nForcas internas [N]\n")
        f.write(str(internal_forces))
        f.write("\n\nTensoes internas [Pa]\n")
        f.write(str(internal_tensions))
