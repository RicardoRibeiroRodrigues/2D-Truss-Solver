"""
Module to use the Solver class to solve the truss problem.
"""

from src.solver import Solver

solver = Solver("data/entrada2.xlsx")
solver.solve("saida.txt")
