"""
Module to use the Solver class to solve the truss problem.
"""

from src.solver import Solver
from argparse import ArgumentParser

parser = ArgumentParser(description="Truss Solver")
parser.add_argument("--input", "-i", help="Input file path", default="data/entrada.xlsx")
parser.add_argument("--output", "-o", help="Output file path", default="saida.txt")
parser.add_argument("--save", "-s", help="Save the plot to a file", action="store_true")
args = parser.parse_args()

solver = Solver(args.input)
solver.solve(args.output, args.save)
