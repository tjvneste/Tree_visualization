#!/usr/bin/env python3
'''
Script for assignment PathoSense:
Takes a Newick file that includes the sample strain and
visualizes a reduced tree based on a input distance in pdf.

Example of how to run the script:
./tree_analysis.py -f full_alignment.fasta.treefile -o pruned_tree.pdf
--sample_strain Porcine_parvovirus_1_1-0004394_BEL_Dec2022_30X Porcine
--type_strains PPV1_IDT_DEU_1964 PPV1_27a_DEU_2001
PPV1_NADL2_USA_1964 PPV1_NADL8_USA_1976
--distance 0.0007

Author: Tristan Vanneste
Date: February 2023
'''

# libraries
from ete3 import Tree, TreeStyle, NodeStyle
import statistics
import argparse


def parse_arguments():
    """Parses and returns arguments from the command line

    Returns:
        parser.parse_args
    """

    parser = argparse.ArgumentParser(description='')

    parser.add_argument("-f",
                        dest="input_file",
                        required=True,
                        type=str,
                        help="input Newick file.")
    parser.add_argument("-o",
                        dest="output_name",
                        required=True,
                        type=str,
                        help="Output name for the pdf.")
    parser.add_argument("--sample_strains",
                        dest="sample_strains",
                        nargs="*",
                        default="",
                        type=str,
                        help="Sample strains.")
    parser.add_argument("--type_strains",
                        dest="type_strains",
                        nargs="*",
                        default="",
                        type=str,
                        help="Type strains.")
    parser.add_argument("--distance",
                        dest="distance",
                        default=0.001,
                        type=float,
                        help="The distance where you want to filter on.")
    # Parse arguments and call main
    return parser.parse_args()

# Helper function used by main


def prune_tree_by_distance(tree: str, dist: float, sample_strains: list, type_strains: list):
     """Prunes the tree by detaching leaves/nodes that have a mean distances to their node
          below a certain threshold (dist).

     Arguments:
          tree: The filename for the Newick tree
          dist: The distance where you want to filter on
          sample_strains: The sample strains
          type_strains: The type strains
     """
     for node in tree.traverse("preorder"):
          distance = []
          for leave in node:
               distance.append(leave.dist)

          if statistics.mean(distance) < dist:
               for leave in node:
                    if leave.name not in (type_strains + sample_strains):
                         leave.delete()

     for node in tree.traverse("preorder"):
          for leave in node:
               if leave.name in type_strains:
                    # Leave is a type_strain
                    style1 = NodeStyle()
                    style1["fgcolor"] = "#000000"
                    style1["shape"] = "circle"
                    style1["size"] = 10
                    style1["vt_line_type"] = 1 # dashed
                    leave.img_style = style1
               elif leave.name in sample_strains:
                    # Leave is a sample_strain
                    style2 = NodeStyle()
                    style2["fgcolor"] = "#14e05c"
                    style2["shape"] = "circle"
                    style2["size"] = 15
                    style2["vt_line_type"] = 2 # dotted
                    leave.img_style = style2
               else:
                    # Leave is no sample_strain or type_strain
                    style3 = NodeStyle()
                    style3["fgcolor"] = "darkred"
                    style3["shape"] = "sphere"
                    style3["size"] = 5
                    style3["vt_line_type"] = 0 # solid
                    leave.img_style = style3

if __name__ == '__main__':
    args = parse_arguments()
    # Load a tree structure from a newick file.
    t = Tree(args.input_file)
    # Pruning the tree
    prune_tree_by_distance(t,
                           args.distance,
                           args.sample_strains,
                           args.type_strains)

    # Saving pruned tree to output
    t.render(args.output_name, tree_style=ts, w=200, units="mm")
