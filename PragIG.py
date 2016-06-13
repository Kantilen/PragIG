#!/usr/bin/python

import build_adjacencies
import build_cBP
import argparse as args
import sys


__author__ = 'klamkiewicz'

parser = args.ArgumentParser(description="Enter two genomes in the input format of Unimog")
parser.add_argument('genomeA', metavar='genomeA', type=str, help="Sequence of the first genome")
parser.add_argument('genomeB', metavar='genomeB', type=str, help="Sequence of the second genome")

arguments = parser.parse_args()

same_content = build_adjacencies.validate_input(arguments.genomeA, arguments.genomeB)

if not same_content[0]:
    print >> sys.stderr, "Content of the genomes is not equal. Check your input!"
    print >> sys.stderr, " ".join(["%s" % x for x in same_content[1]])
    sys.exit(1)

adjacency_setA = build_adjacencies.create_adjacency_set(arguments.genomeA)
adjacency_setB = build_adjacencies.create_adjacency_set(arguments.genomeB)

build_cBP.connect_adjacenies(adjacency_setA, adjacency_setB)