#!/usr/bin/python

import argparse as args
import re
import sys

import adjacency_creation
import build_cBP

__author__ = 'klamkiewicz'

def validate_input(first_genome, second_genome):
    '''
    At the moment no duplications or indel events are allowed. Therefore the set of the genome content
    has to be identical. This is checked here.
    :param first_genome: Sequence of the first genome
    :param second_genome:  Sequence of the second genome
    :return: Boolean variable
    '''

    # remove signs and chromosomes
    first_genome = set(re.sub('[-)|]', '', first_genome))
    second_genome = set(re.sub('[-)|]', '', second_genome))

    # symmetrical difference. Just take elements that are unique in one
    # of the sets
    difference = first_genome ^ second_genome

    return (len(difference) == 0, difference)


parser = args.ArgumentParser(description="Enter two genomes in the input format of Unimog")
parser.add_argument('genomeA', metavar='genomeA', type=str, help="Sequence of the first genome")
parser.add_argument('genomeB', metavar='genomeB', type=str, help="Sequence of the second genome")

arguments = parser.parse_args()

same_content = validate_input(arguments.genomeA, arguments.genomeB)

if not same_content[0]:
    print >> sys.stderr, "Content of the genomes is not equal. Check your input!"
    print >> sys.stderr, " ".join(["%s" % x for x in same_content[1]])
    sys.exit(1)

adjacency_setA = adjacency_creation.create_adjacency_set(arguments.genomeA)
adjacency_setB = adjacency_creation.create_adjacency_set(arguments.genomeB)

build_cBP.connect_adjacencies(adjacency_setA, adjacency_setB)