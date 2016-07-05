#!/usr/bin/python

#################################
# Import section                #
#################################
import argparse as args
import sys
import os
#################################

__author__ = 'klamkiewicz'


class Input:

    def __init__(self, genomes, tree):

        if not (os.path.isfile(genomes)):
            print >> sys.stderr, "Error: Couldn't find genome file. Check your input.\n%s" % (genomes)
            sys.exit(1)

        if not (os.path.isfile(tree)):
            print >> sys.stderr, "Error: Couldn't find tree file. Check your input.\n%s" % (tree)
            sys.exit(1)

        self.genomes = self.read_genomes(genomes)
        self.tree = self.read_tree(tree)


    def read_tree(self,tree):

        num_lines = sum(1 for line in open(tree))
        if num_lines != 1:
            print "Error: Your tree file has more than one line.\nExiting..."
            sys.exit(1)

        with open(tree, 'r') as tree_content:
            newick_tree = tree_content.readline().rstrip("\n")
        return newick_tree


    def read_genomes(self, genomes):
        genome_content = {}
        with open(genomes,'r') as gene_content:
            value = []
            key = ''
            while 1:
                current_genome = gene_content.readline()
                if not current_genome:
                    value.append(')')
                    genome_content.update({key:value})
                    break

                if current_genome.startswith('#'):
                    continue

                if not current_genome.startswith('>'):
                    value.extend(current_genome.rstrip('$\n ').split(' '))

                else:
                    if not key:
                        key = current_genome.strip('\n>')
                        continue
                    value.append(')')
                    genome_content.update({key:value})
                    key = current_genome.strip('\n>')
                    value = []

        return genome_content



    #def read_input(self, genomes, tree):


    #    newick_tree = read_tree(tree)
    #    genome_content = read_genomes(genomes)

    #    return genome_content, newick_tree