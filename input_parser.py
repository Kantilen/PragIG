#!/usr/bin/python

#################################
# Import section                #
#################################
import sys
import os
from Bio import Phylo
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
        return Phylo.read(tree, 'newick')


    def read_genomes(self, genomes):
        genome_content = {}
        with open(genomes,'r') as gene_content:
            value = []
            key = ''
            while 1:
                current_genome = gene_content.readline()
                if not current_genome:
                    #value.append(')')
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
                    #value.append(')')
                    genome_content.update({key:value})
                    key = current_genome.strip('\n>')
                    value = []

        return genome_content

    def find_pairwise_leaves(self, tree):
        pairwise_genomes = []

        for clade in tree.find_clades():
            if clade.count_terminals() == 2:
                leaves = clade.find_clades()
                leaves = [x.name for x in leaves if x.name]
                pairwise_genomes.append(leaves)
        return pairwise_genomes