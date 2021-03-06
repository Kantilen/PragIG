#!/usr/bin/python

#################################
# Import section                #
#################################
import sys
import os
import copy

from Bio import Phylo
#################################

__author__ = 'klamkiewicz'


class Input:
    '''
    We use this class for reading and parsing the input files. There is also a method that
    returns all genomes that are leaves and children of the same ancestor. Each pair of such genomes are
    returned in a list.
    DEPENDENCY: BioPython
    '''
    def __init__(self, genomes, tree, rename):
        '''
        Validation of file paths and reading of the content
        :param genomes: Path to the genome content file
        :param tree: Path to NEWICK tree
        '''
        if not (os.path.isfile(genomes)):
            print >> sys.stderr, "Error: Couldn't find genome file. Check your input.\n%s" % (genomes)
            sys.exit(1)

        if not (os.path.isfile(tree)):
            print >> sys.stderr, "Error: Couldn't find tree file. Check your input.\n%s" % (tree)
            sys.exit(1)

        self.genomes = self.read_genomes(genomes)
        self.tree = self.read_tree(tree, rename)


    def read_tree(self,tree, rename):
        '''
        Reads the NEWICK tree! DEPENDENCY: BioPython
        :param tree: Path to the NEWICK tree
        :param rename: boolean variable, that indicates whether the internal nodes of the tree are renamed
        based on their children
        :return: Some BioPython Object containing information of the tree.
        '''
        newick_tree = Phylo.read(tree, 'newick')

        # if rename is True, each ancestral node gets a new label
        # For this, the two children labels are joined together
        if rename:
            copied_tree = copy.deepcopy(newick_tree)
            for node in copied_tree.get_nonterminals(order="postorder"):
                leaves = node.find_clades()
                leaves = [x for x in leaves if x.name and x != node]

                original_leaves = []
                for leaf in leaves:
                    original_leaf = newick_tree.find_clades(leaf.name)
                    original_leaves.extend([x for x in original_leaf])

                original_node = newick_tree.common_ancestor(original_leaves)
                original_node.name = "".join([x.name for x in leaves if x.name])

                node.name = "".join([x.name for x in leaves if x.name])

                for leaf in leaves:
                    copied_tree.collapse(leaf)

        return (newick_tree, max(newick_tree.depths().values()))


    def read_genomes(self, genomes):
        '''
        This function reads the genome content.
        At the moment only the GRIMM format is supported.
        :param genomes: Path to the genome content file
        :return: dictionary of the genome name as key and content as value
        '''
        genome_content = {}
        with open(genomes,'r') as gene_content:
            value = []
            key = ''
            while 1:
                current_genome = gene_content.readline()
                if not current_genome:
                    genome_content.update({key:value})
                    break
                if current_genome.startswith('#'):
                    continue
                if not current_genome.startswith('>'):
                    value.extend(current_genome.rstrip('\n ').split(' '))
                else:
                    if not key:
                        key = current_genome.strip('\n>')
                        continue
                    genome_content.update({key:value})
                    key = current_genome.strip('\n>')
                    value = []

        return genome_content

    def find_pairwise_leaves(self, tree):
        '''
        Traverses a given tree and returns all pairwise genomes (same ancestor) as a list.
        :param tree:
        :return: list of tuples. Each tuple represents two siblings in the tree.
        '''
        pairwise_genomes = []

        for clade in tree.find_clades():
            if clade.count_terminals() == 2:
                leaves = clade.find_clades()
                leaves = [x.name for x in leaves if x.name and not x == clade]
                distances_to_ancestor = [clade.distance(leaves[0],clade), clade.distance(leaves[1],clade)]
                pairwise_genomes.append((leaves,distances_to_ancestor))
        return pairwise_genomes

    def find_all_leaves(self, tree):
        '''
        Returns a list of all leaves in tree.
        :param tree: NEWICK tree
        :return: list of leaves
        '''
        return tree.get_terminals()