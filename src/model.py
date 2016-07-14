#!/usr/bin/python

__author__ = 'klamkiewicz'

#################################
# Import section                #
#################################
import re
import numpy as np
import sys
#################################

class Adjacency():
    def __init__(self, first_ex, second_ex):
        self.first_ex = first_ex
        self.second_ex = second_ex

    def __repr__(self):
        return "%s%s" % (self.first_ex, self.second_ex)

    def __str__(self):
        return "(%s,%s)" % (self.first_ex, self.second_ex)

    def __eq__(self,other):
        if isinstance(other, Adjacency):
            return (self.first_ex == other.first_ex and self.second_ex == other.second_ex) or \
               (self.first_ex == other.second_ex and self.second_ex == other.first_ex)
        else:
            return False

    def __hash__(self):
        return hash(self.__repr__())

    def get_extremities(self):
        return [self.first_ex, self.second_ex]

    def is_telomere(self):
        return (self.second_ex == None)

    def is_in_list(self, adj_list):
        return (self in adj_list)

    def contains_extremity(self, ext):
        if not ext:
            return False
        return (self.first_ex == ext or self.second_ex == ext)


class Genome():
    def __init__(self,name,content):
        self.name = name
        self.content = content
        self.adjacency_set = self.create_adjacency_set()


    def create_adjacency_set(self):
        '''
        This function reads the genome content and creates the adjacency set of it.
        Note that the genome should be represented in the UniMog input notation.
        The corresponding adjacency set is stored in the global variable.
        '''
        adjacencies = []  # returned value
        current_chromosome = []  # this is needed if there are more chromosomes

        # if not sign is given, linear chromosomes are default
        if not (self.content[-1] == ')' or self.content[-1] == '$'):
            self.content.append('$')

        # main iteration
        for index, gene in enumerate(self.content):
            if not (re.match('[\d*\w*]+', gene) or gene.startswith('-')):  # chromosome sign detected
                # the first_content and last extremity have to be considered separately
                first = current_chromosome[0]
                last = current_chromosome[-1]
                current_chromosome.remove(first)
                current_chromosome.remove(last)

                # connecting adjacencies
                current_chromosome = iter(current_chromosome)
                adjacencies.extend([Adjacency(ex,next(current_chromosome)) for ex in current_chromosome])

                if gene == '$':  # linear chromosome needs telomeres
                    adjacencies.insert(0, Adjacency(first,None))
                    adjacencies.append(Adjacency(last,None))
                if gene == ')':  # circular chromosomes needs another adjacency
                    adjacencies.append(Adjacency(first,last))

                current_chromosome = []  # clear the chromosome
                continue

                # orientation of the gene is considered here
            if gene.startswith('-'):  # negative sign is ignored, but extremities are switched
                current_chromosome.append('%sh' % gene[1:])
                current_chromosome.append('%st' % gene[1:])
            else:
                current_chromosome.append('%st' % gene)
                current_chromosome.append('%sh' % gene)

        return adjacencies

    def create_binary_vector(self, inter_adj):
        '''
        This function creates an numpy array that represents the binary vector of
        adjacencies in the two genomes that have to compared.
        '''
        preprocessing_dict = dict.fromkeys(self.adjacency_set)  # This makes O(k+n) instead of O(k*n)!
        binary = np.zeros([1, len(inter_adj)], dtype=int)

        for index, int_adj in enumerate(inter_adj):
            if int_adj in preprocessing_dict:
                np.put(binary, index, 1)
        return binary

    def length(self):
        return len(self.content)

    def adj_length(self):
        return len(self.adjacency_set)

    def contains(self,adjacency):
        return (adjacency in self.adjacency_set)


    @staticmethod
    def genome_from_adjacencies(name, adjacency_set):
        extremities = {'h': 't', 't':'h'}
        adjacencies = {}
        telomere = []
        content = []

        for adjacency in adjacency_set:
            if adjacency.first_ex in adjacencies.keys():
                print adjacency, adjacencies[adjacency.first_ex]
            if adjacency.second_ex in adjacencies.keys():
                print adjacency, adjacencies[adjacency.second_ex]
            if adjacency.second_ex.startswith('Telo'):
                telomere.append(adjacency.first_ex)
                continue
            if adjacency.first_ex.startswith('Telo'):
                telomere.append(adjacency.second_ex)
                continue
            adjacencies.update({adjacency.first_ex : adjacency.second_ex, adjacency.second_ex : adjacency.first_ex})

        current_gene = telomere[0]

        while adjacencies:

            try:
                adjacencies.pop(current_gene)
            except KeyError:
                pass

            single_extremity = re.search('([\d*]+)([t|h])', current_gene)
            orientation = single_extremity.group(2)
            gene = single_extremity.group(1)
            if orientation == 'h':
                content.append('-%s' % (gene))
            else:
                content.append('%s' % (gene))
            other_extremity = "%s%s" % (gene,extremities[orientation])

            try:
                current_gene = adjacencies[other_extremity]
                adjacencies.pop(other_extremity)
            except KeyError:
                break

        return Genome(name, content)