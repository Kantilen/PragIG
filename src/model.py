#!/usr/bin/python

__author__ = 'klamkiewicz'

#################################
# Import section                #
#################################
import re
import numpy as np
#################################

class Adjacency():
    def __init__(self, first_ex, second_ex):
        self.first_ex = first_ex
        self.second_ex = second_ex

    def __repr__(self):
        return "(%s,%s)" % (self.first_ex, self.second_ex)

    def __eq__(self,other):
        if isinstance(other, Adjacency):
            return (self.first_ex == other.first_ex and self.second_ex == other.second_ex) or \
               (self.first_ex == other.second_ex and self.second_ex == other.first_ex)
        else:
            return False

    def __hash__(self):
        return hash(self.__repr__())

    def is_telomere(self):
        return (self.second_ex == None)

    def is_in_list(self, adj_list):
        return (self in adj_list)

    def contains_extremity(self, ext):
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

    def create_binary_vector(self, inter_adjacency):
        binary = np.zeros([1,len(inter_adjacency)], dtype=int)
        for adjacency in self.adjacency_set:
            try:
                np.put(binary,1,inter_adjacency.index(adjacency))
            except ValueError:
                continue

    def length(self):
        return len(self.content)

    def adj_length(self):
        return len(self.adjacency_set)

    def contains(self,adjacency):
        return (adjacency in self.adjacency_set)