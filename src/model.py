#!/usr/bin/python

__author__ = 'klamkiewicz'

#################################
# Import section                #
#################################
import re
#################################

class Adjacency():
    def __init__(self, first_ex, second_ex):
        self.first_ex = first_ex
        self.second_ex = second_ex

    def __str__(self):
        return "%s%s" % (self.first_ex, self.second_ex)

    def equal(self,adjacency):
        return (self.first_ex == adjacency.first_ex and self.second_ex == adjacency.second_ex) or \
               (self.first_ex == adjacency.second_ex and self.second_ex == adjacency.first_ex)

    def is_telomere(self):
        return (self.first_ex == None) or (self.second_ex == None)


class Genome():
    def __init__(self,name,content):
        self.name = name
        self.content = content
        self.adjacency_set = set()

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

        self.adjacency_set = adjacencies

    def length(self):
        return len(self.content)

    def adj_length(self):
        return len(self.adjacency_set)