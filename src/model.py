#!/usr/bin/python
from __future__ import division

__author__ = 'klamkiewicz'

#################################
# Import section                #
#################################
import re
from collections import Counter
import math

import networkx as nx
#################################

class Adjacency():
    '''
    Class that models adjacencies. Two extremities are needed in order to create an adjacency.
    adj = Adjacency("1h","2t")
    '''
    def __init__(self, first_ex, second_ex):
        self.first_ex = first_ex
        self.second_ex = second_ex

    def __repr__(self):
        return "%s%s" % (self.first_ex, self.second_ex)

    def __str__(self):
        return "(%s,%s)" % (self.first_ex, self.second_ex)

    def __eq__(self,other):
        '''
        New function to evaluate equality between adjacencies.
        1h2t is the same adjacency as 2t1h. Therefore, the or connection in the if-statement is needed
        :param other: Adjacency that has to be compared against
        :return: Ture/False if the adjacencies are the same
        '''
        if isinstance(other, Adjacency):
            return (self.first_ex == other.first_ex and self.second_ex == other.second_ex) or \
               (self.first_ex == other.second_ex and self.second_ex == other.first_ex)
        else:
            return False

    def __hash__(self):
        '''
        Same as __eq__, in order to use the adjacencies as keys in a dictionary, the hash function
        has to be modified.
        :return:
        '''
        return hash(self.first_ex) ^ hash(self.second_ex)

    def get_extremities(self):
        '''
        Returns a list of the two extremities of the adjacency
        :return: list of extremities
        '''
        return [self.first_ex, self.second_ex]

    def is_telomere(self):
        '''
        Evaluates whether the adjacency is a telomere or not
        :return: True/False
        '''
        return (self.second_ex == None)

    def is_in_list(self, adj_list):
        '''
        Evaluates whether the adjacency is part of a component.
        :param adj_list: Cycle/Component/List of adjacencies
        :return: boolean value
        '''
        return (self in adj_list)

    def contains_extremity(self, ext):
        '''
        Evaluates whether the extremity ext is part of this adjacency
        :param ext: some extremity
        :return: boolean value
        '''
        if not ext:
            return False
        return (self.first_ex == ext or self.second_ex == ext)


class Genome():
    '''
    This class models Genomes.
    Only the GRIMM format is currently supported.
    '''
    def __init__(self,name,content):
        self.name = name
        self.content = content
        self.adjacency_set = self.create_adjacency_set()


    def create_adjacency_set(self):
        '''
        This function reads the genome content and creates the adjacency set of it.
        Note that the genome should be represented in the GRIMM notation.
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
                try:
                    first = current_chromosome[0]
                    last = current_chromosome[-1]
                    current_chromosome.remove(first)
                    current_chromosome.remove(last)
                except IndexError:
                    print index
                    print self.content[index-5:index+5]

                # connecting graph
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

    def expected_distance_to_genome(self, other_genome):
        '''
        Calculates the expected DCJ distance between the genome and a second genome other_genome.
        Equation is taken from Biller et al., 2015.
        :param other_genome: Genome Object that saves the second genome
        :return: expected DCJ distance of the two genomes
        '''
        breakpoints = sum(1 for adj in other_genome.adjacency_set if not self.contains(adj))
        observed = sum(1 for adj in self.adjacency_set if not adj.is_telomere())
        all_adj = self.adj_length()

        try:
            upper = math.log(1-((breakpoints*(2*all_adj -1) ) / (observed*(2*all_adj -2))))
            lower = math.log(1 - (1/(all_adj-1)) - 1/all_adj)
        except ValueError as error:
            return error
        distance = upper / lower
        return int(math.floor(distance))

    def length(self):
        '''
        Returns the number of genes of the genome
        :return: length of genome content
        '''
        return len(self.content)

    def adj_length(self):
        '''
        Returns the number of adjacencies in the genome
        :return: length of adjacency set
        '''
        return len(self.adjacency_set)

    def contains(self,adjacency):
        '''
        Determines whether an adjacency is present in the genome or not
        :param adjacency: The adjacency that has queried.
        :return: boolean variable
        '''
        return (adjacency in self.adjacency_set)

    def chr_number(self):
        '''
        Returns the number of chromosomes in the genome.
        :return: Number of chromosomes
        '''
        return Counter(self.content)['$'] + Counter(self.content)[')']

    def linear_chromosomes(self):
        '''
        Returns the number of LINEAR chromosomes in the genome.
        :return: number of linear chromosomes
        '''
        return Counter(self.content)['$']


    @staticmethod
    def genome_from_adjacencies(name, adjacency_set):
        '''
        Creates a new Genome() Object from a given adjacency set.
        For this, the genome content in GRIMM format is reconstructed from the adjacency set.
        :param name: name of the new genome
        :param adjacency_set: the adjacency set that represents the content of the genome
        :return: a new Genome() object
        '''
        extremities = {'h': 't', 't':'h'}
        adjacencies = {}
        telomere = []
        content = []

        for adjacency in adjacency_set:
            if adjacency.first_ex.startswith('Telo') and adjacency.second_ex.startswith('Telo'):
                continue
            if adjacency.second_ex.startswith('Telo'):
                telomere.append(adjacency.first_ex)
                adjacencies.update({adjacency.first_ex : None})
                continue
            if adjacency.first_ex.startswith('Telo'):
                telomere.append(adjacency.second_ex)
                adjacencies.update({adjacency.second_ex : None})
                continue
            adjacencies.update({adjacency.first_ex : adjacency.second_ex, adjacency.second_ex : adjacency.first_ex})

        current_adjacency = None

        while adjacencies:
            if not current_adjacency:
                if telomere:
                    current_adjacency = telomere[0]
                    telomere.remove(current_adjacency)
                else:
                    current_adjacency = adjacencies.keys()[0]

            try:
                adjacencies.pop(current_adjacency)
            except KeyError:
                content.append(')')
                current_adjacency = None
                continue

            single_extremity = re.search('([\d*]+)([t|h])', current_adjacency)
            gene = single_extremity.group(1)
            orientation = single_extremity.group(2)

            content.append('-%s' % gene) if orientation == 'h' else content.append('%s' % gene)

            other_extremity = "%s%s" % (gene, extremities[orientation])

            if other_extremity in telomere:
                telomere.remove(other_extremity)
                adjacencies.pop(other_extremity)
                content.append('$')
                current_adjacency = None
                continue

            current_adjacency = adjacencies[other_extremity]
            adjacencies.pop(other_extremity)

        return Genome(name, content)