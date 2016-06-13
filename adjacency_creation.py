#!/usr/bin/python

#################################
# Import section                #
#################################
import re
#################################

__author__ = 'klamkiewicz'

def create_adjacency_set(genome_content):
    '''
    This function reads the genome content and creates the adjacency set of it.
    Note that the genome should be represented in the UniMog input notation.
    :param genome_content: UniMog representation of the genome. Type=String
    :return: Adjacency set of the genome. Type=List (of strings)
    '''
    adjacencies = [] # returned value
    current_chromosome = [] # this is needed if there are more chromosomes

    # some preprocessing. Force spaces at chromosome characters
    #TODO: Auslagern - den Kram brauch ich oefter, nicht nur hier
    genome_content = re.split('\s+' ,genome_content.replace('|', ' | ').replace(')', ' ) ').strip())

    # if not sign is given, linear chromosomes are default
    if not (genome_content[-1] == ')' or genome_content[-1] == '|'):
        genome_content.append('|')

    # main iteration
    for index, gene in enumerate(genome_content):
        if not (re.match('[\d*\w*]+', gene)  or gene.startswith('-') ): # chromosome sign detected

            # the first and last extremity have to be considered separately
            first = current_chromosome[0]
            last = current_chromosome[-1]
            current_chromosome.remove(first)
            current_chromosome.remove(last)

            # connecting adjacencies
            current_chromosome = iter(current_chromosome)
            adjacencies.extend([ex + next(current_chromosome) for ex in current_chromosome])

            if gene == '|': # linear chromosome needs telomeres
                adjacencies.insert(0,first)
                adjacencies.append(last)
            if gene == ')': # circular chromosomes needs another adjacency
                adjacencies.append('%s%s' % (last,first))

            current_chromosome = [] # clear the chromosome
            continue

        # orientation of the gene is considered here
        if gene.startswith('-'): # negative sign is ignored, but extremities are switched
            current_chromosome.append('%sh' % gene[1:])
            current_chromosome.append('%st' % gene[1:])
        else:
            current_chromosome.append('%st' % gene)
            current_chromosome.append('%sh' % gene)
    return adjacencies