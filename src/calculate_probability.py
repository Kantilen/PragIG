#!/usr/bin/python
from __future__ import division

__author__ = 'klamkiewicz'

#################################
# Import section                #
#################################
import math

import networkx as nx
#################################

def optimal_scenarios(graph):
    '''
    Calculates the number of optimal sorting scenarios between two genomes and their
    circular breakpoint graph.
    Equation is taken from Braga and Stoye, 2010.
    Logarithm is used to deal with huge numbers
    :param graph: circular breakpoint graph of two input genomes
    :return: number of optimal sorting scenarios
    '''
    distances = {}
    for component in nx.connected_component_subgraphs(graph):
        # cycles of length 2 are already sorted
        if len(component) == 2:
            continue
        # otherwise the distance for specific cycle is saved
        distances[component] = int((len(component)/2) - 1)
    upper = 0
    lower = 0
    prod = 1
    for distance in distances.values():
        upper += distance
        lower += math.log10(math.factorial(distance))
        for i in range(distance-1):
            prod += math.log10(distance+1)

    result = math.log10(math.factorial(upper)) - lower + prod
    return result

def all_scenarios(length, distance):
    '''
    Calculates the number of all sorting scenarios of a given genome length and distance.
    For each step, n*(n-1) possible operation exist, there are distance steps in total
    :param length: number of genes in the genomes
    :param distance: number of steps in the sorting scenarios.
    :return: number of all possible sorting scenarios.
    '''
    result = 0
    for i in range(int(distance)):
        result += (math.log10(length) + math.log10(length-1))
    return result

def calculate_prob_ancestor(probability, sorting_scen, all_scen, alpha, distance, tree_distance, expected_distance, lower_bound, upper_bound):
    '''
    Calculates the probability that a sampled genome is the ancestor of its children.
    :param probability: old probability that was calculated in an earlier step
    :param sorting_scen: number of optimal sorting scenarios between two genomes
    :param all_scen: number of all possible sorting scenarios of a given length
    :param alpha: alpha parameter from PragIG. Determines the strictness of the filter
    :param distance: DCJ distance between the two genomes
    :param tree_distance: the DCJ distance implied by the tree between the two genomes
    :param expected_distance: the expected DCJ distance calculated as porposed in Biller et al., 2015
    :param lower_bound: lower bound for the filter step
    :param upper_bound: upper bound for the filter step
    :return: the new probability for the sampled genome
    '''
    prob_sampled_genomes = sorting_scen - all_scen

    if alpha == 1.0 or distance <= tree_distance <= expected_distance:
        probability += prob_sampled_genomes
    else:
        if tree_distance <= distance:
            probability += prob_sampled_genomes + math.log10(tree_distance - lower_bound) - \
                           math.log10(distance - lower_bound)
        else:
            probability += prob_sampled_genomes + math.log10(upper_bound - tree_distance) - \
                           math.log10(expected_distance * (1 - alpha))

    return probability