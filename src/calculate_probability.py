#!/usr/bin/python
from __future__ import division
import sys
import networkx as nx
import math

from collections import defaultdict
__author__ = 'klamkiewicz'

transition_length = {}
form = (0,1)
avoid = (0,0)
cut = (1,0)
keep = (1,1)
identifier = [form, avoid, cut, keep]

def preprocess_transitions(length, size):

    # 0 -> 1
    prob_form_adj = (1 / (size * (size - 1)))
    # 0 -> 0
    prob_avoid_adj = 1 - (prob_form_adj)
    # 1 -> 0
    prob_cut_adj = (2/size)
    # 1 -> 1
    prob_keep_adj = 1 - prob_cut_adj

    probabilities = [{1: math.log10(prob_form_adj)}, {1: math.log10(prob_avoid_adj)}, {1: math.log10(prob_cut_adj)}, {1: math.log10(prob_keep_adj)}]

    first_transition_probabilities = dict(zip(identifier,probabilities))
    transition_length.update(dynamic_table(length, first_transition_probabilities))

def dynamic_table(length, transitions):

    for index in range(2,length+1):
        for event in transitions.keys():
            i,j = event
            transitions[event][index] = ((transitions[(i,0)][index-1]+transitions[(0,j)][1]) + (transitions[(i,1)][index-1]+transitions[(1,j)][1]))
    return transitions


def calculate_probability(binaries, ancestral, distances):
    prob = 0.0
    for index, element in enumerate(ancestral):
        for identifier in binaries.keys():
            prob += transition_length[(binaries[identifier][index],element)][distances[identifier]]
    return prob

def optimal_scenarios(graph):
    distances = {}
    for component in nx.connected_component_subgraphs(graph):
        if len(component) == 2:
            continue
        #TODO: Check what's going on here.... float?
        distances[component] = int((len(component)/2) - 1)
    result = 0
    upper = 0
    lower = 0
    prod = 1
    for distance in distances.values():
        upper += distance
        lower += math.log10(math.factorial(distance))
        for i in range(distance-1):
            prod += math.log10(distance+1)
        #prod *= (distance + 1) ** (distance - 1)

    result = math.log10(math.factorial(upper)) - lower + prod
    return result

def all_scenarios(length, distance):
    result = 0
    for i in range(int(distance)):
        result += (math.log10(length) + math.log10(length-1))
    return result

def calculate_prob_ancestor(probability, sorting_scen, all_scen, alpha, distance, tree_distance, expected_distance, lower_bound, upper_bound):
    prob_sampled_genomes = sorting_scen - all_scen
    # try:
    if alpha == 1.0 or distance <= tree_distance <= expected_distance:
        probability += prob_sampled_genomes
    else:
        if tree_distance <= distance:
            probability += prob_sampled_genomes + math.log10(tree_distance - lower_bound) - \
                           math.log10(distance - lower_bound)
        else:
            probability += prob_sampled_genomes + math.log10(upper_bound - tree_distance) - \
                           math.log10(expected_distance * (1 - alpha))
            # except ValueError:
            #    print >> sys.stderr, identifier, tree_distances[identifier], distance, arguments.alpha, distance*arguments.alpha
            #    sys.exit(0)

    return probability