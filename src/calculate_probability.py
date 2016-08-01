#!/usr/bin/python
from __future__ import division
import sys
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
    prob_cut_adj = float(2/size)
    # 1 -> 1
    prob_keep_adj = 1 - prob_cut_adj

    probabilities = [{0:0, 1:prob_form_adj}, {0:1, 1:prob_avoid_adj}, {0:0, 1:prob_cut_adj}, {0:1, 1:prob_keep_adj}]
    first_transition_probabilities = dict(zip(identifier,probabilities))

    transition_length.update(dynamic_table(length, first_transition_probabilities))
    #for key, value in transition_length.items():
    #    print key
    #    print value
    #sys.exit(0)

def dynamic_table(length, transitions):
    for index in range(2,length+1):
        for event in transitions.keys():
            i,j = event
            transitions[event][index] = max(transitions[(i,0)][index-1]*transitions[(0,j)][1],
                                            transitions[(i,1)][index-1]*transitions[(1,j)][1])
    return transitions


def calculate_probability(first, second, ancestral, distances):
    prob = 1.0
    for index, element in enumerate(ancestral):
        prob *= transition_length[(first[index],element)][distances[0]]*transition_length[second[index],element][distances[1]]
    return prob