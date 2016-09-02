#!/usr/bin/python
from __future__ import division
import sys
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