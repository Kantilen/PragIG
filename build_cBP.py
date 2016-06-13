#!/usr/bin/python

import re

__author__ = 'klamkiewicz'

visited_vertices_in_a = []
visited_vertices_in_b = []


def connect_adjacenies(adjA, adjB):

    for adjacency in adjA:
        single_extremity = re.search('([\d*\w*]+[t|h])([\d*\w*]+[t|h])', adjacency)
        if not single_extremity:
            print adjacency,
        else:
            print single_extremity.group(1), single_extremity.group(2),
        #print single_extremity.group(0)

    print

    for adjacency in adjB:
        single_extremity = re.search('([\d*\w*]+[t|h])([\d*\w*]+[t|h])', adjacency)
        if not single_extremity:
            print adjacency,
        else:
            print single_extremity.group(1), single_extremity.group(2),