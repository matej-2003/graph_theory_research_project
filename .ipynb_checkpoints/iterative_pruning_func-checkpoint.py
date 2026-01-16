import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import math
from graph_functions import *

def fully_infect(G):
    for i in range(len(G.nodes)):
        G.nodes[i]["infected"] = 1

def minimise_infection(G, v=0, start_layer=0):
    layers = list(nx.bfs_layers(G, v))
    # print(layers)
    
    for i in range(start_layer, len(layers), 2):
        layer = layers[i]
        # print(layer)
    
        for x in layer:
            G.nodes[x]["infected"] = 0

def calc_graph_values(G):
    for i, e in enumerate(G):
        if not G.nodes[i].get("infected"):
            infected_count = 0
    
            for y in G.neighbors(e):
                if G.nodes[y].get("infected"):
                    infected_count += 1
    
            G.nodes[i]["infection_count"] = infected_count
    
        # print(e, G.nodes[i])
    
    
    for i, e in enumerate(G.nodes):
        if G.nodes[i].get("infected"):
            c = 0
            for y in G.neighbors(e):
                if G.nodes[y].get("infection_count"):
                    c += 1
    
            G.nodes[i]["infection_influence"] = c
    
        # print(e, G.nodes[i])

    potential_candidates = []
    for i, e in enumerate(G.nodes):
        if G.nodes[i].get("infected") and isinstance(G.nodes[y].get("infection_count"), int) and G.nodes[y].get("infection_count") > 2:
            potential_candidates.append(e)

    # print(potential_candidates)
    return potential_candidates



def display_graph_stage(G, prism=True, figsize=(6, 3), v=0):
    labels = {
        n: (G.nodes[n].get("infection_influence") if G.nodes[n].get("infected") else G.nodes[n].get("infection_count"))
        for n in G.nodes()
    }

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    if prism:
        display_graph2(G, show_title=True, labels=labels, pos=prism_pos(G), ax=axes[0])
    else:
        display_graph2(G, show_title=True, labels=labels, pos=mobius_pos(G), ax=axes[0])
        
    # display_graph2(G, show_title=True, labels=labels, pos=distance_partitions_pos(G, v=v, x_spacing=3.0, y_spacing=3.0), ax=axes[1])
    display_graph2(G, show_title=True, labels=labels, pos=nx.bfs_layout(G, v), ax=axes[1])

    plt.tight_layout()
    plt.show()


def display_graphs_row(*graphs):
    fig, axes = plt.subplots(1, len(graphs), figsize=(6, 3))
    for i, G in enumerate(graphs):
        display_graph2(G, show_title=False, ax=axes[i])

    plt.tight_layout()
    plt.show()


def display_graph_distance_par(G, pos=None, v=0, figsize=(8, 3)):
    labels = {
        n: ("{}, {}".format(G.nodes[n].get("infected_n"), G.nodes[n].get("healthy_n")))
        for n in G.nodes()
    }
    
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    display_graph2(G, labels=labels, pos=pos, ax=axes[0])
    # display_graph2(G, labels=labels, pos=distance_partitions_pos(G, v=v, x_spacing=6.0, y_spacing=3.0), ax=axes[1],)
    display_graph2(G, labels=labels, pos=nx.bfs_layout(G, v)) #, x_spacing=6.0, y_spacing=3.0), ax=axes[1],)

    plt.tight_layout()
    plt.show()


def calc_values(G):
    for i in G.nodes:
        e = G.nodes[i]
        adj = G.adj[i]

        infected_n = 0
        for x in adj:
            z =  G.nodes[x]
            if z.get("infected"):
                infected_n += 1

        G.nodes[i]["infected_n"] = infected_n
        G.nodes[i]["healthy_n"] = G.degree(i) - infected_n
