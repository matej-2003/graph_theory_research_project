#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import math


def save_graph(G, n, a):
	# nx.write_graphml(G, f"./graphs/G_n={n}_a={a}.graphml")
	nx.write_graph6(G, f"./graphs/G_n={n}_a={a}.g6")

def save_graphs(graph_list):
	for G, n, a in graph_list:
		save_graph(G, n, a)

def display_graphs(graphs, figsize=(15, 15), show_titles=True, force_row=False, with_labels=False, pos=None):
	def axes_division(x):
		factor1 = 1
		for i in range(1, int(math.sqrt(x)) + 1):
			if x % i == 0:
				factor1 = i
		return int(factor1), int(x / factor1)

	# layout override
	if force_row:
		rows, cols = 1, len(graphs)
	else:
		rows, cols = axes_division(len(graphs))

	fig, axes = plt.subplots(rows, cols, figsize=figsize)

	# ensure axes is flat array
	axes = np.array(axes).flatten()

	for (G, n, a), ax in zip(graphs, axes):
		if pos == None:
			pos_ = nx.spring_layout(G, seed=100)
		else:
			pos_ = pos(G)
		nx.draw(
			G,
			pos_,
			ax=ax,
			with_labels=with_labels,
			node_size=20,
			font_size=10,
			node_color=['red' if G.nodes[v].get("infected") else 'lightgray'
						for v in G],
		)

		if show_titles:
			ax.set_title(f"Graph C({n}, {a})")
		else:
			ax.set_title("")

		ax.set_aspect("equal")

	plt.tight_layout()
	plt.show()

def display_graph(G, figsize=(6, 6), show_title=True, with_labels=False, pos=None, node_size=60, font_size=10):
	fig, ax = plt.subplots(figsize=figsize)

	if pos == None:
		pos = nx.spring_layout(G, seed=100)

	nx.draw(
		G,
		pos,
		ax=ax,
		with_labels=with_labels,
		node_size=node_size,
		font_size=font_size,
		node_color=[
			'red' if G.nodes[v].get("infected") else 'lightgray'
			for v in G
		],
	)

	if show_title:
		ax.set_title(f"Graph C({len(G.nodes)})")
	else:
		ax.set_title("")

	ax.set_aspect("equal")
	plt.tight_layout()
	plt.show()

def generate_graphs(n):
	n_half = int(n/2)
	graphs = []
	
	def add_candidate_graph(H, a):
		nonlocal graphs
		for (G_i, _, _) in graphs:
			if nx.is_isomorphic(G_i, H):
				return
		graphs.append((H, n, a))
	
	for a in range(n_half):
		gcd = np.gcd(a, n_half)
		condition = (1 <= a < n_half and gcd == 1)
		
		if condition:
			g = nx.circulant_graph(n, [a, -a, n_half])
			add_candidate_graph(g, a)
			# print(f"a={a} n={n}")
	
	print(f"number of graphs for n={n} :", len(graphs))
	return graphs

def is_infectable(G, node):
	node_N = list(nx.all_neighbors(G, node))
	if len(node_N) < 2:
		return False
	
	infected_N_counter = sum(1 for i_neighbor in node_N 
							 if G.nodes[i_neighbor].get("infected"))
	return infected_N_counter >= 2

def fully_endemic(G_):
	G = G_.copy()
	max_itter_number = len(G.nodes)
	infected_nodes = [n for n in G.nodes if G.nodes[n].get("infected")]
	healthy_nodes = [n for n in G.nodes if not G.nodes[n].get("infected")]

	itter_number = 0
	while healthy_nodes:
		# find all new infections in THIS wave
		newly_infected = [node for node in healthy_nodes if is_infectable(G, node)]

		if not newly_infected:
			break

		# infect all of them at once
		for node in newly_infected:
			G.nodes[node]["infected"] = True

		# update tracking sets
		for node in newly_infected:
			infected_nodes.append(node)
			healthy_nodes.remove(node)

		itter_number += 1
		if itter_number > max_itter_number:
			break

	return (G, infected_nodes, healthy_nodes)

def infection_spread(G):
	graph_stages = [G.copy()]
	infected_nodes = [u for u in G.nodes if G.nodes[u].get("infected")]
	healthy_nodes = [u for u in G.nodes if not G.nodes[u].get("infected")]

	itter_number = 0
	max_iters = len(G.nodes)

	while healthy_nodes:
		# determine whole wave
		newly_infected = [node for node in healthy_nodes if is_infectable(G, node)]

		if not newly_infected:
			break

		# infect all of them simultaneously
		for node in newly_infected:
			G.nodes[node]["infected"] = True

		# update node lists
		for node in newly_infected:
			infected_nodes.append(node)
			healthy_nodes.remove(node)

		# add stage AFTER wave completes
		graph_stages.append(G.copy())

		itter_number += 1
		if itter_number > max_iters:
			break

	return [H for H in graph_stages], infected_nodes, healthy_nodes

def show_graph_infection(H):
	print("C=", end="")
	for j in H.nodes:
		if H.nodes[j].get("infected"):
			print(1, end="")
		else:
			print(0, end="")
	print()

def infection_graph_permuation(vertex_list, infection_number, z=0):
	# Base case: no more infections to place
	if infection_number == 0:
		yield vertex_list
		return

	# Try placing the next infection at any allowed index
	for i in range(z, len(vertex_list)):
		if vertex_list[i] == 0:
			new_list = vertex_list[:]  # copy
			new_list[i] = 1
			yield from infection_graph_permuation(new_list, infection_number - 1, i + 1)

def infect_graph(G, infection_mask):
	H = G.copy()
	for i, e in enumerate(infection_mask):
		H.nodes[i]["infected"] = bool(e)
	return H

def generate_graphs2(n):
	n_half = int(n/2)
	graphs = []
	
	def add_candidate_graph(H, a):
		nonlocal graphs
		for (G_i, n_i, a_i) in graphs:
			if nx.is_isomorphic(G_i, H):
				print(f"n={n_i} a={a_i} is isomorpic")
				return
		graphs.append((H, n, a))
	
	for a in range(n_half):
		gcd = np.gcd(a, n_half)
		condition = (1 <= a < n_half and gcd == 1)
		
		if condition:
			g = nx.circulant_graph(n, [a, -a, n_half])
			add_candidate_graph(g, a)
	
	print(f"number of graphs for n={n} :", len(graphs))
	return graphs

def get_cycles(G):
    even_nodes = [v for v in G.nodes if v % 2 == 0]
    odd_nodes  = [v for v in G.nodes if v % 2 == 1]
    
    def cycle_subgraph(G, nodes):
        n = len(nodes)
        edges = [
            (nodes[i], nodes[(i+1) % n])
            for i in range(n)
        ]
        # print(edges)
        return G.edge_subgraph(edges)
    
    
    C_top    = cycle_subgraph(G, even_nodes)
    C_bottom = cycle_subgraph(G, odd_nodes)

    return C_top, C_bottom 

def ordered_cycle_nodes(C):
    nodes = list(C.nodes)
    start = nodes[0]

    order = [start]
    prev = None
    curr = start

    while True:
        nbrs = list(C.neighbors(curr))
        nxt = nbrs[0] if nbrs[0] != prev else nbrs[1]
        if nxt == start:
            break
        order.append(nxt)
        prev, curr = curr, nxt

    return order

def prism_pos(G, r_inner=1.0, r_outer=2.0, shift_top=None):
    C_top, C_bottom = get_cycles(G)

    top_nodes = ordered_cycle_nodes(C_top)
    bottom_nodes = ordered_cycle_nodes(C_bottom)

    if len(top_nodes) != len(bottom_nodes):
        raise ValueError("Top and bottom cycles must have same length")

    n = len(top_nodes)

    # Apply shift to top cycle
    if shift_top is None:
        shift_top = n // 2  # default n/2 counterclockwise

    top_nodes = top_nodes[shift_top:] + top_nodes[:shift_top]

    pos = {}

    for i in range(n):
        angle = 2 * math.pi * i / n

        pos[top_nodes[i]] = (
            r_inner * math.cos(angle),
            r_inner * math.sin(angle)
        )

        pos[bottom_nodes[i]] = (
            r_outer * math.cos(angle),
            r_outer * math.sin(angle)
        )

    return pos





def mobius_pos(G, radius=1.0):
    pos = {}

    nodes = list(G.nodes)
    n = len(nodes)
    angle_off = 2 * math.pi / (n + 1)

    for i, node in enumerate(nodes):
        angle = angle_off * i
        pos[node] = (
            radius * math.cos(angle),
            radius * math.sin(angle)
        )

    return pos




def color_distance_partitions(G, v = 0, figsize=(3, 3), with_labels=False, pos=None):
    # get BFS layers
    layers = list(nx.bfs_layers(G, v))
    
    # assign a color index to each node based on distance
    node_to_layer = {}
    for i, layer in enumerate(layers):
        for node in layer:
            node_to_layer[node] = i
    
    
    print(layers)
    # layout (spring works fine, radial is optional)
    if pos == None:
        pos_ = nx.spring_layout(G, seed=100)
    else:
        pos_ = pos(G)
    
    # colors by layer
    node_colors = [node_to_layer[n] for n in G.nodes()]
    
    plt.figure(figsize=figsize)
    cmap = plt.cm.viridis
    max_layer = max(node_colors) if node_colors else 0
    norm = plt.Normalize(vmin=0, vmax=max_layer if max_layer > 0 else 1)

    edge_colors = [
        'red' if G.nodes[n].get("infected", False) else cmap(norm(node_colors[i]))
        for i, n in enumerate(G.nodes())
    ]
    nx.draw(
        G,
        pos_,
        with_labels=with_labels,
        node_color=node_colors,
        edgecolors=edge_colors,
        linewidths=2,
        cmap=plt.cm.viridis,
        node_size=150,
        font_size=8,
    )
    
    plt.title(f"Distance partitions from v = {v}")
    plt.show()


def visualize_distance_partitions(G, v = 0, figsize=(8, 4), with_labels=False):
    # compute distance layers
    layers = list(nx.bfs_layers(G, v))
    
    # build custom positions
    pos = {}
    x_spacing = 2.0
    y_spacing = 1.0
    
    for i, layer in enumerate(layers):
        x = i * x_spacing
        y_start = -(len(layer) - 1) * y_spacing / 2
    
        for j, node in enumerate(sorted(layer)):
            y = y_start + j * y_spacing
            pos[node] = (x, y)
    
    # color by layer
    node_colors = []
    for node in G.nodes():
        for i, layer in enumerate(layers):
            if node in layer:
                node_colors.append(i)
                break
    print("node_colors", node_colors)
    cmap = plt.cm.viridis
    max_layer = max(node_colors) if node_colors else 0
    norm = plt.Normalize(vmin=0, vmax=max_layer if max_layer > 0 else 1)

    edge_colors = [
        'red' if G.nodes[n].get("infected", False) else cmap(norm(node_colors[i]))
        for i, n in enumerate(G.nodes())
    ]
    
    plt.figure(figsize=figsize)
    nx.draw(
        G,
        pos,
        with_labels=with_labels,
        node_color=node_colors,
        edgecolors=edge_colors,
        linewidths=2,
        cmap=plt.cm.viridis,
        node_size=150,
        font_size=8,
    )
    
    plt.title(f"Distance-layer visualization from v = {v}")
    plt.axis("off")
    plt.show()




def distance_partitions_pos(G, v=0):
    layers = list(nx.bfs_layers(G, v))

    pos = {}
    x_spacing = 2.0
    y_spacing = 1.0

    for i, layer in enumerate(layers):
        x = i * x_spacing
        y_start = -(len(layer) - 1) * y_spacing / 2

        for j, node in enumerate(sorted(layer)):
            pos[node] = (x, y_start + j * y_spacing)

    return pos

def display_graph2(
    G,
    figsize=(6, 6),
    show_title=True,
    with_labels=False,
    pos=None,
    node_size=60,
    font_size=10,
    labels=None,
):
    fig, ax = plt.subplots(figsize=figsize)

    if pos is None:
        pos = nx.spring_layout(G, seed=100)

    nx.draw(
        G,
        pos,
        ax=ax,
        with_labels=False,  # always false, labels handled separately
        node_size=node_size,
        font_size=font_size,
        node_color=[
            'red' if G.nodes[n].get("infected", False) else 'lightgray'
            for n in G.nodes()
        ],
    )

    # draw labels if provided
    if labels is not None:
        nx.draw_networkx_labels(
            G,
            pos,
            labels=labels,
            font_size=font_size,
            ax=ax,
        )
    elif with_labels:
        # fallback: default node labels
        nx.draw_networkx_labels(
            G,
            pos,
            labels={n: n for n in G.nodes()},
            font_size=font_size,
            ax=ax,
        )

    if show_title:
        ax.set_title(f"Graph C({len(G.nodes)})")
    else:
        ax.set_title("")

    ax.set_aspect("equal")
    plt.tight_layout()
    plt.show()





