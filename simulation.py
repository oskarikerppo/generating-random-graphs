#!/usr/bin/python
# -*- coding: latin-1 -*-
#Simulation of Erdős–Rényi random network

#Oskari Kerppo
#9.5.2018

#matplotlib used for visualization


import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math
import random
import numpy as np


#https://gist.github.com/rougier/ebe734dcc6f4ff450abf
def binomial(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke.
    See http://stackoverflow.com/questions/3025162/statistics-combinations-in-python
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0




class network:
	"""Class for a network containing nodes and edges"""
	type = "Network"

	def __init__(self, name="Network"):
		self.name = name
		self.nodes = []
		self.edges = [] 

	def add_node(self, x):
		"""adds a node to network. If argument is list, all items are added"""
		if type(x) is list:
			for item in x:
				self.nodes.append(node(item))
		else:
			self.nodes.append(node(x))

	def list_nodes(self):
		"""prints and returns list of all nodes in network"""
		#print(self.nodes)
		return self.nodes

	def add_edge(self, x):
		"""Adds an edge between all tuples in x"""
		"""Warning !!! Tuples must be of type int"""
		all_nodes = self.list_nodes()
		if type(x) is tuple:
			if type(x[0]) is int and type(x[1]) is int:
				if x[0] < len(all_nodes) and x[1] < len(all_nodes):
					self.edges.append(x)
				else:
					print("Edge index out of range! Skipping this node")
			else:
				print("Invalid tuple type! Skipping this edge")
		elif type(x) is list:
			for item in x:
				if type(item) is tuple:
					if type(item[0]) is int and type(item[1]) is int:
						if item[0] < len(all_nodes) and item[1] < len(all_nodes):
							self.edges.append(item)	
						else:
							print("Edge index out of range! Skipping this node")
							continue
					else:
						print("Invalid tuple type! Skipping this edge")
						continue
				else:
					print("Invalid data type! Skipping this edge")
					continue
		else:
			print("Invalid data! No edges added")


	def list_edges(self):
		"""prints and returns list of all edges in network as tuples"""
		#print(self.edges)
		return self.edges


	def draw_network(self, title):
		"""Draws network using matplotlib"""
		"""Each edge is assigned a random position in unit circle"""
		"""Plot nodes"""
		xCoordinates = []
		yCoordinates = []
		all_nodes = self.list_nodes()
		for node in all_nodes:
			coordinates = node.get_coordinates()
			xCoordinates.append(coordinates[0])
			yCoordinates.append(coordinates[1])
		plt.plot(xCoordinates,yCoordinates,'bo')
		"""Plot edges"""
		all_edges = self.list_edges()
		for edge in all_edges:
			firstNode = all_nodes[edge[0]].get_coordinates()
			secondNode = all_nodes[edge[1]].get_coordinates()
			plt.plot([firstNode[0],secondNode[0]],[firstNode[1],secondNode[1]], 'b-')
		plt.axis([-1.5,1.5,-1.5,1.5])
		plt.title(title)
		#plt.show() 


	def get_adjacency_matrix(self):
		"""Returns adjacency matrix of network as numpy array"""
		all_nodes = self.list_nodes()
		numOfNodes = len(all_nodes)
		all_edges = self.list_edges()
		A = np.zeros((numOfNodes,numOfNodes))
		for edge in all_edges:
			A[edge[0],edge[1]] = 1
			A[edge[1],edge[0]] = 1
		return A

	def get_degree(self, node):
		"""Returns the degree of node in argument"""
		A = self.get_adjacency_matrix()
		degree = np.sum(A[node])
		return degree

	def get_mean_degree(self):
		"""Returns mean degree of network"""
		A = self.get_adjacency_matrix()
		all_nodes = self.list_nodes()
		numOfNodes = len(all_nodes)
		m = (np.sum(A)) / 2
		mean_degree = (2 * m) / numOfNodes
		return mean_degree 

	def get_clustering_coefficient(self):
		"""Returns culstering voefficient of network"""
		A = self.get_adjacency_matrix()
		all_nodes = self.list_nodes()
		numOfNodes = len(all_nodes)
		k = 0.0
		k_2 = 0.0
		for i in range(0, numOfNodes):
			a = np.sum(A[i])
			k += a
			k_2 += math.pow(a, 2)
		k = k / numOfNodes
		k_2 = k_2 / numOfNodes	
		#print("degree:" + str(k))
		#print("second moment: " + str(k_2))
		#print(k_2)	
		clustering_coefficient = math.pow((k_2-k), 2) / (numOfNodes * math.pow(k, 3))
		return clustering_coefficient


	def get_clustering_coefficient_2(self):
		"""Returns culstering coefficient of network"""
		A = self.get_adjacency_matrix()
		A2 = np.dot(A,A)
		A3 = np.dot(A2, A)
		numOfTriangles = np.trace(A3) / 6.0
		denominator = np.sum(A2) - np.trace(A2)
		clustering_coefficient =  6 * numOfTriangles / denominator
		return clustering_coefficient


	def exists_path(self, n, m):
		"""Returns length of shortest path between vertices n and m if exists, and zero otherwise"""
		edges = self.list_edges()
		paths = []
		travelled_length = 0
		paths_from_n = [x[1] if x[1] != n else x[0] for x in edges if n in x]
		for node in paths_from_n:
			paths.append([n,node])
		#print("Paths from n: " + str(paths))
		if len(paths_from_n) == 0:
			return 0
		elif [n,m] in paths:
			return 1
		else:
			while True:
				travelled_length += 1
				new_paths = []
				for path in paths:
					if path == None:
						continue
					start_point = path[0]
					#print("Starting point: " + str(start_point))
					end_point = path[-1]
					#print("End point: " + str(end_point))
					path_connects_to = [x[1] if x[1] != end_point else x[0] for x in edges if end_point in x]
					path_connects_to = [x for x in path_connects_to if x not in path]
					#print("Node " + str(end_point) + " connects to: " + str(path_connects_to))
					if m in path_connects_to:
						return travelled_length + 1
					for node in path_connects_to:
						#print("Node: " + str(node))
						#print("Current path: " + str(path))
						new_path = path
						new_paths.append(new_path + [node])
						#print("new path: " + str(new_path))
					#print("New paths: " + str(new_paths))
				if new_paths == []:
					return 0
				else:
					paths = [x for x in new_paths]
					#print("Paths from n: " + str(paths))







	def diameter(self):
		"""Returns the diameter of network"""
		"""Returns 0 if diameter is infinite"""
		edges = self.list_edges()
		nodes = self.list_nodes()
		pairs_of_nodes = []
		paths = []
		for i in range(len(nodes)):
			for j in range(len(nodes)):
				if i < j:
					pairs_of_nodes.append([i,j])
		#print("Pairs of nodes: " + str(pairs_of_nodes))
		longest_path = 0
		longest_start = 0
		longest_end = 0
		for pair in pairs_of_nodes:
			path = self.exists_path(pair[0],pair[1])
			if path == 0:
				#return 0
				pass
			paths.append(path)
			if path > longest_path:
				longest_start = pair[0]
				longest_end = pair[1]
				longest_path = path

		#print("Longest path starts at " + str(longest_start) + " and ends at " + str(longest_end) + " with length " + str(longest_path))
		plt.text(-1.45, 1.35, "Longest path starts at " + str(longest_start) + " and ends at " + str(longest_end) + " with length " + str(longest_path), fontsize=11)
		if 0 in paths:
			return [0, longest_path, sum(paths) / float(len(paths))]
		else:
			return [max(paths), longest_path, sum(paths) / float(len(paths))]



class node:
	"""Class for nodes in network"""
	"""Each node occupies a point inside unit circle with uniformly random position"""
	type = "Node"

	def __init__(self, name="Some node"):
		self.name = name
		self.coordinates = []
		t = 2*math.pi*random.random()
		u = random.random()
		u += random.random()
		r = u
		if u > 1:
			r = 2 - u
		self.coordinates.append(r*math.cos(t))
		self.coordinates.append(r*math.sin(t))

	def get_coordinates(self):
		return self.coordinates




def Gnm(nodes, edges):
	"""Generates a network of type G(n,m)"""
	n = network()
	nodes_list = [x for x in range(0, nodes)]
	n.add_node(nodes_list)
	allPossibleEdges = []
	for sourceNode in range(0,nodes):
		for targetNode in range(0, nodes):
			if targetNode > sourceNode:
				allPossibleEdges.append([sourceNode,targetNode])
	addedEdges = 0
	numberOfPossibleEdges = len(allPossibleEdges)
	while addedEdges < numberOfPossibleEdges and addedEdges < edges:
		randomChoice = random.randint(0,len(allPossibleEdges)-1)
		n.add_edge((allPossibleEdges[randomChoice][0], allPossibleEdges[randomChoice][1]))
		allPossibleEdges.pop(randomChoice)
		addedEdges += 1
	plt.figure(1)
	ax = plt.gca()
	circle1 = plt.Circle((0, 0), 1, color='b', fill=False)
	ax.add_artist(circle1)
	n.draw_network("G(n,m) model with " + str(nodes) + " nodes and " + str(edges) + " edges")
	#print(n.get_adjacency_matrix())
	mean_degree = n.get_mean_degree()
	plt.text(-1.45, 1.2, "Mean degree: " + str(mean_degree), fontsize=11)		
	#clustering_coefficient = n.get_clustering_coefficient()
	clustering_coefficient_2 = n.get_clustering_coefficient_2()
	plt.text(-1.45, 1.05, "Clustering coefficient: " + str(clustering_coefficient_2), fontsize=11)
	#print("Nodes: " + str(n.list_nodes()))
	#print("Edges: " + str(n.list_edges()))
	diameter = n.diameter()
	longest_path = diameter[1]
	average_path_length = diameter[2]
	diameter = diameter[0]
	#print("Diameter: " + str(diameter))
	plt.text(-1.45,0.9, "Diameter of network: " + str(diameter), fontsize=11)	
	plt.text(-1.45,0.75, "Average path length: " + str(average_path_length), fontsize=11)	
	nodes = n.list_nodes()
	for node in nodes:
		coords = node.get_coordinates()
		plt.text(coords[0]+0.05,coords[1]+0.05, str(node.name), fontsize=15, fontweight='bold')
	edges = n.list_edges()
	degrees = []
	for node in n.list_nodes():
		degrees.append(n.get_degree(node.name))
	return [mean_degree, clustering_coefficient_2, diameter, longest_path, average_path_length, len(edges), degrees]



def Gnp(nodes, p):
	"""Generates a network of type G(n,p)"""
	added = 0.0
	not_added = 0.0
	m = network()
	nodes_list = [x for x in range(0, nodes)]
	m.add_node(nodes_list)
	for sourceNode in range(0,nodes):
		for targetNode in range(0, nodes):
			if targetNode > sourceNode:
				dice = random.random()
				if dice <= p:
					m.add_edge((sourceNode, targetNode))
					added += 1
				else:
					not_added += 1
	#print(added)
	#print(not_added)
	total = added + not_added
	#print(added / total)
	plt.figure(2)
	ax = plt.gca()
	circle1 = plt.Circle((0, 0), 1, color='b', fill=False)
	ax.add_artist(circle1)
	m.draw_network("G(n,p) model with " + str(nodes) + " nodes and  p equal to " + str(p))
	#print(m.get_adjacency_matrix())
	mean_degree = m.get_mean_degree()
	plt.text(-1.45, 1.2, "Mean degree: " + str(mean_degree), fontsize=11)	
	#clustering_coefficient = m.get_clustering_coefficient()
	clustering_coefficient_2 = m.get_clustering_coefficient_2()
	plt.text(-1.45, 1.05, "Clustering coefficient: " + str(clustering_coefficient_2), fontsize=11)	
	#print("Nodes: " + str(m.list_nodes()))
	#print("Edges: " + str(m.list_edges()))
	diameter = m.diameter()
	longest_path = diameter[1]
	average_path_length = diameter[2]
	diameter = diameter[0]
	#print("Diameter: " + str(diameter))
	plt.text(-1.45,0.9, "Diameter of network: " + str(diameter), fontsize=11)	
	plt.text(-1.45,0.75, "Average path length: " + str(average_path_length), fontsize=11)	
	nodes = m.list_nodes()
	for node in nodes:
		coords = node.get_coordinates()
		plt.text(coords[0]+0.05,coords[1]+0.05, str(node.name), fontsize=15, fontweight='bold')
	edges = m.list_edges()
	degrees = []
	for node in m.list_nodes():
		degrees.append(m.get_degree(node.name))
	return [mean_degree, clustering_coefficient_2, diameter, longest_path, average_path_length, len(edges), degrees]




def main(number_of_vertices, number_of_edges, p, number_of_simulations):
	"""Generates many instances of random networks and collects statistics of properties"""
	gnp_average_degree = []
	gnp_average_clustering = []
	gnp_average_diameter = []
	gnp_average_longest_path = []
	gnp_average_path_length = []
	gnp_number_of_edges = []
	gnp_degree_list = []

	for i in range(0,number_of_simulations):
		gn = Gnp(number_of_vertices, p)
		gnp_average_degree.append(gn[0])
		gnp_average_clustering.append(gn[1])
		gnp_average_diameter.append(gn[2])
		gnp_average_longest_path.append(gn[3])
		gnp_average_path_length.append(gn[4])
		gnp_number_of_edges.append(gn[5])
		gnp_degree_list += gn[6]

	gnm_average_degree = []
	gnm_average_clustering = []
	gnm_average_diameter = []
	gnm_average_longest_path = []
	gnm_average_path_length = []
	gnm_number_of_edges = []
	gnm_degree_list = []

	for i in range(0,number_of_simulations):
		gm = Gnm(number_of_vertices, number_of_edges)
		gnm_average_degree.append(gm[0])
		gnm_average_clustering.append(gm[1])
		gnm_average_diameter.append(gm[2])
		gnm_average_longest_path.append(gm[3])
		gnm_average_path_length.append(gm[4])
		gnm_number_of_edges.append(gm[5])
		gnm_degree_list += gm[6]

	print("Gnp average mean degree: " + str(sum(gnp_average_degree) / float(len(gnp_average_degree))))
	print("Standard deviation: " + str(np.std(gnp_average_degree)))
	print("Gnp average clustering coefficient: " + str(sum(gnp_average_clustering) / float(len(gnp_average_clustering))))
	print("Standard deviation: " + str(np.std(gnp_average_clustering)))
	print("Gnp average diameter: " + str(sum(gnp_average_diameter) / float(len(gnp_average_diameter))))
	print("Standard deviation: " + str(np.std(gnp_average_diameter)))
	print("Gnp average longest path: " + str(sum(gnp_average_longest_path) / float(len(gnp_average_longest_path))))
	print("Standard deviation: " + str(np.std(gnp_average_longest_path)))
	print("Gnp average path length: " + str(sum(gnp_average_path_length) / float(len(gnp_average_path_length))))
	print("Standard deviation: " + str(np.std(gnp_average_path_length)))
	print("Gnp average number of edges: " + str(sum(gnp_number_of_edges) / float(len(gnp_number_of_edges))))
	print("Standard deviation: " + str(np.std(gnp_number_of_edges)))

	print("\n-----------------------------------------------------------\n")

	print("Gnm average mean degree: " + str(sum(gnm_average_degree) / float(len(gnm_average_degree))))
	print("Standard deviation: " + str(np.std(gnm_average_degree)))
	print("Gnm average clustering coefficient: " + str(sum(gnm_average_clustering) / float(len(gnm_average_clustering))))
	print("Standard deviation: " + str(np.std(gnm_average_clustering)))
	print("Gnm average diameter: " + str(sum(gnm_average_diameter) / float(len(gnm_average_diameter))))
	print("Standard deviation: " + str(np.std(gnm_average_diameter)))
	print("Gnm average longest path: " + str(sum(gnp_average_longest_path) / float(len(gnp_average_longest_path))))
	print("Standard deviation: " + str(np.std(gnp_average_longest_path)))
	print("Gnm average path length: " + str(sum(gnm_average_path_length) / float(len(gnm_average_path_length))))
	print("Standard deviation: " + str(np.std(gnm_average_path_length)))
	print("Gnm average number of edges: " + str(sum(gnm_number_of_edges) / float(len(gnm_number_of_edges))))
	print("Standard deviation: " + str(np.std(gnm_number_of_edges)))

	print("\n---------------------------MATHS FOR G(n,p)--------------------------------\n")


	c_mean_number_of_edges = binomial(number_of_vertices, 2) * p
	print("Calculated mean number of edges: " + str(c_mean_number_of_edges))
	c_mean_degree = (number_of_vertices - 1) * p
	print("Calculated mean degree: " + str(c_mean_degree))
	c_clustering = c_mean_degree / (number_of_vertices - 1)
	print("Calculated clustering coefficient: " + str(c_clustering))
	c_diameter = math.log(number_of_edges)/math.log(c_mean_degree)
	print("Calculated diameter (minus constant value): " + str(c_diameter))


	#Sample figures and degree distributions as a histogram
	plt.close()
	plt.cla()
	plt.clf()
	gnp = Gnp(number_of_vertices, p)
	gnm = Gnm(number_of_vertices, number_of_edges)
	plt.figure(3)
	hist_1 = plt.hist(gnp[6])
	plt.xlim(xmin=0, xmax=max(gnp[6])+1)
	plt.xticks(range(int(max(gnp[6]))+1))
	plt.title("Degree distribution of sample G(n,p) model with " + str(number_of_vertices) + " and p equal to " + str(p))
	plt.xlabel("Degree")
	plt.figure(4)
	hist_2 = plt.hist(gnm[6])
	plt.xlim(xmin=0, xmax=max(gnm[6])+1)
	plt.xticks(range(int(max(gnm[6]))+1))
	plt.title("Degree distribution of sample G(n,m) model with " + str(number_of_vertices) + " vertices and " + str(number_of_edges) + " number of edges")
	plt.xlabel("Degree")
	plt.show()







if __name__ == "__main__":
	"""Arguments: number of vertices, number of edges in G(n,m), p in G(n,p), nmber of runs of simulation"""
	main(20,57,0.3,10)

