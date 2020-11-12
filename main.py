import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import networkx as nx

c4 = np.array([[0,1,0,1], [1,0,1,0], [0,1,0,1], [1,0,1,0]])
a1 = np.array([[0,1,1], [1,0,0], [1,0,0]])
a3 = np.array([[0,1,1,1,0,0], [1,0,0,0,1,1], [1,0,0,0,1,1], [1,0,0,0,1,1], [0,1,1,1,0,0], [0,1,1,1,0,0]])
c3 = np.array([[0,1,1], [1,1,1], [1,1,1]])
c5 = np.array([[0,1,0,0,1], [1,0,1,0,0], [0,1,0,1,0], [0,0,1,0,1], [1,0,0,1,0]])
nd = np.array([[0,1,1,0,1], 
			   [1,0,1,0,0],
			   [1,1,0,1,1],
			   [0,0,1,0,1],
			   [1,0,1,1,0]])


def check_symmetric(a, tol=1e-8):
    return np.allclose(a, a.T, atol=tol)

def is_regular_graph(matrix):
	#check that the graph is regular
	sum_in_rows = np.sum(matrix, axis=1)
	return all(sum_in_rows == sum_in_rows[0])

def is_distance_regular_graph(matrix):
	if not is_regular_graph(matrix):
		return False
	#check that the graph is distance-regular
	bs, cs = [], []
	for v in range(len(matrix)):
		distances = Dijkstra(v, matrix)
		diam = max(distances)
		for i in range(diam+1):
			ys = []
			for w in range(len(distances)):
				if distances[w] == i:
					ys.append(w)
			#now d(v,y) = i for y in ys

			#bi = number neighbours of y at distance i+1 from v
			#ci = number neighbours of y at distance i-1 from v
			bi, ci, bij, cij = 0, 0, 0, 0
			#find bi ci for first neighbour
			y = ys[0]
			for w in neighbours_of_vertex(matrix, y):
				if distances[w] == i+1:
					bi += 1
				elif distances[w] == i-1:
					ci += 1
			bij, cij = bi, ci
			if v == 0:
				bs.append(bi)
				cs.append(ci)
			#find bi ci for other neighbours and check if they equal
			for y in ys:
				bi, ci = 0, 0
				for w in neighbours_of_vertex(matrix, y):
					if distances[w] == i+1:
						bi += 1
					elif distances[w] == i-1:
						ci += 1
				if bij != bi or cij != ci:
					return False
	return True, (bs[:diam] + cs[1:])


def is_distance_regular_graph1(matrix):
	if not is_regular_graph(matrix):
		return False
	#check that the graph is distance-regular
	bs, cs = [], []
	distances = Dijkstra(0, matrix)
	diam = max(distances)
	for i in range(diam+1):
		ys = []
		for w in range(len(distances)):
			if distances[w] == i:
				ys.append(w)
		#now d(v,y) = i for y in ys

		#bi = number neighbours of y at distance i+1 from v
		#ci = number neighbours of y at distance i-1 from v
		bi, ci, bij, cij = 0, 0, 0, 0
		#find bi ci for first neighbour
		y = ys[0]
		for w in neighbours_of_vertex(matrix, y):
			if distances[w] == i+1:
				bi += 1
			elif distances[w] == i-1:
				ci += 1
		bij, cij = bi, ci
		bs.append(bi)
		cs.append(ci)
		#find bi ci for other neighbours and check if they equal
		for y in ys:
			bi, ci = 0, 0
			for w in neighbours_of_vertex(matrix, y):
				if distances[w] == i+1:
					bi += 1
				elif distances[w] == i-1:
					ci += 1
			if bij != bi or cij != ci:
				return False
	return True, (bs[:diam] + cs[1:])

def get_spectrum(matrix):
	w,v = np.linalg.eig(matrix)
	# count of various lambda's must be equal diam + 1
	# getting rid of errors in calculations
	if len(np.unique(w)) != max(Dijkstra(0,matrix))+1:
		w = w.astype(int)
	return np.unique(w, return_counts = True)

def Dijkstra(S, matrix):
	# S - vertex from which we count the distance to others
	matrix = matrix.copy()
	matrix = np.where(matrix == 0, float('inf'), matrix)
	valid = [True]*len(matrix)        
	weight = [float('inf')]*len(matrix)
	weight[S] = 0
	for i in range(len(matrix)):
		min_weight = float('inf')
		ID_min_weight = -1
		for j in range(len(matrix)):
			if valid[j] and weight[j] < min_weight:
				min_weight = weight[j]
				ID_min_weight = j
		for k in range(len(matrix)):
			if weight[ID_min_weight] + matrix[ID_min_weight][k] < weight[k]:
				weight[k] = weight[ID_min_weight] + matrix[ID_min_weight][k]
		valid[ID_min_weight] = False
	return np.array(weight).astype(int)

def neighbours_of_vertex(matrix, vertex):
	neighbours = []
	i=0
	for n in matrix[vertex]:
		if n==1:
			neighbours.append(i)
		i+=1
	return neighbours

def intersect_array(matrix):
	distances = Dijkstra(0, matrix)
	diam = max(distances)
	bs, cs = [], []
	for i in range(diam+1):
		for w in range(len(distances)):
			if distances[w] == i:
				y = w
				break
		#now d(v,y) = i

		#bi = number neighbours of y at distance i+1 from v
		#ci = number neighbours of y at distance i-1 from v
		bi = 0
		ci = 0
		for w in neighbours_of_vertex(matrix, y):
			if distances[w] == i+1:
				bi += 1
			elif distances[w] == i-1:
				ci += 1
		bs.append(bi)
		cs.append(ci)
	return (bs[:diam] + cs[1:])

def layed_represent(start_point, matrix):
	distances = Dijkstra(start_point, matrix)
	diam = max(distances)
	result_array = []
	for i in range(diam+1):
		temp_array = []
		num = 0
		for dist in distances:
			if dist==i:
				temp_array.append(num)
			num += 1
		result_array.append(temp_array)
	return result_array


def draw_graph(adjacency_matrix):
    rows, cols = np.where(adjacency_matrix == 1)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    gr.add_edges_from(edges)
    nx.draw(gr, node_size=500, with_labels=True)
    plt.show()


#print(c5)
#print(is_distance_regular_graph(c5))
#print(layed_represent(0, c5))
#print(intersect_array(c5))
#print(get_spectrum(c5))

#draw_graph(c4)

def get_edges(adjacency_matrix):
	edges = []
	for i in range(len(adjacency_matrix)):
		for j in range(len(adjacency_matrix)):
			if adjacency_matrix[i][j] != 0 and i > j:
				edges.append([i, j])
	return edges

def draw_layer_represent(start_point, adjacency_matrix):
	nodelist = layed_represent(start_point, adjacency_matrix)
	print(nodelist)

	height, widht = 600, 600
	f_widht = widht/(len(nodelist)+1)
	f_heights = [height/(len(i)+1) for i in nodelist]
	
	#find positions for vertexes
	p_widht, p_height = [], []
	lower_bound_widht = 0
	for i, node in enumerate(nodelist):
		lower_bound_widht = lower_bound_widht + f_widht
		for m in range(len(node)):
			p_widht.append(lower_bound_widht)
		lower_bound_height, upper_bound_height = 0, 600
		for j, point in enumerate(node):
			if j%2==0:
				lower_bound_height = lower_bound_height + f_heights[i]
				p_height.append(lower_bound_height)
			else:
				upper_bound_height = upper_bound_height - f_heights[i]
				p_height.append(upper_bound_height)

	list_of_points = []
	for column in nodelist:
		for vertex in column:
			list_of_points.append(vertex)

	#sort points with positions by names
	z = zip(p_widht, p_height, list_of_points)
	zs = sorted(z, key=lambda tup: tup[2])
	p_widht = [z[0] for z in zs]
	p_height = [z[1] for z in zs]
	list_of_points = [z[2] for z in zs]

	cf = plt.gcf()
	ax = cf.gca()
	ax.set_axis_off()

	edges = get_edges(adjacency_matrix)
	edge_pos = np.asarray([([[p_widht[e[0]], p_height[e[0]]],[p_widht[e[1]], p_height[e[1]]]] ) for e in edges])
	edge_collection = LineCollection(edge_pos)
	ax.add_collection(edge_collection)

	for iteration, label in enumerate(list_of_points):
		plt.text(p_widht[iteration], p_height[iteration], label, size=12, color='w', family="sans-serif",
	  		weight="normal",
	    	alpha=None,
	    	bbox=None,
	    	horizontalalignment="center",
	    	verticalalignment="center",
	    	clip_on=True)
	plt.scatter(p_widht, p_height, s=300, c="#1f78b4")
	plt.show()
	
draw_layer_represent(0, nd)

#draw_graph(nd)
#print(Dijkstra(0, nd))
'''
rows, cols = np.where(tt == 1)
edges = zip(rows.tolist(), cols.tolist())
gr = nx.Graph()
gr.add_edges_from(edges)
print(nx.dijkstra_path(gr, 0, 2))
'''