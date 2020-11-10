import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

c4 = np.array([[0,1,0,1], [1,0,1,0], [0,1,0,1], [1,0,1,0]])
a1 = np.array([[0,1,1], [1,0,0], [1,0,0]])
a3 = np.array([[0,1,1,1,0,0], [1,0,0,0,1,1], [1,0,0,0,1,1], [1,0,0,0,1,1], [0,1,1,1,0,0], [0,1,1,1,0,0]])
c3 = np.array([[0,1,1], [1,1,0], [1,0,1]])
c5 = np.array([[0,1,0,0,1], [1,0,1,0,0], [0,1,0,1,0], [0,0,1,0,1], [1,0,0,1,0]])
nd = np.array([[0,1,0,0,1], [1,0,1,0,0], [1,1,0,1,0], [0,0,1,0,1], [1,0,1,1,0]])
tt = np.array([[1, 0, 0, 0, 0, 0],  [0, 1, 1, 1, 1, 0],  [0, 0, 0, 0, 0, 1],
        [0, 1, 0, 0, 0, 0],  [1, 0, 1, 1, 0, 1],  [0, 0, 0, 0, 1, 0],
        [0, 0, 1, 0, 0, 0],  [1, 1, 0, 0, 1, 1],  [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0, 0],  [1, 1, 0, 0, 1, 1],  [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, 0],  [1, 0, 1, 1, 0, 1],  [0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1], [0, 1, 1, 1, 1, 0], [1, 0, 0, 0, 0, 0]])

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
    nx.draw(gr, node_size=500)
    plt.show()

#print(c5)
print(is_distance_regular_graph(c5))
print(layed_represent(0, c5))
#print(intersect_array(c5))
print(get_spectrum(c5))

#draw_graph(c5)
