import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import networkx as nx


class Graph:
	def __init__(self, adjacency_matrix=[]):
		self.adjacency_matrix = adjacency_matrix
		self.edges = self._get_edges()
		
	def add_vertex(self):
		temp_matrix = np.zeros((len(self.adjacency_matrix)+1, len(self.adjacency_matrix)+1))
		temp_matrix[:len(self.adjacency_matrix), :len(self.adjacency_matrix)] = self.adjacency_matrix
		self.adjacency_matrix = temp_matrix

	def add_edge(self, v1, v2):
		if v1 < len(self.adjacency_matrix) and v2 < len(self.adjacency_matrix):
			self.adjacency_matrix[v1, v2] = 1
			self.adjacency_matrix[v2, v1] = 1
		else: 
			print("You must choose a vertex from the existing ones. From 0 to", len(self.adjacency_matrix)-1)

	def get_adjacency_matrix(self):
		return self.adjacency_matrix

	def _check_symmetric(self, a, tol=1e-8):
	    return np.allclose(a, a.T, atol=tol)

	def _is_regular_graph(self):
		#check that the graph is regular
		sum_in_rows = np.sum(self.adjacency_matrix, axis=1)
		return all(sum_in_rows == sum_in_rows[0])

	def _get_edges(self):
		edges = []
		for i in range(len(self.adjacency_matrix)):
			for j in range(len(self.adjacency_matrix)):
				if self.adjacency_matrix[i][j] != 0 and i > j:
					edges.append([i, j])
		return edges

	def _neighbours_of_vertex(self, vertex):
		neighbours = []
		i=0
		for n in self.adjacency_matrix[vertex]:
			if n==1:
				neighbours.append(i)
			i+=1
		return neighbours

	def is_distance_regular_graph(self):
		if not self._is_regular_graph():
			return False
		#check that the graph is distance-regular
		bs, cs = [], []
		for v in range(len(self.adjacency_matrix)):
			distances = self.Dijkstra(v)
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
				for w in self._neighbours_of_vertex(y):
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
					for w in self._neighbours_of_vertex(y):
						if distances[w] == i+1:
							bi += 1
						elif distances[w] == i-1:
							ci += 1
					if bij != bi or cij != ci:
						return False
		return True, (bs[:diam] + cs[1:])

	def get_spectrum(self):
		from math import ceil
		w,v = np.linalg.eig(self.adjacency_matrix)
		# count of various lambda's must be equal diam + 1
		# getting rid of errors in calculations
		if len(np.unique(w)) != max(self.Dijkstra(0))+1:
			w = np.round(w)
		return np.unique(w, return_counts = True)

	def Dijkstra(self, S):
		# S - vertex from which we count the distance to others
		matrix = self.adjacency_matrix.copy()
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

	def intersect_array(self):
		is_dist_regular = self.is_distance_regular_graph()
		if type(is_dist_regular) is not bool and is_dist_regular[0]:
			return is_dist_regular[1]
		else:
			raise ValueError("Graph isn't distance regular, we can't find intersect array")
		return

	def intersect_matrix(self):
		distances = self.Dijkstra(0)
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
			for w in self._neighbours_of_vertex(y):
				if distances[w] == i+1:
					bi += 1
				elif distances[w] == i-1:
					ci += 1
			bs.append(bi)
			cs.append(ci)
		a = []
		for i in range(diam+1):
			a.append(bs[0]-bs[i]-cs[i])
		mtrx = np.zeros((diam+1, diam+1))
		for i in range(diam+1):
			mtrx[i,i] = a[i]
		for i in range(diam):
			mtrx[i, i+1] = bs[i]
			mtrx[i+1, i] = cs[i+1]
		return mtrx

	def get_spectrum_dist_regular(self):
		distances = self.Dijkstra(0)
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
			for w in self._neighbours_of_vertex(y):
				if distances[w] == i+1:
					bi += 1
				elif distances[w] == i-1:
					ci += 1
			bs.append(bi)
			cs.append(ci)
		a = []
		for i in range(diam+1):
			a.append(bs[0]-bs[i]-cs[i])
		intersect_matrix = np.zeros((diam+1, diam+1))
		for i in range(diam+1):
			intersect_matrix[i,i] = a[i]
		for i in range(diam):
			intersect_matrix[i, i+1] = bs[i]
			intersect_matrix[i+1, i] = cs[i+1]
		w, v = np.linalg.eig(intersect_matrix)
		us = []
		for theta in w:
			u = [1, theta/bs[0]]
			for i in range(1, diam):
				u.append((theta*u[i]-cs[i]*u[i-1]-a[i]*u[i])/bs[i])
			us.append(u)
		ks = [1]
		ms = []
		for i in range(diam):
			ks.append(bs[i]*ks[i]/cs[i+1])
		for vec in us:
			ms.append(len(distances)/sum([ks[i]*vec[i]*vec[i] for i in range(diam+1)]))
		return np.round(w), np.round(ms)

	def layed_represent(self, start_point):
		distances = self.Dijkstra(start_point)
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

	def draw_layer_represent(self, start_point):
		nodelist = self.layed_represent(start_point)
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

		
		edge_pos = np.asarray([([[p_widht[e[0]], p_height[e[0]]],[p_widht[e[1]], p_height[e[1]]]] ) for e in self.edges])
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
		
	def draw(self):
	    rows, cols = np.where(self.adjacency_matrix == 1)
	    edges = zip(rows.tolist(), cols.tolist())
	    gr = nx.Graph()
	    gr.add_edges_from(edges)
	    nx.draw(gr, node_size=500, with_labels=True)
	    plt.show()


c4 = np.array([[0,1,0,1], [1,0,1,0], [0,1,0,1], [1,0,1,0]])
a1 = np.array([[0,1,1], [1,0,0], [1,0,0]])
a3 = np.array([[0,1,1,1,0,0], [1,0,0,0,1,1], [1,0,0,0,1,1], [1,0,0,0,1,1], [0,1,1,1,0,0], [0,1,1,1,0,0]])
c3 = np.array([[0,1,1], [1,0,1], [1,1,0]])
c5 = np.array([[0,1,0,0,1], [1,0,1,0,0], [0,1,0,1,0], [0,0,1,0,1], [1,0,0,1,0]])
nd = np.array([[0,1,1,0,1], 
			   [1,0,1,0,0],
			   [1,1,0,1,1],
			   [0,0,1,0,1],
			   [1,0,1,1,0]])


def main():
	gr = Graph(c5)
	#print(gr.get_adjacency_matrix())
	#print(gr.is_distance_regular_graph())
	#print(gr.intersect_array())
	print(gr.get_spectrum())
	#gr.add_vertex()
	#gr.add_edge(5, 0)
	#gr.add_edge(4, 5)
	#print(gr.intersect_matrix())
	print(gr.get_spectrum_dist_regular())
	print(gr.is_distance_regular_graph())
	#gr.draw()
	#gr.draw_layer_represent(0)


if __name__ == '__main__':
	main()