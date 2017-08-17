#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <functional>
#include "heap.cpp"

#define MAXWEIGHT 10000


class Graph
{
public:
	struct Node
	{
		int vertex;
		double weight;
		
		bool operator<(const Node x) const
		{
			return weight < x.weight;
		}
		
		bool operator>(const Node x) const
		{
			return weight > x.weight;
		}

		bool operator==(Node x)
		{
			if (vertex == x.vertex && weight == x.weight)
				return true;
			else
				return false;
		}
	};

	int vertices, edges;
	std::vector<std::list<Node>> adj;
	bool digraph;
	
	//Constructor, reads edges with weights
	Graph(int v, int e, bool orgraph)
	{
		vertices = v;
		edges = e;
		digraph = orgraph;
		for (int i = 0; i < v; i++)
		{
			std::list<Node> a;
			adj.push_back(a);
		}
		for (int i = 0; i < e; i++)
		{
			int v1, v2;
			double w;
			std::cin >> v1 >> v2 >> w;
			adj[v1].push_back(Node{ v2, w });
			if (!digraph)
				adj[v2].push_back(Node{ v1, w });
		}
	}
};


void dijkstra(double * weights, Graph& graph, int start)
{
	int *indices = new int[graph.vertices]; //indices[v] is the index of v in the heap
	for (int i = 0; i < graph.vertices; i++)
		weights[i] = MAXWEIGHT;
	weights[start] = 0;
	Heap < Graph::Node, std::greater<Graph::Node> > heap(graph.vertices, indices);

	for (int i = 0; i < graph.vertices; i++)
		heap.push(Graph::Node{ i, weights[i] });
	
	while (heap.n > 1)
	{
		Graph::Node u = heap.pop();
		for (std::list<Graph::Node>::const_iterator i = graph.adj[u.vertex].begin();
			i != graph.adj[u.vertex].end(); i++)
			if (weights[i->vertex] > weights[u.vertex] + i->weight)
			{
				//relaxation
				weights[i->vertex] = weights[u.vertex] + i->weight;
				heap.arr[indices[i->vertex]].weight = weights[u.vertex] + i->weight;
				heap.refresh(indices[i->vertex]);
			}
	}
	delete[] indices;
}


//Finds distances from start vertex to all the others in time O(V*E).
//Returns true if successful, false if there are cycles with negative weight.
bool ford_bellman(double * weights, Graph& graph, int start)
{
	//init:
	for (int i = 0; i < graph.vertices; i++)
		weights[i] = MAXWEIGHT;
	weights[start] = 0;
	//main procedure:
	for (int i = 0; i < graph.vertices - 1; i++)
		for (int j = 0; j < graph.vertices; j++)
			for (std::list<Graph::Node>::iterator iter = graph.adj[j].begin();
				iter != graph.adj[j].end(); iter++)
				if (weights[iter->vertex] > weights[j] + iter->weight)
					weights[iter->vertex] = weights[j] + iter->weight; //relaxation
	
	//checking for cycles with negative weight:
	for (int j = 0; j < graph.vertices; j++)
		for (std::list<Graph::Node>::iterator iter = graph.adj[j].begin();
			iter != graph.adj[j].end(); iter++)
			if (weights[iter->vertex] > weights[j] + iter->weight)
				return false;
	return true;
}


//Finds shortest distances between all the vertices.
//Needs adjacency-matrix. Time O(V^3).
//Returns false if there are cycles with negative weight.
bool floyd_warshall(std::vector<std::vector<double> >& weights)
{
	int v = weights.size();
	for (int k = 0; k < v; k++)
		for (int i = 0; i < v; i++)
			for (int j = 0; j < v; j++)
				if (weights[i][j] > weights[i][k] + weights[k][j])
					weights[i][j] = weights[i][k] + weights[k][j];
	//checking for cycles with negative weight (my idea)
	for (int i = 0; i < v; i++)
		if (weights[i][i] != 0)
			return false;
	return true;
}


/*
init graph as adjacency-matrix
std::vector<std::vector<double>> graph(v);
for (int i = 0; i < v; i++)
for (int j = 0; j < v; j++)
if (i == j)
graph[i].push_back(0);
else
graph[i].push_back(MAXWEIGHT);
for (int i = 0; i < e; i++)
{
int u, v;
double w;
std::cin >> u >> v >> w;
graph[u][v] = w;
graph[v][u] = w;//if not directed
}
*/


int main()
{
	int v, e;
	std::cin >> v >> e;
	Graph graph(v, e, false);
		
	double * weights = new double[v];
	dijkstra(weights, graph, 0);
	for (int i = 0; i < v; i++)
		std::cout << weights[i] << " ";

	delete[] weights;

	std::cout << "\n";
	system("pause");
	return 0;
}