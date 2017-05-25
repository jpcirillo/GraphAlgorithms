import java.util.ArrayList;

//Just a general note, whenever an emoji is output, it's because I have to return some int to be printed by the main method, so I return a 0 and make it into a face so it's less awkward for all of us
public class GraphAlgorithms {
	private final int Dijkstra = 1;
	private final int aStar = 2;

	public GraphAlgorithms() { // default constructor
	}

	public int processGraph(int[][] graph, int algorithm, int s, int d) {
		// Precondition: graph is an adjacency array representing a graph,
		// algorithm is an algorithm code (1 for Dijkstra, 2 for A*), s is a
		// valid node, and d is also a valid node
		// Postcondition: the length of the shortest path from s to d in the
		// graph, as found by the specified algorithm, is returned.
		switch (algorithm) {
		case Dijkstra:
			return Dijkstra(graph, s, d);
		case aStar:
			return aStar(graph, s, d);
		default: // Check for invalid algorithm code
			System.out.println("Invalid algorithm, please enter 1 for Dijkstra or 2 for A* ;");
			return 0;
		}
	}

	private int Dijkstra(int[][] graph, int s, int d) {
		// Precondition: graph is an adjacency array representing a graph, s is
		// a valid node, and d is also a valid node
		// Postcondition: the shortest path from s to d, as found by the
		// Dijkstra algorithm, is printed, and the weighted length of that path
		// is returned
		int dim = graph[0].length; // dim is the number of nodes

		if (s >= dim || s < 0) { // checks for invalid source vertex
			System.out.print("The source vertex " + s + " does not exist. o.");
			return 0;
		}

		if (d >= dim || d < 0) { // checks for invalid destination vertex
			System.out.print("The destination vertex " + d + " does not exist. o.");
			return 0;
		}

		// ArrayLists to store visited and unvisited nodes, as is required by
		// the Dijkstra algorithm
		ArrayList<Integer> visited = new ArrayList<>();
		ArrayList<Integer> unvisited = new ArrayList<>();

		// ArrayList which will later be used to store the nodes that make up
		// the shortest path, in order
		ArrayList<Integer> Path = new ArrayList<>();

		int distance[] = new int[dim]; // Integer array which stores, at each
										// index, the shortest distance found so
										// far to get to the index's
										// corresponding node

		int cameFrom[] = new int[dim]; // Integer array which stores, at each
										// index, the node that the index's
										// corresponding node was come to from
		for (int i = 0; i < dim; i++) { // Initialize each cameFrom to -1, as
										// this is an impossible node
			cameFrom[i] = -1;
		}

		for (int i = 0; i < dim; i++) { // Iterates through each node
			unvisited.add(i); // All nodes are originally unvisited
			if (i != s) { // Unless the node is the start node, we initialize
							// its distance to infinity (or 2.147 billion, in
							// this case)
				distance[i] = Integer.MAX_VALUE;
			} else { // The start node's distance is initialized to 0
				distance[i] = 0;
			}
		}

		do {
			int current = minDistNode(unvisited, distance); // Each time that
															// the algorithm
															// runs through, it
															// operates on the
															// unvisited node
															// which has the
															// smallest
															// distance, which
															// minDistNode(unvisited,
															// distance) returns
			visited.add((Integer) current); // the current node is now visited,
											// and no longer unvisited
			unvisited.remove((Integer) current);

			if (unvisited.contains(d)) { // this will be reached when current is
											// the destination node, and in that
											// case, we have already finished
											// constructing our path, and don't
											// need to do any of the stuff in
											// this loop
				for (int i = 0; i < dim; i++) { // iterates through each node
					if (graph[current][i] != 0 && (distance[i] > distance[current] + graph[current][i])) { // if
																											// graph[current][i]
																											// !=
																											// 0,
																											// current
																											// is
																											// adjacent
																											// to
																											// i
						// we check whether the nodes adjacent to current
						// currently have a distance that is worse than the
						// distance to reach it by way of the current node, and
						// if it is we update it with this better distance
						distance[i] = distance[current] + graph[current][i];
						cameFrom[i] = current; // whenever we improve the
												// distance to i, it means we
												// have found the new best way
												// to reach i, so we set the
												// node that we came to it from
												// to the current node
					}
				}
			}
		} while (unvisited.contains(d)); // iterates until the destination has
											// been visited, which would mean we
											// have found the path

		if (distance[d] < 0) { // When there is a disconnect in the graph
								// stopping a path from existing, the final
								// distance to reach the destination ends up
								// being a negative integer, so this checks for
								// that error.
			System.out.print("No path is present from " + s + " to " + d + ", disconnect exists in graph. o.");
			return 0;
		}

		Path = constructPath(cameFrom, s, d, Path); // calls constructPath,
													// which adds the nodes of
													// the path to Path in order
		System.out.print(Path + " ");
		return distance[d]; // returns final distance to reach the destination
	}

	int minDistNode(ArrayList<Integer> unvis, int[] dists) {
		// Precondition: unvis is an ArrayList of integers which represent the
		// nodes that have not yet been visited, and dists is an array filled
		// with the distances to reach each node
		// Postcondition: The number of the node which is unvisited and has the
		// lowest distance is returned.

		int min = dists[unvis.get(0)]; // sets initial values to min and
										// minNode, which are the distance for
										// the first node in unvis and the first
										// node in unvis
		int minNode = unvis.get(0);

		for (int element : unvis) {
			if (dists[element] < min) { // checks if each unvisited element's
										// distance is lower than the current
										// min distance
				min = dists[element]; // if it is, it updates the current min
										// distance and the min node
				minNode = element;
			}
		}
		return minNode; // returns the min node found
	}

	public ArrayList<Integer> constructPath(int[] from, int s, int current, ArrayList<Integer> Path) {
		// Precondition: from is an integer array storing the integer
		// representations of the nodes from which each node was visited, s is
		// the start node, the user should initially pass the destination node
		// as current, and the user should pass an empty ArrayList<Integer> as
		// Path
		// Postcondition: Path is filled with the nodes that form the shortest
		// path from the start node to the destination node in order and then
		// Path is returned.

		if (current == s) { // base case: we have reached the start node
			Path.add(s); // we add the start node to Path, and then return Path
			return Path;
		} else { // recursive case: call constructPath, with the new destination
					// being the node that the current destination was reached
					// through
			Path = constructPath(from, s, from[current], Path);
			Path.add(current); // after recursive call, the current destination
								// is added. Because this is after the recursive
								// call, Path will end up being in the right
								// order. We then return Path
			return Path;
		}
	}

	private int aStar(int[][] graph, int s, int d) {
		// Precondition: graph is an adjacency array representing a graph, s is
		// a valid node, and d is also a valid node
		// Postcondition: the shortest path from s to d, as found by the A*
		// algorithm, is printed, and the weighted length of that path is
		// returned

		// ArrayLists openList and closedList, as are needed for the A*
		// algorithm
		ArrayList<Integer> openList = new ArrayList<>();
		ArrayList<Integer> closedList = new ArrayList<>();

		ArrayList<Integer> Path = new ArrayList<>(); // ArrayList which will
														// store the shortest
														// path to reach the
														// destination

		// dim is the number of nodes
		int dim = graph[0].length;

		if (s >= dim || s < 0) { // checks that source vertex is a valid vertex
			System.out.print("The source vertex " + s + " does not exist. :^");
			return 0;
		}
		if (d >= dim || d < 0) { // checks that the destination vertex is a
									// valid vertex
			System.out.print("The destination vertex " + d + " does not exist. :^");
			return 0;
		}

		// int arrays to store each node's g cost and h cost
		int[] gCosts = new int[dim];
		int[] hCosts = heuristic(graph, d);

		// just like in Dijkstra's method, this array stores what node was used
		// to reach each node, and will be used to reconstruct the path.
		int[] cameFrom = new int[dim];

		int current = s; // we first consider the start node, and also add it to
							// openList
		openList.add(current);
		gCosts[current] = 0; // because current is the start node, its g cost is
								// 0

		while (!openList.isEmpty()) { // if open list is empty, it means that
										// there is no path to reach the
										// destination
			int minF = gCosts[openList.get(0)] + hCosts[openList.get(0)]; // initial
																			// value
																			// for
																			// the
																			// minimum
																			// f
																			// cost,
																			// which
																			// is
																			// the
																			// sum
																			// of
																			// the
																			// g
																			// and
																			// h
																			// cost
			int minNode = openList.get(0); // initial value for the min node,
											// which is the first node in open
											// list
			for (int element : openList) { // iterates through the elements in
											// openList and if any of them has a
											// lower f cost than the current
											// minF, it is assigned as the new
											// minNode
				if (minF > gCosts[element] + hCosts[element]) {
					minF = gCosts[element] + hCosts[element];
					minNode = element;
				}
			}

			current = minNode; // This is the node in openList with the lowest f
								// cost, which is what we always consider in the
								// A* algorithm

			if (current == d) { // if current is the destination, then we have
								// found our path

				Path = constructPath(cameFrom, s, d, Path); // calls
															// constructPath,
															// which stores the
															// shortest path in
															// an ArrayList
															// (Path)
				System.out.print(Path + " ");

				int thisDest = d; // initial values for the current destination
									// and the distance to reach the destination
				int theDist = 0;
				while (thisDest != s) { // iterates backwards through cameFrom
										// and adds the weights between each 2
										// nodes of the path to find the total
										// weight to reach d
					theDist += graph[cameFrom[thisDest]][thisDest];
					thisDest = cameFrom[thisDest];
				}
				return theDist; // this is the final sum of all the weights of
								// the edges we traversed
			}

			openList.remove((Integer) current); // Once we have looked at
												// current, we remove it from
												// openList and add it to
												// closedList, because we are
												// done with it.
			closedList.add(current);

			for (int i = 0; i < dim; i++) { // iterates through all nodes
				if (graph[current][i] != 0 && !closedList.contains(i)) { // checks
																			// if
																			// each
																			// node
																			// is
																			// adjacent
																			// to
																			// the
																			// current
																			// node,
																			// and
																			// if
																			// that
																			// node
																			// also
																			// isn't
																			// on
																			// closedlist
																			// (meaning
																			// we
																			// haven't
																			// yet
																			// looked
																			// at
																			// it)
					gCosts[i] = gCosts[current] + graph[current][i]; // we set
																		// the g
																		// cost
																		// of
																		// all
																		// adjacent
																		// nodes
																		// to
																		// the g
																		// cost
																		// of
																		// current
																		// + the
																		// weight
																		// of
																		// the
																		// edge
																		// between
																		// current
																		// and
																		// that
																		// adjacent
																		// node
					cameFrom[i] = current; // stores the path that we took to
											// reach the adjacent nodes so we
											// can later reconstruct the full
											// path to reach d
					openList.add(i); // we add the adjacent nodes to openlist,
										// as we will look at the smallest one
										// of them the next time through
				}
			}
		}
		System.out.print("No path is present from " + s + " to " + d + ", disconnect exists in graph. :^");
		// if the while list ever actually ends, it means that there is no
		// possible path, as all nodes that were reachable from the start node
		// have been considered (none of them are in openlist), and none of them
		// was the destination node
		return 0;
	}

	private int[] heuristic(int[][] graph, int d) { // method that uses our
													// heuristic to return an
													// int array containing each
													// node's h cost
		// Our heuristic is fairly simple. It starts at d, which it assigns an h
		// cost of 0. It then assigns all adjacent nodes an h cost of 1. It then
		// assigns an h value of 2 to all nodes that were adjacent to those
		// nodes and that did not already have an h value. It continues this
		// pattern, assigning to the nodes adjacent to the previous nodes which
		// do not already have an h value an h value one higher than that of the
		// previous nodes.
		ArrayList<Integer> thisLevel = new ArrayList<>(); // thisLevel is a list
															// of the nodes we
															// are currently
															// finding adjacent
															// nodes to
		ArrayList<Integer> nextLevel = new ArrayList<>(); // nextLevel is a list
															// of the nodes that
															// we will find
															// adjacent nodes to
															// the next time we
															// run through the
															// do while loop

		int dim = graph[0].length; // dim is the number of nodes
		boolean[] finished = new boolean[dim]; // finished is an int array which
												// stores whether or not we have
												// already assigned an h value
												// to each node
		int[] hVals = new int[dim]; // an int array which stores the actual h
									// costs of each node

		for (int i = 0; i < dim; i++) { // initializes all h costs to 0 and sets
										// all nodes as not finished
			hVals[i] = 0;
			finished[i] = false;
		}
		finished[d] = true; // the first time through, we only are finding
							// adjacent nodes to d, so d is added to thisLevel
							// and set to finished.
		thisLevel.add(d);
		int count = 1; // this is the number that we assign as the h value each
						// time we run through, it will be incremented each time
						// we move on to the next set of adjacent nodes.
		do {
			for (int i = 0; i < thisLevel.size(); i++) { // iterates once for
															// each node in
															// thisLevel
				for (int j = 0; j < dim; j++) { // iterates through all nodes
					if (graph[thisLevel.get(i)][j] != 0 && !finished[j]) { // if
																			// a
																			// node
																			// hasn't
																			// already
																			// been
																			// assigned
																			// an
																			// h
																			// value
																			// (meaning
																			// it
																			// is
																			// not
																			// finished)
																			// and
																			// is
																			// adjacent
																			// to
																			// i,
						hVals[j] = count; // we set the h cost of that node to
											// count
						nextLevel.add(j); // and we add it to nextLevel
						finished[j] = true; // and mark it as finished
					}
				}
			}
			count++; // increment count
			thisLevel.clear(); // we reset thisLevel and then add all elements
								// of nextLevel to it, as those are the nodes we
								// will be looking at during the next run
								// through.
			thisLevel.addAll(nextLevel);
			nextLevel.clear(); // we clear nextLevel in preparation for the next
								// run through
		} while (!thisLevel.isEmpty()); // if thisLevel is empty, then there are
										// no more nodes to consider and we are
										// done
		return hVals; // we return the int array containing all of the h costs
	}
}