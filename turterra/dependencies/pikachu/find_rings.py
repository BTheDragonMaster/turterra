#!/usr/bin/env python

from pprint import pprint
from copy import copy
import pikachu

class Graph(Structure):
    def __init__(self, molecule):
        super().__init__(molecule.graph, molecule.bonds, molecule.bond_lookup)
        self.time = 0

    def get_sssr(self):
        pass

    def get_rings(self):
        adjacency_matrix = self.get_component_adjacency_matrix()
        if not adjacency_matrix:
            return None

        connected_components = self.get_graph_components(adjacency_matrix)

        for component in connected_components:



    def get_graph_components(self, adjacency_matrix):
        visited = {}
        components = []
        count = 0

        for atom in self.graph:
            visited[atom] = False

        for atom in self.graph:
            if not visited[atom]:
                component = []
                visited[atom] = True
                component.append(atom)
                count += 1
                self.dfs_components(atom, visited, adjacency_matrix, component)
                if len(component) > 1:
                    components.append(component)

        return components

    def dfs_components(self, atom, visited, adjacency_matrix, component):
        for neighbour in adjacency_matrix[atom]:
            is_adjacent = adjacency_matrix[atom][neighbour]

            if not is_adjacent or visited[neighbour] or atom == neighbour:
                continue

            visited[neighbour] = True
            component.append(neighbour)
            self.dfs_components(neighbour, visited, adjacency_matrix, component)


    def get_component_adjacency_matrix(self):
        adjacency_matrix = {}

        for atom_1 in self.graph:
            adjacency_matrix[atom_1] = {}
            for atom_2 in self.graph:
                adjacency_matrix[atom_1][atom_2] = 0

        for bond_name, bond in self.bonds:
            adjacency_matrix[bond.atom_1][bond.atom_2] = 1
            adjacency_matrix[bond.atom_2][bond.atom_1] = 1


        bridges = self.get_bridges()

        for atom_1, atom_2 in bridges:
            adjacency_matrix[atom_1][atom_2] = 0
            adjacency_matrix[atom_2][atom_1] = 0

        return adjacency_matrix


    def get_bridges(self):
        visited = {}
        disc = {}
        low = {}
        parent = {}
        adjacency_list = self.get_adjacency_list()
        bridges = []
        self.time = 0

        for atom in self.graph:
            visited[atom] = False
            disc[atom] = 0
            parent[atom] = None
            low[atom] = 0

        for atom in self.graph:
            if not visited[atom]:
                dfs_bridges(atom, visited, disc, low, parent, adjacency_list, bridges)

        return bridges

    def dfs_bridges(self, atom, visited, disc, low, parent, adjacency_list, bridges):
        visited[atom] = True
        self.time += 1
        disc[atom] = low[atom] = self.time

        for neighbour in self.graph[atom]:
            if not visited[neighbour]:
                parent[neighbour] = atom
                graph.find_bridges(neighbour, visited, disc, low, parent, adjacency_list, bridges)

                low[atom] = min(low[atom], low[neighbour])

                if low[neighbour] > disc[atom]:
                    bridges.append(self.bond_lookup[atom][neighbour])

                elif neighbour != parent[atom]:
                    low[atom] = min(low[atom], disc[neighbour])


class TestGraph:
    def __init__(self):
        self.graph = {1: [2, 7],
                      2: [1, 3, 6],
                      3: [2, 4],
                      4: [3, 5],
                      5: [4, 6],
                      6: [5, 2, 7],
                      7: [6, 1]}
        self.bond_lookup = {1: {2: 1,
                                7: 7},
                            2: {1: 1,
                                3: 2,
                                6: 8},
                            3: {2: 2,
                                4: 3},
                            4: {3: 3,
                                5: 4},
                            5: {4: 4,
                                6: 5},
                            6: {5: 5,
                                2: 8,
                                7: 7},
                            7: {6: 6,
                                1: 7}}


def find_ring_members(pid_1, pid_2):
    ring_members = []
    for vertex_1, neighbours in pid_1.items():
        for vertex_2, paths in neighbours.items():
            pid_1_path_nr = len(paths)
            pid_2_path_nr = len(pid_2[vertex_1][vertex_2])
            if pid_1_path_nr > 2:
                ring_members.append((vertex_1, vertex_2))
            if pid_1_path_nr == 1 and pid_2_path_nr > 1:
                ring_members.append((vertex_1, vertex_2))

    return ring_members

def find_sssr(graph):
    pass

def compute_distance_matrix_fw(graph):
    """
    Use Floyd-Warshall algorithm to compute the shortest paths between all vertice pairs in a graph

    """

    vertices = list(graph.graph.keys())
    shortest_paths = {}
    all_paths = {}
    pid_1 = {}
    pid_2 = {}
    for vertex_1 in vertices:
        shortest_paths[vertex_1] = {}
        pid_1[vertex_1] = {}
        pid_2[vertex_1] = {}
        all_paths[vertex_1] = {}
        for vertex_2 in vertices:
            pid_2[vertex_1][vertex_2] = []
            all_paths[vertex_1][vertex_2] = {}

            if vertex_1 == vertex_2:
                shortest_path = 0
                pid_1[vertex_1][vertex_2] = []

            elif vertex_1 in graph.graph[vertex_2]:
                shortest_path = 1
                pid_1[vertex_1][vertex_2] = [[graph.bond_lookup[vertex_1][vertex_2]]]
            else:
                shortest_path = float('inf')
                pid_1[vertex_1][vertex_2] = []

            shortest_paths[vertex_1][vertex_2] = shortest_path

    for vertex_k in vertices:
        for vertex_i in vertices:
            for vertex_j in vertices:
                path_length = shortest_paths[vertex_i][vertex_k] + shortest_paths[vertex_k][vertex_j]
                path = pid_1[vertex_i][vertex_k] + pid_1[vertex_k][vertex_j]
                if not path_length in all_paths[vertex_i][vertex_j]:
                    all_paths[vertex_i][vertex_j][path_length] = []

                if not path in all_paths[vertex_i][vertex_j][path_length]:
                    all_paths[vertex_i][vertex_j][path_length].append(path)

                if shortest_paths[vertex_i][vertex_j] > path_length:
                    shortest_paths[vertex_i][vertex_j] = path_length
                    if pid_1[vertex_i][vertex_j] and len(pid_1[vertex_i][vertex_j][0]) == path_length:
                        pid_1[vertex_i][vertex_j].append(path)
                    else:
                        pid_1[vertex_i][vertex_j] = [path]

    for vertex_1, neighbours in shortest_paths.items():
        for vertex_2, distance in neighbours.items():
            try:
                pid_2[vertex_1][vertex_2] = all_paths[vertex_1][vertex_2][distance + 1]
            except KeyError:
                pass


    return shortest_paths, pid_1, pid_2

if __name__ == "__main__":
    test_graph = TestGraph()
    shortest_paths, pid_1, pid_2 = compute_distance_matrix_fw(test_graph)
    pprint(shortest_paths)
    pprint(pid_1)
    pprint(pid_2)



