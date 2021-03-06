#!/usr/bin/env python
import pikachu
import find_cycles
import copy
from pprint import pprint

def get_cyclic_label(cycle_nr):
    if cycle_nr > 9:
        return '%' + str(cycle_nr)
    else:
        return str(cycle_nr)

def get_chiral_permutations(order):
    permutations = [tuple(order)]
    permutations.append((order[0], order[3], order[1], order[2]))
    permutations.append((order[0], order[2], order[3], order[1]))
    permutations.append((order[1], order[0], order[3], order[2]))
    permutations.append((order[1], order[2], order[0], order[3]))
    permutations.append((order[1], order[3], order[2], order[0]))
    permutations.append((order[2], order[0], order[1], order[3]))
    permutations.append((order[2], order[3], order[0], order[1]))
    permutations.append((order[2], order[1], order[3], order[0]))
    permutations.append((order[3], order[0], order[2], order[1]))
    permutations.append((order[3], order[1], order[0], order[2]))
    permutations.append((order[3], order[2], order[1], order[0]))

    return permutations

def determine_chirality(order, chirality):
    original_order = tuple(order)
    chiral_permutations = get_chiral_permutations(original_order)
    order.sort(key=lambda x: x.nr)
    new_order = tuple(order)
    if new_order in chiral_permutations:
        if chirality == 'counterclockwise':
            new_chirality = 'counterclockwise'
        else:
            new_chirality = 'clockwise'
    else:
        if chirality == 'counterclockwise':
            new_chirality = 'clockwise'
        else:
            new_chirality = 'counterclockwise'

    return new_chirality

class GraphToSmiles:
    def __init__(self, structure):
        self.original_structure = structure
        self.structure = copy.deepcopy(structure)


        self.remove_hydrogens()

        self.add_representations()

        self.find_branch_points()
        self.find_terminal_nodes()
        self.find_cycles()

       # print('branch points', self.branch_points)
      #  print('terminal nodes', self.terminal_nodes)
      #  print('cycles', self.cycles)

        self.make_smiles_components()
        self.find_original_atom_indices()
        self.resolve_chiral_centres()
        print(''.join(self.components))

    def is_numerical_component(self, component):
        for character in component:
            if character not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '%']:
                return False

        return True

    def find_original_atom_indices(self):
        for i in range(len(self.components) - 1, -1, -1):
            component = self.components[i]
            if self.is_numerical_component(component):
                atom_to_adjust = None
                for atom, index in self.atom_to_index.items():
                    if index == i:
                        atom_to_adjust = atom

                current_index = i
                while self.is_numerical_component(self.components[current_index]) or component in ('=', '#'):
                    current_index -= 1

                self.atom_to_index[atom_to_adjust] = current_index




    def resolve_chiral_centres(self):
        for atom in self.original_structure.graph:
            if atom.chiral:
                neighbours = atom.neighbours

                indices_and_atoms = []

                cyclic_neighbours = []

                for cycle_nr in self.atom_to_cycle_nr[atom]:
                    cyclic_atoms = self.cycle_nr_to_atoms[cycle_nr]

                    if cyclic_atoms[0] == atom:
                        other_atom = cyclic_atoms[1]
                    elif cyclic_atoms[1] == atom:
                        other_atom = cyclic_atoms[0]
                    else:
                        raise Exception("Cycle should contain the current atom.")

                    cyclic_neighbours.append((other_atom, cycle_nr))

                index = None


                for neighbour in neighbours:
                    if neighbour.type == 'H':
                        index = self.atom_to_index[atom]
                    else:
                        if not neighbour in [atom_and_cycle[0] for atom_and_cycle in cyclic_neighbours]:
                            index = self.atom_to_index[neighbour]
                        else:
                            for atom_and_cycle in cyclic_neighbours:
                                if neighbour == atom_and_cycle[0]:
                                    cycle_nr = atom_and_cycle[1]
                                    break


                            for i, component in enumerate(self.components):

                                if component == cycle_nr and i > self.atom_to_index[atom]:
                                    index = i

                                    break

                    indices_and_atoms.append((index, neighbour))

                atom_order = [atom for _,atom in sorted(indices_and_atoms)]
                chirality = determine_chirality(atom_order, atom.chiral)
                if chirality == 'counterclockwise':
                    chiral_symbol = '@'
                elif chirality == 'clockwise':
                    chiral_symbol = '@@'
                chiral_centre_index = self.atom_to_index[atom]
                self.components[chiral_centre_index] = self.components[chiral_centre_index].replace('X', chiral_symbol)

    def get_branch_levels(self):
        atom_to_branch = []
        index_to_atom = {}
        for atom, index in atom_to_index.items():
            atom_to_branch[atom] = None
            index_to_atom[index] = atom


        branch = 0
        for i, component in enumerate(self.components):
            if component == '(':
                branch += 1
            elif component == ')':
                branch -= 1
            else:
                if i in index_to_atom:
                    atom = index_to_atom[i]


    def make_smiles_components(self):
        self.components = []
        self.atom_to_index = {}
        self.cycle_nr_to_atoms = {}
        self.atom_to_cycle_nr = {}
        self.bonds_to_index = {}

        is_branched = {}
        cycle_nr = 0

        for atom in self.structure.graph:
            self.atom_to_index[atom] = None
            self.atom_to_cycle_nr[atom] = []
            is_branched[atom] = False

        if self.terminal_nodes:
            first_atom = list(self.terminal_nodes)[0]
        elif self.branch_points:
            first_atom = list(self.branch_points)[0]
        else: #this happens when the entire graph is cyclic
            first_atom = list(self.structure.graph.keys())[0]

        self.atom_to_index[first_atom] = 0

        current_atom = first_atom
        self.components.append(self.representations[current_atom])

        working_graph = copy.deepcopy(self.structure)
        atoms_left = set(working_graph.graph.keys())
        atoms_added = set()
        atoms_added.add(current_atom)

        while atoms_left:

            if working_graph.graph[current_atom]:
                cyclic = False
                next_atom = working_graph.graph[current_atom][0]
                if self.atom_to_index[next_atom] != None:
                    cyclic = True
                    cycle_nr += 1
                    cyclic_label = get_cyclic_label(cycle_nr)

                bond_type = working_graph.bond_lookup[current_atom][next_atom].type
                bond_symbol = pikachu.BOND_PROPERTIES.bond_type_to_symbol[bond_type]

                if bond_type == 'single' and current_atom.aromatic and next_atom.aromatic:
                    bond_symbol = '-'

                if cyclic:

                    if bond_symbol:
                        cyclic_label_idx_1 = self.atom_to_index[next_atom]
                        self.add_insert([bond_symbol, cyclic_label], cyclic_label_idx_1)
                        cyclic_label_idx_2 = self.atom_to_index[current_atom]
                        self.add_insert([bond_symbol, cyclic_label], cyclic_label_idx_2)
                        offset = 2
                    else:
                        cyclic_label_idx_1 = self.atom_to_index[next_atom]
                        self.add_insert([cyclic_label], cyclic_label_idx_1)
                        cyclic_label_idx_2 = self.atom_to_index[current_atom]
                        self.add_insert([cyclic_label], cyclic_label_idx_2)
                        offset = 1


                    self.atom_to_index[next_atom] += offset
                    self.atom_to_index[current_atom] += offset

                    self.cycle_nr_to_atoms[cyclic_label] = (current_atom, next_atom)
                    self.atom_to_cycle_nr[current_atom].append(cyclic_label)
                    self.atom_to_cycle_nr[next_atom].append(cyclic_label)

                else:
                    break_point = self.atom_to_index[current_atom]


                    if is_branched[current_atom]:
                        if bond_symbol:
                            insert = ["(", bond_symbol, self.representations[next_atom], ")"]
                            offset = 3
                        else:
                            insert = ["(", self.representations[next_atom], ")"]
                            offset = 2
                    else:
                        if bond_symbol:
                            insert = [bond_symbol, self.representations[next_atom]]
                            offset = 2
                        else:
                            insert = [self.representations[next_atom]]
                            offset = 1

                    self.add_insert(insert, break_point)
                    self.atom_to_index[next_atom] = break_point + offset


                if not cyclic:
                    is_branched[current_atom] = True

                atoms_added.add(next_atom)

                working_graph.remove_bond_between_atoms(current_atom, next_atom)
                current_atom = next_atom

            else:
                del working_graph.graph[current_atom]
                atoms_left.remove(current_atom)
                for atom in working_graph.graph:
                    if atom in atoms_added:
                        current_atom = atom
                        break

        print(''.join(self.components))


    def get_neighbour_nr(self, node):
        neighbours = self.structure.graph[node]
        neighbour_nr = 0
        for neighbour in neighbours:
            if neighbour.type != 'H':
                neighbour_nr += 1

        return neighbour_nr

    def remove_hydrogens(self):
        for atom in list(self.structure.graph.keys()):
            if atom.type == 'H' and atom.charge == 0:
                neighbour = atom.neighbours[0]
                if neighbour.type in ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I'] \
                   and neighbour.charge == 0 and not neighbour.chiral and not neighbour.pyrrole:
                    bond = atom.bonds[0]
                    del self.structure.bonds[bond.nr]
                    del self.structure.bond_lookup[atom][neighbour]
                    del self.structure.bond_lookup[neighbour][atom]
                    
                    self.structure.graph[neighbour].remove(atom)
                    self.structure.graph[atom].remove(neighbour)

                    if not self.structure.graph[neighbour]:
                        del self.structure.graph[neighbour]

                    if not self.structure.graph[atom]:
                        del self.structure.graph[atom]

    def get_bond_symbol(self, last_node, next_node):
        bond_type = self.structure.bond_lookup[last_node][next_node].type
        symbol = pikachu.BOND_PROPERTIES.bond_type_to_symbol[bond_type]
        return symbol

    def remove_explicit_hydrogen(self, hydrogen):
        neighbour = hydrogen.neighbours[0]

        bond = hydrogen.bonds[0]

        del self.structure.bonds[bond.nr]
        del self.structure.bond_lookup[hydrogen][neighbour]
        del self.structure.bond_lookup[neighbour][hydrogen]

        self.structure.graph[neighbour].remove(hydrogen)

        del self.structure.graph[hydrogen]


    def remove_bonds_and_nodes(self, last_node, next_node):
        bond = self.structure.bond_lookup[last_node][next_node]
        
        del self.structure.bonds[bond.nr]
        del self.structure.bond_lookup[last_node][next_node]
        del self.structure.bond_lookup[next_node][last_node]

        if next_node in self.structure.graph[last_node]:
            self.structure.graph[last_node].remove(next_node)
        if last_node in self.structure.graph[next_node]:
            self.structure.graph[next_node].remove(last_node)

        if not self.structure.graph[last_node]:
            del self.structure.graph[last_node]

        if not self.structure.graph[next_node]:
            del self.structure.graph[next_node]

    def make_simplified_graph_from_collapsed(self):
        self.simplified_graph = {}
        for atom, next_atom_and_edges in self.collapsed_graph.items():
            self.simplified_graph[atom] = []
            for next_atom, next_edge, edge_atoms in next_atom_and_edges:
                self.simplified_graph[atom].append(next_atom)

    def make_branch_lookup(self):
        branch_lookup = {}
        for node in self.simplified_graph:
            branch_lookup[node] = 0

        return branch_lookup

    def make_index_dict(self):
        index_dict = {}
        
        for node in self.simplified_graph:
            index_dict[node] = None

        return index_dict

    def adjust_atom_indices(self, insert, break_point):
        index_jump = len(insert)

        for atom, index in self.atom_to_index.items():
            if index != None and index > break_point:
                self.atom_to_index[atom] += index_jump

    def add_insert(self, insert, break_point):
        half_1 = self.components[:break_point + 1]
        half_2 = self.components[break_point + 1:]

        self.adjust_atom_indices(insert, break_point)

        self.components = half_1 + insert + half_2               



    def add_representations(self):
        self.representations = {}
        for atom in self.structure.graph:
            if atom.type != 'H':
                if atom.aromatic:
                    self.representations[atom] = atom.type.lower()
                else:
                    self.representations[atom] = atom.type

        self.add_chiral_placeholders()
        self.add_hydrogen_representations()
        self.add_charge_representations()

    def add_chiral_placeholders(self):
        for atom in self.structure.graph:
            if atom.chiral:
                self.representations[atom] = f'[{self.representations[atom]}X]'

    def add_charge_representations(self):
        for atom in self.structure.graph:
            if atom.charge != 0:
                if charge == 1:
                    charge_string = '+'
                elif charge == -1:
                    charge_string = '-'
                elif charge_string > 1:
                    charge_string = '+%d' % charge
                elif charge_string < -1:
                    charge_string = '-%d' % charge
                    
                representation = self.representations[atom]
                if representation[-1] == ']':
                    self.representations[atom] = representation[:-1] + charge_string + ']'
                else:
                    self.representations[atom] = '[' + representation + charge_string + ']'

    def count_hydrogens(self):
        hydrogen_counts = {}
        atom_to_hydrogens = {}

        for atom, neighbours in self.structure.graph.items():
            hydrogen_counts[atom] = 0
            atom_to_hydrogens[atom] = []
            for neighbour in neighbours:
                if neighbour.type == 'H':
                    hydrogen_counts[atom] += 1
                    atom_to_hydrogens[atom].append(neighbour)

        return hydrogen_counts, atom_to_hydrogens

    def add_hydrogen_representations(self):
        hydrogen_counts, atom_to_hydrogens = self.count_hydrogens()
        
        for atom, count in hydrogen_counts.items():
            if count > 0:
                if count > 1:
                    hydrogen_string = f'H{count}'
                elif count == 1:
                    hydrogen_string = 'H'



                representation = self.representations[atom]
                if representation[-1] == ']':
                    self.representations[atom] = representation[:-1] + hydrogen_string + ']'
                else:
                    self.representations[atom] = '[' + representation + hydrogen_string + ']'

                explicit_hydrogens = atom_to_hydrogens[atom]
                for hydrogen in explicit_hydrogens:
                    self.remove_explicit_hydrogen(hydrogen)

    def find_cycles(self):
        cycles = find_cycles.Cycles(self.structure)

        self.cycles = cycles.find_minimal_cycles()

    def find_branch_points(self):
        self.branch_points = set()
        for node in self.structure.graph:
            if self.is_branch_point(node):
                self.branch_points.add(node)

    def find_terminal_nodes(self):
        self.terminal_nodes = set()
        for node in self.structure.graph:
            if self.is_terminal_node(node):
                self.terminal_nodes.add(node)

    def find_cyclic_pairs(self):
        self.cyclic_pairs = set()
        self.cyclic_nodes = set()
        for cycle in self.cycles:
            for i, atom in enumerate(cycle):
                previous_atom = cycle[i - 1]
                self.cyclic_nodes.add(previous_atom)
                self.cyclic_nodes.add(atom)
                pair = tuple(sorted((previous_atom, atom), key = lambda x: x.nr))
                self.cyclic_pairs.add(pair)

    def is_terminal_node(self, node):
        if self.get_neighbour_nr(node) == 1:
            return True
        else:
            return False

    def is_branch_point(self, node):
        if self.get_neighbour_nr(node) >= 3:
            return True
        else:
            return False




if __name__ == "__main__":
    smiles = 'c1ccc2c(c1)c(c[nH]2)C[C@@H](C(=O)O)N'
 #   smiles = 'C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N'
 #   smiles = 'c1(C)c(C)c(C)c(C)c(C)c1(C)'
 #   smiles = 'c1(c(c(c(c(c1C)C)C)C)C)C'
    smiles = 'c12ccccc1cccc2'
 #   smiles = 'c1ccccc1'
    smiles = 'C1=CC=C(C=C1)C[C@@H](C(=O)O)N'
    smiles = 'CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C'
  #  smiles = 'CCCCCCCCCC(=O)NC(CC1=CNC2=CC=CC=C21)C(=O)NC(CC(=O)N)C(=O)NC(CC(=O)O)C(=O)NC3C(OC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)C(C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C'
    smiles = 'CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O'
  #  smiles = 'CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O'
#    smiles = 'c1ccccc1-c2ccccc2'
  #  smiles = "C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6=C7C8=C1C1=C8C%10=C%10C%11=C2C2=C3C3=C4C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%12C5=C%11C4=C3C3=C5C(=C81)C%10=C23"
    structure = pikachu.Smiles(smiles).smiles_to_structure()
    collapsed_structure = GraphToSmiles(structure)
 #   pprint(collapsed_structure.collapsed_graph)
 #   pprint(collapsed_structure.simplified_graph)
  #  pprint(collapsed_structure.representations)
 #   pprint(collapsed_structure.cycles)
#    pprint(collapsed_structure.simplified_graph)
#    pprint(collapsed_structure.representations)
#    pprint(collapsed_structure.node_to_edge_dict)

    smiles = 'CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C'
    smiles = 'c1ccc2c(c1)c(c[nH]2)C[C@@H](C(=O)O)N'
    structure = pikachu.Smiles(smiles).smiles_to_structure()
    kekule_structure = structure.kekulise()
    GraphToSmiles(kekule_structure)
    GraphToSmiles(structure)
    
    

        
