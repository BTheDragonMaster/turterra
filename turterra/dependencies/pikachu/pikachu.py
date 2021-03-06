import logging
from pprint import pprint
import copy
from collections import OrderedDict, defaultdict
from typing import *
from dataclasses import dataclass
import sys
from math import sqrt

import turterra.dependencies.pikachu.find_cycles as find_cycles
from turterra.dependencies.pikachu.pikachu_errors import SmilesError
from turterra.dependencies.pikachu.math_functions import *



sys.setrecursionlimit(100000)


def read_smiles(smiles_string):
    if smiles_string:
        try:
            smiles = Smiles(smiles_string)
            structure = smiles.smiles_to_structure()
        except SmilesError as e:
            print(f'Error parsing "{smiles_string}": {e.message}')
            return

        return structure

    return

def compare_matches(match_1, match_2):
    matching = True
    for key in match_1:
        if not key in match_2:
            matching = False
            break

        if match_1[key] != match_2[key]:
            matching = False
            break

    return matching

def compare_all_matches(matches):
    matching_pairs = set([])
    
    for i, match_1 in enumerate(matches):
        for j, match_2 in enumerate(matches):
            if i != j:
                if compare_matches(match_1, match_2):
                    matching_pairs.add(tuple(sorted([i, j])))

    matches_to_remove = set([])

    for matching_pair in matching_pairs:
        matches_to_remove.add(matching_pair[1])

    matches_to_remove = sorted(list(matches_to_remove), reverse = True)
    
    for match_to_remove in matches_to_remove:
        del matches[match_to_remove]

def check_aromatic(atom_set):
    aromatic = True
    for atom in atom_set:
        if atom.hybridisation == 'sp2':
            pass
        else:
            aromatic = False
            break

    if aromatic:
        pi_electron_nr = 0
        for atom in atom_set:
            for orbital_name in atom.valence_shell.orbitals:
                orbital = atom.valence_shell.orbitals[orbital_name]
                if orbital.orbital_type == 'p':
                    for electron in orbital.electrons:
                        if electron.atom == atom:
                            pi_electron_nr += 1

        if not pi_electron_nr % 4 == 2:
            aromatic = False

    return aromatic

def check_five_ring(atom_set):
    assert len(atom_set) == 5

    sp2_hybridised = []
    sp3_hybridised_lone_pair = []

    aromatic = False
    heteroatom = None

    for atom in atom_set:
        if atom.hybridisation == 'sp2':
            sp2_hybridised.append(atom)
        elif atom.hybridisation == 'sp3':
            if atom.calc_electron_pair_nr() > 0:
                sp3_hybridised_lone_pair.append(atom)

    if len(sp2_hybridised) == 4 and len(sp3_hybridised_lone_pair) == 1:
        
        pi_electron_nr = 0
        for atom in sp2_hybridised:
            for orbital_name in atom.valence_shell.orbitals:
                orbital = atom.valence_shell.orbitals[orbital_name]
                if orbital.orbital_type == 'p':
                    for electron in orbital.electrons:
                        if electron.atom == atom:
                            pi_electron_nr += 1
            if pi_electron_nr % 4 == 0:
                aromatic = True
                heteroatom = sp3_hybridised_lone_pair[0]

    return aromatic, heteroatom

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


def check_same_chirality(atom_1, atom_2, match):

    equivalent_atom_list = []
    for atom in atom_1.neighbours:
        if atom.type == 'H':
            for atom_2 in atom_2.neighbours:
                if atom_2.type == 'H':
                    equivalent_atom_list.append(atom_2)
                    break
        else:
            equivalent_atom_list.append(match[atom])

    permutation = equivalent_atom_list[:]

    if len(equivalent_atom_list) != 4:
        lone_pairs = atom_2.lone_pairs

        try:
            assert len(equivalent_atom_list) + len(lone_pairs) == 4
        except AssertionError:
            raise SmilesError('chiral centre')

        permutation += lone_pairs

    chiral_permutations = get_chiral_permutations(permutation)


    if atom_1.chiral == atom_2.chiral:
        if tuple(permutation) in chiral_permutations:
            return True
        else:
            return False
    else:
        if tuple(permutation) in chiral_permutations:
            return False
        else:
            return True
        
class BondProperties:
    """
    A class storing various properties of bonds

    Attributes
    ----------
    bond_type_to_weight: dict of {bond type: bond_weight, ->}
        with bond type (str) the type of covalent bond, and bond weight (int)
        the number of times that bond is counted in determining if an atom is
        in its excited state or not
    bond_type_to_symbol: dict of {bond type: SMILES symbol, ->}
        with bond type (str) the type of covalent bond, and SMILES symbol the
        text symbol used in SMILES strings to represent that bond
    bond_type_to_p_orbitals: dict of {bond type: p orbitals, ->}
        with bond type (str) the type of covalent bond, and p orbitals (int)
        the number of p orbitals involved in the formation of that bond
    """

    type_to_dash2d_input = {'single': 1,
                            'double': 2,
                            'triple': 3,
                            'quadruple': 4}

    bond_type_to_weight = {'single': 1,
                           'double': 2,
                           'triple': 3,
                           'quadruple': 4,
                           'aromatic': 1}

    bond_type_to_symbol = {'single': '',
                           'double': '=',
                           'triple': '#',
                           'aromatic': ''}
    
    bond_type_to_p_orbitals = {'single': 0,
                               'double': 1,
                               'triple': 2,
                               'quadruple': 3,
                               'aromatic': 1}

BOND_PROPERTIES = BondProperties()        

class AtomProperties:
    """
    A class storing various properties of atoms

    Attributes
    ----------
    element_to_valences: dict of {element: [valence, ->], ->}
        with element (str) an abbreviation of an element of the periodic table,
        and valence (int) a possible number of covalent bonds that element can
        make
    element_to_atomic_nr: dict of {element: atomic nr, ->}
        with element (str) an abbreviation of an element of the periodic table,
        and atomic nr (int) its atomic number
    element_to_radius: dict of {element: radius, ->}
        with element (str) an abbreviation of an element of the periodic table,
        and radius (float) the atom radius of that element in Angstrom.
    element_to_valence_electrons: dict of {element: valence electrons}
        with element (str) an abbreviation of an element of the periodic table,
        and valence electrons (int) the number of electrons in the valence
        shell of that element
    element_to_shell_nr: dict of {element: shell nr, ->}
        with element (str) an abbreviation of an element of the periodic table,
        and shell nr (int) the number of electron shells of that element
    shell_to_electron_nr: dict of {shell: electron nr, ->}
        with shell (int) the index of an electron shell counting from the
        nucleus, and electron nr (int) the number of electrons in that shell
    shell_to_orbitals: dict of {shell: [orbital, ->], ->}
        with shell (int) the index of an electron shell counting from the
        nucleus, and orbital (str) an electron orbital in that shell
    orbital_order: tuple of (orbital, ->)
        with the orbitals (str) placed in order of them being filled with
        electrons in nature
    orbital_type_to_orbital_number: dict of {orbital type: orbital number, ->}
        with orbital type (str) the type of orbital (s, p, d, f, g), and
        orbital number the number of orbitals of that type an atom has
    
    """
    element_to_valences = {'C': [4],
                           'O': [2],
                           'N': [3],
                           'P': [3, 5],
                           'S': [2, 4, 6],
                           'Cl': [1, 7],
                           'Br': [1, 7],
                           'F': [1],
                           'H': [1],
                           'I': [1, 7],
                           'Na': [1],
                           'Fe': [2, 3],
                           '*': [1]}

    element_to_atomic_nr = {'H': 1,
                            'He': 2,
                            'Li': 3,
                            'Be': 4,
                            'B': 5,
                            'C': 6,
                            'N': 7,
                            'O': 8,
                            'F': 9,
                            'Ne': 10,
                            'Na': 11,
                            'Mg': 12,
                            'Al': 13,
                            'Si': 14,
                            'P': 15,
                            'S': 16,
                            'Cl': 17,
                            'Ar': 18,
                            'K': 19,
                            'Ca': 20,
                            'Sc': 21,
                            'Ti': 22,
                            'V': 23,
                            'Cr': 24,
                            'Mn': 25,
                            'Fe': 26,
                            'Co': 27,
                            'Ni': 28,
                            'Cu': 29,
                            'Zn': 30,
                            'Ga': 31,
                            'Ge': 32,
                            'As': 33,
                            'Se': 34,
                            'Br': 35,
                            'Kr': 36,
                            'I': 53}

    element_to_radius = {'H': 0.37,
                        'He': 0.32,
                        'Li': 1.34,
                        'Be': 0.90,
                        'B': 0.82,
                        'C': 0.77,
                        'N': 0.75,
                        'O': 0.73,
                        'F': 0.71,
                        'Ne': 0.69,
                        'Na': 1.54,
                        'Mg': 1.30,
                        'Al': 1.18,
                        'Si': 1.11,
                        'P': 1.06,
                        'S': 1.02,
                        'Cl': 0.99,
                        'Ar': 0.97,
                        'K': 1.96,
                        'Ca': 1.74,
                        'Sc': 1.44,
                        'Ti': 1.36,
                        'V': 1.25,
                        'Cr': 1.27,
                        'Mn': 1.39,
                        'Fe': 1.25,
                        'Co': 1.26,
                        'Ni': 1.21,
                        'Cu': 1.38,
                        'Zn': 1.31,
                        'Ga': 1.26,
                        'Ge': 1.22,
                        'As': 1.19,
                        'Se': 1.16,
                        'Br': 1.14,
                        'Kr': 1.10,
                        'I': 1.33}


    element_to_valence_electrons = {'H': 1,
                                    'C': 4,
                                    'O': 6,
                                    'N': 5,
                                    'P': 5,
                                    'S': 6,
                                    'Cl': 7,
                                    'Br': 7,
                                    'F': 7,
                                    'I': 7,
                                    'Na': 1}

    element_to_shell_nr = {'H': 1,
                           'C': 2,
                           'N': 2,
                           'O': 2,
                           'F': 2,
                           'P': 3,
                           'S': 3,
                           'Cl': 3,
                           'Br': 4,
                           'I': 5,
                           'Na': 3,
                           'Fe': 4,
                           '*': 1}

    shell_to_electron_nr = {1: 2,
                            2: 8,
                            3: 18,
                            4: 32,
                            5: 32}
                           

    shell_to_orbitals = {1: ['1s'],
                         2: ['2s', '2p1', '2p2', '2p3'],
                         3: ['3s', '3p1', '3p2', '3p3', '3d1', '3d2', '3d3', '3d4', '3d5'],
                         4: ['4s', '4p1', '4p2', '4p3', '4d1', '4d2', '4d3', '4d4', '4d5', '4f1',
                             '4f2', '4f3', '4f4', '4f5', '4f6', '4f7'],
                         5: ['5s', '5p1', '5p2', '5p3', '5d1', '5d2', '5d3', '5d4', '5d5', '5f1',
                             '5f2', '5f3', '5f4', '5f5', '5f6', '5f7']}

    orbital_order = ('1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s', '4d', '5p', '6s', '4f', '5d', '6p', '7s', '5f', '6d', '7p', '8s', '5g', '6f', '7d', '8p', '9s')

    
    orbital_type_to_orbital_number = {'s': 1,
                                      'p': 3,
                                      'd': 5,
                                      'f': 7,
                                      'g': 9}
    
ATOM_PROPERTIES = AtomProperties()

class Structure:
    """
    Class to store a molecular structure

    Attributes
    ----------
    graph: dict of {atom: [atom, ->], ->}
        with each atom (Atom object) pointing to a list of atoms
        covalently bonded to that atom
    bonds: dict of {bond nr: bond, ->}
        with bond nr (int) an arbitrarily chosen and unique number for a bond,
        and bond (Bond object) the bond associated with that number
    bond_lookup: dict of {atom: {atom: bond, ->}, ->}
        with each pair of atoms (Atom object) pointing to the bond (Bond)
        between those atoms
    """

    def __init__(self, graph = None, bonds = None, bond_lookup = None):
        if graph:
            self.graph = graph
        else:
            self.graph = {}
        if bonds:
            self.bonds = bonds
        else:
            self.bonds = {}

        if bond_lookup:
            self.bond_lookup = bond_lookup
        else:
            self.bond_lookup = {}

    def get_next_in_ring(self, ring, current_atom, previous_atom):
        neighbours = self.graph[current_atom]

        for neighbour in neighbours:
            for member in ring.members:
                if neighbour == member:
                    if previous_atom != neighbour:
                        return neighbour

        return None

    def to_dash_molecule2d_input(self):
        nodes = []
        links = []

        kekulised_structure = self.kekulise()

        for atom in kekulised_structure.graph:
            atom_dict = {}
            atom_dict['id'] = atom.nr
            atom_dict['atom'] = atom.type
            nodes.append(atom_dict)

        for bond_nr, bond in kekulised_structure.bonds.items():
            assert bond.type != 'aromatic'
            bond_dict = {}

            bond_dict['id'] = bond_nr
            bond_dict['source'] = bond.atom_1.nr
            bond_dict['target'] = bond.atom_2.nr
            bond_dict['bond'] = BOND_PROPERTIES.bond_type_to_weight[bond.type]
            links.append(bond_dict)

        dash_molecule2d_input = {'nodes': nodes, 'links': links}
        return dash_molecule2d_input

    def set_double_bond_chirality(self):

        for bond_nr, bond in self.bonds.items():

            #iterate over all double bonds

            if bond.type == 'double' or bond.type == 'triple':

                #double bonds can only be chiral if they have exactly three bonds

                if len(bond.atom_1.bonds) == 3 and len(bond.atom_2.bonds) == 3:


                    #define atoms adjacent to the atoms involved in the double bond
                    #also keep track of the chiral symbol that defines these bonds

                    atom_1_1 = None
                    atom_1_2 = None
                    chiral_1_1 = None
                    chiral_1_2 = None

                    atom_2_1 = None
                    atom_2_2 = None
                    chiral_2_1 = None
                    chiral_2_2 = None

                    for bond_1 in bond.atom_1.bonds:


                        if bond_1.type == 'single':
                            if bond.atom_1 == bond_1.atom_1:
                                if bond_1.chiral_symbol == '/':
                                    direction = 'up'
                                elif bond_1.chiral_symbol == '\\':
                                    direction = 'down'
                                else:
                                    direction = None

                                if not atom_1_1:
                                    atom_1_1 = bond_1.atom_2
                                    chiral_1_1 = direction
                                else:
                                    atom_1_2 = bond_1.atom_2
                                    chiral_1_2 = direction

                            elif bond.atom_1 == bond_1.atom_2:

                                if bond_1.chiral_symbol == '/':
                                    direction = 'down'
                                elif bond_1.chiral_symbol == '\\':
                                    direction = 'up'
                                else:
                                    direction = None

                                if not atom_1_1:
                                    atom_1_1 = bond_1.atom_1
                                    chiral_1_1 = direction
                                else:
                                    atom_1_2 = bond_1.atom_1
                                    chiral_1_2 = direction


                    for bond_2 in bond.atom_2.bonds:

                        if bond_2.type == 'single':
                            if bond.atom_2 == bond_2.atom_1:
                                if bond_2.chiral_symbol == '/':
                                    direction = 'up'
                                elif bond_2.chiral_symbol == '\\':
                                    direction = 'down'
                                else:
                                    direction = None

                                if not atom_2_1:
                                    atom_2_1 = bond_2.atom_2
                                    chiral_2_1 = direction
                                else:
                                    atom_2_2 = bond_2.atom_2
                                    chiral_2_2 = direction

                            elif bond.atom_2 == bond_2.atom_2:
                                if bond_2.chiral_symbol == '/':
                                    direction = 'down'
                                elif bond_2.chiral_symbol == '\\':
                                    direction = 'up'
                                else:
                                    direction = None

                                if not atom_2_1:
                                    atom_2_1 = bond_2.atom_1
                                    chiral_2_1 = direction
                                else:
                                    atom_2_2 = bond_2.atom_1
                                    chiral_2_2 = direction



                    chiral_1 = False
                    chiral_2 = False

                    if chiral_1_1 or chiral_1_2:
                        chiral_1 = True

                    if chiral_2_1 or chiral_2_2:
                        chiral_2 = True

                    if chiral_1 and chiral_2:
                        if chiral_1_1 == chiral_1_2:

                            raise SmilesError('chiral double bond')
                        if chiral_2_2 == chiral_2_1:

                            raise SmilesError('chiral double bond')

                        first_atom = None
                        second_atom = None

                        if chiral_1_1:
                            first_atom = atom_1_1
                            first_other_atom = atom_1_2
                            first_chiral_symbol = chiral_1_1

                        else:
                            first_atom = atom_1_2
                            first_other_atom = atom_1_1
                            first_chiral_symbol = chiral_1_2

                        if chiral_2_1:
                            second_atom = atom_2_1
                            second_other_atom = atom_2_2
                            second_chiral_symbol = chiral_2_1
                        else:
                            second_atom = atom_2_2
                            second_other_atom = atom_2_1
                            second_chiral_symbol = chiral_2_2

                        bond.chiral_dict[first_atom] = {}
                        bond.chiral_dict[first_other_atom] = {}
                        bond.chiral_dict[second_atom] = {}
                        bond.chiral_dict[second_other_atom] = {}

                        if first_chiral_symbol == second_chiral_symbol:
                            bond.chiral_dict[first_atom][second_atom] = 'cis'
                            bond.chiral_dict[second_atom][first_atom] = 'cis'

                            bond.chiral_dict[first_other_atom][second_other_atom] = 'cis'
                            bond.chiral_dict[second_other_atom][first_other_atom] = 'cis'

                            bond.chiral_dict[first_atom][second_other_atom] = 'trans'
                            bond.chiral_dict[second_other_atom][first_atom] = 'trans'

                            bond.chiral_dict[first_other_atom][second_atom] = 'trans'
                            bond.chiral_dict[second_atom][first_other_atom] = 'trans'

                        else:
                            bond.chiral_dict[first_atom][second_atom] = 'trans'
                            bond.chiral_dict[second_atom][first_atom] = 'trans'

                            bond.chiral_dict[first_other_atom][second_other_atom] = 'trans'
                            bond.chiral_dict[second_other_atom][first_other_atom] = 'trans'

                            bond.chiral_dict[first_atom][second_other_atom] = 'cis'
                            bond.chiral_dict[second_other_atom][first_atom] = 'cis'

                            bond.chiral_dict[first_other_atom][second_atom] = 'cis'
                            bond.chiral_dict[second_atom][first_other_atom] = 'cis'

                        bond.chiral = True


    def get_atom_representations(self):
        atoms = sorted(self.graph.keys())
        atom_to_repr = {}
        for atom in atoms:
            atom_repr = atom.type
            neighbours = self.graph[atom]
            if atom.chiral:
                if atom.get_hydrogen_nr == 1:
                    neighbours = sorted(neighbours, key = lambda x: x.nr)
                    if neighbours[1].type == 'H':
                        if atom.chiral == 'counterclockwise':
                            atom_repr = f'[{atom.type}@H]'
                        elif atom.chiral == 'clockwise':
                            atom_repr = f'[{atom.type}@@H]'
                        else:
                            raise Exception("Can't have a chirality that is not clockwise or counterclockwise.")
                    else:
                        raise Exception("Hydrogen in incorrect position.")
                elif atom.get_hydrogen_nr == 0:
                    if atom.chiral == 'counterclockwise':
                        atom_repr = f'[{atom.type}@]'
                    elif atom.chiral == 'clockwise':
                        atom_repr = f'[{atom.type}@@]'
                    else:
                        raise Exception("Can't have a chirality that is not clockwise or counterclockwise.")

        return atom_to_repr

    def structure_to_smiles(self):
        pass

    def find_next_atom_nr(self):
        """
        Return the next available integer to label an atom

        Output
        ------
        next_atom_nr: int
        """
        atom_nrs = []

        for atom in self.graph:
            atom_nrs.append(atom.nr)

        atom_nrs.sort()

        next_atom_nr = None

        for i, atom_nr in enumerate(atom_nrs):
            if i != atom_nr:
                next_atom_nr = i
                break
        
        if next_atom_nr == None:
            next_atom_nr = atom_nrs[-1] + 1

        return next_atom_nr

    def get_steric_atoms(self):
        """
        Return list of atoms that are chiral centres

        Output
        ------
        steric atoms: list of [atom, ->]
            with each atom (Atom object) a chiral centre in the structure
        """
        steric_atoms = []
        for atom in self.graph:
            if atom.chiral:
                steric_atoms.append(atom)

        return steric_atoms

    def find_next_bond_nr(self):
        """
        Return the next available integer to label a bond

        Output
        ------
        next_bond_nr: int
        """
        bond_nrs = list(self.bonds.keys())
        bond_nrs.sort()

        next_bond_nr = None

        for i, bond_nr in enumerate(bond_nrs):
            if i != bond_nr:
                next_bond_nr = i
                break
        
        if next_bond_nr == None:
            next_bond_nr = bond_nrs[-1] + 1

        return next_bond_nr

    def make_bond_lookup(self):
        """
        Update bond_lookup from current collection of bonds
        """
        self.bond_lookup = {}
        for bond in self.bonds:
            atom_1 = self.bonds[bond].atom_1
            atom_2 = self.bonds[bond].atom_2
            if not atom_1 in self.bond_lookup:
                self.bond_lookup[atom_1] = {}
            if not atom_2 in self.bond_lookup:
                self.bond_lookup[atom_2] = {}
            self.bond_lookup[atom_1][atom_2] = self.bonds[bond]
            self.bond_lookup[atom_2][atom_1] = self.bonds[bond]
            

    def make_atom_index(self):
        """
        Return dict of {atom nr: atom}

        Output
        ------
        atom_nr_to_atom: dict of {atom nr: atom}
            with atom nr (int) the unique number associated with an atom based
            on the order they occurred in in the input SMILES, and atom
            (Atom object) the atom associated with that number
        """
        atom_nr_to_atom = {}
        for atom in self.graph:
            atom_nr_to_atom[atom.nr] = atom
        return atom_nr_to_atom

    def add_atom(self, atom_type, neighbours, chiral = None, charge = 0, aromatic = False):
        """
        Add a new atom to an already existing structure

        Input
        -----
        atom_type: str, an abbreviation of an element of the periodic table
        neighbours: list of Atom objects
        chiral: str or None, str indicating the chirality (clockwise or 
            counterclockwise) if the new atom is a chiral centre, None
            if the new atom is not a chiral centre
        charge: int, charge of the new atom
        """

        
        next_atom_nr = self.find_next_atom_nr()
        atom = Atom(atom_type, next_atom_nr, chiral, charge, aromatic)
        atom.add_shell_layout()

        
        for i, neighbour in enumerate(neighbours):
            next_bond_nr = self.find_next_bond_nr()
            self.make_bond(atom, neighbour, next_bond_nr)

        return atom

    def find_cycles(self):
        """
        Find all cycles in a structure and store them
        """
        self.cycles = find_cycles.Cycles(self)

    def promote_electrons_in_five_rings(self):
        """
        Promote electrons in aromatic five rings to p-orbitals
        """
        five_rings = self.cycles.find_five_membered()
        for five_ring in five_rings:
            aromatic, heteroatom = check_five_ring(five_ring)
            if aromatic:
                heteroatom.promote_lone_pair_to_p_orbital()
                if heteroatom.type == 'N':
                    heteroatom.pyrrole = True

    def find_aromatic_cycles(self):
        """
        Returns cycles that are aromatic

        Output
        ------
        aromatic_cycles: list of [[Atom, ->], ->], with each list of
            atoms the atoms that comprise an aromatic cycle

        """
        aromatic_cycles = []
        for cycle in self.cycles.unique_cycles:
            if check_aromatic(cycle):
                aromatic_cycles.append(cycle)

        return aromatic_cycles
                
    def find_aromatic_systems(self):
        """
        Return ring systems that are aromatic

        Output
        ------
        aromatic_ring_systems: list of [[Atom, ->], ->], with each list of
            atoms the atoms that comprise an aromatic ring system
        """
        
        ring_systems = self.cycles.find_cyclic_systems()
        aromatic_ring_systems = []
        for ring_system in ring_systems:
            if check_aromatic(ring_system):
                aromatic_ring_systems.append(ring_system)

        return aromatic_ring_systems

    def set_bonds_to_aromatic(self, aromatic_systems):
        for aromatic_system in aromatic_systems:
            for atom_1 in aromatic_system:
                atom_1.aromatic = True
                for atom_2 in aromatic_system:
                    if atom_1 in self.graph[atom_2]:
                        bond = self.bond_lookup[atom_1][atom_2]
                        if bond.type != 'aromatic':
                            bond.make_aromatic()


    def refine_structure(self):
        """
        """
        
        self.add_shells()
        self.add_hydrogens()
        self.add_shells()
        
        self.sort_by_nr()

        self.refine_p_bonds()
        self.hybridise_atoms()
        
        self.refine_s_bonds()
        self.drop_electrons()
        self.set_atom_neighbours()
        self.set_connectivities()
        self.make_lone_pairs()
        self.make_bond_lookup()
        self.find_cycles()
        self.promote_electrons_in_five_rings()
        aromatic_systems = self.find_aromatic_systems()
        aromatic_cycles = self.find_aromatic_cycles()
        self.set_bonds_to_aromatic(aromatic_systems)
        self.set_bonds_to_aromatic(aromatic_cycles)

    def make_lone_pairs(self):
        for atom in self.graph:
            atom.make_lone_pairs()

    def remove_bond_between_atoms(self, atom_1, atom_2):
        bond = self.bond_lookup[atom_1][atom_2]

        del self.bond_lookup[atom_1][atom_2]
        del self.bond_lookup[atom_2][atom_1]
        del self.bonds[bond.nr]


        if atom_2 in self.graph[atom_1]:
            self.graph[atom_1].remove(atom_2)
        if atom_1 in self.graph[atom_2]:
            self.graph[atom_2].remove(atom_1)


    def break_bond(self, bond):
        """
        Break a bond in the structure

        Input
        -----
        bond: Bond object
        """
        
        atom_1 = bond.atom_1
        atom_2 = bond.atom_2
        
        bond.break_bond()
        
        del self.bond_lookup[atom_1][atom_2]
        del self.bond_lookup[atom_2][atom_1]
        del self.bonds[bond.nr]

        if atom_2 in self.graph[atom_1]:
            self.graph[atom_1].remove(atom_2)
        if atom_1 in self.graph[atom_2]:
            self.graph[atom_2].remove(atom_1)
        

    def break_bond_by_nr(self, bond):
        """
        Break a bond in the structure based on bond number

        Input
        -----
        bond: Bond object or int, with int the number of a bond in the graph
        """
        if type(bond) == int:
            bond_nr = bond
            bond = self.bonds[bond]

        self.break_bond(bond)



    def break_bond_between_atoms(self, atom_1, atom_2):
        """
        Break a bond in the structure based on the numbers of neighbouring
        atoms

        Input
        -----
        atom_1: Atom object or int, with int the number of an atom in the graph
        atom_2: Atom object or int, with int the number of an atom in the graph

        """
        atom_nr_to_atom = self.make_atom_index()
        
        if atom_1.type == int:
            atom_1 = atom_nr_to_atom[atom_1]
        if atom_2.type == int:
            atom_2 = atom_nr_to_atom[atom_2]

        bond = self.bond_lookup[atom_1][atom_2]
        self.break_bond(bond)

    def remove_atom(self, atom_to_remove):
        """
        Remove an atom from the structure

        Input
        -----
        atom_to_remove: Atom object, atom that is to be removed from the
            structure
        """

        for bond in atom_to_remove.bonds:
            self.break_bond_by_nr(bond.nr)
            
        for atom in self.graph:
            if atom_to_remove in self.graph[atom]:
                self.graph[atom].remove(atom_to_remove)

        del self.graph[atom_to_remove]

    def set_connectivities(self):
        for atom in self.graph:
            if atom.type != 'H':
                atom.set_connectivity()

    def set_atom_neighbours(self):
        for atom in self.graph:
            atom.set_neighbours(self)

    def get_connectivities(self):
        connectivities = {}
        for atom in self.graph:
            if atom.type != 'H':
                connectivity = atom.connectivity
                if not connectivity in connectivities:
                    connectivities[connectivity] = []
                
                    connectivities[connectivity].append(atom)

        
        return connectivities

    def get_chiral_double_bonds(self):
        chiral_bonds = []
        for bond_nr in self.bonds:
            bond = self.bonds[bond_nr]
            if bond.chiral:
                chiral_bonds.append(bond)

        return chiral_bonds
                

    def check_chiral_double_bonds(self, child, match):
        chirality_matches = True
        chiral_bonds = child.get_chiral_double_bonds()

        for chiral_bond in chiral_bonds:
            neighbour_1, neighbour_2 = chiral_bond.neighbours

            parent_neighbour_1 = match[neighbour_1]
            parent_neighbour_2 = match[neighbour_2]
            parent_bond = self.bond_lookup[parent_neighbour_1][parent_neighbour_2]

            if not parent_bond.chiral:
                chirality_matches = False
                break
            else:
                matching_chirality = chiral_bond.check_same_chirality(parent_bond, match)
                if not matching_chirality:
                    chirality_matches = False
                    break

        return chirality_matches

    def check_chiral_centres(self, child, match):
        chirality_matches = True
        chiral_centres = child.get_steric_atoms()
        
        for chiral_centre in chiral_centres:
            parent_atom = match[chiral_centre]
            if parent_atom.chiral:
                child_neighbours = chiral_centre.neighbours
                parent_neighbours = parent_atom.neighbours
                
                chirality_matches = check_same_chirality(chiral_centre, parent_atom, match)
                if not chirality_matches:
                    break
            else:
                chirality_matches = False
                break
            
        return chirality_matches

    def find_substructures(self, substructure, check_chiral_centres = True, check_chiral_double_bonds = True):
        matches = []
        if self.is_substructure_atom_composition(substructure):
            if self.is_substructure_atom_connectivity(substructure):
                matches = self.find_substructure_in_structure(substructure)

        if check_chiral_centres:
            final_matches = []

            for match in matches:
                if self.check_chiral_centres(substructure, match):
                    final_matches.append(match)
            matches = final_matches

        if check_chiral_double_bonds:
            final_matches = []
            for match in matches:
                if self.check_chiral_double_bonds(substructure, match):
                    final_matches.append(match)
            matches = final_matches

        return matches
                                                      

    def is_substructure_atom_composition(self, child):
        
        atom_counts_self = self.get_atom_counts()
        atom_counts_child = child.get_atom_counts()

        can_be_substructure = True

        for atom_type in atom_counts_child:
            try:
                atom_nr_self = atom_counts_self[atom_type]
                atom_nr_child = atom_counts_child[atom_type]
                if atom_nr_child > atom_nr_self:
                    can_be_substructure = False
                    break
            except KeyError:
                can_be_substructure = False
                break

        return can_be_substructure

    def get_connectivity_counts(self):
        connectivities = {}
        for atom in self.graph:
            if atom.type != 'H':
                if not atom.type in connectivities:
                    connectivities[atom.type] = {}
                connectivity = atom.connectivity
                if not connectivity in connectivities[atom.type]:
                    connectivities[atom.type][connectivity] = 0
                connectivities[atom.type][connectivity] += 1

        return connectivities

    def get_substructure_connectivity_counts(self, atom_connectivities_child):
        substructure_connectivity_counts = {}
        for atom_type in atom_connectivities_child:
            substructure_connectivity_counts[atom_type] = {}
            for connectivity in atom_connectivities_child[atom_type]:
                substructure_connectivity_counts[atom_type][connectivity] = 0
                for atom in self.graph:
                    if atom.type == atom_type and atom.potential_same_connectivity(connectivity):
                            substructure_connectivity_counts[atom_type][connectivity] += 1

        return substructure_connectivity_counts

    def get_substructure_connectivities(self, atom_connectivities_child):
        substructure_connectivities = {}
        for substructure_connectivity in atom_connectivities_child:
            substructure_connectivities[substructure_connectivity] = []
            for atom in self.graph:
                if atom.type != 'H' and atom.potential_same_connectivity(substructure_connectivity):
                        substructure_connectivities[substructure_connectivity].append(atom)

        
        return substructure_connectivities

    def is_substructure_atom_connectivity(self, child):
        
        atom_connectivities_child = child.get_connectivity_counts()
        atom_connectivities_self = self.get_substructure_connectivity_counts(atom_connectivities_child)

        can_be_substructure = True

        for atom_type in atom_connectivities_child:
            for connectivity in atom_connectivities_child[atom_type]:
                connectivity_nr_self = atom_connectivities_self[atom_type][connectivity]
                connectivity_nr_child = atom_connectivities_child[atom_type][connectivity]
                if connectivity_nr_child > connectivity_nr_self:
                    can_be_substructure = False
                    break

        return can_be_substructure

    def make_bond_dict(self):
        bond_dict = {}

        for atom in self.graph:
            if atom.type != 'H':
                bond_dict[atom] = 0
                for neighbour in atom.neighbours:
                    if neighbour.type != 'H':
                        bond_dict[atom] += 1
                

        return bond_dict

    def make_match_dict(self):
        match_dict = {}
        for atom in self.graph:
            if atom.type != 'H':
                match_dict[atom] = None

        return match_dict

    def compress_graph(self):
        compressed_graph = {}
        for atom in self.graph:
            if atom.type != 'H':
                new_atom = copy.deepcopy(atom)
                compressed_graph[new_atom] = copy.deepcopy(self.graph[atom])

    def reset_visited(self):
        for atom in self.graph:
            atom.visited = False
        for bond_nr in self.bonds:
            bond = self.bonds[bond_nr]
            bond.visited = False

    def remove_visited(self):
        for atom in self.graph:
            delattr(atom, 'visited')
        for bond_nr in self.bonds:
            bond = self.bonds[bond_nr]
            delattr(bond, 'visited')

    def traceback(self, placed_atoms_parent, match_dict, reverse_match_dict, parent_options, atoms_to_place, child_bond_dict):

        new_parent_candidate = None
        new_child_candidate = None
        previous_child_candidate = None


        for i in range(len(placed_atoms_parent) - 1, -1, -1):
            current_atom = placed_atoms_parent[i]

            try:
                for option in parent_options[current_atom]:
                    bond = self.bond_lookup[current_atom][option]
                        
                    if not bond.visited:
                        new_parent_candidate = current_atom
                        new_child_candidate = reverse_match_dict[current_atom]
                        break
            except KeyError:
                pass

            if new_parent_candidate:
                break
            else:
                del placed_atoms_parent[i]
                child_current = reverse_match_dict[current_atom]
                atoms_to_place.add(child_current)
                    
                match_dict[child_current] = None

                child_bond_dict[child_current] += 1
                
                for bond in child_current.bonds:
                    for atom in bond.neighbours:
                        if atom != child_current:
                            if atom in match_dict:
                                child_bond_dict[atom] += 1
                    
                    
                del reverse_match_dict[current_atom]

        return new_child_candidate, new_parent_candidate

    def find_substructure_in_structure(self, child):
       
        atom_connectivities_child = child.get_connectivities()
        atom_connectivities_parent = self.get_substructure_connectivities(atom_connectivities_child)
        #Sort based on the complexity of the connectivity
        
        connectivities = sorted(list(atom_connectivities_child.keys()),
                                key = lambda x: len(set(x)), reverse = True)

        starting_connectivity = connectivities[0]
        starting_atom = atom_connectivities_child[starting_connectivity][0]
        
        seeds = atom_connectivities_parent[starting_connectivity]
        matches = []

        for seed in seeds:

            self.reset_visited()

            child_bond_dict = child.make_bond_dict()
            
            
            match_active = True
            match = child.make_match_dict()
            reverse_match = {}
            match[starting_atom] = seed
            reverse_match[seed] = starting_atom
            
            atoms_to_place = set(copy.copy(list(child.graph.keys())))
            for atom in copy.deepcopy(atoms_to_place):
                if atom.type == 'H':
                    atoms_to_place.remove(atom)
                    
            atoms_to_place.remove(starting_atom)
            
            placed_atoms_child = set([starting_atom])
            placed_atoms_parent = [seed]
            parent_options = {}
            
            current_atom_child = starting_atom
            current_atom_parent = seed

            while atoms_to_place and match_active:

                child_neighbours = current_atom_child.neighbours
                next_atom_child = None
                
                for neighbour in child_neighbours:
                    if neighbour in atoms_to_place and neighbour.type != 'H':
                        next_atom_child = neighbour
                        
                        break

                if next_atom_child:
                    child_bond = child.bond_lookup[current_atom_child][next_atom_child]
                    next_atom_parent_options = []
                    parent_neighbours = current_atom_parent.neighbours
                    for neighbour in parent_neighbours:
                        parent_bond = self.bond_lookup[current_atom_parent][neighbour]
                        
                        
                        if neighbour.type == next_atom_child.type and child_bond.type == parent_bond.type:
                            if not self.bond_lookup[current_atom_parent][neighbour].visited:
                                if neighbour.potential_same_connectivity(next_atom_child.connectivity):
                                    next_atom_parent_options.append(neighbour)

                    parent_options[current_atom_parent] = []

                    if next_atom_parent_options:
                        
                        next_atom_parent_options = sorted(next_atom_parent_options, key = lambda x: len(set(x.connectivity)), reverse = True)
                        for option in next_atom_parent_options:
                            bond = self.bond_lookup[current_atom_parent][option]
                            parent_options[current_atom_parent].append(bond)


                        next_atom_parent = next_atom_parent_options[0]
                        self.bond_lookup[current_atom_parent][next_atom_parent].visited = True
                        
                        match[next_atom_child] = next_atom_parent
                        reverse_match[next_atom_parent] = next_atom_child
                        placed_atoms_parent.append(next_atom_parent)
                        
                        child_bond_dict[current_atom_child] -= 1
                        child_bond_dict[next_atom_child] -= 1
                        
                        current_atom_child = next_atom_child
                        current_atom_parent = next_atom_parent
                        atoms_to_place.remove(current_atom_child)


                    else:
                        new_child_candidate, new_parent_candidate = self.traceback(placed_atoms_parent, match, reverse_match, parent_options, atoms_to_place, child_bond_dict)
                        if not new_child_candidate:
                            match_active = False
                        else:
                            current_atom_child = new_child_candidate
                            current_atom_parent = new_parent_candidate


                    

                else:
                    if atoms_to_place:

                        for atom in match.keys():
                            if atom != current_atom_child:
                                if child_bond_dict[atom] > 0:
                                    current_atom_child = atom
                                    current_atom_parent = match[atom]
                                    break
                    else:
                        pass

            if match_active:
                matches.append(match)

        compare_all_matches(matches)

        self.remove_visited()

        return matches

    def drop_electrons(self):
        for atom in self.graph:
            atom.drop_electrons()

    def add_shells(self):
        for atom in self.graph:
            if not atom.shells:
                atom.add_shell_layout()

    def hybridise_atoms(self, atoms = None):

        if not atoms:
            for atom in self.graph:
                atom.hybridise()

        else:
            for atom in atoms:
                atom.hybridise()

    def add_hydrogens(self):
        max_atom_nr = max([atom.nr for atom in self.graph])
        max_bond_nr = max(list(self.bonds.keys()))
        
        for atom in list(self.graph.keys()):
            hydrogens_to_add = atom.calc_hydrogens()
            for i in range(hydrogens_to_add):
                max_atom_nr += 1
                max_bond_nr += 1
                hydrogen = Atom('H', max_atom_nr, None, 0, False)
                self.add_bond(atom, hydrogen, 'single', max_bond_nr)
                

    def refine_p_bonds(self):
        for bond_nr in self.bonds:
            bond = self.bonds[bond_nr]
            if bond.type != 'single':

                bond.combine_p_orbitals()


    def refine_s_bonds(self):
        for bond_nr, bond in self.bonds.items():
            bond.combine_hybrid_orbitals()
        
    def get_atom_counts(self):
        atom_counts = {}
        for atom in self.graph:
            if atom.type != 'H':
                if atom.type not in atom_counts:
                    atom_counts[atom.type] = 0
                atom_counts[atom.type] += 1

        return atom_counts
            

    def add_disconnected_atom(self, atom):
        self.graph[atom] = []

    def sort_by_nr(self):
        for atom in self.graph:
            self.graph[atom].sort(key = lambda x: x.nr)

    def make_dummy_bond(self, atom_1, atom_2, bond_nr, dummy = False):
        if dummy:
            bond_type = 'dummy'
        else:
            bond_type = 'single'
            
        if atom_1 in self.graph:
            self.graph[atom_1].append(atom_2)
        else:
            self.graph[atom_1] = [atom_2]

        if atom_2 in self.graph:
            self.graph[atom_2].append(atom_1)
        else:
            self.graph[atom_2] = [atom_1]

        bond = Bond(atom_1, atom_2, bond_type, bond_nr)

        
        atom_1.add_bond(bond)
        atom_2.add_bond(bond)
        
        self.bonds[bond_nr] = bond

        if not atom_1 in self.bond_lookup:
            self.bond_lookup[atom_1] = {}
        if not atom_2 in self.bond_lookup:
            self.bond_lookup[atom_2] = {}

        self.bond_lookup[atom_1][atom_2] = bond
        self.bond_lookup[atom_2][atom_1] = bond

    def make_bond(self, atom_1, atom_2, bond_nr):
        if atom_1 in self.graph:
            self.graph[atom_1].append(atom_2)
        else:
            self.graph[atom_1] = [atom_2]

        if atom_2 in self.graph:
            self.graph[atom_2].append(atom_1)
        else:
            self.graph[atom_2] = [atom_1]

        bond = Bond(atom_1, atom_2, 'single', bond_nr)
        
        atom_1.add_bond(bond)
        atom_2.add_bond(bond)
        
        self.bonds[bond_nr] = bond

        electron_1 = None
        electron_2 = None

        orbital_1 = None
        orbital_2 = None

        for orbital in atom_1.valence_shell.orbitals:
            if atom_1.valence_shell.orbitals[orbital].electron_nr == 1 and not atom_1.valence_shell.orbitals[orbital].electrons[0].aromatic:
                orbital_1 = atom_1.valence_shell.orbitals[orbital]
                electron_1 = orbital_1.electrons[0]
                break
                

        for orbital in atom_2.valence_shell.orbitals:
            if atom_2.valence_shell.orbitals[orbital].electron_nr == 1 and not atom_2.valence_shell.orbitals[orbital].electrons[0].aromatic:
                orbital_2 = atom_2.valence_shell.orbitals[orbital]
                electron_2 = orbital_2.electrons[0]
                break

        orbital_1.add_electron(electron_2)
        orbital_2.add_electron(electron_1)

        orbital_1.set_bond(bond, 'sigma')
        orbital_2.set_bond(bond, 'sigma')

        if not atom_1 in self.bond_lookup:
            self.bond_lookup[atom_1] = {}
        if not atom_2 in self.bond_lookup:
            self.bond_lookup[atom_2] = {}

        self.bond_lookup[atom_1][atom_2] = bond
        self.bond_lookup[atom_2][atom_1] = bond

    def add_bond(self, atom_1, atom_2, bond_type, bond_nr, chiral_symbol=None):
        if atom_1 in self.graph:
            self.graph[atom_1].append(atom_2)
        else:
            self.graph[atom_1] = [atom_2]

        if atom_2 in self.graph:
            self.graph[atom_2].append(atom_1)
        else:
            self.graph[atom_2] = [atom_1]

        bond = Bond(atom_1, atom_2, bond_type, bond_nr)

        bond.chiral_symbol = chiral_symbol
        
        atom_1.add_bond(bond)
        atom_2.add_bond(bond)
        
        self.bonds[bond_nr] = bond

        if not atom_1 in self.bond_lookup:
            self.bond_lookup[atom_1] = {}
        if not atom_2 in self.bond_lookup:
            self.bond_lookup[atom_2] = {}

        self.bond_lookup[atom_1][atom_2] = bond
        self.bond_lookup[atom_2][atom_1] = bond


    def print_graph(self):
        pprint(self.graph)

    def print_bonds(self):
        pprint(self.bonds)

    def find_aromatic_bonds(self):
        aromatic_bonds = []
        for bond in self.bonds:
            if self.bonds[bond].type == 'aromatic':
                aromatic_bonds.append(self.bonds[bond])

        return aromatic_bonds

    def find_aromatic_atoms(self):
        aromatic_bonds = self.find_aromatic_bonds()
        aromatic_atoms = []
        for atom in self.graph:
            for bond in atom.bonds:
                if bond in aromatic_bonds:
                    aromatic_atoms.append(atom)

        return aromatic_atoms

    def find_aromatic_graphs(self):
        aromatic_atoms = self.find_aromatic_atoms()
        aromatic_bonds = self.find_aromatic_bonds()
        aromatic_graph = {}
        
        for atom in copy.deepcopy(self.graph):
            if atom in aromatic_atoms:
                aromatic_graph[copy.deepcopy(atom)] = copy.deepcopy(self.graph[atom])

        for atom in aromatic_graph:
            for bond in atom.bonds:
                if bond.type != 'aromatic':
                    atom.bonds.remove(bond)

        for atom_1 in aromatic_graph:
            
            for atom_2 in aromatic_graph[atom_1]:
                aromatic = False
                for bond in atom_1.bonds:
                    if bond.atom_1 == atom_2 or bond.atom_2 == atom_2:
                        aromatic = True

                if not aromatic:
                    aromatic_graph[atom_1].remove(atom_2)

        aromatic_structure = Structure(aromatic_graph, aromatic_bonds)
        
                
        aromatic_structures = aromatic_structure.split_disconnected_structures()

        return aromatic_structures
        
    def find_pi_subgraph(self, prune = True):
        pi_subgraph = {}

        for bond_nr, bond in self.bonds.items():
            if bond.type == 'aromatic':


                # prune the subgraph as kekulisation can only occur in atoms
                # #that have an unpaired electron

                unpaired_electrons_1 = 0
                unpaired_electrons_2 = 0

                for orbital_name, orbital in bond.atom_1.valence_shell.orbitals.items():
                    if len(orbital.electrons) == 1:
                        unpaired_electrons_1 += 1

                for orbital_name, orbital in bond.atom_2.valence_shell.orbitals.items():
                    if len(orbital.electrons) == 1:
                        unpaired_electrons_2 += 1

                if unpaired_electrons_1 and unpaired_electrons_2:

                    if not bond.atom_1 in pi_subgraph:
                        pi_subgraph[bond.atom_1] = []
                    if not bond.atom_2 in pi_subgraph:
                        pi_subgraph[bond.atom_2] = []

                    pi_subgraph[bond.atom_1].append(bond.atom_2)
                    pi_subgraph[bond.atom_2].append(bond.atom_1)

                elif not prune:

                    if not bond.atom_1 in pi_subgraph:
                        pi_subgraph[bond.atom_1] = []
                    if not bond.atom_2 in pi_subgraph:
                        pi_subgraph[bond.atom_2] = []

                    pi_subgraph[bond.atom_1].append(bond.atom_2)
                    pi_subgraph[bond.atom_2].append(bond.atom_1)


        return pi_subgraph


    def kekulise(self):

        pruned = self.find_pi_subgraph(prune = True)
        unpruned = self.find_pi_subgraph(prune = False)

        aromatic_unmatched = set(unpruned.keys()) - set(pruned.keys())

        matching = Match.from_structure(self)
        unmatched_nodes = matching.unmatched_nodes()
        if unmatched_nodes != 0:
            raise Exception("This structure cannot be kekulised!")
        else:
            double_bond_pairs = set()
            single_bond_pairs = set()

            for node in matching.nodes:
                double_bond_pair = tuple(sorted([node.atom, node.mate.atom], key = lambda x: x.nr))
                if not double_bond_pair in double_bond_pairs:
                    double_bond_pairs.add(double_bond_pair)

                for neighbour in node.neighbors:
                    if neighbour.index != node.mate.index:
                        single_bond_pair = tuple(sorted([node.atom, neighbour.atom], key = lambda x: x.nr))
                        if not single_bond_pair in single_bond_pairs:
                            single_bond_pairs.add(single_bond_pair)

            for atom in aromatic_unmatched:
                for neighbour in atom.neighbours:
                    if neighbour in pruned:
                        single_bond_pair = tuple(sorted([atom, neighbour], key = lambda x: x.nr))
                        if not single_bond_pair in single_bond_pairs:
                            single_bond_pairs.add(single_bond_pair)



        kekule_structure = copy.deepcopy(self)

        for pair in double_bond_pairs:
            bond = kekule_structure.bond_lookup[pair[0]][pair[1]]
            bond.type = 'double'
            bond.aromatic = False
            bond.atom_1.aromatic = False
            bond.atom_2.aromatic = False

        for pair in single_bond_pairs:
            bond = kekule_structure.bond_lookup[pair[0]][pair[1]]
            bond.type = 'single'
            bond.aromatic = False
            bond.atom_1.aromatic = False
            bond.atom_2.aromatic = False
            bond.atom_1.pyrrole = False
            bond.atom_2.pyrrole = False

        return kekule_structure

        

    def break_any_bond(self, atom_1, atom_2):
        """Remove edges from structure to break any bond

        atom_1: tuple of (str, int), with str atom type and int atom number
        atom_2: tuple of (str, int), with str atom type and int atom number

        """
        bonds_removed = 0
        while atom_2 in self.structure[atom_1]:
            self.structure[atom_1].remove(atom_2)
            bonds_removed += 1
        while atom_1 in self.structure[atom_2]:
            self.structure[atom_2].remove(atom_1)
            
        self.bonds[atom_1] -= bonds_removed
        self.bonds[atom_2] -= bonds_removed

    def remove_an_atom(self, atom):
        """Remove a node from structure to remove an atom

        atom_1: tuple of (str, int), with str atom type and int atom number
        """
        del self.strucrure[atom]

    def remove_isolated_atoms(self):
        """Remove nodes without outgoing edges from structure

        atom_1: tuple of (str, int), with str atom type and int atom number
        """
        for atom in list(self.graph.keys()):
            if not self.graph[atom]:
                del self.graph[atom]

    def split_disconnected_structures(self):
        """Return list of unconnected structures from structure


        Output:
        new_graphs: list of dicts of {node: [node, ->], ->}, with each dict a
            bidirectional graph that is not connected to the other graphs in
            the list
        """
        working_graph = copy.deepcopy(self)
        new_graphs = []
        working_graph.make_bond_nr_dict()
        

        working_graph.remove_connectors()

        start_node = list(working_graph.graph.keys())[0]

        paths_collection = []
        paths = []

        while start_node:
            path = working_graph.find_a_path(start_node)
            paths.append(path)

            potential_start_nodes = working_graph.find_start_nodes(paths)

            try:
                start_node = potential_start_nodes[0]
                
            except IndexError:
                paths_collection.append(paths)
                paths = []
                potential_start_nodes = working_graph.find_new_start_node()
                
                
                try:
                    start_node = potential_start_nodes[0]
                    
                except IndexError:
                    paths_collection.append(paths)
                    start_node = None
                
        for paths in paths_collection:
            if paths:
                new_graph = working_graph.put_paths_in_graph(paths)
                new_graphs.append(new_graph)

        #add back connectors

        for new_graph in new_graphs:
            for node in new_graph:
                new_graph[node] = self.graph[node]

        #Add lone atoms
        if working_graph.graph:
            for atom in working_graph.graph:
                new_graph = {atom: []}
                if new_graph not in new_graphs:
                    new_graphs.append(new_graph)

        new_structures = []

        for new_graph in new_graphs:
            new_structures.append(Structure(new_graph))

        for new_structure in new_structures:
            new_structure.infer_bonds()
            new_structure.make_bond_lookup()

        return new_structures

    def infer_bonds(self):
        self.bonds = {}
        for atom in self.graph:
            for bond in atom.bonds:
                if not bond.nr in self.bonds:
                    self.bonds[bond.nr] = bond

    def find_lowest_atom_nr(self):
        lowest_atom_nr = 4242424242424242

        for atom in self.structure:
            if atom.nr < lowest_atom_nr:
                lowest_atom_nr = atom.nr

        return lowest_atom_nr

    def sort_substructures_by_nr(self):
        disconnected_structures = self.split_disconnected_structures()

        index_dict = {}

        for i, structure in enumerate(disconnected_structures):
            lowest_atom_nr = structure.find_lowest_atom_nr()
            index_dict[lowest_atom_nr] = i
            
        new_disconnected_structures = []

        sorted_indices = list(index_dict.keys())
        sorted_indices.sort()
        
        for index in sorted_indices:
            new_disconnected_structures += [disconnected_structures[index_dict[index]]]

        return new_disconnected_structures

    def make_sorted_copy(self):
        copy = Structure(self.structure)
        for atom in copy.structure:
            copy.structure[atom].sort()

        return copy

    def compare_structures(self, structure_2):
        sorted_copy_group_1 = self.make_sorted_copy()
        sorted_copy_group_2 = structure_2.make_sorted_copy()
        
        if sorted_copy_group_1.structure == sorted_copy_group_2.structure:
            return True
        else:
            return False

    def label_smiles_nrs(self, K_to_S_dict):
        for atom in self.structure:
            atom.set_smiles_nr(K_to_S_dict)

    def get_true_neighbour(self, atom_1, atom_2):
        bond = self.bond_lookup[atom_1][atom_2]
        true_neighbour = None
        for atom in self.graph[atom_1]:
            if atom in self.graph[atom_2] and self.bond_lookup[atom_1][atom].type != 'dummy' and  self.bond_lookup[atom_2][atom].type != 'dummy':
                true_neighbour = atom
                break

        return true_neighbour

    def structure_to_smiles(self):
        valence_dict = self.make_valence_dict()
        strings = []
        for atom in self.structure:
            atom_string = []
            H_string = self.make_H_string(atom, valence_dict)
            atom_string += [H_string]
            
            bonds = self.bond_record[atom]
            for bond in bonds:
                nr_string = self.make_nr_string(bond)
                atom_string += [nr_string]

            atom_string += ['.']
            

            strings += [''.join(atom_string)]

        string = ''.join(strings)[:-1]
        mol = Chem.MolFromSmiles(string)
        for i in range(mol.GetNumAtoms()):
            mol.GetAtomWithIdx(i).SetAtomMapNum(i)
            

        return Smiles(Chem.MolToSmiles(mol))
                    
                
        
            
        
        
                

    #========================================================================
    #Auxillary functions
    #========================================================================

    def make_H_string(self, atom, valence_dict):
        H_nr = valence_dict[atom]['max'] - valence_dict[atom]['current']
        if H_nr > 1:
            H_string = "[%sH%d]" % (atom.type, H_nr)
        elif H_nr == 1:
            H_string = "[%sH]" % (atom.type)
        else:
            H_string = "[%s]" % (atom.type)
            
        return H_string

    def make_nr_string(self, bond):
        bond_nr, bond_type = bond
        if bond_type == 1:
            string = ''
        elif bond_type == 2:
            string = '='
        elif bond_type == 3:
            string = '#'

        if bond_nr > 9:
            bond_string = '%%%d' % bond_nr
        else:
            bond_string = str(bond_nr)

        return ''.join([string, bond_string])
        
            
    def make_valence_dict(self):
        valence_dict = {}
        for atom in self.structure:
            valence_dict[atom] = {}
            valence_dict[atom]['current'] = len(self.structure[atom])
            valence_dict[atom]['max'] = atom.get_valence(len(self.structure[atom]))

        return valence_dict

    def record_bonds(self):
        bond_record = {}

        counter = 1

        self.bond_record = {}
        for atom_1 in self.structure:
            self.bond_record[atom_1] = []
            
            for atom_2 in self.structure[atom_1]:
                bond_type = self.structure[atom_1].count(atom_2)
                
                if atom_1.nr > atom_2.nr:
                    pair = (atom_2, atom_1)
                else:
                    pair = (atom_1, atom_2)

                if not pair in bond_record:
                    bond_nr = counter
                    bond_record[pair] = bond_nr
                    counter += 1
                else:
                    bond_nr = bond_record[pair]
                    
                if not (bond_nr, bond_type) in self.bond_record[atom_1]:
                    self.bond_record[atom_1] += [(bond_nr, bond_type)]
            self.bond_record[atom_1].sort()
                

    def put_paths_in_graph(self, paths):
        """Return single structure from bond paths

        Input:
        paths: list of [atom, ->], with each atom a tuple of (str, int), with
            str atom type and int atom number

        Output:
        rest_group_graph: dict of {atom: [atom, ->], ->}, with each atom a tuple
            of (str, int), str representing atom type, and int representing atom
            number. Example: ('O', 1). This graph represents the side chain of
            an amino acid
        """
        rest_group_graph = {}
        for path in paths:
            current_atom = path[0]
            if path[1:]:
                for atom in path[1:]:

                    next_atom = atom
                    if current_atom in rest_group_graph:
                        rest_group_graph[current_atom] += [next_atom]
                    else:
                        rest_group_graph[current_atom] = [next_atom]

                    if next_atom in rest_group_graph:
                        rest_group_graph[next_atom] += [current_atom]
                    else:
                        rest_group_graph[next_atom] = [current_atom]

                    current_atom = atom
            else:
                rest_group_graph = {current_atom: []}

        return rest_group_graph

    def get_bond_nr(self, atom_1, atom_2):
        bond_nr = self.structure[atom_1].count(atom_2)
        return bond_nr

    def make_bond_nr_dict(self):
        """
        """
        self.bond_nr_dict = {}
        for atom in self.graph:
            self.bond_nr_dict[atom] = len(atom.bonds)


    def find_a_path(self, start_atom):
        """Return a list of linked atoms from a structure

        Input:
        structure: dict of {atom: [atom, ->], ->}, with each atom a tuple of
            (str, int), str representing atom type, and int representing atom
            number. Example: ('O', 1)
        start_atom: tuple of (str, int), with str atom type and int atom number
        bond_dict: dict of {atom: remaining_bonds, ->}, with atom tuple of
            (str, int), with str atom type and int atom number, and remaining bonds
            int

        Output:
        path: list of [atom, ->], where adjacent atoms are connected, and each atom
            is a tuple of (str, int), with str atom type and int atom number
        """

        nodes = list(self.graph.keys())
        
        current_atom = start_atom
        path = [current_atom]

        if len(self.graph[current_atom]) == 0:
            path = [current_atom]
            return path


        # keep trying to extend the path until there are no bonds to traverse
        while True:
            try:

                next_atom = self.graph[current_atom][0]

                path.append(next_atom)

                # remove traversed bond from structure
                self.graph[current_atom].remove(next_atom)
                if not next_atom == current_atom:
                    self.graph[next_atom].remove(current_atom)

                self.bond_nr_dict[current_atom] -= 1
                self.bond_nr_dict[next_atom] -= 1

                # remove atom from structure if no more untraversed bonds come
                # from it
                if not self.graph[current_atom]:
                    del self.graph[current_atom]

                current_atom = next_atom

                if not self.graph[current_atom]:
                    del self.graph[current_atom]
                
            except KeyError:
                break
            
        return path

    

    def find_start_nodes(self, paths):
        """Return atoms that still have outgoing bonds within an existing path

        Input:
        paths: list of [atom, ->], with each atom a tuple of (str, int), with
            str atom type and int atom number
        bond_dict: dict of {atom: remaining_bonds, ->}, with atom tuple of
            (str, int), with str atom type and int atom number, and remaining
            bonds int
            
        Output:
        start_atoms: list of [atom, ->], with each atom a tuple of (str, int),
            with str atom type and int atom number


        """
        
        start_atoms = []
        for path in paths:
            for atom in path:
                if self.bond_nr_dict[atom] > 0:
                    start_atoms.append(atom)

        return start_atoms

    def remove_connectors(self):
        """Remove nodes that only have incoming edges from graph

        Input:
        working_graph: dict of {node: [node, ->], ->}, representing a graph

        """

        for node in self.graph:
            for next_node in self.graph[node]:
                if next_node not in self.graph:
                    self.graph[node].remove(next_node)

    def find_new_start_node(self):
        """Return list of nodes that still have outgoing edges

        Input:
        edges_dict: dict of {node: int, ->}, with int representing the number of
            outgoing edges of that node

        Output:
        start_nodes: list of [node, ->], with each node an immutable object
        """
        start_nodes = []
        for atom in self.graph:
            if self.bond_nr_dict[atom] > 0:
                start_nodes.append(atom)

        return start_nodes





class WorkingStructure(Structure):
    def __init__(self, graph = {}, bonds = {}):
        Structure.__init__(self, graph)

class RingSystem():
    def __init__(self, atoms, bonds):
        self.atoms = atoms
        self.bonds = bonds

    def get_pi_electrons(self):
        self.pi_electrons = 0
        for bond in self.bonds:
            if bond.type == 'aromatic':
                self.pi_electrons += 1
            if bond.type == 'double':
                self.pi_electrons += 2


class Bond:
    bond_types = {'single', 'double', 'triple', 'quadruple', 'aromatic', 'ionic', 'dummy'}

    
                 
    def __init__(self, atom_1, atom_2, bond_type, bond_nr):
        atoms = [atom_1, atom_2]
        atoms.sort(key = lambda a: a.nr)

        self.atom_1 = atoms[0]
        self.atom_2 = atoms[1]
        self.neighbours = atoms

        self.chiral_symbol = None

        try:

            assert bond_type in self.bond_types

        except AssertionError:
            print(bond_type)
        
        self.type = bond_type
        self.nr = bond_nr
        self.aromatic = False
        if bond_type == 'aromatic':
            self.aromatic = True
        self.electrons = []

        self.chiral = False
        self.chiral_dict = {}

        if self.type == 'dummy':
            self.cbond = 0.3

        else:
            self.cbond = 0.1

        self.draw = BondDrawProperties()

    def __eq__(self, bond):
        return self.nr == bond.nr

    def __hash__(self):
        return self.nr

    def __repr__(self):
        return f'{self.type}_{self.nr}:{self.atom_1}_{self.atom_2}'

    def make_aromatic(self):
        self.type = 'aromatic'

        self.atom_1.aromatic = True
        self.atom_2.aromatic = True

        for orbital_name, orbital in self.atom_1.valence_shell.orbitals.items():
            if orbital.orbital_type == 'p':
                for electron in orbital.electrons:
                    if electron.atom != self.atom_1:
                        orbital.remove_electron(electron)
                    else:
                        electron.set_aromatic()

        for orbital_name, orbital in self.atom_2.valence_shell.orbitals.items():
            if orbital.orbital_type == 'p':
                for electron in orbital.electrons:
                    if electron.atom != self.atom_2:
                        orbital.remove_electron(electron)
                    else:
                        electron.set_aromatic()
            
        

        
    def check_same_chirality(self, parent_bond, match):
        same_chirality = True
        for atom in self.chiral_dict:
            if atom.type != 'H':
                parent_atom = match[atom]
                for atom_2 in self.chiral_dict[atom]:
                    if atom_2.type != 'H':
                        parent_atom_2 = match[atom_2]
                        orientation = self.chiral_dict[atom][atom_2]
                        parent_orientation = parent_bond.chiral_dict[parent_atom][parent_atom_2]
                        if orientation != parent_orientation:
                            same_chirality = False
                            break
            if not same_chirality:
                break

        return same_chirality

    def break_bond(self):
        assert self.type == 'single'

        electron_1, electron_2 = self.electrons
        orbital_1 = electron_1.orbital
        orbital_2 = electron_2.orbital

        orbital_1.remove_electron(electron_2)
        orbital_2.remove_electron(electron_1)
        
        orbital_1.remove_bond()
        orbital_2.remove_bond()
        
        self.atom_1.remove_neighbour(self.atom_2)
        self.atom_2.remove_neighbour(self.atom_1)
        
        self.atom_1.remove_bond(self)
        self.atom_2.remove_bond(self)

        
    def combine_hybrid_orbitals(self):

        s_bonding_orbital_1 = None
        s_bonding_orbital_2 = None

        for orbital_name in self.atom_1.valence_shell.orbitals:
            orbital = self.atom_1.valence_shell.orbitals[orbital_name]
            if 's' in orbital.orbital_type and orbital.electron_nr == 1:
                s_bonding_orbital_1 = orbital

        if not s_bonding_orbital_1:
            promotable = self.atom_1.is_promotable()
            if promotable:
                self.atom_1.promote_pi_bond_to_d_orbital()
                for orbital_name in self.atom_1.valence_shell.orbitals:
                    orbital = self.atom_1.valence_shell.orbitals[orbital_name]
                    if 's' in orbital.orbital_type and orbital.electron_nr == 1:
                        s_bonding_orbital_1 = orbital
                    

        for orbital_name in self.atom_2.valence_shell.orbitals:
            orbital = self.atom_2.valence_shell.orbitals[orbital_name]
            if 's' in orbital.orbital_type and orbital.electron_nr == 1:
                s_bonding_orbital_2 = orbital

        if not s_bonding_orbital_2:
            promotable = self.atom_2.is_promotable()
            if promotable:
                self.atom_2.promote_pi_bond_to_d_orbital()
                for orbital_name in self.atom_2.valence_shell.orbitals:
                    orbital = self.atom_2.valence_shell.orbitals[orbital_name]
                    if 's' in orbital.orbital_type and orbital.electron_nr == 1:
                        s_bonding_orbital_2 = orbital

        electron_1 = s_bonding_orbital_1.electrons[0]
        electron_2 = s_bonding_orbital_2.electrons[0]

        self.electrons.append(electron_1)
        self.electrons.append(electron_2)

        s_bonding_orbital_1.add_electron(electron_2)
        s_bonding_orbital_2.add_electron(electron_1)

        s_bonding_orbital_1.set_bond(self, 'sigma')
        s_bonding_orbital_2.set_bond(self, 'sigma')

    def combine_p_orbitals(self):
        assert self.type != 'single'

        

        if self.atom_1.pyrrole or self.atom_2.pyrrole:
            pass
        else:
            p_bonding_orbitals_1 = []
            electrons_found = 0
           
            for orbital_name, orbital in self.atom_1.valence_shell.orbitals.items():
                if orbital.orbital_type == 'p' and orbital.electron_nr == 1:
                    if (orbital.electrons[0].aromatic and self.type == 'aromatic') or not orbital.electrons[0].aromatic:
                        electrons_found += 1
                        p_bonding_orbitals_1.append(orbital)

                        if electrons_found == BOND_PROPERTIES.bond_type_to_p_orbitals[self.type]:
                            break


            p_bonding_orbitals_2 = []
            electrons_found = 0
                    
            for orbital_name, orbital in self.atom_2.valence_shell.orbitals.items():
                if orbital.orbital_type == 'p' and orbital.electron_nr == 1:
                    if (orbital.electrons[0].aromatic and self.type == 'aromatic') or not orbital.electrons[0].aromatic:
                        electrons_found += 1
                        p_bonding_orbitals_2.append(orbital)
                        
                        if electrons_found == BOND_PROPERTIES.bond_type_to_p_orbitals[self.type]:
                            break


            assert len(p_bonding_orbitals_1) == len(p_bonding_orbitals_2)

            if self.type == 'aromatic':
                assert len(p_bonding_orbitals_1) == len(p_bonding_orbitals_2) == 1
                electron_1 = p_bonding_orbitals_1[0].electrons[0]
                electron_2 = p_bonding_orbitals_2[0].electrons[0]
                
                electron_1.set_aromatic()
                electron_2.set_aromatic()

                self.electrons.append(electron_1)
                self.electrons.append(electron_2)
                
                p_bonding_orbitals_1[0].set_bond(self, 'pi')
                p_bonding_orbitals_2[0].set_bond(self, 'pi')
            

            else:

                for i in range(len(p_bonding_orbitals_1)):
                    electron_1 = p_bonding_orbitals_1[i].electrons[0]
                    electron_2 = p_bonding_orbitals_2[i].electrons[0]
                    

                    p_bonding_orbitals_1[i].add_electron(electron_2)
                    p_bonding_orbitals_2[i].add_electron(electron_1)

                    self.electrons.append(electron_1)
                    self.electrons.append(electron_2)
                    
                    p_bonding_orbitals_1[i].set_bond(self, 'pi')
                    p_bonding_orbitals_2[i].set_bond(self, 'pi')



    def calc_bond_spring_force(self):
        #calculate the distance between a bond's two neighbours
        r = self.calc_bond_length_model()
        bond_spring_force = self.cbond * (self.blength - r)
        return bond_spring_force

    def calc_angle_force(self, structure):
        r = self.calc_bond_length_model()
        
        true_neighbour = structure.get_true_neighbour(self.atom_1, self.atom_2)
        bond_angle = true_neighbour.get_bond_angle(structure)
        double_repulsion = False

        if bond_angle == 120:
            angle_force = 0.3 * (sqrt(3) * self.blength - r)
        elif bond_angle == 180:
            angle_force = 0.3 * (0.65 * blength - r)
        elif bond_angle == 90:
            angle_force = 0
            double_repulsion = True
        else:
            angle_force = 0

        return angle_force, double_repulsion
            
        
    def calc_bond_length_model(self):
        r = self.atom_1.calc_distance(self.atom_2)
        return r

    def calc_bond_length(self):
        atom_1_radius = ATOM_PROPERTIES.element_to_radius[atom_1.type]
        atom_2_radius = ATOM_PROPERTIES.element_to_radius[atom_2.type]
        bond_length = atom_1_radius + atom_2_radius
        return bond_length

class Shell:
    
    def __init__(self, atom, shell_nr):
        self.shell_nr = shell_nr
        self.orbital_sets = {}
        self.orbitals = {}
        self.bonding_orbitals = []
        self.atom = atom

        
        self.define_orbitals()

        self.find_bonding_orbitals()

    def define_orbitals(self):
        self.orbitals = {}
        self.orbital_sets[f'{self.shell_nr}s'] = OrbitalSet(self.atom, self.shell_nr, 's')
        if self.shell_nr >= 2:
            self.orbital_sets[f'{self.shell_nr}p'] = OrbitalSet(self.atom, self.shell_nr, 'p')
        if self.shell_nr >= 3:
            self.orbital_sets[f'{self.shell_nr}d'] = OrbitalSet(self.atom, self.shell_nr, 'd')
        if self.shell_nr >= 4:
            self.orbital_sets[f'{self.shell_nr}f'] = OrbitalSet(self.atom, self.shell_nr, 'f')

        for orbital_set in self.orbital_sets:
            for orbital in self.orbital_sets[orbital_set].orbitals:
                self.orbitals[orbital.__hash__()] = orbital

    def __hash__(self):
        return f'{self.atom.nr}_{self.shell_nr}'

    def __repr__(self):
        return f'{self.atom.nr}_{self.shell_nr}'

    def hybridise(self, hybridisation):
        if hybridisation == 'sp3':
            self.sp_hybridise(3)
        elif hybridisation == 'sp2':
            self.sp_hybridise(2)
        elif hybridisation == 'sp':
            self.sp_hybridise(1)
        elif hybridisation == 'sp3d':
            self.spd_hybridise(1)
        elif hybridisation == 'sp3d2':
            self.spd_hybridise(2)
        elif hybridisation == None:
            pass

        for orbital_name in self.orbitals:
            orbital = self.orbitals[orbital_name]
            for electron in orbital.electrons:
                if self.atom == electron.atom:
                    electron.set_orbital(orbital)

    def dehybridise(self):
        for orbital_set in self.orbital_sets:
            for i, orbital in enumerate(self.orbital_sets[orbital_set].orbitals):
                if orbital.orbital_type not in {'s', 'p', 'd', 'f'}:
                    new_orbital_type = self.orbital_sets[orbital_set].orbital_type
                    if new_orbital_type != 's':
                        new_orbital_nr = i + 1
                    else:
                        new_orbital_nr = None

                    orbital.orbital_type = new_orbital_type
                    orbital.orbital_nr = new_orbital_nr

        for orbital_name in self.orbitals:
            orbital = self.orbitals[orbital_name]
            for electron in orbital.electrons:
                if self.atom == electron.atom:
                    electron.set_orbital(orbital)


    def sp_hybridise(self, p_nr):
        #print("P nr:", p_nr, self.atom)
        hybridised_p = 0

        orbital_nr = 1
        
        for orbital_name in self.orbitals:
            orbital = self.orbitals[orbital_name]
            if orbital.orbital_type == 's':
                orbital.orbital_nr = orbital_nr
                orbital.orbital_type = f'sp{p_nr}'
                orbital_nr += 1
            elif orbital.orbital_type == 'p':
                if not orbital.bond or orbital.bonding_orbital == 'sigma':
                    if hybridised_p < p_nr:
                        orbital.orbital_type = f'sp{p_nr}'
                        orbital.orbital_nr = orbital_nr
                        hybridised_p += 1
                        orbital_nr += 1
            

                        
    def spd_hybridise(self, d_nr):
        hybridised_d = 0

        orbital_nr = 1

        for orbital_name in self.orbitals:
            orbital = self.orbitals[orbital_name]
            if orbital.orbital_type == 's':
                orbital.orbital_type = f'sp3d{d_nr}'
                orbital.orbital_nr = orbital_nr
                orbital_nr += 1
            if orbital.orbital_type == 'p':
                orbital.orbital_type = f'sp3d{d_nr}'
                orbital.orbital_nr = orbital_nr
                orbital_nr += 1
            elif orbital.orbital_type == 'd':
                if not orbital.bond or orbital.bonding_orbital == 'sigma':
                    if hybridised_d < d_nr:
                        orbital.orbital_type = f'sp3d{d_nr}'
                        hybridised_d += 1
                        orbital.orbital_nr = orbital_nr
                        orbital_nr += 1

    def excite(self):
        try:
            assert self.is_excitable()
        except AssertionError:
            raise SmilesError('violated_bonding_laws')
        electron_nr = self.count_electrons()
        
        for orbital_set in ATOM_PROPERTIES.orbital_order:
            if orbital_set in self.orbital_sets:
                for orbital in self.orbital_sets[orbital_set].orbitals:
                    for i in range(orbital.electron_nr):
                        orbital.empty_orbital()
                    
                    if electron_nr > 0:
                        orbital.fill_orbital()
                        electron_nr -= 1
                
                        
    def get_lone_pair_nr(self):
        lone_pair_nr = 0
        for orbital_name in self.orbitals:
            orbital = self.orbitals[orbital_name]
            if orbital.electron_nr == 2:
                if orbital.electrons[0].atom == orbital.electrons[1].atom:
                    lone_pair_nr += 1
        return lone_pair_nr

    def get_lone_electrons(self):
        lone_electrons = 0
        for orbital_name in self.orbitals:
            orbital = self.orbitals[orbital_name]
            if orbital.electron_nr == 1:
                lone_electrons += 1

        return lone_electrons

    def find_bonding_orbitals(self):
        self.bonding_orbitals = []
        
        for orbital in self.orbitals:
            if self.orbitals[orbital].electron_nr == 1:
                self.bonding_orbitals.append(orbital)
                
    def count_electrons(self):
        electron_nr = 0
        for orbital_name in self.orbitals:
            orbital = self.orbitals[orbital_name]
            electron_nr += orbital.electron_nr

        return electron_nr

    def count_orbitals(self):
        orbital_nr = len(list(self.orbitals.keys()))

        return orbital_nr
        

    def is_excitable(self):
        electron_nr = self.count_electrons()
        orbital_nr = self.count_orbitals()

        if orbital_nr >= electron_nr:
            return True
        else:
            return False

    def drop_electrons(self):

        
        lone_orbitals = []


        for orbital_set in ATOM_PROPERTIES.orbital_order:
            if orbital_set in self.orbital_sets:
                for orbital in self.orbital_sets[orbital_set].orbitals:
                    if orbital.electron_nr == 1:
                        if not orbital.electrons[0].aromatic:
                            lone_orbitals.append(orbital)


        while len(lone_orbitals) > 1 and (lone_orbitals[0].orbital_type != lone_orbitals[-1].orbital_type\
                or lone_orbitals[0].orbital_nr != lone_orbitals[-1].orbital_nr):

            receiver_orbital = lone_orbitals[0]
            donor_orbital = lone_orbitals[-1]

            moved_electron = donor_orbital.electrons[0]

            donor_orbital.remove_electron(moved_electron)
            receiver_orbital.add_electron(moved_electron)
            receiver_orbital.electrons[1].set_orbital(receiver_orbital)
            

            del lone_orbitals[0]
            del lone_orbitals[-1]


    def print_shell(self):
        for orbital in self.orbitals:
            print(self.orbitals[orbital])

            
class OrbitalSet():
    def __init__(self, atom, shell_nr, orbital_type):
        self.atom = atom
        self.shell_nr = shell_nr
        
        self.orbital_type = orbital_type
        self.orbitals = []
        self.define_orbitals()
        self.capacity = len(self.orbitals) * 2
        

    def __repr__(self):
        return f'{self.shell_nr}{orbital_type}'


    def define_orbitals(self):
        if self.orbital_type == 's':
            self.append_s_orbital()
        if self.orbital_type == 'p':
            self.append_p_orbitals()
        if self.orbital_type == 'd':
            self.append_d_orbitals()
        if self.orbital_type == 'f':
            self.append_f_orbitals()
        
    def append_s_orbital(self):
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 's'))
    
    def append_p_orbitals(self):
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'p', 1))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'p', 2))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'p', 3))

    def append_d_orbitals(self):
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'd', 1))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'd', 2))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'd', 3))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'd', 4))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'd', 5))
            
    def append_f_orbitals(self):
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'f', 1))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'f', 2))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'f', 3))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'f', 4))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'f', 5))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'f', 6))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'f', 7))

    def fill_orbitals(self, electrons):
        while electrons > 0:
            for orbital in self.orbitals:
                if electrons > 0:
                    orbital.fill_orbital()
                    electrons -= 1
                else:
                    break
        
class Orbital():
    subtype_dict = {'p': {1: 'x',
                          2: 'y',
                          3: 'z'},
                    'd': {1: 'z^2',
                          2: 'zx',
                          3: 'yz',
                          4: 'xy',
                          5: 'x^2-y^2'},
                    'f': {1: 'z^3-3/5zr^2',
                          2: 'x^3-3/5xr^2',
                          3: 'y^3-3/5yr^2',
                          4: 'xyz',
                          5: 'y(x^2-z^2)',
                          6: 'x(z^2-y^2)',
                          7: 'z(x^2-y^2)'}
                    }
    
    def __init__(self, atom, shell_nr, orbital_type, orbital_nr = None):
        self.shell_nr = shell_nr
        self.orbital_type = orbital_type
        self.orbital_nr = orbital_nr
        self.electron_nr = 0
        self.electrons = []
        self.atom = atom
        self.bond = None
        self.bonding_orbital = None


    def __hash__(self):
        if self.orbital_nr:                               
                                        
            return f'{self.shell_nr}{self.orbital_type}{self.orbital_nr}'
        
        else:
            return f'{self.shell_nr}{self.orbital_type}'

    def __repr__(self):
        if self.orbital_nr:                               
                                        
            return f'{self.shell_nr}{self.orbital_type}{self.orbital_nr}'
        
        else:
            return f'{self.shell_nr}{self.orbital_type}'

    def set_electron_nr(self):
        self.electron_nr = len(self.electrons)

    def set_bond(self, bond, bonding_orbital):
        self.bond = bond
        self.bonding_orbital = bonding_orbital

    def remove_bond(self):
        self.bond = None
        self.bonding_orbital = None

    def fill_orbital(self):
        """
        """
        assert self.electron_nr < 2

        self.electrons.append(Electron(self.shell_nr, self.orbital_type,
                                       self.orbital_nr, 0.5, self.atom))
        self.set_electron_nr()
        
        if self.electron_nr == 2:
            self.electrons[0].pair(self.electrons[1])

    def empty_orbital(self):
        """
        """
        assert self.electron_nr > 0
        
        del self.electrons[-1]
        self.set_electron_nr()

        if self.electron_nr == 1:
            self.electrons[0].unpair()
            
        
    def add_electron(self, electron):
        assert self.electron_nr < 2

        self.electrons.append(electron)
        self.set_electron_nr()
        

        if self.electron_nr == 2:
            self.electrons[0].pair(self.electrons[1])
            

    def remove_electron(self, electron):
        assert electron in self.electrons
        
        self.electrons.remove(electron)
        self.set_electron_nr()

        if self.electron_nr == 1:

            self.electrons[0].unpair()
 

class Atom:

    def __new__(cls, atom_type, atom_nr, chiral, charge, aromatic):
        self = super().__new__(cls)  # Must explicitly create the new object
        # Aside from explicit construction and return, rest of __new__
        # is same as __init__
        self.type = atom_type
        self.nr = atom_nr
        self.chiral = chiral
        self.charge = charge
        self.aromatic = aromatic
        self.shells = {}
        
        return self  # __new__ returns the new object

    def __getnewargs__(self):
        # Return the arguments that *must* be passed to __new__
        return (self.type, self.nr, self.chiral, self.charge, self.aromatic)

    def __init__(self, atom_type, atom_nr, chiral, charge, aromatic):
        
        self.type = atom_type
        self.nr = atom_nr
        self.aromatic = aromatic
 #       self.atomic_nr = ATOM_PROPERTIES.element_to_atomic_nr[self.type]
        self.bonds = []

        self.chiral = chiral
        # remove?
        self.charge = charge
        self.pyrrole = False
        self.shells = {}
        self.lone_pairs = []
        self.draw = AtomDrawProperties()
        
    def __eq__(self, atom):
        if type(atom) == Atom:
            return self.nr == atom.nr
        else:
            return False

    def __hash__(self):
        return self.nr

    def __repr__(self):
        if self.charge == 0:
            charge_string = ''
        elif self.charge > 0:
            if self.charge == 1:
                charge_string = '+'
            else:
                charge_string = str(self.charge) + '+'
        elif self.charge < 0:
            if self.charge == -1:
                charge_string = '-'
            else:
                charge_string = str(abs(self.charge)) + '-'
        
        return f'{self.type}{charge_string}_{self.nr}'

    def make_lone_pairs(self):

        lone_pair_nr = self.valence_shell.get_lone_pair_nr()
        for i in range(lone_pair_nr):
            self.lone_pairs.append(LonePair(self, i + 10000))



    def set_neighbours(self, structure):
        self.neighbours = structure.graph[self]

    def remove_neighbour(self, neighbour):
        self.neighbours.remove(neighbour)

    def set_connectivity(self):
        self.connectivity = self.get_connectivity()

    def get_connectivity(self):
        connectivity = []

        for bond in self.bonds:
            for atom in bond.neighbours:
                if atom.type != 'H' and atom != self:
                    bond_type = bond.type
                    connectivity.append(f'{atom.type}_{bond_type}')

        connectivity = tuple(sorted(connectivity))
        return connectivity

    def same_connectivity(self, atom):
        if self.type == atom.type:
            if len(self.connectivity) == len(atom.connectivity):
                if self.chiral == atom.chiral == None:
                    if set(self.connectivity) == set(atom.connectivity):
                        return True
                    else:
                        return False
                else:
                    pass
            else:
                return False
                

        else:
            return False

    def potential_same_connectivity(self, substructure_connectivity):
        parent_connectivity_copy = list(copy.copy(self.connectivity))
        substructure_connectivity_copy = list(copy.copy(substructure_connectivity))

        same_connectivity = True

        for atom in substructure_connectivity_copy:
            if atom in parent_connectivity_copy:
                parent_connectivity_copy.remove(atom)
            else:
                same_connectivity = False

        return same_connectivity

    def set_order(self):
        self.order = 0
        for neighbour in self.neighbours:
            if neighbour.type != 'H':
                self.order += 1

    def get_drawn_neighbours(self):
        drawn_neighbours = []
        for neighbour in self.neighbours:
            if neighbour.draw.is_drawn:
                drawn_neighbours.append(neighbour)

        return drawn_neighbours

            
            

    def add_shell_layout(self):
        self.shell_nr = ATOM_PROPERTIES.element_to_shell_nr[self.type]
        self.make_shells()
        self.fill_shells()

        if self.type == 'Fe':

            self.valence_shell.print_shell()
            for orbital_nr, orbital in self.valence_shell.orbitals.items():
                print(orbital.electrons)




        self.excitable = self.valence_shell.is_excitable()

        if self.type == 'C' and self.charge == 0:
            self.excite()
        else:
            bond_weights = [BOND_PROPERTIES.bond_type_to_weight[bond.type] for bond in self.bonds]
            aromatic_count = 0
            for bond in self.bonds:
                if bond.type == 'aromatic':
                    aromatic_count += 1


            nr_of_nonH_bonds = sum(bond_weights) + int(aromatic_count / 2)

            if self.pyrrole:
                nr_of_nonH_bonds -= 1

            if nr_of_nonH_bonds > ATOM_PROPERTIES.element_to_valences[self.type][0]:
                if self.excitable:
                    self.excite()
                else:
                    raise SmilesError('violated_bonding_laws')

    def make_shells(self):
        for i in range(self.shell_nr):
            current_shell = i + 1
            self.shells[current_shell] = Shell(self, current_shell)

        self.valence_shell = self.shells[self.shell_nr]
            

    def fill_shells(self):
        electrons_assigned = 0

        electrons_remaining = ATOM_PROPERTIES.element_to_atomic_nr[self.type] - self.charge
        if self.type == 'Fe':
            print('charge', self.charge)
            print(electrons_remaining)
        

        #Iterate over the orbitals in order of them being filled

        
        for orbital in ATOM_PROPERTIES.orbital_order:
            if electrons_remaining > 0:
                shell = int(orbital[0])
                orbital_set = self.shells[shell].orbital_sets[orbital]
                electrons_to_dump = min([electrons_remaining, orbital_set.capacity])
                orbital_set.fill_orbitals(electrons_to_dump)
                electrons_remaining -= electrons_to_dump
            else:
                break



    def excite(self):
        assert self.excitable
        
        self.valence_shell.excite()

    def remove_bond(self, bond):
        self.bonds.remove(bond)

    def calc_electron_pair_nr(self):

        bond_nr = self.calc_bond_nr()
        bonds_accounted_for = 0
        
        electron_nr = 0
        orbital_nr = len(list(self.valence_shell.orbitals.keys()))
        
        for orbital_name, orbital in self.valence_shell.orbitals.items():
            if orbital.electron_nr == 1:
                electron_nr += 1
            elif orbital.electron_nr == 2:
                if orbital.electrons[0].atom == orbital.electrons[1].atom:
                    electron_nr += 2
                else:
                    bonds_accounted_for += 1

        bonds_to_make = bond_nr - bonds_accounted_for

        unbonded_electrons = electron_nr - bonds_to_make

        if unbonded_electrons % 2 != 0:
            print("Warning! Rogue electron.")

        electron_pair_nr = int(unbonded_electrons / 2)
        
        return electron_pair_nr

    def drop_electrons(self):
        if self.valence_shell.get_lone_electrons() > 1:
            self.valence_shell.drop_electrons()
        

    def calc_bond_nr(self):

        bond_nr = 0
        aromatic_bond_nr = 0

        for bond in self.bonds:
            if bond.type == 'single':
                bond_nr += 1
            elif bond.type == 'double':
                bond_nr += 2
            elif bond.type == 'triple':
                bond_nr += 3
            elif bond.type == 'quadruple':
                bond_nr += 4
            elif bond.type == 'aromatic':
                aromatic_bond_nr += 1

        if aromatic_bond_nr == 2:
            if not self.pyrrole:
                bond_nr += 3
            else:
                bond_nr += 2
        elif aromatic_bond_nr == 3:
            bond_nr += 4

        return bond_nr

    def is_promotable(self):
        promotable = False
        for orbital_set in self.valence_shell.orbital_sets:
            if 'd' in orbital_set:
                promotable = True

        return promotable

    def promote_lone_pair_to_p_orbital(self):
        assert self.calc_electron_pair_nr() > 0
        assert self.hybridisation == 'sp3'

        self.valence_shell.dehybridise()
        #self.valence_shell.hybridise('sp2')

        p_orbitals = []
        sp2_orbitals = []

        for orbital_name, orbital in self.valence_shell.orbitals.items():
            if orbital.electron_nr == 2:
                if orbital.electrons[0].atom != orbital.electrons[1].atom:
                    sp2_orbitals.append(orbital)
                elif orbital.electrons[0].atom == orbital.electrons[1].atom == self:
                    p_orbitals.append(orbital)

        if len(p_orbitals) > 1:
            for i in range(1, len(p_orbitals)):
                sp2_orbitals.append(p_orbitals[i])

        p_orbital = p_orbitals[0]

        p_orbital.orbital_type = 'p'

        for orbital in sp2_orbitals:
            orbital.orbital_type = 'sp2'

        self.hybridisation = 'sp2'

        for orbital_name, orbital in self.valence_shell.orbitals.items():
            for electron in orbital.electrons:
                if electron.atom == self:
                    electron.set_orbital(orbital)
                    if orbital.orbital_type == 'p':
                        electron.set_aromatic()
            

    def promote_pi_bond_to_d_orbital(self):
        assert self.is_promotable()
            
        donor_orbitals = []
        receiver_orbitals = []
        for orbital_name in self.valence_shell.orbitals:
            orbital = self.valence_shell.orbitals[orbital_name]
            if orbital.orbital_type == 'p' and orbital.electron_nr == 2:
                if orbital.electrons[0].atom != orbital.electrons[1].atom:
                    donor_orbitals.append(orbital)

            elif orbital.orbital_type == 'd' and orbital.electron_nr == 1:
                receiver_orbitals.append(orbital)

        donor_orbital = donor_orbitals[0]
        receiver_orbital = receiver_orbitals[0]

        moved_electron = None

        for electron in donor_orbital.electrons:
            if electron.atom != self:
                moved_electron = electron


        donor_orbital.remove_electron(moved_electron)
        receiver_orbital.add_electron(moved_electron)

        receiver_orbital.set_bond(donor_orbital.bond, 'pi')
        donor_orbital.remove_bond()

        self.valence_shell.dehybridise()

        self.hybridise()


    def calc_hydrogens(self):
        hydrogens = 0
        if self.type in ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']:

            bond_nr = self.calc_bond_nr()
            if bond_nr in ATOM_PROPERTIES.element_to_valences[self.type]:
                hydrogens = 0
            else:
                max_bonds = self.valence_shell.get_lone_electrons()
                hydrogens = max_bonds - bond_nr

        return hydrogens

    def add_bond(self, bond):
        self.bonds.append(bond)

    def calc_repulsion_force(self, atom, crep, repulsion_threshold):
        self_coords = self.get_coords()
        atom_coords = atom.get_coords()
        
        distance = squared_distance(self_coords, atom_coords)

        if repulsion_threshold > distance > 0.1:
            vector = []
            
            for i, coord in enumerate(self_coords):
                atom_coord = atom_coords[i]
                diff = atom_coord - coord
                vector.append(crep * diff / distance)

        elif distance <= 0.1:
            vector = force_to_vector(1.0, atom, self)
 #           print('stom_1', self, self.get_coords())
#            print('atom_2', atom, atom.get_coords())
#            print(vector)

        else:
            vector = [0.0, 0.0, 0.0]

        return vector

    def hybridise(self):
        hybridisation = self.get_hybridisation()
        self.valence_shell.hybridise(hybridisation)
        self.set_hybridisation()

    def set_hybridisation(self):
        self.hybridisation = 's'
        for orbital_name in self.valence_shell.orbitals:
            orbital = self.valence_shell.orbitals[orbital_name]
            if orbital.orbital_type in {'sp', 'sp2', 'sp3', 'sp3d', 'sp3d2'}:
                self.hybridisation = orbital.orbital_type
                break

    def get_hybridisation(self):
        steric_number = self.get_steric_number()
        #Make dict

        if steric_number == 1:
            hybridisation = 's'
        elif steric_number == 2:
            hybridisation = 'sp'
        elif steric_number == 3:
            hybridisation = 'sp2'
        elif steric_number == 4:
            hybridisation = 'sp3'
        elif steric_number == 5:
            hybridisation = 'sp3d'
        elif steric_number == 6:
            hybridisation = 'sp3d2'
        elif steric_number == 0:
            hybridisation = None
        

        return hybridisation

    def get_steric_number(self):
        return self.calc_electron_pair_nr() + len(self.bonds) 

    def get_valence(self):
        try:
            valence = ATOM_PROPERTIES.element_to_valences[self.type][0]
            return valence

        except KeyError:
            print("Unknown atom: %s" % self.type)
            return "Unknown atom: %s" % self.type

    def get_coords(self):
        return [self.x, self.y, self.z]

    def get_bond_angle(self, structure):
        non_dummy = []
        for atom in structure.graph[self]:
            if structure.bond_lookup[self][atom].type != 'dummy':
                non_dummy.append(atom)

        if len(non_dummy) == 3:
            angle = 120
        elif len(non_dummy) == 4:
            angle = 90
        elif len(non_dummy) == 2:
            angle = 180
        elif len(non_dummy) == 1:
            angle = 0
        else:
            angle = None

        return angle

    def move_atom(self, vector):
        current_position = self.get_coords()
        new_position = [self.x + vector[0],
                        self.y + vector[1],
                        self.z + vector[2]]
        distance = euclidean(current_position, new_position)
        if distance > 1.5:
            x_displacement = vector[0] * 1.5 / distance
            y_displacement = vector[1] * 1.5 / distance
            z_displacement = vector[2] * 1.5 / distance
        else:
            x_displacement = vector[0]
            y_displacement = vector[1]
            z_displacement = vector[2]

        self.x += x_displacement
        self.y += y_displacement
        self.z += z_displacement

    def get_hydrogen_nr(self, structure):
        hydrogen_count = 0
        for atom in structure.graph[self]:
            if atom.type == 'H':
                hydrogen_count += 1

        return hydrogen_count
            
                

    def calc_distance(self, atom):
        distance = squared_distance(self.get_coords(), atom.get_coords())
        return distance

class LonePair:
    def __init__(self, atom, nr):
        self.nr = nr
        self.parent = atom

class Electron:
    def __init__(self, shell_nr, orbital_type, orbital_nr, spin, atom):
        self.shell_nr = shell_nr
        self.orbital_type = orbital_type
        self.orbital_nr = orbital_nr
        self.spin = spin
        self.atom = atom
        self.paired = False
        self.partner = None
        self.aromatic = False

    def __repr__(self):

        if self.aromatic:
            aromatic_string = '*'
        else:
            aromatic_string = ''
        if self.orbital_nr:
            return f'{self.atom}_{self.shell_nr}{self.orbital_type}{self.orbital_nr}_{self.spin}{aromatic_string}'
        else:
            return f'{self.atom}_{self.shell_nr}{self.orbital_type}_{self.spin}{aromatic_string}'

    def __hash__(self):
        if self.orbital_nr:
            return f'{self.atom}_{self.shell_nr}{self.orbital_type}{self.orbital_nr}_{self.spin}'
        else:
            return f'{self.atom}_{self.shell_nr}{self.orbital_type}_{self.spin}'
        
    def set_orbital(self, orbital):
        self.orbital_type = orbital.orbital_type
        self.orbital_nr = orbital.orbital_nr
        self.orbital = orbital

    def set_paired(self):
        self.paired = True

    def set_unpaired(self):
        self.paired = False

    def set_aromatic(self):
        self.aromatic = True

    def set_unaromatic(self):
        self.aromatic = False

    def pair(self, electron):
        self.spin = 0.5
        electron.spin = -0.5
        
        self.set_paired()
        electron.set_paired()
        self.partner = electron
        electron.partner = self

    def unpair(self):
        self.set_unpaired()
        self.partner = None
        self.spin = 0.5
                                                    

def calc_charge(sign, value):
    if sign == '+':
        return value
    elif sign == '-':
        return value * -1
    else:
        raise Exception("Wrong character to indicate charge!")

def parse_explicit(component):
    skip = False
    informative = component[1:-1]
    
    charges = []
    hydrogen = None
    numbers = []
    element = []
    chirals = []

    hydrogens = 0
    chiral = None
    charge = 0


    for i, character in enumerate(informative):
        if skip:
            skip = False
            continue
        if character.isupper():
            
            if character == 'H':
                hydrogen = i
            else:
                try:
                    if informative[i + 1].islower():
                        element.append(i)
                        element.append(i + 1)
                        skip = True
                    else:
                        element.append(i)
                except IndexError:
                    element.append(i)
        elif character.islower():
            element.append(i)
        elif character.isdigit():
            numbers.append(i)
        elif character == '+' or character == '-':
            charges.append(i)
        elif character == '@':
            chirals.append(i)

    element = ''.join([informative[x] for x in element])

    #Parsing the charge

    if len(charges) == 1:
        index = charges[0]
        charge_type = informative[index]
        try:
            if (index + 1) in numbers:
                charge_value = int(informative[index + 1])
                charge = calc_charge(charge_type, charge_value)
            else:
                charge = calc_charge(charge_type, 1)
                
        except IndexError:
            charge = calc_charge(charge_type, 1)

    elif len(charges) == 0:
        charge = 0
    else:
        charge_type = informative[charges[0]]
        charge_value = len(charges)
        charge = calc_charge(charge_type, charge_value)

    #Parsing an explicit hydrogen

    if hydrogen != None:
        try:
            if (hydrogen + 1) in numbers:
                hydrogens = int(informative[hydrogen + 1])
            else:
                hydrogens = 1
        except IndexError:
            hydrogens = 1
    else:
        hydrogens = 0


    #Parsing chirality

    if len(chirals) == 1:
        chiral = 'counterclockwise'
    elif len(chirals) == 2:
        chiral = 'clockwise'
    else:
        chiral = None

    #dealing with lone explicit hydrogens

    if not element and hydrogen == 0:
        element = 'H'
        hydrogens = 0

    return element, chiral, charge, hydrogens


def make_character_dict() -> Dict[str, str]:
    """Create dict of {character: label, ->} to label smiles characters
    """
    character_dict = {}
    atoms = ["C", "O", "N", "S", "B", "P", "F", "I", "c", "n", "o", '*',
             'Cl', 'Br', 'p', 'b', 'p', 's']
    cyclic = list(range(1,100))
    

    for atom in atoms:
        character_dict[atom] = "atom"
    for number in cyclic:
        character_dict[str(number)] = "cyclic"

    character_dict["="] = "double_bond"
    character_dict["("] = "branch_start"
    character_dict[")"] = "branch_end"
    character_dict['\\'] = 'chiral_double_bond'
    character_dict['/'] = 'chiral_double_bond'
    character_dict['#'] = 'triple_bond'
    character_dict['$'] = 'quadruple_bond'
    character_dict['.'] = 'split'
    character_dict['-'] = 'single_bond'
    character_dict[':'] = 'aromatic_bond'
    
    return character_dict

class Smiles():
    character_dict = make_character_dict()
    two_atom_dict = {'B': {'r'}, 'C': {'l'}}
    
            
    def __init__(self, string: str) -> None:
        self.smiles = string
        self.get_components()
        

    def get_components(self):
        self.components = []

        skip = False
        double_digits = False
        square_brackets = False
        
        for i, character in enumerate(self.smiles):
            if skip:
                skip = False
            elif square_brackets:
                component += character
                if character == ']':
                    square_brackets = False
                    self.components.append(component)
                    component = ''
            elif double_digits:
                try:
                    next_character = self.smiles[i + 1]
                    
                    if next_character not in {'0', '1', '2', '3', '4',
                                              '5', '6', '7', '8', '9'}:
                        double_digits = False
                        self.components.append(component + character)
                        component = ''
                    else:
                        component += character
                except IndexError:
                    self.components.append(component + character)
            else:
                
                if character in self.two_atom_dict:
                    try:
                        next_character = self.smiles[i + 1]
                        if next_character in self.two_atom_dict[character]:
                            self.components.append(character + next_character)
                            skip = True
                        else:
                            self.components.append(character)
                            
                        
                    except IndexError:
                        self.components.append(character)
                elif character == '[':
                    square_brackets = True
                    component = character
                elif character == '%':
                    double_digits = True
                    component = ''
                else:
                    self.components.append(character)
            

    def smiles_to_structure(self):

        bond_labels = {'single_bond',
                       'double_bond',
                       'triple_bond',
                       'quadruple_bond',
                       'chiral_double_bond',
                       'aromatic_bond'}

        structure = Structure()

        #Keeps track of the layer of brackets the atom is in

        branch_level = 0

        #Keeps track of which cyclic atoms have been encountered
        cyclic_dict = {}

        #Keeps track of the last atom encountered per branch level
        last_atoms_dict = {0: None}


        #Keeps track of the nature of the bond
        bond_type = 'single'
        bond_chiral_symbol = None

        #Keeps track of chiral centres
        chiral_dict = {}
        

        explicit = False
        pyrrole = False

        atom_nr = -1
        bond_nr = -1
        
        for i, component in enumerate(self.components):
            if component[0] == '[':
                label = 'atom'
                explicit = True
            else:
                label = self.character_dict[component]

            if label in bond_labels:
                try:
                    next_component = self.components[i + 1]
                    if next_component[0] == '[':
                        next_label = 'atom'
                    else:
                        next_label = self.character_dict[next_component]

                    if next_label != 'atom' and next_label != 'cyclic':
                        raise SmilesError('bond')
                except IndexError:
                    raise SmilesError('bond')

            if label == 'split':
                #Starts disconnected structure; set everything back to default
                branch_level = 0
                cyclic_dict = {}
                last_atoms_dict = {0: None}
                double = False
                triple = False
                chiral_dict = {}

            elif label == "atom":

                

                if not explicit:
                    element = component
                    chiral = None
                    charge = 0
                    hydrogens = 0
                else:
                    element, chiral, charge, hydrogens = parse_explicit(component)
                    if element == 'n' and hydrogens == 1:
                        pyrrole = True
                    explicit = False

                if element.islower():
                    aromatic = True
                    element = element.upper()
                else:
                    aromatic = False

                #Determine atom
                
                atom_nr += 1

                #Turn atom into an atom object
                atom_2 = Atom(element, atom_nr, chiral, charge, aromatic)
                if pyrrole:
                    atom_2.pyrrole = True
                    pyrrole = False

                    
                for i in range(hydrogens):
                    atom_nr += 1
                    bond_nr += 1
                    hydrogen = Atom('H', atom_nr, None, 0, False)
                    structure.add_bond(atom_2, hydrogen, 'single', bond_nr)

                
                        
                    
                atom_1 = self.get_last_atom(branch_level, last_atoms_dict)

                previous_atom_branch_level = branch_level

                #Go back through the branches to identify the atom that the
                #current atom is attached to; call this atom_1
                
                while previous_atom_branch_level > 0 and not atom_1:
                    previous_atom_branch_level -= 1
                    atom_1 = last_atoms_dict[previous_atom_branch_level]

                if atom_1:
                    bond_nr += 1
                    
                    if atom_1.aromatic and atom_2.aromatic:
                        if not bond_type == 'explicit_single':
                            bond_type = 'aromatic'
                        else:
                            bond_type = 'single'

                    if bond_type == 'single_chiral':
                        bond_type = 'single'
                        
                    structure.add_bond(atom_1, atom_2, bond_type, bond_nr, bond_chiral_symbol)
                    

                    if atom_1.chiral:
                        #Add current atom to list of residues
                        chiral_dict[atom_1].append(atom_2)                          
                            
                    if atom_2.chiral:
                        chiral_dict[atom_2] = [atom_1]
                        
                        if hydrogens == 1:
                            chiral_dict[atom_2].append(hydrogen)

                    bond_type = 'single'
                    bond_chiral_symbol = None

                            
                            

                else:
                    
                    #This happens only if atom_2 is completely disconnected
                    #of any atoms already in the graph
                    structure.add_disconnected_atom(atom_2)
                    if atom_2.chiral:
                        chiral_dict[atom_2] = []
                        if hydrogens == 1:
                            chiral_dict[atom_2].append(hydrogen)

                #Set atom_2 as the last atom in the current branch level
                self.track_last_atoms_per_branch(atom_2, branch_level,
                                                 last_atoms_dict)

                
            elif label == "single_bond":
                bond_type = 'explicit_single'

            elif label == 'aromatic_bond':
                bond_type = 'aromatic'
                
            elif label == "double_bond":
                bond_type = 'double'
                
            elif label == "triple_bond":
                bond_type = 'triple'
                
            elif label == "quadruple_bond":
                bond_type = 'quadruple'

            #If there are brackets: go up or down a branch level
            elif label == "branch_start":
                branch_level += 1
                last_atoms_dict[branch_level] = None

            elif label == "branch_end":
                last_atoms_dict[branch_level] = None
                branch_level -= 1

            elif label == "cyclic":

                cycle_nr = int(component)
                atom = self.get_last_atom(branch_level, last_atoms_dict)
                

                #If we encounter a number that hasn't been encountered before
                
                if self.is_new_cycle(cyclic_dict, cycle_nr):
                    self.start_cycle(cycle_nr, atom, cyclic_dict)
                    
                    if atom in chiral_dict:
                        self.add_cycle_placeholder(chiral_dict, atom, cycle_nr)

                #Otherwise look up the atom that the cycle closes on
                else:
                    bond_nr += 1
                    
                        
                    atom_1, atom_2 = self.end_cycle(cycle_nr, atom,
                                                    cyclic_dict)
                    
                    if atom_1.aromatic and atom_2.aromatic:
                        if not bond_type == 'explicit_single':
                            bond_type = 'aromatic'
                        else:
                            bond_type = 'single'

                    if bond_type == 'single_chiral':
                        bond_type = 'single'
                        if bond_chiral_symbol == '/':
                            bond_chiral_symbol = '\\'
                        elif bond_chiral_symbol == '\\':
                            bond_chiral_symbol = '/'

                        structure.add_bond(atom_1, atom_2, bond_type, bond_nr, bond_chiral_symbol)
                        bond_chiral_symbol = None

                    else:
                        
                        structure.add_bond(atom_1, atom_2, bond_type, bond_nr)
                    
                    
                    if atom_1 in chiral_dict:
                        self.replace_cycle_placeholder(chiral_dict, atom_1, atom_2, cycle_nr)
                    if atom_2 in chiral_dict:
                        chiral_dict[atom_2].append(atom_1)

                    bond_type = 'single'

            elif label == 'chiral_double_bond':
                bond_type = 'single_chiral'
                bond_chiral_symbol = component




        structure.refine_structure()
        structure.set_double_bond_chirality()

        for atom in chiral_dict:
            order = chiral_dict[atom]

            #Save chirality in the atom object
            current_chirality = atom.chiral
            new_chirality = self.determine_chirality(order, current_chirality, atom)


            atom.chiral = new_chirality

        return structure



    def kekulize(self) -> str:
        """Return kekulized smiles from smiles.

        Input:
        smiles: str, smiles string

        Output:
        K_smiles: str, smiles string, kekulized version of smiles
        """
        
        return K_smiles

 

    def draw_smiles(self, img_dir: str, highlight_map: Dict[int, List[float]], bond_list: List[int], bond_colour: List[float], ID: str) -> None:
        """Create image from smiles string

        Input:
        img_dir: str, directory of output image
        highlight_map: dict of {atom_nr: [float_1, float_2, float_3], ->}, with
            atom_nr int and the combination of floats representing an RGB colour
        bond_list: list of int, with each int a bond index of the molecule
        bond_colour: list of [float_1, float_2, float_3], representing an
            RGB colour
        ID: str, name of the molecule
        """

        mol = Chem.MolFromSmiles(self.smiles)
        mol = Chem.Mol(mol.ToBinary())

        if not mol.GetNumConformers():
            rdDepictor.Compute2DCoords(mol)

        colormap = dict(zip(bond_list, [bond_colour] * len(bond_list)))
        
        print(colormap)        
        
        size = (500, 500)
 #       drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
#        drawer.DrawMolecule(mol, highlightAtoms = list(highlight_map.keys()),
#                            highlightAtomColors = highlight_map,
#                            highlightBonds = bond_list,
#                            highlightBondColors=colormap)
#        drawer.FinishDrawing()
        
 #       img = Chem.Draw.MolToImage(mol, size=size, kekulize=True,
#                                   wedge_bonds = True, fitImage = True,
#                                   highlightMap = highlight_map,
#                                   highlightBonds = bond_list,
#                                   highlightColor = bond_colour)

        img = Chem.Draw.MolToFile(mol, "%s/%s.svg" % (img_dir, ID),
                                   size=size, kekulize=True,
                                   wedge_bonds = True, fitImage = True,
                                   highlightMap = highlight_map,
                                   highlightColor = bond_colour,
                                   highlightBonds = bond_list)

        
 #       img.save("%s/%s.png" % (img_dir, ID))
#        svg = drawer.GetDrawingText()
 #       with open("%s/%s.svg" % (img_dir, ID), 'w') as f:
#            f.write(svg)
        

    #=========================================================================
    #Auxillary functions to smiles_to_structure
    #=========================================================================
        

    def add_chiral_atom(self, chiral_dict: Dict['Atom', Dict[str, Any]], last_atom: 'Atom', current_atom: 'Atom') -> None:
        """Place current_atom in one of the four bond slots of last_atom

        Input:
        chiral_dict: dict of {atom: {'direction': direction, 'order':
            [atom_1, atom_2, atom_3, atom_4]}}, with atom Atom Object,
            direction int, and order a list of Atom Object or int or None
        last_atom: Atom Object, chiral atom
        current_atom: Atom Object, neighbour of chiral atom

        """


        chiral_dict[last_atom].append(current_atom)

    def add_cycle_placeholder(self, chiral_dict, atom, cycle_nr):
        chiral_dict[atom].append(cycle_nr)

    def replace_cycle_placeholder(self, chiral_dict, chiral_atom, current_atom, cycle_nr):
        for i, atom in enumerate(chiral_dict[chiral_atom]):
            if type(atom) == int:
                if atom == cycle_nr:
                    chiral_dict[chiral_atom][i] = current_atom
        

    def get_chiral_permutations(self, order):
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
        
        
        
    def determine_chirality(self, order, chirality, atom):
        if len(order) != 4:
            lone_pairs = atom.lone_pairs

            try:
                assert len(order) + len(lone_pairs) == 4
            except AssertionError:
                raise SmilesError('chiral centre')

            original_order = tuple(order + lone_pairs)

        else:
            original_order = tuple(order)

        chiral_permutations = self.get_chiral_permutations(original_order)
        order.sort(key = lambda x: x.nr)
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



    def is_new_cycle(self, cyclic_dict: Dict[int, 'Atom'], cycle_nr: int) -> bool:
        """Return bool, True if a new cycle is recorded, False if not

        Input:
        cyclic_dict: dict of {cycle_nr: atom, ->}, with cycle_nr int and atom
            an Atom Object
        cycle_nr: int, nr of the current cycle

        Output:
        bool: True if an atom with cycle_nr at position 0 does not yet exist in
            cyclic_atoms, False if it does
        """
        if cycle_nr in cyclic_dict:
            return False
        else:
            return True

    def start_cycle(self, cycle_nr: int, atom: 'Atom', cyclic_dict: Dict[int, 'Atom']) -> None:
        """Add a new atom and corresponding cycle number to cyclic dict

        Input:
        cycle_nr: int, nr of the current cycle
        atom: Atom Object
        cyclic_dict: dict of {cycle_nr: atom, ->}, with cycle_nr int and atom
            an Atom Object

        """
        cyclic_dict[cycle_nr] = atom

    def end_cycle(self, cycle_nr: int, atom: 'Atom', cyclic_dict: Dict[int, 'Atom']) -> Tuple['Atom']:
        """Return pair of atoms that close a cycle

        Input:
        cycle_nr: int, nr of the current cycle
        atom: Atom Object
        cyclic_dict: dict of {cycle_nr: atom, ->}, with cycle_nr int and atom
            an Atom Object

        Output:
        atom_pair: tuple of two atoms, with each atom an Atom Object
        
        """
        atom_old = cyclic_dict[cycle_nr]
        atom_pair = (atom_old, atom)
        
        del cyclic_dict[cycle_nr]

        #Implement atom pair here?
        
        return atom_pair

    def track_last_atoms_per_branch(self, new_atom: 'Atom', current_level: int, last_atoms_dict: Dict[int, 'Atom']) -> None:
        """Update which atom was the last to occur at the current branch level

        Input:
        new_atom: Atom Object
        current_level: int, current branch level, level to which new_atom is
            added
        last_atoms_dict: dict of {int: atom, ->}, with int representing a branch
            level, and atom representing the last atom that occurs in that
            branch. 
        """
        last_atoms_dict[current_level] = new_atom

    def get_last_atom(self, current_level: int, last_atoms_dict: Dict[int, 'Atom']) -> 'Atom':
        """Return the last atom in the current branch level

        Input:
        current_level: int, current branch level, level from which the last atom
            is to be extracted
        last_atoms_dict: dict of {int: atom, ->}, with int representing a branch
            level, and atom representing the last atom that occurs in that
            branch.

        Output:
        last_atom: Atom Object, last atom that was encountered in the
            current_level branch
        """
        last_atom = last_atoms_dict[current_level]
        return last_atom
        
    def check_character(self, character: str) -> str:
        """Return the label for a character to determine its type

        Input:
        character: str, substring of smiles string of length 1

        Output:
        label: str, descriptor of character
        """

        try:
            #Consider making character_dict a class attribute defined before
            #__init__
            label = character_dict[character]
        except KeyError:
            label = "Unknown"
        return label

    def update_structure(self, atom_1: 'Atom', atom_2: 'Atom', structure_graph: Dict['Atom', List['Atom']], bond_type: str) -> None:
        """Add an atom to the structure graph

        Input:
        atom_1: Atom Object
        atom_2: Atom Object
        structure_graph: dict of {atom: [atom, ->], ->}, with each atom an
        Atom Object
        """
        if atom_1 in structure_graph:
            structure_graph[atom_1] += [atom_2]
        else:
            structure_graph[atom_1] = [atom_2]

        if atom_2 in structure_graph:
            structure_graph[atom_2] += [atom_1]
        else:
            structure_graph[atom_2] = [atom_1]


#=============================================================================
#Code adapted from https://github.com/yorkyer/edmonds-blossom to find perfect
# matchings for kekulisation
#=============================================================================


class Node:
    index = 0

    def __init__(self):
        self.atom = None
        self.neighbors = []
        self.is_visited = False
        self.parent = None
        self.mate = None
        self.index = Node.index
        Node.index += 1

    def __repr__(self):
        return str(self.atom)

    def set_atom(self, atom):
        self.atom = atom


class SuperNode(Node):

    def __init__(self):
        Node.__init__(self)
        self.subnodes = []
        self.original_edges = []

    def circle(self, node):
        for i, v in enumerate(self.subnodes):
            if v == node:
                break
        assert i < len(self.subnodes)

        if i > 0 and self.subnodes[i].mate == self.subnodes[i - 1] or i == 0 and self.subnodes[i].mate == self.subnodes[
            -1]:
            return self.subnodes[i::-1] + self.subnodes[:i:-1]
        else:
            return self.subnodes[i::] + self.subnodes[:i]


class Path:

    def __init__(self):
        self.nodes = []

    def head(self):
        return self.nodes[0]

    def tail(self):
        return self.nodes[-1]

    def replace_head(self):
        assert isinstance(self.nodes[0], SuperNode)
        snode = self.nodes.pop(0)
        for node in snode.subnodes:
            if self.nodes[0] in node.neighbors:
                if node.mate == None:
                    self.nodes.insert(0, node)
                else:
                    for v in snode.circle(node):
                        self.nodes.insert(0, v)
                return
        logging.error("cannot replace head node.")

    def replace_tail(self):
        assert isinstance(self.nodes[-1], SuperNode)
        snode = self.nodes.pop()
        for node in snode.subnodes:
            if self.nodes[-1] in node.neighbors:
                if node.mate == None:
                    self.nodes.append(node)
                else:
                    for v in snode.circle(node):
                        self.nodes.append(v)
                return
        logging.error("cannot replace tail node.")

    def __repr__(self):
        return str(self.nodes)


class Match:

    def __init__(self, nodes):
        self.nodes = nodes
        self.freenodes = []
        for node in nodes:
            self.freenodes.append(node)
        self.supernodes = []

    def from_structure(structure):
        pi_subgraph = structure.find_pi_subgraph()
        node_nr = len(pi_subgraph.keys())

        atom_to_node = {}
        nodes = [Node() for i in range(node_nr)]

        for i, atom in enumerate(pi_subgraph):
            atom_to_node[atom] = i

        for atom, neighbours in pi_subgraph.items():
            index_1 = atom_to_node[atom]
            nodes[index_1].set_atom(atom)

            for neighbour in neighbours:
                index_2 = atom_to_node[neighbour]
                nodes[index_1].neighbors.append(nodes[index_2])


        return Match(nodes)

    def from_edges(N, edges):
        nodes = [Node() for i in range(N)]
        for i, j in edges:
            nodes[i].neighbors.append(nodes[j])
        return Match(nodes)

    def clear_nodes(self):
        for node in self.nodes:
            node.is_visited = False
            node.parent = None

    def find_augmenting_path(self, root):
        self.clear_nodes()
        queue = [root]
        while len(queue) > 0:
            cur_node = queue.pop(0)
            cur_node.is_visited = True
            for node in cur_node.neighbors:
                if node == cur_node.parent:
                    continue

                elif node.is_visited:
                    cycle = self.find_cycles(node, cur_node)
                    if len(cycle) % 2 == 1:
                        logging.debug('blossom: {}'.format(cycle))
                        snode = self.shrink_blossom(cycle)
                        self.supernodes.append(snode)
                        for v in cycle:
                            if v in queue:
                                queue.remove(v)
                            if v.is_visited:
                                snode.is_visited = True
                                snode.parent = v.parent
                        queue.append(snode)
                        break

                elif node.mate == None:
                    node.parent = cur_node
                    return self.construct_augmenting_path(node)

                else:
                    node.is_visited = True
                    node.mate.is_visited = True
                    node.parent = cur_node
                    node.mate.parent = node
                    queue.append(node.mate)
        raise Exception('cannot find an augmenting path')

    def unmatched_nodes(self):
        self.maximum_matching()

        count = 0
        for node in self.nodes:
            if node.mate != None:
                count += 1

        return len(self.nodes) - count

    def maximum_matching(self):
        while len(self.freenodes) > 0:
            logging.debug('freenodes: {}'.format(self.freenodes))

            for node in self.freenodes:
                try:
                    path = self.find_augmenting_path(node)
                    logging.debug('augmenting path: {}'.format(path.nodes))
                    self.invert_path(path)
                    self.freenodes.remove(path.nodes[0])
                    self.freenodes.remove(path.nodes[-1])
                    break
                except Exception as e:
                    logging.info(e)
            else:
                logging.info('Tried all free nodes, no more augmenting path.')

                break


    def invert_path(self, path):
        assert len(path.nodes) % 2 == 0
        for i in range(0, len(path.nodes), 2):
            path.nodes[i].mate = path.nodes[i + 1]
            path.nodes[i + 1].mate = path.nodes[i]

    def construct_augmenting_path(self, node):
        path = Path()
        path.nodes.append(node)
        node = node.parent
        path.nodes.append(node)
        while node.mate != None:
            node = node.parent
            path.nodes.append(node)

        while len(self.supernodes) > 0:
            snode = self.supernodes.pop()
            self.expand_supernode(snode)
            if snode == path.head():
                path.replace_head()
            elif snode == path.tail():
                path.replace_tail()

        while path.nodes[0].mate != None:
            path.nodes.insert(path.nodes[0].parent, 0)

        while path.nodes[-1].mate != None:
            path.nodes.append(path.nodes[-1].parent)

        return path

    def find_cycles(self, node1, node2):
        def find_ancestors(node):
            ancestors = [node]
            while node.parent != None:
                node = node.parent
                ancestors.append(node)
            return ancestors

        ancestors1 = find_ancestors(node1)
        ancestors2 = find_ancestors(node2)
        i = len(ancestors1) - 1
        j = len(ancestors2) - 1
        while ancestors1[i] == ancestors2[j]:
            i -= 1
            j -= 1

        cycle = ancestors1[:i + 1] + ancestors2[j + 1::-1]
        return cycle

    def shrink_blossom(self, blossom):
        snode = SuperNode()
        for node in blossom:
            snode.subnodes.append(node)
            for adj_node in node.neighbors:
                if adj_node not in blossom:
                    snode.original_edges.append((node, adj_node))

        for node1, node2 in snode.original_edges:
            node1.neighbors.remove(node2)
            node2.neighbors.remove(node1)
            node2.neighbors.append(snode)
            snode.neighbors.append(node2)

        return snode

    def expand_supernode(self, snode):
        assert isinstance(snode, SuperNode)
        for node1, node2 in snode.original_edges:
            node1.neighbors.append(node2)
            node2.neighbors.append(node1)
            node2.neighbors.remove(snode)
            snode.neighbors.remove(node2)


#=============================================================================

class AtomDrawProperties:
    def __init__(self, x=0, y=0):
        self.rings = []
        self.original_rings = []
        self.ring_anchors = set()
        self.is_bridge_atom = False
        self.is_bridge = False
        self.bridged_ring = None
        self.is_drawn = True
        self.has_hydrogen = False
        self.positioned = False
        self.previous_position = Vector(0, 0)
        self.position = Vector(x, y)
        self.angle = None
        self.force_positioned = False
        self.connected_to_ring = False
        self.draw_explicit = False
        self.previous_atom = None

    def set_position(self, vector):
        self.position = vector

    def set_previous_position(self, previous_atom):
        self.previous_position = previous_atom.previous_position
        self.previous_atom = previous_atom

    def get_angle(self, reference_vector=None):
        if not reference_vector:
            vector = Vector.subtract_vectors(self.position, self.previous_position)
        else:
            vector = Vector.subtract_vectors(self.position, reference_vector)

        return vector.angle()

class BondDrawProperties:
    def __init__(self):
        self.center = False



if __name__ == "__main__":

    #structure = read_smiles("[C-]#[N+]/C=C\\C1=CNC2=CC=CC=C21")
    #for bond_nr, bond in structure.bonds.items():
     #   if bond.chiral:
     #       pprint(bond.chiral_dict)

    structure = read_smiles(r'C=C(/C=C/C1=CC=CC2=C1NC=C2CC(N3O[Fe](ON4C(CC5=CNC6=C5C=CC=C6/C=C/C(C)=C)=C(N=C(C4=O)C)OC)ON7C(CC8=CNC9=C8C=CC=C9/C=C/C(C)=C)=C(N=C(C7=O)C)OC)=C(N=C(C3=O)C)OC)C')







 #   graph = build_graph()
#    print(graph.nodes)
#    print(graph.edges)
#    print(graph)

#    nx.draw(graph)
#    plt.show()


    
        
