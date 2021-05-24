#!/usr/bin/env python


from sys import argv
from typing import (
    Iterable, Callable, List, Optional, Generator, Dict, Union, Tuple, Type,
    DefaultDict,
)

import json
import numpy as np
import re

from dataclasses import dataclass, field, InitVar, asdict
from collections import defaultdict
from copy import deepcopy

TreeDict = Dict[str, Union[str, int, float, List[Optional['TreeDict']]]]


class Placement:
    def __init__(self, name, edge_number, likelihood, like_weight_ratio, distal_length, pendant_length):
        self.name = name
        self.edge_number = int(edge_number)
        self.likelihood = float(likelihood)
        self.like_weight_ratio = float(like_weight_ratio)
        self.distal_length = float(distal_length)
        self.pendant_length = float(pendant_length)

    def __repr__(self):
        return self.name

def parse_jplace_epa(jplace_file):
    with open(jplace_file, 'r') as jplace:
        jplace_data = json.load(jplace)

    tree = jplace_data['tree']

    pqueries = jplace_data['placements']
    placement_names = []
    placement_data = []

    for pquery in pqueries:
        placement_name = pquery['n'][0]
        placement_info = pquery['p'][0]
        placement_names.append(placement_name)
        placement_data.append(placement_info)


    placement_list = []

    for i, name in enumerate(placement_names):
        edge_number, likelihood, like_weight_ratio, distal_length, pendant_length = placement_data[i]
        placement = Placement(name, edge_number, likelihood, like_weight_ratio, distal_length, pendant_length)
        placement_list.append(placement)

    print(tree)
    tree = JPlaceTree.from_jplace(tree)

    return tree, placement_list

@dataclass
class JPlaceTree:
    """[summary]
    Returns:
        [type]: [description]
    """

    name: Optional[str] = None
    edge_number: Optional[str] = None
    length: float = 0.0
    children: Optional[List['Tree']] = field(default_factory=list)

    ID: InitVar[Optional[Union[int, str]]] = None
    depth: InitVar[Optional[int]] = None
    parent: InitVar[Optional['Tree']] = None
    cumulative_length: InitVar[float] = 0.0

    def __post_init__(self, ID, *args, **kwargs):
        """[summary]
        Args:
            ID ([type]): [description]
        """
        self.ID = ID

    @property
    def loc(self) -> 'Tree':
        """Name based index
        Example:
            >>> from picea import Tree
            >>> newick = '(((a,b),(c,d)),e);'
            >>> tree = Tree.from_newick(newick)
            >>> tree.loc['a']
            Tree(name='a', length=0.0, children=[])
        Returns:
            Tree: tree node matching name
        Raises:
            IndexError
        """
        return TreeIndex(
            iterator=self.depth_first,
            eq_func=lambda node, name: node.name == name
        )

    @property
    def iloc(self) -> 'Tree':
        """Index based index
        Example:
            >>> from picea import Tree
            >>> newick = '(((a,b),(c,d)),e);'
            >>> tree = Tree.from_newick(newick)
            >>> tree.iloc[2]
            Tree(name='', length=0.0, children=[Tree(name='a', length=0.0, \
children=[]), Tree(name='b', length=0.0, children=[])])
        Returns:
            Tree: tree node matching index
        """
        return TreeIndex(
            iterator=self.depth_first,
            eq_func=lambda node, index: node.ID == index
        )

    @property
    def root(self) -> 'Tree':
        """Root node of the (sub)tree
        Returns:
            Tree: Root node
        """
        root = self
        while root.parent:
            root = root.parent
        return root

    @property
    def nodes(self) -> List['Tree']:
        """A list of all tree nodes in breadth-first order
        Returns:
            list: A list of all tree nodes
        """
        return list(self.breadth_first())

    @property
    def leaves(self) -> List['Tree']:
        """A list of leaf nodes only
        Returns:
            list: A list of leaf nodes only
        """
        return [n for n in self.nodes if not n.children]

    @property
    def links(self) -> List[Tuple['Tree', 'Tree']]:
        """A list of all (parent, child) combinations
        Returns:
            list: All (parent,child) combinations
        """
        _links = []
        for node in self.nodes:
            if node.children:
                for child in node.children:
                    _links.append((node, child))
        return _links

    @classmethod
    def from_jplace(
            cls,
            string: Optional[str] = None,
            filename: Optional[str] = None
    ) -> 'Tree':
        """Parse a newick formatted string into a Tree object
        Arguments:
            newick_string (string): Newick formatted tree string
        Returns:
            Tree: Tree object
        """
        assert filename or string
        assert not (filename and string)
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()
        tokens = re.split(r'\s*(;|\(|\)|,|:|\{|\})\s*', string)
        print(tokens)
        ID = 0
        tree = cls(ID=ID, length=0.0, cumulative_length=0.0)
        ancestors = list()

        for i, token in enumerate(tokens):
            if token == '(':
                ID += 1
                subtree = cls(ID=ID)
                tree.children = [subtree]
                ancestors.append(tree)
                tree = subtree
            elif token == ',':
                ID += 1
                subtree = cls(ID=ID)
                ancestors[-1].children.append(subtree)
                tree = subtree
            elif token == ')':
                tree = ancestors.pop()
            elif token == '{':
                pass
            elif token == '}':
                pass

            else:
                previous_token = tokens[i - 1]
                if previous_token in ('(', ')', ','):
                    tree.name = token
                elif previous_token == ':':
                    tree.length = float(token)
                elif previous_token == '{':
                    tree.edge_number = int(token)
        tree.depth = 0
        queue = [tree]
        while queue:
            node = queue.pop(0)
            for child in node.children:
                child.parent = node
                child.depth = node.depth + 1
                child.cumulative_length = node.cumulative_length \
                                          + abs(child.length)
            queue += node.children

        return tree

    def to_newick(
        self,
        branch_lengths: bool = True
    ) -> str:
        """Make a Newick formatted string
        Args:
            branch_lengths (bool, optional): Whether to include branch lengths\
             in the Newick string. Defaults to True.
        Returns:
            String: Newick formatted tree string
        """
        if self.name:
            name = str(self.name)
        else:
            name = ''

        if self.children:
            subtree_string = ','.join([
                c.to_newick(branch_lengths=branch_lengths)
                for c in self.children
            ])
            newick = f'({subtree_string}){name}'
        else:
            newick = name

        if branch_lengths and self.ID != 0:
            length = self.length
            if length == 0:
                length = int(0)
            newick += f':{length}'

        if self == self.root:
            newick += ';'

        return newick

    def to_dict(self) -> TreeDict:
        """[summary]
        Returns:
            TreeDict: [description]
        """
        return asdict(self)

    def breadth_first(self) -> Generator['Tree', None, None]:
        """Generator implementing breadth first search starting at root node
        """
        queue = [self]
        while queue:
            node = queue.pop(0)
            queue += node.children
            yield node

    def depth_first(
            self,
            post_order: bool = True
    ) -> Generator['Tree', None, None]:
        """Generator implementing depth first search in either post- or
        pre-order traversel
        Keyword Arguments:
            post_order (bool, optional): Depth first search in post-order
            traversal or not. Defaults to True
        """
        if not post_order:
            yield self
        for child in self.children:
            yield from child.depth_first(post_order=post_order)
        if post_order:
            yield self

    def rename_leaves(
            self,
            rename_func: Callable,
            inplace: bool = True
    ) -> Optional['Tree']:
        """[summary]
        """
        tree = self if inplace else deepcopy(self)
        for leaf in tree.leaves:
            leaf.name = rename_func(leaf.name)

def parse_jplace(jplace_file):
    jplace_data = json.load(jplace_file)

    tree = jplace_data['tree']
    fields = jplace_data['fields']

    placements = jplace_data['placements']

    name_to_placement = {}

    placement_names = placements['n']
    placement_data = placements['p']

    for i, name in enumerate(placement_names):
        associated_data = placement_data[i]
        for j, field in enumerate(fields):
            name_to_placement[name][field] = associated_data[j]

    return tree, name_to_placement

def place_branch(tree, placement):
    for node in tree.breadth_first():
        if node.edge_number == placement.edge_number:
            new_internal_node = Tree(ID='placed_internal_node')
            new_internal_node.name = str(placement.like_weight_ratio)
            new_internal_node.length = placement.distal_length
            new_internal_node.parent = node.parent
            new_internal_node.cumulative_length = new_internal_node.parent.length + new_internal_node.length
            new_internal_node.children = [node]
            new_internal_node.depth = node.depth

            new_leaf_node = Tree(ID='placed_leaf_node')
            new_leaf_node.name = placement.name
            new_leaf_node.length = placement.pendant_length
            new_leaf_node.parent = new_internal_node
            new_leaf_node.cumulative_length = new_leaf_node.parent.length + new_leaf_node.length
            new_internal_node.children.append(new_leaf_node)

            node.parent.children.remove(node)
            node.parent = new_internal_node
            new_internal_node.parent.children.append(new_internal_node)

            node.length -= placement.distal_length

            for downstream in node.breadth_first():
                downstream.depth += 1


if __name__ == "__main__":
    tree, placement_list = parse_jplace_epa(argv[1])
    print(Tree(**tree.to_dict()))
    print(tree.to_newick())
    for placement in placement_list:
        place_branch(tree, placement)
    print(Tree(**tree.to_dict()))
    print(tree.to_newick())






