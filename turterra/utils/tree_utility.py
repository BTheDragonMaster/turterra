from dataclasses import dataclass, field, InitVar, asdict
import typing
from pathlib import Path
from Bio import Phylo
from Bio.Phylo import Newick
import numpy as np
import subprocess
import os
from turterra.utils.sequence_utility import get_sequences_from_fasta

from typing import (
    Iterable,
    Callable,
    List,
    Optional,
    Generator,
    Dict,
    Union,
    Tuple,
    Type,
    DefaultDict,
)

import json
import re
from collections import defaultdict
from copy import deepcopy
from sys import argv
from shutil import copyfile

TreeDict = Dict[str, Union[str, int, float, List[Optional["TreeDict"]]]]


@dataclass
class Tree:
    nodes: list
    edges: list
    col_positions: dict
    row_positions: dict
    tree: Newick.Tree

    @classmethod
    def from_newick(
        cls,
        newick_tree_file: typing.Union[Path, str],
        column_width=80,
        xlen=30,
        ylen=10,
        grabbable=False,
    ):
        """
        Generate positions of nodes from a phylogenetic tree, as well as edges between nodes
        Parameters
        ----------
        newick_tree_file
            tree file in Newick format
        column_width
        xlen
            horizontal tree size
        ylen
            vertical tree size
        grabbable
            if True - nodes can be dragged and moved

        Returns
        -------
        Tree object
        """
        tree = Phylo.read(newick_tree_file, "newick")
        taxa = tree.get_terminals()

        # Some constants for the drawing calculations
        max_label_width = max(len(str(taxon)) for taxon in taxa)
        drawing_width = column_width - max_label_width - 1

        # Create a mapping of each clade to its column position.
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        # Potential drawing overflow due to rounding -- 1 char per tree layer
        fudge_margin = int(np.ceil(np.log2(len(taxa))))
        cols_per_branch_unit = (drawing_width - fudge_margin) / float(
            max(depths.values())
        )

        col_positions = dict(
            (clade, int(blen * cols_per_branch_unit + 1.0))
            for clade, blen in depths.items()
        )

        row_positions = dict((taxon, 2 * idx) for idx, taxon in enumerate(taxa))

        def calc_row(clade):
            for subclade in clade:
                if subclade not in row_positions:
                    calc_row(subclade)
            row_positions[clade] = (
                row_positions[clade.clades[0]] + row_positions[clade.clades[-1]]
            ) // 2

        calc_row(tree.root)
        tree_class = cls([], [], col_positions, row_positions, tree)
        tree_class.add_to_elements(tree.clade, "r", xlen, ylen, grabbable)
        return tree_class

    def add_to_elements(self, clade, clade_id, xlen=30, ylen=10, grabbable=False):
        children = clade.clades

        pos_x = self.col_positions[clade] * xlen
        pos_y = self.row_positions[clade] * ylen

        cy_source = {
            "data": {"id": clade_id},
            "position": {"x": pos_x, "y": pos_y},
            "classes": "nonterminal",
            "grabbable": grabbable,
        }
        self.nodes.append(cy_source)

        if clade.is_terminal():
            cy_source["data"]["name"] = clade.name
            cy_source["classes"] = "terminal"

        for n, child in enumerate(children):
            # The "support" node is on the same column as the parent clade,
            # and on the same row as the child clade. It is used to create the
            # 90 degree angle between the parent and the children.
            # Edge config: parent -> support -> child

            support_id = clade_id + "s" + str(n)
            child_id = clade_id + "c" + str(n)
            pos_y_child = self.row_positions[child] * ylen

            cy_support_node = {
                "data": {"id": support_id},
                "position": {"x": pos_x, "y": pos_y_child},
                "grabbable": grabbable,
                "classes": "support",
            }

            cy_support_edge = {
                "data": {
                    "source": clade_id,
                    "target": support_id,
                    "sourceCladeId": clade_id,
                },
            }

            cy_edge = {
                "data": {
                    "source": support_id,
                    "target": child_id,
                    "length": clade.branch_length,
                    "sourceCladeId": clade_id,
                },
            }

            # if clade.confidence and clade.confidence.value:
            #     cy_source['data']['confidence'] = clade.confidence.value

            self.nodes.append(cy_support_node)
            self.edges.extend([cy_support_edge, cy_edge])

            self.add_to_elements(child, child_id, xlen, ylen, grabbable)


SELECTED_HIGHLIGHT_STYLE = {
            "label": "data(name)",
            "width": 25,
            "height": 25,
            "shape": "star",
            "text-valign": "center",
            "text-halign": "right",
            "background-color": "#FF0000",
        }
TREE_STYLESHEET = [
    {
        "selector": ".nonterminal",
        "style": {
            "label": "data(confidence)",
            "background-opacity": 0,
            "text-halign": "left",
            "text-valign": "top",
            "selectable": "yes",
        },
    },
    {"selector": ".support", "style": {"background-opacity": 0}},
    {
        "selector": "edge",
        "style": {
            "source-endpoint": "inside-to-node",
            "target-endpoint": "inside-to-node",
        },
    },
    {
        "selector": ".terminal",
        "style": {
            "label": "data(name)",
            "width": 10,
            "height": 10,
            "text-valign": "center",
            "text-halign": "right",
            "background-color": "#222222",
            "selectable": "yes",
        },
    },
    {  # change 'classes' of selected nodes to 'selected'
        "selector": ".selected",
        "style": SELECTED_HIGHLIGHT_STYLE
    },
]


class Placement:
    def __init__(
        self,
        name,
        edge_number,
        likelihood,
        like_weight_ratio,
        distal_length,
        pendant_length,
    ):
        self.name = name
        self.edge_number = int(edge_number)
        self.likelihood = float(likelihood)
        self.like_weight_ratio = float(like_weight_ratio)
        self.distal_length = float(distal_length)
        self.pendant_length = float(pendant_length)

    def __repr__(self):
        return self.name


def parse_jplace(jplace_file):
    with open(jplace_file, "r") as jplace:
        jplace_data = json.load(jplace)

    tree = jplace_data["tree"]

    pqueries = jplace_data["placements"]
    placement_names = []
    placement_data = []

    for pquery in pqueries:
        placement_name = pquery["n"][0]
        placement_info = pquery["p"][0]
        placement_names.append(placement_name)
        placement_data.append(placement_info)

    placement_list = []

    for i, name in enumerate(placement_names):
        (
            edge_number,
            likelihood,
            like_weight_ratio,
            distal_length,
            pendant_length,
        ) = placement_data[i]
        placement = Placement(
            name,
            edge_number,
            likelihood,
            like_weight_ratio,
            distal_length,
            pendant_length,
        )
        placement_list.append(placement)

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
    children: Optional[List["Tree"]] = field(default_factory=list)

    ID: InitVar[Optional[Union[int, str]]] = None
    depth: InitVar[Optional[int]] = None
    parent: InitVar[Optional["Tree"]] = None
    cumulative_length: InitVar[float] = 0.0

    def __post_init__(self, ID, *args, **kwargs):
        """[summary]
        Args:
            ID ([type]): [description]
        """
        self.ID = ID

    @property
    def loc(self) -> "Tree":
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
            iterator=self.depth_first, eq_func=lambda node, name: node.name == name
        )

    @property
    def iloc(self) -> "Tree":
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
            iterator=self.depth_first, eq_func=lambda node, index: node.ID == index
        )

    @property
    def root(self) -> "Tree":
        """Root node of the (sub)tree
        Returns:
            Tree: Root node
        """
        root = self
        while root.parent:
            root = root.parent
        return root

    @property
    def nodes(self) -> List["Tree"]:
        """A list of all tree nodes in breadth-first order
        Returns:
            list: A list of all tree nodes
        """
        return list(self.breadth_first())

    @property
    def leaves(self) -> List["Tree"]:
        """A list of leaf nodes only
        Returns:
            list: A list of leaf nodes only
        """
        return [n for n in self.nodes if not n.children]

    @property
    def links(self) -> List[Tuple["Tree", "Tree"]]:
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
        cls, string: Optional[str] = None, filename: Optional[str] = None
    ) -> "Tree":
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
        tokens = re.split(r"\s*(;|\(|\)|,|:|\{|\})\s*", string)
        ID = 0
        tree = cls(ID=ID, length=0.0, cumulative_length=0.0)
        ancestors = list()

        for i, token in enumerate(tokens):
            if token == "(":
                ID += 1
                subtree = cls(ID=ID)
                tree.children = [subtree]
                ancestors.append(tree)
                tree = subtree
            elif token == ",":
                ID += 1
                subtree = cls(ID=ID)
                ancestors[-1].children.append(subtree)
                tree = subtree
            elif token == ")":
                tree = ancestors.pop()
            elif token == "{":
                pass
            elif token == "}":
                pass

            else:
                previous_token = tokens[i - 1]
                if previous_token in ("(", ")", ","):
                    tree.name = token
                elif previous_token == ":":
                    tree.length = float(token)
                elif previous_token == "{":
                    tree.edge_number = int(token)
        tree.depth = 0
        queue = [tree]
        while queue:
            node = queue.pop(0)
            for child in node.children:
                child.parent = node
                child.depth = node.depth + 1
                child.cumulative_length = node.cumulative_length + abs(child.length)
            queue += node.children

        return tree

    def to_newick(self, branch_lengths: bool = True) -> str:
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
            name = ""

        if self.children:
            subtree_string = ",".join(
                [c.to_newick(branch_lengths=branch_lengths) for c in self.children]
            )
            newick = f"({subtree_string}){name}"
        else:
            newick = name

        if branch_lengths and self.ID != 0:
            length = self.length
            if length == 0:
                length = int(0)
            newick += f":{length}"

        if self == self.root:
            newick += ";"

        return newick

    def to_dict(self) -> TreeDict:
        """[summary]
        Returns:
            TreeDict: [description]
        """
        return asdict(self)

    def breadth_first(self) -> Generator["Tree", None, None]:
        """Generator implementing breadth first search starting at root node
        """
        queue = [self]
        while queue:
            node = queue.pop(0)
            queue += node.children
            yield node

    def depth_first(self, post_order: bool = True) -> Generator["Tree", None, None]:
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
        self, rename_func: Callable, inplace: bool = True
    ) -> Optional["Tree"]:
        """[summary]
        """
        tree = self if inplace else deepcopy(self)
        for leaf in tree.leaves:
            leaf.name = rename_func(leaf.name)


def place_branch(tree, placement):
    for node in tree.breadth_first():
        if node.edge_number == placement.edge_number:
            new_internal_node = JPlaceTree(ID="placed_internal_node")
            new_internal_node.name = str(placement.like_weight_ratio)
            new_internal_node.length = placement.distal_length
            new_internal_node.parent = node.parent
            new_internal_node.cumulative_length = (
                new_internal_node.parent.length + new_internal_node.length
            )
            new_internal_node.children = [node]
            new_internal_node.depth = node.depth

            new_leaf_node = JPlaceTree(ID="placed_leaf_node")
            new_leaf_node.name = placement.name
            new_leaf_node.length = placement.pendant_length
            new_leaf_node.parent = new_internal_node
            new_leaf_node.cumulative_length = (
                new_leaf_node.parent.length + new_leaf_node.length
            )
            new_internal_node.children.append(new_leaf_node)

            node.parent.children.remove(node)
            node.parent = new_internal_node
            new_internal_node.parent.children.append(new_internal_node)

            node.length -= placement.distal_length

            for downstream in node.breadth_first():
                downstream.depth += 1


def make_directory(parent, dir_name):
    dir = os.path.join(parent, dir_name)
    if not os.path.exists(dir):
        os.mkdir(dir)

    return dir


def add_branch(placement_nr, ref_msa, ref_nwk, query_faa, out_dir, tmp_tree):
    command = [
        "epa-ng",
        "--ref-msa",
        ref_msa,
        "--tree",
        ref_nwk,
        "--query",
        query_faa,
        "--model",
        "JTT",
        "--outdir",
        out_dir,
        "--redo",
    ]
    subprocess.check_call(command, shell=False)

    tree, placements = parse_jplace(os.path.join(out_dir, "epa_result.jplace"))
    placement = placements[0]

    place_branch(tree, placement)
    updated_tree = os.path.join(tmp_tree, f"updated_tree_{placement_nr}.nwk")

    with open(updated_tree, "w") as out_file:
        out_file.write(tree.to_newick())

    return updated_tree

def add_branches(ref_msa, query_fasta, query_folder, nwk_tree, epa_output_dir, tmp_tree, data_folder):

    query_ids = list(get_sequences_from_fasta(query_fasta).keys())

    updated_tree = nwk_tree
    ref_to_seq = get_sequences_from_fasta(ref_msa)

    for i, query_id in enumerate(query_ids):
        aligned_seq = ref_to_seq[query_id]
        query_out_dir = os.path.join(query_folder, f"{query_id}.faa")
        with open(query_out_dir, "w") as out_file:
            out_file.write(f">{query_id}\n{aligned_seq}")

        epa_dir = make_directory(epa_output_dir, query_id)

        updated_tree = add_branch(
            i, ref_msa, updated_tree, query_out_dir, epa_dir, tmp_tree
        )


    copyfile(updated_tree, os.path.join(data_folder, 'updated_tree.nwk'))
