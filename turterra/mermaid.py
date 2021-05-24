import json
import re
from collections import defaultdict
from dataclasses import dataclass
import typing as ty
from prompt_toolkit import prompt
from prompt_toolkit.shortcuts import confirm
from prompt_toolkit.completion import WordCompleter
import typer
from pathlib import Path

RE_FUNCTION = "@app\.callback\((?P<parameters>[\s\S]*?)(?=\)\n    def )"
RE_OUTPUT = 'Output\("(?P<component>[a-z-_A-Z]*)", "(?P<property>["a-z-_A-Z]*)"\)'
RE_INPUT = '(?:Input|Trigger)\("(?P<component>[a-z-_A-Z]*)", "(?P<property>["a-z-_A-Z]*)"\)'
RE_STATE = 'State\("(?P<component>[a-z-_A-Z]*)", "(?P<property>["a-z-_A-Z]*)"\)'
PALETTE = ('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933', '#CC6677', '#882255', '#AA4499', '#DDDDDD')


@dataclass(eq=True, frozen=True)
class Node:
    name: str
    prop: str

    @property
    def parts(self):
        return self.name.split("-")

    @property
    def display_name(self):
        return " ".join(f"{n[0].upper()}{n[1:]}" for n in self.parts[1:-1])

    @property
    def panel_name(self):
        return self.parts[0]

    @property
    def node_type(self):
        return self.parts[-1]

    @property
    def node_name(self):
        return self.name.replace("style", "style1").replace("graph", "graph1").replace("upload", "upload1")

    @property
    def prop_name(self):
        return self.prop.replace("style", "style1").replace("graph", "graph1").replace("upload", "upload1")

    @property
    def node_prop_name(self):
        return f"{self.node_name}:{self.prop_name}"


@dataclass
class Edge:
    input: Node
    states: ty.List[Node]
    output: Node

    def format_edge(self, arrow_type, property_as_edge: bool, state_on_edge: bool):
        input_node = self.input.node_name
        output_node = self.output.node_name
        edge = ""
        if property_as_edge:
            edge = f"{self.input.prop}<br/>"
        if state_on_edge:
            edge += "States: " + "<br/>".join(f"{s.display_name}:{s.prop_name}" for s in self.states)
        if len(edge):
            return f"{input_node}{arrow_type}|{edge}|{output_node}"
        else:
            return f"{input_node}{arrow_type}{output_node}"


def get_nodes_and_edges(callbacks):
    nodes = set()
    edges = []
    for group in re.finditer(RE_FUNCTION, callbacks):
        parameters = group.group("parameters")
        outputs = re.findall(RE_OUTPUT, parameters)
        inputs = re.findall(RE_INPUT, parameters)
        states = re.findall(RE_STATE, parameters)
        for component, prop in inputs + outputs + states:
            nodes.add(Node(component, prop))
        states = [Node(sc, sp) for sc, sp in states]
        for ic, ip in inputs:
            for oc, op in outputs:
                edges.append(Edge(Node(ic, ip),
                                  states,
                                  Node(oc, op)))
    return nodes, edges


def get_contrasting_text_color(hex_str):
    hex_str = hex_str.lower()[1:]
    (r, g, b) = (hex_str[:2], hex_str[2:4], hex_str[4:])
    return '#000' if 1 - (int(r, 16) * 0.299 + int(g, 16) * 0.587 + int(b, 16) * 0.114) / 255 < 0.5 else '#fff'


def get_node_classes(node_types, node_classes: dict):
    for i, node_type in enumerate(node_types):
        if node_type in node_classes:
            continue
        node_shape = prompt(f"Node shape for {node_type}: ",
                            completer=WordCompleter(["Round", "Stadium", "Subroutine",
                                                     "Cylinder", "Circle", "Asymmetric",
                                                     "Rhombus", "Hexagon", "Parallelogram",
                                                     "Parallelogram alt", "Trapezoid", "Trapezoid alt"]), )
        node_color = PALETTE[i % len(PALETTE)]
        font_color = get_contrasting_text_color(node_color)
        classdef = f"classDef {node_type} fill:{node_color},color:{font_color}"
        bracket_open, bracket_close = get_brackets(node_shape)
        node_classes[node_type] = (classdef, bracket_open, bracket_close)
    return node_classes


def get_brackets(node_shape):
    if node_shape == "Round":
        return "(", ")"
    elif node_shape == "Stadium":
        return "([", "])"
    elif node_shape == "Subroutine":
        return "[[", "]]"
    elif node_shape == "Cylinder":
        return "[(", ")]"
    elif node_shape == "Circle":
        return "((", "))"
    elif node_shape == "Asymmetric":
        return ">", "]"
    elif node_shape == "Rhombus":
        return "{", "}"
    elif node_shape == "Hexagon":
        return "{{", "}}"
    elif node_shape == "Parallelogram":
        return "[/", "/]"
    elif node_shape == "Parallelogram alt":
        return "[\\", "\\]"
    elif node_shape == "Trapezoid":
        return "[\\", "/]"
    elif node_shape == "Trapezoid alt":
        return "[/", "\\]"
    else:
        raise ValueError(f"Unknown node shape: {node_shape}")


def main(filename: Path,
         output_file: Path,
         node_file: ty.Optional[Path] = typer.Argument(None)):
    if output_file.exists():
        raise typer.BadParameter(
            f"Output file {output_file} already exists, cowardly refusing to overwrite. Please delete it and try again"
        )
    with open(filename) as f:
        nodes, edges = get_nodes_and_edges(f.read())
    if node_file:
        with open(node_file) as f:
            node_classes = json.load(f)
    else:
        node_classes = {}
    node_types = set(n.node_type for n in nodes)
    node_classes = get_node_classes(node_types, node_classes)
    node_panels = defaultdict(list)
    nodes = list(nodes)
    for i, node in enumerate(nodes):
        node_panels[node.panel_name].append(i)
    arrow_type = prompt("Choose arrow type: ", completer=WordCompleter(["-->", "==>", "-.->"]))
    states_on_edge = confirm("Show States on edge?")

    with open(output_file, "w") as f:
        f.write("```mermaid\n")
        f.write("graph LR\n")
        for node_type in node_classes:
            classdef, _, _ = node_classes[node_type]
            f.write(f"{classdef}\n")
        f.write(f"style legend fill:#FFF,color:#000,stroke-width:0px;\n")
        f.write("subgraph legend [Legend]\n")
        for node_type in node_classes:
            _, bracket_open, bracket_close = node_classes[node_type]
            f.write(
                f"legend-{node_type}{bracket_open}{node_type[0].upper() + node_type[1:]}{bracket_close}:::{node_type}\n")
        f.write("end\n\n")
        for p in node_panels:
            f.write(f"style {p} fill:#FFF,color:#000,stroke:#000,stroke-width:2px;\n")
            f.write(f"subgraph {p} [{p[0].upper() + p[1:]} Panel]\n")
            for index in node_panels[p]:
                node = nodes[index]
                classdef, bracket_open, bracket_close = node_classes[node.node_type]
                f.write(f"{node.node_name}{bracket_open}{node.display_name}{bracket_close}:::{node.node_type};\n")
            f.write("end\n\n")
        for edge in edges:
            f.write(f"{edge.format_edge(arrow_type, True, states_on_edge)}\n")
        f.write("```\n")
    if confirm("Save node classes to file?"):
        output_node_file = prompt("Filename: ", default="node_classes.json")
        with open(output_node_file, "w") as f:
            json.dump(node_classes, f)


if __name__ == "__main__":
    typer.run(main)
