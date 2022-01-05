import dash
import dash_bio
import dash_core_components as dcc
import dash_html_components as html
from dash_extensions.enrich import Trigger, ServersideOutput, Output, Input, State
import os
from zipfile import ZipFile, ZIP_DEFLATED
import zlib
import io
from dash_extensions.snippets import send_file
from turterra.app_layout import BUTTON_STYLE
from turterra.turterra import TurterraData
from turterra.utils import sequence_utility, structure_utility
import base64
from caretta.helper import group_indices
from turterra.utils.tree_utility import TREE_STYLESHEET, SELECTED_HIGHLIGHT_STYLE


def group_into_ranges(coverage):
    indices = []
    for stretch in coverage:
        start = stretch["start"]
        end = stretch["end"]
        for index in range(start, end + 1):
            indices.append(index)

    indices = list(set(indices))
    indices.sort()

    ranges = []
    current_range = [indices[0]]
    for index in indices[1:]:
        if index == current_range[-1] + 1:
            current_range.append(index)
        else:
            ranges.append(current_range)
            current_range = [index]

    ranges.append(current_range)
    new_coverage = []

    for stretch in ranges:
        new_coverage.append(
            {"start": stretch[0], "end": stretch[-1], "bgcolor": "lightblue"}
        )

    return new_coverage


def register_callbacks(app):
    @app.callback(
        [
            ServersideOutput("input-turterra-data", "data"),
            Output("phylogeny-tree-viewer", "elements"),
        ],
        [
            Trigger("input-load-button", "n_clicks"),
            Trigger("upload-add-to-analysis-button", "n_clicks"),
        ],
        [State("input-turterra-data", "data"),],
    )
    def load_and_add_turterra_data(turterra_data):
        ctx = dash.callback_context
        trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if trigger_id == "input-load-button":
            turterra_data = TurterraData.from_folder("test_data")
        if trigger_id == "upload-add-to-analysis-button":
            turterra_data.add_data_points("uploaded_data", "test_data")
        if turterra_data:
            tree_elements = {
                "nodes": turterra_data.tree.nodes,
                "edges": turterra_data.tree.edges,
            }
            return [turterra_data, tree_elements]
        return [None, []]

    # Filter accessions
    @app.callback(
        [
            # species list
            Output("filter-species-dropdown", "options"),
            # compounds list
            Output("filter-compound-dropdown", "options"),
            # accessions list
            Output("filter-accession-dropdown", "options"),
        ],
        [
            Input("input-turterra-data", "data"),
            Trigger("filter-filter-button", "n_clicks"),
        ],
        [
            # species to filter by
            State("filter-species-dropdown", "value"),
            # compounds to filter by
            State("filter-compound-dropdown", "value"),
            # turterra data
            # State("input-turterra-data", "data"),
        ],
        prevent_initial_call=True,
    )
    def filter_accessions(
        turterra_data, species_filter_dropdown, compound_filter_dropdown
    ):
        ctx = dash.callback_context
        trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if turterra_data:
            species = [
                {"label": x, "value": x} for x in set(turterra_data.species.values())
            ]
            compounds = [{"label": c, "value": c} for c in turterra_data.compounds_dict]
            accessions = [{"label": x, "value": x} for x in turterra_data.accessions]
            if trigger_id == "filter-filter-button":
                print("Filtering...")
                if species_filter_dropdown:
                    compounds = set(
                        c
                        for x in turterra_data.accessions
                        for c in turterra_data.compounds_mapping[x]
                        if turterra_data.species[x] in species_filter_dropdown
                    )
                    compounds = [{"label": c, "value": c} for c in compounds]
                    accessions = [
                        x
                        for x in accessions
                        if turterra_data.species[x["label"]] in species_filter_dropdown
                    ]
                if compound_filter_dropdown:
                    accessions = [
                        x
                        for x in accessions
                        if any(
                            c in compound_filter_dropdown
                            for c in turterra_data.compounds_mapping[x["label"]]
                        )
                    ]
                print(
                    f"After filtering: {len(species)} species, {len(compounds)} compounds, {len(accessions)} accessions"
                )
            return species, compounds, accessions
        else:
            return [], [], []

    # select-all
    @app.callback(
        Output("filter-accession-dropdown", "value"),
        [Trigger("filter-select-all-button", "n_clicks")],
        [State("filter-accession-dropdown", "options")],
    )
    def select_all_accessions(accessions_dropdown):
        return [a["value"] for a in accessions_dropdown]

    # highlight selected accessions on tree
    @app.callback(
        Output("phylogeny-tree-viewer", "stylesheet"),
        [Input("filter-accession-dropdown", "value")],
        [
            # tree nodes
            State("phylogeny-tree-viewer", "elements"),
        ],
    )
    def highlight_nodes_in_tree(selected_accessions, elements):
        if not selected_accessions:
            return TREE_STYLESHEET
        selected_accessions = set(selected_accessions)
        selected_nodes = []
        for element in elements["nodes"]:
            if (
                "name" in element["data"]
                and element["data"]["name"] in selected_accessions
            ):
                selected_nodes.append(element["data"]["id"])
        stylesheet = [
            {"selector": f'node[id = "{node}"]', "style": SELECTED_HIGHLIGHT_STYLE}
            for node in selected_nodes
        ]
        return TREE_STYLESHEET + stylesheet

    # add entries from tree
    @app.callback(
        [
            Output("filter-accession-dropdown", "value"),
            Output("filter-accession-dropdown", "options"),
        ],
        [Input("phylogeny-tree-viewer", "tapNode")],
        [
            State("filter-accession-dropdown", "value"),
            State("filter-accession-dropdown", "options"),
            State("input-turterra-data", "data"),
        ],
    )
    def add_accessions_from_tree(
        node_data, selected_accessions, possible_accessions, turterra_data
    ):
        if node_data and turterra_data:
            accessions = set(x["value"] for x in possible_accessions)
            tree_accessions = []
            if not selected_accessions:
                selected_accessions = []
            if "name" in node_data["data"]:  # only terminal nodes have a name
                tree_accessions.append(node_data["data"]["name"].split("\t")[0])
            # else:
            #     mapping = defaultdict(list)
            #     for node in tree_elements["nodes"]:
            #         mapping[node["data"]["sourceCladeId"]].append(node["data"]["name"])
            #     for edge in node_data["edgesData"]:
            #         tree_accessions += [node.split("\t")[0] for node in mapping[edge["sourceCladeId"]]]
            for accession in tree_accessions:
                if accession not in selected_accessions:
                    selected_accessions.append(accession)
                if accession not in accessions:
                    possible_accessions.append(dict(label=accession, value=accession))
        return selected_accessions, possible_accessions

    # update accession dropdown list in all panels
    @app.callback(
        [
            # selected sequence key options
            Output("sequence-accession-select-dropdown", "options"),
            # selected structure key options
            Output("structure-accession-select-dropdown", "options"),
            # selected compound key options
            Output("compound-accession-select-dropdown", "options"),
            # color scheme options
            Output("structure-color-scheme-select-dropdown", "options"),
            # color scheme value
            Output("structure-color-scheme-select-dropdown", "value"),
        ],
        [
            Input("filter-accession-dropdown", "value"),
            # Whether the show button has been clicked
            Trigger("filter-show-button", "n_clicks"),
        ],
        [State("input-turterra-data", "data")],
    )
    def update_accessions(accession_filter_dropdown, turterra_data):
        if not accession_filter_dropdown:
            return [], [], [], [], ""
        sequence_key_select_dropdown_options = [
            dict(label=a, value=a) for a in accession_filter_dropdown
        ]
        compound_key_select_dropdown_options = sequence_key_select_dropdown_options  # TODO: maybe some accessions don't have compounds
        structure_key_select_dropdown_options = [
            {"label": a, "value": a}
            for a in accession_filter_dropdown
            if a in turterra_data.structures
        ]
        if len(accession_filter_dropdown) > 1:
            color_schemes = [
                "sequence conservation",
                "structure conservation",
                "atom",
                "residue",
                "residue_type",
            ]
            color_scheme_select_dropdown_value = "structure conservation"
        elif len(accession_filter_dropdown) == 1:
            color_schemes = ["atom", "residue", "residue_type"]
            color_scheme_select_dropdown_value = "residue"
        else:
            color_schemes = []
            color_scheme_select_dropdown_value = ""
        color_scheme_select_dropdown_options = [
            {"label": x, "value": x} for x in color_schemes
        ]
        return (
            sequence_key_select_dropdown_options,
            structure_key_select_dropdown_options,
            compound_key_select_dropdown_options,
            color_scheme_select_dropdown_options,
            color_scheme_select_dropdown_value,
        )

    # load/update sequence panel
    @app.callback(
        # Outputs
        # sequence panel
        Output("sequence-sequence-viewer", "children"),
        # Inputs
        [
            # Accessions
            Input("sequence-accession-select-dropdown", "options"),
            # Selected accession
            Input("sequence-accession-select-dropdown", "value"),
            # Type of view
            Input("sequence-alignment-radio-select-button", "value"),
        ],
        # States
        [
            # turterra data
            State("input-turterra-data", "data"),
        ],
        prevent_initial_call=True,
    )
    def update_sequence_panel(
        accessions, selected_accession, alignment_select_radio_buttons, turterra_data,
    ):
        if selected_accession:
            accessions = [a["label"] for a in accessions]
            sequence_alignment = turterra_data.get_sub_alignment(accessions, "sequence")
            structure_alignment = turterra_data.get_sub_alignment(
                accessions, "structure"
            )
            if alignment_select_radio_buttons == "Single Sequence":
                return [
                    dash_bio.SequenceViewer(
                        id="sequence-single-sequence-viewer",
                        sequence=sequence_alignment[selected_accession]
                        .replace("-", "")
                        .upper(),
                        title=selected_accession,
                        charsPerLine=50,
                        coverage=[],
                    ),
                    get_sequence_viewer_information(selected_accession, turterra_data),
                    html.Button(
                        "Reset selection",
                        id="structure-reset-selection-button",
                        className="pure-button pure-button-primary",
                        style=BUTTON_STYLE,
                    ),
                ]
            if len(accessions) > 1:
                if alignment_select_radio_buttons == "Sequence Alignment":
                    alignment = sequence_utility.alignment_to_fasta(sequence_alignment)
                else:
                    alignment = sequence_utility.alignment_to_fasta(structure_alignment)
                return dash_bio.AlignmentChart(
                    id="sequence-sequence-alignment-viewer",
                    data=alignment,
                    showgap=False,
                    height=500,
                    showconservation=True,
                    showconsensus=False,
                )
            elif len(accessions) == 1:
                accession = accessions[0]
                return [
                    dash_bio.SequenceViewer(
                        id="sequence-single-sequence-viewer",
                        sequence=sequence_alignment[accession].replace("-", "").upper(),
                        title=accession,
                        charsPerLine=100,
                        coverage=[],
                    ),
                    get_sequence_viewer_information(accession, turterra_data),
                    html.Button(
                        "Reset selection",
                        id="structure-reset-selection-button",
                        className="pure-button pure-button-primary",
                        style=BUTTON_STYLE,
                    ),
                ]
        return ""

    # make structure information
    @app.callback(
        # Outputs
        # structure data
        ServersideOutput("structure-structure-data", "data"),
        # Inputs
        [
            # selected accession for displaying structure
            Input("structure-accession-select-dropdown", "value"),
            # selected visualization scheme
            Input("structure-visualisation-select-dropdown", "value"),
            # selected color scheme
            Input("structure-color-scheme-select-dropdown", "value"),
            # Whether the show button has been clicked
            Trigger("filter-show-button", "n_clicks"),
            # whether the refresh button has been clicked
            Trigger("structure-refresh-button", "n_clicks"),
        ],
        # States
        [
            # accessions
            State("structure-accession-select-dropdown", "options"),
            # turterra data
            State("input-turterra-data", "data"),
        ],
    )
    def get_structure(
        selected_accession,
        selected_visualisation,
        selected_color_scheme,
        accessions,
        turterra_data,
    ):
        if accessions and selected_accession:
            accessions = [a["label"] for a in accessions]
            sequence_alignment = turterra_data.get_sub_alignment(accessions, "sequence")
            structure_alignment = turterra_data.get_sub_alignment(
                accessions, "structure"
            )
            structure_data = turterra_data.structures[selected_accession].plot_data
            styles_data = turterra_data.structures[selected_accession].get_style_data(
                selected_visualisation,
                selected_color_scheme,
                sequence_alignment,
                structure_alignment,
            )
            return {"styles": styles_data, "structure": structure_data}
        return None

    # make 3d structure viewer
    @app.callback(
        # panel with 3D viewer
        [
            Output("structure-structure-container", "children"),
            Output("structure-model-information", "children"),
        ],
        [
            # current selected structure data
            Input("structure-structure-data", "data")
        ],
        [
            State("structure-accession-select-dropdown", "value"),
            # turterra data
            State("input-turterra-data", "data"),
        ],
    )
    def show_structure(structure_data, selected_accession, turterra_data):
        if structure_data and turterra_data:
            return (
                dash_bio.Molecule3dViewer(
                    id="structure-structure-viewer",
                    modelData=structure_data["structure"],
                    styles=structure_data["styles"],
                    selectionType="residue",
                ),
                get_structure_information(selected_accession, turterra_data),
            )
        else:
            return "", ""

    # sync selected accession across panels
    @app.callback(
        [
            Output("sequence-accession-select-dropdown", "value"),
            Output("structure-accession-select-dropdown", "value"),
            Output("compound-accession-select-dropdown", "value"),
        ],
        [
            Input("sequence-accession-select-dropdown", "value"),
            Input("structure-accession-select-dropdown", "value"),
            Input("compound-accession-select-dropdown", "value"),
        ],
    )
    def sync_selected_accession(
        sequence_accession, structure_accession, compound_accession
    ):
        ctx = dash.callback_context
        trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if trigger_id == "sequence-accession-select-dropdown":
            return sequence_accession, sequence_accession, sequence_accession
        elif trigger_id == "structure-accession-select-dropdown":
            return structure_accession, structure_accession, structure_accession
        elif trigger_id == "compound-accession-select-dropdown":
            return compound_accession, compound_accession, compound_accession
        else:
            return sequence_accession, structure_accession, compound_accession

    # reset selected residues in sequence and structure panels
    @app.callback(
        # Outputs
        [
            Output("sequence-single-sequence-viewer", "coverage"),
            Output("sequence-single-sequence-viewer", "mouseSelection"),
            Output("structure-selected-residue-indices-data", "data"),
            Output("structure-refresh-button", "n_clicks"),
        ],
        # Inputs
        Trigger("structure-reset-selection-button", "n_clicks"),
        # States
        State("structure-refresh-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_selected_residues(n_clicks):
        coverage = []
        mouse_selection = None
        if n_clicks != None:
            n_clicks += 1
        else:
            n_clicks = 1
        return coverage, mouse_selection, [], n_clicks

    # show selected residues in structure panel
    @app.callback(
        # Outputs
        Output("structure-selected-residue-indices-data", "data"),
        # Inputs
        [
            # selected residues in structure viewer
            Input("sequence-single-sequence-viewer", "mouseSelection"),
            Input("structure-structure-viewer", "selectedAtomIds"),
        ],
        # States
        [
            State("structure-structure-data", "data"),
            State("structure-accession-select-dropdown", "value"),
            State("input-turterra-data", "data"),
            State("structure-refresh-button", "n_clicks"),
            State("structure-selected-residue-indices-data", "data"),
        ],
    )
    def show_selected_residues_structure(
        single_sequence_mouse_selection,
        atom_ids_current,
        structure_data,
        selected_structure_key,
        turterra_data,
        structure_loaded,
        atom_ids,
    ):
        if not structure_loaded:
            return None
        selected_residue_ids = []
        if not atom_ids:
            atom_ids = []
        if atom_ids_current:
            atom_ids += atom_ids_current
        if single_sequence_mouse_selection is not None:
            if single_sequence_mouse_selection["selection"]:
                selected_residue_ids += list(
                    range(
                        single_sequence_mouse_selection["start"],
                        single_sequence_mouse_selection["end"] + 1,
                    )
                )

        atom_ids = set(atom_ids)

        if not selected_residue_ids:
            return list(atom_ids)
        start_index, end_index = structure_utility.map_structure_to_sequence(
            structure_data["structure"], turterra_data.sequences[selected_structure_key]
        )
        # map residues to atom indices
        indices = group_indices(
            [
                structure_data["structure"]["atoms"][i]["residue_index"]
                for i in range(len(structure_data["structure"]["atoms"]))
            ]
        )
        if selected_residue_ids:
            for residue_id in selected_residue_ids:
                if start_index <= residue_id <= end_index:
                    # add all atom ids of selected residue
                    atom_ids = atom_ids.union(indices[residue_id - start_index - 1])
        return list(atom_ids)

    # highlight selected residues in sequence
    @app.callback(
        # Outputs
        [
            Output("structure-structure-data", "data"),
            Output("structure-residue-information", "children"),
            Output("sequence-single-sequence-viewer", "coverage"),
        ],
        # Inputs
        [
            # change structure viewer's selection
            Input("structure-selected-residue-indices-data", "data"),
            Input('structure-residue-color-picker', 'value'),
            Input('structure-residue-visualisation-select-dropdown', 'value'),
        ],
        # States
        [
            State("sequence-single-sequence-viewer", "mouseSelection"),
            State("structure-structure-data", "data"),
            State("structure-accession-select-dropdown", "value"),
            State("input-turterra-data", "data"),
            State("sequence-single-sequence-viewer", "coverage"),
        ],
    )
    def show_selected_residues_structure_information(
        selected_atom_ids,
        residue_color,
        residue_visualization_type,
        single_sequence_mouse_selection,
        structure_data,
        selected_structure_key,
        turterra_data,
        coverage,
    ):

        information = "No atom has been selected. Click somewhere on the molecular structure to select an atom."

        if (
            not selected_atom_ids and not single_sequence_mouse_selection
        ) or not structure_data:
            return structure_data, information, coverage
        start_index, _ = structure_utility.map_structure_to_sequence(
            structure_data["structure"], turterra_data.sequences[selected_structure_key]
        )
        residues = {}
        if selected_atom_ids:
            for atom_id in selected_atom_ids:
                residues[
                    structure_data["structure"]["atoms"][atom_id]["residue_index"]
                ] = structure_data["structure"]["atoms"][atom_id]
                structure_data["styles"][atom_id]["color"] = residue_color["hex"]
                structure_data["styles"][atom_id]["visualization_type"] = residue_visualization_type
        residue_names = []
        for residue in sorted(residues.keys()):
            sequence_id = residue + start_index
            coverage.append(
                {"start": sequence_id - 1, "end": sequence_id, "bgcolor": "lightblue"}
            )
            residue_names.append(residues[residue]["residue_name"])

        if single_sequence_mouse_selection is not None:
            coverage.append(
                {
                    "start": single_sequence_mouse_selection["start"],
                    "end": single_sequence_mouse_selection["end"],
                    "bgcolor": "lightblue",
                }
            )
        if len(coverage):
            coverage = group_into_ranges(coverage)

        if residue_names:
            information = get_selected_residues_information(residue_names)

        return structure_data, information, coverage

    # select compound to view (if multiple exist)
    @app.callback(
        Output("compound-radio-select-button-container", "children"),
        [Input("compound-accession-select-dropdown", "value")],
        [
            # turterra data
            State("input-turterra-data", "data"),
        ],
    )
    def compound_selector(selected_accession, turterra_data):
        if selected_accession:
            return get_compound_radio_buttons(selected_accession, turterra_data)
        else:
            return []

    # update compound viewer
    @app.callback(
        Output("compound-compound-container", "children"),
        [Input("compound-radio-select-button", "value")],
        [State("input-turterra-data", "data")],
    )
    def update_compound(selected_compound, turterra_data):
        if turterra_data:
            return [
                html.Label(f"Currently viewing {selected_compound}"),
                # TODO: make this zoomable
                dash_bio.Molecule2dViewer(
                    id="compound-compound-viewer",
                    modelData=turterra_data.compounds_dict[selected_compound].plot_data,
                    height=300,
                    width=300,
                ),
            ]
        else:
            return ""

    # Export FASTA file for selected accession
    @app.callback(
        # Outputs
        # Link to download file
        Output("sequence-single-download", "data"),
        # Inputs
        [
            # Whether the FASTA download button has been clicked
            Trigger("sequence-single-download-button", "n_clicks")
        ],
        # States
        [
            # Which accession selected
            State("sequence-accession-select-dropdown", "value"),
            State("input-turterra-data", "data"),
        ],
    )
    def download_single_sequence(
        accession, turterra_data
    ):
        if accession and turterra_data:
            return dict(
                content=sequence_utility.format_as_fasta({accession: turterra_data.sequences[accession]}),
                filename=f"{accession}.fasta",
            )
        else:
            return None

    # Export FASTA file for selected accessions
    @app.callback(
        # Outputs
        # Link to download file
        Output("sequence-multiple-download", "data"),
        # Inputs
        [
            # Whether the FASTA download button has been clicked
            Trigger("sequence-multiple-download-button", "n_clicks")
        ],
        # States
        [
            State("sequence-accession-select-dropdown", "options"),
            State("input-turterra-data", "data"),
        ],
    )
    def download_multiple_sequences(
            accessions, turterra_data
    ):
        if accessions and turterra_data:
            accessions = [a["label"] for a in accessions]
            sequences = {k: turterra_data.sequences[k] for k in accessions}
            return dict(
                content=sequence_utility.format_as_fasta(sequences),
                filename="sequences.fasta",
            )
        else:
            return None

    # Export FASTA alignment file
    @app.callback(
        # Outputs
        # Link to download file
        Output("sequence-alignment-download", "data"),
        # Inputs
        [
            # Whether the FASTA download button has been clicked
            Trigger("sequence-alignment-download-button", "n_clicks")
        ],
        # States
        [
            # Which alignment selected
            State("sequence-alignment-radio-select-button", "value"),
            # Filtered accessions
            State("sequence-accession-select-dropdown", "options"),
            State("input-turterra-data", "data"),
        ],
    )
    def download_alignment(
            alignment_select_value, accessions, turterra_data
    ):
        if accessions and alignment_select_value and turterra_data:
            accessions = [a["label"] for a in accessions]
            if (
                    alignment_select_value == "Single Sequence"
                    or alignment_select_value == "Sequence Alignment"
            ):
                alignment = turterra_data.get_sub_alignment(accessions, "sequence")
            else:
                alignment = turterra_data.get_sub_alignment(accessions, "structure")
            if not alignment:
                return ""
            return dict(
                content=sequence_utility.format_as_fasta(alignment),
                filename="alignment.fasta",
            )
        else:
            return None

    # Export PDB file for selected accession
    @app.callback(
        # Outputs
        # Link to download file
        Output("structure-single-download", "data"),
        # Inputs
        [
            # Whether the PDB download button has been clicked
            Trigger("structure-single-download-button", "n_clicks")
        ],
        # States
        [
            # Which accession selected
            State("structure-accession-select-dropdown", "value"),
            State("input-turterra-data", "data"),
        ],
    )
    def download_single_structure(
        accession, turterra_data
    ):
        if accession and turterra_data:
            with open(turterra_data.structures[accession].filename) as f:
                return dict(
                    content=f.read(),
                    filename=f"{accession}.pdb",
                )
        else:
            return None

    # Export PDB files (zipped) for selected accessions
    @app.callback(
        # Outputs
        # Link to download file
        Output("structure-multiple-download", "data"),
        # Inputs
        [
            # Whether the PDB download button has been clicked
            Trigger("structure-multiple-download-button", "n_clicks")
        ],
        # States
        [
            State("structure-accession-select-dropdown", "options"),
            State("input-turterra-data", "data"),
        ],
    )
    def download_multiple_structures(
            accessions, turterra_data
    ):
        if accessions and turterra_data:
            accessions = [a["label"] for a in accessions]
            output_filename = "pdbs.zip"
            pdb_zip_file = ZipFile(output_filename, mode="w", compression=ZIP_DEFLATED)
            for accession in accessions:
                with open(turterra_data.structures[accession].filename) as f:
                    pdb_zip_file.writestr(f"pdbs/{accession}.pdb", f.read())
            pdb_zip_file.close()
            return send_file("pdbs.zip")
        else:
            return None

    # Function to export tree alignment file
    @app.callback(
        # Outputs
        # Link to download file
        Output("phylogeny-newick-download", "data"),
        # Inputs
        [
            # Whether the tree download button has been clicked
            Trigger("phylogeny-newick-download-button", "n_clicks")
        ],
        # States
        # [
        #
        # ],
    )
    def download_tree():
        return None
        # TODO: @BTheDragonMaster
        # return dict(content=tree_newick_string,
        #             filename="tree.newick")

    # Function to upload new sequences
    @app.callback(
        [
            Output("upload-sequences-list", "children"),
            Output("upload-add-to-analysis-button", "style"),
        ],
        [
            Input("upload-sequences-upload", "filename"),
            Input("upload-sequences-upload", "contents"),
        ],
    )
    def upload_sequences(file_name, content):
        if file_name is not None and content is not None:
            content_type, content_string = content.split(",")
            decoded = base64.b64decode(content_string)
            fasta_dir = os.path.join("uploaded_data", "uploaded_sequences.fasta")
            with open(fasta_dir, "wb") as fasta:
                fasta.write(decoded)

            fasta_ids = list(
                sequence_utility.get_sequences_from_fasta(fasta_dir).keys()
            )

            return (
                [html.P("Uploaded sequences:"), html.Br()] + [html.Li(fasta_id) for fasta_id in fasta_ids],
                dict(style=BUTTON_STYLE),
            )
        return None, dict(display="none")

    # Function to upload new sequences
    @app.callback(
        Output("upload-structures-list", "children"),
        [
            Input("upload-structures-upload", "filename"),
            Input("upload-structures-upload", "contents"),
        ],
    )
    def upload_structures(file_names, contents):
        if file_names is not None and contents is not None:
            for file_name, content in zip(file_names, contents):
                content_type, content_string = content.split(",")
                decoded = base64.b64decode(content_string)
                with open(
                    os.path.join("uploaded_data/structures", file_name), "wb"
                ) as pdb:
                    pdb.write(decoded)

            return [html.P("Uploaded structures:"), html.Br()] + [html.Li(file_name)
                for file_name in os.listdir("uploaded_data/structures")]
        return None

    @app.callback(
        [
            ServersideOutput("input-turterra-data", "data"),
            Output("phylogeny-tree-viewer", "elements"),
        ],
        [
            Input("load-zip", "filename"),
            Input("load-zip", "contents"),
        ],
    )
    def upload_zipped_data(file_name, zip_contents):
        if file_name and zip_contents:
            if file_name.endswith('.zip'):
                content_type, content_string = zip_contents.split(',')
                decoded = base64.b64decode(content_string)
                zip_str = io.BytesIO(decoded)
                folder = '.'.join(file_name.split('.')[:-1])

                with ZipFile(zip_str) as z:
                    z.extractall("data")

                directory = os.path.join("data", folder)
                turterra_data = TurterraData.from_folder(directory)

                if turterra_data:
                    tree_elements = {
                        "nodes": turterra_data.tree.nodes,
                        "edges": turterra_data.tree.edges,
                    }
                    return [turterra_data, tree_elements]

        return [None, []]


def get_sequence_viewer_information(accession: str, turterra_data: TurterraData):
    return html.Div(
        children=[
            html.Div(f"Accession: {accession}"),
            html.Div(f"Organism: {turterra_data.species[accession]}"),
        ]
        + [
            html.Div(f"{header}: {turterra_data.extra_data[accession][header]}")
            for header in turterra_data.extra_headers
            if accession in turterra_data.extra_data
        ]
    )


def get_selected_residues_information(residues):
    return html.Div(
        [html.Label("Selected Residues:"), html.Div("{}".format(", ".join(residues))),]
    )


def get_structure_information(accession: str, turterra_data: TurterraData):
    if turterra_data.structures[accession].is_model:
        templates_str = ", ".join(turterra_data.structures[accession].templates)
        return html.Label(
            f"Note: This structure is constructed through homology modelling using MODELLER 9.23 with {templates_str} as templates",
            style={"color": "red", "textAlign": "center"},
        )
    else:
        return html.Label(
            f"Structure PDB ID: {turterra_data.structures[accession].pdb_id}",
            style={"color": "green", "textAlign": "center"},
        )


def get_compound_radio_buttons(accession: str, turterra_data: TurterraData):
    compounds = turterra_data.compounds_mapping[accession]
    return dcc.RadioItems(
        id="compound-radio-select-button",
        options=[{"label": c, "value": c,} for c in compounds],
        value=compounds[0],
    )
