import dash_core_components as dcc
import dash_daq as daq
import dash_cytoscape
import dash_html_components as html
from dash_extensions import Download

from turterra.utils import tree_utility

INTRODUCTION_TEXT = """Folder needs to have:
* sequences.fasta
* alignment.fasta
* smiles.tsv - mapping name of compound to compound smiles
* tree.txt - newick representation of tree
* structures - folder with PDB files (*_model.pdb represent homology models)
* data.txt - tab-separated file with Accession, Species, Compounds as must-have columns, and other information as optional columns
"""


def get_layout():
    return html.Div(
        children=[
            get_introduction_panel(),
            get_loading_panel(),
            html.Div(
                className="pure-g", children=[get_input_panel(), get_tree_panel(), ]
            ),
            html.Div(
                className="pure-g",
                children=[get_sequence_panel(), get_structure_panel(), ],
            ),
            html.Div(
                className="pure-g", children=[get_compound_panel(), get_upload_panel(), ]
            ),
        ]
    )


BOX_STYLE = {
    "border-width": "0.5%",
    "border-style": "solid",
    "border-color": "rgba(0,0,0,0.5)",
    "border-radius": "10px",
    "background-color": "#f9f7f7",
    "width": "45.5%",
    "padding": "2%",
}

BUTTON_STYLE = {
    "text-align": "center",
    "margin": "0 auto",
    "display": "flex",
    "justify-content": "center",
}


def get_introduction_panel():
    return html.Div(
        html.Div(
            [
                html.H1("Turterra", style={"text-align": "center"}),
                html.H5(
                    """A web portal to explore protein families""",
                    style={"text-align": "center"},
                ),
                html.Br(),
                # html.P(dcc.Markdown(INTRODUCTION_TEXT), style={"text-align": "left"}),
                html.Div(
                    [
                        html.Button(
                            "Load data",
                            id="input-load-button",
                            className="pure-button pure-button-primary",
                            style=BUTTON_STYLE,
                        ),
                        dcc.Loading(dcc.Store(id="input-turterra-data"), type="dot"),
                        html.Br(),
                    ],
                    className="row",
                ),
            ],
        ),
    )


def get_input_panel():
    return html.Div(
        id="input-panel",
        children=[
            html.H3("Filter", style={"text-align": "center"}),
            html.Br(),
            html.Div(
                html.P(
                    """Filter and select accessions below, or by clicking on nodes in the phylogenetic tree"""
                ),
            ),
            html.Div(
                [
                    dcc.Dropdown(
                        id="filter-species-dropdown",
                        options=[],
                        value=[],
                        placeholder="Filter by species",
                        multi=True,
                    ),
                    dcc.Dropdown(
                        id="filter-compound-dropdown",
                        options=[],
                        value=[],
                        placeholder="Filter by specificity",
                        multi=True,
                    ),
                    html.Br(),
                    html.Button(
                        "Filter",
                        id="filter-filter-button",
                        className="pure-button pure-button-primary",
                        style=BUTTON_STYLE,
                    ),
                    html.Br(),
                ],
            ),
            html.Div(
                [
                    dcc.Dropdown(
                        id="filter-accession-dropdown",
                        options=[],
                        value=[],
                        placeholder="Select accessions",
                        multi=True,
                        className="eight columns",
                        style={'max-height': '30vh', 'overflow': 'auto'}
                    ),
                    html.Br(),
                    html.Button(
                        children="Select all",
                        id="filter-select-all-button",
                        className="pure-button pure-button-primary",
                        style=BUTTON_STYLE,
                    ),
                    html.Br(),
                    html.Div(
                        html.Button(
                            "Show results",
                            id="filter-show-button",
                            className="pure-button pure-button-primary",
                            style=BUTTON_STYLE,
                        ),
                    ),
                    html.Br(),
                ],
            ),
            html.Br(),
        ],
        className="pure-u-1 pure-u-md-1-2",
        style=BOX_STYLE,
    )


def get_tree_panel():
    return html.Div(
        id="phylogeny-panel",
        children=[
            html.H3("Phylogeny", style={"text-align": "center"}),
            html.Br(),
            html.Div(
                html.P(
                    """Inspect phylogeny below by dragging the tree and scrolling to zoom. Any selected entries 
              after running a search above are highlighted with a red star. Additionally, click any node in the tree 
              to add this accession to the selected entries."""
                ),
            ),
            html.Br(),
            html.Button(
                "Download tree",
                id="phylogeny-download-button",
                className="pure-button pure-button-primary",
                style=BUTTON_STYLE,
            ),
            html.Br(),
            Download(id="phylogeny-download-download"),
            html.Br(),
            html.Div(
                dash_cytoscape.Cytoscape(
                    id="phylogeny-tree-viewer",
                    elements=[],
                    stylesheet=tree_utility.TREE_STYLESHEET,
                    layout={"name": "preset"},
                    minZoom=0.1,
                    maxZoom=1,
                    style={"height": "70vh", "width": "100%"},
                ),
                id="phylogeny-tree-container",
            ),
            html.Br(),
            html.Br(),
        ],
        className="pure-u-1 pure-u-md-1-2",
        style=BOX_STYLE,
    )


def get_sequence_panel():
    return html.Div(
        id="sequence-panel",
        children=[
            html.H3("Sequence and Alignments", style={"text-align": "center"}),
            html.Br(),
            html.Div(
                [
                    dcc.RadioItems(
                        id="sequence-alignment-radio-select-button",
                        options=[
                            dict(label="Single Sequence", value="Single Sequence"),
                            dict(
                                label="Sequence Alignment", value="Sequence Alignment"
                            ),
                            dict(
                                label="Structure Alignment", value="Structure Alignment"
                            ),
                        ],
                        labelStyle={"display": "inline-block"},
                        value="Single Sequence",
                    ),
                ],
            ),
            html.Br(),
            html.Div(
                [
                    dcc.Dropdown(
                        id="sequence-accession-select-dropdown",
                        options=[],
                        value=None,
                        clearable=False,
                        className="four columns",
                    )
                ],
            ),
            html.Br(),
            html.Br(),
            html.Div(
                children=[], id="sequence-sequence-viewer"
            ),
            html.Br(),
            html.Br(),
            html.Div(
                className="pure-g",
                children=[
                    html.Button(
                        "Download sequence",
                        id="sequence-single-download-button",
                        className="pure-button pure-button-primary pure-u-1-3",
                        style=BUTTON_STYLE,
                    ),
                    Download(id="sequence-single-download"),
                    html.Button(
                        "Download sequences",
                        id="sequence-multiple-download-button",
                        className="pure-button pure-button-primary pure-u-1-3",
                        style=BUTTON_STYLE,
                    ),
                    Download(id="sequence-multiple-download"),
                    html.Button(
                        "Download alignment",
                        id="sequence-alignment-download-button",
                        className="pure-button pure-button-primary pure-u-1-3",
                        style=BUTTON_STYLE,
                    ),
                    Download(id="sequence-alignment-download")])
        ],
        className="pure-u-1 pure-u-md-1-2",
        style=BOX_STYLE,
    )


def get_structure_panel():
    return html.Div(
        id="structure-panel",
        children=[
            html.H3("Structure", style={"text-align": "center"}),
            html.Br(),
            html.Div(
                [
                    html.Div(
                        [
                            dcc.Dropdown(
                                id="structure-accession-select-dropdown",
                                options=[],
                                value=None,
                                clearable=False,
                                placeholder="Select an accession"
                            ),
                            dcc.Dropdown(
                                id="structure-visualisation-select-dropdown",
                                options=[
                                    {"label": x, "value": x}
                                    for x in ["cartoon", "stick", "sphere"]
                                ],
                                value="cartoon",
                                clearable=False,
                                placeholder="Structure visualisation type",
                            ),
                            dcc.Dropdown(
                                id="structure-color-scheme-select-dropdown",
                                options=[],
                                value=None,
                                clearable=False,
                                placeholder="Structure color scheme",
                            ),
                        ],
                    ),
                    html.Br(),
                    html.Button(
                        "Load structure",
                        id="structure-refresh-button",
                        className="pure-button pure-button-primary",
                        style=BUTTON_STYLE,
                    ),
                    html.Br(),
                    dcc.Loading(dcc.Store(id="structure-structure-data"), type="dot"),
                    dcc.Store(id="structure-selected-residue-indices-data"),
                    html.Br(),
                    html.Div(
                        [
                            html.Div(id="structure-structure-container"),
                            html.P("", id="structure-model-information"),
                        ],
                    ),
                    html.Br(),
                    html.Div(
                        className="pure-g",
                        children=[
                            daq.ColorPicker(
                                id='structure-residue-color-picker',
                                label='Selected residue color',
                                value=dict(hex='#119DFF'),
                                className="pure-u-1-2"
                            ),
                            html.Div(className="pure-u-1-2", children=[dcc.Dropdown(
                                id="structure-residue-visualisation-select-dropdown",
                                options=[
                                    {"label": x, "value": x}
                                    for x in ["cartoon", "stick", "sphere"]
                                ],
                                value="cartoon",
                                clearable=False,
                                placeholder="Selected residue visualisation type",
                            )]),
                        ],
                    ),
                    html.Br(),
                    html.Div(id="structure-residue-information", style={'max-height': '300px', 'overflow': 'auto'}),
                    html.Br(),
                    html.Div(
                        className="pure-g",
                        children=[
                            html.Button(
                                "Download structure",
                                id="structure-single-download-button",
                                className="pure-button pure-button-primary pure-u-1-2",
                                style=BUTTON_STYLE,
                            ),
                            Download(id="structure-single-download"),
                            html.Button(
                                "Download structures",
                                id="structure-multiple-download-button",
                                className="pure-button pure-button-primary pure-u-1-2",
                                style=BUTTON_STYLE,
                            ),
                            Download(id="structure-multiple-download"),
                        ])
                ],
            ),
        ],
        className="pure-u-1 pure-u-md-1-2",
        style=BOX_STYLE,
    )


def get_compound_panel():
    return html.Div(
        id="compound-panel",
        children=[
            html.H3("Compound structure", style={"text-align": "center"}),
            html.Br(),
            html.Div(
                [
                    html.P("Select an accession to see its linked compounds"),
                    html.Br(),
                    dcc.Dropdown(
                        id="compound-accession-select-dropdown",
                        options=[],
                        value=None,
                        clearable=False,
                    ),
                    html.Br(),
                    html.Div(id="compound-radio-select-button-container"),
                    html.Br(),
                    html.Div(
                        id="compound-compound-container", style={"textAlign": "center"},
                    ),
                    html.Br(),
                ],
            ),
        ],
        className="pure-u-1 pure-u-md-1-2",
        style=BOX_STYLE,
    )


def get_loading_panel():
    return html.Div(
        id="loading-panel",
        children=[
            html.Div(
                html.P(
                    """Upload a .zip folder generated by turterra-build or manually constructed according to
                    guidelines described in the wiki."""
                ),
            ),
            html.Div(
                [
                    dcc.Upload(
                        id="load-zip",
                        children=html.Button(
                            "Load data",
                            id="load-zip-button",
                            className="pure-button pure-button-primary",
                            style=BUTTON_STYLE,
                        ),
                    ),
                ],
            ),

        ],
        className="pure-u-1 pure-u-md-1-2",
        style=BOX_STYLE,
    )


def get_upload_panel():
    return html.Div(
        id="upload-panel",
        children=[
            html.H3("Upload new sequences", style={"text-align": "center"}),
            html.Br(),
            html.Div(
                html.P(
                    """Upload sequences as fasta and structural models as pdb. Structural models should be of the 
                    format 'seqid_model.pdb', with seqid matching the fasta identifier corresponding to the modelled
                    sequence."""
                ),
            ),
            html.Div(
                [
                    dcc.Upload(
                        id="upload-sequences-upload",
                        children=html.Button(
                            "Upload sequences",
                            id="upload-sequences-button",
                            className="pure-button pure-button-primary",
                            style=BUTTON_STYLE,
                        ),
                    ),
                    html.Br(),
                    dcc.Upload(
                        id="upload-structures-upload",
                        children=html.Button(
                            "Upload structures",
                            id="upload-structures-button",
                            className="pure-button pure-button-primary",
                            style=BUTTON_STYLE,
                        ),
                        multiple=True,
                    ),
                ],
            ),
            html.Br(),
            html.Div(id="upload-sequences-list"),
            html.Br(),
            html.Div(id="upload-structures-list"),
            html.Div(
                id="upload-add-to-analysis",
                children=html.Button(
                    "Add to analysis",
                    id="upload-add-to-analysis-button",
                    style=dict(display="none"),
                    className="pure-button pure-button-primary",
                ),
            ),
            html.Br(),
        ],
        className="pure-u-1 pure-u-md-1-2",
        style=BOX_STYLE,
    )
