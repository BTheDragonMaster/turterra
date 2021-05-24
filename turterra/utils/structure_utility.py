from Bio.Data.IUPACData import protein_letters_3to1
from dataclasses import dataclass, field
import typing
from pathlib import Path
import json
from dash_bio_utils import styles_parser
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import re
import parmed as pmd

from turterra.utils import sequence_utility


@dataclass
class Structure:
    name: str
    filename: typing.Union[Path, str]
    is_model: bool
    pdb_id: typing.Union[str, None] = None
    templates: typing.List[str] = field(default_factory=list)
    ndope_score: typing.Union[float, None] = None

    @classmethod
    def from_file(cls, filename, pdb_id=None):
        filename = Path(filename)
        name = filename.stem
        if name.endswith("_model"):
            is_model = True
            name = name[: -len("_model")]
            templates = get_templates_from_model_file(filename)
            ndope_score = get_ndope_score_from_model_file(filename)
        else:
            is_model = False
            templates = []
            ndope_score = None
        return cls(name, filename, is_model, pdb_id, templates, ndope_score)

    @property
    def plot_data(self):
        return json.loads(create_data(str(self.filename)))

    def get_style_data(
        self,
        visualization_scheme,
        color_scheme,
        sequence_alignment,
        structure_alignment,
    ):
        if color_scheme in ["sequence conservation", "structure conservation"]:
            if color_scheme == "sequence conservation":
                conservations = sequence_utility.alignment_conservation(
                    sequence_alignment
                )
                start_index, end_index = map_structure_to_sequence(
                    self.plot_data, sequence_alignment[self.name].replace("-", "")
                )
            else:
                conservations = sequence_utility.alignment_conservation(
                    structure_alignment
                )
                start_index, end_index = map_structure_to_sequence(
                    self.plot_data, structure_alignment[self.name].replace("-", "")
                )
            conservations = conservations[start_index:end_index]
            viridis_r_cmap = plt.get_cmap("viridis_r", 50)

            styles_data = []
            current_residue = 0
            for x in self.plot_data["atoms"]:
                if x["residue_index"] != current_residue:
                    current_residue = x["residue_index"]
                color = to_hex(
                    viridis_r_cmap(conservations[current_residue - 1])
                )  # res_index starts at one
                styles_data.append(
                    {"color": color, "visualization_type": visualization_scheme}
                )

        else:
            styles_data = json.loads(
                styles_parser.create_style(
                    pdb_path=self.filename,
                    style=visualization_scheme,
                    mol_color=color_scheme,
                )
            )
        return styles_data


def get_templates_from_model_file(model_file):
    templates = []
    started = False
    with open(model_file) as f:
        for line in f:
            if line.startswith("REMARK"):
                if "TEMPLATE: " in line:
                    started = True
                    template = line.split("TEMPLATE: ")[1].split()[0]
                    templates.append(template)
            elif started:
                break
    return templates


def get_ndope_score_from_model_file(model_file):
    with open(model_file) as f:
        for line in f:
            if "Normalized DOPE score" in line:
                score = line.split(":")[-1].strip()
                return float(score)


def map_structure_to_sequence(structure_data, sequence):
    """
    Find the start and end indexes of the structure residue positions in the sequence

    Parameters
    ----------
    structure_data
        data of the structure
    sequence

    Returns
    -------
    start index, end index
    """
    residue_list = []
    for atom in structure_data["atoms"]:
        if atom["residue_name"] not in residue_list:
            residue_list.append(atom["residue_name"])

    structure_sequence = "".join(
        [
            protein_letters_3to1[residue[0].upper() + residue[1:3].lower()]
            for residue in residue_list
        ]
    )

    mapping_start = sequence.find(structure_sequence)
    mapping_end = mapping_start + len(structure_sequence)
    return mapping_start, mapping_end


def create_data(pdb_path):
    """
    Parse the protein data bank (PDB) file to generate
    input modelData

    @param pdb_path
    Name of the biomolecular structure file in PDB format

    """

    top = pmd.load_file(pdb_path)

    # Read PDB file to create atom/bond information
    with open(pdb_path, 'r') as infile:
        # store only non-empty lines
        lines = [l.strip() for l in infile if l.strip()]

    # Initialize all variables
    var_nchains = []
    serial = []
    atm_name = []
    res_name = []
    chain = []
    res_id = []
    positions = []
    occupancy = []
    temp_factor = []
    atom_type = []
    ct = 0

    datb = {
        'atoms': [],
        'bonds': []
    }

    # Variables that store the character positions of different
    # parameters from the molecule PDB file
    serialpos = [6, 11]
    atm_namepos = [12, 16]
    r_namepos = [17, 20]
    chainpos = [21, 22]
    r_idpos = [22, 26]
    xpos = [30, 38]
    ypos = [38, 46]
    zpos = [46, 54]
    occupos = [54, 60]
    bfacpos = [60, 66]
    atm_typepos = [77, 79]

    for l in lines:
        line = l.split()
        if "ATOM" in line[0] or "HETATM" in line[0]:
            serial.append(int(l[serialpos[0]:serialpos[1]]))
            atm_name.append(l[atm_namepos[0]:atm_namepos[1]].strip())
            val_r_name = l[r_namepos[0]:r_namepos[1]].strip()
            res_name.append(val_r_name)
            chain_val = l[chainpos[0]:chainpos[1]].strip()
            chain.append(chain_val)
            if chain_val not in var_nchains:
                var_nchains.append(chain_val)
            val_r_id = int(l[r_idpos[0]:r_idpos[1]])
            res_id.append(val_r_id)
            x = float(l[xpos[0]:xpos[1]])
            y = float(l[ypos[0]:ypos[1]])
            z = float(l[zpos[0]:zpos[1]])
            positions.append([x, y, z])
            occupancy.append(l[occupos[0]:occupos[1]].strip())
            temp_factor.append(l[bfacpos[0]:bfacpos[1]].strip())
            atom_type.append(l[atm_typepos[0]:atm_typepos[1]].strip())
            ct += 1

    # Create list of atoms
    tmp_res = res_id[0]
    resct = 1
    for i in range(len(chain)):  # pylint: disable=consider-using-enumerate
        if tmp_res != res_id[i]:
            tmp_res = res_id[i]
            resct += 1
        datb['atoms'].append({
            "name": atm_name[i],
            "chain": chain[i],
            "positions": positions[i],
            "residue_index": resct,
            "element": atom_type[i],
            "residue_name": res_name[i] + str(res_id[i]),
            "serial": i,
        })

    # Create list of bonds using the parmed module
    for i in range(len(top.bonds)):
        bondpair = top.bonds[i].__dict__
        atom1 = re.findall(r"\[(\d+)\]", str(bondpair['atom1']))
        atom2 = re.findall(r"\[(\d+)\]", str(bondpair['atom2']))
        datb['bonds'].append({
            'atom2_index': int(atom1[0]),
            'atom1_index': int(atom2[0])
        })

    return json.dumps(datb)
