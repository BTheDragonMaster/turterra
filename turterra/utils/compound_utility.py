from dataclasses import dataclass
import turterra.dependencies.pikachu.pikachu as pikachu


@dataclass
class Compound:
    name: str
    plot_data: dict

    @classmethod
    def from_smiles(cls, name: str, smiles_string: str):
        return cls(name, smile_to_dict(smiles_string))


def smile_to_dict(smiles):
    """ Convert a smiles string into a nodes and links dictionary
  :param smile: Smiles string
  :return: Dict of: - nodes: list of dicts containing id and atom type
                    - links: list of dicts containing id, source id, target id and bond strength
  """
    if smiles:
        mol = pikachu.read_smiles(smiles)
        dash_input = mol.to_dash_molecule2d_input()
        return dash_input

    else:
        return {'nodes': [], 'links': []}
