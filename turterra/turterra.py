from dataclasses import dataclass, field
import typing
import os
from pathlib import Path
import pandas as pnd
from turterra.utils import (
    structure_utility,
    tree_utility,
    compound_utility,
    sequence_utility,
)
import typer
from caretta import multiple_alignment
from shutil import copyfile
import pickle


@dataclass
class TurterraData:
    accessions: typing.List[str]
    sequences: typing.Dict[str, str]
    compounds_dict: typing.Dict[str, compound_utility.Compound]
    compounds_mapping: typing.Dict[str, typing.List[str]]
    structures: typing.Dict[str, structure_utility.Structure]
    species: typing.Dict[str, str]
    tree: tree_utility.Tree
    sequence_alignment: typing.Dict[str, str]
    structure_alignment: typing.Dict[str, str]
    structure_alignment_object: multiple_alignment.StructureMultiple
    extra_headers: typing.List[str] = field(default_factory=list)
    extra_data: typing.Dict[str, typing.Dict[str, typing.Any]] = field(
        default_factory=dict
    )

    @classmethod
    def from_folder(cls, folder: typing.Union[str, Path]):
        folder = Path(folder)
        sequences_file = folder / "sequences.fasta"
        sequences = sequence_utility.get_sequences_from_fasta(sequences_file)
        typer.echo("Loaded sequences")
        smiles_file = folder / "smiles.tsv"
        smiles_df = pnd.read_csv(smiles_file, sep="\t")
        compounds_dict = dict()
        with typer.progressbar(
            zip(smiles_df["Name"], smiles_df["SMILES"]), label="Loading SMILES"
        ) as progress:
            for (name, smiles) in progress:
                compounds_dict[name] = compound_utility.Compound.from_smiles(
                    name, smiles
                )
        tree_file = folder / "tree.txt"
        tree = tree_utility.Tree.from_newick(tree_file)
        typer.echo("Loaded tree")

        alignment_file = folder / "sequence_alignment.fasta"
        sequence_alignment = sequence_utility.get_sequences_from_fasta(alignment_file)
        typer.echo("Loaded sequence alignment")

        data_file = folder / "full_data.txt"
        dataframe = pnd.read_csv(data_file, sep="\t")
        necessary_columns = {"Accession", "Species", "Compounds"}
        assert all(c in dataframe.columns for c in necessary_columns)

        species = dict(zip(dataframe["Accession"], dataframe["Species"]))
        compounds = {}
        with typer.progressbar(
            list(zip(dataframe["Accession"], dataframe["Compounds"])),
            label="Mapping compounds",
        ) as progress:
            for name, compound_names in progress:
                compound_names = compound_names.replace('"', "").split(", ")
                compounds[name] = [c for c in compound_names if c in compounds_dict]

        structure_folder = folder / "structures"
        structure_alignment_folder = folder / "caretta_results"
        structure_alignment_file = folder / "structure_alignment.fasta"
        structure_alignment_object_file = os.path.join(folder, 'structure_alignment_object.pickle')

        if not structure_alignment_file.exists():
            multiple_alignment.trigger_numba_compilation()
            structure_alignment_object = multiple_alignment.StructureMultiple.align_from_pdb_files(
                structure_folder,
                output_folder=structure_alignment_folder,
                write_fasta=True,
                write_pdb=True,
                write_features=False,
                only_dssp=True,
                write_class=False,
                write_matrix=True,
                verbose=True,
            )

            with open(structure_alignment_object_file, 'wb') as pickled_alignment:
                pickle.dump(structure_alignment_object, pickled_alignment)

            copyfile(
                structure_alignment_folder / "result.fasta",
                folder / "structure_alignment.fasta",
            )
            typer.echo("Made structure alignment")

        with open(structure_alignment_object_file, 'rb') as pickled_alignment:
            structure_alignment_object = pickle.load(pickled_alignment)

        structures = dict()
        with typer.progressbar(
            list((structure_alignment_folder / "superposed_pdbs").glob("*.pdb")),
            label="Loading structures",
        ) as progress:
            for structure_file in progress:
                structure = structure_utility.Structure.from_file(structure_file)
                structures[structure.name] = structure
        structure_alignment = sequence_utility.get_sequences_from_fasta(
            structure_alignment_file
        )
        structure_alignment = {
            k.replace("_model", ""): v for k, v in structure_alignment.items()
        }
        typer.echo("Loaded structure alignment")

        extra_headers = [c for c in dataframe.columns if c not in necessary_columns]
        extra_data = {
            h: dict(zip(dataframe["Accession"], dataframe[h])) for h in extra_headers
        }
        typer.echo("Loaded extra data")

        return cls(
            list(dataframe["Accession"]),
            sequences,
            compounds_dict,
            compounds,
            structures,
            species,
            tree,
            sequence_alignment,
            structure_alignment,
            structure_alignment_object,
            extra_headers,
            extra_data,
        )

    def add_data_points(self, upload_folder, data_folder):
        new_sequences_file = os.path.join(upload_folder, "uploaded_sequences.fasta")

        for id, sequence in sequence_utility.get_sequences_from_fasta(new_sequences_file).items():
            self.sequences[id] = sequence
            self.compounds_mapping[id] = ['Unknown']
            self.species[id] = 'Unknown'
            self.accessions.append(id)

        self.compounds_dict['Unknown'] = compound_utility.Compound.from_smiles('Unknown', '')

        original_alignment = os.path.join(data_folder, "sequence_alignment.fasta")

        new_alignment_file = os.path.join(upload_folder, "updated_alignment.faa")
        sequence_utility.add_sequences_to_alignment(new_sequences_file, original_alignment, new_alignment_file)

        self.sequence_alignment = sequence_utility.get_sequences_from_fasta(new_alignment_file)

        new_structure_dir = os.path.join(upload_folder, "structures")
        new_superposed_structure_dir = os.path.join(new_structure_dir, "superposed")

        msa_class_new = multiple_alignment.StructureMultiple.from_pdb_files(new_structure_dir,
                                                                            output_folder=Path(data_folder) / "caretta_results")
        new_structure_alignment = self.structure_alignment_object.get_profile_alignment(msa_class_new, 1, 0.01)

        self.structure_alignment_object.structures += msa_class_new.structures
        self.structure_alignment_object.sequences = {**self.structure_alignment_object.sequences, **msa_class_new.sequences}

        self.structure_alignment_object.reference_structure_index = 0

        #self.structure_alignment_object.reference_structure_index = np.argmin(
       #     np.median(self.structure_alignment_object.pairwise_distance_matrix, axis=0)
       # )

        self.structure_alignment_object.write_superposed_pdbs(new_superposed_structure_dir, multiple_alignment.alignment_to_numpy(new_structure_alignment))

        self.structure_alignment = {
            k.replace("_model", ""): v for k, v in new_structure_alignment.items()
        }

        for pdb in os.listdir(new_superposed_structure_dir):
            pdb_dir = os.path.join(new_superposed_structure_dir, pdb)
            if pdb_dir[-4:] == '.pdb':
                structure = structure_utility.Structure.from_file(pdb_dir)
                self.structures[structure.name] = structure

        original_tree = os.path.join(data_folder, "tree.txt")
        epa_dir = os.path.join(upload_folder, 'epa-ng')
        tmp_dir = os.path.join(epa_dir, 'temp')
        query_alignments_dir = os.path.join(epa_dir, "query_alignments")
        new_tree = os.path.join(upload_folder, 'updated_tree.nwk')

        tree_utility.add_branches(
            new_alignment_file,
            new_sequences_file,
            query_alignments_dir,
            original_tree,
            epa_dir,
            tmp_dir,
            upload_folder)

        self.tree = tree_utility.Tree.from_newick(new_tree)

        print("Sequences added to analysis.")

    def get_sub_alignment(self, accessions, which):
        if which == "sequence":
            return sequence_utility.get_alignment_subselection(
                self.sequence_alignment, accessions, error=False
            )
        else:
            return sequence_utility.get_alignment_subselection(
                self.structure_alignment, accessions, error=False
            )









