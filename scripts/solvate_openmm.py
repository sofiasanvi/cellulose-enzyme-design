#!/usr/bin/env python3
"""Prepare and solvate an OpenMM system.

This script only performs model preparation:
  input PDB -> add hydrogens -> add solvent/ions -> output solvated PDB

For systems with small molecules, pass one or more SDF files with --ligand-sdf.
The ligands must already be present in the input PDB. The SDF files only provide
chemical identity/parameters; they do not place ligands into the protein.
"""

import argparse
from pathlib import Path

from openmm.app import ForceField, Modeller, PDBFile
from openmm.unit import molar, nanometers


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add hydrogens, water, and ions to a PDB file for OpenMM."
    )
    parser.add_argument(
        "--input-pdb",
        required=True,
        help="Input PDB file, e.g. msact_complex_with_ligands.pdb",
    )
    parser.add_argument(
        "--output-pdb",
        default="msact_complex_solvated.pdb",
        help="Output solvated PDB file. Default: msact_complex_solvated.pdb",
    )
    parser.add_argument(
        "--ligand-sdf",
        action="append",
        default=[],
        help=(
            "Small-molecule SDF file for OpenFF parameters. "
            "Use once per ligand, e.g. --ligand-sdf methyl_acetate.sdf "
            "--ligand-sdf beta_glucose.sdf"
        ),
    )
    parser.add_argument(
        "--protein-forcefield",
        default="amber14-all.xml",
        help="Protein force field XML. Default: amber14-all.xml",
    )
    parser.add_argument(
        "--water-forcefield",
        default="amber14/tip3p.xml",
        help="Water force field XML. Default: amber14/tip3p.xml",
    )
    parser.add_argument(
        "--small-molecule-forcefield",
        default="openff-2.2.1",
        help="OpenFF force field for ligands. Default: openff-2.2.1",
    )
    parser.add_argument(
        "--water-model",
        default="tip3p",
        help="Water model for Modeller.addSolvent(). Default: tip3p",
    )
    parser.add_argument(
        "--padding-nm",
        type=float,
        default=1.0,
        help="Water-box padding in nm. Default: 1.0",
    )
    parser.add_argument(
        "--ionic-strength-molar",
        type=float,
        default=0.15,
        help="Ionic strength in molar. Default: 0.15",
    )
    parser.add_argument(
        "--ph",
        type=float,
        default=7.0,
        help="pH used when adding hydrogens. Default: 7.0",
    )
    return parser.parse_args()


def build_forcefield(args: argparse.Namespace) -> ForceField:
    forcefield = ForceField(args.protein_forcefield, args.water_forcefield)

    if args.ligand_sdf:
        from openff.toolkit import Molecule
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator

        molecules = []
        for sdf in args.ligand_sdf:
            print(f"Reading ligand SDF: {sdf}")
            molecules.append(Molecule.from_file(sdf))

        smirnoff = SMIRNOFFTemplateGenerator(
            molecules=molecules,
            forcefield=args.small_molecule_forcefield,
        )
        forcefield.registerTemplateGenerator(smirnoff.generator)

    return forcefield


def main() -> None:
    args = parse_args()

    input_pdb = Path(args.input_pdb)
    if not input_pdb.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb}")

    print("Loading input PDB...")
    pdb = PDBFile(str(input_pdb))

    print("Building force field...")
    forcefield = build_forcefield(args)

    print("Preparing model...")
    modeller = Modeller(pdb.topology, pdb.positions)

    print(f"Adding hydrogens at pH {args.ph}...")
    modeller.addHydrogens(forcefield, pH=args.ph)

    print("Adding solvent and ions...")
    modeller.addSolvent(
        forcefield,
        model=args.water_model,
        padding=args.padding_nm * nanometers,
        ionicStrength=args.ionic_strength_molar * molar,
    )

    print(f"Writing solvated PDB: {args.output_pdb}")
    with open(args.output_pdb, "w") as handle:
        PDBFile.writeFile(modeller.topology, modeller.positions, handle)

    print("Done.")


if __name__ == "__main__":
    main()
