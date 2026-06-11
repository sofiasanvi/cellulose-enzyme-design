#!/usr/bin/env python3
"""Run an OpenMM simulation from a prepared/solvated PDB file.

This script only performs simulation:
  solvated PDB -> build system -> minimize -> equilibrate -> production MD

For systems with small molecules, pass the same SDF files that were used during
solvation with --ligand-sdf. The SDF files provide ligand parameters.
"""

import argparse
from pathlib import Path

from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform
from openmm.app import (
    DCDReporter,
    ForceField,
    HBonds,
    PME,
    PDBFile,
    Simulation,
    StateDataReporter,
)
from openmm.unit import atmospheres, kelvin, nanometers, picoseconds


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run an OpenMM MD simulation from a solvated PDB file."
    )
    parser.add_argument(
        "--input-pdb",
        default="msact_complex_solvated.pdb",
        help="Prepared/solvated input PDB. Default: msact_complex_solvated.pdb",
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
    parser.add_argument("--steps", type=int, default=100000, help="Production MD steps.")
    parser.add_argument(
        "--equilibration-steps",
        type=int,
        default=10000,
        help="Equilibration steps before production.",
    )
    parser.add_argument(
        "--minimize-iterations",
        type=int,
        default=1000,
        help="Maximum energy minimization iterations.",
    )
    parser.add_argument("--dt-ps", type=float, default=0.002, help="Timestep in ps.")
    parser.add_argument(
        "--temperature-k",
        type=float,
        default=300.0,
        help="Simulation temperature in K.",
    )
    parser.add_argument(
        "--friction-per-ps",
        type=float,
        default=1.0,
        help="Langevin friction coefficient in 1/ps.",
    )
    parser.add_argument(
        "--pressure-atm",
        type=float,
        default=1.0,
        help="Pressure in atm for Monte Carlo barostat.",
    )
    parser.add_argument(
        "--barostat-interval",
        type=int,
        default=25,
        help="Barostat attempt interval in steps.",
    )
    parser.add_argument(
        "--nonbonded-cutoff-nm",
        type=float,
        default=1.0,
        help="Nonbonded cutoff in nm.",
    )
    parser.add_argument(
        "--platform",
        default="CPU",
        help="OpenMM platform name, e.g. CPU, CUDA, OpenCL. Default: CPU",
    )
    parser.add_argument(
        "--traj-output",
        default="msact_traj.dcd",
        help="Output DCD trajectory file. Default: msact_traj.dcd",
    )
    parser.add_argument(
        "--data-output",
        default="msact_data.csv",
        help="Output state data file. Default: msact_data.csv",
    )
    parser.add_argument(
        "--report-interval",
        type=int,
        default=1000,
        help="Reporter interval in steps. Default: 1000",
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

    print("Loading prepared PDB...")
    pdb = PDBFile(str(input_pdb))

    print("Building force field...")
    forcefield = build_forcefield(args)

    temperature = args.temperature_k * kelvin
    pressure = args.pressure_atm * atmospheres

    print("Building OpenMM system...")
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=args.nonbonded_cutoff_nm * nanometers,
        constraints=HBonds,
        rigidWater=True,
    )
    system.addForce(MonteCarloBarostat(pressure, temperature, args.barostat_interval))

    integrator = LangevinMiddleIntegrator(
        temperature,
        args.friction_per_ps / picoseconds,
        args.dt_ps * picoseconds,
    )

    platform = Platform.getPlatformByName(args.platform)
    simulation = Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)

    print("Performing energy minimization...")
    simulation.minimizeEnergy(maxIterations=args.minimize_iterations)

    print("Equilibrating...")
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.step(args.equilibration_steps)

    print("Running production simulation...")
    simulation.reporters.append(DCDReporter(args.traj_output, args.report_interval))
    simulation.reporters.append(
        StateDataReporter(
            args.data_output,
            args.report_interval,
            totalSteps=args.steps,
            step=True,
            time=True,
            potentialEnergy=True,
            temperature=True,
            density=True,
            speed=True,
            separator="\t",
        )
    )

    simulation.currentStep = 0
    simulation.step(args.steps)
    print("Done.")


if __name__ == "__main__":
    main()
