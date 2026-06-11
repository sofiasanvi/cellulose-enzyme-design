#!/usr/bin/env python3

import argparse
from pathlib import Path

from openmm import LangevinMiddleIntegrator, Platform
from openmm.app import (
    GromacsGroFile,
    GromacsTopFile,
    Simulation,
    PME,
    HBonds,
    DCDReporter,
    StateDataReporter,
)
from openmm.unit import kelvin, picoseconds, nanometers


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run an OpenMM simulation from GROMACS .gro and .top files."
    )

    parser.add_argument(
        "--input-gro",
        default="cellulose_water.gro",
        help="Input GROMACS .gro file. Default: cellulose_water.gro"
    )

    parser.add_argument(
        "--input-top",
        default="topol_cellulose_water.top",
        help="Input GROMACS .top topology file. Default: topol_cellulose_water.top"
    )

    parser.add_argument(
        "--include-dir",
        default="charmm36-jul2021",
        help="Directory containing included force-field files. Default: charmm36-jul2021"
    )

    parser.add_argument(
        "--steps",
        type=int,
        default=100000,
        help="Number of production MD steps. Default: 100000"
    )

    parser.add_argument(
        "--equilibration-steps",
        type=int,
        default=10000,
        help="Number of equilibration steps. Default: 10000"
    )

    parser.add_argument(
        "--minimize-iterations",
        type=int,
        default=1000,
        help="Maximum minimization iterations. Default: 1000"
    )

    parser.add_argument(
        "--temperature-k",
        type=float,
        default=300.0,
        help="Temperature in Kelvin. Default: 300"
    )

    parser.add_argument(
        "--friction-per-ps",
        type=float,
        default=1.0,
        help="Langevin friction coefficient in 1/ps. Default: 1.0"
    )

    parser.add_argument(
        "--dt-ps",
        type=float,
        default=0.002,
        help="Timestep in picoseconds. Default: 0.002"
    )

    parser.add_argument(
        "--nonbonded-cutoff-nm",
        type=float,
        default=1.0,
        help="Nonbonded cutoff in nm. Default: 1.0"
    )

    parser.add_argument(
        "--platform",
        default="CPU",
        help="OpenMM platform: CPU, CUDA, or OpenCL. Default: CPU"
    )

    parser.add_argument(
        "--traj-output",
        default="cellulose_traj.dcd",
        help="Output trajectory DCD file. Default: cellulose_traj.dcd"
    )

    parser.add_argument(
        "--data-output",
        default="cellulose_data.csv",
        help="Output data CSV file. Default: cellulose_data.csv"
    )

    parser.add_argument(
        "--report-interval",
        type=int,
        default=1000,
        help="Reporter interval in steps. Default: 1000"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    input_gro = Path(args.input_gro)
    input_top = Path(args.input_top)
    include_dir = Path(args.include_dir)

    if not input_gro.exists():
        raise FileNotFoundError(f"Could not find GRO file: {input_gro}")

    if not input_top.exists():
        raise FileNotFoundError(f"Could not find TOP file: {input_top}")

    if not include_dir.exists():
        raise FileNotFoundError(f"Could not find include directory: {include_dir}")

    print(f"Loading GRO file: {input_gro}")
    gro = GromacsGroFile(str(input_gro))

    print(f"Loading topology file: {input_top}")
    print(f"Using include directory: {include_dir}")

    top = GromacsTopFile(
        str(input_top),
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir=str(include_dir),
    )

    print("Building OpenMM system...")
    system = top.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=args.nonbonded_cutoff_nm * nanometers,
        constraints=HBonds,
    )

    temperature = args.temperature_k * kelvin

    integrator = LangevinMiddleIntegrator(
        temperature,
        args.friction_per_ps / picoseconds,
        args.dt_ps * picoseconds,
    )

    print(f"Using platform: {args.platform}")
    platform = Platform.getPlatformByName(args.platform)

    simulation = Simulation(
        top.topology,
        system,
        integrator,
        platform,
    )

    simulation.context.setPositions(gro.positions)

    print("Minimizing...")
    simulation.minimizeEnergy(maxIterations=args.minimize_iterations)

    print("Equilibrating...")
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.step(args.equilibration_steps)

    print("Adding reporters...")

    simulation.reporters.append(
        DCDReporter(args.traj_output, args.report_interval, enforcePeriodicBox=False)
    )

    simulation.reporters.append(
        StateDataReporter(
            args.data_output,
            args.report_interval,
            step=True,
            time=True,
            potentialEnergy=True,
            temperature=True,
            speed=True,
            separator="\t",
        )
    )

    print("Running production MD...")
    simulation.step(args.steps)

    print("Done.")
    print(f"Trajectory written to: {args.traj_output}")
    print(f"Data written to: {args.data_output}")


if __name__ == "__main__":
    main()