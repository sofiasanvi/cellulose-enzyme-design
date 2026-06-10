import argparse
import mdtraj as md


def main():
    parser = argparse.ArgumentParser(
        description="Convert a GROMACS .gro file to a .pdb file."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input .gro file, e.g. npt_new.gro"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output .pdb file, e.g. npt_new.pdb"
    )

    args = parser.parse_args()

    print(f"Loading GRO file: {args.input}")
    traj = md.load(args.input)

    print(f"Number of atoms: {traj.n_atoms}")
    print(f"Number of frames: {traj.n_frames}")

    print(f"Saving PDB file: {args.output}")
    traj.save_pdb(args.output)

    print("Done.")


if __name__ == "__main__":
    main()