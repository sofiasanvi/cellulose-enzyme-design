#!/usr/bin/env python3

import argparse


def filter_gro(input_file, output_file, keep_resnames):
    keep_resnames = set(keep_resnames)

    with open(input_file, "r") as f:
        title = f.readline()
        n_atoms = int(f.readline().strip())

        atom_lines = [f.readline() for _ in range(n_atoms)]
        box_line = f.readline()

    kept_lines = []

    for line in atom_lines:
        # GRO-format:
        # columns 1-5   residue number
        # columns 6-10  residue name
        # columns 11-15 atom name
        # columns 16-20 atom number
        resname = line[5:10].strip()

        if resname in keep_resnames:
            kept_lines.append(line)

    # Renumber atom numbers so the output GRO is consistent
    renumbered_lines = []
    for i, line in enumerate(kept_lines, start=1):
        atom_number = i % 100000
        new_line = line[:15] + f"{atom_number:5d}" + line[20:]
        renumbered_lines.append(new_line)

    with open(output_file, "w") as f:
        f.write(title)
        f.write(f"{len(renumbered_lines):5d}\n")
        for line in renumbered_lines:
            f.write(line)
        f.write(box_line)

    print(f"Wrote: {output_file}")
    print(f"Kept residue names: {', '.join(sorted(keep_resnames))}")
    print(f"Kept atoms: {len(renumbered_lines)}")


def main():
    parser = argparse.ArgumentParser(
        description="Filter a GROMACS .gro file by residue name."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input GRO file, e.g. npt_new.gro"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output GRO file, e.g. cellulose_water.gro"
    )

    parser.add_argument(
        "--keep",
        nargs="+",
        required=True,
        help="Residue names to keep, e.g. BGLC SOL"
    )

    args = parser.parse_args()

    filter_gro(
        input_file=args.input,
        output_file=args.output,
        keep_resnames=args.keep
    )


if __name__ == "__main__":
    main()