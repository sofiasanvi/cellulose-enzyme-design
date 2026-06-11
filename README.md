# cellulose-enzyme-design

## Project Structure

```text
cellulose-enzyme-design/
├── config/                  # Files controlling the simulations (GROMACS/OpenMM)
│   ├── charmm36-jul2021/    # Force field directory
│   ├── mdout.mdp
│   └── npt_new.mdp
│
├── data/                    # Input structures, topologies, and ligands
│   ├── structures/          # PDB and GRO files (.pdb, .gro)
│   ├── topology/            # Topology files (.top, .itp, .prm)
│   └── ligands/             # Small molecules (.sdf)
│
├── scripts/                 # Python scripts for running and analyzing simulations
│   ├── convert_gro_to_pdb.py
│   ├── filter_gro_residues.py
│   ├── fix_pbc.py
│   └── running_openmm.py
│
├── trajectories/            # Large simulation trajectory files
│   ├── dcd/                 # DCD trajectory files
│   └── xtc/                 # XTC trajectory files
│
├── results/                 # Analysis data, plots, and logs
│   ├── cellulose_data.csv
│   └── msact_data.csv
│
├── movies/                  # Rendered videos and slow-motion movies
│
├── .gitignore
├── LICENSE
└── README.md
```

