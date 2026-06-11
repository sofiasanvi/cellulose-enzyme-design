from openmm import *
from openmm.app import *
from openmm.unit import *

from openff.toolkit import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

# -------------------------
# Input files
# -------------------------

pdb = PDBFile("2Q0Q_chainA_clean_MET.pdb")
methyl_acetate = Molecule.from_file("Conformer3D_COMPOUND_CID_6584.sdf")
beta_glucose = Molecule.from_file("Conformer3D_COMPOUND_CID_64689.sdf")

# Small-molecule force field for methyl acetate + glucose
smirnoff = SMIRNOFFTemplateGenerator(
    molecules=[methyl_acetate, beta_glucose],
    forcefield="openff-2.2.1"   # use the OpenFF version installed in your env
)

# Protein + water force field
forcefield = ForceField(
    "amber14-all.xml",
    "amber14/tip3p.xml"
)

# Register ligand parameter generator
forcefield.registerTemplateGenerator(smirnoff.generator)

# -------------------------
# Model preparation
# -------------------------

print("Preparing model...")
modeller = Modeller(pdb.topology, pdb.positions)

# Add hydrogens to protein/ligands where possible
modeller.addHydrogens(forcefield, pH=7.0)

# Add water box and ions
modeller.addSolvent(
    forcefield,
    model="tip3p",
    padding=1.0*nanometers,
    ionicStrength=0.15*molar
)

# Save prepared system for inspection
with open("msact_complex_solvated.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

# -------------------------
# System configuration
# -------------------------

nonbondedMethod = PME
nonbondedCutoff = 1.0*nanometers
constraints = HBonds
rigidWater = True

dt = 0.002*picoseconds
temperature = 300*kelvin
friction = 1.0/picosecond
pressure = 1.0*atmospheres
barostatInterval = 25

# Start short for testing
steps = 100000          # 100k steps = 200 ps with 2 fs timestep
equilibrationSteps = 10000

platform = Platform.getPlatformByName("CPU")

# -------------------------
# Build system
# -------------------------

print("Building system...")
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=nonbondedMethod,
    nonbondedCutoff=nonbondedCutoff,
    constraints=constraints,
    rigidWater=rigidWater
)

system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))

integrator = LangevinMiddleIntegrator(temperature, friction, dt)

simulation = Simulation(
    modeller.topology,
    system,
    integrator,
    platform
)

simulation.context.setPositions(modeller.positions)

# -------------------------
# Minimize and equilibrate
# -------------------------

print("Performing energy minimization...")
simulation.minimizeEnergy(maxIterations=1000)

print("Equilibrating...")
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(equilibrationSteps)

# -------------------------
# Production simulation
# -------------------------

print("Simulating...")

simulation.reporters.append(DCDReporter("msact_traj.dcd", 1000))
simulation.reporters.append(StateDataReporter(
    "msact_data.csv",
    1000,
    totalSteps=steps,
    step=True,
    time=True,
    potentialEnergy=True,
    temperature=True,
    density=True,
    speed=True,
    separator="\t"
))

simulation.step(steps)

print("Done.")