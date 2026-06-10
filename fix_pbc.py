import mdtraj as md

traj = md.load("cellulose_traj.dcd", top="npt_new.gro")
traj.save_xtc("cellulose_traj.xtc")