import pytraj as pt
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


fep_dir = Path("/home/gridsan/ywang3/Project/Capping/AmberTI/eIF4E/FEP.run03/M7G_OME/run.01/fep")
prep_dir = Path("/home/gridsan/ywang3/Project/Capping/AmberTI/eIF4E/FEP.run03/M7G_OME/run.01/prepration")
tasks = ["unified"]
pdb_output = fep_dir.joinpath("pdb_traj")
if not os.path.exists(pdb_output):
    os.mkdir(pdb_output)

for system in ["ligands", "complex"]:
    system_dir = fep_dir.joinpath(system)
    for task in tasks:
        task_dir = system_dir.joinpath(task)
        lambda_dirs = task_dir.glob("[0-9]*[0-9]")
        for lmd_dir in lambda_dirs:
            print(lmd_dir)
            traj = pt.iterload(str(lmd_dir.joinpath("prod.nc")), str(prep_dir.joinpath(f"{system}_{task}.parm7")))
            print("loaded.")
            traj.autoimage()
            print("autoimaged")
            if system == "complex":
                traj = pt.align(traj, mask="@CA", ref=-1)
            else:
                traj = pt.align(traj, mask=":1,2", ref=-1)
            # pt.write_traj(str(lmd_dir.joinpath("prod.new.pdb")), traj["!:WAT,Na+,Cl-"], options='model', frame_indices=[i*100 for i in range(10)], overwrite=True)
            pt.write_traj(str(pdb_output.joinpath(f"prod.{system}.{task}.{lmd_dir.name}.pdb")), traj["!:WAT,Na+,Cl-"], options='model', overwrite=True)

