from pathlib import Path
import os


base_path = Path(os.getcwd())
ligand_path = base_path.joinpath("ligands")
ligands = ligand_path.glob("*.sdf")
submit_path = base_path.joinpath("submit")
if not os.path.exists(submit_path):
    os.mkdir(submit_path)
    
submit = True
cwd = os.getcwd()


for ligand in ligands:
    script = ["#!/bin/bash", "source /etc/profile",  "source $HOME/env/amber22.env", "set -e"]
    this_ligand = ligand_path.joinpath(ligand.stem)
    if not os.path.exists(this_ligand):
        os.mkdir(this_ligand)
    script += [f"cd {this_ligand.resolve()}"]
    script += ["cp maketop.py ./"]
    script += [f"python maketop.py {ligand.resolve()} {ligand.stem}"]
    submit_script = submit_path.joinpath(f"{ligand.stem}.sub.sh")
    with open(submit_script, "w") as fp:
        fp.write("\n".join(script))
    if submit:
        os.chdir(submit_path)
        os.system(f"LLsub {submit_script.resolve()} -s 2")
        os.chdir(cwd)

