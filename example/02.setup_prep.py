from pathlib import Path
import os
import json
import copy

base_path = Path(os.getcwd())
protein_path = base_path.joinpath("dcp2_protein_only.pdb")
ligand_path = base_path.joinpath("ligands")
ppair_file = base_path.joinpath("perturbation_pair.json")
submit_path = base_path.joinpath("submit")
fep_dir = base_path.joinpath("FEP.run01")

if not os.path.exists(submit_path):
    os.mkdir(submit_path)

if not os.path.exists(fep_dir):
    os.mkdir(fep_dir)

with open(ppair_file, "r") as fp:
    perturbation_pair = json.load(fp)

submit = True
cwd = os.getcwd()

with open(base_path.joinpath("configs/config.aces.json"), "r") as fp:
    config_tmpl = json.load(fp)


for key, value in perturbation_pair.items():
    mol_names = list(value.keys())
    pair_name = "_".join(mol_names)
    pair_dir = fep_dir.joinpath(pair_name)
    if not os.path.exists(pair_dir):
        os.mkdir(pair_dir)
    
    config = copy.deepcopy(config_tmpl)
    
    config["prep"]["timask1"] = f":{mol_names[0]}"
    config["prep"]["timask2"] = f":{mol_names[1]}"
    config["prep"]["scmask1"] = f":{mol_names[0]}@{value[mol_names[0]]}"
    config["prep"]["scmask2"] = f":{mol_names[1]}@{value[mol_names[1]]}"
    if config['protocol'] == 'split':
        config["decharge"]["crgmask"] = f":2@{value[mol_names[0]]}"
        config["recharge"]["crgmask"] = f":1@{value[mol_names[1]]}"
        config["vdw_bonded"]["scmask1"] = f":1@{value[mol_names[0]]}"
        config["vdw_bonded"]["scmask2"] = f":2@{value[mol_names[1]]}"
        config["vdw_bonded"]["crgmask"] = f":1@{value[mol_names[0]]} | :2@{value[mol_names[1]]}"
    elif config['protocol'] == 'unified':
        config["unified"]["scmask1"] = f":1@{value[mol_names[0]]}"
        config["unified"]["scmask2"] = f":2@{value[mol_names[1]]}"


    with open(pair_dir.joinpath("config.json"), "w") as fp:
        json.dump(config, fp, indent=4)
    
    script = ["#!/bin/bash", "source /etc/profile", "module load cuda/11.2",  "source $HOME/env/amber22.env", "set -e"]
    script += [f"cd {pair_dir.resolve()}"]
    script += ["cp /home/gridsan/ywang3/Project/Capping/AmberTI/eIF4E/prep.py ./"]
    script += [f"python prep.py " \
               f" --lpath1 {ligand_path.joinpath(mol_names[0]).resolve()}" \
               f" --lpath2 {ligand_path.joinpath(mol_names[1]).resolve()}" \
               f" --lname1 {mol_names[0]}" \
               f" --lname2 {mol_names[1]}" \
               f" --ppdb {protein_path.resolve()}" \
               f" --config {pair_dir.joinpath('config.json').resolve()}" \
               f" --workdir {pair_dir.joinpath('run.01')}" \
               ]
    script += ["touch done"]
    submit_script = submit_path.joinpath(f"{pair_name}.sub.sh")
    with open(submit_script, "w") as fp:
        fp.write("\n".join(script))
    if submit:
        # if not os.path.exists(pair_dir.joinpath("done")):
        os.chdir(submit_path)
        # os.system(f"LLsub {submit_script.resolve()} -s 8 -g volta:1")
        os.system(f"bash {submit_script.resolve()}")
        os.chdir(cwd)

