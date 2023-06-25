import os
import shutil
from pathlib import Path
from amberti.md import em, heat, pressurize, production
from amberti.logger import getLogger
from amberti.utils import set_directory, set_tag, print_list
from amberti.op import (
    create_simulation_box, 
    make_charge_transform, 
    concatenate_pdb, 
    extract_conf
)


logger = getLogger()


def equilibrate_systems(
        config,
        top_dir   
    ):
    systems = ["ligands", "complex"]

    cwd = Path(os.getcwd())
    for system in systems:
        system_dir = cwd.joinpath(system)
        if not os.path.exists(system_dir):
            os.mkdir(system_dir)
        prmtop = top_dir.joinpath(f"{system}_vdw_bonded.parm7")
        with set_directory(system_dir):
            em(
                "min",
                prmtop=prmtop, 
                conf=top_dir.joinpath(f"{system}_vdw_bonded.rst7"),
                ref=top_dir.joinpath(f"{system}_vdw_bonded.rst7"),
                maxcyc=config["prep"]["em"]["maxcyc"],
                resstraint_wt=config["prep"]["em"]["resstraint_wt"],
                fep=True,
                clambda=0.5,
                scalpha=0.5,
                scbeta=12.0,
                ifsc=1,
                timask1=config["prep"]["timask1"],
                timask2=config["prep"]["timask2"],
                scmask1=config["prep"]["scmask1"],
                scmask2=config["prep"]["scmask2"],
                fname="min.in",
                run=True
            )

            logger.info("Pre-heating ...")
            heat(
                defname="pre_heating",
                prmtop=prmtop, 
                conf="min.rst7",
                ref="min.rst7",
                nsteps=20,
                dt=0.0002,
                temp=5.0,
                resstraint_wt=config["prep"]["heat"]["resstraint_wt"],
                fep=True,
                clambda=0.5,
                scalpha=0.5,
                scbeta=12.0,
                ifsc=1,
                timask1=config["prep"]["timask1"],
                timask2=config["prep"]["timask2"],
                scmask1=config["prep"]["scmask1"],
                scmask2=config["prep"]["scmask2"],
                tempi=4.0,
                ofreq=1,
                fname="pre_heat.in",
                run=True
            )

            logger.info("Second Emergy minimizing ...")
            em(
                "min2",
                prmtop=prmtop,
                conf="pre_heating.rst7",
                ref="pre_heating.rst7",
                maxcyc=config["prep"]["em"]["maxcyc"],
                resstraint_wt=config["prep"]["em"]["resstraint_wt"],
                fep=True,
                clambda=0.5,
                scalpha=0.5,
                scbeta=12.0,
                ifsc=1,
                timask1=config["prep"]["timask1"],
                timask2=config["prep"]["timask2"],
                scmask1=config["prep"]["scmask1"],
                scmask2=config["prep"]["scmask2"],
                fname="min2.in",
                run=True
            )

            logger.info(f"Heating from 5K to {config['prep']['temp']}K ...")
            heat(
                defname="heat",
                prmtop=prmtop, 
                conf="min2.rst7",
                ref="min2.rst7",
                nsteps=config["prep"]["heat"]["nsteps"],
                temp=config['prep']['temp'],
                resstraint_wt=config["prep"]["heat"]["resstraint_wt"],
                fep=True,
                clambda=0.5,
                scalpha=0.5,
                scbeta=12.0,
                ifsc=1,
                timask1=config["prep"]["timask1"],
                timask2=config["prep"]["timask2"],
                scmask1=config["prep"]["scmask1"],
                scmask2=config["prep"]["scmask2"],
                tempi=5.0,
                ofreq=config["prep"]["heat"]["ofreq"],
                fname="heat.in",
                run=True
            )

            logger.info("Pre-Pressurising ...")
            pressurize(
                defname="pre_press",
                prmtop=prmtop, 
                conf="heat.rst7",
                ref="heat.rst7",
                nsteps=1000,
                dt=0.002,
                temp=config['prep']['temp'],
                resstraint_wt=config['prep']['pressurize']['resstraint_wt'],
                irest=1, ntx=5,
                fep=True,
                clambda=0.5,
                scalpha=0.5,
                scbeta=12.0,
                ifsc=1,
                timask1=config["prep"]["timask1"],
                timask2=config["prep"]["timask2"],
                scmask1=config["prep"]["scmask1"],
                scmask2=config["prep"]["scmask2"],
                ofreq=10,
                fname="pre_press.in",
                run=True
            )

            logger.info("Pressurising ...")
            pressurize(
                defname="pressurize",
                prmtop=prmtop, 
                conf="pre_press.rst7",
                ref="pre_press.rst7",
                nsteps=config["prep"]["pressurize"]["nsteps"],
                dt=0.002,
                temp=config["prep"]["temp"],
                resstraint_wt=config["prep"]["pressurize"]["resstraint_wt"],
                irest=1, ntx=5,
                fep=True,
                clambda=0.5,
                scalpha=0.5,
                scbeta=12.0,
                ifsc=1,
                timask1=config["prep"]["timask1"],
                timask2=config["prep"]["timask2"],
                scmask1=config["prep"]["scmask1"],
                scmask2=config["prep"]["scmask2"],
                ofreq=config["prep"]["pressurize"]["ofreq"],
                fname="pressurize.in",
                run=True
            )

            logger.info("Equilibrium finished.")

    return


def prep(ppdb, lpath1, lname1, lpath2, lname2, config):

    cwd = Path(os.getcwd())
    lpath1 = Path(lpath1)
    lpath2 = Path(lpath2)
    ligand_forcefield = config["ligand_forcefield"]
    protein_forcefield = config["protein_forcefield"]

    lib1 = lpath1.joinpath(f"{lname1}.lib")
    frcmod1 = lpath1.joinpath(f"{lname1}.{ligand_forcefield}.frcmod")
    lpdb1 = lpath1.joinpath(f"{lname1}.pdb")
    lib2 = lpath2.joinpath(f"{lname2}.lib")
    frcmod2 = lpath2.joinpath(f"{lname2}.{ligand_forcefield}.frcmod")
    lpdb2 = lpath2.joinpath(f"{lname2}.pdb")
    llpdb = cwd.joinpath(f"{lname1}_{lname2}.pdb")
    logger.info("Concatenating the two ligand PDBs.")
    concatenate_pdb(lpdb1, lpdb2, llpdb)
    logger.info("Creating the simulation boxes.")
    create_simulation_box(
                          lib1, lib2, 
                          frcmod1, frcmod2, 
                          lpdb1, lpdb2, ppdb, llpdb,
                          ligand_forcefield,
                          protein_forcefield,
                          water='tip3p',
                          size_ligand=config["ligand_box_size"], 
                          size_complex=config["complex_box_size"], 
                          resize=config["resize"]
    )
    logger.info("Start equilibrating.")
    equilibrate_systems(config=config, top_dir=cwd)

    for system in ["ligands", "complex"]:
        rst7 = f"{system}_vdw_bonded.rst7"
        if os.path.exists(rst7):
            os.rename(rst7, rst7+".leap")
        shutil.copy(f"{system}/pressurize.rst7", f"{system}_vdw_bonded.rst7")
        extract_conf(f"{system}_vdw_bonded.rst7", f"{system}_vdw_bonded.parm7", 
                     f"{system}_solvated.pdb", select="!(:1,2)", overwrite=True)
        extract_conf(f"{system}_vdw_bonded.rst7", f"{system}_vdw_bonded.parm7", 
                     f"{system}_{lname1}.pdb", select=":1", overwrite=True)
        extract_conf(f"{system}_vdw_bonded.rst7", f"{system}_vdw_bonded.parm7", 
                     f"{system}_{lname2}.pdb", select=":2", overwrite=True)
    
    make_charge_transform(
            lib1, lib2, 
            frcmod1, frcmod2, 
            "ligands_solvated.pdb", f"ligands_{lname1}.pdb", f"ligands_{lname2}.pdb",
            "complex_solvated.pdb", f"complex_{lname1}.pdb", f"complex_{lname2}.pdb",
            ligand_forcefield,
            protein_forcefield,
            water='tip3p',
    )

    return 


def setTI(
        config,
        top_dir,
    ):
    """
        To shedule the lambda list, a dict is desired here.
    """
    systems = ["ligands", "complex"]
    steps = ["decharge", "vdw_bonded", "recharge"]
    top_dir = Path(top_dir)
    tasks = {}
    
    cwd = Path(os.getcwd())
    for system in systems:
        system_dir = cwd.joinpath(system)
        if not os.path.exists(system_dir):
            os.mkdir(system_dir)

        for step in steps:
            step_dir = system_dir.joinpath(step)
            if not os.path.exists(step_dir):
                os.mkdir(step_dir)
            
            for lmb in config[step]["lambdas"]:
                lmb_dir = step_dir.joinpath(str(lmb))
                if not os.path.exists(lmb_dir):
                    os.mkdir(lmb_dir)
                prmtop = top_dir.joinpath(f"{system}_{step}.parm7")
                tasks[f"{system}_{step}_{lmb}"] = {}
                tasks[f"{system}_{step}_{lmb}"]["path"] = lmb_dir.resolve()
                tasks[f"{system}_{step}_{lmb}"]["command"] = []

                with set_directory(lmb_dir):
                    tasks[f"{system}_{step}_{lmb}"]["command"].append( em(
                        "min",
                        prmtop=prmtop, 
                        conf=top_dir.joinpath(f"{system}_{step}.rst7"),
                        ref=top_dir.joinpath(f"{system}_{step}.rst7"),
                        maxcyc=config["TI"]["em"]["maxcyc"],
                        resstraint_wt=config["TI"]["em"]["resstraint_wt"],
                        fep=True,
                        clambda=lmb,
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        fname="min.in",
                        run=False
                    ))

                    logger.info("Pre-heating ...")
                    tasks[f"{system}_{step}_{lmb}"]["command"].append(heat(
                        defname="pre_heating",
                        prmtop=prmtop, 
                        conf="min.rst7",
                        ref="min.rst7",
                        nsteps=20,
                        dt=0.0002,
                        temp=5.0,
                        resstraint_wt=config["TI"]["heat"]["resstraint_wt"],
                        fep=True,
                        clambda=lmb,
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        tempi=4.0,
                        ofreq=1,
                        fname="pre_heat.in",
                        run=False
                    ))

                    logger.info("Second Emergy minimizing ...")
                    tasks[f"{system}_{step}_{lmb}"]["command"].append( em(
                        "min2",
                        prmtop=prmtop,
                        conf="pre_heating.rst7",
                        ref="pre_heating.rst7",
                        maxcyc=config["TI"]["em"]["maxcyc"],
                        resstraint_wt=config["TI"]["em"]["resstraint_wt"],
                        fep=True,
                        clambda=lmb,
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        fname="min2.in",
                        run=False
                    ))

                    logger.info(f"Heating from 5K to {config['TI']['temp']}K ...")
                    tasks[f"{system}_{step}_{lmb}"]["command"].append( heat(
                        defname="heat",
                        prmtop=prmtop, 
                        conf="min2.rst7",
                        ref="min2.rst7",
                        nsteps=config["TI"]["heat"]["nsteps"],
                        temp=config['TI']['temp'],
                        resstraint_wt=config["TI"]["heat"]["resstraint_wt"],
                        fep=True,
                        clambda=lmb,
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        tempi=5.0,
                        ofreq=config["TI"]["heat"]["ofreq"],
                        fname="heat.in",
                        run=False
                    ))

                    logger.info("Pre-Pressurising ...")
                    tasks[f"{system}_{step}_{lmb}"]["command"].append( pressurize(
                        defname="pre_pressing",
                        prmtop=prmtop, 
                        conf="heat.rst7",
                        ref="heat.rst7",
                        nsteps=1000,
                        dt=0.002,
                        temp=config['TI']['temp'],
                        resstraint_wt=config['TI']['pressurize']['resstraint_wt'],
                        irest=1, ntx=5,
                        fep=True,
                        clambda=lmb,
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        ofreq=10,
                        fname="pre_press.in",
                        run=False
                    ))

                    logger.info("Pressurising ...")
                    tasks[f"{system}_{step}_{lmb}"]["command"].append( pressurize(
                        defname="pressurize",
                        prmtop=prmtop, 
                        conf="pre_press.rst7",
                        ref="pre_press.rst7",
                        nsteps=config["TI"]["pressurize"]["nsteps"],
                        dt=0.002,
                        temp=config["TI"]["temp"],
                        resstraint_wt=config["TI"]["pressurize"]["resstraint_wt"],
                        irest=1, ntx=5,
                        fep=True,
                        clambda=lmb,
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        ofreq=config["TI"]["pressurize"]["ofreq"],
                        fname="pressurize.in",
                        run=False
                    ))

                    logger.info("Equilibrium finished.")

                    logger.info("Production run ...")
                    tasks[f"{system}_{step}_{lmb}"]["command"].append( production(
                        defname="prod",
                        prmtop=prmtop, 
                        conf="pressurize.rst7",
                        ref="pressurize.rst7",
                        nsteps=config["TI"]["production"]["nsteps"],
                        dt=0.002,
                        temp=config["TI"]["temp"],
                        resstraint_wt=config["TI"]["production"]["resstraint_wt"],
                        irest=1, ntx=5,
                        fep=True,
                        clambda=lmb,
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifmbar=1,
                        mbar_states=len(config[step]["lambdas"]),
                        mbar_lambda=print_list(config[step]["lambdas"], by=", "),
                        ifsc=config[step]["ifsc"],
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        ntwe=config["TI"]["production"]["ntwe"],
                        ntwx=config["TI"]["production"]["ntwx"],
                        ntpr=config["TI"]["production"]["ntpr"],
                        ntwr=config["TI"]["production"]["ntwr"],
                        fname="prod.in",
                        run=False
                    ))
    return tasks


def submit_jobs(tasks, env, groups=-1, submit=True):
    cwd = Path(os.getcwd())
    submit_scripts_dir = cwd.joinpath("submit")
    if not os.path.exists(submit_scripts_dir):
        os.mkdir(submit_scripts_dir)
    
    scripts = {}
    for task, info in tasks.items():
        script = ["#!/bin/bash"]
        script += env
        script.append(f"cd {info['path']}")
        script += info["command"]
        script += ["touch done.tag"]
        scripts[task] = script
    
    if groups > 0:
        ngroups = len(scripts.keys())
        task_list = list(scripts.keys())
        RuntimeError("Not implemented.")
    else:
        for task, script in scripts.items():
            sub = submit_scripts_dir.joinpath(f"{task}.sh")
            with open(sub, "w") as fp:
                fp.write("\n".join(script))

            if submit:
                if not os.path.exists(Path(tasks[task]['path']).joinpath("done.tag")):
                    with set_directory(submit_scripts_dir):
                        os.system(f"LLsub {sub.name} -s 8 -g volta:1")


def run(
        ppdb, 
        lpath1, lname1, 
        lpath2, lname2,
        config
    ):
    cwd = Path(os.getcwd())
    tag = "done.tag"

    prep_dir = cwd.joinpath("prepration")
    if not os.path.exists(prep_dir.joinpath(tag)):
        with set_directory(prep_dir):
            prep(ppdb, lpath1, lname1, lpath2, lname2, config)
            set_tag(tag)
    else:
        logger.info("preparation has already existed. skip through it.")

    fep_run_dir = cwd.joinpath("fep")
    if not os.path.exists(fep_run_dir.joinpath(tag)):
        with set_directory(fep_run_dir):
            tasks = setTI(config, top_dir=prep_dir.resolve())
            submit_jobs(tasks, config['slurm_env'], submit=False)
            set_tag(tag)
    else:
        logger.info("feprun has already existed. skip through it.")

