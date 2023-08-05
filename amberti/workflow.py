import os
import json
import shutil
import random
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
            logger.info("Pre-heating ...")
            # just for debug
            heat(
                defname="debug",
                prmtop=prmtop, 
                conf=top_dir.joinpath(f"{system}_vdw_bonded.rst7"),
                ref=top_dir.joinpath(f"{system}_vdw_bonded.rst7"),
                nsteps=5,
                dt=0.00001,
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
                fname="debug.in",
                run=True
            )

            em(
                "min",
                prmtop=prmtop, 
                conf="debug.rst7",
                ref="debug.rst7",
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

            logger.info("Pre-Pressurising1 ...")
            pressurize(
                defname="pre_press",
                prmtop=prmtop, 
                conf="heat.rst7",
                ref="heat.rst7",
                nsteps=3000,
                dt=0.002,
                temp=config['prep']['temp'],
                resstraint_wt=config['prep']['pressurize_res']['resstraint_wt'],
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

            logger.info("Pre-Pressurising2 ...")
            pressurize(
                defname="pre_press2",
                prmtop=prmtop, 
                conf="pre_press.rst7",
                ref="pre_press.rst7",
                nsteps=3000,
                dt=0.002,
                temp=config['prep']['temp'],
                resstraint_wt=config['prep']['pressurize_res']['resstraint_wt'],
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
                fname="pre_press2.in",
                run=True
            )

            logger.info("Pressurising ...")
            pressurize(
                defname="pressurize",
                prmtop=prmtop, 
                conf="pre_press2.rst7",
                ref="pre_press2.rst7",
                nsteps=config["prep"]["pressurize_res"]["nsteps"],
                dt=0.002,
                temp=config["prep"]["temp"],
                resstraint_wt=config["prep"]["pressurize_res"]["resstraint_wt"],
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
                ofreq=config["prep"]["pressurize_res"]["ofreq"],
                fname="pressurize.in",
                run=True
            )

            logger.info("Pressurising ...")
    return


def prep(ppdb, lpath1, lname1, lpath2, lname2, config):

    cwd = Path(os.getcwd())
    lpath1 = Path(lpath1)
    lpath2 = Path(lpath2)
    ligand_forcefield = config["ligand_forcefield"]
    protein_forcefield = config["protein_forcefield"]

    lib1 = lpath1.joinpath(f"{lname1}.{ligand_forcefield}.lib")
    frcmod1 = lpath1.joinpath(f"{lname1}.{ligand_forcefield}.frcmod")
    lpdb1 = lpath1.joinpath(f"{lname1}.pdb")
    lib2 = lpath2.joinpath(f"{lname2}.{ligand_forcefield}.lib")
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
                        resize=config["resize"],
                        ligand_cation=config["ligand_cation"],
                        ligand_anion=config["ligand_anion"],
                        complex_cation=config["complex_cation"],
                        complex_anion=config["complex_anion"]
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
        # for unified protocol
        shutil.copy(f"{system}/pressurize.rst7", f"{system}_unified.rst7")
        shutil.copy(f"{system}_vdw_bonded.parm7", f"{system}_unified.parm7")

    
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
        protocol="split",
        remd=False,
    ):

    systems = ["ligands", "complex"]

    if protocol == "split":
        steps = ["decharge", "vdw_bonded", "recharge"] 
    elif protocol == "unified":
        steps = ["unified"]

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
                tasks[f"{system}_{step}_{lmb}"]["path"] = str(lmb_dir.resolve())
                tasks[f"{system}_{step}_{lmb}"]["command"] = []

                with set_directory(lmb_dir):

                    logger.info("Pre-heating ...")
                    tasks[f"{system}_{step}_{lmb}"]["command"].append(heat(
                        defname="debug",
                        prmtop=prmtop, 
                        conf=top_dir.joinpath(f"{system}_{step}.rst7"),
                        ref=top_dir.joinpath(f"{system}_{step}.rst7"),
                        nsteps=5,
                        dt=0.00005,
                        temp=5.0,
                        resstraint_wt=config["TI"]["heat"]["resstraint_wt"],
                        fep=True,
                        clambda=lmb,
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        aces=config[step].get("aces", 0),
                        aces_setting=config.get("aces_setting", None),
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        crgmask=config[step].get("crgmask", None),
                        tempi=4.0,
                        ofreq=1,
                        fname="debug.in",
                        run=False
                    ))


                    tasks[f"{system}_{step}_{lmb}"]["command"].append( em(
                        "min",
                        prmtop=prmtop, 
                        conf="debug.rst7",
                        ref="debug.rst7",
                        maxcyc=config["TI"]["em"]["maxcyc"],
                        resstraint_wt=config["TI"]["em"]["resstraint_wt"],
                        fep=True,
                        clambda=lmb,
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        aces=config[step].get("aces", 0),
                        aces_setting=config.get("aces_setting", None),
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        crgmask=config[step].get("crgmask", None),
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
                        aces=config[step].get("aces", 0),
                        aces_setting=config.get("aces_setting", None),
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        crgmask=config[step].get("crgmask", None),
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
                        aces=config[step].get("aces", 0),
                        aces_setting=config.get("aces_setting", None),
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        crgmask=config[step].get("crgmask", None),
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
                        aces=config[step].get("aces", 0),
                        aces_setting=config.get("aces_setting", None),
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        crgmask=config[step].get("crgmask", None),
                        tempi=5.0,
                        ofreq=config["TI"]["heat"]["ofreq"],
                        fname="heat.in",
                        run=False
                    ))

                    logger.info("Pre-Pressurising ...")
                    tasks[f"{system}_{step}_{lmb}"]["command"].append( pressurize(
                        defname="pre_pressing1",
                        prmtop=prmtop, 
                        conf="heat.rst7",
                        ref="heat.rst7",
                        nsteps=5000,
                        dt=0.002,
                        temp=config['TI']['temp'],
                        resstraint_wt=config['TI']['pressurize_res']['resstraint_wt'],
                        irest=1, ntx=5,
                        fep=True,
                        clambda=lmb,
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        aces=config[step].get("aces", 0),
                        aces_setting=config.get("aces_setting", None),
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        crgmask=config[step].get("crgmask", None),
                        ofreq=100,
                        fname="pre_pressing1.in",
                        run=False
                    ))

                    logger.info("Pre-Pressurising2 ...")
                    tasks[f"{system}_{step}_{lmb}"]["command"].append( pressurize(
                        defname="pre_pressing2",
                        prmtop=prmtop, 
                        conf="pre_pressing1.rst7",
                        ref="pre_pressing1.rst7",
                        nsteps=5000,
                        dt=0.002,
                        temp=config['TI']['temp'],
                        resstraint_wt=config['TI']['pressurize_res']['resstraint_wt'],
                        irest=1, ntx=5,
                        fep=True,
                        clambda=lmb,
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        aces=config[step].get("aces", 0),
                        aces_setting=config.get("aces_setting", None),
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        crgmask=config[step].get("crgmask", None),
                        ofreq=100,
                        fname="pre_pressing2.in",
                        run=False
                    ))

                    logger.info("Pre-Pressurising3 ...")
                    tasks[f"{system}_{step}_{lmb}"]["command"].append( pressurize(
                        defname="pre_pressing3",
                        prmtop=prmtop, 
                        conf="pre_pressing2.rst7",
                        ref="pre_pressing2.rst7",
                        nsteps=5000,
                        dt=0.002,
                        temp=config['TI']['temp'],
                        resstraint_wt=config['TI']['pressurize_res']['resstraint_wt']*0.6,
                        irest=1, ntx=5,
                        fep=True,
                        clambda=lmb,
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        aces=config[step].get("aces", 0),
                        aces_setting=config.get("aces_setting", None),
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        crgmask=config[step].get("crgmask", None),
                        ofreq=100,
                        fname="pre_pressing3.in",
                        run=False
                    ))

                    logger.info("Pressurising ...")
                    tasks[f"{system}_{step}_{lmb}"]["command"].append( pressurize(
                        defname="pressurize_res",
                        prmtop=prmtop, 
                        conf="pre_pressing3.rst7",
                        ref="pre_pressing3.rst7",
                        nsteps=config["TI"]["pressurize_res"]["nsteps"],
                        dt=0.002,
                        temp=config["TI"]["temp"],
                        resstraint_wt=config["TI"]["pressurize_res"]["resstraint_wt"],
                        irest=1, ntx=5,
                        fep=True,
                        clambda=lmb,
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        aces=config[step].get("aces", 0),
                        aces_setting=config.get("aces_setting", None),
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        crgmask=config[step].get("crgmask", None),
                        ofreq=config["TI"]["pressurize_res"]["ofreq"],
                        fname="pressurize_res.in",
                        run=False
                    ))

                    logger.info("Pressurising ...")
                    tasks[f"{system}_{step}_{lmb}"]["command"].append( pressurize(
                        defname="pressurize",
                        prmtop=prmtop, 
                        conf="pressurize_res.rst7",
                        ref="pressurize_res.rst7",
                        nsteps=config["TI"]["pressurize"]["nsteps"],
                        dt=0.002,
                        temp=config["TI"]["temp"],
                        resstraint_wt=None,
                        irest=1, ntx=5,
                        fep=True,
                        clambda=lmb,
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        aces=config[step].get("aces", 0),
                        aces_setting=config.get("aces_setting", None),
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        crgmask=config[step].get("crgmask", None),
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
                        resstraint_wt=None,
                        irest=1, ntx=5,
                        iremd=1 if config['remd'] else 0,
                        numexchg=config["TI"]["production"]['numexchg'],
                        gremd_acyc=1 if config['remd'] else None,
                        fep=True,
                        clambda=lmb,
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifmbar=1,
                        mbar_states=len(config[step]["lambdas"]),
                        mbar_lambda=print_list(config[step]["lambdas"], by=", "),
                        ifsc=config[step]["ifsc"],
                        aces=config[step].get("aces", 0),
                        aces_setting=config.get("aces_setting", None),
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step].get("scmask1", None),
                        scmask2=config[step].get("scmask2", None),
                        crgmask=config[step].get("crgmask", None),
                        ntwe=config["TI"]["production"]["ntwe"],
                        ntwx=config["TI"]["production"]["ntwx"],
                        ntpr=config["TI"]["production"]["ntpr"],
                        ntwr=config["TI"]["production"]["ntwr"],
                        fname="prod.in",
                        run=False
                    ))
                    with open("run.sh", "w") as fp:
                        fp.write("\n".join(tasks[f"{system}_{step}_{lmb}"]["command"]))
                        fp.write("\n")
                        fp.write("touch done.tag")
                        fp.write("\n")
    if remd:
        # write the groupfile
        n_task = len(tasks[f"{systems[0]}_{steps[0]}_{config[step]['lambdas'][0]}"]["command"])
        step = steps[0]
        run_files = []
        for system in systems:
            task_list = []
            system_dir = cwd.joinpath(system)
            groupfile_dir = system_dir.joinpath("groupfile")
            output_dir = system_dir.joinpath("output")
            for idir in [groupfile_dir, output_dir]: 
                if not os.path.exists(idir):
                    os.mkdir(idir)
            
            for task in range(n_task):
                if task in [1, 3]:
                    for lmd in config[step]['lambdas']:
                        old_command = tasks[f"{system}_{step}_{lmd}"]["command"][task]
                        new_command_line = [config['serial_execute']] + [
                            ii if "-" in ii else str(Path(tasks[f"{system}_{step}_{lmd}"]["path"]).joinpath(ii)) for ii in old_command.split()[1:]
                        ]
                        task_list.append(f"{' '.join(new_command_line)}")
                else:
                    task_file = groupfile_dir.joinpath(f"task.{task}.groupfile")
                    with open(task_file, "w") as fp:
                        for lmd in config[step]['lambdas']:
                            index_i = tasks[f"{system}_{step}_{lmd}"]["command"][task].index("-i")
                            old_command = tasks[f"{system}_{step}_{lmd}"]["command"][task][index_i:]
                            new_command_line = [
                                ii if "-" in ii else str(Path(tasks[f"{system}_{step}_{lmd}"]["path"]).joinpath(ii)) for ii in old_command.split() 
                            ]
                            fp.write(f"{' '.join(new_command_line)}\n")

                    n_lambdas = len(config[step]['lambdas'])
                    if task < n_task - 1:
                        command = f"mpirun  --oversubscribe -np {config['np']} {config['mpi_execute']} -ng {config['ng']} -groupfile {str(task_file)}"
                    else:
                        command = f"mpirun  --oversubscribe -np {config['np']} {config['mpi_execute']} -rem 3" \
                                    f" -remlog {str(output_dir)}/remt.{task}.log -ng {config['ng']} -groupfile {str(task_file)}"
                    task_list.append(command)

            with open(groupfile_dir.joinpath(f"run.sh"), "w") as fp:
                fp.write("\n".join(task_list))
                # fp.write("\ntouch done.tag\n")
            run_files.append(groupfile_dir.joinpath(f"run.sh"))
        return run_files

    else:
        with open("tasks.json", "w") as fp:
            json.dump(tasks, fp, indent=2)

        return tasks


def submit_jobs(tasks, env, ngroup=-2, submit=False, n_lambda=-1):
    cwd = Path(os.getcwd())
    submit_scripts_dir = cwd.joinpath("submit")
    if not os.path.exists(submit_scripts_dir):
        os.mkdir(submit_scripts_dir)

    scripts = {}
    script_head = ["#!/bin/bash"] + env + ["set -e"]
    if ngroup == -2:
        logger.info("REMD is enabled.")
        for idx, task in enumerate(tasks):
            script = [f"cd {task.parent}"]
            script += [f"[ ! -f done.tag ] && bash run.sh && touch done.tag"]
            sub = submit_scripts_dir.joinpath(f"unified.{idx}.sh")
            with open(sub, "w") as fp:
                fp.write("\n".join(script_head))
                fp.write("\n")
                fp.write("\n".join(script))
        return 

    for task, info in tasks.items():
        script = [f"cd {info['path']}"]
        script += info["command"]
        script += ["touch done.tag"]
        scripts[task] = script
    
    if ngroup > 0:
        group_size = len(scripts) // ngroup
        keys = list(scripts.keys())
        random.shuffle(keys)
        for igroup in range(ngroup):
            if igroup < len(scripts) % ngroup:
                igroup_size = group_size + 1
            else:
                igroup_size = group_size

            sub = submit_scripts_dir.joinpath(f"grouped.{igroup}.sh")
            with open(sub, "w") as fp:
                fp.write("\n".join(script_head))
                fp.write("\n")
                for _ in range(igroup_size):
                    # if not os.path.exists(Path(tasks[task]['path']).joinpath("done.tag")):
                    fp.write("\n".join(scripts[keys.pop()]))
                    fp.write("\n")
            if submit:
                with set_directory(submit_scripts_dir):
                    os.system(f"LLsub {sub.name} -s 6 -g volta:1")
        assert len(keys) == 0
    else:
        for task, script in scripts.items():
            sub = submit_scripts_dir.joinpath(f"{task}.sh")
            with open(sub, "w") as fp:
                fp.write("\n".join(script_head))
                fp.write("\n")
                fp.write("\n".join(script))

            if submit:
                if not os.path.exists(Path(tasks[task]['path']).joinpath("done.tag")):
                    with set_directory(submit_scripts_dir):
                        os.system(f"LLsub {sub.name} -s 8 -g volta:1")


def run(
        ppdb, 
        lpath1, lname1, 
        lpath2, lname2,
        config,
        submit=False,
        ngroup=8,
        fep_dir_name="fep",
        overwrite=True
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

    fep_run_dir = cwd.joinpath(fep_dir_name)
    if overwrite:  # not os.path.exists(fep_run_dir.joinpath(tag)) or 
        with set_directory(fep_run_dir):
            tasks = setTI(config, top_dir=prep_dir.resolve(), protocol=config['protocol'], remd=config['remd'])
            if config['remd']:
                ngroup = -2
            submit_jobs(tasks, config['slurm_env'], submit=submit, ngroup=ngroup, n_lambda=len(config['unified']['lambdas']))
            set_tag(tag)
    else:
        logger.info("feprun has already existed. skip through it.")

