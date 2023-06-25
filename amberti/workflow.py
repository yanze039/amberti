import os
from pathlib import Path
from amberti.md import em, heat, pressurize, production
from amberti.logger import getLogger
from amberti.utils import set_directory


logger = getLogger()


def prep(
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
                conf=top_dir.joinpath(f"{system}_vdw_bonded.parm7"),
                ref=top_dir.joinpath(f"{system}_vdw_bonded.parm7"),
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
                run=False
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
                run=False
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
                run=False
            )

            logger.info(f"Heating from 5K to {config['prep']['temp']}K ...")
            heat(
                defname="heating",
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
                run=False
            )

            logger.info("Pre-Pressurising ...")
            pressurize(
                defname="pre_pressing",
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
                run=False
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
                run=False
            )

            logger.info("Equilibrium finished.")

    return


def setTI(
        config,
        top_dir,
        lambdas: dict
    ):
    """
        To shedule the lambda list, a dict is desired here.
        lambdas =  {
            "decharge": [0.0, 0.1], 
            "vdw_transfrom": [0.0, 1.0],
            "recharge": [0.0, 1.0]
        }
    """
    systems = ["ligands", "complex"]
    steps = ["decharge", "vdw_bonded", "recharge"]
    top_dir = Path(top_dir)
    
    cwd = Path(os.getcwd())
    for system in systems:
        system_dir = cwd.joinpath(system)
        if not os.path.exists(system_dir):
            os.mkdir(system_dir)

        for step in steps:
            step_dir = system_dir.joinpath(step)
            if not os.path.exists(step_dir):
                os.mkdir(step_dir)
            
            for lmb in lambdas[step]:
                lmb_dir = step_dir.joinpath(lmb)
                if not os.path.exists(lmb_dir):
                    os.mkdir(lmb_dir)
                prmtop = top_dir.joinpath(f"{system}_{step}.parm7")
                with set_directory(lmb_dir):
                    em(
                        "min",
                        prmtop=prmtop, 
                        conf=top_dir.joinpath(f"{system}_{step}.rst7"),
                        ref=top_dir.joinpath(f"{system}_{step}.rst7"),
                        maxcyc=config["TI"]["em"]["maxcyc"],
                        resstraint_wt=config["TI"]["em"]["resstraint_wt"],
                        fep=True,
                        clambda=config[step]["clambda"],
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step]["scmask1"],
                        scmask2=config[step]["scmask2"],
                        fname="min.in",
                        run=False
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
                        resstraint_wt=config["TI"]["heat"]["resstraint_wt"],
                        fep=True,
                        clambda=config[step]["clambda"],
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step]["scmask1"],
                        scmask2=config[step]["scmask2"],
                        tempi=4.0,
                        ofreq=1,
                        fname="pre_heat.in",
                        run=False
                    )

                    logger.info("Second Emergy minimizing ...")
                    em(
                        "min2",
                        prmtop=prmtop,
                        conf="pre_heating.rst7",
                        ref="pre_heating.rst7",
                        maxcyc=config["TI"]["em"]["maxcyc"],
                        resstraint_wt=config["TI"]["em"]["resstraint_wt"],
                        fep=True,
                        clambda=config[step]["clambda"],
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step]["scmask1"],
                        scmask2=config[step]["scmask2"],
                        fname="min2.in",
                        run=False
                    )

                    logger.info(f"Heating from 5K to {config['TI']['temp']}K ...")
                    heat(
                        defname="heating",
                        prmtop=prmtop, 
                        conf="min2.rst7",
                        ref="min2.rst7",
                        nsteps=config["TI"]["heat"]["nsteps"],
                        temp=config['TI']['temp'],
                        resstraint_wt=config["TI"]["heat"]["resstraint_wt"],
                        fep=True,
                        clambda=config[step]["clambda"],
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step]["scmask1"],
                        scmask2=config[step]["scmask2"],
                        tempi=5.0,
                        ofreq=config["TI"]["heat"]["ofreq"],
                        fname="heat.in",
                        run=False
                    )

                    logger.info("Pre-Pressurising ...")
                    pressurize(
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
                        clambda=config[step]["clambda"],
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step]["scmask1"],
                        scmask2=config[step]["scmask2"],
                        ofreq=10,
                        fname="pre_press.in",
                        run=False
                    )

                    logger.info("Pressurising ...")
                    pressurize(
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
                        clambda=config[step]["clambda"],
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifsc=config[step]["ifsc"],
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step]["scmask1"],
                        scmask2=config[step]["scmask2"],
                        ofreq=config["TI"]["pressurize"]["ofreq"],
                        fname="pressurize.in",
                        run=False
                    )

                    logger.info("Equilibrium finished.")

                    logger.info("Production run ...")
                    production(
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
                        clambda=config[step]["clambda"],
                        scalpha=config[step]["scalpha"],
                        scbeta=config[step]["scbeta"],
                        ifmbar=1,
                        mbar_states=len(config[step]["lambdas"]),
                        mbar_lambda=", ".join(config[step]["lambdas"]),
                        ifsc=config[step]["ifsc"],
                        timask1=config[step]["timask1"],
                        timask2=config[step]["timask2"],
                        scmask1=config[step]["scmask1"],
                        scmask2=config[step]["scmask2"],
                        ntwe=config["TI"]["production"]["ntwe"],
                        ntwx=config["TI"]["production"]["ntwx"],
                        ntpr=config["TI"]["production"]["ntpr"],
                        ntwr=config["TI"]["production"]["ntwr"],
                        fname="prod.in",
                        run=False
                    )

