from amberti.utils import run_command
from amberti.logger import getLogger
from alchemlyb.postprocessors.units import R_kJmol, kJ2kcal
from alchemlyb.parsing.amber import SectionParser
import pandas as pd
import numpy as np
import os


k_b = R_kJmol * kJ2kcal

logger = getLogger()

def antechamber(mol, ncharge, output, name=None, charge_method="bcc", forcefield="gaff", mol_type="sdf", fo="mol2"):
    if name is None:
        mol_name = ".".join(mol.split(".")[:-1])
        name = mol_name[:3]
    command = [
        "antechamber", 
        "-i", mol, "-fi", mol_type,
        "-o", output, "-fo", fo, 
        "-rn", name, 
        "-at", forcefield, 
        "-c", charge_method, "-nc", str(ncharge),
        "-an", "yes", "-dr", "no", "-pf", "yes", 
    ]
    # antechamber -i M7G.sdf -fi sdf -o M7G.gaff.mol2 -fo mol2 -rn M7G -at gaff -an yes -dr no -pf yes -c bcc -nc -2
    logger.info(f'Running "{" ".join(command)}"')
    run_command(command)
    

def parmchk2(mol, output, forcefield="gaff", mol_type="mol2"):
    # parmchk2 -i M7G.gaff.mol2 -f mol2 -o M7G.gaff.frcmod -s gaff -a yes
    mol_name = ".".join(mol.split(".")[:-1])
    command = [
        "parmchk2", 
        "-i", mol, "-f", mol_type,
        "-o", output,
        "-s", forcefield, "-a", "yes"
    ]
    logger.info(f'Running "{" ".join(command)}"')
    run_command(command)


def tleap(script, fname="tleap.in"):
    with open(fname, "w") as fp:
        fp.write(script)
    command = ["tleap", "-f", fname]
    logger.info(f'Running "{" ".join(command)}"')
    return_code, _, _ = run_command(command)
    assert return_code == 0


def cpptraj(script, prmtop, fname="cpptraj.in"):
    with open(fname, "w") as fp:
        fp.write(script)
    command = ["cpptraj", "-p", prmtop, fname]
    logger.info(f'Running "{" ".join(command)}"')
    run_command(command)
    return


def pmemd(defname, md_in, prmtop, conf, outtraj=True, ref=None, cuda=True, mpi=False, run=True):

    if cuda:
        pmemd_binary = "pmemd.cuda"
    elif mpi:
        pmemd_binary = "pmemd.mpi"
    else:
        pmemd_binary = "pmemd"

    command = [
        pmemd_binary,
        "-i", str(md_in), 
        "-p", str(prmtop),
        "-c", str(conf),
    ]

    if ref is not None:
        command += ["-ref", str(ref)]
    
    command += [ "-O", 
        "-o", f"{defname}.out",
        "-e", f"{defname}.en",
        "-inf", f"{defname}.info",
        "-r", f"{defname}.rst7",
        "-l", f"{defname}.log"
    ]

    if outtraj:
        command += [
            "-x", f"{defname}.nc"
        ]
    
    if run:
        logger.info(f'Running "{" ".join(command)}"')
        return_code, _, _ = run_command(command)
        assert return_code == 0
        return return_code
    else:
        return " ".join(command)
    


def make_ligand_topology(
          mol, 
          ncharge, 
          name=None, 
          charge_method="bcc", 
          forcefield="gaff", 
          mol_type="sdf", 
          fo="mol2"
    ):
    if name is None:
        mol_name = ".".join(mol.split(".")[:-1])
        name = mol_name[:3]
    
    # these files are supposed to be created.
    mol2 = f"{name}.{forcefield}.{fo}"
    frcmod = f"{name}.{forcefield}.frcmod"
    lib = f"{name}.{forcefield}.lib"
    pdb = f"{name}.pdb"
    antechamber(
        mol, ncharge, output=mol2, name=name, charge_method=charge_method, 
        forcefield=forcefield, mol_type=mol_type, fo=fo
    )
    parmchk2(mol2, frcmod, mol_type=fo, forcefield=forcefield)

    tleap_in = [
        f"source leaprc.{forcefield}",
        f"loadamberparams {frcmod}",
        f"{name} = loadmol2 {mol2} ",
        f"saveoff {name} {lib}",
        f"savepdb {name} {pdb}",
        "quit"
    ]
    tleap_in_file = f"tleap.{name}.in"
    tleap("\n".join(tleap_in), fname=tleap_in_file)

    # check validity
    assert os.path.exists(mol2)
    assert os.path.exists(lib)
    assert os.path.exists(frcmod)
    assert os.path.exists(pdb)

    logger.info(f'{mol2}, {frcmod} are created successfully.')



def extract_unk(outfile, T, mbar_lambdas):

    beta = 1 / (k_b * T)

    finished = False

    all_mbar_energy = []
    with SectionParser(outfile) as secp:
        line = secp.skip_lines(5)

        for line in secp:
            if "      A V E R A G E S   O V E R" in line:
                _ = secp.skip_after("^|=========================================")
            elif line.startswith("MBAR Energy analysis"):

                mbar = secp.extract_section(
                    "^MBAR", "^ ---", mbar_lambdas, extra=line
                )
                all_mbar_energy.append(mbar)


                if None in mbar:
                    msg = "Something strange parsing the following MBAR section."
                    msg += "\nMaybe the mbar_lambda values are incorrect?"
                    logger.error("{}\n{}", msg, mbar)
                    raise ValueError(msg)
                
                # print(mbar)

        #         reference_energy = mbar[file_datum.mbar_lambda_idx]
        #         for lmbda, energy in enumerate(mbar):
        #             if energy > 0.0:
        #                 high_E_cnt += 1

        #             file_datum.mbar_energies[lmbda].append(
        #                 beta * (energy - reference_energy)
        #             )
            elif line == "   5.  TIMINGS\n":
                finished = True

        # mbar_df = pd.DataFrame(
        #     file_datum.mbar_energies,
        #     index=np.array(file_datum.mbar_lambdas, dtype=np.float64),
        #     columns=pd.MultiIndex.from_arrays(
        #         [mbar_time, np.repeat(file_datum.clambda, len(mbar_time))],
        #         names=["time", "lambdas"],
        #     ),
        # ).T


    return np.array(all_mbar_energy)

