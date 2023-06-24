from amberti.utils import run_command
from amberti.logger import getLogger
import os

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


def build_lib(forcefield, frcmod, mol2, name=None):
    if name is None:
        mol_name = ".".join(mol2.split(".")[:-1])
        name = mol_name[:3]
    tleap_in = [
        f"source leaprc.{forcefield}",
        f"loadamberparams {frcmod}",
        f"{name} = loadmol2 {mol2} ",
        f"saveoff {name} {name}.lib",
        f"savepdb {name} {name}.pdb",
        "quit"
    ]
    tleap_in_file = f"tleap.{name}.in"
    with open(tleap_in_file, "w") as fp:
        fp.write("\n".join(tleap_in))
    command = ["tleap", "-f", tleap_in_file]
    logger.info(f'Running "{" ".join(command)}"')
    run_command(command)


def maketop(
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
    frcmod = f"{name}.frcmod"
    lib = f"{name}.lib"
    pdb = f"{name}.pdb"
    antechamber(
        mol, ncharge, output=mol2, name=name, charge_method=charge_method, 
        forcefield=forcefield, mol_type=mol_type, fo=fo
    )
    parmchk2(mol2, frcmod, mol_type=fo, forcefield=forcefield)
    build_lib(forcefield, frcmod, mol2, name=name)

    # check validity
    assert os.path.exists(mol2)
    assert os.path.exists(lib)
    assert os.path.exists(frcmod)
    assert os.path.exists(pdb)

    logger.info(f'{mol2}, {frcmod} are created successfully.')
    