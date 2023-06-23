from amberti.utils import run_command
from amberti.logger import getLogger

logger = getLogger()

def antechamber(mol, ncharge, name=None, charge_method="bcc", forcefield="gaff", mol_type="sdf", fo="mol2"):
    if name is None:
        mol_name = ".".join(mol.split(".")[:-1])
        name = mol_name[:3]
    command = [
        "antechamber", 
        "-i", mol, "-fi", mol_type,
        "-o", f"{name}.{forcefield}.mol2", "-fo", fo, 
        "-rn", name, 
        "-at", forcefield, 
        "-c", charge_method, "-nc", ncharge,
        "-an", "yes", "-dr", "no", "-pf", "yes", 
    ]
    # antechamber -i M7G.sdf -fi sdf -o M7G.gaff.mol2 -fo mol2 -rn M7G -at gaff -an yes -dr no -pf yes -c bcc -nc -2
    logger.info(f'Running "{" ".join(command)}"')
    run_command(command)
    

def parmchk2(mol, forcefield="gaff", mol_type="mol2", fo="mol2"):
    # parmchk2 -i M7G.gaff.mol2 -f mol2 -o M7G.gaff.frcmod -s gaff -a yes
    mol_name = ".".join(mol.split(".")[:-1])
    command = [
        "parmchk2", 
        "-i", mol, "-f", mol_type,
        "-o", f"{mol_name}.frcmod"
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
        f"MOL = loadmol2 {mol2} ",
        f"saveoff MOL {name}.lib",
        f"savepdb MOL {name}.pdb",
        "quit"
    ]
    tleap_in_file = f"tleap.{name}.in"
    with open(tleap_in_file, "w") as fp:
        fp.writelines(tleap_in)
    run_command(["tleap", "-f", tleap_in_file])


    """tleap -f - <<_EOF

quit
_EOF"""
    