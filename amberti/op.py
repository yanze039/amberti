from amberti.logger import getLogger
from amberti.amber import tleap
import parmed as pmd
import pytraj as pt
from pathlib import Path


logger = getLogger()


def create_complex_box(
                        protein,
                        protein_name,
                        ligand,
                        ligand_name,
                        lib,
                        frcmod,
                        protein_forcefield,
                        ligand_forcefield,
                        water='tip3p',
                        box_size=15.0,
                        ion_num=12
    ):
    """Create a simulation box with protein and water."""
    if not isinstance(protein, Path):
        protein = Path(protein)
    else:
        protein = protein
    
    if not isinstance(ligand, Path):
        ligand = Path(ligand)
    else:
        ligand = ligand

    system_name = f"{protein_name}_{ligand_name}"
    protein = protein.absolute()
    ligand = ligand.absolute()

    if water == "tip3p":
        waterbox = "TIP3PBOX"
        ionparm = "ionsjc_tip3p"
    else:
        logger.error("Water box style can only support tip3p. Not implemented.")
        RuntimeError("Not implemented.")
    
    scripts = [
        f"source leaprc.water.{water}",
        f"source leaprc.protein.{protein_forcefield}",
        f"source leaprc.{ligand_forcefield}",
        f"loadAmberParams frcmod.{ionparm}",
        f"loadoff {lib}",
        f"loadamberparams {frcmod}",

        # load the coordinates and create the complex
        f"protein = loadpdb {protein}",
        f"ligands = loadMol2 {ligand}",
        # create proteins in solution
        "complex = combine { protein ligands }",
        f"solvatebox complex {waterbox} {box_size}",

        # "addions ligands Na+ 0",
        "addionsrand complex Na+ 0",
        "addionsrand complex Cl- 0",
        f"addionsrand complex Na+ {ion_num}",
        f"addionsrand complex Cl- {ion_num}",
        f"savepdb complex {system_name}.tleap.pdb",
        f"saveamberparm complex {system_name}.tleap.parm7 {system_name}.tleap.rst7",
        "quit"
        
        ]
    fname = f"tleap.{system_name}.in"
    tleap("\n".join(scripts), fname=fname)


def create_ligand_box(
                        ligand,
                        ligand_name,
                        lib,
                        frcmod,
                        ligand_forcefield,
                        water='tip3p',
                        box_size=15.0,
                        ion_num=12
            ):
    """Create a simulation box with protein and water."""
    
    if not isinstance(ligand, Path):
        ligand = Path(ligand)
    else:
        ligand = ligand

    system_name = ligand_name
    ligand = ligand.absolute()

    if water == "tip3p":
        waterbox = "TIP3PBOX"
        ionparm = "ionsjc_tip3p"
    else:
        logger.error("Water box style can only support tip3p. Not implemented.")
        RuntimeError("Not implemented.")
    
    scripts = [
        f"source leaprc.water.{water}",
        f"source leaprc.{ligand_forcefield}",
        f"loadAmberParams frcmod.{ionparm}",
        f"loadoff {lib}",
        f"loadamberparams {frcmod}",

        # load the coordinates and create the complex
        f"ligands = loadMol2 {ligand}",
        # create proteins in solution
        f"solvatebox ligands {waterbox} {box_size}",

        # "addions ligands Na+ 0",
        "addionsrand ligands Na+ 0",
        "addionsrand ligands Cl- 0",
        f"addionsrand ligands Na+ {ion_num}",
        f"addionsrand ligands Cl- {ion_num}",
        f"savepdb ligands {system_name}.tleap.pdb",
        f"saveamberparm ligands {system_name}.tleap.parm7 {system_name}.tleap.rst7",
        "quit"
        ]
    fname = f"tleap.{system_name}.in"
    tleap("\n".join(scripts), fname=fname)


def create_simulation_box(
                          lib1, lib2, 
                          frcmod1, frcmod2, 
                          lpdb1, lpdb2, ppdb, llpdb,
                          ligand_forcefield,
                          protein_forcefield,
                          water='tip3p',
                          size_ligand=15.0, size_complex=12.0,
                          resize=0.75,
                          ligand_cation=12, 
                          ligand_anion=10,
                          complex_cation=51, 
                          complex_anion=50,
    ):
    if water == "tip3p":
        waterbox = "TIP3PBOX"
        ionparm = "ionsjc_tip3p"
    else:
        logger.error("Water box style can only support tip3p. Not implemented.")
        RuntimeError("Not implemented.")

    scripts = [
        f"source leaprc.{ligand_forcefield}",
        f"source leaprc.water.{water}",
        f"source leaprc.protein.{protein_forcefield}",
        f"loadAmberParams frcmod.{ionparm}",
        f"loadoff {lib1}",
        f"loadoff {lib2}",
        f"loadamberparams {frcmod1}",
        f"loadamberparams {frcmod2}",

        # load the coordinates and create the complex
        f"mol1 = loadpdb {lpdb1}",
        f"mol2 = loadpdb {lpdb2}",
        f"protein = loadpdb {ppdb}",
        f"ligands = loadpdb {llpdb}",
        "complex1 = combine {mol1 protein}",
        "complex2 = combine {mol2 protein}",
        "complex3 = combine {ligands protein}",
        
        # create ligands in solution for vdw+bonded transformation
        f"solvatebox ligands {waterbox} {size_ligand} {resize}",
        f"addionsrand ligands Na+ {ligand_cation} Cl- {ligand_anion}",
        "savepdb ligands ligands_vdw_bonded.pdb",
        "saveamberparm ligands ligands_vdw_bonded.parm7 ligands_vdw_bonded.rst7",

        # create complex in solution for vdw+bonded transformation
        f"solvatebox complex3 {waterbox} {size_complex} {resize}",
        f"addionsrand complex3 Na+ {complex_cation} Cl- {complex_anion}",
        "savepdb complex3 complex_vdw_bonded.pdb",
        "saveamberparm complex3 complex_vdw_bonded.parm7 complex_vdw_bonded.rst7",
        "quit"
        
        ]
    fname = "tleap.buildtop.in"
    tleap("\n".join(scripts), fname=fname)



def make_charge_transform(
            lib1, lib2, 
            frcmod1, frcmod2, 
            lsolv, lmol1, lmol2,
            csolv, cmol1, cmol2,
            ligand_forcefield,
            protein_forcefield,
            water='tip3p',
    ):
    if water == "tip3p":
        ionparm = "ionsjc_tip3p"
    else:
        logger.error("Water box style can only support tip3p. Not implemented.")
        RuntimeError("Not implemented.")

    scripts = [
        f"source leaprc.{ligand_forcefield}",
        f"source leaprc.water.{water}",
        f"source leaprc.protein.{protein_forcefield}",
        f"loadAmberParams frcmod.{ionparm}",
        f"loadoff {lib1}",
        f"loadoff {lib2}",
        f"loadamberparams {frcmod1}",
        f"loadamberparams {frcmod2}",

        # load the coordinates and create the complex
        f"lsolv = loadpdb {lsolv}",
        f"lmol1 = loadpdb {lmol1}",
        f"lmol2 = loadpdb {lmol2}",

        f"csolv = loadpdb {csolv}",
        f"cmol1 = loadpdb {cmol1}",
        f"cmol2 = loadpdb {cmol2}",

        # decharge transformation
        'decharge = combine {lmol1 lmol1 lsolv}',
        'setbox decharge vdw',
        'savepdb decharge ligands_decharge.pdb',
        'saveamberparm decharge ligands_decharge.parm7 ligands_decharge.rst7',

        'decharge = combine {cmol1 cmol1 csolv}',
        'setbox decharge vdw',
        'savepdb decharge complex_decharge.pdb',
        'saveamberparm decharge complex_decharge.parm7 complex_decharge.rst7',

        # recharge transformation
        'recharge = combine {lmol2 lmol2 lsolv}',
        'setbox recharge vdw',
        'savepdb recharge ligands_recharge.pdb',
        'saveamberparm recharge ligands_recharge.parm7 ligands_recharge.rst7',

        'recharge = combine {cmol2 cmol2 csolv}',
        'setbox recharge vdw',
        'savepdb recharge complex_recharge.pdb',
        'saveamberparm recharge complex_recharge.parm7 complex_recharge.rst7',

        "quit"
        ]
    
    fname = "tleap.charge.in"
    logger.info("Create the decharging and recharging topology.")
    tleap("\n".join(scripts), fname=fname)


def concatenate_pdb(pdb1, pdb2, out):
    with open(pdb1, "r") as fp:
        mol1 = fp.read().strip().split("\n")
    with open(pdb2, "r") as fp:
        mol2 = fp.read().strip().split("\n")
    mol1 = mol1[:-1]  # remove the `END` section.
    mol = mol1 + mol2
    with open(out, "w") as fp:
        fp.write("\n".join(mol))
        fp.write("\n")


def extract_conf(rst7, prmtop, out, select, overwrite=False):
    traj = pt.iterload(rst7, prmtop)
    pt.write_traj(out, traj[str(select)], overwrite=overwrite)


def get_atom_name_from_atom_idx(mol2file, atom_idx):
    if isinstance(atom_idx, int):
        atom_idx = [atom_idx]
    mol2 = pmd.load_file(str(mol2file))
    atom_names = []
    for atom in atom_idx:
        atom_names.append(mol2[atom].name)
    return atom_names