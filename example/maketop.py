import sys,os
kit_path = os.path.abspath("/home/gridsan/ywang3/Project/amberti/")
sys.path.insert(0, kit_path)
import amberti
from amberti.amber import make_ligand_topology
from pathlib import Path
from amberti.utils import set_directory
import os
import argparse

parser = argparse.ArgumentParser(
        description="AmberTI",
        # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n"
    )

parser.add_argument("MOL")
parser.add_argument("NAME")
parsed_args = parser.parse_args()

make_ligand_topology(
    parsed_args.MOL, 
    ncharge=-2, 
    name=parsed_args.NAME, 
    charge_method="bcc", 
    forcefield="gaff2",
    mol_type="sdf", 
    fo="mol2"
    )
