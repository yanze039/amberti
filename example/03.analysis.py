import sys,os
kit_path = os.path.abspath("/home/gridsan/ywang3/Project/amberti/")
sys.path.insert(0, kit_path)
import amberti
from amberti.amber import make_ligand_topology
from pathlib import Path
from amberti.utils import set_directory
import json
from amberti.analysis import analyze
import os
import argparse

parser = argparse.ArgumentParser(
        description="AmberTI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n"
    )

parser.add_argument("WORKPATH")
parser.add_argument("OUTPUT")
parser.add_argument("--filename", default="prod.out")
parser.add_argument("--mode", default="mbar")
parser.add_argument("--suffix", default=None)
parser.add_argument("--protocol", default="unified")
parsed_args = parser.parse_args()


analyze(
    parsed_args.WORKPATH,
    parsed_args.OUTPUT,
    filename=parsed_args.filename,
    suffix = parsed_args.suffix,
)

