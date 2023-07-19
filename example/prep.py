import sys,os
kit_path = os.path.abspath("../amberti/")
sys.path.insert(0, kit_path)
import amberti
from amberti.amber import make_ligand_topology
from pathlib import Path
from amberti.utils import set_directory
import json
from amberti.workflow import run
import os
import argparse

parser = argparse.ArgumentParser(
        description="AmberTI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n"
    )

parser.add_argument("--lpath1")
parser.add_argument("--lpath2")
parser.add_argument("--lname1")
parser.add_argument("--lname2")
parser.add_argument("--ppdb")
parser.add_argument("--config")
parser.add_argument("--workdir")
parsed_args = parser.parse_args()


ppdb = parsed_args.ppdb
lpath1 = parsed_args.lpath1
lname1 = parsed_args.lname1
lpath2 = parsed_args.lpath2
lname2 = parsed_args.lname2


with open(parsed_args.config, "r") as fp:
    data = json.load(fp)


workdir = Path(parsed_args.workdir)

if not os.path.exists(workdir):
    os.mkdir(workdir)

with set_directory(workdir):
    run(
        ppdb, 
        lpath1, lname1, 
        lpath2, lname2,
        data
        )

