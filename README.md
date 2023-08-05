# Python Package to run Amber TI

## Release Update

### July 31th

* Ligands alignment from SMILES.

### July 24th

* Add supports for H-REMD.
* Reproduce ACES technique.

## Example

### Make topology

```python
from amberti.amber import make_ligand_topology
from pathlib import Path
from amberti.utils import set_directory
import os

lig_path = Path("test_toy/lig")
with set_directory(lig_path):
    make_ligand_topology(
    "M7G.sdf", 
    ncharge=-2, 
    name="M7G", 
    charge_method="bcc", 
    forcefield="gaff2", 
    mol_type="sdf", 
    fo="mol2"
    )
```

### Run TI simulation

```python

import json
from pathlib import Path
from amberti.utils import set_directory
from amberti.workflow import run
import os

ppdb = cwd.joinpath("1ipb_protein_only.pdb")
lpath1 = Path("/lig/M7G")
lname1 = "M7G"
lpath2 = Path("/lig/OME")
lname2 = "OME"


with open("config.json", "r") as fp:
    data = json.load(fp)


workdir = cwd.joinpath("test_run")
if not os.path.exists(workdir):
    os.mkdir(workdir)

with set_directory(workdir):
    run(
        ppdb, 
        lpath1, lname1, 
        lpath2, lname2,
        data
        )
```

### An example of `config.json` file

By changing the TI/softcore mask part, this JSON should run well.

```JSON
{
    "ligand_forcefield": "gaff2",
    "protein_forcefield": "ff14SB",
    "ligand_box_size": 15.0,
    "complex_box_size": 12.0,
    "resize": 0.75,
    "prep": {
        "temp": 300,
        "timask1": ":M7G",
        "timask2": ":BNC",
        "scmask1": ":M7G@H26",
        "scmask2": ":BNC@C1,C4,C9,C12,C17,C22,H27,H28,H29,H30,H31",
        "em": {
            "maxcyc": 2000,
            "resstraint_wt": 5.0
        },
        "heat": {
            "nsteps": 5000,
            "resstraint_wt": 5.0,
            "ofreq": 1000
        },
        "pressurize": {
            "nsteps": 5000,
            "resstraint_wt": 5.0,
            "ofreq": 1000
        }
    },
    "TI": {
        "temp": 300,
        "em": {
            "maxcyc": 5000,
            "resstraint_wt": 5.0
        },
        "heat": {
            "nsteps": 50000,
            "resstraint_wt": 5.0,
            "ofreq": 2500
        },
        "pressurize": {
            "nsteps": 500000,
            "resstraint_wt": 5.0,
            "ofreq": 2500
        },
        "production": {
            "nsteps": 3000000,
            "resstraint_wt": 5.0,
            "ntwe": 3000,
            "ntwx": 3000,
            "ntwr": 3000,
            "ntpr": 30000
        }
    },
    "decharge": {
        "lambdas": [
            0.0,
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.6,
            0.7,
            0.8,
            0.9,
            1.0
        ],
        "scalpha": 0.5,
        "scbeta": 12.0,
        "ifsc": 0,
        "timask1": ":1",
        "timask2": ":2",
        "crgmask": ":2@H26"
    },
    "vdw_bonded": {
        "lambdas": [
            0.0,
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.6,
            0.7,
            0.8,
            0.9,
            1.0
        ],
        "scalpha": 0.5,
        "scbeta": 12.0,
        "ifsc": 1,
        "timask1": ":1",
        "timask2": ":2",
        "scmask1": ":1@H26",
        "scmask2": ":2@C1,C4,C9,C12,C17,C22,H27,H28,H29,H30,H31",
        "crgmask": ":1@H26 | :2@C1,C4,C9,C12,C17,C22,H27,H28,H29,H30,H31"
    },
    "recharge": {
        "lambdas": [
            0.0,
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.6,
            0.7,
            0.8,
            0.9,
            1.0
        ],
        "scalpha": 0.5,
        "scbeta": 12.0,
        "ifsc": 0,
        "timask1": ":1",
        "timask2": ":2",
        "crgmask": ":1@C1,C4,C9,C12,C17,C22,H27,H28,H29,H30,H31"
    },
    "slurm_env": [
        "source /etc/profile",
        "module load cuda/11.2",
        "source $HOME/env/amber22.env"
    ]
}

```

## Note & Before you start

* Set up the simulation parameters in `config.json`. Be careful about the simulation steps as they may affect the performance significantly.

* Check your system trajectory out that the system has been well equilibrated.

## Trouble Shoot
### Input Error

* Don't use decimal as the input for some parameters. For example, `ifsc=0` works but `ifsc=0.0` doesn't work. The error message is like `Fortran runtime error: Cannot match namelist object name .0crgmask.`

* `md.in` file does end with an extra blank line.

### Simulation Error

## ERROR: Calculation halted.  Periodic box dimensions have changed too much from their initial values.

This may happen along with the following message:
```text
 Your system density has likely changed by a large amount, probably from
  starting the simulation from a structure a long way from equilibrium.

  [Although this error can also occur if the simulation has blown up for some reason]

  The GPU code does not automatically reorganize grid cells and thus you
  will need to restart the calculation from the previous restart file.
  This will generate new grid cells and allow the calculation to continue.
  It may be necessary to repeat this restarting multiple times if your system
  is a long way from an equilibrated density.

  Alternatively you can run with the CPU code until the density has converged
  and then switch back to the GPU code.
```

A simple solution is to gradually pressurize your system in multiple samll steps. When approching to the target pressure, run a relatively long simulation (like $1 \ \mathrm{ns}$) to further equilibrate the system.


## 
