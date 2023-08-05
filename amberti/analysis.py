import os
import re
import json
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

import alchemlyb
from alchemlyb.estimators import MBAR
from alchemlyb.preprocessing import  decorrelate_u_nk
from alchemlyb.parsing.amber import file_validation, SectionParser, convert_to_pandas
from alchemlyb.convergence import fwdrev_cumavg_Rc
from alchemlyb.visualisation import plot_convergence
from alchemlyb.convergence import forward_backward_convergence
from alchemlyb.visualisation import plot_mbar_overlap_matrix
from loguru import logger
from alchemlyb.parsing import _init_attrs_dict
from alchemlyb.parsing.util import anyopen
from alchemlyb.postprocessors.units import R_kJmol, kJ2kcal


k_b = R_kJmol * kJ2kcal

_FP_RE = r"[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?"


@_init_attrs_dict
def extract(outfile, T):

    beta = 1 / (k_b * T)

    file_datum = file_validation(outfile)

    if not np.isclose(T, file_datum.T, atol=0.01):
        msg = f"The temperature read from the input file ({file_datum.T:.2f} K)"
        msg += f" is different from the temperature passed as parameter ({T:.2f} K)"
        logger.error(msg)
        raise ValueError(msg)

    finished = False
    with SectionParser(outfile) as secp:
        line = secp.skip_lines(5)
        high_E_cnt = 0
        nensec = 0
        old_nstep = -1
        for line in secp:
            if "      A V E R A G E S   O V E R" in line:
                _ = secp.skip_after("^|=========================================")
            elif line.startswith(" NSTEP"):
                nstep, dvdl = secp.extract_section(
                    "^ NSTEP", "^ ---", ["NSTEP", "DV/DL"], extra=line
                )
                if nstep != old_nstep and dvdl is not None and nstep is not None:
                    if finished:
                        raise ValueError(
                            "TI Energy detected after the TIMINGS section. Did you concatenate the output file?"
                        )
                    file_datum.gradients.append(dvdl)
                    nensec += 1
                    old_nstep = nstep
            elif line.startswith("MBAR Energy analysis") and file_datum.have_mbar:
                if finished:
                    raise ValueError(
                        "MBAR Energy detected after the TIMINGS section. Did you concatenate the output file?"
                    )
                mbar = secp.extract_section(
                    "^MBAR", "^ ---", file_datum.mbar_lambdas, extra=line
                )

                if None in mbar:
                    msg = "Something strange parsing the following MBAR section."
                    msg += "\nMaybe the mbar_lambda values are incorrect?"
                    logger.error("{}\n{}", msg, mbar)
                    raise ValueError(msg)

                reference_energy = mbar[file_datum.mbar_lambda_idx]
                for lmbda, energy in enumerate(mbar):
                    if energy == float("inf"):
                        print("WARNING: MBAR energy is infinite")

                    if energy > 0.0:
                        high_E_cnt += 1

                    file_datum.mbar_energies[lmbda].append(
                        beta * (energy - reference_energy)
                    )
            elif line == "   5.  TIMINGS\n":
                finished = True

        if high_E_cnt:
            logger.warning(
                "{} MBAR energ{} > 0.0 kcal/mol",
                high_E_cnt,
                "ies are" if high_E_cnt > 1 else "y is",
            )

    if not finished:
        logger.warning("WARNING: file {} is a prematurely terminated run", outfile)

    if file_datum.have_mbar:
        mbar_time = [
            file_datum.t0 + (frame_index + 1) * file_datum.dt * file_datum.ntpr
            for frame_index in range(len(file_datum.mbar_energies[0]))
        ]

        mbar_df = pd.DataFrame(
            file_datum.mbar_energies,
            index=np.array(file_datum.mbar_lambdas, dtype=np.float64),
            columns=pd.MultiIndex.from_arrays(
                [mbar_time, np.repeat(file_datum.clambda, len(mbar_time))],
                names=["time", "lambdas"],
            ),
        ).T
    else:
        logger.info('WARNING: No MBAR energies found! "u_nk" entry will be None')
        mbar_df = None

    if not nensec:
        logger.warning("WARNING: File {} does not contain any dV/dl data", outfile)
        dHdl_df = None
    else:
        logger.info("Read {} dV/dl data points in file {}", nensec, outfile)
        dHdl_df = convert_to_pandas(file_datum)
        dHdl_df["dHdl"] *= beta
    
    max_value = np.nanmax(mbar_df[mbar_df != np.inf])
    mbar_df.replace([np.inf, -np.inf], max_value, inplace=True)
    return {"u_nk": mbar_df, "dHdl": dHdl_df}


def extract_u_nk(outfile, T):
    extracted = extract(outfile, T)
    return extracted["u_nk"]



def analyze(workpath, output_dir, filename, temp=300, protocol='unified', test=True, suffix=None):
    base_dir = Path(workpath)
    output_dir = Path(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    result = {}
    if protocol == 'unified':
        tasks = ["unified"]
    elif protocol == 'split':
        tasks = [ "decharge", "vdw_bonded", "recharge"]
    
    if test:
        systems = ["ligands"]
    else:
        systems = ["ligands", "complex"]

    for cl in systems:
        wkdir = base_dir.joinpath(cl)
        for job in tasks:
            if suffix is None:
                job = f"{job}"
            else:
                job = f"{job}.{suffix}"
            jobdir = wkdir.joinpath(job)
            lmd_lists = list(jobdir.glob("[0-9]*.*[0-9]"))
            lmd_lists.sort()
            data_list = []
            
            for lmd in lmd_lists:
                print(lmd)
                # data2 = extract_u_nk(lmd.joinpath(filename), T=temp)
                data2 = alchemlyb.parsing.amber.extract_u_nk(lmd.joinpath(filename), T=temp)
                # exit(0)
                # data2[data2<1e5].dropna()
                decorrelated_u_nk = decorrelate_u_nk(data2, method='dE', remove_burnin=True)
                # data_list.append(data2)
                data_list.append(decorrelated_u_nk)
            u_nk = alchemlyb.concat(data_list)
            u_nk = u_nk.sort_index(level=u_nk.index.names[1:])

            mbar = MBAR()
            mbar.fit(u_nk)
            fe = mbar.delta_f_
            ax2 = plot_mbar_overlap_matrix(mbar.overlap_matrix)
            ax2.figure.savefig(output_dir.joinpath(f'overlap_{cl}_{job}.png'), bbox_inches='tight', pad_inches=0.0)
        
            result[f"{cl}-{job}"] = fe.loc[0.00, 1.00]
            # R_c, running_average = fwdrev_cumavg_Rc(data_list, tol=2)
            df = forward_backward_convergence(data_list)
            ax = plot_convergence(df)
            ax.figure.savefig(output_dir.joinpath(f'dF_t_{cl}_{job}.png'))


    if test:
        with open(output_dir.joinpath('result.json'), "w") as fp:
            json.dump(result, fp, indent=4)
        return
    
    if protocol == 'unified':
        ddG = result["complex-unified"] - result["ligands-unified"]
    elif protocol == 'split':
        ddG = result["complex-vdw_bonded"] + result["complex-recharge"] + result["complex-decharge"] - \
            (result["ligands-vdw_bonded"] + result["ligands-recharge"] + result["ligands-decharge"])

    result["ddG"] = ddG

    with open(output_dir.joinpath('result.json'), "w") as fp:
        json.dump(result, fp, indent=4)
    
    return result
    