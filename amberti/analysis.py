import alchemlyb
from alchemlyb.estimators import TI, MBAR
from alchemlyb.preprocessing import decorrelate_dhdl
from alchemlyb.parsing.amber import extract_dHdl, extract_u_nk
from alchemlyb.convergence import fwdrev_cumavg_Rc
from alchemlyb.visualisation import plot_convergence
from alchemlyb.convergence import forward_backward_convergence
from alchemlyb.visualisation import plot_mbar_overlap_matrix
import os
from pathlib import Path
import json
import matplotlib.pyplot as plt


def analyze(workpath, output_dir, filename, temp=300, mode='mbar'):
    base_dir = Path(workpath)
    output_dir = Path(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    result = {}
    for cl in ["ligands", "complex"]:
        wkdir = base_dir.joinpath(cl)
        for job in [ "decharge", "vdw_bonded", "recharge"]:
            jobdir = wkdir.joinpath(job)
            lmd_lists = list(jobdir.glob("[0-9]*.*[0-9]"))
            lmd_lists.sort()
            data_list = []
            if mode == 'mbar':
                for lmd in lmd_lists:
                    print(lmd)
                    data2 = extract_u_nk(lmd.joinpath(filename), T=temp)
                    decorrelated_dhdl = decorrelate_dhdl(data2, remove_burnin=True)
                    # data_list.append(data2)
                    data_list.append(decorrelated_dhdl)
                u_nk = alchemlyb.concat(data_list)
                u_nk = u_nk.sort_index(level=u_nk.index.names[1:])

                mbar = MBAR()
                mbar.fit(u_nk)
                fe = mbar.delta_f_
                ax2 = plot_mbar_overlap_matrix(mbar.overlap_matrix)
                ax2.figure.savefig(output_dir.joinpath(f'overlap_{cl}_{job}.png'), bbox_inches='tight', pad_inches=0.0)
            else:
                raise RuntimeError()
            result[f"{cl}-{job}"] = fe.loc[0.00, 1.00]
            # R_c, running_average = fwdrev_cumavg_Rc(data_list, tol=2)
            df = forward_backward_convergence(data_list)
            ax = plot_convergence(df)
            ax.figure.savefig(output_dir.joinpath(f'dF_t_{cl}_{job}.png'))
            

    ddG = result["complex-vdw_bonded"] + result["complex-recharge"] + result["complex-decharge"] - \
        (result["ligands-vdw_bonded"] + result["ligands-recharge"] + result["ligands-decharge"])

    result["ddG"] = ddG

    with open(output_dir.joinpath('result.json'), "w") as fp:
        json.dump(result, fp, indent=4)
    
    return result
    