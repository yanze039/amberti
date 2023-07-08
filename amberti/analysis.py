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


def analyze(workpath, output_dir, filename, temp=300, mode='ti'):
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
            # if mode == "ti":
            #     for lmd in lmd_lists:
            #         data1 = extract_dHdl(lmd.joinpath(filename), T=temp)
            #         # data2 = extract_u_nk(lmd.joinpath(filename), T=temp)
            #         print(data1)
            #         # decorrelated_dhdl = decorrelate_dhdl(data1, remove_burnin=True)
            #         data_list.append(data1)
            #         # data_list.append(data2)
            #     dHdl = alchemlyb.concat(data_list)
            #     # u_nk = alchemlyb.concat(data_list)
            #     ti = TI()
            #     ti.fit(dHdl)
            #     fe = ti.delta_f_
            #     # mbar = MBAR()
            #     # mbar.fit(u_nk)
            #     # fe = mbar.delta_f_
            if mode == 'mbar':
                for lmd in lmd_lists:
                    print(lmd)
                    data2 = extract_u_nk(lmd.joinpath(filename), T=temp)
                    # decorrelated_dhdl = decorrelate_dhdl(data1, remove_burnin=True)
                    data_list.append(data2)
                u_nk = alchemlyb.concat(data_list)
                u_nk = u_nk.sort_index(level=u_nk.index.names[1:])
                print(u_nk)
                exit(0)
                mbar = MBAR()
                mbar.fit(u_nk)
                fe = mbar.delta_f_
                ax2 = plot_mbar_overlap_matrix(mbar.overlap_matrix)
                ax2.figure.savefig(output_dir.joinpath(f'overlap_{cl}_{job}.png'), bbox_inches='tight', pad_inches=0.0)
            else:
                raise RuntimeError()
            result[f"{cl}-{job}"] = fe.loc[0.00, 1.00]
            # R_c, running_average = fwdrev_cumavg_Rc(data_list, tol=2)
            df = forward_backward_convergence(data_list, mode)
            ax = plot_convergence(df)
            ax.figure.savefig(output_dir.joinpath(f'dF_t_{cl}_{job}.png'))
            
            plt.figure()
            query_lambda = 0.5
            for ii in range(len(data_list)):
                print(data_list[ii][query_lambda])
                plt.hist(data_list[ii][query_lambda], alpha=0.7, label=f"{data_list[ii].index[0][1]:.1f}")
            plt.xlabel('Energy (kcal/mol)')
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='$\lambda$')
            plt.title(f'{cl}-{job}-Hamiltonian under $\lambda={query_lambda}$')
            plt.ylabel('Frequency')
            plt.tight_layout()    
            plt.ylim(0, 60)
            plt.savefig(output_dir.joinpath(f'overlap-{cl}-{job}-Hamiltonian-{query_lambda}.png'))
        break
    print(fe)
    return

    ddG = result["complex-vdw_bonded"] + result["complex-recharge"] + result["complex-decharge"] - \
        (result["ligands-vdw_bonded"] + result["ligands-recharge"] + result["ligands-decharge"])

    result["ddG"] = ddG

    with open(output_dir.joinpath('result.json'), "w") as fp:
        json.dump(result, fp, indent=4)
    
    return result
    