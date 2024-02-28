'''
Author: Edoardo Caldarelli
Affiliation: Institut de Robòtica i Informàtica Industrial, CSIC-UPC
email: ecaldarelli@iri.upc.edu
January 2024
'''


import numpy as np
import matplotlib.pyplot as plt
import pathlib

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})

# ms = np.logspace(1, 2.6, num=20, dtype=int)  # np.arange(10, 200, 15)
ms = np.around(np.logspace(1, 2.3, num=20)).astype(int)  # np.arange(10, 200, 15)

system = 'duffing'
path_to_data = pathlib.Path(f"./{system}")
path_to_data.mkdir(exist_ok=True)
plt.figure(figsize=[8, 4])
orders = [30, 15, 0]
styles = ['-', '-.', '--']
for i, kapprox in enumerate(['nystrom', 'splines', 'eigfuns']):
    if kapprox == 'eigfuns':
        rmses = np.loadtxt(f"{path_to_data}/all_rmses_{kapprox}.csv", delimiter=',')
    else:
        rmses = np.loadtxt(f"{path_to_data}/all_rmses_{kapprox}_double_dataset.csv")

    plt.plot(ms, np.median(rmses, axis=0), color=f'C{i}', linewidth=3, zorder=orders[i], linestyle=styles[i])
    plt.fill_between(ms, np.percentile(rmses, axis=0, q=15),
                     np.percentile(rmses, axis=0, q=85), alpha=0.3, color=f'C{i}', label='_nolegend_', zorder=orders[i])
plt.ylabel(r'$\mathrm{RMSE}_{\%}$')
plt.xlabel('$m$')
plt.legend(["Nys.\ Matérn-5/2", "Splines", "EigFuns"], bbox_to_anchor=(0.0, 1.02, 1.0, 0.2), loc='lower left',
           mode='expand',
           borderaxespad=0, ncol=3, handlelength=1.0)
plt.xscale('log')
plt.yscale('log')
#plt.ylim(0.0, 60.0)
plt.grid(visible=True, which='both')
plt.tight_layout()
plt.savefig(f"{path_to_data}/rmse_open_loop_three_algos.png", dpi=300, bbox_inches='tight',
            pad_inches=0)
plt.show()