'''
Author: Edoardo Caldarelli
Affiliation: Institut de Robòtica i Informàtica Industrial, CSIC-UPC
email: ecaldarelli@iri.upc.edu
July 2024
'''


import numpy as np
import matplotlib.pyplot as plt
import pathlib

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})

# ms = np.logspace(1, 2.6, num=20, dtype=int)  # np.arange(10, 200, 15)
ms = np.around(np.logspace(1, 2.3, num=20)).astype(int)  # np.arange(10, 200, 15)

system = 'hjb'
path_to_data = pathlib.Path(f"./{system}")
path_to_data.mkdir(exist_ok=True)
plt.figure(figsize=[10, 4])
for i, kapprox in enumerate(['nystrom']):

    rmses = np.loadtxt(f"{path_to_data}/control_rmse.csv").T

    plt.plot(ms[:-3], np.median(rmses, axis=0), color=f'C{i}', linewidth=3)
    plt.fill_between(ms[:-3], np.percentile(rmses, axis=0, q=15),
                     np.percentile(rmses, axis=0, q=85), alpha=0.3, color=f'C{i}', label='_nolegend_')
plt.axhline(13.79886699903511, linewidth=3, color='C1', linestyle='--')
plt.ylabel(r'$\mathrm{RMSE}_{\%}^{\mathbf u}$')
plt.xlabel('$m$')
plt.legend(["Nys.\ Matérn-5/2", "Exact kernel"], bbox_to_anchor=(0.0, 1.02, 1.0, 0.2), loc='lower left',
           mode='expand',
           borderaxespad=0, ncol=3, handlelength=1.0)
# plt.xscale('log')
#plt.yscale('log')
#plt.ylim(0.0, 60.0)
plt.grid(visible=True, which='both')
plt.tight_layout()
plt.savefig(f"{path_to_data}/control_rmse_hjb.png", dpi=300, bbox_inches='tight',
            pad_inches=0)
plt.show()