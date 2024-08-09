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
times = np.arange(0.0, 5.05, 0.05)
plt.figure(figsize=[8, 4])
path_to_experiment = pathlib.Path(f"./8x8_cloth_swing_xyz")
n_seeds = 50
n_states = 192
labels=['Nyström', 'Splines']
for i, kapprox in enumerate(['nys', 'splines']):
    all_rmses = []
    for seed in range(0, n_seeds):
        label = labels[i] if seed == 0 else '_'
        rmses = np.loadtxt(f"{path_to_experiment}/sim_results/REBUTTAL_rmse_{kapprox}_seed_{seed}.csv", delimiter=',') / np.sqrt(n_states)
        plt.scatter(np.argmin(rmses) * 0.05, np.amin(rmses), color=f'C{i}', alpha=0.7, label=label)
    # plt.plot(times, np.array(all_rmses).T, color=f'C{i}', linewidth=0.5, alpha=0.7)
plt.ylabel(r'Minimum $\mathrm{RMSE}$ [m]')
plt.xlabel('Time to reach min. $\mathrm{RMSE}$ [s]')
plt.legend(bbox_to_anchor=(0.0, 1.02, 1.0, 0.2), loc='lower left',
           mode='expand',
           borderaxespad=0, ncol=3, handlelength=1.0)
#plt.yscale('log')
#plt.xscale('log')

plt.grid(visible=True, which='both')
plt.tight_layout()
plt.savefig(f"rebuttal_cloth_swing_regulation_rmse_scatter.png", dpi=300, bbox_inches='tight',
            pad_inches=0)
plt.show()