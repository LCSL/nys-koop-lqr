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
times = np.arange(0.0, 5.05, 0.05)
plt.figure(figsize=[10, 4])
path_to_experiment = pathlib.Path(f"./8x8_cloth_swing_xyz")
n_seeds = 50
n_states = 192
for i, kapprox in enumerate(['nys', 'splines']):
    all_rmses = []
    for seed in range(0, n_seeds):
        rmses = np.loadtxt(f"{path_to_experiment}/sim_results/rmse_{kapprox}_seed_{seed}.csv", delimiter=',') / np.sqrt(n_states)
        all_rmses.append(rmses)
    plt.plot(times[:50], np.median(all_rmses, axis=0)[:50], color=f'C{i}', linewidth=2)
    plt.fill_between(times[:50], np.percentile(all_rmses, axis=0, q=15)[:50],
                     np.percentile(all_rmses, axis=0, q=85)[:50], alpha=0.3, color=f'C{i}', label='_nolegend_')
plt.ylabel(r'$\mathrm{RMSE}$ [m]')
plt.xlabel('$\mathrm{t}$ [s]')
plt.legend(["Nyström RBF", "Splines"], bbox_to_anchor=(0.0, 1.02, 1.0, 0.2), loc='lower left',
           mode='expand',
           borderaxespad=0, ncol=3, handlelength=1.0)
#plt.yscale('log')
#plt.xscale('log')

plt.grid(visible=True, which='both')
plt.tight_layout()
plt.savefig(f"cloth_swing_regulation_rmse.png", dpi=300, bbox_inches='tight',
            pad_inches=0)
plt.show()