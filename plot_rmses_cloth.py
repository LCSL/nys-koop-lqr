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

ms = np.logspace(1, 2.6, num=20, dtype=int) # np.arange(10, 200, 15)
plt.figure(figsize=[10, 4])
path_to_experiment = pathlib.Path(f"./8x8_cloth_swing_xyz")

for i, kapprox in enumerate(['nystrom', 'splines']):
    pathdata = pathlib.Path(f"{path_to_experiment}/sim_results/{kapprox}/data")
    rmses = np.loadtxt(f"{pathdata}/all_rmses_{kapprox}_cloth_swing_angle.csv")
    plt.plot(ms, np.median(rmses, axis=0), color=f'C{i}', linewidth=2)
    plt.fill_between(ms, np.percentile(rmses, axis=0, q=15),
                     np.percentile(rmses, axis=0, q=85), alpha=0.3, color=f'C{i}', label='_nolegend_')
plt.ylabel(r'$\mathrm{RMSE}$ [m]')
plt.xlabel('$m$')
plt.legend(["Nyström RBF", "Splines"], bbox_to_anchor=(0.0, 1.02, 1.0, 0.2), loc='lower left',
           mode='expand',
           borderaxespad=0, ncol=3, handlelength=1.0)
plt.yscale('log')
plt.xscale('log')

# plt.ylim(0.03, 0.4)
plt.grid(visible=True, which='both')
plt.tight_layout()
plt.savefig(f"cloth_swing_rmse.png", dpi=300, bbox_inches='tight',
             pad_inches=0)
plt.show()