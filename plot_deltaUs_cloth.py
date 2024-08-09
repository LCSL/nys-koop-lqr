'''
Author: Edoardo Caldarelli
Affiliation: Institut de Robòtica i Informàtica Industrial, CSIC-UPC
email: ecaldarelli@iri.upc.edu
July 2024
'''

import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})
times = np.arange(0.0, 5.05, 0.05)
for kapprox in ['nystrom', 'splines']:
    plt.figure(figsize=[10, 4])
    deltas = np.loadtxt(f"8x8_cloth_swing_xyz/sim_results/REBUTTAL_delta_us{kapprox}_seed_0.csv", delimiter=',')
    styles = ['-', '-.']
    widths = [3, 3]
    for i in range(0, 6):
        plt.plot(times, deltas[i, :].T, linewidth=widths[i % 2], linestyle=styles[i % 2])
    plt.xlabel("t [s]")
    plt.ylabel("Control components [m]")

    plt.legend(["$u_0$", "$u_1$", "$u_2$", "$u_3$", "$u_4$", "$u_5$"], bbox_to_anchor=(0.0, 1.02, 1.0, 0.2), loc='lower left',
               mode='expand',
               borderaxespad=0, ncol=6, handlelength=1.0)
    plt.grid(visible=True, which='both')
    plt.ylim(-0.15, 0.25)
    plt.tight_layout()
    plt.savefig(f"oscillating_control_{kapprox}.png", dpi=300, bbox_inches='tight',
                pad_inches=0)
    plt.show()