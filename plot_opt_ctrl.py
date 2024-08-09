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
simulation_horizon = 2
Ts = 0.01
plt.figure(figsize=[10, 4])
u_exact_kernel = np.loadtxt("hjb/control_exact_kernel.csv")
us = np.loadtxt("hjb/control_nystrom_kernel.csv")
u_opt_true= np.loadtxt("hjb/control_optimal.csv")

plt.plot(np.arange(0, 10 * simulation_horizon, 0.01), u_opt_true, linewidth=2, color='C2')
plt.plot(np.arange(0, 10 * simulation_horizon, 0.01), u_exact_kernel, linewidth=3, color='C1', linestyle = '-')
plt.plot(np.arange(0, 10 * simulation_horizon, 0.01), us.squeeze(), linewidth=3, linestyle='-.', color='C0')
plt.legend(["True opt.\ control", "Exact kernel", "Nys.\ Matérn-5/2"], bbox_to_anchor=(0.0, 1.02, 1.0, 0.2),
           loc='lower left',
           mode='expand',
           borderaxespad=0, ncol=3, handlelength=1.0)
plt.xlabel("t [s]")
plt.grid(visible=True, which='both')
plt.tight_layout()
plt.savefig(f"hjb/control_comparison_hjb.png", dpi=300, bbox_inches='tight',
            pad_inches=0)
plt.show()