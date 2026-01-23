import h5py 
import numpy as np 
import matplotlib.pyplot as plt 
import os 


script_dir = os.path.dirname(os.path.abspath(__file__)) 
subdir = "00" 
file_name = "sim_rho0100_rhoA0100_Dt0.1000_alphaA00.0100_sigmaA0.1111_kappa1.0000_etaAA-0.0000_etaAB3.0000_etaBA-3.0000_etaBB-0.0000_Lx128_Ly128_N3276800_n20000_20260120_154823.h5" 
file_path = os.path.join(script_dir, subdir, file_name) 

with h5py.File(file_path, mode = 'r') as file: 
    print(file.keys()) 
    parameters_int = file['parameters_int'] 
    parameters_float = file['parameters_float'] 
    t_vec = file['time_vector'][:] 
    R1A = file['results_R1A'][:] 
    R2A = file['results_R2A'][:] 
    R3A = file['results_R3A'][:] 
    R4A = file['results_R4A'][:] 
    R5A = file['results_R5A'][:] 
    R6A = file['results_R6A'][:] 
    R7A = file['results_R7A'][:] 
    R8A = file['results_R8A'][:] 
    TA = file['results_TA'][:] 
    R1B = file['results_R1B'][:] 
    R2B = file['results_R2B'][:] 
    R3B = file['results_R3B'][:] 
    R4B = file['results_R4B'][:] 
    R5B = file['results_R5B'][:] 
    R6B = file['results_R6B'][:] 
    R7B = file['results_R7B'][:] 
    R8B = file['results_R8B'][:] 
    TB = file['results_TB'][:] 
    par_i = dict(zip(
        [x.decode('utf-8') for x in parameters_int['name'][:]], 
        parameters_int['value'][:] 
    )) 
    par_f = dict(zip(
        [x.decode('utf-8') for x in parameters_float['name'][:]], 
        parameters_float['value'][:] 
    )) 

print("Completed steps: ", t_vec[R1A.shape[0] - 1]) 
# print(par_f) 
rho0 = par_f['rho0'] 
rhoA0 = par_f['rhoA0'] 
Lx = par_i['Lx'] 
Ly = par_i['Ly'] 
n = par_i['n'] 
kappa = par_f['kappa'] 
Dt = par_f['Dt'] 
alphaA0 = par_f['alphaA0'] 
alphaB0 = par_f['alphaB0'] 
etaAA = par_f['etaAA'] 
etaAB = par_f['etaAB'] 
etaBA = par_f['etaBA'] 
etaBB = par_f['etaBB'] 
### 
nn = R1A.shape[0]#len(t_vec) 
rhoA = R1A + R2A + R3A + R4A + R5A + R6A + R7A + R8A + TA 
# mean_rhoA = np.mean(rhoA, axis = 2) 
# mA = R1A - R2A 
# mA = (R1A - R2A) / rho0 #/ (R1A + R2A) ### #!!!!! 
# mA = (R1A.astype(np.int32) - R2A.astype(np.int32)) / rho0 #/ (R1A + R2A) ### 
# mean_mA = np.mean(mA, axis = 2) 
rhoB = R1B + R2B + R3B + R4B + R5B + R6B + R7B + R8B + TB 
# mean_rhoB = np.mean(rhoB, axis = 2) 
# mB = (R1B - R2B) / rho0 #/ (R1B + R2B) ### #!!!!! 
# mB = (R1B.astype(np.int32) - R2B.astype(np.int32)) / rho0 #/ (R1B + R2B) ### 
# mean_mB = np.mean(mB, axis = 2) 
# rhoAB = rhoA + rhoB 
# mean_rhoAB = np.mean(rhoAB, axis = 2) 
# rhoA__B = (rhoA - rhoB) / (rhoA + rhoB) #!!!!! 
# rhoA__B = (rhoA.astype(np.int32) - rhoB.astype(np.int32)) / (rhoA + rhoB) 


# max_value_rho = np.maximum(np.max(mean_rhoA), np.max(mean_rhoB)) 
# min_value_rho = np.minimum(np.min(mean_rhoA), np.min(mean_rhoB)) 
# max_value_m = np.maximum(np.max(mean_mA), np.max(mean_mB)) 
# min_value_m = np.minimum(np.min(mean_mA), np.min(mean_mB)) 
# max_value_rhoAB = np.max(mean_rhoAB) 
# min_value_rhoAB = np.min(mean_rhoAB) 
# print(np.mean(mean_mA[0])) 
# print(np.mean(mean_mB[0])) 

plt.ion() 
x = np.arange(Lx) 
fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, figsize = (5, 12)) 
plt.subplots_adjust(hspace = 0.1) 

# f_rhoA = ax1.imshow(rhoA[0, :, :].T, cmap = 'viridis', origin = 'lower', extent = [0, Lx-1, 0, Ly-1], aspect = 'equal', interpolation = 'nearest') 
f_rhoA = ax1.imshow(rhoA[0, :, :], cmap = 'viridis', origin = 'lower', aspect = 'auto', interpolation = 'nearest') 
fig.colorbar(f_rhoA, ax = ax1, fraction = 0.046, pad = 0.04) 
# line_rhoA, = ax2.plot(x, mean_rhoA[0, :], color = 'black', label = r'$\langle \rho_A \rangle _y$') 
# line_rhoA1, = ax2.plot(x, np.mean(R1B[0], axis = 1), color = 'red', label = r'$B ^+$') 
# line_rhoA2, = ax2.plot(x, np.mean(R2B[0], axis = 1), color = 'blue', label = r'$B ^-$') 
# line_rhoAT, = ax2.plot(x, np.mean(TB[0], axis = 1), color = 'green')#, label = r'$B ^T$') 
# ax2.grid(True) 
# # ax2.set_ylim(min_value_rho - 5, max_value_rho + 5) 
# ax2.set_ylim(0, max_value_rho + 5) 
# ax2.legend(loc = 'upper right', frameon = False) 
# f_rhoB = ax2.imshow(rhoB[0, :, :].T, cmap = 'viridis', origin = 'lower', extent = [0, Lx-1, 0, Ly-1], aspect = 'equal', interpolation = 'nearest') 
f_rhoB = ax2.imshow(rhoB[0, :, :], cmap = 'viridis', origin = 'lower', aspect = 'auto', interpolation = 'nearest') 
fig.colorbar(f_rhoB, ax = ax2, fraction = 0.046, pad = 0.04) 
ax2.set_xlabel("x") 
# line_rhoB, = ax4.plot(x, mean_rhoB[0, :], color = 'black', label = r'$\langle \rho_B \rangle _y$') 
# line_rhoB1, = ax4.plot(x, np.mean(R1A[0], axis = 1), color = 'red', label = r'$A ^+$') 
# line_rhoB2, = ax4.plot(x, np.mean(R2A[0], axis = 1), color = 'blue', label = r'$A ^-$') 
# line_rhoBT, = ax4.plot(x, np.mean(TA[0], axis = 1), color = 'green')#, label = r'$A ^T$') 
# ax4.grid(True) 
# # ax4.set_ylim(min_value_rho - 5, max_value_rho + 5) 
# ax4.set_ylim(0, max_value_rho + 5) 
# ax4.legend(loc = 'upper right', frameon = False) 
# f_rhoA__B = ax5.imshow(rhoA__B[0, :, :].T, cmap = 'bwr', origin = 'lower', vmin = -1, vmax = 1, extent = [0, Lx-1, 0, Ly-1], aspect = 'auto', interpolation = 'nearest') 
# cbar = fig.colorbar(f_rhoA__B, ax = ax5, fraction = 0.046, pad = 0.04) 
# # cbar.set_label("rhoA → Red   |   rhoB → Blue", fontsize=12) 
# line_rhoAB, = ax6.plot(x, mean_rhoAB[0, :], color = 'black', label = r'$\langle \rho_A \rangle _y + \langle \rho_B \rangle _y$')  
# ax6.grid(True) 
# ax6.set_ylim(min_value_rhoAB - 5, max_value_rhoAB + 5) 
# ax6.legend(loc = 'upper right', frameon = False) 
# line_mA, = ax7.plot(x, mean_mA[0, :], color = 'red', label = r'$\langle m_A \rangle _y$') 
# line_mB, = ax7.plot(x, mean_mB[0, :], color = 'blue', label = r'$\langle m_B \rangle _y$') 
# ax7.set_xlabel("x") 
# ax7.grid(True) 
# ax7.set_ylim(min_value_m*1.1, max_value_m*1.1) 
# ax7.legend(loc = 'upper right', frameon = False) 


fig.suptitle(fr'$\kappa = {kappa:.3f}, \rho_0 = {rho0}, \rho_{{A0}} = {rhoA0}, D_t = {Dt:.3f}, \alpha_A^0 = {alphaA0:.4f}, \eta_{{AA}} = {etaAA:.3f}, \eta_{{AB}} = {etaAB:.3f}, \eta_{{BA}} = {etaBA:.3f}, \eta_{{BB}} = {etaBB:.3f}$') 
for j in range(nn): 
    print(f"current time: t = {t_vec[j]}") 

    f_rhoA.set_data(rhoA[j, :, :]) 
    # line_rhoA.set_ydata(mean_rhoA[j, :]) 
    # line_rhoA1.set_ydata(np.mean(R1B[j], axis = 1)) 
    # line_rhoA2.set_ydata(np.mean(R2B[j], axis = 1)) 
    # line_rhoAT.set_ydata(np.mean(TB[j], axis = 1)) 
    f_rhoB.set_data(rhoB[j, :, :]) 
    # line_rhoB.set_ydata(mean_rhoB[j, :]) 
    # line_rhoB1.set_ydata(np.mean(R1A[j], axis = 1)) 
    # line_rhoB2.set_ydata(np.mean(R2A[j], axis = 1)) 
    # line_rhoBT.set_ydata(np.mean(TA[j], axis = 1)) 
    # f_rhoA__B.set_data(rhoA__B[j, :, :].T) 
    # line_rhoAB.set_ydata(mean_rhoAB[j, :]) 
    # line_mA.set_ydata(mean_mA[j, :]) 
    # line_mB.set_ydata(mean_mB[j, :]) 
    ax1.set_title(f"t = {t_vec[j]}") 
    # ax1.set_title(f"t = {1*n + t_vec[j]}") 
    fig.canvas.draw() 
    fig.canvas.flush_events() 
    # print(np.mean(mean_mA[j])) 
    # print(np.mean(mean_mB[j])) 
    print(np.sum(rhoA[j, :, :]) + np.sum(rhoB[j, :, :])) 
    print(np.sum(rhoA[j, :, :])) 
    print(np.sum(rhoB[j, :, :])) 
    input("Press Enter to continue...") 


