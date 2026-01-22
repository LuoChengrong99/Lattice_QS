import h5py 
import numpy as np 
import matplotlib.pyplot as plt 
import os 


script_dir = os.path.dirname(os.path.abspath(__file__)) 
subdir = "11" 
file_name = "sim_rho0100_rhoA0100_Dt0.1000_alphaA00.0100_sigmaA0.1111_kappa1.0000_etaAA-0.0000_etaAB3.0000_etaBA-3.0000_etaBB-0.0000_Lx128_Ly128_N3276800_n500_t19999_20260122_144744.h5" 
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
    TA = file['results_TA'][:] 
    R1B = file['results_R1B'][:] 
    R2B = file['results_R2B'][:] 
    R3B = file['results_R3B'][:] 
    R4B = file['results_R4B'][:] 
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
rhoB0 = par_f['rhoB0'] 
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
rhoA = R1A + R2A + R3A + R4A + TA 
rhoB = R1B + R2B + R3B + R4B + TB 


plt.ion() 
x = np.arange(Lx) 
# fig, axes = plt.subplots(2, 6, figsize = (14, 4)) 
fig, axes = plt.subplots(2, 6, sharex = True, sharey = True, figsize = (14, 4)) 
ax1 = axes[0, 0]
ax2 = axes[0, 1]
ax3 = axes[0, 2]
ax4 = axes[0, 3]
ax5 = axes[0, 4]
ax6 = axes[0, 5]
ax7 = axes[1, 0]
ax8 = axes[1, 1]
ax9 = axes[1, 2]
ax10 = axes[1, 3]
ax11 = axes[1, 4]
ax12 = axes[1, 5]
# plt.subplots_adjust(hspace = 0.01, wspace = 0.5) 
plt.subplots_adjust(wspace = 0.3) 

# f_rhoA = ax1.imshow(rhoA[0, :, :], cmap = 'viridis', origin = 'lower', aspect = 'equal', interpolation = 'nearest') 
f_rhoA = ax1.imshow(rhoA[0, :, :], cmap = 'viridis', origin = 'lower', aspect = 'auto', interpolation = 'nearest') 
fig.colorbar(f_rhoA, ax = ax1, fraction = 0.046, pad = 0.04) 
f_A1 = ax2.imshow(R1A[0, :, :], cmap = 'viridis', origin = 'lower', aspect = 'auto', interpolation = 'nearest') 
fig.colorbar(f_A1, ax = ax2, fraction = 0.046, pad = 0.04) 
f_A2 = ax3.imshow(R2A[0, :, :], cmap = 'viridis', origin = 'lower', aspect = 'auto', interpolation = 'nearest') 
fig.colorbar(f_A2, ax = ax3, fraction = 0.046, pad = 0.04) 
f_A3 = ax4.imshow(R3A[0, :, :], cmap = 'viridis', origin = 'lower', aspect = 'auto', interpolation = 'nearest') 
fig.colorbar(f_A3, ax = ax4, fraction = 0.046, pad = 0.04) 
f_A4 = ax5.imshow(R4A[0, :, :], cmap = 'viridis', origin = 'lower', aspect = 'auto', interpolation = 'nearest') 
fig.colorbar(f_A4, ax = ax5, fraction = 0.046, pad = 0.04) 
f_AT = ax6.imshow(TA[0, :, :], cmap = 'viridis', origin = 'lower', aspect = 'auto', interpolation = 'nearest', vmin = 0) 
fig.colorbar(f_AT, ax = ax6, fraction = 0.046, pad = 0.04) 

f_rhoB = ax7.imshow(rhoB[0, :, :], cmap = 'viridis', origin = 'lower', aspect = 'auto', interpolation = 'nearest') 
fig.colorbar(f_rhoB, ax = ax7, fraction = 0.046, pad = 0.04) 
# ax7.set_xlabel("x") 
f_B1 = ax8.imshow(R1B[0, :, :], cmap = 'viridis', origin = 'lower', aspect = 'auto', interpolation = 'nearest') 
fig.colorbar(f_B1, ax = ax8, fraction = 0.046, pad = 0.04) 
f_B2 = ax9.imshow(R2B[0, :, :], cmap = 'viridis', origin = 'lower', aspect = 'auto', interpolation = 'nearest') 
fig.colorbar(f_B2, ax = ax9, fraction = 0.046, pad = 0.04) 
f_B3 = ax10.imshow(R3B[0, :, :], cmap = 'viridis', origin = 'lower', aspect = 'auto', interpolation = 'nearest') 
fig.colorbar(f_B3, ax = ax10, fraction = 0.046, pad = 0.04) 
f_B4 = ax11.imshow(R4B[0, :, :], cmap = 'viridis', origin = 'lower', aspect = 'auto', interpolation = 'nearest') 
fig.colorbar(f_B4, ax = ax11, fraction = 0.046, pad = 0.04) 
f_BT = ax12.imshow(TB[0, :, :], cmap = 'viridis', origin = 'lower', aspect = 'auto', interpolation = 'nearest', vmin = 0) 
fig.colorbar(f_BT, ax = ax12, fraction = 0.046, pad = 0.04) 




# fig.suptitle(fr'$\kappa = {kappa:.3f}, \rho_0 = {rho0}, \rho_{{A0}} = {rhoA0}, D_t = {Dt:.3f}, \alpha_A^0 = {alphaA0:.4f}, \eta_{{AA}} = {etaAA:.3f}, \eta_{{AB}} = {etaAB:.3f}, \eta_{{BA}} = {etaBA:.3f}, \eta_{{BB}} = {etaBB:.3f}$') 
for j in range(nn): 
    print(f"current time: t = {t_vec[j]}") 
    fig.suptitle(fr'$\kappa = {kappa:.3f}, \rho_0 = {rho0}, \rho_{{A0}} = {rhoA0}, D_t = {Dt:.3f}, \alpha_A^0 = {alphaA0:.4f}, \eta_{{AA}} = {etaAA:.3f}, \eta_{{AB}} = {etaAB:.3f}, \eta_{{BA}} = {etaBA:.3f}, \eta_{{BB}} = {etaBB:.3f}. t = {t_vec[j]}$') 

    f_rhoA.set_data(rhoA[j, :, :]) 
    f_A1.set_data(R1A[j, :, :]) 
    f_A2.set_data(R2A[j, :, :]) 
    f_A3.set_data(R3A[j, :, :]) 
    f_A4.set_data(R4A[j, :, :]) 
    f_AT.set_data(TA[j, :, :]) 

    f_rhoB.set_data(rhoB[j, :, :]) 
    f_B1.set_data(R1B[j, :, :]) 
    f_B2.set_data(R2B[j, :, :]) 
    f_B3.set_data(R3B[j, :, :]) 
    f_B4.set_data(R4B[j, :, :]) 
    f_BT.set_data(TB[j, :, :]) 

    # ax1.set_title(f"t = {t_vec[j]}") 
    # ax1.set_title(f"t = {1*n + t_vec[j]}") 
    fig.canvas.draw() 
    fig.canvas.flush_events() 
    print(np.sum(rhoA[j, :, :]) + np.sum(rhoB[j, :, :])) 
    print(np.sum(rhoA[j, :, :])) 
    print(np.sum(rhoB[j, :, :])) 
    input("Press Enter to continue...") 

