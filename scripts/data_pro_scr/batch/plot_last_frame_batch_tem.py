import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data_dir = r"D:\AAAdata\VS code\aaaAnaconda\202511\202511RTtwoSpecies\20251210\2S_rt_1d"
figure_save_dir = r"D:\AAAdata\VS code\aaaAnaconda\202511\202511RTtwoSpecies\20251210\2S_rt_1d\figure"
os.makedirs(figure_save_dir, exist_ok = True)

h5_files = glob.glob(os.path.join(data_dir, "*.h5"))
# print("Found files:")
# for f in h5_files:
#     print(f)

for file_path in h5_files:
    print(f"\n Processing: {file_path}")

    df1 = pd.read_hdf(file_path, key = "parameters")
    df2 = pd.read_hdf(file_path, key = "results")
    rho0 = df1["rho0"].item()
    Lx = df1["Lx"].item()
    Dt = df1["Dt"].item()
    # alphaA0 = df1["alpphaA0"].item()
    # alphaB0 = df1["alpphaB0"].item()
    if 'alphaA0' in df1.columns:
        alphaA0 = df1['alphaA0'].item()
    elif 'alpphaA0' in df1.columns:
        alphaA0 = df1['alpphaA0'].item()
    else:
        raise KeyError("Neither 'alphaA0' nor 'alpphaA0' found.")
    if 'alphaB0' in df1.columns:
        alphaB0 = df1['alphaB0'].item()
    elif 'alpphaB0' in df1.columns:
        alphaB0 = df1['alpphaB0'].item()
    else:
        raise KeyError("Neither 'alphaB0' nor 'alpphaB0' found.")
    etaAA = df1["etaAA"].item()
    etaAB = df1["etaAB"].item()
    etaBA = df1["etaBA"].item()
    etaBB = df1["etaBB"].item()
    record_times = df1.loc[0, "record_times"]
    t_vec = np.array(record_times.split(","), dtype = np.int32)
    nn = len(t_vec)

    R1A = df2.pivot(index = "time", columns = "x", values = "R1A").values
    R2A = df2.pivot(index = "time", columns = "x", values = "R2A").values
    TA = df2.pivot(index = "time", columns = "x", values = "TA").values
    R1B = df2.pivot(index = "time", columns = "x", values = "R1B").values
    R2B = df2.pivot(index = "time", columns = "x", values = "R2B").values
    TB = df2.pivot(index = "time", columns = "x", values = "TB").values

    rhoA = R1A + R2A + TA
    mA = (R1A - R2A) / rho0
    rhoB = R1B + R2B + TB
    mB = (R1B - R2B) / rho0
    rhoAB = rhoA + rhoB
    max_value_rho = np.maximum(np.max(rhoA), np.max(rhoB))
    min_value_rho = np.minimum(np.min(rhoA), np.min(rhoB))
    max_value_m = np.maximum(np.max(mA), np.max(mB))
    min_value_m = np.minimum(np.min(mA), np.min(mB))
    max_value_rhoAB = np.max(rhoAB)
    min_value_rhoAB = np.min(rhoAB)

    j = nn - 1
    x = np.arange(Lx)
    fig, axes = plt.subplots(5, 1, sharex = True, figsize = (6, 6))
    ax1, ax2, ax3, ax4, ax5 = axes
    plt.subplots_adjust(hspace = 0.05)

    ax1.plot(x, rhoA[j, :], color = "red", label = r"$\rho_A$")
    ax1.grid(True)
    ax1.set_title(f"t = {t_vec[j]}")
    ax1.set_ylabel(r"$\rho_A$")
    ax1.set_ylim(min_value_rho - 5, max_value_rho + 5)

    ax2.plot(x, rhoB[j, :], color = "blue", label = r"$\rho_B$")
    ax2.grid(True)
    ax2.set_ylabel(r"$\rho_B$")
    ax2.set_ylim(min_value_rho - 5, max_value_rho + 5)

    ax3.plot(x, rhoAB[j, :], color = 'black', label = r'$\rho_A + \rho_B$') 
    ax3.grid(True) 
    ax3.set_ylabel(r'$\rho_A + \rho_B$')
    ax3.set_ylim(min_value_rhoAB - 5, max_value_rhoAB + 5) 

    ax4.plot(x, rhoA[j, :], color = 'red', label = r'$\rho_A$') 
    ax4.plot(x, rhoB[j, :], color = 'blue', label = r'$\rho_B$') 
    ax4.grid(True) 
    ax4.set_ylim(min_value_rho - 5, max_value_rho + 5) 
    ax4.legend() 

    ax5.plot(x, mA[j, :], color = 'red', label = r'$m_A$') 
    ax5.plot(x, mB[j, :], color = 'blue', label = r'$m_B$') 
    ax5.set_xlabel("x") 
    ax5.grid(True) 
    ax5.set_ylim(min_value_m*1.1, max_value_m*1.1) 
    ax5.legend()

    fig.suptitle(fr'$\rho_0 = {rho0}, D_t = {Dt:.3f}, \alpha_A^0 = {alphaA0:.4f}, \eta_{{AA}} = {etaAA:.3f}, \eta_{{AB}} = {etaAB:.3f}, \eta_{{BA}} = {etaBA:.3f}, \eta_{{BB}} = {etaBB:.3f}$') 

    figure_name = os.path.splitext(os.path.basename(file_path))[0]
    figure_save_path = os.path.join(figure_save_dir, f"{figure_name}.png")
    fig.savefig(figure_save_path, dpi = 300, bbox_inches = "tight")
    plt.close(fig)
    print(f"Saved: {figure_save_path}")

print("\n All files processed.")


