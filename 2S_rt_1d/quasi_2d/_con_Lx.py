import h5py 
import numpy as np 
import time 
import os 


script_dir = os.path.dirname(os.path.abspath(__file__)) 
subdir = "rhoA0_rhoB0" 
file_name_last_str = f"sim_t49999_rho050_rhoA040_Dt0.1000_alphaA00.0100_sigmaA0.1111_kappa1.0000_etaAA-0.0000_etaAB3.0000_etaBA3.0000_etaBB-0.0000_Lx2048_Ly8_N1310720_n50000_20260108_202713" 
file_name_last = f"{file_name_last_str}.h5" 
file_path_last = os.path.join(script_dir, subdir, file_name_last) 
timestamp = time.strftime("%Y%m%d_%H%M%S")  
file_name = f"con_Lx_{file_name_last_str}_{timestamp}.h5" 
file_path = os.path.join(script_dir, subdir, file_name) 
datasets_to_copy = [
    'parameters_int', 'parameters_float', 'time_vector' 
]
with h5py.File(file_path_last, mode = 'r') as file_last: 
    with h5py.File(file_path, mode = 'w') as file: 
        for name in datasets_to_copy: 
            file_last.copy(file_last[name], file, name = name) 
    parameters_int = file_last['parameters_int'] 
    R1AL = file_last['results_R1A'][-1] 
    R2AL = file_last['results_R2A'][-1] 
    TAL = file_last['results_TA'][-1] 
    R1BL = file_last['results_R1B'][-1] 
    R2BL = file_last['results_R2B'][-1] 
    TBL = file_last['results_TB'][-1] 
    # R1AR = file_last['results_R1A'][-1] #->-> 
    # R2AR = file_last['results_R2A'][-1] 
    # TAR = file_last['results_TA'][-1] 
    # R1BR = file_last['results_R1B'][-1] 
    # R2BR = file_last['results_R2B'][-1] 
    # TBR = file_last['results_TB'][-1] #->-> 
    R1AR = file_last['results_R2A'][-1] #-><- 
    R2AR = file_last['results_R1A'][-1] 
    TAR = file_last['results_TA'][-1] 
    R1BR = file_last['results_R2B'][-1] 
    R2BR = file_last['results_R1B'][-1] 
    TBR = file_last['results_TB'][-1] #-><- 
    par_i = dict(zip(
        [x.decode('utf-8') for x in parameters_int['name'][:]], 
        parameters_int['value'][:] 
    )) 

Lx0 = par_i['Lx'] 
Ly = par_i['Ly'] 
Lx = 2*Lx0 
R1A = np.vstack((R1AL, R1AR)) 
R2A = np.vstack((R2AL, R2AR)) 
TA = np.vstack((TAL, TAR)) 
R1B = np.vstack((R1BL, R1BR)) 
R2B = np.vstack((R2BL, R2BR)) 
TB = np.vstack((TBL, TBR)) 
file = h5py.File(file_path, 'a') 
res_R1A = file.create_dataset('results_R1A', shape = (1, Lx, Ly), maxshape = None, dtype = np.int32, chunks = (1, Lx, Ly)) 
res_R2A = file.create_dataset('results_R2A', shape = (1, Lx, Ly), maxshape = None, dtype = np.int32, chunks = (1, Lx, Ly)) 
res_TA = file.create_dataset('results_TA', shape = (1, Lx, Ly), maxshape = None, dtype = np.int32, chunks = (1, Lx, Ly)) 
res_R1B = file.create_dataset('results_R1B', shape = (1, Lx, Ly), maxshape = None, dtype = np.int32, chunks = (1, Lx, Ly)) 
res_R2B = file.create_dataset('results_R2B', shape = (1, Lx, Ly), maxshape = None, dtype = np.int32, chunks = (1, Lx, Ly)) 
res_TB = file.create_dataset('results_TB', shape = (1, Lx, Ly), maxshape = None, dtype = np.int32, chunks = (1, Lx, Ly)) 
res_R1A[0] = R1A 
res_R2A[0] = R2A 
res_TA[0] = TA 
res_R1B[0] = R1B 
res_R2B[0] = R2B 
res_TB[0] = TB 
file.close() 