#include <iostream> 
#include <iomanip> 
#include <sstream> 
#include <filesystem> 
#include <string> 
#include <cstring> 
#include <chrono> 
#include <ctime> 
#include <vector> 
#include <cmath> 
#include <unordered_map> 
#include "H5Cpp.h" 
//
#include "rand.h" 
#include "particle_counting.h" 
#include "dynamics.h" 
#include "parameters_vec.h" 
#include "io.h" 


// 
int main()
{
    auto start_time = std::chrono::high_resolution_clock::now(); 
    Ran myran(2026); 

    H5::StrType str_type(H5::PredType::C_S1, 64); 
    str_type.setCset(H5T_CSET_UTF8); 
    str_type.setStrpad(H5T_STR_NULLTERM); 
    H5::CompType comp_type_i(sizeof(Par_i)); 
    H5::CompType comp_type_f(sizeof(Par_f)); 
    comp_type_i.insertMember("name", HOFFSET(Par_i, name), str_type); 
    comp_type_i.insertMember("value", HOFFSET(Par_i, value), H5::PredType::NATIVE_UINT); 
    comp_type_f.insertMember("name", HOFFSET(Par_f, name), str_type); 
    comp_type_f.insertMember("value", HOFFSET(Par_f, value), H5::PredType::NATIVE_DOUBLE); 
    // read the last 
    auto dir_last = std::filesystem::current_path(); 
    std::filesystem::path subdir_last("11"); ///// 
    ////// 
    const H5std_string file_name_last("sim_rho0100_rhoA0100_Dt0.1000_alphaA00.0100_sigmaA0.1111_kappa1.0000_etaAA-0.0000_etaAB3.0000_etaBA-3.0000_etaBB-0.0000_Lx128_Ly128_N3276800_n20000_20260122_111047.h5"); ///// 
    H5::H5File file_last((dir_last / subdir_last / file_name_last).string(), H5F_ACC_RDONLY); 
    H5::DataSet parameters_int = file_last.openDataSet("parameters_int"); 
    H5::DataSet parameters_float = file_last.openDataSet("parameters_float"); 
    H5::DataSet time_vector_last = file_last.openDataSet("time_vector"); 
    H5::DataSet res_R1A_last = file_last.openDataSet("results_R1A"); 
    H5::DataSet res_R2A_last = file_last.openDataSet("results_R2A"); 
    H5::DataSet res_R3A_last = file_last.openDataSet("results_R3A"); 
    H5::DataSet res_R4A_last = file_last.openDataSet("results_R4A");  
    H5::DataSet res_TA_last = file_last.openDataSet("results_TA"); 
    H5::DataSet res_R1B_last = file_last.openDataSet("results_R1B"); 
    H5::DataSet res_R2B_last = file_last.openDataSet("results_R2B"); 
    H5::DataSet res_R3B_last = file_last.openDataSet("results_R3B"); 
    H5::DataSet res_R4B_last = file_last.openDataSet("results_R4B"); 
    H5::DataSet res_TB_last = file_last.openDataSet("results_TB"); 
    H5::DataSpace dataspace_i_last = parameters_int.getSpace(); 
    hsize_t dims_i_last[1]; 
    dataspace_i_last.getSimpleExtentDims(dims_i_last, nullptr); 
    std::vector<Par_i> par_vec_i_last(dims_i_last[0]); 
    parameters_int.read(par_vec_i_last.data(), comp_type_i); 
    H5::DataSpace dataspace_f_last = parameters_float.getSpace(); 
    hsize_t dims_f_last[1]; 
    dataspace_f_last.getSimpleExtentDims(dims_f_last, nullptr); 
    std::vector<Par_f> par_vec_f_last(dims_f_last[0]); 
    parameters_float.read(par_vec_f_last.data(), comp_type_f); 
    H5::DataSpace dataspace_t_vec_last = time_vector_last.getSpace(); 
    hsize_t dims_time_last[1]; 
    dataspace_t_vec_last.getSimpleExtentDims(dims_time_last, nullptr); 
    hsize_t offset_t_vec_last[1] = {dims_time_last[0] - 1}; 
    hsize_t count_t_vec_last[1] = {1}; 
    dataspace_t_vec_last.selectHyperslab(H5S_SELECT_SET, count_t_vec_last, offset_t_vec_last); 
    H5::DataSpace memspace_t_vec_last(1, count_t_vec_last); 
    unsigned int t_last; 
    time_vector_last.read(&t_last, H5::PredType::NATIVE_UINT, memspace_t_vec_last, dataspace_t_vec_last); 
    H5::DataSpace dataspace_res_last = res_R1A_last.getSpace(); 
    hsize_t dims_res_last[3]; 
    dataspace_res_last.getSimpleExtentDims(dims_res_last, nullptr); 
    hsize_t offset_res_last[3] = {dims_res_last[0] - 1, 0, 0}; 
    hsize_t count_res_last[3] = {1, dims_res_last[1], dims_res_last[2]}; 
    dataspace_res_last.selectHyperslab(H5S_SELECT_SET, count_res_last, offset_res_last); 
    H5::DataSpace memspace_res_last(3, count_res_last); 
    const std::size_t Lx0Ly0 = dims_res_last[1] * dims_res_last[2]; 
    std::vector<unsigned int> R1A(Lx0Ly0); 
    std::vector<unsigned int> R2A(Lx0Ly0); 
    std::vector<unsigned int> R3A(Lx0Ly0); 
    std::vector<unsigned int> R4A(Lx0Ly0); 
    std::vector<unsigned int> TA(Lx0Ly0); 
    std::vector<unsigned int> R1B(Lx0Ly0); 
    std::vector<unsigned int> R2B(Lx0Ly0); 
    std::vector<unsigned int> R3B(Lx0Ly0); 
    std::vector<unsigned int> R4B(Lx0Ly0); 
    std::vector<unsigned int> TB(Lx0Ly0); 
    res_R1A_last.read(R1A.data(), H5::PredType::NATIVE_UINT, memspace_res_last, dataspace_res_last); 
    res_R2A_last.read(R2A.data(), H5::PredType::NATIVE_UINT, memspace_res_last, dataspace_res_last); 
    res_R3A_last.read(R3A.data(), H5::PredType::NATIVE_UINT, memspace_res_last, dataspace_res_last); 
    res_R4A_last.read(R4A.data(), H5::PredType::NATIVE_UINT, memspace_res_last, dataspace_res_last); 
    res_TA_last.read(TA.data(), H5::PredType::NATIVE_UINT, memspace_res_last, dataspace_res_last); 
    res_R1B_last.read(R1B.data(), H5::PredType::NATIVE_UINT, memspace_res_last, dataspace_res_last); 
    res_R2B_last.read(R2B.data(), H5::PredType::NATIVE_UINT, memspace_res_last, dataspace_res_last); 
    res_R3B_last.read(R3B.data(), H5::PredType::NATIVE_UINT, memspace_res_last, dataspace_res_last); 
    res_R4B_last.read(R4B.data(), H5::PredType::NATIVE_UINT, memspace_res_last, dataspace_res_last); 
    res_TB_last.read(TB.data(), H5::PredType::NATIVE_UINT, memspace_res_last, dataspace_res_last); 
    file_last.close(); 
    // 
    std::unordered_map<std::string, unsigned int> par_i; 
    for (const auto& p : par_vec_i_last) { par_i.emplace(p.name, p.value); } 
    std::unordered_map<std::string, double> par_f; 
    for (const auto& p : par_vec_f_last) { par_f.emplace(p.name, p.value); } 
    //
 
    //
    const double rho0 = par_f.at("rho0"); 
    const double rhoA0 = par_f.at("rhoA0"); 
    const double rhoB0 = par_f.at("rhoB0"); 
    const unsigned int Lx = par_i.at("Lx"); // 
    const unsigned int Ly = par_i.at("Ly"); 
    const unsigned int N = par_i.at("N"); 
    const unsigned int NA = par_i.at("NA"); 
    const unsigned int NB = par_i.at("NB"); // 
    // const unsigned int Lx = static_cast<unsigned int>(2*par_i.at("Lx")); //_con_LxLy 
    // const unsigned int Ly = static_cast<unsigned int>(2*par_i.at("Ly"));
    // const unsigned int N = static_cast<unsigned int>(4*par_i.at("N")); 
    // const unsigned int NA = static_cast<unsigned int>(4*par_i.at("NA")); 
    // const unsigned int NB = static_cast<unsigned int>(4*par_i.at("NB")); //_con_LxLy 

    const double omegaAA = par_f.at("omegaAA"); 
    const double omegaAB = par_f.at("omegaAB"); 
    const double omegaBA = par_f.at("omegaBA"); 
    const double omegaBB = par_f.at("omegaBB"); 
    // 
    const double kappa = par_f.at("kappa"); 
    const double v = par_f.at("v"); 
    // const double v_dia = v / std::sqrt(2); 
    const double Dt = par_f.at("Dt"); 
    // const double Dt_dia = Dt / std::sqrt(2); 
    const double sigma_A = par_f.at("sigmaA"); 
    const double alphaA0 = par_f.at("alphaA0"); 
    const double betaA0 = par_f.at("betaA0"); 
    const double xiAA = par_f.at("xiAA"); 
    const double xiAB = par_f.at("xiAB"); 
    const double sigma_B = par_f.at("sigmaB"); 
    const double alphaB0 = par_f.at("alphaB0"); 
    const double betaB0 = par_f.at("betaB0"); 
    const double xiBA = par_f.at("xiBA"); 
    const double xiBB = par_f.at("xiBB"); 
    // 
    const double etaAA = par_f.at("etaAA"); 
    const double etaAB = par_f.at("etaAB"); 
    const double etaBA = par_f.at("etaBA"); 
    const double etaBB = par_f.at("etaBB"); 
    const double zetaAA = par_f.at("zetaAA"); 
    const double zetaAB = par_f.at("zetaAB"); 
    const double zetaBA = par_f.at("zetaBA"); 
    const double zetaBB = par_f.at("zetaBB"); 
    // 
    // const unsigned int n = par_i.at("n"); ///// 
    // const unsigned int save_step = 10; ///// 
    const unsigned int n = 500; ///// 
    const unsigned int save_step = 2; ///// 
    std::vector<unsigned int> t_vec; 
    for (std::size_t i = 0; i < n; i = i + save_step) { t_vec.push_back(i); } 
    t_vec.push_back(n-1); 
    const unsigned int nn = t_vec.size(); 

    std::vector<ParticleStates> sigma(N); 
    std::vector<int> x(N); 
    std::vector<int> y(N); 
    // 
    // initialize 
    std::size_t ii = 0; 
    for (std::size_t jj = 0; jj < Lx; ++jj) {
        for (std::size_t kk = 0; kk < Ly; ++kk) {
            for (std::size_t i = ii; i < (ii + R1A[jj + kk * Lx]); ++i) {
                x[i] = jj; 
                y[i] = kk; 
                sigma[i] = ParticleStates::Running_1_state; 
            }
            ii += R1A[jj + kk * Lx]; 
        }
    } 
    for (std::size_t jj = 0; jj < Lx; ++jj) {
        for (std::size_t kk = 0; kk < Ly; ++kk) {
            for (std::size_t i = ii; i < (ii + R2A[jj + kk * Lx]); ++i) {
                x[i] = jj; 
                y[i] = kk; 
                sigma[i] = ParticleStates::Running_2_state; 
            } 
            ii += R2A[jj + kk * Lx]; 
        }
    } 
    for (std::size_t jj = 0; jj < Lx; ++jj) {
        for (std::size_t kk = 0; kk < Ly; ++kk) {
            for (std::size_t i = ii; i < (ii + R3A[jj + kk * Lx]); ++i) {
                x[i] = jj; 
                y[i] = kk; 
                sigma[i] = ParticleStates::Running_3_state; 
            } 
            ii += R3A[jj + kk * Lx]; 
        }
    } 
    for (std::size_t jj = 0; jj < Lx; ++jj) {
        for (std::size_t kk = 0; kk < Ly; ++kk) {
            for (std::size_t i = ii; i < (ii + R4A[jj + kk * Lx]); ++i) {
                x[i] = jj; 
                y[i] = kk; 
                sigma[i] = ParticleStates::Running_4_state; 
            } 
            ii += R4A[jj + kk * Lx]; 
        }
    } 
    for (std::size_t jj = 0; jj < Lx; ++jj) {
        for (std::size_t kk = 0; kk < Ly; ++kk) {
            for (std::size_t i = ii; i < (ii + TA[jj + kk * Lx]); ++i) {
                x[i] = jj; 
                y[i] = kk; 
                sigma[i] = ParticleStates::Tumbling_0_state; 
            }
            ii += TA[jj + kk * Lx]; 
        }
    } 
    for (std::size_t jj = 0; jj < Lx; ++jj) {
        for (std::size_t kk = 0; kk < Ly; ++kk) {
            for (std::size_t i = ii; i < (ii + R1B[jj + kk * Lx]); ++i) {
                x[i] = jj; 
                y[i] = kk; 
                sigma[i] = ParticleStates::Running_1_state; 
            }
            ii += R1B[jj + kk * Lx]; 
        }
    } 
    for (std::size_t jj = 0; jj < Lx; ++jj) {
        for (std::size_t kk = 0; kk < Ly; ++kk) {
            for (std::size_t i = ii; i < (ii + R2B[jj + kk * Lx]); ++i) {
                x[i] = jj; 
                y[i] = kk; 
                sigma[i] = ParticleStates::Running_2_state; 
            }
            ii += R2B[jj + kk * Lx]; 
        }
    } 
    for (std::size_t jj = 0; jj < Lx; ++jj) {
        for (std::size_t kk = 0; kk < Ly; ++kk) {
            for (std::size_t i = ii; i < (ii + R3B[jj + kk * Lx]); ++i) {
                x[i] = jj; 
                y[i] = kk; 
                sigma[i] = ParticleStates::Running_3_state; 
            }
            ii += R3B[jj + kk * Lx]; 
        }
    } 
    for (std::size_t jj = 0; jj < Lx; ++jj) {
        for (std::size_t kk = 0; kk < Ly; ++kk) {
            for (std::size_t i = ii; i < (ii + R4B[jj + kk * Lx]); ++i) {
                x[i] = jj; 
                y[i] = kk; 
                sigma[i] = ParticleStates::Running_4_state; 
            }
            ii += R4B[jj + kk * Lx]; 
        }
    } 
    for (std::size_t jj = 0; jj < Lx; ++jj) {
        for (std::size_t kk = 0; kk < Ly; ++kk) {
            for (std::size_t i = ii; i < (ii + TB[jj + kk * Lx]); ++i) {
                x[i] = jj; 
                y[i] = kk; 
                sigma[i] = ParticleStates::Tumbling_0_state; 
            }
            ii += TB[jj + kk * Lx]; 
        }
    } 
    std::cout << "ii: " << ii << std::endl; 

    // 
    std::vector<Par_i> par_vec_i; 
    parameters_vec_i(Lx, Ly, N, NA, NB, n, par_vec_i); 
    std::vector<Par_f> par_vec_f; 
    parameters_vec_f(
        rho0, rhoA0, rhoB0, 
        omegaAA, omegaAB, omegaBA, omegaBB, kappa, v, Dt, 
        sigma_A, alphaA0, zetaAA, zetaAB, betaA0, xiAA, xiAB, 
        sigma_B, alphaB0, zetaBA, zetaBB, betaB0, xiBA, xiBB, 
        etaAA, etaAB, etaBA, etaBB, par_vec_f
    ); 
    // File.h5 
    auto now = std::chrono::system_clock::now(); 
    std::time_t t_now = std::chrono::system_clock::to_time_t(now); 
    std::tm tm = *std::localtime(&t_now); 
    std::ostringstream timestamp; 
    timestamp << std::put_time(&tm, "%Y%m%d_%H%M%S"); 
    std::ostringstream file_name; 
    file_name << "sim_rho0" << rho0 << "_rhoA0" << rhoA0 << "_Dt" << std::fixed << std::setprecision(4) << Dt << "_alphaA0" << alphaA0 << "_sigmaA" << sigma_A << "_kappa" << kappa << "_etaAA" << etaAA << "_etaAB" << etaAB << "_etaBA" << etaBA << "_etaBB" << etaBB << std::defaultfloat << "_Lx" << Lx << "_Ly" << Ly << "_N" << N << "_n" << n << "_t" << t_last << "_" << timestamp.str() << ".h5"; 
    auto dir = std::filesystem::current_path(); 
    std::filesystem::path subdir("11"); ////// 
    initial_setting_wri(
        dir, subdir, file_name, Lx, Ly, par_vec_i, par_vec_f, t_vec, 
        &R1A[0], &R2A[0], &R3A[0], &R4A[0], &TA[0], 
        &R1B[0], &R2B[0], &R3B[0], &R4B[0], &TB[0]
    ); 
 
    // 
    std::vector<unsigned int> idx(N); 
    for (std::size_t i = 0; i < N; ++i) { idx[i] = i; } 
    unsigned int j = 1; 
    for (std::size_t t = 1; t < n; ++t) 
    {
        shuffle_idx(idx, myran); 
        for (std::size_t i = 0; i < N; ++i) 
        {
            double r = myran.doub() * ((1 + kappa)*(1 + kappa)*std::max(std::max(alphaA0, betaA0), std::max(alphaB0, betaB0)) + v + 4*Dt); 
            int xx = x[idx[i]]; 
            int yy = y[idx[i]]; 
            int jj = y[idx[i]] * Lx + x[idx[i]]; 
            ParticleStates s = sigma[idx[i]]; 
            unsigned int rhoA_xy = R1A[jj] + R2A[jj] + R3A[jj] + R4A[jj] + TA[jj]; 
            unsigned int rhoB_xy = R1B[jj] + R2B[jj] + R3B[jj] + R4B[jj] + TB[jj]; 
            if (idx[i] < NA) 
            {
                if (s != ParticleStates::Tumbling_0_state) 
                {
                    running_state(Lx, Ly, kappa, v, Dt, alphaA0, zetaAA, zetaAB, rho0, rhoA_xy, rhoB_xy, r, 
                        x[idx[i]], y[idx[i]], sigma[idx[i]]
                    ); 
                    counting_running_A(
                        Lx, Ly, 
                        jj, x[idx[i]], y[idx[i]], 
                        s, sigma[idx[i]], 
                        R1A, R2A, R3A, R4A, TA 
                    ); 
                }
                else 
                {
                    tumbling_state(Lx, Ly, kappa, Dt, betaA0, xiAA, xiAB, rho0, rhoA_xy, rhoB_xy, r, 
                        x[idx[i]], y[idx[i]], sigma[idx[i]]
                    ); 
                    counting_tumbling_A(
                        Lx, Ly, 
                        jj, x[idx[i]], y[idx[i]], 
                        s, sigma[idx[i]], 
                        R1A, R2A, R3A, R4A, TA 
                    ); 
                }
            }
            else 
            {
                if (s != ParticleStates::Tumbling_0_state) 
                {
                    running_state(Lx, Ly, kappa, v, Dt, alphaB0, zetaBA, zetaBB, rho0, rhoA_xy, rhoB_xy, r, 
                        x[idx[i]], y[idx[i]], sigma[idx[i]]
                    ); 
                    counting_running_B(
                        Lx, Ly, 
                        jj, x[idx[i]], y[idx[i]], 
                        s, sigma[idx[i]], 
                        R1B, R2B, R3B, R4B, TB 
                    ); 
                }
                else 
                {
                    tumbling_state(Lx, Ly, kappa, Dt, betaB0, xiBA, xiBB, rho0, rhoA_xy, rhoB_xy, r, 
                        x[idx[i]], y[idx[i]], sigma[idx[i]]
                    ); 
                    counting_tumbling_B(
                        Lx, Ly, 
                        jj, x[idx[i]], y[idx[i]], 
                        s, sigma[idx[i]], 
                        R1B, R2B, R3B, R4B, TB 
                    ); 
                }
            }
        }
        if (j < nn && t == t_vec[j]) 
        {
            results_wri(
                dir, subdir, file_name, j, Lx, Ly, 
                &R1A[0], &R2A[0], &R3A[0], &R4A[0], &TA[0], 
                &R1B[0], &R2B[0], &R3B[0], &R4B[0], &TB[0] 
            ); 
            std::cout << "j = " << j << "\t"; 
            auto end_time = std::chrono::high_resolution_clock::now(); 
            std::chrono::duration<double> running_time = end_time - start_time; 
            std::cout << std::fixed << std::setprecision(4) << "Running time: " << running_time.count() << "s" << std::endl; 
            j += 1; 
            // for (std::size_t i = 0; i < Lx*Ly; i++) { std::cout << R1A[i] << "\t"; } std::cout << "\n"; 
        }
    }


    //////
    std::cout << std::endl; 
    unsigned int A = 0, A11 = 0, A22 = 0, A33 = 0, A44 = 0, AT = 0; 
    unsigned int B = 0, B11 = 0, B22 = 0, B33 = 0, B44 = 0, BT = 0; 
    for (std::size_t i = 0; i < (Lx * Ly); ++i) 
    {
        A = A + R1A[i] + R2A[i] + R3A[i] + R4A[i] + TA[i]; 
        A11 = A11 + R1A[i]; 
        A22 = A22 + R2A[i]; 
        A33 = A33 + R3A[i]; 
        A44 = A44 + R4A[i]; 
        AT = AT + TA[i]; 
        B = B + R1B[i] + R2B[i] + R3B[i] + R4B[i] + TB[i]; 
        B11 = B11 + R1B[i]; 
        B22 = B22 + R2B[i]; 
        B33 = B33 + R3B[i]; 
        B44 = B44 + R4B[i]; 
        BT = BT + TB[i]; 
    }
    std::cout << "A:" << A << "\t" << "A1:" << A11 << "\t" << "A2:" << A22 << "\t" << "A3:" << A33 << "\t" << "A4:" << A44 << "\t" << "A tumble:" << AT << std::endl; 
    std::cout << "B:" << B << "\t" << "B1:" << B11 << "\t" << "B2:" << B22 << "\t" << "B3:" << B33 << "\t" << "B4:" << B44 << "\t" << "B tumble:" << BT << std::endl; 
    //////
    auto end_time = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> running_time = end_time - start_time; 
    double total_time = running_time.count(); 
    if (total_time < 60) {
        std::cout << std::fixed << std::setprecision(4) << "Total running time: " << total_time << "s" << std::endl; 
    } 
    else if (total_time < 3600) {
        total_time = total_time / 60; 
        std::cout << std::fixed << std::setprecision(4) << "Total running time: " << total_time << "m" << std::endl; 
    }
    else {
        total_time = total_time / 3600; 
        std::cout << std::fixed << std::setprecision(4) << "Total running time: " << total_time << "h" << std::endl; 
    }

    return 0; 
}