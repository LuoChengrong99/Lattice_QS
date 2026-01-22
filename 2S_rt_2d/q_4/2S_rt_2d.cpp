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
#include "H5Cpp.h" 
//
#include "rand.h" 
#include "particle_counting.h" 
#include "dynamics.h" 
#include "parameters_vec.h" 
#include "io.h" 


// 
int main(int argc, char* argv[]) //std::stoi(), int; std::stod(), double; 
{
    auto start_time = std::chrono::high_resolution_clock::now(); 
    Ran myran(2025678); //std::stoi(argv[1])
    
    const double rho0 = 50*2; 
    const double rhoA0 = 50*2; 
    const double rhoB0 = rhoA0; 
    const unsigned int Lx = 64; 
    const unsigned int Ly = 32; 
    const unsigned int N = static_cast<unsigned int>((rhoA0 + rhoB0)*Lx*Ly); 
    const unsigned int NA = static_cast<unsigned int>(rhoA0*Lx*Ly); 
    const unsigned int NB = static_cast<unsigned int>(rhoB0*Lx*Ly); 
    // 
    const double omegaAA = 9.0/10; 
    const double omegaAB = -3.0; 
    const double omegaBA = -omegaAB; 
    const double omegaBB = omegaAA; 
    // 
    const double kappa = 1.0; 
    const double v = 1.0; 
    // const double v_dia = v / std::sqrt(2); 
    const double Dt = 0.1; 
    // const double Dt_dia = Dt / std::sqrt(2); 
    const double sigma_A = 1.0/9; 
    const double alphaA0 = 0.01; 
    const double betaA0 = alphaA0 / sigma_A; 
    const double xiAA = 0.0; 
    const double xiAB = 0.0; 
    const double sigma_B = sigma_A; 
    const double alphaB0 = alphaA0; 
    const double betaB0 = alphaB0 / sigma_B; 
    const double xiBA = 0.0; 
    const double xiBB = 0.0; 
    // 
    const double etaAA = 1/(1 + sigma_A) - omegaAA; 
    const double etaAB = -omegaAB; 
    const double etaBA = -omegaBA; 
    const double etaBB = 1/(1 + sigma_B) - omegaBB; 
    const double zetaAA = (etaAA * (alphaA0 + betaA0) * (alphaA0 + betaA0) + xiAA * alphaA0) / betaA0; 
    const double zetaAB = (etaAB * (alphaA0 + betaA0) * (alphaA0 + betaA0) + xiAB * alphaA0) / betaA0; 
    const double zetaBA = (etaBA * (alphaB0 + betaB0) * (alphaB0 + betaB0) + xiBA * alphaB0) / betaB0; 
    const double zetaBB = (etaBB * (alphaB0 + betaB0) * (alphaB0 + betaB0) + xiBB * alphaB0) / betaB0; 
    // 
    const unsigned int n = 20000; 
    const unsigned int save_step = 100; 
    std::vector<unsigned int> t_vec; 
    for (std::size_t i = 0; i < n; i = i + save_step) { t_vec.push_back(i); } 
    // for (std::size_t i = 0; i < (n - 200); i = i + save_step) { t_vec.push_back(i); } // 
    // for (std::size_t i = (n - 200); i < n; i = i + 2) { t_vec.push_back(i); } // 
    t_vec.push_back(n-1); 
    const unsigned int nn = t_vec.size(); 
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

    //
    std::vector<ParticleStates> sigma(N); 
    std::vector<int> x(N); 
    std::vector<int> y(N); 
    std::vector<unsigned int> R1A(Lx*Ly, 0); 
    std::vector<unsigned int> R2A(Lx*Ly, 0); 
    std::vector<unsigned int> R3A(Lx*Ly, 0); 
    std::vector<unsigned int> R4A(Lx*Ly, 0);  
    std::vector<unsigned int> TA(Lx*Ly, 0); 
    std::vector<unsigned int> R1B(Lx*Ly, 0); 
    std::vector<unsigned int> R2B(Lx*Ly, 0); 
    std::vector<unsigned int> R3B(Lx*Ly, 0); 
    std::vector<unsigned int> R4B(Lx*Ly, 0); 
    std::vector<unsigned int> TB(Lx*Ly, 0); 
    // 
    // initialize 
    initial_setting(
        Lx, Ly, N, NA, sigma, x, y, 
        R1A, R2A, R3A, R4A, TA, 
        R1B, R2B, R3B, R4B, TB, 
        myran
    ); 
    // 
    // File.h5 
    auto now = std::chrono::system_clock::now(); 
    std::time_t t_now = std::chrono::system_clock::to_time_t(now); 
    std::tm tm = *std::localtime(&t_now); 
    std::ostringstream timestamp; 
    timestamp << std::put_time(&tm, "%Y%m%d_%H%M%S"); 
    std::ostringstream file_name; 
    file_name << "sim_rho0" << rho0 << "_rhoA0" << rhoA0 << "_Dt" << std::fixed << std::setprecision(4) << Dt << "_alphaA0" << alphaA0 << "_sigmaA" << sigma_A << "_kappa" << kappa << "_etaAA" << etaAA << "_etaAB" << etaAB << "_etaBA" << etaBA << "_etaBB" << etaBB << std::defaultfloat << "_Lx" << Lx << "_Ly" << Ly << "_N" << N << "_n" << n << "_" << timestamp.str() << ".h5"; 
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