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
#include "dynamics.h" 
#include "parameters_vec.h" 
#include "io.h" 


// 
int main(int argc, char* argv[])
{
    auto start_time = std::chrono::high_resolution_clock::now(); 
    Ran myran(2025); 
    
    const double rho0 = 50; 
    // const double rhoA0 = 45.55; ////// 
    const double rhoA0 = atof(argv[1]); ////// 
    const double rhoB0 = rhoA0; 
    const unsigned int Lx = 1024*4*4; 
    const unsigned int Ly = 8; 
    const unsigned int N = static_cast<unsigned int>((rhoA0 + rhoB0)*Lx*Ly); 
    const unsigned int NA = static_cast<unsigned int>(rhoA0*Lx*Ly); 
    const unsigned int NB = static_cast<unsigned int>(rhoB0*Lx*Ly); 
    // 
    const double omegaAA = 9.0/10; 
    const double omegaAB = -3.0; 
    const double omegaBA = omegaAB; 
    const double omegaBB = omegaAA; 
    // 
    const double kappa = 1; 
    const double v = 1.0; 
    const double Dt = 0.1; 
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
    const unsigned int n = 50000; 
    const unsigned int save_step = 100; 
    std::vector<unsigned int> t_vec; 
    for (std::size_t i = 0; i < n; i = i + save_step) { t_vec.push_back(i); } 
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
    unsigned int R1A[Lx][Ly] = {0}; 
    unsigned int R2A[Lx][Ly] = {0}; 
    unsigned int TA[Lx][Ly] = {0}; 
    unsigned int R1B[Lx][Ly] = {0}; 
    unsigned int R2B[Lx][Ly] = {0}; 
    unsigned int TB[Lx][Ly] = {0}; 
    // 
    // initialize 
    for (std::size_t i = 0; i < N; ++i) 
    {
        double r = myran.doub(); 
        x[i] = static_cast<int> (myran.doub() * Lx); 
        y[i] = static_cast<int> (myran.doub() * Ly); 
        if (i < NA) 
        {
            if (r < 0.25) 
            {
                R1A[x[i]][y[i]] += 1; 
                sigma[i] = ParticleStates::Running_right_state; 
            }
            else if (r < 0.5) 
            {
                R2A[x[i]][y[i]] += 1; 
                sigma[i] = ParticleStates::Running_left_state; 
            }
            else 
            {
                TA[x[i]][y[i]] += 1; 
                sigma[i] = ParticleStates::Tumbling_0_state; 
            }
        }
        else 
        {
            if (r < 0.25) 
            {
                R1B[x[i]][y[i]] += 1; 
                sigma[i] = ParticleStates::Running_right_state; 
            }
            else if (r < 0.5) 
            {
                R2B[x[i]][y[i]] += 1; 
                sigma[i] = ParticleStates::Running_left_state; 
            }
            else 
            {
                TB[x[i]][y[i]] += 1; 
                sigma[i] = ParticleStates::Tumbling_0_state; 
            }
        }
    }
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
    std::filesystem::path subdir("rhoA0_rhoB0"); ////// 
    initial_setting_wri(
        dir, subdir, file_name, Lx, Ly, par_vec_i, par_vec_f, t_vec, 
        &R1A[0][0], &R2A[0][0], &TA[0][0], &R1B[0][0], &R2B[0][0], &TB[0][0]
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
            ParticleStates s = sigma[idx[i]]; 
            unsigned int rhoA_xy = R1A[xx][yy] + R2A[xx][yy] + TA[xx][yy]; 
            unsigned int rhoB_xy = R1B[xx][yy] + R2B[xx][yy] + TB[xx][yy]; 
            if (idx[i] < NA) 
            {
                if (s != ParticleStates::Tumbling_0_state) 
                {
                    running_state(Lx, Ly, kappa, v, Dt, alphaA0, zetaAA, zetaAB, rho0, rhoA_xy, rhoB_xy, r, x[idx[i]], y[idx[i]], sigma[idx[i]]); 
                    if (s == ParticleStates::Running_right_state) 
                    {
                        if (sigma[idx[i]] == ParticleStates::Tumbling_0_state) 
                        {
                            R1A[xx][yy] -= 1; 
                            TA[x[idx[i]]][y[idx[i]]] += 1; 
                        }
                        else 
                        {
                            R1A[xx][yy] -= 1; 
                            R1A[x[idx[i]]][y[idx[i]]] += 1; 
                        }
                    }
                    else 
                    {
                        if (sigma[idx[i]] == ParticleStates::Tumbling_0_state) 
                        {
                            R2A[xx][yy] -= 1; 
                            TA[x[idx[i]]][y[idx[i]]] += 1; 
                        }
                        else 
                        {
                            R2A[xx][yy] -= 1; 
                            R2A[x[idx[i]]][y[idx[i]]] += 1; 
                        }
                    }
                }
                else 
                {
                    double r0 = myran.doub(); 
                    tumbling_state(Lx, Ly, kappa, Dt, betaA0, xiAA, xiAB, rho0, rhoA_xy, rhoB_xy, r, r0, x[idx[i]], y[idx[i]], sigma[idx[i]]); 
                    if (sigma[idx[i]] == ParticleStates::Running_right_state) 
                    {
                        TA[xx][yy] -= 1; 
                        R1A[x[idx[i]]][y[idx[i]]] += 1; 
                    }
                    else if (sigma[idx[i]] == ParticleStates::Running_left_state) 
                    {
                        TA[xx][yy] -= 1; 
                        R2A[x[idx[i]]][y[idx[i]]] += 1; 
                    }
                    else 
                    {
                        TA[xx][yy] -= 1; 
                        TA[x[idx[i]]][y[idx[i]]] += 1; 
                    }
                }
            }
            else 
            {
                if (s != ParticleStates::Tumbling_0_state) 
                {
                    running_state(Lx, Ly, kappa, v, Dt, alphaB0, zetaBA, zetaBB, rho0, rhoA_xy, rhoB_xy, r, x[idx[i]], y[idx[i]], sigma[idx[i]]); 
                    if (s == ParticleStates::Running_right_state) 
                    {
                        if (sigma[idx[i]] == ParticleStates::Tumbling_0_state) 
                        {
                            R1B[xx][yy] -= 1; 
                            TB[x[idx[i]]][y[idx[i]]] += 1; 
                        }
                        else 
                        {
                            R1B[xx][yy] -= 1; 
                            R1B[x[idx[i]]][y[idx[i]]] += 1; 
                        }
                    }
                    else 
                    {
                        if (sigma[idx[i]] == ParticleStates::Tumbling_0_state) 
                        {
                            R2B[xx][yy] -= 1; 
                            TB[x[idx[i]]][y[idx[i]]] += 1; 
                        }
                        else 
                        {
                            R2B[xx][yy] -= 1; 
                            R2B[x[idx[i]]][y[idx[i]]] += 1; 
                        }
                    }
                }
                else 
                {
                    double r0 = myran.doub(); 
                    tumbling_state(Lx, Ly, kappa, Dt, betaB0, xiBA, xiBB, rho0, rhoA_xy, rhoB_xy, r, r0, x[idx[i]], y[idx[i]], sigma[idx[i]]); 
                    if (sigma[idx[i]] == ParticleStates::Running_right_state) 
                    {
                        TB[xx][yy] -= 1; 
                        R1B[x[idx[i]]][y[idx[i]]] += 1; 
                    }
                    else if (sigma[idx[i]] == ParticleStates::Running_left_state) 
                    {
                        TB[xx][yy] -= 1; 
                        R2B[x[idx[i]]][y[idx[i]]] += 1; 
                    }
                    else 
                    {
                        TB[xx][yy] -= 1; 
                        TB[x[idx[i]]][y[idx[i]]] += 1; 
                    }
                }
            }
        }
        if (j < nn && t == t_vec[j]) 
        {
            results_wri(
                dir, subdir, file_name, j, Lx, Ly, 
                &R1A[0][0], &R2A[0][0], &TA[0][0], &R1B[0][0], &R2B[0][0], &TB[0][0] 
            ); 
            std::cout << "j = " << j << "\t"; 
            auto end_time = std::chrono::high_resolution_clock::now(); 
            std::chrono::duration<double> running_time = end_time - start_time; 
            std::cout << std::fixed << std::setprecision(4) << "Running time: " << running_time.count() << "s" << std::endl; 
            j += 1; 
        }
    }


    //////
    std::cout << std::endl; 
    unsigned int A = 0, A11 = 0, A22 = 0, AT = 0, B = 0, B11 = 0, B22 = 0, BT = 0; 
    for (std::size_t i = 0; i < Lx; ++i) 
    {
        for (std::size_t j = 0; j < Ly; ++j) 
        {
            A = A + R1A[i][j] + R2A[i][j] + TA[i][j]; 
            A11 = A11 + R1A[i][j]; 
            A22 = A22 + R2A[i][j]; 
            AT = AT + TA[i][j]; 
            B = B + R1B[i][j] + R2B[i][j] + TB[i][j]; 
            B11 = B11 + R1B[i][j]; 
            B22 = B22 + R2B[i][j]; 
            BT = BT + TB[i][j]; 
        }
    }
    std::cout << "A:" << A << "\t" << "A+:" << A11 << "\t" << "A-:" << A22 << "\t" << "A tumble:" << AT << std::endl; 
    std::cout << "B:" << B << "\t" << "B+:" << B11 << "\t" << "B-:" << B22 << "\t" << "B tumble:" << BT << std::endl; 
    //////
    auto end_time = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> running_time = end_time - start_time; 
    std::cout << std::fixed << std::setprecision(4) << "Total running time: " << running_time.count() << "s" << std::endl; 

    return 0; 
}