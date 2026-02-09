#include <iostream> 
#include "domain2D.h"
#include "cellList2D.h"
#include "2S_rt_2d_2pi.h"
#include "exporter2D.h"

int main(int argc, char* argv[]) {//std::stoi(argv[1]), std::stod(argv[1])
    double Lx = 10. * 2;
    double Ly = 5. * 4;
    double phiA = 100. * 1;
    double phiB = 100. * 1;
    double rho0 = 100. * 1;
    double h = 0.01;
    unsigned int n_step = 500000;
    //
    double omegaAA = 9.0/10;// + 3;
    double omegaAB = -3.; //omegaAB = 0, only particle A
    double omegaBA = -omegaAB;
    double omegaBB = omegaAA;
    double kappa = 1;
    double Dt = 0.0;
    double alphaA0 = 0.05;
    unsigned int seed = 2026;
    std::string ini_mode = "new"; // should be "new" or "restart"
    //
    double v0 = 1;
    double sigma_A = 1.0/9;
    double betaA0 = alphaA0 / sigma_A;
    double xiAA = 0.0; // xi = 0 means beta = beta0
    double xiAB = 0.0;
    double sigma_B = sigma_A;
    double alphaB0 = alphaA0;
    double betaB0 = alphaB0 / sigma_B;
    double xiBA = 0.0;
    double xiBB = 0.0;
    //
    double etaAA = 1/(1 + sigma_A) - omegaAA; 
    double etaAB = -omegaAB; 
    double etaBA = -omegaBA; 
    double etaBB = 1/(1 + sigma_B) - omegaBB; 
    double zetaAA = (etaAA * (alphaA0 + betaA0) * (alphaA0 + betaA0) + xiAA * alphaA0) / betaA0; 
    double zetaAB = (etaAB * (alphaA0 + betaA0) * (alphaA0 + betaA0) + xiAB * alphaA0) / betaA0; 
    double zetaBA = (etaBA * (alphaB0 + betaB0) * (alphaB0 + betaB0) + xiBA * alphaB0) / betaB0; 
    double zetaBB = (etaBB * (alphaB0 + betaB0) * (alphaB0 + betaB0) + xiBB * alphaB0) / betaB0; 
    unsigned int n_par_A = static_cast<unsigned int>(std::round(Lx * Ly * phiA));
    unsigned int n_par_A_half = n_par_A / 2;
    unsigned int n_par_B = static_cast<unsigned int>(std::round(Lx * Ly * phiB));
    unsigned int n_par_B_half = n_par_B / 2;
    unsigned int n_par = n_par_A + n_par_B;

    typedef BiNode<RTP_QS_2pi> node_t;
    Ranq2 myran(static_cast<unsigned int>(seed));
    Vec_2<double> gl_l(Lx, Ly);
    double r_cut = 1;
    Grid_2 grid(gl_l, r_cut);
    PeriodicDomain_2 pdm(gl_l);
    CellListNode_2<node_t> cl(pdm, grid);
    std::vector<node_t> p_arr;

    // ini integrator
    EM_RTP_QS_2pi integrator(h, v0, Dt);
    // ini rate calculator
    Trans_Rates rates(kappa, rho0,
        alphaA0, alphaB0, zetaAA, zetaAB, zetaBA, zetaBB,
        betaA0, betaB0, xiAA, xiAB, xiBA, xiBB);
    // cal force
    LinearDensityKernal kernal(r_cut);
    auto f1 = [&kernal](node_t* p1, node_t* p2) {
        kernal(*p1, *p2);
    };
    auto f2 = [&kernal, &pdm](node_t* p1, node_t* p2) {
        kernal(*p1, *p2, pdm);
    };

    //set output
    char basename[255];
    char log_file[255];
    char gsd_file[255];
    // char folder[] = "D:/Adata/VS code/C/202511RTtwoSpecies/2S_rt_2d_2pi/00/";
    char folder[] = "/mnt/d/Adata/VS code/C/202511RTtwoSpecies/2S_rt_2d_2pi/data11/";
    std::snprintf(basename, 255, "L%g_%g_p%g_%g_r%g_omAA%.3f_AB%.3f_BA%.3f_BB%.3f_k%g_Dt%g_al%g_s%d",
        Lx, Ly, phiA, phiB, rho0,
        omegaAA, omegaAB, omegaBA, omegaBB,
        kappa, Dt, alphaA0, seed);
    std::snprintf(log_file, 255, "%slog_%s.dat", folder, basename);
    std::snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);

    unsigned int snap_interval = 1000;
    unsigned int log_interval = 1000;
    exporter::LogExporter log(log_file, 0, n_step, log_interval, n_par);
    exporter::Snap_GSD_2 gsd(gsd_file, 0, n_step, snap_interval, gl_l, ini_mode);

    //ini particles
    p_arr.reserve(n_par);
    if (ini_mode == "new") {
        for (unsigned int i = 0; i < n_par; i++) {//or starting from i = n_par_A_half
            p_arr.emplace_back(myran, gl_l, Vec_2<double>());
            if (i < n_par_A_half) {
                p_arr[i].type_id = 0; //A species, R
            } else if (i < n_par_A) {
                p_arr[i].type_id = 2; //A species, T !!!
            } else if (i < n_par_A + n_par_B_half) {
                p_arr[i].type_id = 1; //B species, R !!!
            } else {
                p_arr[i].type_id = 3; //B species, T
            }
        }
    } else {
        gsd.read_last_frame(p_arr);
    }

    cl.create(p_arr);

    for (unsigned int t = 0; t < n_step; t++) {
        cl.for_each_pair(f1, f2);
        kernal.normalize(p_arr);
        for (unsigned int i = 0; i < n_par; i++) {
            double prob = myran.doub();
            if (p_arr[i].type_id < 2) {
                double alpha = rates.cal_alpha_S(p_arr[i]);
                if (prob < alpha * h) {//alpha * h ??? as h -> 0
                    update_type_id(p_arr[i], 0., 0., 0.); //parameters prob, alpha(beta) and h are not used in this case
                }
                integrator.update_pos(p_arr[i], pdm, myran);
            } else {
                double beta = rates.cal_beta_S(p_arr[i]);
                if (prob < beta * h) {
                    update_type_id(p_arr[i], prob, beta, h);
                }
                integrator.update_pos(p_arr[i], pdm, myran);
            }
        }
        // std::cout << p_arr[0].rho_local[0] << std::endl;
        // std::cout << p_arr[0].rho_local[1] << std::endl;
        kernal.reset_local_density(p_arr);
        cl.recreate(p_arr);
        gsd.dump(t, p_arr);
        log.record(t);
    }

    //
    // std::cout << p_arr[0].pos.x << ", " << p_arr[0].pos.y << std::endl;
    // std::cout << p_arr[0].u.x << ", " << p_arr[0].u.y << std::endl;
    // std::cout << p_arr[0].type_id << std::endl;

    return 0; 
}