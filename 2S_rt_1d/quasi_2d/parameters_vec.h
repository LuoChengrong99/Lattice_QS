struct Par_i {
    char name[64]; 
    unsigned int value; 
}; 
struct Par_f {
    char name[64]; 
    double value; 
}; 
void parameters_vec_i(
    unsigned int Lx, unsigned int Ly, 
    unsigned int N, unsigned int NA, unsigned int NB, unsigned int n, 
    std::vector<Par_i>& par_vec_i
) 
{
    par_vec_i.emplace_back(); 
    std::strncpy(par_vec_i.back().name, "Lx", 63); 
    par_vec_i.back().value = Lx; 
    par_vec_i.emplace_back(); 
    std::strncpy(par_vec_i.back().name, "Ly", 63); 
    par_vec_i.back().value = Ly; 
    par_vec_i.emplace_back(); 
    std::strncpy(par_vec_i.back().name, "N", 63); 
    par_vec_i.back().value = N; 
    par_vec_i.emplace_back(); 
    std::strncpy(par_vec_i.back().name, "NA", 63); 
    par_vec_i.back().value = NA; 
    par_vec_i.emplace_back(); 
    std::strncpy(par_vec_i.back().name, "NB", 63); 
    par_vec_i.back().value = NB; 
    par_vec_i.emplace_back(); 
    std::strncpy(par_vec_i.back().name, "n", 63); 
    par_vec_i.back().value = n; 
}
void parameters_vec_f(
    double rho0, double rhoA0, double rhoB0, 
    double omegaAA, double omegaAB, double omegaBA, double omegaBB, double kappa, double v, double Dt, 
    double sigma_A, double alphaA0, double zetaAA, double zetaAB, double betaA0, double xiAA, double xiAB, 
    double sigma_B, double alphaB0, double zetaBA, double zetaBB, double betaB0, double xiBA, double xiBB, 
    double etaAA, double etaAB, double etaBA, double etaBB, 
    std::vector<Par_f>& par_vec_f
) 
{
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "rho0", 63); 
    par_vec_f.back().value = rho0; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "rhoA0", 63); 
    par_vec_f.back().value = rhoA0; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "rhoB0", 63); 
    par_vec_f.back().value = rhoB0; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "omegaAA", 63); 
    par_vec_f.back().value = omegaAA; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "omegaAB", 63); 
    par_vec_f.back().value = omegaAB; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "omegaBA", 63); 
    par_vec_f.back().value = omegaBA; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "omegaBB", 63); 
    par_vec_f.back().value = omegaBB; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "kappa", 63); 
    par_vec_f.back().value = kappa; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "v", 63); 
    par_vec_f.back().value = v; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "Dt", 63); 
    par_vec_f.back().value = Dt; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "sigmaA", 63); 
    par_vec_f.back().value = sigma_A; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "alphaA0", 63); 
    par_vec_f.back().value = alphaA0; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "zetaAA", 63); 
    par_vec_f.back().value = zetaAA; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "zetaAB", 63); 
    par_vec_f.back().value = zetaAB; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "betaA0", 63); 
    par_vec_f.back().value = betaA0; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "xiAA", 63); 
    par_vec_f.back().value = xiAA; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "xiAB", 63); 
    par_vec_f.back().value = xiAB; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "sigmaB", 63); 
    par_vec_f.back().value = sigma_B; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "alphaB0", 63); 
    par_vec_f.back().value = alphaB0; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "zetaBA", 63); 
    par_vec_f.back().value = zetaBA; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "zetaBB", 63); 
    par_vec_f.back().value = zetaBB; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "betaB0", 63); 
    par_vec_f.back().value = betaB0; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "xiBA", 63); 
    par_vec_f.back().value = xiBA; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "xiBB", 63); 
    par_vec_f.back().value = xiBB; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "etaAA", 63); 
    par_vec_f.back().value = etaAA; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "etaAB", 63); 
    par_vec_f.back().value = etaAB; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "etaBA", 63); 
    par_vec_f.back().value = etaBA; 
    par_vec_f.emplace_back(); 
    std::strncpy(par_vec_f.back().name, "etaBB", 63); 
    par_vec_f.back().value = etaBB; 
}