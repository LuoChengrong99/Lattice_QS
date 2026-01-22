#include "particle_state.h" 
void pbc(const unsigned int Lx, const unsigned int Ly, int &x, int &y) 
{
    if (x < 0) {
        x += Lx; 
    } else if (x >= Lx) {
        x -= Lx; 
    }
    if (y < 0) {
        y += Ly; 
    } else if (y >= Ly) {
        y -= Ly; 
    }
}
double cal_alpha_S(const double kappa, const double alpha_S0, const double zeta_SA, const double zeta_SB, 
    const unsigned int rho0, const unsigned int rho_A, const unsigned int rho_B
) 
{
    double d_rho0 = static_cast<double>(rho0); 
    double d_rho_A = static_cast<double>(rho_A); 
    double d_rho_B = static_cast<double>(rho_B); 
    double alpha_S = alpha_S0 * (1 + kappa * std::tanh(zeta_SA * (d_rho_A - d_rho0) / (kappa * alpha_S0 * d_rho0))) * (1 + kappa * std::tanh(zeta_SB * (d_rho_B - d_rho0) / (kappa * alpha_S0 * d_rho0))); 
    return alpha_S; 
}
double cal_beta_S(const double kappa, const double beta_S0, const double xi_SA, const double xi_SB, 
    const unsigned int rho0, const unsigned int rho_A, const unsigned int rho_B
) 
{
    double d_rho0 = static_cast<double>(rho0); 
    double d_rho_A = static_cast<double>(rho_A); 
    double d_rho_B = static_cast<double>(rho_B); 
    double beta_S = beta_S0 * (1 + kappa * std::tanh(xi_SA * (d_rho_A - d_rho0) / (kappa * beta_S0 * d_rho0))) * (1 + kappa * std::tanh(xi_SB * (d_rho_B - d_rho0) / (kappa * beta_S0 * d_rho0))); 
    return beta_S; 
}
void running_state(const unsigned int Lx, const unsigned int Ly, const double kappa, const double v, const double Dt, 
    const double alpha_S0, const double zeta_SA, const double zeta_SB, 
    const unsigned int rho0, const unsigned int rhoA_xy, const unsigned int rhoB_xy, const double r, 
    int &x, int &y, ParticleStates &s
) 
{
    double alpha_S = cal_alpha_S(kappa, alpha_S0, zeta_SA, zeta_SB, 
        rho0, rhoA_xy, rhoB_xy
    ); 
    if (r < alpha_S) 
    {
        x = x; 
        y = y; 
        s = ParticleStates::Tumbling_0_state; 
    }
    else if (r < (alpha_S + Dt)) //R 
    {
        x = x + 1; 
        y = y; 
        pbc(Lx, Ly, x, y); 
        s = s; 
    }
    else if (r < (alpha_S + 2*Dt)) //U 
    {
        x = x; 
        y = y + 1; 
        pbc(Lx, Ly, x, y); 
        s = s; 
    }
    else if (r < (alpha_S + 3*Dt)) //L 
    {
        x = x - 1; 
        y = y; 
        pbc(Lx, Ly, x, y); 
        s = s; 
    }
    else if (r < (alpha_S + 4*Dt)) //D 
    {
        x = x; 
        y = y - 1; 
        pbc(Lx, Ly, x, y); 
        s = s; 
    }
    else if (r < (alpha_S + 4*Dt + v)) //self-propelling 
    {
        if (s == ParticleStates::Running_1_state) //R 
        {
            x = x + 1; 
            y = y; 
            pbc(Lx, Ly, x, y); 
            s = s; 
        }
        else if (s == ParticleStates::Running_2_state) //U 
        {
            x = x; 
            y = y + 1; 
            pbc(Lx, Ly, x, y); 
            s = s; 
        }
        else if (s == ParticleStates::Running_3_state) //L 
        {
            x = x - 1; 
            y = y; 
            pbc(Lx, Ly, x, y); 
            s = s; 
        }
        else //D 
        {
            x = x; 
            y = y - 1; 
            pbc(Lx, Ly, x, y); 
            s = s; 
        }
    }
    else 
    {
        x = x; 
        y = y; 
        s = s; 
    }

    return; 
}
void tumbling_state(const unsigned int Lx, const unsigned int Ly, const double kappa, const double Dt, 
    const double beta_S0, const double xi_SA, const double xi_SB, 
    const unsigned int rho0, const unsigned int rhoA_xy, const unsigned int rhoB_xy, const double r, 
    int &x, int &y, ParticleStates &s
) 
{ 
    double beta_S = cal_beta_S(kappa, beta_S0, xi_SA, xi_SB, 
        rho0, rhoA_xy, rhoB_xy
    ); 
    if (r < beta_S) 
    {
        x = x; 
        y = y; 
        s = static_cast<ParticleStates>(int(4 * r / beta_S) + 1); 
    }
    else if (r < (beta_S + Dt)) //R 
    {
        x = x + 1; 
        y = y; 
        pbc(Lx, Ly, x, y); 
        s = s; 
    }
    else if (r < (beta_S + 2*Dt)) //U 
    {
        x = x; 
        y = y + 1; 
        pbc(Lx, Ly, x, y); 
        s = s; 
    }
    else if (r < (beta_S + 3*Dt)) //L 
    {
        x = x - 1; 
        y = y; 
        pbc(Lx, Ly, x, y); 
        s = s; 
    }
    else if (r < (beta_S + 4*Dt)) //D 
    {
        x = x; 
        y = y - 1; 
        pbc(Lx, Ly, x, y); 
        s = s; 
    }
    else 
    {
        x = x; 
        y = y; 
        s = s; 
    }

    return; 
}