enum class ParticleStates : int {
    Running_left_state = -1, 
    Running_right_state = 1, 
    Tumbling_0_state = 0 
}; 
void pbc(unsigned int Lx, unsigned int Ly, int &x, int &y) 
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
double cal_alpha_S(double kappa, double alpha_S0, double zeta_SA, double zeta_SB, unsigned int rho0, unsigned int rho_A, unsigned int rho_B) 
{
    double d_rho0 = static_cast<double>(rho0); 
    double d_rho_A = static_cast<double>(rho_A); 
    double d_rho_B = static_cast<double>(rho_B); 
    double alpha_S = alpha_S0 * (1 + kappa * tanh(zeta_SA * (d_rho_A - d_rho0) / (kappa * alpha_S0 * d_rho0))) * (1 + kappa * tanh(zeta_SB * (d_rho_B - d_rho0) / (kappa * alpha_S0 * d_rho0))); 
    return alpha_S; 
}
double cal_beta_S(double kappa, double beta_S0, double xi_SA, double xi_SB, unsigned int rho0, unsigned int rho_A, unsigned int rho_B) 
{
    double d_rho0 = static_cast<double>(rho0); 
    double d_rho_A = static_cast<double>(rho_A); 
    double d_rho_B = static_cast<double>(rho_B); 
    double beta_S = beta_S0 * (1 + kappa * tanh(xi_SA * (d_rho_A - d_rho0) / (kappa * beta_S0 * d_rho0))) * (1 + kappa * tanh(xi_SB * (d_rho_B - d_rho0) / (kappa * beta_S0 * d_rho0))); 
    return beta_S; 
}
void running_state(unsigned int Lx, unsigned int Ly, double kappa, double v, double Dt, double alpha_S0, double zeta_SA, double zeta_SB, unsigned int rho0, unsigned int rhoA_xy, unsigned int rhoB_xy, double r, int &x, int &y, ParticleStates &s) 
{
    double alpha_S = cal_alpha_S(kappa, alpha_S0, zeta_SA, zeta_SB, rho0, rhoA_xy, rhoB_xy); 
    if (r < alpha_S) 
    {
        x = x; 
        y = y; 
        s = ParticleStates::Tumbling_0_state; 
    }
    else if (r < (alpha_S + v + Dt)) //self-propelling 
    {
        if (s == ParticleStates::Running_right_state) 
        {
            x = x + 1; 
            y = y; 
            pbc(Lx, Ly, x, y); 
            s = ParticleStates::Running_right_state; 
        }
        else 
        {
            x = x - 1; 
            y = y; 
            pbc(Lx, Ly, x, y); 
            s = ParticleStates::Running_left_state; 
        }
    }
    else if (r < (alpha_S + v + 2*Dt)) 
    {
        if (s == ParticleStates::Running_right_state) 
        {
            x = x - 1; 
            y = y; 
            pbc(Lx, Ly, x, y); 
            s = ParticleStates::Running_right_state; 
        }
        else 
        {
            x = x + 1; 
            y = y; 
            pbc(Lx, Ly, x, y); 
            s = ParticleStates::Running_left_state; 
        }
    }
    else if (r < (alpha_S + v + 3*Dt)) 
    {
        x = x; 
        y = y + 1; 
        pbc(Lx, Ly, x, y); 
        s = s; 
    }
    else if (r < (alpha_S + v + 4*Dt)) 
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
void tumbling_state(unsigned int Lx, unsigned int Ly, double kappa, double Dt, double beta_S0, double xi_SA, double xi_SB, unsigned int rho0, unsigned int rhoA_xy, unsigned int rhoB_xy, double r, double r0, int &x, int &y, ParticleStates &s) 
{ 
    double beta_S = cal_beta_S(kappa, beta_S0, xi_SA, xi_SB, rho0, rhoA_xy, rhoB_xy); 
    if (r < beta_S) 
    {
        if (r0 < 0.50) 
        {
            x = x; 
            y = y; 
            s = ParticleStates::Running_left_state; 
        }
        else 
        {
            x = x; 
            y = y; 
            s = ParticleStates::Running_right_state; 
        }
    }
    else if (r < (beta_S + Dt)) 
    {
        x = x + 1; 
        y = y; 
        pbc(Lx, Ly, x, y); 
        s = s; 
    }
    else if (r < (beta_S + 2*Dt)) 
    {
        x = x - 1; 
        y = y; 
        pbc(Lx, Ly, x, y); 
        s = s; 
    }
    else if (r < (beta_S + 3*Dt)) 
    {
        x = x; 
        y = y + 1; 
        pbc(Lx, Ly, x, y); 
        s = s; 
    }
    else if (r < (beta_S + 4*Dt)) 
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