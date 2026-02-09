#pragma once
#include <cstdint>
#include <cmath>
#include "vect.h"
#define M_PI 3.14159265358979323846

class RTP_QS_2pi {
    public:
    RTP_QS_2pi() : pos(), u(), type_id(0) {}
    RTP_QS_2pi(const Vec_2<double>& pos0) : pos(pos0), u(), type_id(0) {}
    RTP_QS_2pi(const Vec_2<double>& pos0, const Vec_2<double>& u0)
        : pos(pos0), u(u0), type_id(0) {}
    template <typename TRan>
    RTP_QS_2pi(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o);

    double get_theta() const { return std::atan2(u.y, u.x); }

    Vec_2<double> pos;
    Vec_2<double> u;
    uint32_t type_id;
    float rho_local[2]{}; //index are calculated as p.type_id % 2: 0 for A species, 1 for B species
};

template <typename TRan>
RTP_QS_2pi::RTP_QS_2pi(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o)
: pos(o.x + myran.doub() * l.x, o.y + myran.doub() * l.y) {
    double theta = myran.doub() * M_PI * 2.0;
    // double theta;//theta = 0 or pi
    // if (myran.doub() < 0.50) {
    //     theta = 0.0;
    // } else {
    //     theta = M_PI;
    // }//
    u.x = std::cos(theta);
    u.y = std::sin(theta);
    type_id = 0;
}

class EM_RTP_QS_2pi {
    public:
    EM_RTP_QS_2pi(double h, double v0, double Dt)
     : h_(h), v0_(v0), Dt_(std::sqrt(24 * Dt *h)), trans_noise_on_(Dt > 0.0) {}
    
    template <typename TPar, class BdyCondi, class TRan>
    void update_pos(TPar& p, const BdyCondi& bc, TRan& myran) const;

    protected:
    double h_;
    double v0_;
    double Dt_;
    bool trans_noise_on_;
};

template <typename TPar, class BdyCondi, class TRan>
void EM_RTP_QS_2pi::update_pos(TPar& p, const BdyCondi& bc, TRan& myran) const {
    if (p.type_id == 2 || p.type_id == 3) {
        p.pos.x = p.pos.x;
        p.pos.y = p.pos.y;
    } else {
        double displace = h_ * v0_;
        p.pos.x += p.u.x * displace;
        p.pos.y += p.u.y * displace;
    }
    if (trans_noise_on_) {
        p.pos.x += (myran.doub() - 0.5) * Dt_;
        p.pos.y += (myran.doub() - 0.5) * Dt_;
    }
    bc.tangle(p.pos);
}

struct Trans_Rates {
    double kappa_, rho0_;
    double alpha_A0_, alpha_B0_, zeta_AA_, zeta_AB_, zeta_BA_, zeta_BB_;
    double beta_A0_, beta_B0_, xi_AA_, xi_AB_, xi_BA_, xi_BB_;
    bool beta_dependent_on_density_;
    Trans_Rates(
        double kappa, double rho0,
        double alpha_A0, double alpha_B0, double zeta_AA, double zeta_AB, double zeta_BA, double zeta_BB,
        double beta_A0, double beta_B0, double xi_AA, double xi_AB, double xi_BA, double xi_BB
    ) : kappa_(kappa), rho0_(rho0),
    alpha_A0_(alpha_A0), alpha_B0_(alpha_B0), zeta_AA_(zeta_AA), zeta_AB_(zeta_AB), zeta_BA_(zeta_BA), zeta_BB_(zeta_BB),
    beta_A0_(beta_A0), beta_B0_(beta_B0), xi_AA_(xi_AA), xi_AB_(xi_AB), xi_BA_(xi_BA), xi_BB_(xi_BB),
    beta_dependent_on_density_(xi_AA != 0.0 || xi_AB != 0.0 || xi_BA != 0.0 || xi_BB != 0.0) {}

    template <typename TPar>
    double cal_alpha_S(const TPar& p) const;
    template <typename TPar>
    double cal_beta_S(const TPar& p) const;
};

template<typename TPar>
double Trans_Rates::cal_alpha_S(const TPar& p) const {
    double alpha_S;
    double d_rhoA = (p.rho_local[0] - rho0_) / (rho0_ * kappa_);
    double d_rhoB = (p.rho_local[1] - rho0_) / (rho0_ * kappa_);
    if (p.type_id == 0) {
        alpha_S = alpha_A0_ * (1 + kappa_ * std::tanh(zeta_AA_ * d_rhoA / alpha_A0_)) * (1 + kappa_ * std::tanh(zeta_AB_ * d_rhoB / alpha_A0_));
    } else {//type_id == 1
        alpha_S = alpha_B0_ * (1 + kappa_ * std::tanh(zeta_BA_ * d_rhoA / alpha_B0_)) * (1 + kappa_ * std::tanh(zeta_BB_ * d_rhoB / alpha_B0_));
    }
    return alpha_S;
}

template <typename TPar>
double Trans_Rates::cal_beta_S(const TPar& p) const {
    double beta_S;
    if (beta_dependent_on_density_) {
        double d_rhoA = (p.rho_local[0] - rho0_) / (rho0_ * kappa_);
        double d_rhoB = (p.rho_local[1] - rho0_) / (rho0_ * kappa_);
        if (p.type_id == 2) {
            beta_S = beta_A0_ * (1 + kappa_ * std::tanh(xi_AA_ * d_rhoA / beta_A0_)) * (1 + kappa_ * std::tanh(xi_AB_ * d_rhoB / beta_A0_));
        } else {//type_id == 3
            beta_S = beta_B0_ * (1 + kappa_ * std::tanh(xi_BA_ * d_rhoA / beta_B0_)) * (1 + kappa_ * std::tanh(xi_BB_ * d_rhoB / beta_B0_));
        }
    } else {
        if (p.type_id == 2) {
            beta_S = beta_A0_;
        } else {//type_id == 3
            beta_S = beta_B0_;
        }
    }
    return beta_S;
}

template <typename TPar>
void update_type_id(TPar& p, const double prob, const double beta, const double h) {
    if (p.type_id < 2) {
        p.type_id += 2; //from R to T
    } else {
        double theta = (prob / (beta * h)) * M_PI * 2.0;
        // double theta;//theta = 0 or pi
        // if ((prob / (beta * h)) < 0.50) {
        //     theta = 0;
        // } else {
        //     theta = M_PI;
        // }//
        p.u.x = std::cos(theta);
        p.u.y = std::sin(theta);
        p.type_id -= 2; //from T to R
    }
}

class LinearDensityKernal {
public:
  LinearDensityKernal(double r_cut) : r_cut_(r_cut), r_cut_square_(r_cut* r_cut) {}

  template <typename TPar>
  void accum_density(double r12, TPar& p1, TPar& p2) const;

  template <typename TPar>
  void operator ()(TPar& p1, TPar& p2) const;

  template <typename TPar, typename BdyCondi>
  void operator ()(TPar& p1, TPar& p2, const BdyCondi& bc) const;

  template <typename TPar>
  void normalize(std::vector<TPar>& p_arr) const;

  template <typename TPar>
  void reset_local_density(std::vector<TPar>& p_arr) const;

protected:
  double r_cut_;
  double r_cut_square_;
};

template<typename TPar>
void LinearDensityKernal::accum_density(double r12, TPar& p1, TPar& p2) const {
  double density_weighted = r_cut_ - r12;
//   double density_weighted = 1;//const weighted
//   p1.rho_local[p2.type_id] += density_weighted;
//   p2.rho_local[p1.type_id] += density_weighted;
  p1.rho_local[p2.type_id % 2] += density_weighted; //type_id 0,2 for A species; 1,3 for B species
  p2.rho_local[p1.type_id % 2] += density_weighted;

  //!TODO need divided by a nomalized factor
}

template <typename TPar>
void LinearDensityKernal::operator ()(TPar& p1, TPar& p2) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  const double r12_square = r12.square();
  if (r12_square < r_cut_square_) {
    double r12_module = sqrt(r12_square);
    accum_density(r12_module, p1, p2);
  }
}

template <typename TPar, typename BdyCondi>
void LinearDensityKernal::operator ()(TPar& p1, TPar& p2, const BdyCondi& bc) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  bc.untangle(r12);
  const double r12_square = r12.square();
  if (r12_square < r_cut_square_) {
    double r12_module = sqrt(r12_square);
    accum_density(r12_module, p1, p2);
  }
}

template <typename TPar>
void LinearDensityKernal::normalize(std::vector<TPar>& p_arr) const {
  double norm = 3. / M_PI;
//   double norm = 1. / M_PI;
  for (auto& p : p_arr) {
    //p.rho_local[p.type_id] += 1;
    p.rho_local[0] *= norm;
    p.rho_local[1] *= norm;
  }
}

template <typename TPar>
void LinearDensityKernal::reset_local_density(std::vector<TPar>& p_arr) const {
  for (auto& p : p_arr) {
    p.rho_local[0] = 0.;
    p.rho_local[1] = 0.;
  }
}
