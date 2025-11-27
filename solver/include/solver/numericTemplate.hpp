#pragma once

#include <cmath>

#include "pre/parameter.h"
#include "solver/variableDef.h"
namespace solver {

inline PRIM_VAR ConvertFromConv(const CONS_VAR &var, double cp, double gamma) {
  double gam1 = gamma - 1.0;
  double rgas = gam1 * cp / gamma;
  double g1cp = gam1 * cp;
  PRIM_VAR result;
  double rhoq = var.xmom * var.xmom + var.ymom * var.ymom;
  result.dens = var.dens;
  result.velx = var.xmom / var.dens;
  result.vely = var.ymom / var.dens;
  result.press = gam1 * (var.ener - 0.5 * rhoq / var.dens);
  return result;
}

inline VISC_VAR ConvertToVISC(const CONS_VAR &cv, const double &gamma,
                              const double &cp, const double &prandtl,
                              const double &refVisc) {
  double gam1 = gamma - 1.0;
  double rgas = gam1 * cp / gamma;
  double g1cp = gam1 * cp;
  double s1 = 110.0;
  double s2 = 288.16;
  double s12 = 1.0 + s1 / s2;
  double cppr = cp / prandtl;
  double rhoq = cv.xmom * cv.xmom + cv.ymom * cv.ymom;
  DEPEND_VAR result;
  result.press = gam1 * (cv.ener - 0.5 * rhoq / cv.dens);
  result.temp = result.press / (rgas * cv.dens);
  result.cs = std::sqrt(g1cp * result.temp);
  result.gamma = gamma;
  result.Cp = cp;
  double rat = std::sqrt(result.temp / s2) * s12 / (1.0 + s1 / result.temp);
  double mu = refVisc * rat;
  double lambda = mu * cppr;
  return {mu, lambda};
}

struct ROENumeric {
  inline static void ComputeResidual(const preprocess::parameter &param,
                                     double nx, double ny, double ds,
                                     PRIM_VAR vari, PRIM_VAR varj,
                                     double gammai, double gammaj,
                                     CONS_VAR &resi, CONS_VAR &resj) {
    double rl = vari.dens;
    double ul = vari.velx;
    double vl = vari.vely;
    double pl = vari.press;
    double hl = gammai / (gammai - 1.0) * pl / rl + 0.5 * (ul * ul + vl * vl);

    double rr = varj.dens;
    double ur = varj.velx;
    double vr = varj.vely;
    double pr = varj.press;
    double hr = gammaj / (gammaj - 1.0) * pr / rr + 0.5 * (ur * ur + vr * vr);

    double rav = std::sqrt(rl * rr);
    double gam1 = 0.5 * (gammai + gammaj) - 1.0;
    double dd = rav / rl;
    double dd1 = 1.0 / (1.0 + dd);
    double uav = (ul + dd * ur) * dd1;
    double vav = (vl + dd * vr) * dd1;
    double hav = (hl + dd * hr) * dd1;
    double q2a = 0.5 * (uav * uav + vav * vav);
    double c2a = gam1 * (hav - q2a);
    double cav = std::sqrt(c2a);
    double uv = uav * nx + vav * ny;
    double du = (ur - ul) * nx + (vr - vl) * ny;

    double h1 = std::abs(uv - cav);
    double h2 = std::abs(uv);
    double h4 = std::abs(uv + cav);
    double delta = param.entropyCorreRoe * h4;

    auto entropyCorr = [](const double &z, const double &d) {
      if (z > d) {
        return z;
      } else {
        return 0.5 * (z * z + d * d) / d;
      }
    };

    double eabs1 = entropyCorr(h1, delta);
    double eabs2 = entropyCorr(h2, delta);
    double eabs4 = entropyCorr(h4, delta);

    h1 = rav * cav * du;
    h2 = eabs1 * (pr - pl - h1) / (2.0 * c2a);
    double h3 = eabs2 * (rr - rl - (pr - pl) / c2a);
    h4 = eabs2 * rav;
    double h5 = eabs4 * (pr - pl + h1) / (2.0 * c2a);

    double fd[4];

    fd[0] = h2 + h3 + h5;
    fd[1] = h2 * (uav - cav * nx) + h3 * uav + h4 * (ur - ul - du * nx) +
            h5 * (uav + cav * nx);
    fd[2] = h2 * (vav - cav * ny) + h3 * vav + h4 * (vr - vl - du * ny) +
            h5 * (vav + cav * ny);
    fd[3] = h2 * (hav - cav * uv) + h3 * q2a +
            h4 * (uav * (ur - ul) + vav * (vr - vl) - uv * du) +
            h5 * (hav + cav * uv);

    double qsrl = (ul * nx * ds + vl * ny * ds) * rl;
    double qsrr = (ur * nx * ds + vr * ny * ds) * rr;

    double pav = 0.5 * (pl + pr);

    double fc[4];

    fc[0] = 0.5 * (qsrl + qsrr);
    fc[1] = 0.5 * (qsrl * ul + qsrr * ur) + nx * ds * pav;
    fc[2] = 0.5 * (qsrl * vl + qsrr * vr) + ny * ds * pav;
    fc[3] = 0.5 * (qsrl * hl + qsrr * hr);

    resi.dens += (fc[0] - 0.5 * fd[0] * ds);
    resi.xmom += (fc[1] - 0.5 * fd[1] * ds);
    resi.ymom += (fc[2] - 0.5 * fd[2] * ds);
    resi.ener += (fc[3] - 0.5 * fd[3] * ds);

    resj.dens -= (fc[0] - 0.5 * fd[0] * ds);
    resj.xmom -= (fc[1] - 0.5 * fd[1] * ds);
    resj.ymom -= (fc[2] - 0.5 * fd[2] * ds);
    resj.ener -= (fc[3] - 0.5 * fd[3] * ds);
  }
};

struct SLAU2Numeric {
  inline static double pressureFuncLeft(const double &m) {
    if (std::abs(m) >= 1.0) {
      return 0.5 * (1 + (m >= 0 ? 1 : -1));
    } else {
      return 0.25 * (m + 1) * (m + 1) * (2 - m);
    }
  }

  inline static double pressureFuncRight(const double &m) {
    if (std::abs(m) >= 1.0) {
      return 0.5 * (1 - (m >= 0 ? 1 : -1));
    } else {
      return 0.25 * (m - 1) * (m - 1) * (2 + m);
    }
  }

  inline static void ComputeResidual(const preprocess::parameter &param,
                                     double nx, double ny, double ds,
                                     PRIM_VAR vari, PRIM_VAR varj,
                                     double gammai, double gammaj,
                                     CONS_VAR &resi, CONS_VAR &resj) {
    double rl = vari.dens;
    double ul = vari.velx;
    double vl = vari.vely;
    double pl = vari.press;
    double hl = gammai / (gammai - 1.0) * pl / rl + 0.5 * (ul * ul + vl * vl);

    double rr = varj.dens;
    double ur = varj.velx;
    double vr = varj.vely;
    double pr = varj.press;
    double hr = gammaj / (gammaj - 1.0) * pr / rr + 0.5 * (ur * ur + vr * vr);

    double cl = (gammai - 1.0) * (hl - 0.5 * (ul * ul + vl * vl));
    double cr = (gammaj - 1.0) * (hr - 0.5 * (ur * ur + vr * vr));

    double veln_left = ul * nx + vl * ny;
    double veln_right = ur * nx + vr * ny;
    double veln_bar_abs =
        (rl * std::abs(veln_left) + rr * std::abs(veln_right)) / (rl + rr);
    double c_interface = 0.5 * (cl + cr);
    double m_left = veln_left / c_interface;
    double m_right = veln_right / c_interface;
    double g = -1.0 * std::max(std::min(m_left, 0.0), -1.0) *
               std::min(std::max(m_right, 0.0), 1.0);
    double veln_bar_abs_p = (1.0 - g) * veln_bar_abs + g * std::abs(veln_left);
    double veln_bar_abs_n = (1.0 - g) * veln_bar_abs + g * std::abs(veln_right);

    double m_hat = std::min(
        1.0, 1 / c_interface *
                 std::sqrt(0.5 * (ul * ul + vl * vl + ur * ur + vr * vr)));
    double X = std::pow((1.0 - m_hat), 2);

    double fp_left = pressureFuncLeft(m_left);
    double fp_right = pressureFuncRight(m_right);
    double p_interface =
        0.5 * (pl + pr) + 0.5 * (fp_left - fp_right) * (pl - pr) +
        std::sqrt(0.5 * (ul * ul + vl * vl + ur * ur + vr * vr)) *
            (fp_left + fp_right - 1) * c_interface * 0.5 * (rr + rl);
    double mass_interface = 0.5 * (rl * (veln_left + veln_bar_abs_p) +
                                   rr * (veln_right - veln_bar_abs_n) -
                                   X / c_interface * (pr - pl));
    double mass_flux_left = 0.5 * (mass_interface + std::abs(mass_interface));
    double mass_flux_right = 0.5 * (mass_interface - std::abs(mass_interface));

    double fc[4];

    fc[0] = (mass_flux_left + mass_flux_right) * ds;
    fc[1] =
        (mass_flux_left * ul + mass_flux_right * ur + p_interface * nx) * ds;
    fc[2] =
        (mass_flux_left * vl + mass_flux_right * vr + p_interface * ny) * ds;
    fc[3] = (mass_flux_left * hl + mass_flux_right * hr) * ds;

    resi.dens += fc[0];
    resi.xmom += fc[1];
    resi.ymom += fc[2];
    resi.ener += fc[3];

    resj.dens -= fc[0];
    resj.xmom -= fc[1];
    resj.ymom -= fc[2];
    resj.ener -= fc[3];
  }
};

struct AUSMNumeric {
  inline static void ComputeResidual(const preprocess::parameter &param,
                                     double nx, double ny, double ds,
                                     PRIM_VAR vari, PRIM_VAR varj,
                                     double gammai, double gammaj,
                                     CONS_VAR &resi, CONS_VAR &resj) {
    double rl = vari.dens;
    double ul = vari.velx;
    double vl = vari.vely;
    double pl = vari.press;
    double hl = gammai / (gammai - 1.0) * pl / rl + 0.5 * (ul * ul + vl * vl);

    double rr = varj.dens;
    double ur = varj.velx;
    double vr = varj.vely;
    double pr = varj.press;
    double hr = gammaj / (gammaj - 1.0) * pr / rr + 0.5 * (ur * ur + vr * vr);
    double cl2 = (gammai - 1.0) * (hl - 0.5 * (ul * ul + vl * vl));
    double cl = std::sqrt(cl2);
    double cr2 = (gammaj - 1.0) * (hr - 0.5 * (ur * ur + vr * vr));
    double cr = std::sqrt(cr2);

    double veln_left = ul * nx + vl * ny;
    double veln_right = ur * nx + vr * ny;

    double m_left = veln_left / cl;
    double m_right = veln_right / cr;

    double m_left_p;
    double p_left_p;
    if (m_left >= 1.0) {
      m_left_p = m_left;
      p_left_p = pl;
    } else if (std::abs(m_left) < 1.0) {
      m_left_p = 0.25 * (m_left + 1.0) * (m_left + 1.0);
      p_left_p = 0.25 * pl * (m_left + 1.0) * (m_left + 1.0) * (2.0 - m_left);
    } else if (m_left <= -1.0) {
      m_left_p = 0.0;
      p_left_p = 0.0;
    }

    double m_right_n;
    double p_right_n;
    if (m_right >= 1.0) {
      m_right_n = 0.0;
      p_right_n = 0.0;
    } else if (std::abs(m_right) < 1.0) {
      m_right_n = -0.25 * (m_right - 1.0) * (m_right - 1.0);
      p_right_n =
          0.25 * pr * (m_right - 1.0) * (m_right - 1.0) * (2.0 + m_right);
    } else if (m_right <= -1.0) {
      m_right_n = m_right;
      p_right_n = pr;
    }

    double m_interface = m_left_p + m_right_n;
    double p_interface = p_left_p + p_right_n;
    double m_interface_abs;
    double delta = 0.1;
    if (std::abs(m_interface) > delta) {
      m_interface_abs = std::abs(m_interface);
    } else {
      m_interface_abs =
          (m_interface * m_interface + delta * delta) / (2.0 * delta);
    }

    double fc[4];

    fc[0] = 0.5 * m_interface * (rl * cl + rr * cr) * ds -
            0.5 * m_interface_abs * (rr * cr - rl * cl) * ds;
    fc[1] = 0.5 * m_interface * (rl * cl * ul + rr * cr * ur) * ds -
            0.5 * m_interface_abs * (rr * cr * ur - rl * cl * ul) * ds +
            p_interface * nx * ds;
    fc[2] = 0.5 * m_interface * (rl * cl * vl + rr * cr * vr) * ds -
            0.5 * m_interface_abs * (rr * cr * vr - rl * cl * vl) * ds +
            p_interface * ny * ds;
    fc[3] = 0.5 * m_interface * (rl * cl * hl + rr * cr * hr) * ds -
            0.5 * m_interface_abs * (rr * cr * hr - rl * cl * hl) * ds;

    resi.dens += fc[0];
    resi.xmom += fc[1];
    resi.ymom += fc[2];
    resi.ener += fc[3];

    resj.dens -= fc[0];
    resj.xmom -= fc[1];
    resj.ymom -= fc[2];
    resj.ener -= fc[3];
  }
};

struct AUSMUP2Numeric {
  template <preprocess::flowType Type>
  inline static void ComputeResidual(const preprocess::parameter &param,
                                     double nx, double ny, double ds,
                                     PRIM_VAR vari, PRIM_VAR varj,
                                     double gammai, double gammaj,
                                     CONS_VAR &resi, CONS_VAR &resj) {
    double rl = vari.dens;
    double ul = vari.velx;
    double vl = vari.vely;
    double pl = vari.press;
    double hl = gammai / (gammai - 1.0) * pl / rl + 0.5 * (ul * ul + vl * vl);

    double rr = varj.dens;
    double ur = varj.velx;
    double vr = varj.vely;
    double pr = varj.press;
    double hr = gammaj / (gammaj - 1.0) * pr / rr + 0.5 * (ur * ur + vr * vr);

    double cl2 = 2.0 * (gammai - 1.0) / (gammai + 1.0) * hl;
    double cl = std::sqrt(cl2);
    double cr2 = 2.0 * (gammaj - 1.0) / (gammaj + 1.0) * hr;
    double cr = std::sqrt(cr2);

    double veln_left = ul * nx + vl * ny;
    double veln_right = ur * nx + vr * ny;
    double c_bar_left = cl2 / (std::max(cl, veln_left));
    double c_bar_right = cr2 / (std::max(cr, -veln_right));

    double c_interface = std::min(c_bar_left, c_bar_right);
    double m_left = veln_left / c_interface;
    double m_right = veln_right / c_interface;
    double m_bar_2 = (veln_left * veln_left + veln_right * veln_right) /
                     (2.0 * c_interface * c_interface);
    double m_infty;
    if constexpr (Type == preprocess::flowType::External) {
      m_infty = param.MaInfinity;
    } else if constexpr (Type == preprocess::flowType::Internal) {
      double mach = std::sqrt(param.refMach2);
      m_infty = 0.2 * std::min(mach, 1.0);
    }
    double m0_2 = std::min(1.0, std::max(m_bar_2, m_infty * m_infty));
    double m0 = std::sqrt(m0_2);
    double fa = m0 * (2.0 - m0);
    double mp = -0.25 / fa * std::max(1.0 - m_bar_2, 0.0) * (pr - pl) /
                (0.5 * (rr + rl) * c_interface * c_interface);

    double fp_left;
    double fp_right;
    double fm_left;
    double fm_right;
    if (std::abs(m_left) >= 1.0) {
      fm_left = 0.5 * (m_left + std::abs(m_left));
      fp_left = 0.5 * (1.0 + (m_left >= 0.0 ? 1.0 : -1.0));
    } else {
      fm_left = 0.25 * (m_left + 1.0) * (m_left + 1.0) +
                0.125 * (m_left * m_left - 1.0) * (m_left * m_left - 1.0);
      fp_left = 0.25 * (m_left + 1.0) * (m_left + 1.0) * (2.0 - m_left);
    }
    if (std::abs(m_right) >= 1.0) {
      fm_right = 0.5 * (m_right - std::abs(m_right));
      fp_right = 0.5 * (1.0 - (m_right >= 0.0 ? 1.0 : -1.0));
    } else {
      fm_right = -0.25 * (m_right - 1.0) * (m_right - 1.0) -
                 0.125 * (m_right * m_right - 1.0) * (m_right * m_right - 1.0);
      fp_right = 0.25 * (m_right - 1.0) * (m_right - 1.0) * (2.0 + m_right);
    }
    double m_interface = fm_left + fm_right + mp;

    double p_interface =
        0.5 * (pl + pr) + 0.5 * (fp_left - fp_right) * (pl - pr) +
        std::sqrt(0.5 * (ul * ul + vl * vl + ur * ur + vr * vr)) *
            (fp_left + fp_right - 1) * c_interface * 0.5 * (rr + rl);

    double mass_dot;
    if (m_interface > 0.0) {
      mass_dot = m_interface * c_interface * rl;
    } else {
      mass_dot = m_interface * c_interface * rr;
    }

    double mass_flux_left = 0.5 * (mass_dot + std::abs(mass_dot));
    double mass_flux_right = 0.5 * (mass_dot - std::abs(mass_dot));

    double fc[4];

    fc[0] = (mass_flux_left + mass_flux_right) * ds;
    fc[1] =
        (mass_flux_left * ul + mass_flux_right * ur + p_interface * nx) * ds;
    fc[2] =
        (mass_flux_left * vl + mass_flux_right * vr + p_interface * ny) * ds;
    fc[3] = (mass_flux_left * hl + mass_flux_right * hr) * ds;

    resi.dens += fc[0];
    resi.xmom += fc[1];
    resi.ymom += fc[2];
    resi.ener += fc[3];

    resj.dens -= fc[0];
    resj.xmom -= fc[1];
    resj.ymom -= fc[2];
    resj.ener -= fc[3];
  };
};

struct ViscousNumericBound {
  inline static void ComputeResidual(CONS_VAR vari, CONS_VAR varj,
                                     VISC_VAR visci, VISC_VAR visicj,
                                     PRIM_VAR gradx, PRIM_VAR grady,
                                     double gradTx, double gradTy, double sx,
                                     double sy, CONS_VAR &resi) {
    double uav = 0.5 * (vari.xmom / vari.dens + varj.xmom / varj.dens);
    double vav = 0.5 * (vari.ymom / vari.dens + varj.ymom / varj.dens);
    double mav = 0.5 * (visci.mu + visicj.mu);
    double kav = 0.5 * (visci.lambda + visicj.lambda);

    double tauxx = 2.0 / 3.0 * mav * (2.0 * gradx.velx - grady.vely);
    double tauyy = 2.0 / 3.0 * mav * (2.0 * grady.vely - gradx.velx);
    double tauxy = mav * (grady.velx + gradx.vely);
    double phix = uav * tauxx + vav * tauxy + kav * gradTx;
    double phiy = uav * tauxy + vav * tauyy + kav * gradTy;

    double fv[3];
    fv[0] = sx * tauxx + sy * tauxy;
    fv[1] = sx * tauxy + sy * tauyy;
    fv[2] = sx * phix + sy * phiy;

    resi.xmom -= fv[0];
    resi.ymom -= fv[1];
    resi.ener -= fv[2];
  }
};

}  // namespace solver