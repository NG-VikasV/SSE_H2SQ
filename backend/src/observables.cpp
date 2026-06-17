#include "../include/variables.hpp"
#include "../include/observables.hpp"

#include <iostream>
#include <ostream>

#include <cmath>


// ======================================================
// Accumulate one SSE configuration
// ======================================================
void Observables::accumulate(const SSEConfig& cfg,
                             const Parameters& prm,
                             const Geometry& geom)
{


    // --------------------------------------------------
    // Energy (SSE exact)
    // E = -<n>/beta + Shift
    // Shift = Sum of constants added to diagonal weights
    // --------------------------------------------------
    double shift = 0.0;

    // J0 shift
    shift += double(geom.blackPlaquettes.size()) * std::abs(prm.J0);

    // J1 shift
    shift += double(geom.blackPlaquettes.size()) * 2.0 * std::abs(prm.J1);

    // J2 shift
    double C_J2 = 0.0;
    for(const auto& b : geom.j2bonds) {
        C_J2 += std::abs(b.J_val * prm.J2);
    }
    C_J2 += (geom.j2_constant_term * prm.J2);
    shift += C_J2; 

    // J3 shift
    shift += double(geom.j3bonds.size()) * std::abs(prm.J3);


    // hx shift
    shift += double(prm.Ns) * std::abs(prm.hx); 

    double e = ((-double(cfg.n_ops) / prm.beta) + shift) / double(prm.Ns);

    // SSE-correct estimators for E^2, E^4:
    // In SSE, <H'^k> = <n(n-1)...(n-k+1)> / beta^k  (falling factorial / beta^k)
    // Expanding H = C - H'  via binomial theorem:
    //
    //   <H^2> = C^2 - 2C * n/beta            + n(n-1)/beta^2
    //   <H^4> = C^4 - 4C^3 * n/beta
    //                + 6C^2 * n(n-1)/beta^2
    //                - 4C   * n(n-1)(n-2)/beta^3
    //                +        n(n-1)(n-2)(n-3)/beta^4
    double n     = double(cfg.n_ops);
    double beta  = prm.beta;
    double Ns    = double(prm.Ns);
    double C     = shift;                 // total constant (not per-site)
    double ib    = 1.0 / beta;           // 1/beta
    double ib2   = ib * ib;
    double ib3   = ib2 * ib;
    double ib4   = ib2 * ib2;

    // Falling factorials
    double ff1 = n;
    double ff2 = n * (n - 1.0);
    double ff3 = ff2 * (n - 2.0);
    double ff4 = ff3 * (n - 3.0);

    // Correct estimator for total <H^2>
    double e2_total = C*C
                    - 2.0*C * ff1 * ib
                    + ff2 * ib2;

    // Correct estimator for total <H^4>
    double e4_total = C*C*C*C
                    - 4.0*C*C*C * ff1 * ib
                    + 6.0*C*C   * ff2 * ib2
                    - 4.0*C     * ff3 * ib3
                    +             ff4 * ib4;

    double e2 = e2_total / (Ns * Ns);
    double e4 = e4_total / (Ns * Ns * Ns * Ns);

    E  += e;
    E2 += e2;
    E4 += e4;

    // --------------------------------------------------
    // Transverse energy contribution (optional bookkeeping)
    //   E_trans = -hx * <σ^x>
    //   where <σ^x> = N_hx / (β N)
    // --------------------------------------------------
    int n_hx = 0;
    for (const auto& op : cfg.op_string) {
        if (op.type == OpType::OffDiagonal &&
            op.subtype == OpSubtype::TransverseX_offdiag) {
            n_hx++;
        }
    }

    double mx_inst = 0.0;
    if (std::abs(prm.hx) > 1e-15) {
        mx_inst = double(n_hx) / (prm.beta * prm.Ns * prm.hx);
    }
    double e_trans_inst = -prm.hx * mx_inst;

    // Transverse energy column
    E_trans += e_trans_inst;

    // --------------------------------------------------
    // Transverse magnetization
    // --------------------------------------------------
    mx  += mx_inst;
    mx2 += mx_inst * mx_inst;
    mx4 += mx_inst * mx_inst * mx_inst * mx_inst;

    // --------------------------------------------------
    // Longitudinal magnetization (Staggered)
    // --------------------------------------------------
    double mz_inst = 0.0;
    for (int i = 0; i < prm.Ns; ++i) {
        int ix = i % prm.Lx;
        int iy = (i / prm.Lx) % prm.Ly;
        int iz = i / (prm.Lx * prm.Ly);
        int phase = ((ix + iy + iz) % 2 == 0) ? 1 : -1;
        mz_inst += 0.5 * cfg.spins.spin[i] * phase;
    }

    mz_inst /= double(prm.Ns);

    mz  += mz_inst;
    mz2 += mz_inst * mz_inst;
    mz4 += mz_inst * mz_inst * mz_inst * mz_inst;

    // --------------------------------------------------
    // Plaquette / order parameter
    // --------------------------------------------------
    double op_inst = 0.0;
    for (const auto& p : geom.blackPlaquettes) {
        int prod = 1;
        for (int s : p.sites)
            prod *= cfg.spins[s];
        op_inst += prod;
    }

    op_inst /= double(geom.blackPlaquettes.size());

    OP  += op_inst;
    OP2 += op_inst * op_inst;
    OP4 += op_inst * op_inst * op_inst * op_inst;

    // --------------------------------------------------
    // Ice-rule indicator (example: |prod| == 1 → ice)
    // You can refine this later without touching main()
    // --------------------------------------------------
    double ice_inst = 0.0;
    for (const auto& p : geom.blackPlaquettes) {
        int sum = 0;
        for (int s : p.sites)
            sum += cfg.spins[s];
        if (std::abs(sum) == 0)
            ice_inst += 1.0;
    }
    ice_inst /= double(geom.blackPlaquettes.size());

    ice  += ice_inst;
    ice2 += ice_inst * ice_inst;

    // --------------------------------------------------
    // Stiffness placeholders (set via winding later)
    // --------------------------------------------------
    // NOTE: Proper SSE stiffness comes from winding numbers.
    // Here we only accumulate placeholders so the pipeline
    // remains consistent.
    rho_x += 0.0;
    rho_y += 0.0;

    n_samples++;
}

// ======================================================
// Normalize (NO derived physics here)
// ======================================================
void Observables::normalize()
{
    if (n_samples == 0) return;

    double inv = 1.0 / double(n_samples);

    E  *= inv;  E2 *= inv;  E4 *= inv;
    E_trans *= inv;

    mz *= inv;  mz2 *= inv;  mz4 *= inv;
    mx *= inv;  mx2 *= inv;  mx4 *= inv;

    OP *= inv;  OP2 *= inv;  OP4 *= inv;

    ice *= inv; ice2 *= inv;

    rho_x *= inv;
    rho_y *= inv;
}
