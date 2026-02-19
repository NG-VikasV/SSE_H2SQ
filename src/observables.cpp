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

    E  += e;
    E2 += e * e;
    E4 += e * e * e * e;

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
    // Longitudinal magnetization
    // --------------------------------------------------
    double mz_inst = 0.0;
    for (int s : cfg.spins.spin)
        mz_inst += s;

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
