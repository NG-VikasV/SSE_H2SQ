#ifndef SLMC_UPDATES_HPP
#define SLMC_UPDATES_HPP

#include <cmath>
#include "variables.hpp"

// =======================================
// SLMC effective model for SSE
// =======================================

struct SLMC_EffectiveModel {

    // Learned effective transverse field
    double hx_eff = 0.0;

    // Accumulators
    long long n_samples = 0;
    long long n_hx_ops  = 0;

    // -----------------------------------
    // Reset accumulators
    // -----------------------------------
    void reset() {
        hx_eff = 0.0;
        n_samples = 0;
        n_hx_ops  = 0;
    }

    // -----------------------------------
    // Accumulate SSE data
    // -----------------------------------
    void accumulate(const SSEConfig& cfg) {

        for (const auto& op : cfg.op_string) {
            if (op.type == OpType::OffDiagonal &&
                op.subtype == OpSubtype::TransverseX) {
                n_hx_ops++;
            }
        }

        n_samples++;
    }

    // -----------------------------------
    // Fit effective transverse field
    // -----------------------------------
    void fit(const Parameters& prm) {
        if (n_samples == 0) return;

        double avg_nhx = double(n_hx_ops) / double(n_samples);

        // ⟨σ^x⟩ = ⟨N_hx⟩ / (β h_x)
        double mx = avg_nhx / (prm.beta * prm.hx);

        // Effective classical coupling
        // hx_eff reproduces same ⟨σ^x⟩
        hx_eff = prm.hx * mx;
    }
};

#endif // SLMC_UPDATES_HPP
