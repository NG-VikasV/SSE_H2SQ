#ifndef SSE_OBSERVABLES_HPP
#define SSE_OBSERVABLES_HPP

#include <cmath>
#include <vector>
#include <ostream>

// Forward declarations
struct SSEConfig;
struct Parameters;
struct Geometry;

struct Observables {

    // ===============================
    // Counters
    // ===============================
    long long n_samples = 0;

    // ===============================
    // Energy-related (SSE exact)
    // ===============================
    double E      = 0.0;
    double E2     = 0.0;
    double E4     = 0.0;

    // Transverse-field contribution
    double E_trans = 0.0;

    // ===============================
    // Magnetization
    // ===============================
    double mz     = 0.0;
    double mz2    = 0.0;
    double mz4    = 0.0;  // Add this

    double mx     = 0.0;
    double mx2    = 0.0;
    double mx4    = 0.0;  // Add this


    // ===============================
    // Order parameter (plaquette / ice)
    // ===============================
    double OP     = 0.0;
    double OP2    = 0.0;
    double OP4    = 0.0;

    double ice    = 0.0;
    double ice2   = 0.0;

    // ===============================
    // Stiffness / winding (optional)
    // ===============================
    double rho_x  = 0.0;
    double rho_y  = 0.0;

    // ===============================
    // Control
    // ===============================
    void reset() {
        n_samples = 0;
        E = E2 = E4 = 0.0;
        E_trans = 0.0;
        mz = mz2 = mz4 = 0.0;  // Update this
        mx = mx2 = mx4 = 0.0;  // Update this
        OP = OP2 = OP4 = 0.0;
        ice = ice2 = 0.0;
        rho_x = rho_y = 0.0;
    }

    // ===============================
    // Accumulate one SSE measurement
    // ===============================
    void accumulate(const SSEConfig& cfg,
                    const Parameters& prm,
                    const Geometry& geom);

    // ===============================
    // Normalize averages
    // ===============================
    void normalize();

    // ===============================
    // Binder cumulants
    // ===============================
    double binder_energy() const {
        return 1.0 - E4 / (3.0 * E2 * E2);
    }

    double binder_op() const {
        return 1.0 - OP4 / (3.0 * OP2 * OP2);
    }

    double binder_mz() const {  // Add this
        return 1.0 - mz4 / (3.0 * mz2 * mz2);
    }    
};

#endif
