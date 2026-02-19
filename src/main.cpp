// ======================================================
// Standard headers
// ======================================================
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>

// ======================================================
// Project headers
// ======================================================
#include "../include/variables.hpp"     // Parameters, SSEConfig, Geometry
#include "../include/updates.hpp"       // SSEUpdates
#include "../include/observables.hpp"   // Observables (exact fields only)
#include "../include/writeresults.hpp"  // ResultWriter
#include "../include/read_input.hpp"    // read_parameters

// ======================================================
// Main
// ======================================================
int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: ./sse_sim input.in\n";
        return EXIT_FAILURE;
    }

    using clock = std::chrono::high_resolution_clock;
    auto t_start = clock::now();

    // --------------------------------------------------
    // Read parameters
    // --------------------------------------------------
    Parameters prm;
    read_parameters(argv[1], prm);
    prm.print_parameters();
    
    // --------------------------------------------------
    // Geometry
    // --------------------------------------------------
    Geometry geom;
    geom.build(prm.Lx, prm.Ly, prm.Lz);
    geom.print_black_plaquettes();
    geom.print_j2_bonds();


    // --------------------------------------------------
    // Initialize SSE configuration
    // --------------------------------------------------
    SSEConfig cfg;
    // --------------------------------------------------
    // RNG and update engine
    // --------------------------------------------------
    RNG rng;
    SSEUpdates updates;

    // Initial estimate for max_ops 
    // <n> ~ beta * |E|. For J0=1, E ~ Ns. So <n> ~ beta * Ns.
    // We add a safety buffer.
    cfg.max_ops = static_cast<int>(2.0 * prm.beta * prm.Ns) + 1000;
    
    // Safety check
    if (cfg.max_ops < 1000) cfg.max_ops = 1000;

    cfg.n_ops   = 0;
    cfg.op_string.resize(cfg.max_ops);


    // Random start
    cfg.spins.spin.resize(prm.Ns);
    for (int i = 0; i < prm.Ns; ++i) {
        cfg.spins.spin[i] = (rng.uniform() < 0.5) ? -1 : 1;
    }

    // cfg.spins.spin.resize(prm.Ns, -1);   // cold start


    // --------------------------------------------------
    // Equilibration
    // --------------------------------------------------
    for (int step = 0; step < prm.n_therm; ++step) {
        // std::cout << "Max Ops:  " << cfg.max_ops << " No of ops: "<< cfg.n_ops << std::endl;
        updates.sweep(cfg, prm, geom, rng);
        updates.adjust_operator_string(cfg);
    }

    // --------------------------------------------------
    // Measurements
    // --------------------------------------------------
    Observables obs;
    obs.reset();

    for (int step = 0; step < prm.n_measure; ++step) {
        updates.sweep(cfg, prm, geom, rng);
        obs.accumulate(cfg, prm, geom);
    }

    obs.normalize();   // normalize using n_samples inside obs

    // --------------------------------------------------
    // Output
    // --------------------------------------------------
    ResultWriter writer("SSE_data.txt", /*append=*/true);

    // Metadata (strongly recommended)
    // writer.write_comment("SSE simulation");
    // writer.write_comment(
    //     "Lx=" + std::to_string(prm.Lx) +
    //     " Ly=" + std::to_string(prm.Ly) +
    //     " beta=" + std::to_string(prm.beta)
    // );

    // Headers must match EXACT observable fields
    // Headers must match EXACT observable fields
    std::vector<std::string> headers = {
        "beta", "J0", "J1", "J2", "J3", "hx",
        "E", "E2", "E4",
        "E_trans",
        "mz", "mz2", "mz4",
        "mx", "mx2", "mx4",
        "OP", "OP2", "OP4",
        "ice", "ice2",
        "rho_x", "rho_y"
    };

    std::vector<double> values = {
        prm.beta, prm.J0, prm.J1, prm.J2, prm.J3, prm.hx,
        obs.E, obs.E2, obs.E4,
        obs.E_trans,
        obs.mz, obs.mz2, obs.mz4,
        obs.mx, obs.mx2, obs.mx4,
        obs.OP, obs.OP2, obs.OP4,
        obs.ice, obs.ice2,
        obs.rho_x, obs.rho_y
    };

    writer.write_row(headers, values);

    // --------------------------------------------------
    // Timing
    // --------------------------------------------------
    auto t_end = clock::now();
    double elapsed =
        std::chrono::duration<double>(t_end - t_start).count();

    std::cout << "SSE simulation finished in "
              << std::fixed << std::setprecision(3)
              << elapsed << " seconds\n";

    return 0;
}
