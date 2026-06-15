// ======================================================
// Standard headers
// ======================================================
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>
#include <fstream>

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
    // Start small so it scales automatically using Sandvik style L = n + n/3
    cfg.max_ops = 20;

    cfg.n_ops   = 0;
    cfg.op_string.resize(cfg.max_ops);


    // Random start
    cfg.spins.spin.resize(prm.Ns);
    for (int i = 0; i < prm.Ns; ++i) {
        cfg.spins.spin[i] = (rng.uniform() < 0.5) ? -1 : 1;
    }

    // cfg.spins.spin.resize(prm.Ns, -1);   // cold start


    std::cout << "[SSE-CHECK] OK | Initialization complete" << std::endl;
    auto last_update = std::chrono::steady_clock::now();
    auto phase_start = std::chrono::steady_clock::now();
    
    // ────────────────────────────────────────────────────────────────
    // ETA estimator: mirrors tqdm's algorithm exactly.
    //   • Track `avg_spt` = EMA of seconds-per-step  (alpha = 0.3, tqdm default)
    //   • ETA = avg_spt * remaining_steps
    //   • Near completion (progress > 80%) we blend in the global average
    //     so the countdown converges smoothly to 0 (prevents late-run drift).
    //   • Reset state on phase transition (therm → measure).
    // ────────────────────────────────────────────────────────────────
    const double ALPHA = 0.3;          // tqdm default smoothing factor
    double avg_spt  = 0.0;             // seconds per step (EMA)
    int    last_step = 0;              // step index at last emit
    double t_phase   = 0.0;           // total elapsed since phase start (sec)
    std::string cur_phase;

    auto emit_progress = [&](const std::string& phase, int current, int total) {
        auto now = std::chrono::steady_clock::now();
        double ms_since_last = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_update).count();

        // ── Phase transition: reset all estimator state ──
        if (phase != cur_phase) {
            cur_phase  = phase;
            avg_spt    = 0.0;
            last_step  = 0;
            t_phase    = 0.0;
            last_update = now;
            phase_start = now;
            ms_since_last = 0.0;
        }

        if (ms_since_last > 400.0 || current == total) {
            double dt_sec   = ms_since_last / 1000.0;
            int    dn       = current - last_step;     // steps since last emit

            if (dn > 0 && dt_sec > 0.0) {
                double inst_spt = dt_sec / double(dn);  // seconds per step (instant)

                // tqdm EMA: warm-start on first sample, blend thereafter
                avg_spt = (avg_spt <= 0.0)
                    ? inst_spt
                    : ALPHA * inst_spt + (1.0 - ALPHA) * avg_spt;

                last_step = current;
            }

            // Global average (elapsed / steps done) — most accurate near end
            t_phase = std::chrono::duration_cast<std::chrono::milliseconds>(now - phase_start).count() / 1000.0;
            double global_spt = (current > 0) ? t_phase / double(current) : 0.0;

            // Blend: near completion, weight global avg more (it converges to 0)
            double progress = double(current) / double(total);
            double ema_weight = 1.0 - progress;  // 1.0 at start → 0.0 at end
            double blended_spt = ema_weight * avg_spt + (1.0 - ema_weight) * global_spt;

            double eta = (blended_spt > 0.0) ? blended_spt * double(total - current) : 0.0;

            std::cout << "[PROGRESS] " << phase << " " << current << " " << total
                      << " " << std::fixed << std::setprecision(1) << eta << std::endl;
            last_update = now;
        }
    };

    for (int step = 0; step < prm.n_therm; ++step) {
        updates.sweep(cfg, prm, geom, rng);
        updates.adjust_operator_string(cfg);
        emit_progress("therm", step + 1, prm.n_therm);
    }
    std::cout << "[SSE-CHECK] OK | Equilibration finished" << std::endl;

    // --------------------------------------------------
    // Measurements
    // --------------------------------------------------
    Observables obs;
    obs.reset();
    phase_start = std::chrono::steady_clock::now();

    std::ofstream spin_file("spin_configs.json");
    spin_file << "[\n";
    int dump_interval = std::max(1, prm.n_measure / 3); // 4 snapshots: 0, 33%, 66%, 100%

    bool first_dump = true;
    for (int step = 0; step < prm.n_measure; ++step) {
        updates.sweep(cfg, prm, geom, rng);
        obs.accumulate(cfg, prm, geom);
        emit_progress("measure", step + 1, prm.n_measure);

        if (step == 0 || step % dump_interval == 0 || step == prm.n_measure - 1) {
            if (!first_dump) spin_file << ",\n";
            first_dump = false;
            spin_file << "  {\"step\": " << step << ", \"spins\": [";
            for (size_t i = 0; i < cfg.spins.spin.size(); ++i) {
                spin_file << cfg.spins.spin[i] << (i + 1 == cfg.spins.spin.size() ? "" : ", ");
            }
            spin_file << "]}";
        }
    }
    spin_file << "\n]\n";
    spin_file.close();

    std::cout << "[SSE-CHECK] OK | Measurement finished" << std::endl;

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
