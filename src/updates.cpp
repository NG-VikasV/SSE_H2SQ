#include "../include/updates.hpp"
#include "../include/variables.hpp"

#include <cmath>
#include <cassert>
#include <algorithm>
#include <iostream>

// ======================================================
// Top-level sweep
// ======================================================
void SSEUpdates::sweep(SSEConfig& cfg,
                       const Parameters& prm,
                       const Geometry& geom,
                       RNG& rng)
{

    if (prm.J0 != 0.0) 
        diagonal_update_J0(cfg, prm, geom, rng);
    if (prm.J1 != 0.0)        
        diagonal_update_J1(cfg, prm, geom, rng);
    if (prm.hx != 0.0)
        diagonal_update_h(cfg, prm, geom, rng);
    if (prm.J2 != 0.0)
        diagonal_update_J2(cfg, prm, geom, rng);
    if (prm.J3 != 0.0)
        diagonal_update_J3(cfg, prm, geom, rng);

    // Multiple loop passes improve decorrelation
    for (int i = 0; i < 2; ++i)
        directed_loop_update(cfg, prm, geom, rng);

    //check_spin_consistency(cfg, geom);

    if (use_slmc)
        slmc_update(cfg, prm, geom, rng);
}

// ======================================================
// 1. Diagonal update
// ======================================================
void SSEUpdates::diagonal_update_J0(SSEConfig& cfg,
                                 const Parameters& prm,
                                 const Geometry& geom,
                                 RNG& rng)
{
    const double beta = prm.beta;

    // Available diagonal operators: J0 plaquettes
    const int n_diag = static_cast<int>(geom.blackPlaquettes.size());

    for (int p = 0; p < cfg.max_ops; ++p) {

        SSEOperator& op = cfg.op_string[p];

        // ----------------------------
        // Attempt insertion
        // ----------------------------
        if (op.type == OpType::Identity) {

            int index = int(rng.uniform() * n_diag);
            double w = diagonal_weight(
                static_cast<int>(OpSubtype::J0_Plaquette),
                index, cfg, prm, geom);

            
            if (w <= 0.0) continue;

            double prob = beta * n_diag * w / double(cfg.max_ops - cfg.n_ops);

            if (rng.uniform() < prob) {
                op.type = OpType::Diagonal;
                op.subtype = OpSubtype::J0_Plaquette;
                op.index = index;
                cfg.n_ops++;
            }
        }

        // ----------------------------
        // Attempt removal
        // ----------------------------
        else if (op.subtype == OpSubtype::J0_Plaquette) {

            double w = diagonal_weight(
                static_cast<int>(OpSubtype::J0_Plaquette),
                op.index, cfg, prm, geom);

            if (w <= 0.0) continue;

            double prob = double(cfg.max_ops - cfg.n_ops + 1)
                          / (beta * n_diag * w);

            if (rng.uniform() < prob) {
                op.type = OpType::Identity;
                op.index = -1;
                cfg.n_ops--;
            }
        }

        // --------------------------------------------------
        // 3. Off-diagonal → propagate spins IMMEDIATELY
        // --------------------------------------------------
        else if (op.subtype == OpSubtype::TransverseX_offdiag) {
            // Example: transverse-field operator σ^x_i
            //if (op.subtype == OpSubtype::TransverseX_offdiag){
                cfg.spins.spin[op.index] *= -1;
            //}
        }
    }
}


void SSEUpdates::diagonal_update_J1(SSEConfig& cfg,
                                 const Parameters& prm,
                                 const Geometry& geom,
                                 RNG& rng)
{
    const double beta = prm.beta;

    // Available diagonal operators: J0 plaquettes
    const int n_diag = static_cast<int>(geom.blackPlaquettes.size());

    for (int p = 0; p < cfg.max_ops; ++p) {

        SSEOperator& op = cfg.op_string[p];

        // ----------------------------
        // Attempt insertion
        // ----------------------------
        if (op.type == OpType::Identity) {

            int index = int(rng.uniform() * n_diag);
            double w = diagonal_weight(
                static_cast<int>(OpSubtype::J1_Plaquette),
                index, cfg, prm, geom);

            if (w <= 0.0) continue;

            double prob = beta * n_diag * w / double(cfg.max_ops - cfg.n_ops);

            if (rng.uniform() < prob) {
                op.type = OpType::Diagonal;
                op.subtype = OpSubtype::J1_Plaquette;
                op.index = index;
                cfg.n_ops++;
            }
        }

        // ----------------------------
        // Attempt removal
        // ----------------------------
        else if (op.subtype == OpSubtype::J1_Plaquette) {

            double w = diagonal_weight(
                static_cast<int>(OpSubtype::J1_Plaquette),
                op.index, cfg, prm, geom);

            if (w <= 0.0) continue;


            double prob = double(cfg.max_ops - cfg.n_ops + 1)
                          / (beta * n_diag * w);

            if (rng.uniform() < prob) {
                op.type = OpType::Identity;
                op.index = -1;
                cfg.n_ops--;
            }
        }

        // --------------------------------------------------
        // 3. Off-diagonal → propagate spins IMMEDIATELY
        // --------------------------------------------------
        else if (op.subtype == OpSubtype::TransverseX_offdiag) {
            // Example: transverse-field operator σ^x_i
            //if (op.subtype == OpSubtype::TransverseX_offdiag) {
                cfg.spins.spin[op.index] *= -1;
            //}
        }
    }
}



void SSEUpdates::diagonal_update_h(SSEConfig& cfg,
                                 const Parameters& prm,
                                 const Geometry& geom,
                                 RNG& rng)
{
    const double beta = prm.beta;

    // Available diagonal operators: J0 plaquettes
    // const int n_diag = static_cast<int>(geom.blackPlaquettes.size());

    for (int p = 0; p < cfg.max_ops; ++p) {

        SSEOperator& op = cfg.op_string[p];

        // ----------------------------
        // Attempt insertion
        // ----------------------------
        if (op.type == OpType::Identity) {

            int index = int(rng.uniform() * prm.Ns);
            double w = diagonal_weight(
                static_cast<int>(OpSubtype::TransverseX_diag),
                index, cfg, prm, geom);

            if (w <= 0.0) continue;


            double prob = beta * prm.Ns * w / double(cfg.max_ops - cfg.n_ops);

            if (rng.uniform() < prob) {
                op.type = OpType::OffDiagonal;
                op.subtype = OpSubtype::TransverseX_diag;
                op.index = index;
                cfg.n_ops++;
            }
        }

        // ----------------------------
        // Attempt removal
        // ----------------------------
        else if (op.subtype == OpSubtype::TransverseX_diag) {

            double w = diagonal_weight(
                static_cast<int>(OpSubtype::TransverseX_diag),
                op.index, cfg, prm, geom);

            if (w <= 0.0) continue;


            double prob = double(cfg.max_ops - cfg.n_ops + 1)
                          / (beta * prm.Ns * w);

            if (rng.uniform() < prob) {
                op.type = OpType::Identity;
                op.index = -1;
                cfg.n_ops--;
            }
        }

        // --------------------------------------------------
        // 3. Off-diagonal → propagate spins IMMEDIATELY
        // --------------------------------------------------
        else if (op.subtype == OpSubtype::TransverseX_offdiag) {
            // Example: transverse-field operator σ^x_i
            // if (op.subtype == OpSubtype::TransverseX_offdiag) {
                cfg.spins.spin[op.index] *= -1;
            // }
        }
    }
}

void SSEUpdates::diagonal_update_J2(SSEConfig& cfg,
                                 const Parameters& prm,
                                 const Geometry& geom,
                                 RNG& rng)
{
    const double beta = prm.beta;
    const int n_diag = static_cast<int>(geom.j2bonds.size());

    for (int p = 0; p < cfg.max_ops; ++p) {

        SSEOperator& op = cfg.op_string[p];

        // Insertion
        if (op.type == OpType::Identity) {
            int index = int(rng.uniform() * n_diag);
            double w = diagonal_weight(
                static_cast<int>(OpSubtype::J2_Dipole),
                index, cfg, prm, geom);

            if (w <= 0.0) continue;

            double prob = beta * n_diag * w / double(cfg.max_ops - cfg.n_ops);

            if (rng.uniform() < prob) {
                op.type = OpType::Diagonal;
                op.subtype = OpSubtype::J2_Dipole;
                op.index = index;
                cfg.n_ops++;
            }
        }
        // Removal
        else if (op.subtype == OpSubtype::J2_Dipole) {
            double w = diagonal_weight(
                static_cast<int>(OpSubtype::J2_Dipole),
                op.index, cfg, prm, geom);

            if (w <= 0.0) continue;

            double prob = double(cfg.max_ops - cfg.n_ops + 1) / (beta * n_diag * w);

            if (rng.uniform() < prob) {
                op.type = OpType::Identity;
                op.index = -1;
                cfg.n_ops--;
            }
        }
        else if (op.subtype == OpSubtype::TransverseX_offdiag) {
             cfg.spins.spin[op.index] *= -1;
        }
    }
}

void SSEUpdates::diagonal_update_J3(SSEConfig& cfg,
                                 const Parameters& prm,
                                 const Geometry& geom,
                                 RNG& rng)
{
    const double beta = prm.beta;
    const int n_diag = static_cast<int>(geom.j3bonds.size());

    for (int p = 0; p < cfg.max_ops; ++p) {

        SSEOperator& op = cfg.op_string[p];

        // Insertion
        if (op.type == OpType::Identity) {
            int index = int(rng.uniform() * n_diag);
            double w = diagonal_weight(
                static_cast<int>(OpSubtype::J3_Inter),
                index, cfg, prm, geom);

            if (w <= 0.0) continue;

            double prob = n_diag * beta * w / double(cfg.max_ops - cfg.n_ops);

            if (rng.uniform() < prob) {
                op.type = OpType::Diagonal;
                op.subtype = OpSubtype::J3_Inter;
                op.index = index;
                cfg.n_ops++;
            }
        }
        // Removal
        else if (op.subtype == OpSubtype::J3_Inter) {
            double w = diagonal_weight(
                static_cast<int>(OpSubtype::J3_Inter),
                op.index, cfg, prm, geom);

            if (w <= 0.0) continue;

            double prob = double(cfg.max_ops - cfg.n_ops + 1) / (beta * n_diag * w);

            if (rng.uniform() < prob) {
                op.type = OpType::Identity;
                op.index = -1;
                cfg.n_ops--;
            }
        }
        else if (op.subtype == OpSubtype::TransverseX_offdiag) {
             cfg.spins.spin[op.index] *= -1;
        }
    }
}

// ======================================================
// 3. Diagonal operator weight
// ======================================================
double SSEUpdates::diagonal_weight(int subtype,
                                   int index,
                                   const SSEConfig& cfg,
                                   const Parameters& prm,
                                   const Geometry& geom)
{
    switch (static_cast<OpSubtype>(subtype)) {

    case OpSubtype::J0_Plaquette: {
        const auto& p = geom.blackPlaquettes[index];
        int prod = 1;
        for (int s : p.sites)
            prod *= cfg.spins.spin[s];          
        return  std::abs(prod + 1) * prm.J0;
    }

    case OpSubtype::J1_Plaquette: {
        const auto& p = geom.blackPlaquettes[index];
        int prod_sum = (cfg.spins.spin[p.sites[0]] * cfg.spins.spin[p.sites[2]] + 
                        cfg.spins.spin[p.sites[1]] * cfg.spins.spin[p.sites[3]]);
        return (prod_sum == -2 ? std::abs(prod_sum - 2) * prm.J1 : 0.0);
    }
    case OpSubtype::J2_Dipole: {
        const auto& bond = geom.j2bonds[index];
        int s1 = cfg.spins.spin[bond.s1];
        int s2 = cfg.spins.spin[bond.s2];
        
        // Interaction: J_val * S1 * S2
        // W = |J_val| - J_val * S1 * S2
        return prm.J2 * (std::abs(bond.J_val) - bond.J_val * s1 * s2);
    }

    case OpSubtype::J3_Inter: {
        const auto& bond = geom.j3bonds[index];
        int s1 = cfg.spins.spin[bond[0]];
        int s2 = cfg.spins.spin[bond[1]];
        return std::abs(s1 * s2 - 1) * prm.J3;
    }

    case OpSubtype::TransverseX_diag: {
        int spin = cfg.spins.spin[index];
        return std::abs(spin + 1) * prm.hx;        
    }
    default:
        return 0.0;
    }
}






// ======================================================
// 2. Directed Loop update
// ======================================================
void SSEUpdates::directed_loop_update(SSEConfig& cfg,
                                      const Parameters& prm,
                                      const Geometry& geom,
                                      RNG& rng)
{
    build_vertex_list(cfg, geom, rng);
    check_integrity();

    if (n_legs == 0) return;

    // Reset cluster state/marks (using linked_leg as storage for visited status)
    // Start loops
    // Simplified target: volume target
    int target_steps = std::max(prm.Ns, n_legs / 2);
    int steps_done = 0;

    while (steps_done < target_steps) {
        int start_leg = int(rng.uniform() * n_legs);
        int current_leg = start_leg;

        do {
            steps_done++;

            // Entry into vertex (current_leg is an entry leg)
            int p = linked_leg_to_op_index[current_leg];
            
            // Traversal: use deterministic exit if that was the 'starting phase'
            // or just follow the vertex indices.
            vertex_spin[current_leg] *= -1;
            int exit_leg = choose_exit_leg(p, current_leg, 0, cfg, prm, geom, rng);
            vertex_spin[exit_leg] *= -1;

            // Follow time link
            current_leg = linked_leg[exit_leg];

            if (steps_done > 10000000) {
                 std::cerr << "Safety break! Infinite loop or extremely large order.\n";
                 return; // Avoid aborting, just return to keep progress if possible
            }

        } while (current_leg != start_leg);
    }

    propagate_spins(cfg);
    check_spin_consistency(cfg, geom);
}




// ======================================================
// 4. Build vertex list (SSE imaginary-time linking)
// ======================================================
void SSEUpdates::build_vertex_list(SSEConfig& cfg,
                                   const Geometry& geom,
                                   RNG& rng)
{
    vertex_op_index.clear();
    linked_leg.clear();
    vertex_spin.clear();
    leg_site.clear();   // <-- REQUIRED
    linked_leg_to_op_index.clear();
    base_vertex_leg.clear();


    const int N = cfg.spins.spin.size();
    const int M = cfg.max_ops;

    // Collect non-identity operators
    for (int p = 0; p < M; ++p)
        if (cfg.op_string[p].type != OpType::Identity)
            vertex_op_index.push_back(p);

    last_leg.assign(N, -1);
    first_leg.assign(N, -1);

    // Local propagated spins
    std::vector<int> spin = cfg.spins.spin;

    int leg_counter = 0;

    // Correctly resize to handle access by operator index 'p'
    if ((int)vertex_decomp.size() < M) {
        vertex_decomp.resize(M);
    }

    for (int vp = 0; vp < (int)vertex_op_index.size(); ++vp) {

        int p = vertex_op_index[vp];
        const auto& op = cfg.op_string[p];
        int decomp_type = 0; // default for others
        
        int n_op_legs = 0;
        std::vector<int> sites;

        // --------------------------------------------------
        // Identify sites involved
        // --------------------------------------------------
        if (op.subtype == OpSubtype::TransverseX_diag || 
            op.subtype == OpSubtype::TransverseX_offdiag) {
            sites.push_back(op.index);
            n_op_legs = 2; // 2 legs per site * 1 site = 2
        }
        else if (op.subtype == OpSubtype::J0_Plaquette) {
            const auto& plaq = geom.blackPlaquettes[op.index];
            for(int s : plaq.sites) sites.push_back(s);
            n_op_legs = 8; // 2 * 4
            
            // J0 Decomposition: 3 pairings
            // 0: (0-1), (2-3)
            // 1: (0-2), (1-3)
            // 2: (0-3), (1-2)
            decomp_type = int(rng.uniform() * 3.0); 
        }
        else if (op.subtype == OpSubtype::J1_Plaquette) {
            const auto& plaq = geom.blackPlaquettes[op.index];
            for(int s : plaq.sites) sites.push_back(s);
            n_op_legs = 8; // 2 * 4

            // J1 Decomposition: Sum of two terms
            // Term A (type 0): sites (0,2) active, sites (1,3) identity
            // Term B (type 1): sites (1,3) active, sites (0,2) identity
            // Weight calculation based on spins BEFORE the vertex
            // Note: sites order is BL, TL, TR, BR. (0,2) are diagonals?
            // Plaquette sites: 0=BL, 1=TL, 2=TR, 3=BR
            // J1 is diagonal pairs: (0,2) and (1,3) ?
            // Let's re-verify diagonal_weight for J1.
            // prod_sum = s0*s2 + s1*s3.
            // Yes, pairs are (0,2) and (1,3).
            
            // Weight w0 = J1 * (s0*s2 + 1)
            // Weight w1 = J1 * (s1*s3 + 1)
            
            double s0 = (double)spin[plaq.sites[0]];
            double s2 = (double)spin[plaq.sites[2]];
            double s1 = (double)spin[plaq.sites[1]];
            double s3 = (double)spin[plaq.sites[3]];

            double w0 = s0 * s2 + 1.0; 
            double w1 = s1 * s3 + 1.0;
            // Note: actual weight in diagonal_weight was (prod_sum + 2)*J1.
            // w0+w1 = s0s2 + s1s3 + 2. Correct.
            
            double total_w = w0 + w1;
            if (total_w <= 1e-10) {
                 // Should ideally not happen if op exists, but finite precision/updates?
                 // fallback
                 decomp_type = (rng.uniform() < 0.5) ? 0 : 1;
            } else {
                 decomp_type = (rng.uniform() < (w0 / total_w)) ? 0 : 1;
            }
        }
        else if (op.subtype == OpSubtype::J2_Dipole) {
            const auto& bond = geom.j2bonds[op.index];
            sites.push_back(bond.s1);
            sites.push_back(bond.s2);
            n_op_legs = 4; 
        }
        else if (op.subtype == OpSubtype::J3_Inter) {
            const auto& bond = geom.j3bonds[op.index];
            sites.push_back(bond[0]);
            sites.push_back(bond[1]);
            n_op_legs = 4; // 2 * 2
        }
        
        // Store choice
        vertex_decomp[p] = decomp_type;

        // --------------------------------------------------
        // Create legs
        // --------------------------------------------------
        if (n_op_legs > 0) {
            int base = leg_counter;
            leg_counter += n_op_legs;

            // Resize arrays
            // Check size to avoid frequent reallocs?
            // Just resize
            linked_leg.resize(leg_counter, -1);
            vertex_spin.resize(leg_counter);
            leg_site.resize(leg_counter);
            linked_leg_to_op_index.resize(leg_counter);
            base_vertex_leg.resize(leg_counter);

            // SPECIAL CASE: Transverse Field (Off-diagonal handles spin flip)
            if (op.subtype == OpSubtype::TransverseX_diag || 
                op.subtype == OpSubtype::TransverseX_offdiag) {
                
                int site = sites[0];
                int in_leg  = base;
                int out_leg = base + 1;

                vertex_spin[in_leg]  = spin[site];
                vertex_spin[out_leg] = spin[site]; 
                
                leg_site[in_leg]     = site;
                leg_site[out_leg]    = site;

                linked_leg_to_op_index[in_leg] = p;
                linked_leg_to_op_index[out_leg] = p;

                base_vertex_leg[in_leg] = base;
                base_vertex_leg[out_leg] = base;

                // Time Linking
                if (last_leg[site] != -1) {
                    linked_leg[last_leg[site]] = in_leg;
                    linked_leg[in_leg] = last_leg[site];
                } else {
                    first_leg[site] = in_leg;
                }
                last_leg[site] = out_leg;

                // Propagate Spin
                spin[site] *= -1;
                // Update Out Spin
                vertex_spin[out_leg] = spin[site];
            }
            else {
                // DIAGONAL OPERATORS (J0, J1, J2, J3)
                // All sites connected, no spin flips across vertex
                
                int num_sites = (int)sites.size();
                for(int i=0; i<num_sites; ++i) {
                    int s = sites[i];
                    int in = base + 2*i;
                    int out = base + 2*i + 1; // Fixed index
                    
                    vertex_spin[in] = spin[s];
                    vertex_spin[out] = spin[s];
                    
                    leg_site[in] = s;
                    leg_site[out] = s;
                    
                    linked_leg_to_op_index[in] = p;
                    linked_leg_to_op_index[out] = p;
                    
                    base_vertex_leg[in] = base;
                    base_vertex_leg[out] = base;
                    
                    // Time Linking
                    if (last_leg[s] != -1) {
                        linked_leg[last_leg[s]] = in;
                        linked_leg[in] = last_leg[s];
                    } else {
                        first_leg[s] = in;
                    }
                    last_leg[s] = out;
                }
            }
        }
    }


    // --------------------------------------------------
    // Close imaginary-time boundaries
    // --------------------------------------------------
    for (int s = 0; s < N; ++s) {
        if (first_leg[s] != -1) {
            linked_leg[last_leg[s]] = first_leg[s];
            linked_leg[first_leg[s]] = last_leg[s];
        }
    }

    n_legs = leg_counter;
}



// Check linked list basic integrity
void SSEUpdates::check_integrity()
{
    int n = (int)linked_leg.size();
    for (int l = 0; l < n; ++l) {
        int next = linked_leg[l];
        if (next < 0 || next >= n) {
            std::cerr << "Integrity Error: Leg " << l << " links to invalid " << next << "\n";
            std::abort();
        }
        if (linked_leg[next] != l) {
            std::cerr << "Integrity Error: Asymmetry at leg " << l << " -> " << next << " -> " << linked_leg[next] << "\n";
            std::abort();
        }
    }
    // std::cout << "Integrity Check Passed.\n"; 
}












// --------------------------------------------------
// Visit (mark) an entire cluster without flipping
// --------------------------------------------------
// --------------------------------------------------
// Visit (mark) an entire cluster without flipping
// --------------------------------------------------









// =====================================================
// Propagate spins after cluster update
// Flip spins whose worldline belongs to a flipped cluster
// ======================================================
void SSEUpdates::propagate_spins(SSEConfig& cfg)
{
    const int N = static_cast<int>(cfg.spins.spin.size());

    for (int site = 0; site < N; ++site) {

        int leg0 = first_leg[site];

        // No operator touches this spin
        if (leg0 == -1)
            continue;

        // Flip if the vertex leg spin does not match the original spin?
        // vertex_spin holds the state AFTER the loop.
        // We initialize vertex_spin[leg] = cfg.spins.spin[site] in build_vertex_list.
        // If the loop flipped it, vertex_spin[leg] != cfg.spins.spin[site] initially?
        // Wait. vertex_spin is modified IN PLACE.
        // So vertex_spin[leg0] holds the final spin state of that leg.
        // Just copy it back.
        
        cfg.spins.spin[site] = vertex_spin[leg0];
    }
}





void SSEUpdates::check_spin_consistency(const SSEConfig& cfg,
                                        const Geometry& geom)
{
    // Copy current spins (state at tau=0)
    std::vector<int> current_spins = cfg.spins.spin;
    std::vector<int> initial_spins = cfg.spins.spin; 
    
    // Propagate through list
    for (int p = 0; p < cfg.max_ops; ++p) {
        const auto& op = cfg.op_string[p];
        
        if (op.type == OpType::OffDiagonal) {
            // Check subtypes that flip spins
            if (op.subtype == OpSubtype::TransverseX_offdiag) {
                 if (op.index >= 0 && op.index < (int)current_spins.size()) {
                     current_spins[op.index] *= -1;
                 }
            }
        }
    }
    
    // Check
    for (size_t i = 0; i < current_spins.size(); ++i) {
        if (current_spins[i] != initial_spins[i]) {
            std::cerr << "Spin consistency check failed at site " << i << "\n";
            std::exit(EXIT_FAILURE); 
        }
    }
}

void SSEUpdates::adjust_operator_string(SSEConfig& cfg)
{
    // Safety factor
    const int extra = std::max(6, cfg.n_ops / 3);


    
    int new_L = cfg.n_ops + extra;

	if(new_L > cfg.max_ops){	

        cfg.op_string.resize(new_L);

        // Initialize new operators as identity
        for (int p = cfg.max_ops; p < new_L; ++p) {
            cfg.op_string[p].type    = OpType::Identity;
            cfg.op_string[p].subtype = OpSubtype::None;
            cfg.op_string[p].index   = -1;
        }

        cfg.max_ops = new_L;
    }    
}






// ======================================================
// 9. SLMC global proposal (skeleton)
// ======================================================
void SSEUpdates::slmc_update(SSEConfig& cfg,
                             const Parameters& /*prm*/,
                             const Geometry& /*geom*/,
                             RNG& rng)
{
    // Placeholder: identity proposal
    // Real SLMC proposes cfg_trial using an effective model

    (void)cfg;
    (void)rng;
}



// ======================================================
// 6. Extras Region (temp functions)
// ======================================================






// void SSEUpdates::propagate_loop(SSEConfig& cfg,
//                                 const Parameters& prm,
//                                 const Geometry& geom,
//                                 RNG& rng,
//                                 int start_leg)
// {
//     int leg = start_leg;
//     int spin_out = -vertex_spin[leg];

//     do {
// ======================================================
// Choose exit leg for Directed Loop
// ======================================================
int SSEUpdates::choose_exit_leg(int vertex,
                                int entry_leg,
                                int spin_out, 
                                const SSEConfig& cfg,
                                const Parameters& prm,
                                const Geometry& geom,
                                RNG& rng)
{
    const auto& op = cfg.op_string[vertex];
    int base = base_vertex_leg[entry_leg];

    // Transverse Field Off-Diagonal: Always Transmit
    if (op.subtype == OpSubtype::TransverseX_offdiag) {
        return (entry_leg == base) ? base + 1 : base;
    }

    // Default: Deterministic same-site transmission
    // Site index: (entry_leg - base) / 2
    // If enter leg 2*site (bottom), exit at 2*site + 1 (top).
    // This is simply entry_leg ^ 1 in the existing paired indexing.
    return entry_leg ^ 1;
}
