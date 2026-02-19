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
                op.type = OpType::Diagonal;
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
        return  std::abs(prod + 1) * std::abs(prm.J0);
    }

    case OpSubtype::J1_Plaquette: {
        const auto& p = geom.blackPlaquettes[index];
        int s0 = cfg.spins.spin[p.sites[0]];
        int s1 = cfg.spins.spin[p.sites[1]];
        int s2 = cfg.spins.spin[p.sites[2]];
        int s3 = cfg.spins.spin[p.sites[3]];
        
        // Correct ED match: Diagonal (0-2) and (1-3)
        int prod_sum = s0*s2 + s1*s3;
        
        // Weight = C - H = 2|J| - J * prod_sum
        return (prod_sum == 2.0 ? 4.0 * std::abs(prm.J1) : 0.0);
    }
    case OpSubtype::J2_Dipole: {
        const auto& bond = geom.j2bonds[index];
        int s1 = cfg.spins.spin[bond.s1];
        int s2 = cfg.spins.spin[bond.s2];
        
        double J_local = bond.J_val * prm.J2;
        // Standard SSE weight: C - H_bond = |J| - J * s1 * s2
        // If AntiFerro (J>0): Antip(s1*s2=-1) -> J - (-J) = 2J. Parallel -> 0.
        return std::abs(J_local) - J_local * s1 * s2;
    }

    case OpSubtype::J3_Inter: {
        const auto& bond = geom.j3bonds[index];
        int s1 = cfg.spins.spin[bond[0]];
        int s2 = cfg.spins.spin[bond[1]];
        
        return (s1*s2 + 1 == 0 ? 2.0 * std::abs(prm.J3) : 0.0);
    }

    case OpSubtype::TransverseX_diag: {
        return std::abs(prm.hx);        
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
    // Rename to Cluster Update logically, keeping function name for compatibility if needed.
    // Ideally rename, but let's stick to the existing hook.
    build_vertex_list(cfg, geom, rng);
    check_integrity();

    if (n_legs == 0) {
        // Handle trivial case: No operators. All sites are free spins?
        // Yes. But let's check legs.
        // If n_legs=0, we still might have free sites.
        // Iterate all sites below.
    }

    std::vector<int> cluster_id(n_legs, -1);
    std::vector<bool> cluster_flipped; // Store if cluster ID k is flipped
    int n_clusters = 0;
    
    // BFS Queue
    std::vector<int> queue;
    queue.reserve(n_legs);

    // 1. Identify Clusters
    for(int i=0; i<n_legs; ++i) {
        if (cluster_id[i] == -1) {
            int current_cid = n_clusters++;
            bool do_flip = (rng.uniform() < 0.5);
            cluster_flipped.push_back(do_flip);

            // Start BFS
            cluster_id[i] = current_cid;
            queue.clear();
            queue.push_back(i);

            size_t head = 0;
            while(head < queue.size()){
                int u = queue[head++];
                
                // Neighbor 1: Linked Leg (Time)
                int v1 = linked_leg[u];
                if (v1 != -1 && cluster_id[v1] == -1) {
                    cluster_id[v1] = current_cid;
                    queue.push_back(v1);
                }

                // Neighbor 2: Vertex Partner (Space)
                int v2 = vertex_partner[u];
                if (v2 != -1 && cluster_id[v2] == -1) {
                    cluster_id[v2] = current_cid;
                    queue.push_back(v2);
                }
            }
        }
    }

    // 2. Update Spin Configuration (t=0 boundary)
    // Iterate all sites. If site has legs (first_leg != -1), check its cluster.
    // If site has no legs, it's an independent cluster.
    const int N = (int)cfg.spins.spin.size();
    for(int s=0; s<N; ++s) {
        if (first_leg[s] != -1) {
            int leg = first_leg[s];
            if (cluster_flipped[cluster_id[leg]]) {
                cfg.spins.spin[s] *= -1;
            }
        } else {
            // Free site cluster
            if (rng.uniform() < 0.5) {
                cfg.spins.spin[s] *= -1;
            }
        }
    }

    // 3. Update Operator Types (Transverse Field Boundary Checks)
    // Iterate all operators, look for Transverse Field
    // We need to find their legs.
    // We can iterate vertex_op_index.
    int leg_ptr = 0; // track legs linearly? No, build_vertex_list order is preserved.
    // But leg_counter in build_vertex was linear.
    // Better to use linked_leg_to_op_index implicitly?
    // Let's iterate vertices again logic.
    // Re-construct leg mapping is tedious.
    // Easier: Iterate all legs?
    // We know Transverse Field has 2 legs.
    
    // Loop over non-identity ops using vertex_op_index which matches build_vertex_list order
    // We also know n_op_legs.
    // We need to track the base leg index.
    int current_base = 0;
    for (int vp = 0; vp < (int)vertex_op_index.size(); ++vp) {
        int p = vertex_op_index[vp];
        auto& op = cfg.op_string[p];
        
        int n_this_op_legs = 0;
        
        if (op.subtype == OpSubtype::TransverseX_diag || 
            op.subtype == OpSubtype::TransverseX_offdiag) {
            n_this_op_legs = 2;
            
            // Check boundary flip
            int l_in = current_base;
            int l_out = current_base + 1;
            
            bool flip_in = cluster_flipped[cluster_id[l_in]];
            bool flip_out = cluster_flipped[cluster_id[l_out]];
            
            if (flip_in != flip_out) {
                // Toggle Type
                if (op.subtype == OpSubtype::TransverseX_diag) {
                    op.subtype = OpSubtype::TransverseX_offdiag;
                    op.type = OpType::OffDiagonal;
                } else {
                    op.subtype = OpSubtype::TransverseX_diag;
                    op.type = OpType::Diagonal;
                }
            }
        } 
        else if (op.subtype == OpSubtype::J0_Plaquette || op.subtype == OpSubtype::J1_Plaquette) {
            n_this_op_legs = 8;
        }
        else if (op.subtype == OpSubtype::J2_Dipole || op.subtype == OpSubtype::J3_Inter) {
            n_this_op_legs = 4;
        }
        
        current_base += n_this_op_legs;
    }

    // Done. No propagate_spins needed (global update).
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
    vertex_partner.clear();


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
            n_op_legs = 2;
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
            n_op_legs = 8; 

            // J1 Decomposition: None. (Deterministic Diagonal Pairs 0-2, 1-3)
            // Just register sites.
            // No random choice needed.
            decomp_type = 1;
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
            vertex_partner.resize(leg_counter, -1); // Separate spatial links
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
                
                leg_site[in_leg]     = site;
                leg_site[out_leg]    = site;

                linked_leg_to_op_index[in_leg] = p;
                linked_leg_to_op_index[out_leg] = p;

                base_vertex_leg[in_leg] = base;
                base_vertex_leg[out_leg] = base;

                // Time Linking (Standard)
                if (last_leg[site] != -1) {
                    linked_leg[last_leg[site]] = in_leg;
                    linked_leg[in_leg] = last_leg[site];
                } else {
                    first_leg[site] = in_leg;
                }
                last_leg[site] = out_leg;

                // Vertex Linking (Cluster Construction)
                // Transverse Field: Do NOT link In and Out. 
                // This acts as a 'cut' in the cluster.
                // If the cluster connected to In flips, and Out doesn't, we create a domain wall (Kink).

                // Propagate Spin: Only Off-Diagonal (kinks) flip the worldline (removed)
                // if (op.type == OpType::OffDiagonal) {
                //     spin[site] *= -1;
                // }
                
                // Assign Out Leg Spin (Must be the NEW state) (removed)
                // vertex_spin[out_leg] = spin[site]; 
                
                // last_leg[site] = out_leg;
                // vertex_spin[out_leg] = spin[site]; // Spin tracking removed
            }
            else {
                // J0, J1, J2, J3 (Diagonal Interactions)
                // Link legs to force coupled flips.
                // Iterate sites in pairs/groups based on bond structure.

                // Iterate sites participating in the operator
                for(size_t i=0; i<sites.size(); ++i) {
                     int site = sites[i];
                     int l_in = base + 2*i;
                     int l_out = base + 2*i + 1;
                     
                     // vertex_spin[l_in] = spin[s]; // Spin tracking removed
                     // vertex_spin[l_out] = spin[s]; // Spin tracking removed
                     
                     leg_site[l_in] = site;
                     leg_site[l_out] = site;
                     
                     linked_leg_to_op_index[l_in] = p;
                     linked_leg_to_op_index[l_out] = p;
                     
                     base_vertex_leg[l_in] = base;
                     base_vertex_leg[l_out] = base;
                     
                     // Time Linking
                     if (last_leg[site] != -1) {
                         linked_leg[last_leg[site]] = l_in;
                         linked_leg[l_in] = last_leg[site];
                     } else {
                         first_leg[site] = l_in;
                     }
                     last_leg[site] = l_out;
                }

                // Internal Vertex Linking (Horizontal)
                // Connect Site A to Site B to ensure they flip together.
                
                if (op.subtype == OpSubtype::J0_Plaquette) {
                    // J0: 3 Decompositions (0-1/2-3), (0-2/1-3), (0-3/1-2)
                    int type = vertex_decomp[p]; // Randomized in build_vertex_list
                    
                    int pair1_A = -1, pair1_B = -1;
                    int pair2_A = -1, pair2_B = -1;
                    
                    if (type == 0) {
                         // Pair (0,1) and (2,3)
                         pair1_A = 0; pair1_B = 1;
                         pair2_A = 2; pair2_B = 3;
                    } else if (type == 1) {
                         // Pair (0,2) and (1,3)
                         pair1_A = 0; pair1_B = 2;
                         pair2_A = 1; pair2_B = 3;
                    } else {
                         // Pair (0,3) and (1,2)
                         pair1_A = 0; pair1_B = 3;
                         pair2_A = 1; pair2_B = 2;
                    }

                    // Link First Pair
                    {
                        int l_in_A = base + 2*pair1_A;
                        int l_in_B = base + 2*pair1_B;
                        int l_out_A = base + 2*pair1_A + 1;
                        int l_out_B = base + 2*pair1_B + 1;
                        
                        vertex_partner[l_in_A] = l_in_B; vertex_partner[l_in_B] = l_in_A;
                        vertex_partner[l_out_A] = l_out_B; vertex_partner[l_out_B] = l_out_A;
                    }
                    
                    // Link Second Pair
                    {
                        int l_in_A = base + 2*pair2_A;
                        int l_in_B = base + 2*pair2_B;
                        int l_out_A = base + 2*pair2_A + 1;
                        int l_out_B = base + 2*pair2_B + 1;
                        
                        vertex_partner[l_in_A] = l_in_B; vertex_partner[l_in_B] = l_in_A;
                        vertex_partner[l_out_A] = l_out_B; vertex_partner[l_out_B] = l_out_A;
                    }
                }
                else if (op.subtype == OpSubtype::J1_Plaquette) {
                    // J1: Sum of two bonds (0-2) and (1-3) [Diagonal] to match ED.
                    {
                        int k1=0, k2=2;
                        int l_in_A = base + 2*k1; int l_in_B = base + 2*k2;
                        int l_out_A = base + 2*k1 + 1; int l_out_B = base + 2*k2 + 1;
                        
                        vertex_partner[l_in_A] = l_in_B; vertex_partner[l_in_B] = l_in_A;
                        vertex_partner[l_out_A] = l_out_B; vertex_partner[l_out_B] = l_out_A;
                    }
                    // Link (1,3)
                    {
                        int k1=1, k2=3;
                        int l_in_A = base + 2*k1; int l_in_B = base + 2*k2;
                        int l_out_A = base + 2*k1 + 1; int l_out_B = base + 2*k2 + 1;
                        
                        vertex_partner[l_in_A] = l_in_B; vertex_partner[l_in_B] = l_in_A;
                        vertex_partner[l_out_A] = l_out_B; vertex_partner[l_out_B] = l_out_A;
                    }
                }
                else {
                    // J2/J3 (Dipole). 2 sites. Simple Pair.
                    int l_in_A = base;
                    int l_in_B = base + 2;
                    int l_out_A = base + 1;
                    int l_out_B = base + 3;
                    
                    vertex_partner[l_in_A] = l_in_B; vertex_partner[l_in_B] = l_in_A;
                    vertex_partner[l_out_A] = l_out_B; vertex_partner[l_out_B] = l_out_A;
                }
            }
        }
    }

    // Close Periodic Boundary Conditions
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

    // Transverse Field (hx): 50/50 Probabilistic Transmit/Bounce
    if (op.subtype == OpSubtype::TransverseX_diag || 
        op.subtype == OpSubtype::TransverseX_offdiag) {
        return (rng.uniform() < 0.5) ? (entry_leg ^ 1) : entry_leg;
    }

    // Identify active legs
    int base = base_vertex_leg[entry_leg];
    int rel = entry_leg - base;
    
    int partner_site_idx = -1;
    
    // Logic from previous step to find partner
    if (op.subtype == OpSubtype::J0_Plaquette) {
        int type = vertex_decomp[vertex];
        int site_idx = rel / 2;
        if (type == 0) { 
             if (site_idx == 0) partner_site_idx = 1;
             else if (site_idx == 1) partner_site_idx = 0;
             else if (site_idx == 2) partner_site_idx = 3;
             else if (site_idx == 3) partner_site_idx = 2;
        } else if (type == 1) { 
             if (site_idx == 0) partner_site_idx = 2;
             else if (site_idx == 2) partner_site_idx = 0;
             else if (site_idx == 1) partner_site_idx = 3;
             else if (site_idx == 3) partner_site_idx = 1;
        } else { 
             if (site_idx == 0) partner_site_idx = 3;
             else if (site_idx == 3) partner_site_idx = 0;
             else if (site_idx == 1) partner_site_idx = 2;
             else if (site_idx == 2) partner_site_idx = 1;
        }
    }
    else if (op.subtype == OpSubtype::J1_Plaquette) {
        // No decomposition. Deterministic Diagonal Pairing.
        // 0-2 and 1-3.
        int site_idx = rel / 2;
        int partner_site = -1;
        
        if (site_idx == 0) partner_site = 2;
        else if (site_idx == 2) partner_site = 0;
        else if (site_idx == 1) partner_site = 3;
        else if (site_idx == 3) partner_site = 1;
        
        if (partner_site != -1) {
             // In->Out check?
             // target_rel = partner_site * 2 + (rel % 2 == 0 ? 1 : 0); // "Switch" logic
             // Wait, previous code had manual target_rel calc.
             // Let's use standard partners.
             return base + partner_site * 2 + (rel % 2 == 0 ? 1 : 0);
        }
    }
    else if (op.subtype == OpSubtype::J2_Dipole || op.subtype == OpSubtype::J3_Inter) {
        int site_idx = rel / 2; 
        partner_site_idx = (site_idx == 0) ? 1 : 0;
    }

    // If no partner (inactive leg), forced Vertical
    if (partner_site_idx == -1) return entry_leg ^ 1;



    // Assign canonical indices
    // active_in/out vs partner_in/out
    // Let's use 4 legs:
    // 0: entry_leg's site IN
    // 1: entry_leg's site OUT
    // 2: partner's site IN
    // 3: partner's site OUT
    
    // But we need to map actual legs.
    int site_base = (rel / 2) * 2;
    int p_base = partner_site_idx * 2;
    
    int l0 = base + site_base;     // Site A In
    //int l1 = base + site_base + 1; // Site A Out (Unused)
    int l2 = base + p_base;        // Site B In
    int l3 = base + p_base + 1;    // Site B Out
    
    // Current spins on these legs
    // CAUTION: vertex_spin[entry_leg] is ALREADY flipped by caller!
    // We want to evaluate weights of potential OUTCOMES.
    
    // Candidate 1: Vertical (Exit = entry_leg ^ 1)
    // We flip vertex_spin[entry_leg ^ 1].
    // Result: A_in flipped, A_out flipped. B unchanged.
    // State: (-Original_A, Original_B).
    // Let's evaluate this using current vertex_spins.
    // vertex_spin[l0] (if l0 is entry) is -sA. vertex_spin[l1] is sA.
    // If we flip l1, we get -sA.
    // So both A legs are -sA. B legs are sB.
    // Valid "Diagonal" state with spins (-S_A, S_B).
    
    // Reconstruct for Horizontal choice
    // Horizontal flips B_in (l2).
    // A_in is -sA. B_in becomes -sB.
    // A_out is sA. B_out is sB.
    // State: In(-sA, -sB) -> Out(sA, sB).
    // This is valid Off-Diagonal.
    
    // Evaluate weights
    // S1, S2 for Vertical: -sA, sB.
    int sA_v = vertex_spin[l0]; // Already flipped entry
    int sB_v = vertex_spin[l2]; // Unchanged
    // Actually, if entry was l0, sA_v is -Original_A.
    // W_vert = Weight_Diag(sA_v, sB_v).
    
    // S1, S2 for Horizontal:
    // Transition |-sA, -sB> to |sA, sB>.
    // Weight is W_off.
    
    double w_vert = 0.0;
    double w_horiz = 0.0;
    
    // Calculate J_val and Constant Shift C
    double J_val = 0.0;
    double C = 0.0;
    
    if (op.subtype == OpSubtype::J0_Plaquette) {
        // J0 Spectators determine sign.
        int s_idx1 = site_base/2;
        int s_idx2 = p_base/2;
        
        int spec1 = -1, spec2 = -1;
        for(int k=0; k<4; ++k) {
             if (k!=s_idx1 && k!=s_idx2) {
                 if (spec1==-1) spec1 = k;
                 else spec2 = k;
             }
        }
        
        int l_spec1 = base + spec1*2;
        int l_spec2 = base + spec2*2;
        
        int s_spec1 = vertex_spin[l_spec1];
        int s_spec2 = vertex_spin[l_spec2];
        
        double spectator_sign = (double)(s_spec1 * s_spec2);
        J_val = prm.J0 * spectator_sign; 
        C = 2.0 * std::abs(prm.J0);
    }
    else if (op.subtype == OpSubtype::J1_Plaquette) {
        J_val = prm.J1;
        C = 2.0 * std::abs(prm.J1);
    }
    else if (op.subtype == OpSubtype::J2_Dipole) {
        J_val = geom.j2bonds[op.index].J_val * prm.J2; 
        C = 2.0 * std::abs(prm.J2);
    }
    else if (op.subtype == OpSubtype::J3_Inter) {
        J_val = prm.J3;
        C = 2.0 * std::abs(prm.J3);
    }
    
    // Vertical (Diagonal) Weight
    // W' = C - J_val * prod_v.
    double prod_v = sA_v * sB_v;
    w_vert = 0.0;

    if (C > 0.0) {
        w_vert = C - J_val * prod_v;
    } else {
        // Fallback (should not be reached if shifted)
        double J_abs_local = std::abs(J_val);
        w_vert = J_abs_local - J_val * prod_v; 
    }
    
    if (w_vert < 1e-9) w_vert = 0.0;
    
    // Horizontal (Off-Diagonal) Weight
    // J are Strictly Diagonal (Ising-like). 
    // w_horiz = 0 implies we NEVER switch legs.
    if (op.subtype == OpSubtype::J2_Dipole || op.subtype == OpSubtype::J3_Inter || 
        op.subtype == OpSubtype::J1_Plaquette || op.subtype == OpSubtype::J0_Plaquette) {
        w_horiz = 0.0;
    } else {
        // Default behavior for other operators (if any)
        w_horiz = std::abs(J_val); 
    }

    // Determine candidate horizontal leg (even if weight is 0, define it)
    int cand_horiz = -1;
    // Simple partner logic: 0<->2 (In-In), 1<->3 (Out-Out)
    // base legs: 0, 1 (Site A); 2, 3 (Site B).
    // rel 0 -> 2. rel 1 -> 3. rel 2 -> 0. rel 3 -> 1.
    // This assumes 4-leg vertex with 2 sites.
    // For J1/J0 plaquette, we have 4 sites but decomposing into 2-site bonds.
    // Site mapping was done at start of function (l0, l1, l2, l3).
    // l0=entry, l1=entry^1.
    // l2=partner_In, l3=partner_Out.
    // If entry is Even (In), partner is l2 (In).
    // If entry is Odd (Out), partner is l3 (Out).
    if (rel % 2 == 0) cand_horiz = l2; 
    else cand_horiz = l3;

    // Compute probabilities
    double total = w_vert + w_horiz;
    if (total <= 1e-14) return entry_leg ^ 1; // Fallback
    
    if (rng.uniform() < (w_horiz / total)) {
         return cand_horiz;
    } else {
         return entry_leg ^ 1;
    }
}
