#include "../include/updates.hpp"
#include "../include/variables.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

// ======================================================
// Top-level sweep
// ======================================================
void SSEUpdates::sweep(SSEConfig &cfg, const Parameters &prm,
                       const Geometry &geom, RNG &rng) {
  diagonal_update(cfg, prm, geom, rng);

  for (int i = 0; i < 2; ++i)
      directed_loop_update(cfg, prm, geom, rng);

  // Gauge-sector update restores ergodicity: the BFS cluster collapses all
  // sites into one giant cluster (they're connected through shared plaquette
  // vertices), so the cluster update alone can only perform a global Z2 flip.
  // Row/column flips can change topological sector (e.g. uniform -> stripe)
  // while always preserving every J0 plaquette constraint.
  gauge_update(cfg, prm, geom, rng);

  if (use_slmc)
    slmc_update(cfg, prm, geom, rng);
}

// ======================================================
// 1. Diagonal update
// ======================================================
void SSEUpdates::diagonal_update(SSEConfig &cfg, const Parameters &prm,
                                 const Geometry &geom, RNG &rng) {
  // Rebuild active-subtype cache — member vectors avoid per-sweep allocation
  m_active_subtypes.clear();
  m_n_diags.clear();
  std::fill(std::begin(m_n_diag_by_sub), std::end(m_n_diag_by_sub), 0);

  auto add_type = [&](OpSubtype sub, int nd) {
    m_active_subtypes.push_back(sub);
    m_n_diags.push_back(nd);
    m_n_diag_by_sub[static_cast<int>(sub)] = nd;
  };

  if (prm.J0 != 0.0)
    add_type(OpSubtype::J0_Plaquette,
             static_cast<int>(geom.blackPlaquettes.size()));
  if (prm.J1 != 0.0)
    add_type(OpSubtype::J1_Plaquette,
             static_cast<int>(geom.blackPlaquettes.size()));
  if (prm.J2 != 0.0)
    add_type(OpSubtype::J2_Dipole, static_cast<int>(geom.j2bonds.size()));
  if (prm.J3 != 0.0)
    add_type(OpSubtype::J3_Inter, static_cast<int>(geom.j3bonds.size()));
  if (prm.hx != 0.0)
    add_type(OpSubtype::TransverseX_diag, prm.Ns);

  m_N_active = static_cast<int>(m_active_subtypes.size());

  for (int p = 0; p < cfg.max_ops; ++p) {
    SSEOperator &op = cfg.op_string[p];

    if (op.type == OpType::Identity) {
      if (m_N_active > 0) {
        int type_idx = int(rng.uniform() * m_N_active);
        OpSubtype sub = m_active_subtypes[type_idx];
        int n_diag = m_n_diags[type_idx];
        int index = int(rng.uniform() * n_diag);

        double w =
            diagonal_weight(static_cast<int>(sub), index, cfg, prm, geom);
        if (w > 0.0) {
          double prob = m_N_active * prm.beta * n_diag * w /
                        double(cfg.max_ops - cfg.n_ops);
          if (rng.uniform() < prob) {
            op.type = OpType::Diagonal;
            op.subtype = sub;
            op.index = index;
            ++cfg.n_ops;
          }
        }
      }
    } else if (op.type == OpType::Diagonal) {
      OpSubtype sub = op.subtype;
      int n_diag = m_n_diag_by_sub[static_cast<int>(sub)]; // O(1) array lookup

      double w =
          diagonal_weight(static_cast<int>(sub), op.index, cfg, prm, geom);
      if (w > 0.0) {
        double prob = double(cfg.max_ops - cfg.n_ops + 1) /
                      (m_N_active * prm.beta * n_diag * w);
        if (rng.uniform() < prob) {
          op.type = OpType::Identity;
          op.subtype = OpSubtype::None;
          op.index = -1;
          --cfg.n_ops;
        }
      }
    } else if (op.subtype == OpSubtype::TransverseX_offdiag) {
      cfg.spins.spin[op.index] *= -1;
    }
  }
}

// ======================================================
// 3. Diagonal operator weight
// ======================================================
double SSEUpdates::diagonal_weight(int subtype, int index, const SSEConfig &cfg,
                                   const Parameters &prm,
                                   const Geometry &geom) {
  switch (static_cast<OpSubtype>(subtype)) {

  case OpSubtype::J0_Plaquette: {
    const auto &p = geom.blackPlaquettes[index];
    int prod = 1;
    for (int s : p.sites)
      prod *= cfg.spins.spin[s];
    return std::abs(prod + 1) * std::abs(prm.J0);
  }

  case OpSubtype::J1_Plaquette: {
    const auto &p = geom.blackPlaquettes[index];
    int prod_sum = cfg.spins.spin[p.sites[0]] * cfg.spins.spin[p.sites[2]] +
                   cfg.spins.spin[p.sites[1]] * cfg.spins.spin[p.sites[3]];
    return (prod_sum == 2 ? 4.0 * std::abs(prm.J1) : 0.0);
  }

  case OpSubtype::J2_Dipole: {
    const auto &bond = geom.j2bonds[index];
    double J_local = bond.J_val * prm.J2;
    return std::abs(J_local) -
           J_local * cfg.spins.spin[bond.s1] * cfg.spins.spin[bond.s2];
  }

  case OpSubtype::J3_Inter: {
    const auto &bond = geom.j3bonds[index];
    return (cfg.spins.spin[bond[0]] * cfg.spins.spin[bond[1]] + 1 == 0)
               ? 2.0 * std::abs(prm.J3)
               : 0.0;
  }

  case OpSubtype::TransverseX_diag:
    return std::abs(prm.hx);

  default:
    return 0.0;
  }
}

// ======================================================
// 2. Directed Loop update (cluster variant)
// ======================================================
void SSEUpdates::directed_loop_update(SSEConfig &cfg, const Parameters &prm,
                                      const Geometry &geom, RNG &rng) {
  build_vertex_list(cfg, geom, rng);
#ifndef NDEBUG
  check_integrity();
#endif

  // --- 1. Identify clusters via BFS (member scratch — no per-call allocation)
  // ---
  m_cluster_id.assign(n_legs, -1);
  m_cluster_flipped.clear();
  m_cluster_flipped.reserve(n_legs);
  int n_clusters = 0;

  for (int i = 0; i < n_legs; ++i) {
    if (m_cluster_id[i] != -1)
      continue;

    const int cid = n_clusters++;
    const bool do_flip = (rng.uniform() < 0.5);
    m_cluster_flipped.push_back(do_flip);

    m_cluster_id[i] = cid;
    m_bfs_queue.clear();
    m_bfs_queue.push_back(i);

    for (size_t head = 0; head < m_bfs_queue.size(); ++head) {
      const int u = m_bfs_queue[head];

      const int v1 = linked_leg[u];
      if (v1 != -1 && m_cluster_id[v1] == -1) {
        m_cluster_id[v1] = cid;
        m_bfs_queue.push_back(v1);
      }
      const int v2 = vertex_partner[u];
      if (v2 != -1 && m_cluster_id[v2] == -1) {
        m_cluster_id[v2] = cid;
        m_bfs_queue.push_back(v2);
      }
    }
  }

  // --- 2. Update t=0 boundary spins ---
  const int N = static_cast<int>(cfg.spins.spin.size());
  for (int s = 0; s < N; ++s) {
    if (first_leg[s] != -1) {
      if (m_cluster_flipped[m_cluster_id[first_leg[s]]])
        cfg.spins.spin[s] *= -1;
    } else {
      if (rng.uniform() < 0.5)
        cfg.spins.spin[s] *= -1;
    }
  }

  // --- 3. Toggle TransverseX operators that straddle a cluster boundary ---
  // vertex_base_vec[vp] pre-computed in build_vertex_list: no base-tracking
  // needed here.
  for (int vp = 0; vp < static_cast<int>(vertex_op_index.size()); ++vp) {
    const int p = vertex_op_index[vp];
    auto &op = cfg.op_string[p];

    if (op.subtype != OpSubtype::TransverseX_diag &&
        op.subtype != OpSubtype::TransverseX_offdiag)
      continue;

    const int base = vertex_base_vec[vp];
    const bool flip_in = m_cluster_flipped[m_cluster_id[base]];
    const bool flip_out = m_cluster_flipped[m_cluster_id[base + 1]];

    if (flip_in != flip_out) {
      if (op.subtype == OpSubtype::TransverseX_diag) {
        op.subtype = OpSubtype::TransverseX_offdiag;
        op.type = OpType::OffDiagonal;
      } else {
        op.subtype = OpSubtype::TransverseX_diag;
        op.type = OpType::Diagonal;
      }
    }
  }

#ifndef NDEBUG
  check_spin_consistency(cfg, geom);
#endif
}

// ======================================================
// 4. Build vertex list (SSE imaginary-time linking)
// ======================================================
void SSEUpdates::build_vertex_list(SSEConfig &cfg, const Geometry &geom,
                                   RNG &rng) {
  vertex_op_index.clear();
  linked_leg.clear();
  vertex_spin.clear();
  leg_site.clear();
  linked_leg_to_op_index.clear();
  base_vertex_leg.clear();
  vertex_partner.clear();
  vertex_base_vec.clear();

  const int N = static_cast<int>(cfg.spins.spin.size());
  const int M = cfg.max_ops;

  for (int p = 0; p < M; ++p)
    if (cfg.op_string[p].type != OpType::Identity)
      vertex_op_index.push_back(p);

  const int max_legs = cfg.n_ops * 8;
  linked_leg.assign(max_legs, -1);
  vertex_partner.assign(max_legs, -1);
  vertex_spin.resize(max_legs);
  leg_site.resize(max_legs);
  linked_leg_to_op_index.resize(max_legs);
  base_vertex_leg.resize(max_legs);

  last_leg.assign(N, -1);
  first_leg.assign(N, -1);

  vertex_base_vec.reserve(vertex_op_index.size());

  std::vector<int> spin = cfg.spins.spin;

  int leg_counter = 0;

  if (static_cast<int>(vertex_decomp.size()) < M)
    vertex_decomp.resize(M);

  for (int vp = 0; vp < static_cast<int>(vertex_op_index.size()); ++vp) {
    const int p = vertex_op_index[vp];
    const auto &op = cfg.op_string[p];
    int decomp_type = 0;

    // Stack-allocated buffer — avoids one heap allocation per vertex in the hot
    // loop
    int sites_buf[4];
    int n_sites = 0;
    int n_op_legs = 0;

    if (op.subtype == OpSubtype::TransverseX_diag ||
        op.subtype == OpSubtype::TransverseX_offdiag) {
      sites_buf[0] = op.index;
      n_sites = 1;
      n_op_legs = 2;
    } else if (op.subtype == OpSubtype::J0_Plaquette) {
      const auto &plaq = geom.blackPlaquettes[op.index];
      sites_buf[0] = plaq.sites[0];
      sites_buf[1] = plaq.sites[1];
      sites_buf[2] = plaq.sites[2];
      sites_buf[3] = plaq.sites[3];
      n_sites = 4;
      n_op_legs = 8;
      decomp_type = static_cast<int>(rng.uniform() * 3.0);
    } else if (op.subtype == OpSubtype::J1_Plaquette) {
      const auto &plaq = geom.blackPlaquettes[op.index];
      sites_buf[0] = plaq.sites[0];
      sites_buf[1] = plaq.sites[1];
      sites_buf[2] = plaq.sites[2];
      sites_buf[3] = plaq.sites[3];
      n_sites = 4;
      n_op_legs = 8;
      decomp_type = 1;
    } else if (op.subtype == OpSubtype::J2_Dipole) {
      const auto &bond = geom.j2bonds[op.index];
      sites_buf[0] = bond.s1;
      sites_buf[1] = bond.s2;
      n_sites = 2;
      n_op_legs = 4;
    } else if (op.subtype == OpSubtype::J3_Inter) {
      const auto &bond = geom.j3bonds[op.index];
      sites_buf[0] = bond[0];
      sites_buf[1] = bond[1];
      n_sites = 2;
      n_op_legs = 4;
    }

    vertex_decomp[p] = decomp_type;
    vertex_base_vec.push_back(
        leg_counter); // store base for directed_loop_update

    if (n_op_legs > 0) {
      const int base = leg_counter;
      leg_counter += n_op_legs;

      if (op.subtype == OpSubtype::TransverseX_diag ||
          op.subtype == OpSubtype::TransverseX_offdiag) {

        const int site = sites_buf[0];
        const int in_leg = base;
        const int out_leg = base + 1;

        vertex_spin[in_leg] = spin[site];

        leg_site[in_leg] = site;
        leg_site[out_leg] = site;

        linked_leg_to_op_index[in_leg] = p;
        linked_leg_to_op_index[out_leg] = p;

        base_vertex_leg[in_leg] = base;
        base_vertex_leg[out_leg] = base;

        if (last_leg[site] != -1) {
          linked_leg[last_leg[site]] = in_leg;
          linked_leg[in_leg] = last_leg[site];
        } else {
          first_leg[site] = in_leg;
        }
        last_leg[site] = out_leg;
      } else {
        // J0, J1, J2, J3: time-link all site legs
        for (int i = 0; i < n_sites; ++i) {
          const int site = sites_buf[i];
          const int l_in = base + 2 * i;
          const int l_out = l_in + 1;

          leg_site[l_in] = site;
          leg_site[l_out] = site;

          linked_leg_to_op_index[l_in] = p;
          linked_leg_to_op_index[l_out] = p;

          base_vertex_leg[l_in] = base;
          base_vertex_leg[l_out] = base;

          if (last_leg[site] != -1) {
            linked_leg[last_leg[site]] = l_in;
            linked_leg[l_in] = last_leg[site];
          } else {
            first_leg[site] = l_in;
          }
          last_leg[site] = l_out;
        }

        // Internal vertex partner links (force paired-site cluster flips)
        if (op.subtype == OpSubtype::J0_Plaquette) {
          const int type = vertex_decomp[p];
          int p1A, p1B, p2A, p2B;
          if (type == 0) {
            p1A = 0;
            p1B = 1;
            p2A = 2;
            p2B = 3;
          } else if (type == 1) {
            p1A = 0;
            p1B = 2;
            p2A = 1;
            p2B = 3;
          } else {
            p1A = 0;
            p1B = 3;
            p2A = 1;
            p2B = 2;
          }
          {
            int liA = base + 2 * p1A, liB = base + 2 * p1B;
            vertex_partner[liA] = liB;
            vertex_partner[liB] = liA;
            vertex_partner[liA + 1] = liB + 1;
            vertex_partner[liB + 1] = liA + 1;
          }
          {
            int liA = base + 2 * p2A, liB = base + 2 * p2B;
            vertex_partner[liA] = liB;
            vertex_partner[liB] = liA;
            vertex_partner[liA + 1] = liB + 1;
            vertex_partner[liB + 1] = liA + 1;
          }
        } else if (op.subtype == OpSubtype::J1_Plaquette) {
          // Deterministic diagonal pairing: sites (0,2) and (1,3)
          {
            int liA = base + 0, liB = base + 4; // 2*0, 2*2
            vertex_partner[liA] = liB;
            vertex_partner[liB] = liA;
            vertex_partner[liA + 1] = liB + 1;
            vertex_partner[liB + 1] = liA + 1;
          }
          {
            int liA = base + 2, liB = base + 6; // 2*1, 2*3
            vertex_partner[liA] = liB;
            vertex_partner[liB] = liA;
            vertex_partner[liA + 1] = liB + 1;
            vertex_partner[liB + 1] = liA + 1;
          }
        } else {
          // J2/J3: two-site bond — pair sites 0 and 1
          vertex_partner[base + 0] = base + 2;
          vertex_partner[base + 2] = base + 0;
          vertex_partner[base + 1] = base + 3;
          vertex_partner[base + 3] = base + 1;
        }
      }
    }
  }

  // Trim to actual usage
  linked_leg.resize(leg_counter);
  vertex_partner.resize(leg_counter);
  vertex_spin.resize(leg_counter);
  leg_site.resize(leg_counter);
  linked_leg_to_op_index.resize(leg_counter);
  base_vertex_leg.resize(leg_counter);

  // Close periodic boundary conditions
  for (int s = 0; s < N; ++s) {
    if (first_leg[s] != -1) {
      linked_leg[last_leg[s]] = first_leg[s];
      linked_leg[first_leg[s]] = last_leg[s];
    }
  }
  n_legs = leg_counter;
}

// Check linked list basic integrity
void SSEUpdates::check_integrity() {
  const int n = static_cast<int>(linked_leg.size());
  for (int l = 0; l < n; ++l) {
    int next = linked_leg[l];
    if (next < 0 || next >= n) {
      std::cerr << "Integrity Error: Leg " << l << " links to invalid " << next
                << "\n";
      std::abort();
    }
    if (linked_leg[next] != l) {
      std::cerr << "Integrity Error: Asymmetry at leg " << l << " -> " << next
                << " -> " << linked_leg[next] << "\n";
      std::abort();
    }
  }
}

// =====================================================
// Propagate spins after cluster update
// ======================================================
void SSEUpdates::propagate_spins(SSEConfig &cfg) {
  const int N = static_cast<int>(cfg.spins.spin.size());
  for (int site = 0; site < N; ++site) {
    const int leg0 = first_leg[site];
    if (leg0 == -1)
      continue;
    cfg.spins.spin[site] = vertex_spin[leg0];
  }
}

void SSEUpdates::check_spin_consistency(const SSEConfig &cfg,
                                        const Geometry & /*geom*/) {
  std::vector<int> current_spins = cfg.spins.spin;
  const std::vector<int> &initial_spins = cfg.spins.spin;

  for (int p = 0; p < cfg.max_ops; ++p) {
    const auto &op = cfg.op_string[p];
    if (op.type == OpType::OffDiagonal &&
        op.subtype == OpSubtype::TransverseX_offdiag && op.index >= 0 &&
        op.index < static_cast<int>(current_spins.size())) {
      current_spins[op.index] *= -1;
    }
  }

  for (size_t i = 0; i < current_spins.size(); ++i) {
    if (current_spins[i] != initial_spins[i]) {
      std::cerr << "Spin consistency check failed at site " << i << "\n";
      std::exit(EXIT_FAILURE);
    }
  }
}

void SSEUpdates::adjust_operator_string(SSEConfig &cfg) {
  // Sandvik scaling: new_L = n + n/3
  int new_L = cfg.n_ops + cfg.n_ops / 3;
  if (new_L < 20)
    new_L = 20;

  if (new_L > cfg.max_ops) {
    cfg.op_string.resize(new_L);
    for (int p = cfg.max_ops; p < new_L; ++p) {
      cfg.op_string[p].type = OpType::Identity;
      cfg.op_string[p].subtype = OpSubtype::None;
      cfg.op_string[p].index = -1;
    }
    cfg.max_ops = new_L;
  }
}

// ======================================================
// 9. SLMC global proposal (skeleton)
// ======================================================
void SSEUpdates::slmc_update(SSEConfig &cfg, const Parameters & /*prm*/,
                             const Geometry & /*geom*/, RNG &rng) {
  (void)cfg;
  (void)rng;
}

// ======================================================
// 10. Gauge-sector update (row / column flips)
// ======================================================
// Flipping every spin in row r touches exactly the 2 sites of each black
// plaquette that straddles row r (its bottom-row sites and its top-row sites
// each come in pairs within that plaquette). Flipping both → prod unchanged.
// Same argument holds for column flips.  Both moves always satisfy detailed
// balance (w_new = w_old) so they are accepted unconditionally.
// This restores ergodicity that the BFS cluster update loses: at large β the
// cluster collapses to one giant cluster, allowing only the global Z2 flip.
void SSEUpdates::gauge_update(SSEConfig &cfg, const Parameters &prm,
                              const Geometry &geom, RNG &rng) {
  // Only valid when J0 is the only interaction (J2/J3/hx add bonds that
  // cross rows/columns and would change weight on a single-site flip).
  if (prm.J2 != 0.0 || prm.J3 != 0.0 || prm.hx != 0.0)
    return;

  const int Lx = prm.Lx;
  const int Ly = prm.Ly;

  // Attempt an independent row flip for each row
  for (int r = 0; r < Ly; ++r) {
    if (rng.uniform() < 0.5) {
      for (int x = 0; x < Lx; ++x)
        cfg.spins.spin[x + Lx * r] *= -1;
    }
  }

  // Attempt an independent column flip for each column
  for (int c = 0; c < Lx; ++c) {
    if (rng.uniform() < 0.5) {
      for (int y = 0; y < Ly; ++y)
        cfg.spins.spin[c + Lx * y] *= -1;
    }
  }
}

// ======================================================
// Choose exit leg for Directed Loop
// ======================================================
int SSEUpdates::choose_exit_leg(int vertex, int entry_leg, int spin_out,
                                const SSEConfig &cfg, const Parameters &prm,
                                const Geometry &geom, RNG &rng) {
  (void)spin_out;

  const auto &op = cfg.op_string[vertex];

  if (op.subtype == OpSubtype::TransverseX_diag ||
      op.subtype == OpSubtype::TransverseX_offdiag) {
    return (rng.uniform() < 0.5) ? (entry_leg ^ 1) : entry_leg;
  }

  const int base = base_vertex_leg[entry_leg];
  const int rel = entry_leg - base;
  int partner_site_idx = -1;

  if (op.subtype == OpSubtype::J0_Plaquette) {
    // Lookup table replaces 12-branch if-else
    static constexpr int J0_partners[3][4] = {
        {1, 0, 3, 2}, // type 0: 0<->1, 2<->3
        {2, 3, 0, 1}, // type 1: 0<->2, 1<->3
        {3, 2, 1, 0}, // type 2: 0<->3, 1<->2
    };
    partner_site_idx = J0_partners[vertex_decomp[vertex]][rel / 2];
  } else if (op.subtype == OpSubtype::J1_Plaquette) {
    static constexpr int J1_partners[4] = {2, 3, 0, 1}; // 0<->2, 1<->3
    const int partner_site = J1_partners[rel / 2];
    return base + partner_site * 2 + (rel % 2 == 0 ? 1 : 0);
  } else if (op.subtype == OpSubtype::J2_Dipole ||
             op.subtype == OpSubtype::J3_Inter) {
    partner_site_idx = (rel / 2 == 0) ? 1 : 0;
  }

  if (partner_site_idx == -1)
    return entry_leg ^ 1;

  // All current Ising-diagonal operators have w_horiz = 0, so evaluate
  // w_vert and fall through to vertical exit unless weights change.
  const int p_base = partner_site_idx * 2;
  const int l0 = base + (rel / 2) * 2;
  const int l2 = base + p_base;
  const int l3 = base + p_base + 1;

  double J_val = 0.0, C = 0.0;

  if (op.subtype == OpSubtype::J0_Plaquette) {
    const int s_idx1 = (rel / 2), s_idx2 = partner_site_idx;
    int spec1 = -1, spec2 = -1;
    for (int k = 0; k < 4; ++k) {
      if (k != s_idx1 && k != s_idx2) {
        if (spec1 == -1)
          spec1 = k;
        else
          spec2 = k;
      }
    }
    const double sign = static_cast<double>(vertex_spin[base + spec1 * 2] *
                                            vertex_spin[base + spec2 * 2]);
    J_val = prm.J0 * sign;
    C = 2.0 * std::abs(prm.J0);
  } else if (op.subtype == OpSubtype::J1_Plaquette) {
    J_val = prm.J1;
    C = 2.0 * std::abs(prm.J1);
  } else if (op.subtype == OpSubtype::J2_Dipole) {
    J_val = geom.j2bonds[op.index].J_val * prm.J2;
    C = 2.0 * std::abs(prm.J2);
  } else if (op.subtype == OpSubtype::J3_Inter) {
    J_val = prm.J3;
    C = 2.0 * std::abs(prm.J3);
  }

  double w_vert =
      C > 0.0 ? C - J_val * (vertex_spin[l0] * vertex_spin[l2])
              : std::abs(J_val) - J_val * (vertex_spin[l0] * vertex_spin[l2]);
  if (w_vert < 1e-9)
    w_vert = 0.0;

  // w_horiz = 0 for all Ising-diagonal operators; horizontal exit is never
  // chosen. Return the candidate horizontal leg here for completeness if
  // w_horiz ever becomes non-zero.
  const int cand_horiz = (rel % 2 == 0) ? l2 : l3;
  (void)cand_horiz;

  if (w_vert <= 1e-14)
    return entry_leg ^ 1; // fallback
  return entry_leg ^ 1;
}
