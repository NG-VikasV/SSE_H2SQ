#ifndef SSE_UPDATES_HPP
#define SSE_UPDATES_HPP

#include <vector>
#include <stack>

// Forward declarations
struct SSEConfig;
struct Parameters;
struct Geometry;
class RNG;

// ======================================================
// SSE Update Engine
// ======================================================

class SSEUpdates {
public:
    SSEUpdates() = default;

    // --------------------------------------------------
    // Top-level Monte Carlo sweep
    // --------------------------------------------------
    void sweep(SSEConfig& cfg,
               const Parameters& prm,
               const Geometry& geom,
               RNG& rng);

    // --------------------------------------------------
    // Core SSE updates
    // --------------------------------------------------
    void adjust_operator_string(SSEConfig& cfg);

    void check_spin_consistency(const SSEConfig& cfg,
                                const Geometry& geom);




    // Diagonal insertion / removal
    void diagonal_update_J0(SSEConfig& cfg,
                         const Parameters& prm,
                         const Geometry& geom,
                         RNG& rng);

    void diagonal_update_J1(SSEConfig& cfg,
                         const Parameters& prm,
                         const Geometry& geom,
                         RNG& rng);

    void diagonal_update_h(SSEConfig& cfg,
                         const Parameters& prm,
                         const Geometry& geom,
                         RNG& rng);                         
                         
    void diagonal_update_J2(SSEConfig& cfg,
                         const Parameters& prm,
                         const Geometry& geom,
                         RNG& rng);

    void diagonal_update_J3(SSEConfig& cfg,
                         const Parameters& prm,
                         const Geometry& geom,
                         RNG& rng);                     

    // Directed-loop (quantum cluster) update
    void directed_loop_update(SSEConfig& cfg,
                              const Parameters& prm,
                              const Geometry& geom,
                              RNG& rng);

    // --------------------------------------------------
    // Optional: Self-learning Monte Carlo (global move)
    // --------------------------------------------------
    void slmc_update(SSEConfig& cfg,
                     const Parameters& prm,
                     const Geometry& geom,
                     RNG& rng);

private:
    // ==================================================
    // Internal helpers (hot path)
    // ==================================================

    // Number of legs per vertex
    static constexpr int legs_per_vertex = 4;

    // --------------------------------------------------
    // Vertex / leg bookkeeping
    // --------------------------------------------------



    void build_vertex_list(SSEConfig& cfg,
                        const Geometry& geom,
                        RNG& rng);


    void check_integrity();




    // Choose exit leg (directed-loop probabilities)
    // Choose exit leg (directed-loop probabilities)
    int choose_exit_leg(int vertex,
                        int entry_leg,
                        int spin_out,
                        const SSEConfig& cfg,
                        const Parameters& prm,
                        const Geometry& geom,
                        RNG& rng);

                        

    // Propagate changes back to |α⟩
    void propagate_spins(SSEConfig& cfg);


    // Flip a single vertex leg (global leg index)
    inline void flip_leg(int global_leg, int new_spin);
        
    // --------------------------------------------------
    // Directed-loop propagation
    // --------------------------------------------------
    void propagate_loop(SSEConfig& cfg,
                        const Parameters& prm,
                        const Geometry& geom,
                        RNG& rng,
                        int start_leg);


    // --------------------------------------------------
    // Diagonal operator weights
    // --------------------------------------------------
    double diagonal_weight(int subtype,
                           int index,
                           const SSEConfig& cfg,
                           const Parameters& prm,
                           const Geometry& geom);

    // ==================================================
    // Internal storage (reused every sweep)
    // ==================================================

    // List of non-identity operators as vertices
    std::vector<int> vertex_op_index;

    // Decomposition choice for each vertex (J0, J1 support)
    std::vector<int> vertex_decomp;

    // Linked list of legs (size = n_legs)
    std::vector<int> linked_leg;

    std::vector<int> linked_leg_to_op_index;

    // Partner leg within the vertex (Space Link)
    std::vector<int> vertex_partner;

    // Spin value on each leg
    std::vector<int> vertex_spin;

    std::vector<int> leg_site;  
    
    std::vector<int> base_vertex_leg;
    
    // cluster bookkeeping
    std::vector<int> cluster_state; 
    //  >=0 : unvisited
    //  -1  : visited (not flipped)
    //  -2  : flipped

    std::stack<int> stack;


    // imaginary-time boundaries
    std::vector<int> first_leg;          // per site
    std::vector<int> last_leg;           // per site



    int n_legs = 0;

    // --------------------------------------------------
    // SLMC toggle
    // --------------------------------------------------
    bool use_slmc = false;
};

#endif // SSE_UPDATES_HPP
