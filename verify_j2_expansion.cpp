#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <iomanip>

// Simple mock of Geometry build for testing J2 expansion
int main() {
    int Lx = 4;
    int Ly = 4;
    int Lz = 1;
    int Ns = Lx * Ly * Lz;

    struct Plaq {
        int id;
        int sites[4];
    };
    std::vector<Plaq> blackPlaquettes;
    
    // Build Grid & Plaquettes (Same as variables.cpp)
    for (int y = 0; y < Ly; ++y) {
        for (int x = 0; x < Lx; ++x) {
            int i = x + Lx * y;
            // Sites BL, TL, TR, BR (indices)
            int s0 = i;
            int s1 = x + Lx * ((y + 1) % Ly);
            int s2 = ((x + 1) % Lx) + Lx * ((y + 1) % Ly);
            int s3 = ((x + 1) % Lx) + Lx * y;

            if ((x + y) % 2 == 0) { // Black
                blackPlaquettes.push_back({(int)blackPlaquettes.size(), {s0, s1, s2, s3}});
            }
        }
    }

    std::cout << "Lattice " << Lx << "x" << Ly << ". Black Plaquettes: " << blackPlaquettes.size() << "\n";

    // J2 Pairs Verification
    // Logic: Diagonal neighbors. (x,y) and (x+1, y+1) etc.
    // Let's brute force distances or use the grid logic.
    // Here we assume variables.cpp logic finds "Neighbors".
    // Neighbors of a black plaquette at (x,y) [even sum]
    // are at (x+1, y+1), (x-1, y+1), etc.
    // These are ALSO black plaquettes (even sum + even shift = even).
    // They share ONE site (corner).
    
    // Let's iterate all pairs and check site overlaps.
    // If overlap count == 1, they are diagonal neighbors?
    // If overlap count == 2, they are Nearest (share edge) -> No, Checkerboard doesn't share edges between same color.
    // Correct. Black plaquettes only share corners with other Black plaquettes.
    
    int total_connected_pairs = 0;
    double total_constant = 0.0;
    
    // Coefficients (from variables.cpp)
    double cx[4] = { 0.25, 0.25, -0.25, -0.25 };
    double cy[4] = { 0.25, -0.25, -0.25, 0.25 };

    std::cout << "\nScanning all Black-Black Pairs for Overlaps...\n";

    for(size_t i=0; i<blackPlaquettes.size(); ++i) {
        for(size_t j=i+1; j<blackPlaquettes.size(); ++j) {
            // Check overlapping sites
            std::vector<int> overlap;
            for(int k=0; k<4; ++k) {
                for(int m=0; m<4; ++m) {
                    if(blackPlaquettes[i].sites[k] == blackPlaquettes[j].sites[m]) {
                        overlap.push_back(blackPlaquettes[i].sites[k]);
                        
                        // Calculate Constant Contribution
                        // Term = J2 * (cx_k * cx_m + cy_k * cy_m) * S * S
                        // S*S = 1.
                        double term = (cx[k] * cx[m] + cy[k] * cy[m]);
                        total_constant += term;
                        
                        // std::cout << "  Overlap Pair (" << i << "," << j << ") at Site " << blackPlaquettes[i].sites[k] 
                        //           << ". Term=" << term << "\n";
                    }
                }
            }
            
            if(!overlap.empty()) {
                total_connected_pairs++;
            }
        }
    }

    std::cout << "Total Interacting Pairs: " << total_connected_pairs << "\n";
    std::cout << "Computed total_constant (sum of coefficients for self-terms): " << total_constant << "\n";
    
    // Expected?
    // Each pair shares 1 site usually?
    // 8 plaquettes on 4x4 torus.
    // Neighbors per plaquette: 4 diagonal neighbors?
    // Total pairs = 8 * 4 / 2 = 16 pairs.
    // Each overlap contributes some constant.
    
    return 0;
}
