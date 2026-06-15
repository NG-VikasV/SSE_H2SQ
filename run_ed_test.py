import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'backend', 'Analysis', 'ED_analysis'))
from exact_diagonalzation_H2SQ import H2SQ_ED, generate_black_plaquette_indices_4_4

Lx, Ly, Nl = 4, 4, 1
J0, J1, J2, h = 1.0, 0.3, 0.025, 0.0
bplaqs = generate_black_plaquette_indices_4_4(Lx, Ly)
W, V = H2SQ_ED(Lx, Ly, Nl, J0, J1, J2, h, bplaqs, 1)
print("Ground state energy:", W[0])
print("Energy per site:", W[0]/(Lx*Ly*Nl))
