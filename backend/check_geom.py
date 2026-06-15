import json
with open('c:/Users/VikasVijigiri/Documents/SSE_H2SQ/backend/src/geom.json', 'r') as f:
    d = json.load(f)
print("Black Plaquettes:", len(d['black_plaquettes']))
print("J2 bonds:", len(d['j2_bonds']))
print("J3 bonds:", len(d['j3_bonds']))
print("Ns:", d['Lx'] * d['Ly'] * d['Lz'])
