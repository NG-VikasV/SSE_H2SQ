try:
    with open('output_bond_verification.txt', 'r', encoding='latin-1', errors='ignore') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if i >= 26 and i < 80:
                print(line.strip())
except Exception as e:
    print(f"Error: {e}")
