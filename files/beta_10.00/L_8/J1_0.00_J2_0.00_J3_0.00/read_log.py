try:
    # Try utf-16 based on BOM
    with open('output_log.txt', 'r', encoding='utf-16', errors='ignore') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if i > 500: break
            print(f"{i}: {line.strip()}")
except Exception as e:
    print(f"Error: {e}")
