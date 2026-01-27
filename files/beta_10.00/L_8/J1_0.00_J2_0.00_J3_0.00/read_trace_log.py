try:
    with open('output_winding_trace.txt', 'r', encoding='latin-1', errors='ignore') as f:
        lines = f.readlines()
        
        trace = []
        current_step = {}
        
        for line in lines:
            line = line.strip()
            if "Step: " in line:
                # Parse: Step: 1 Cur: 6 Start: ...
                parts = line.split()
                try:
                    s_idx = parts.index("Step:")
                    step_num = int(parts[s_idx+1])
                    current_step['step'] = step_num
                    
                    c_idx = parts.index("Cur:")
                    current_step['cur'] = int(parts[c_idx+1])
                except: pass
                
            if "NextV: " in line:
                # Parse: NextV: 8 Op: 1
                parts = line.split()
                try:
                    n_idx = parts.index("NextV:")
                    current_step['next'] = int(parts[n_idx+1])
                    if "Op:" in parts:
                        o_idx = parts.index("Op:")
                        current_step['op'] = int(parts[o_idx+1])
                    
                    if 'step' in current_step:
                         trace.append(current_step.copy())
                         current_step = {}
                except: pass
                
        # Sort and print
        trace.sort(key=lambda x: x.get('step', 0))
        for t in trace: # Print all found
            op_str = f" (Op {t.get('op')})" if 'op' in t else ""
            print(f"Step {t.get('step')}: Leg {t.get('cur')} -> VertexLeg {t.get('next')}{op_str}")
                
except Exception as e:
    print(f"Error: {e}")
