#!/usr/bin/env python3
import sys
import os

def print_help():
    help_text = """
Usage:
    python update_coords.py <relax.out> <target_1.in> [target_2.in ...]

Description:
    Reads a Quantum ESPRESSO relaxation output file (e.g. relax.out) to find the final
    atomic coordinates and cell parameters. Then updates the specified target input files
    with these new values.

    The script preserves atom-specific constraints (e.g. '0 0 0') in the target files
    if the number of atoms matches.
"""
    print(help_text)

def parse_atomic_positions(lines):
    block = []
    started = False
    for line in lines:
        stripped = line.strip()
        if not stripped:
            if started and len(block) > 1:
                break
            continue
        
        if not started:
            if "ATOMIC_POSITIONS" in stripped.upper():
                started = True
                block.append(line)
        else:
            parts = stripped.split()
            if len(parts) >= 4:
                try:
                    float(parts[1])
                    float(parts[2])
                    float(parts[3])
                    block.append(line)
                except ValueError:
                    break
            else:
                break
    return block

def parse_cell_parameters(lines):
    block = []
    started = False
    for line in lines:
        stripped = line.strip()
        if not stripped:
            if started and len(block) > 1:
                break
            continue
        
        if not started:
            if "CELL_PARAMETERS" in stripped.upper():
                started = True
                block.append(line)
        else:
            parts = stripped.split()
            if len(parts) >= 3:
                try:
                    float(parts[0])
                    float(parts[1])
                    float(parts[2])
                    block.append(line)
                    if len(block) == 4:  # header + 3 lattice vector lines
                        break
                except ValueError:
                    break
            else:
                break
    return block

def extract_final_blocks(relax_out_path):
    if not os.path.exists(relax_out_path):
        print(f"Error: Output file '{relax_out_path}' not found.")
        sys.exit(1)
        
    with open(relax_out_path, 'r') as f:
        lines = f.readlines()
        
    final_lines = []
    in_final = False
    for line in lines:
        if "Begin final coordinates" in line:
            in_final = True
            final_lines = []
        elif "End final coordinates" in line:
            in_final = False
        elif in_final:
            final_lines.append(line)
            
    atomic_positions = []
    cell_parameters = []
    
    if final_lines:
        # Parse from the final coordinates section
        atomic_positions = parse_atomic_positions(final_lines)
        cell_parameters = parse_cell_parameters(final_lines)
        print("Extracted final coordinates from 'Begin final coordinates' block.")
    else:
        # Fallback to scanning the whole file backwards
        print("Warning: 'Begin final coordinates' block not found. Scanning file backwards for last coordinates.")
        for i in range(len(lines) - 1, -1, -1):
            if "ATOMIC_POSITIONS" in lines[i].upper() and not atomic_positions:
                atomic_positions = parse_atomic_positions(lines[i:])
            if "CELL_PARAMETERS" in lines[i].upper() and not cell_parameters:
                cell_parameters = parse_cell_parameters(lines[i:])
            if atomic_positions and cell_parameters:
                break
                
    if not atomic_positions:
        print("Error: Could not find any ATOMIC_POSITIONS in output file.")
        sys.exit(1)
        
    return atomic_positions, cell_parameters

def merge_atomic_positions(old_block, new_block):
    old_coords = old_block[1:]
    new_coords = new_block[1:]
    
    if len(old_coords) != len(new_coords):
        print(f"  Note: Atom counts differ between target ({len(old_coords)}) and output ({len(new_coords)}). Direct block overwrite will be performed.")
        return new_block
        
    merged = [new_block[0]]
    for old_line, new_line in zip(old_coords, new_coords):
        old_parts = old_line.strip().split()
        new_parts = new_line.strip().split()
        
        # Keep constraints (any columns after the first 4: atom, x, y, z)
        if len(old_parts) > 4:
            constraints = old_parts[4:]
            formatted_line = f"  {new_parts[0]:<3} {new_parts[1]:>15} {new_parts[2]:>15} {new_parts[3]:>15}  " + "  ".join(constraints)
            merged.append(formatted_line + "\n")
        else:
            formatted_line = f"  {new_parts[0]:<3} {new_parts[1]:>15} {new_parts[2]:>15} {new_parts[3]:>15}"
            merged.append(formatted_line + "\n")
    return merged

def replace_block_in_content(content_lines, keyword, new_block):
    if not new_block:
        return content_lines, False
        
    start_idx = -1
    for i, line in enumerate(content_lines):
        if keyword in line.upper():
            start_idx = i
            break
            
    if start_idx == -1:
        return content_lines, False
        
    end_idx = len(content_lines)
    for j in range(start_idx + 1, len(content_lines)):
        line_upper = content_lines[j].upper().strip()
        if not line_upper:
            end_idx = j
            break
        if line_upper.startswith("&") or line_upper == "/":
            end_idx = j
            break
            
        is_other_card = False
        for card in ["ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS", "OCCUPATIONS", "CONSTRAINTS"]:
            if line_upper.startswith(card):
                is_other_card = True
                break
        if is_other_card:
            end_idx = j
            break
            
    old_block = content_lines[start_idx:end_idx]
    
    if keyword == "ATOMIC_POSITIONS":
        final_block = merge_atomic_positions(old_block, new_block)
    else:
        final_block = new_block
        
    formatted_block = [line if line.endswith('\n') else line + '\n' for line in final_block]
    new_content_lines = content_lines[:start_idx] + formatted_block + content_lines[end_idx:]
    return new_content_lines, True

def main():
    if len(sys.argv) < 3 or sys.argv[1] in ('-h', '--help'):
        print_help()
        sys.exit(1)
        
    relax_out = sys.argv[1]
    targets = sys.argv[2:]
    
    # Extract blocks
    atomic_positions, cell_parameters = extract_final_blocks(relax_out)
    
    print("\n--- Extracted ATOMIC_POSITIONS ---")
    sys.stdout.write("".join(atomic_positions))
    if cell_parameters:
        print("\n--- Extracted CELL_PARAMETERS ---")
        sys.stdout.write("".join(cell_parameters))
    print("----------------------------------\n")
    
    for target in targets:
        if not os.path.exists(target):
            print(f"Error: Target file '{target}' not found. Skipping.")
            continue
            
        with open(target, 'r') as f:
            content_lines = f.readlines()
            
        # Replace atomic positions
        content_lines, ap_replaced = replace_block_in_content(content_lines, "ATOMIC_POSITIONS", atomic_positions)
        
        # Replace cell parameters (only if they were found in the output)
        cp_replaced = False
        if cell_parameters:
            content_lines, cp_replaced = replace_block_in_content(content_lines, "CELL_PARAMETERS", cell_parameters)
            
        if ap_replaced:
            print(f"Updated ATOMIC_POSITIONS in '{target}'")
        else:
            print(f"Warning: ATOMIC_POSITIONS not found in '{target}'")
            
        if cell_parameters:
            if cp_replaced:
                print(f"Updated CELL_PARAMETERS in '{target}'")
            else:
                print(f"Note: CELL_PARAMETERS not found in '{target}', skipping cell update.")
                
        # Write back to target
        with open(target, 'w') as f:
            f.writelines(content_lines)
            
    print("\nDone!")

if __name__ == "__main__":
    main()
