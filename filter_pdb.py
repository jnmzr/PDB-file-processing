import os
import argparse

def filter_pdb_by_residues(pdb_file, output_file, residue_file=None, renumber=False):
    # Read residue number file if provided
    target_residues = None
    if residue_file:
        with open(residue_file, 'r') as f:
            target_residues = set(int(line.strip()) for line in f if line.strip())
        print(f"Target residues to keep: {sorted(target_residues)}")
    
    # Read PDB file and process
    with open(pdb_file, 'r') as f:
        pdb_lines = f.readlines()
    
    output_lines = []
    in_model = False
    model_atoms = []
    original_residues = []  # Store original residue IDs in order
    model_count = 0
    
    for line in pdb_lines:
        # Skip all REMARK, TITLE, CRYST1 lines regardless of position
        if line.startswith(('REMARK', 'TITLE', 'CRYST1')):
            continue
            
        if line.startswith('MODEL'):
            # Start new MODEL
            in_model = True
            model_atoms = []
            model_count += 1
            output_lines.append(line)
        elif line.startswith('ENDMDL'):
            # Process current MODEL atoms
            in_model = False
            
            # Collect all CA atoms and corresponding residue IDs for current MODEL
            ca_atoms = []
            current_model_residues = []
            
            for atom_line in model_atoms:
                if atom_line.startswith('ATOM') and atom_line[12:16].strip() == 'CA':
                    residue_num = int(atom_line[22:26].strip())
                    # If target_residues is provided, only keep specified residues
                    # Otherwise, keep all CA atoms
                    if target_residues is None or residue_num in target_residues:
                        ca_atoms.append(atom_line)
                        current_model_residues.append(residue_num)
            
            # If first MODEL, record original residue ID order
            if model_count == 1:
                original_residues = current_model_residues.copy()
                print(f"Original residue ID order: {original_residues}")
                print(f"Number of residues: {len(original_residues)}")
            
            # Renumber atoms and residues
            filtered_atoms = []
            atom_counter = 1
            
            for atom_line in ca_atoms:
                if atom_line.startswith('ATOM'):
                    # Parse original residue number
                    orig_residue_num = int(atom_line[22:26].strip())
                    
                    # Update atom number (columns 7-11)
                    new_line = atom_line[:6] + f"{atom_counter:>5}" + atom_line[11:]
                    
                    # Update residue number
                    if renumber:
                        # Renumber from 1
                        new_residue_num = original_residues.index(orig_residue_num) + 1
                    else:
                        # Keep original residue numbers
                        new_residue_num = orig_residue_num
                    
                    new_line = new_line[:22] + f"{new_residue_num:>4}" + new_line[26:]
                    
                    filtered_atoms.append(new_line)
                    atom_counter += 1
                else:
                    filtered_atoms.append(atom_line)
            
            # Add TER record
            if filtered_atoms:
                # Get last atom information to create TER record
                last_atom = filtered_atoms[-1]
                if last_atom.startswith('ATOM'):
                    ter_line = f"TER   {atom_counter:>5}      {last_atom[17:20]} {last_atom[21:26]}\n"
                    filtered_atoms.append(ter_line)
            
            output_lines.extend(filtered_atoms)
            output_lines.append(line)
        elif in_model:
            # Inside MODEL, collect atom information
            if line.startswith('ATOM') or line.startswith('TER'):
                model_atoms.append(line)
            else:
                # Skip other lines inside MODEL (including REMARK, TITLE, CRYST1)
                if not line.startswith(('REMARK', 'TITLE', 'CRYST1')):
                    output_lines.append(line)
        else:
            # Lines outside MODEL are kept directly (except REMARK, TITLE, CRYST1)
            if not line.startswith(('REMARK', 'TITLE', 'CRYST1')):
                output_lines.append(line)
    
    # Write output file
    with open(output_file, 'w') as f:
        f.writelines(output_lines)
    
    # Write original residue ID record file
    ori_resid_file = "ori_resid.log"
    with open(ori_resid_file, 'w') as f:
        for resid in original_residues:
            f.write(f"{resid}\n")
    
    print(f"Processing completed! Results saved to: {output_file}")
    print(f"Original residue IDs recorded in: {ori_resid_file}")
    if renumber:
        print(f"Residue numbering mapping: {dict(zip(range(1, len(original_residues) + 1), original_residues))}")
    else:
        print("Residue numbers kept as original")

def main():
    parser = argparse.ArgumentParser(description='Filter PDB file by residue numbers')
    parser.add_argument('-filter', help='File containing residue numbers to keep (optional)')
    parser.add_argument('-renum', action='store_true', help='Renumber residues starting from 1')
    parser.add_argument('-i', '--input', required=True, help='Input PDB file')
    parser.add_argument('-o', '--output', required=True, help='Output PDB file')
    
    args = parser.parse_args()
    
    # Check if files exist
    if not os.path.exists(args.input):
        print(f"Error: PDB file '{args.input}' does not exist")
        return
    
    if args.filter and not os.path.exists(args.filter):
        print(f"Error: Residue file '{args.filter}' does not exist")
        return
    
    filter_pdb_by_residues(args.input, args.output, args.filter, args.renum)

if __name__ == "__main__":
    main()