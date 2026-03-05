import argparse

def write_chimerax_pairings(dot_bracket, filename="pairings.txt"):
    stack = []
    pairs = []
    
    # 1. Extract all individual base pairs
    for i, char in enumerate(dot_bracket):
        if char == '(':
            stack.append(i + 1)
        elif char == ')':
            if stack:
                start = stack.pop()
                end = i + 1
                pairs.append((start, end))
    
    # Sort pairs by the start index to find contiguous stems
    pairs.sort()

    # 2. Group individual pairs into "stems" (triples: start, end, length)
    stems = []
    if not pairs:
        return "No pairs found."

    current_start, current_end = pairs[0]
    length = 1

    for i in range(1, len(pairs)):
        next_start, next_end = pairs[i]
        # In a standard stem, as start increases by 1, end decreases by 1
        if next_start == current_start + length and next_end == current_end - length:
            length += 1
        else:
            stems.append((current_start, current_end, length))
            current_start, current_end = next_start, next_end
            length = 1
    
    stems.append((current_start, current_end, length))

    # 3. Write to file in 3-column format
    with open(filename, "w") as f:
        for s, e, l in stems:
            f.write(f"{s}\t{e}\t{l}\n")
    
    print(f"File '{filename}' created successfully with {len(stems)} stems.")
    return stems

def get_parser():
    parser = argparse.ArgumentParser(description="Convert dot-bracket notation to ChimeraX pairs format.")
    parser.add_argument("input", type=str, help="Path to the input file containing the dot-bracket notation.")
    parser.add_argument("output", type=str, help="Path to the output file where the ChimeraX pairs will be saved.")
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    
    # Read the dot-bracket notation from the input file
    with open(args.input, 'r') as f:
        dot_bracket = f.read().strip()
    
    # Convert to ChimeraX pairs format
    write_chimerax_pairings(dot_bracket, args.output)
    
if __name__ == "__main__":
    main()