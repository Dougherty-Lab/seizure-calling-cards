import sys

# Input and output files
input_file = "peaks.bed"
output_file = "peaks_shifted_20kb.bed"

# Open input and output files
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        fields = line.strip().split('\t')
        
        # Ensure we have at least 3 columns (chr, start, end)
        if len(fields) >= 3:
            # Convert start and end to integers
            start = int(fields[1])
            end = int(fields[2])
            
            # Shift right by 10kb (10000 bp)
            new_start = start + 20000
            new_end = end + 20000
            
            # Update fields
            fields[1] = str(new_start)
            fields[2] = str(new_end)
            
            # Write updated line to output file
            outfile.write('\t'.join(fields) + '\n')
        else:
            # Copy line as is if it doesn't have enough columns
            outfile.write(line)

print(f"Created shifted BED file: {output_file}")