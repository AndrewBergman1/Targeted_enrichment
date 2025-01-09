import sys

# Get command-line arguments
baits_file = sys.argv[1]
problematic_baits_file = sys.argv[2]
filtered_baits_file = sys.argv[3]

# Open the files
with open(baits_file, "r") as baits, \
     open(problematic_baits_file, "r") as problematic_baits, \
     open(filtered_baits_file, "w") as filtered_baits:

    # Read problematic baits into a set for faster lookup
    problematic_baits_set = {line.strip() for line in problematic_baits}

    # Process each bait
    for line in baits:
        if line.startswith(">"):
            # Extract the bait name (assumes the bait name is on the same line)
            bait_name = line[1:].strip()

            print("Processing", bait_name)

            # Check if the bait is in the problematic baits set
            if bait_name not in problematic_baits_set:
                # Write the header to filtered_baits
                filtered_baits.write(line)

                # Read the next line (the sequence) and write it to filtered_baits
                sequence = next(baits)
                filtered_baits.write(sequence)

print("Filtering complete. Filtered baits saved to", filtered_baits_file)
