import sys

if len(sys.argv) < 3:
    print("Usage: python script.py <resultDB> <problematic_probes_file>")
    sys.exit(1)

resultDB = sys.argv[1]
problematic_probes_file = sys.argv[2]

threshold = 0.8 
problematic_probes = set()

try:
    with open(resultDB, "r") as file:
        for line in file:
            elements = line.strip().split("\t")
            if len(elements) > 3:  # Ensure there are enough elements
                try:
                    identity = float(elements[3])
                    if identity >= threshold:
                        print(f"Probe {elements[0]} exceeds the identity threshold of {threshold}")
                        problematic_probes.add(elements[0])
                except ValueError:
                    print(f"Invalid identity value in line: {line.strip()}")

    print(problematic_probes)

    with open(problematic_probes_file, "w") as saveFile:
        for probe in problematic_probes:
            saveFile.write(f"{probe}\n")

    print("Problematic probes are saved to", problematic_probes_file)

except FileNotFoundError:
    print("The input file was not found.")
except Exception as e:
    print(f"An error occurred: {e}")
