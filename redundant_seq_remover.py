import argparse
from Bio import SeqIO
import sys
import os

def remove_redundant_sequences(input_fastq, output_fastq):
    if not os.path.isfile(input_fastq):
        print(f"Error: The input file '{input_fastq}' does not exist.")
        sys.exit(1)

    unique_sequences = set()
    try:
        with open(output_fastq, "w") as output_handle:
            for record in SeqIO.parse(input_fastq, "fastq"):
                sequence = str(record.seq)
                if sequence not in unique_sequences:
                    unique_sequences.add(sequence)
                    SeqIO.write(record, output_handle, "fastq")
                    
        print(f"Redundant sequences removed. Output saved to {output_fastq}")

    except Exception as e:
        print(f"An error occurred while processing the file: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Remove redundant sequences from a FASTQ file.")
    parser.add_argument("-i", "--input", help="Input FASTQ file")
    parser.add_argument("-o", "--output", help="Output FASTQ file with unique sequences")
    args = parser.parse_args()
    
    # Check if both input and output files are provided
    if not args.input:
        print("Error: No input file specified. Use the -i option to specify the input FASTQ file.")
        sys.exit(1)
    if not args.output:
        print("Error: No output file specified. Use the -o option to specify the output FASTQ file.")
        sys.exit(1)
    
    # Check if the input file has .fastq or .fq extension
    if not (args.input.endswith(".fastq") or args.input.endswith(".fq")):
        print("Error: The input file must have a .fastq or .fq extension.")
        sys.exit(1)
    
    # Run the sequence removal function
    remove_redundant_sequences(args.input, args.output)

if __name__ == "__main__":
    main()
