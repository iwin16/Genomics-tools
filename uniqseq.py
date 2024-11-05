from Bio import SeqIO
import argparse
import sys
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

def process_chunk(records):
    unique_sequences = set()
    output_records = []

    for record in records:
        sequence = str(record.seq)
        if sequence not in unique_sequences:
            unique_sequences.add(sequence)
            output_records.append(record)

    return output_records

def remove_redundant_sequences(input_fastq, output_fastq, chunk_size=10000, num_threads=1):
    if not os.path.isfile(input_fastq):
        sys.exit(f"Error: The input file '{input_fastq}' does not exist.")

    all_unique_records = []
    total_input_sequences = 0
    all_chunks = []
    start_time = time.time()

    # Read and split the FASTQ file into chunks
    with open(input_fastq, "r") as handle:
        records = SeqIO.parse(handle, "fastq")
        chunk = []
        for record in records:
            chunk.append(record)
            total_input_sequences += 1
            if len(chunk) >= chunk_size:
                all_chunks.append(chunk)
                chunk = []
        if chunk:
            all_chunks.append(chunk)

    # Process each chunk in parallel and collect results
    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        future_to_chunk = {executor.submit(process_chunk, chunk): chunk for chunk in all_chunks}
        for future in as_completed(future_to_chunk):
            all_unique_records.extend(future.result())

    # Write final unique sequences to output file
    with open(output_fastq, "w") as output_handle:
        SeqIO.write(all_unique_records, output_handle, "fastq")

    total_time_taken = time.time() - start_time
    total_unique_sequences = len(set(str(record.seq) for record in all_unique_records))
    percent_unique_sequences = (total_unique_sequences / total_input_sequences * 100) if total_input_sequences else 0

    print("\nUnique sequences stats:")
    print(f"Total input sequences: {total_input_sequences}")
    print(f"Total unique sequences: {total_unique_sequences}")
    print(f"Percent unique sequences: {percent_unique_sequences:.2f}%")
    print(f"Total time taken: {total_time_taken:.2f} seconds")

    return total_unique_sequences

def main():
    parser = argparse.ArgumentParser(description="Remove redundant sequences from a FASTQ file.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTQ file with unique sequences")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use (default: 1)")
    args = parser.parse_args()

    # Validate input file extension
    if not args.input.lower().endswith((".fastq", ".fq")):
        sys.exit("Error: The input file must have a .fastq or .fq extension.")

    # Ensure output directory exists
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Run the sequence removal
    remove_redundant_sequences(args.input, args.output, num_threads=args.threads)

if __name__ == "__main__":
    main()
