import argparse
import pysam
import csv
import os
from concurrent.futures import ProcessPoolExecutor

def process_chromosome(bam_file, chrom, output_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        for read in bam.fetch(chrom):
            if read.has_tag("CB") and read.has_tag("xf") and read.has_tag("UB") and read.has_tag("GX"):
                cb = read.get_tag("CB")
                xf = read.get_tag("xf")
                ub = read.get_tag("UB")
                gx = read.get_tag("GX")
            else:
                cb = None
                xf = None
                ub = None
                gx = None
            writer.writerow([cb, xf, ub, gx])
        #print(f"done {chrom}.")
    bam.close()

def main(bam_file, output_prefix, num_threads):
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_prefix):
        os.makedirs(output_prefix)

    bam = pysam.AlignmentFile(bam_file, "rb")
    chromosomes = bam.references
    bam.close()

    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        for chrom in chromosomes:
            output_file = os.path.join(output_prefix, f"{chrom}.tsv")
            futures.append(executor.submit(process_chromosome, bam_file, chrom, output_file))
        for future in futures:
            future.result()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract tags from BAM file and save to individual chromosome TSV files.")
    parser.add_argument('bam_file', type=str, help='Path to the BAM file.')
    parser.add_argument('output_prefix', type=str, help='Output directory prefix.')
    parser.add_argument('--num_threads', type=int, default=1, help='Number of threads to use. Default is 1.')

    args = parser.parse_args()

    main(args.bam_file, args.output_prefix, args.num_threads)
    print("Tags have been extracted and saved to individual chromosome TSV files.")
