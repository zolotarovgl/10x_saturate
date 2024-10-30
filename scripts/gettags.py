import argparse
import pysam
import csv
import os
from concurrent.futures import ProcessPoolExecutor
import subprocess

def checkbam(input_bam,n_check = 10000000):
    with pysam.AlignmentFile(input_bam, "rb") as in_bam:
        print(f'Checking whether {input_bam} contains any secondary alignments or unmapped reads.' )
        for i,read in enumerate(in_bam):
            if i < n_check:
                # Keep only primary alignments (those without the secondary alignment flag 0x100)
                if read.is_secondary or read.is_unmapped:
                    return(True)
            else:
                return(False)
        else:
            return(False)
        
def keep_primary(input_bam,output_bam,num_threads = 1):
    #with pysam.AlignmentFile(input_bam, "rb") as in_bam, pysam.AlignmentFile(output_bam, "wb", header=in_bam.header) as out_bam:
    #    for read in in_bam:
    #        # Keep only primary alignments (those without the secondary alignment flag 0x100)
    #        if not read.is_secondary:
    #            out_bam.write(read)
    # Define the command as a list
    cmd = ["samtools", "view","-@", str(num_threads) , "-h", "-F", "0x100", input_bam, "-o", output_bam]
    # Run the command
    print(f'Running: {" ".join(cmd)}')
    subprocess.run(cmd, check=True)
    cmd = ["samtools", "index", output_bam]
    # Run the command
    print(f'Running: {" ".join(cmd)}')
    subprocess.run(cmd, check=True)

def process_chromosome(bam_file, chrom, output_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        for read in bam.fetch(chrom):
            # cell barcode
            if read.has_tag("CB"):
                cb = read.get_tag("CB")
            else:
                cb = None
            # UMI
            if read.has_tag("UB"):
                ub = read.get_tag("UB")
            elif read.has_tag("pN"):
                ub = read.get_tag("pN")
            else:
                ub = None
            # tag
            if read.has_tag("xf"):
                xf = read.get_tag("xf")
            else:
                xf = None
            # gene
            if read.has_tag("gx"):
                gx = read.get_tag("gx")
            else:
                gx = None
            writer.writerow([cb, xf, ub, gx])
        #print(f"done {chrom}.")
    bam.close()

def main(bam_file, output_prefix, num_threads):
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_prefix):
        os.makedirs(output_prefix)

    # Check if the .bam contains only primary alignments - this is important 
        
    clean_bam = '%s/clean.bam' % output_prefix
    status = checkbam(bam_file)
    if(status):
        print(f'{bam_file} contains secondary alignments or unmapped reads. Filtering ... => {clean_bam}')
        keep_primary(bam_file,clean_bam,num_threads)
        print(f'Created: {clean_bam}')
        bam_file = clean_bam
    else:
        print(f'{bam_file} contains only primary alignments.')

    
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
