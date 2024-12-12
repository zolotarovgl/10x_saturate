import argparse
import pysam
import csv
import os
from concurrent.futures import ProcessPoolExecutor
import subprocess

def checkbam(input_bam,n_check = 10000000):
    # this script checks if there are unmapped reads in the bam file or secondary aligments 
    with pysam.AlignmentFile(input_bam, "rb") as in_bam:
        print(f'Checking whether {input_bam} contains any secondary alignments or unmapped reads.' )
        for i,read in enumerate(in_bam):
            if i < n_check:
                # Keep only primary alignments (those without the secondary alignment flag 0x100)
                if read.is_secondary:
                    print(f'Found secondary alignment in {input_bam}!')
                    return(True)
                if read.is_unmapped:
                    print(f'Found umapped read in {input_bam}!')
                    return(True)
            else:
                return(False)
        else:
            return(False)



def count_bam(input_bam, num_threads):
    # Keep only the primary alignments:
    cmd = ["samtools", "view", "-@", str(num_threads), "-c", input_bam]
    # Run the command and capture the output
    print(f'Running: {" ".join(cmd)}')
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    # Process the output
    output = result.stdout.strip()
    return int(output)

def get_flagstat(input_bam,n_check = 1000000):
    print(f'Counting read stats in {input_bam} using first {n_check} reads')
    stats = {
        "total_reads": 0,
        "mapped_reads": 0,
        "paired_reads": 0,
        "properly_paired_reads": 0,
        "unmapped_reads": 0,
        "secondary_alignments": 0,
        "supplementary_alignments": 0,
        "duplicate_reads": 0
    }

    # Open the BAM file
    with pysam.AlignmentFile(input_bam, "rb") as bamfile:
        for read in bamfile:
            if stats["total_reads"] < n_check:
                stats["total_reads"] += 1
                if read.is_unmapped:
                    stats["unmapped_reads"] += 1
                else:
                    stats["mapped_reads"] += 1

                if read.is_paired:
                    stats["paired_reads"] += 1

                if read.is_proper_pair:
                    stats["properly_paired_reads"] += 1

                if read.is_secondary:
                    stats["secondary_alignments"] += 1

                if read.is_supplementary:
                    stats["supplementary_alignments"] += 1

                if read.is_duplicate:
                    stats["duplicate_reads"] += 1
            else:
                break
    return stats



def filter_bam(input_bam,output_bam,num_threads = 1):
    # keep only the primary alignments:
    cmd = ["samtools", "view","-@", str(num_threads) , "-h", "-F", "260","-F","0x400", input_bam, "-o", output_bam]
    # Run the command
    print(f'Running: {" ".join(cmd)}')
    subprocess.run(cmd, check=True)
    cmd = ["samtools", "index","-@",str(num_threads), output_bam]
    # Run the command
    print(f'Running: {" ".join(cmd)}')
    subprocess.run(cmd, check=True)
    n_input = count_bam(input_bam,num_threads)
    n_output = count_bam(output_bam,num_threads)
    print(f'{input_bam}: {n_input}\n{output_bam}: {n_output}\nRemoved {n_input - n_output} reads.')


def index_bam(bam_file,num_threads = 1):
    if not os.path.exists(bam_file + '.bai'):
        cmd = ["samtools", "index","-@",str(num_threads), bam_file]
        # Run the command
        print(f'Running: {" ".join(cmd)}')
        subprocess.run(cmd, check=True)
    else:
        print(f'{bam_file} is indexed')

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
            # changed this to only output cb or ub
            #writer.writerow([cb, xf, ub, gx])
            writer.writerow([cb, ub])
        #print(f"done {chrom}.")
    bam.close()


def main(bam_file, output_prefix, num_threads):
    # Check if the bam file is indexed:
    index_bam(bam_file,num_threads)

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_prefix):
        os.makedirs(output_prefix)

    n_input = count_bam(bam_file,num_threads)
    stats = get_flagstat(bam_file,n_check = n_input)
    print(f'\n{bam_file} flagstat:')
    print("\n".join([f'{k}: {v}' for k,v in stats.items()]))
    print("\n")

    # BAM cleaning 
    # Check if the .bam contains only primary alignments - this is important 

    #status = checkbam(bam_file)
    status = False
    if stats['duplicate_reads']>0:
        print(f'Found duplicate reads in {bam_file}!')
        status = True
    if stats['unmapped_reads']>0:
        print(f'Found unmapped reads in {bam_file}!')
        status = True
    if stats["secondary_alignments"]>0:
        print(f'Found secondary_alignments in {bam_file}!')
        status = True
    if stats["supplementary_alignments"]>0:
        print(f'Found supplementary_alignments in {bam_file}!')
        status = True

    if(status):
        clean_bam = '%s/clean.bam' % output_prefix
        filter_bam(bam_file,clean_bam,num_threads)
        print(f'Created: {clean_bam}')
        bam_file = clean_bam

        # report the stats for 
        stats = get_flagstat(bam_file,n_check = n_input)
        print(f'\n{bam_file} flagstat:')
        print("\n".join([f'{k}: {v}' for k,v in stats.items()]))
        print('\n')

    else:
        print(f'{bam_file} contains only primary alignments.')


    #############################################################
    # Count tags
    #############################################################
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
