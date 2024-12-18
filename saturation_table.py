#!/usr/bin/env python3
import argparse
import sys
import os
import subprocess

# Argument parsing
parser = argparse.ArgumentParser(description="Compute sequencing saturation curve")
parser.add_argument('-b', '--bam', required=True, help='Input BAM file')
parser.add_argument('-o', '--output', required=True, help='Output file')
parser.add_argument('-n', '--ncells', type=int, required=True, help='Number of cells')
parser.add_argument('-r', '--mapping_rate', type=float, required=True, help='Mapping rate')
parser.add_argument('-s', '--nsteps', type=int, required=False,default = 20, help='Number of steps')
parser.add_argument('-c', '--ncpu', type=int, required=False,default = 1, help='Number of CPUs to use')
parser.add_argument('-t', '--temp', required=False,default = '_tmp', help='Temporary directory')
parser.add_argument('-f', '--tabfile', required=False,default = 'tags.tab', help='File to store tags')
args = parser.parse_args()

BAM = args.bam
OUTPUT = args.output
NCELLS = args.ncells
MAPPING_RATE = args.mapping_rate
NSTEPS = args.nsteps
NCPU = args.ncpu
TEMP = args.temp
TABFILE = os.path.join(args.temp,'tags.tab')



CODE_DIR = os.path.dirname(os.path.abspath(__file__))
gettags_script = os.path.join(CODE_DIR,'scripts/gettags.py')
compute_saturation_script = os.path.join(CODE_DIR,'scripts/compute_saturation.r')

remove_temp = False
remove_tab = False
# Remove temporary directory if it exists
if os.path.exists(TEMP):
    print(f'Be careful! Found temp directory {TEMP}!')
    quit()
else:
    print(f'Creating temp directory {TEMP}')
    os.mkdir(TEMP)


if not os.path.exists(TABFILE):
    print(f"Alignment: {BAM}\nOutput: {OUTPUT}\nN cells: {NCELLS}\nMapping rate: {MAPPING_RATE}\nN threads: {NCPU}")
    gettags_command = ['python', gettags_script,  BAM, TEMP, '--num_threads', str(NCPU)]
    print(f"Executing command: {' '.join(gettags_command)}")
    subprocess.run(gettags_command)


    # Concatenate tag files into TABFILE
    with open(TABFILE, 'wb') as outfile:
        for filename in [x for x in os.listdir(TEMP) if x.endswith('.tsv')]:
            with open(os.path.join(TEMP, filename), 'rb') as infile:
                outfile.write(infile.read())
    print(f"Counting tags done: {TABFILE}")
else:
    print(f'Found tags {TABFILE}! Skipping counting tags')
# Remove temporary directory
if remove_temp:
    subprocess.run(['rm', '-rf', TEMP])
    print(f"Removed temp directory: {TEMP}")
# Execute compute_saturation.r script
subprocess.run(['Rscript', compute_saturation_script, '--ncpu', str(NCPU), '--mapping_rate', str(MAPPING_RATE), TABFILE, OUTPUT, str(NCELLS), '--nstep', str(NSTEPS)])


if remove_tab:
    subprocess.run(['rm', TABFILE])
print(f"Saturation done: {OUTPUT}")

