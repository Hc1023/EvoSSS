# =================================================
# virodecode
# Version: 1.0 using Singularity
# Singularity installed
# 
# Author: Sisi Huang
# Made at: Jiang Lab
# Github: https://github.com/Hc1023/ViroDecode
# Docker: https://hub.docker.com/r/sisih/virodecode
# 
# =================================================

__author__ = 'sisihuang'
import sys, os
import argparse

## Set up arguments

parser = argparse.ArgumentParser()

# python3 run_virodecode_singularity.py -i virodecode_1.sif -1 fq1 -2 fq2
# Required Arguments
parser.add_argument("-i", "--image", type=str, required=True,
                    help="Image path")
parser.add_argument("-1", "--read1", type=str, required=True, 
                    help="Provide read1 file path")
parser.add_argument("-2", "--read2", type=str, required=True,
                    help="Provide read2 file path")


# NOT Required Arguments
parser.add_argument("-o", "--outDir", type=str, default=".", required=False,
                    help="Provide path to a directory to place the output")
parser.add_argument("-r", "--ref", type=str, default="", required=False,
                    help="Reference genome")
parser.add_argument("-t", "--thread", type=str, default="8", required=False,
                    help="Thread")
parser.add_argument("-a", "--align", type=str, default="1", required=False,
                    help="Alignment tool  1: BOWTIE2 (default) 2: BWA MEM")
parser.add_argument("-c", "--coverage", type=str, default="5000", required=False,
                    help="Subsample reads if the mean coverage is higher than xxx")

parser.add_argument("-g", "--gff", type=str, default="0", required=False,
                    help="Gff file, for Icarus visualization")
parser.add_argument("-G", "--genome_version", type=str, default="0", required=False,
                    help="Genome version (e.g. NC_045512.2) for SNP annotation by snpEff")
parser.add_argument("-p", "--amplicon", type=str, default="0", required=False,
                    help="Set to 1 if amplicon sequencing")
parser.add_argument("-b", "--bedpe", type=str, default="0", required=False,
                    help="BEDPE file of primer pair locations")
parser.add_argument("-A", "--assembly", type=str, default="0", required=False,
                    help="Set to 1 if needing assembly results (fasta)")
parser.add_argument("--min_reads2", type=str, default="10", required=False,
                    help="VarScan parameter: Was supported by at least the minimum number of supporting reads")
parser.add_argument("--min_var_freq", type=str, default="0.05", required=False,
                    help="VarScan parameter: Meets the minimum allele frequency threshold")
parser.add_argument("--strand_filter", type=str, default="0", required=False,
                    help="VarScan parameter: Passes a basic strand-bias filter if set to 1")
parser.add_argument("--clean", type=str, default="1", required=False,
                    help="Set to 0 if needing complete processing files")

args = parser.parse_args()


# requirement: read1 and read2 in the same directory
image = os.path.abspath(args.image)
InDir = os.path.dirname(os.path.abspath(args.read1))
OutDir = os.path.abspath(args.outDir)
r1Dir = os.path.abspath(args.read1)
r2Dir = os.path.abspath(args.read2)

fq1 = os.path.basename(args.read1)
fq2 = os.path.basename(args.read2)
ref = args.ref
thread=args.thread
align=args.align
coverage=args.coverage
gff=args.gff
genome_version=args.genome_version
amplicon=args.amplicon
bedpe=args.bedpe
assembly=args.assembly
min_reads2=args.min_reads2
min_var_freq=args.min_var_freq
strand_filter=args.strand_filter
clean=args.clean

dir_bind = f"{r1Dir}:{r1Dir},{r2Dir}:{r2Dir},{OutDir}:{OutDir}"
if bedpe != "0":
    dir_bind += f",{os.path.abspath(bedpe)}:{os.path.abspath(bedpe)}"


if args.ref == "":
    ref = "/ref/NC_045512.2.fna"
    gff = "/ref/sequence.gff3"
    genome_version = "NC_045512.2"
    
    os.system(
        f"singularity exec \
    -B {dir_bind} {image} bash /workfolder/main.sh \
    {InDir} {OutDir} {fq1} {fq2} {ref} {thread} \
    {align} {coverage} {gff} {genome_version} {amplicon} {bedpe} {assembly} \
    {min_reads2} {min_var_freq} {strand_filter} {clean}"
    )

else:

    os.system(
        f"singularity exec \
    -B {dir_bind},{ref}:{ref} {image} bash /workfolder/main.sh \
    {InDir} {OutDir} {fq1} {fq2} {ref} {thread} \
    {align} {coverage} {gff} {genome_version} {amplicon} {bedpe} {assembly} \
    {min_reads2} {min_var_freq} {strand_filter} {clean}"
    )
    