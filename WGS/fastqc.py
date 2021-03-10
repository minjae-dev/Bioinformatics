#!/usr/bin/env python3

import click
from cancerscan_fastq2bam.fastqc import FastQC

@click.command()
@click.option('-r1', 'R1', required=True, type=click.Path(exists=True, dir_okay=False), help='FASTQ Read1 file.')
@click.option('-r2', 'R2', required=True, type=click.Path(exists=True, dir_okay=False), help='FASTQ Read2 file.')
@click.option('-t', 'threads', type=click.INT, default=3, help='Number of threads to run.')
@click.option('-o', 'outdir', required=True, type=click.Path(exists=True, file_okay=False, writable=True), help='Output file directory path.')
def main(R1, R2, threads, outdir):
    fastqc = FastQC(threads)
    fastqc(R1, outdir)
    fastqc(R2, outdir)

if __name__=='__main__':
    main()
