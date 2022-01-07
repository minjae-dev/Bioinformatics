#!/usr/bin/env python3

import os
import re
import click
import shlex
import datetime
import subprocess
from statistics import mean

import numpy as np
import pandas as pd



def log(msg):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %T")
    print(f"{timestamp}: {msg}")


class MpileupPipeline:
    def __init__(self, configs):
        self.mpileup = Mpileup.from_configs(configs)

    def run_sample(self, sample, directory_name = "MPILEUP", overwrite=False):
        mpileup_file = self.__call__(
            outdir = sample.outdir
            , prefix = sample.prefix
            , recal_bam = sample.recal_bam
            , directory_name = directory_name
            , overwrite = overwrite)

        sample.mpileup_result = mpileup_file

    def __call__(self, outdir, prefix, recal_bam, directory_name = "MPILEUP", overwrite = False):
        output_directory = os.path.join(outdir, directory_name)
        print(f'Creating directory: {output_directory}')
        
        if (os.path.exists(output_directory)) and (overwrite is True):
            print(f'Overwriting {output_directory}..')
            rmtree(output_directory)
        elif (os.path.exists(output_directory)) and (overwrite is False):
            raise FileExistsError("Output directory already exists! Remove directory youself or use overwite parameter.")

        # Create output directory and output_prefix
        os.mkdir(output_directory)
        print(f'All outputs will be saved in {output_directory}')

        output_prefix = os.path.join(output_directory, prefix)

        mpileup_file = os.path.join(f'{output_prefix}.mpileup')
        self.mpileup(input_bam = recal_bam, output = mpileup_file)

        return mpileup_file


class Mpileup:
    """
    Samtool's Mpileup pipeline class.  

    Parameters
    ----------
    ref_hg : str
        Path to the reference fasta file (hg19.fa)
    target_bed : 
        Path to BED file for running mpileup on the regions specified in the BED file.
        
    Methods
    -------
    __call__(recal_bam = '/path/to/recal.sorted.bam', output = '/path/to/output.txt')
        Run mpileup with the `input_bam` file and output the result to `output`
    
    """
    __version__ = '1.0.0'

    dependencies = {'exe': ['samtools']}
    
    def __init__(self, ref_hg, target_bed):
        self.ref_hg = ref_hg
        self.target_bed = target_bed
    
    def __call__(self, input_bam, output):
        """
        Runs Samtool's mpileup.
        Takes about ~20 minutes to run.
        
        Parameters
        ----------
        input_bam : str
            Path to the input.bam to run mpileup on.
        output : str
            Absolute path of the output mpileup file to be created.
        """
        cmd = ('samtools mpileup -a ' # -a option for make output including zero coverage region
               f'--fasta-ref {self.ref_hg} '
               f'{input_bam} '
#                '--adjust-MQ 50 '
                '--min-MQ 1 '
               f'--positions {self.target_bed} '
               f'--output {output}')
        log('Run Pileup')

        process = subprocess.run(shlex.split(cmd),
                                stderr=subprocess.PIPE)
        if process.returncode!=0:
            print(cmd)
            raise SystemError(process.stderr.decode("utf-8"))
        log('Mpileup done...')


class MpileupParser:
    """
    A convenient mpileup file parser. The instance of the class is used to 
    conveniently extract some useful statistics about the mpileup file such as 
    average depth, region depth, and uniformity.

    Parameters
    ----------
    path : str
        Path to the mpileup file.        
    
    Attributes
    ----------
    uniformity : float
        Calculated uniformity value from the mpileup file
        
    avg_depth : float
        Average depth of the entire mpileup file
        
    Methods
    -------
    calc_region_depth(chrom='chrX', 3000, 4000)
        Calculates the average depth of the specified region
        
    Example
    -------
    >>> parsed_mpileup = MpileupParser('/path/to/sample.mpileup.txt')
    >>> parsed_mpileup.uniformity
    """

    __version__ = '1.0.0'

    def __init__(self, path):
        self.path = path

        self.mpileup = pd.read_csv(self.path, header=None, sep='\t', 
             names=['chrom', 'pos', 'ref', 'depth', 'seq', 'qual'])
        self.mpileup['ref'] = self.mpileup['ref'].str.upper()
        self.mpileup['seq'] = [clean_seq(seq) for seq in self.mpileup['seq']]
        self.mpileup = self.mpileup[self.mpileup['qual'].notnull()].copy()
        processed_seq = list()
        processed_qual = list()
        for seq, qual in self.mpileup[['seq', 'qual']].itertuples(index=False):
            cleaned_seq = clean_seq(seq)
            if len(cleaned_seq) != len(qual):
                print (cleaned_seq, qual)
            new_seq, new_qual = remove_asterisk(cleaned_seq, qual)
            processed_seq.append(new_seq)
            processed_qual.append(new_qual)
        self.mpileup['seq'] = processed_seq
        self.mpileup['qual'] = processed_qual
            
    def __str__(self):
        print(self.mpileup)
        return ''
    
    def __repr__(self):
        return f'MpileupParser("{self.path}")'
    
    @property
    def uniformity(self):

        no_mitochondria_depth = self.mpileup.loc[
                                    self.mpileup['chrom']!='chrM'
                                    , 'depth'
                                    ]

        mean_depth = mean(no_mitochondria_depth) / 2
        ratio = [0 if d < mean_depth else 1 for d in self.mpileup['depth']]
        return mean(ratio)
    
    @property
    def avg_depth(self):
        return round(mean(self.mpileup['depth']), 2)
    
    def calc_region_depth(self, chrom, start, end):
        """
        Calculates the average depth of the specified region.

        Parameters
        ----------    
        chrom : str
            Chromosome number, X, or Y 
        start : int
            Starting base pair position 
        end : int
            Ending base pair position 

        Returns
        -------
        float 
            Average depth of the specified region.

        Example
        -------
        >>> parsed_mpileup.calc_region_depth('chrX', 3000, 4000)
        >>> parsed_mpileup.calc_region_depth('chr22', 3000, 4000)
        
        """
        cond1 = self.mpileup['chrom']==chrom
        cond2 = start<=self.mpileup['pos']
        cond3 = self.mpileup['pos']<=end

        depths = self.mpileup.loc[(cond1) & (cond2) & (cond3), 'depth']

        if len(depths)==0:
            return 0.0
        else:
            return round(mean(depths), 2)


def clean_seq(seq):
    out_seq = _remove_indel(seq)
    out_seq = _remove_low_quality(out_seq)
    return out_seq


def _remove_indel(seq):
    found = re.search(r'[+-][0-9]+', seq)
    while bool(found) is True:
        # Indel 시작과 끝 지점을 찾아서
        start = found.start()
        end = (found.end() + abs(int(found.group())))
        indel_string = seq[start:end]

        # 삭제
        seq = seq.replace(indel_string, '')
        # Indel이 또 있는지 확인
        found = re.search(r'[+-][0-9]+', seq)

    return seq


def _remove_low_quality(seq):
    # https://en.wikipedia.org/wiki/Pileup_format
    # http://samtools.sourceforge.net/pileup.shtml
    unwanted_symbols = r'\^.'
    unwanted_symbols_2 = r'[0-9\@\'\&\!\?\_\$\:\;F\"\+\#\-\=\%\)\(\/\~\}\[]+' # Not entirely sure, but these reads might be low mapping quality reads so pileup was unsure what allele this read contained  
    # TODO: Split up the regexes and find out the meaning of all these regexes

    out_seq = re.sub(unwanted_symbols, '', seq)
    out_seq = re.sub(unwanted_symbols_2, '', out_seq)

    return out_seq


def remove_asterisk(seq, qual):
    ''' 
    Removes the deletion placeholder symbol (*) in the mpileup output. 
    This removal needs to be dealt separately because the asterisk(*) 
    unlike other symbols have a corresponding base-quality score 
    See: 
    1. http://seqanswers.com/forums/showthread.php?t=3388
    2. http://seqanswers.com/forums/showthread.php?t=5431
    '''
    assert len(seq)==len(qual), f"Unequal seq and qual string length {seq}, {qual}"
    
    seq_list = list(seq)
    qual_list = list(qual)

    while '*' in seq_list:
        pop_idx = seq_list.index('*')
        seq_list.pop(pop_idx)
        qual_list.pop(pop_idx)
    out_seq = ''.join(seq_list)
    out_qual = ''.join(qual_list)
    
    assert len(out_seq)==len(out_qual), f"Unequal out seq and out qual string length {out_seq}, {out_qual}"

    return out_seq, out_qual
