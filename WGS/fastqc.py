#!/usr/bin/env python3

import os
import re
import click
import shlex
import zipfile
import subprocess
from collections import namedtuple

import pandas as pd


class FastQC:
    """
    Runs FastQC on a FASTQ file and unzips.

    Example
    -------
    >>> fastqc = FastQC(3)
    >>> fastqc('/path/to/sample.fastq.gz', '/output/directory')
    """
    dependencies = {'exe':['perl', 'java', 'fastqc', 'unzip']}

    def __init__(self, threads = 2):
        self.threads = threads

    def __call__(self, fastq, outdir):
        zipped_fastqc, data_file = self._get_fastqc_paths(fastq, outdir)

        cmd = f"fastqc {fastq} -q -o {outdir} -f fastq -t {self.threads}"

        print(cmd)
        process = subprocess.run(shlex.split(cmd), stderr=subprocess.PIPE)
        if process.returncode!=0:
            raise SystemError(process.stderr.decode("utf-8"))

        with zipfile.ZipFile(zipped_fastqc, 'r') as zip_ref:
            zip_ref.extractall(outdir)
        os.remove(zipped_fastqc)

        return data_file

    def __repr__(self):
        return f"FastQC({self.threads})"

    @staticmethod
    def _get_fastqc_paths(fastq, outdir):
        """
        Get the prospective zip and data file path when you run FastQC
        >>> fastq = '/example_fastq/D_19_0001.fastq.gz'
        >>> outdir = '/example_outdir'
        >>> FastQC._get_fastqc_paths(fastq, outdir)
        ('/example_outdir/D_19_0001_fastqc.zip', '/example_outdir/D_19_0001_fastqc/fastqc_data.txt')
        """
        file_name = re.sub(r'\.fastq.*', '', os.path.basename(fastq))
        zip_file = f'{file_name}_fastqc.zip'
        unzipped_dir = f'{file_name}_fastqc'

        zip_file_abspath = os.path.join(outdir, zip_file)

        fastqc_data_abspath = os.path.join(outdir, unzipped_dir, 'fastqc_data.txt')

        return zip_file_abspath, fastqc_data_abspath

class FastQCParser:
    """
    A FastQC Parser.
    Parses fastqc_data.txt file to organize and calculates some convenient data
    from the input text file.

    Parameters
    ----------

    fastqc_data : str
        Path to fastqc_data.txt file.

    Attributes
    ----------

    basic_stats : pd.DataFrame
        Displays a tablse with some basic stats

    q30 : float
        Q30 value for the input fastqc_data

    gc : int
        GC-content value for the input fastqc_data

    total_read : int

    overall_quality : str
        "PASS" if Q30 is above 80. Else "WARNING"

    bsq : pd.DataFrame
        Per base sequence quality

    bsc : pd.DataFrame
        Per base sequence content

    sqc : pd.DataFrame
        Per sequence quality scores

    sgc : pd.DataFrame
        Per sequence GC content


    Example
    -------
    >>> fastqc_data = '/path/to/fastqc_data.txt'
    >>> parsed_fastqc = FastQCParser(fastqc_data)
    >>> parsed_fastqc.basic_stats # See basic stats
    >>> parsed_fastqc.q30 # See q30

    """

    def __init__(self, fastqc_data: str):


        self.fastqc_data = fastqc_data

        with open(self.fastqc_data, 'r') as f:
            file_content = f.readlines()

        Section = namedtuple('Section', 'line_number description')

        self.separator_list = list()
        for indx, line in enumerate(file_content):
            if bool(re.match(r'^\>\>', line)):
                self.separator_list.append(Section(indx, line))
        self.last_line_num = self.separator_list[-1].line_number

        self.basic_stats = self._make_data_section(regex = r'^\>\>Basic Statistics')
        self.bsq = self._make_data_section(regex = r'^\>\>Per base sequence quality')
        self.bsc = self._make_data_section(regex = r'^\>\>Per base sequence content')
        self.sqc = self._make_data_section(regex = r'^\>\>Per sequence quality scores')
        self.sgc = self._make_data_section(regex = r'^\>\>Per sequence GC content')

        q30_seq_count = sum(self.sqc.loc[self.sqc['#Quality']>=30, 'Count'])
        total_seq_count = self.basic_stats.loc[self.basic_stats['#Measure']=='Total Sequences', 'Value']

        self.q30 = round((q30_seq_count / float(total_seq_count)) * 100, 2)

        self.gc = int(self.basic_stats.loc[self.basic_stats['#Measure']=='%GC', 'Value'])

    def __repr__(self):
        return f"FastQCParser('{self.fastqc_data}')"

    @property
    def overall_quality(self):
        if not isinstance(self.q30, float):
            return None
        else:
            assert self.q30<100, 'Q30 was calculated as above 100%'
            if self.q30>80:
                return 'PASS'
            else:
                return 'WARNING'

    @property
    def total_read(self):
        if not isinstance(self.basic_stats, pd.DataFrame):
            return None
        total_seq_count = int(self.basic_stats.loc[self.basic_stats['#Measure']=='Total Sequences', 'Value'])
        return total_seq_count

    def _distinguish_data_sections(self, fastqc_data):
        with open(fastqc_data, 'r') as f:
            file_content = f.readlines()

        Section = namedtuple('Section', 'line_number description')

        self.separator_list = list()
        for indx, line in enumerate(file_content):
            if bool(re.match(r'^\>\>', line)):
                self.separator_list.append(Section(indx, line))

    def _get_data_section(self, regex):
        for i in range(0, len(self.separator_list), 2):
            if bool(re.match(regex, self.separator_list[i].description)):
                skiprows = self.separator_list[i].line_number + 1
                skipfooters = self.separator_list[i+1].line_number - 1
                break

        return skiprows, skipfooters

    def _make_data_section(self, regex):
        skiprows, skipfooters = self._get_data_section(regex)

        data = pd.read_csv(self.fastqc_data, skiprows=skiprows, skipfooter=self.last_line_num - skipfooters,
                        header = 0, sep = '\t', engine='python')
        # engine = 'python' because c engine does not support the skipfooter option
        return data


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
