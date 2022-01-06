#!/usr/bin/env python3

import click
import pandas as pd

from 


class CalculateGeneDepth:

    def __init__(self, genelist):
        genelist = pd.read_csv(genelist, sep="\t")[["Symbol","Chrom","Start","End"]]

    def __call__(self, mpileup):
        mpileup = MpileupParser(mpileup)
        depth_table = self.make_depth_table(mpileup)
        return depth_table

    def make_depth_table(self, mpileup: pd.DataFrame) -> pd.DataFrame:
        """
        Parameter
        ---------
        mpileup : pandas.DataFrame
            Index:
                RnageIndex
            Columns:
                Name: chrom
                Name: pos
                Name: ref
                Name: depth
                Name: seq
                Name: qual
            The parsed pileup format of samtools-mpileup.

        Returns
        -------
        pandas.DataFrame
            Index:
                RnageIndex
            Columns:
                Name: HUGO_SYMBOL, dtype: str
                Name: DEPTH, dtype: float64
                Name: VARIANT_GROUP_ID, dtype: str
            The depth of echo gene.
        """

        data = list()
        for gene, chrom, start, end  in self.gene_df[["Symbol","Chrom","Start","End"]].itertuples(index=False):
            depth = mpileup.calc_region_depth(chrom, start, end)
            data.append([gene, depth, 'GENE'])
        return pd.DataFrame(data, columns = ['HUGO_SYMBOL', 'DEPTH', 'VARIANT_GROUP_ID'])
