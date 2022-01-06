
import pandas as pd





class CalcGeneDepth:
    genes = pd.read_csv(gene_table, 
                        sep='\t', 
                        index_col='SYMBOL',
                        names=['ENSEMBL_ID', 'SYMBOL', 'HGNC symbol', 'CHROM', 'START', 'END'], 
                        header=0)
    
    def __init__(self, version):
        if version != None:
            list_path, _ = select_gene_list_version(version)
            self.gene_df = pd.read_csv(list_path, sep="\t")[["Entrez_accession","Ensembl_ID"]]
        else:
            raise Exception(f'version value: {version}')
        
    def __call__(self, mpileup):
        mpileup = MpileupParser(mpileup)
        depth_table = self.make_depth_table(mpileup)
        depth_table['VARIANT_GROUP_ID'] = 'GENE'
        return depth_table
    
    @classmethod
    def from_configs(cls, configs):
        return cls(version = configs.version)
    
    def run_sample(self, sample):
        return self.__call__(sample.mpileup_result)
    
    def make_depth_table(self, mpileup):
        data = list()
        for gene, ensembl_id in self.gene_df.itertuples(index=False):
            chrom, start, end = self._get_genomic_position(gene, ensembl_id)
            depth = mpileup.calc_region_depth(chrom, start, end)
            data.append([gene, depth])
        return pd.DataFrame(data, columns = ['HUGO_SYMBOL', 'DEPTH'])
    
    def _get_genomic_position(self, symbol, ensembl_id):
        """
        Parameter
        ---------
        symbol : str
            The gene symbol to look up its genomic position.
        
        Returns
        -------
        tuple (str, int, int)
            The genomic position of the queried gene symbol.
        """
        try:
            chrom, start, end = self.genes.loc[self.genes['ENSEMBL_ID']==ensembl_id, ['CHROM', 'START', 'END']].loc[symbol,:]
        except KeyError as e:
            print(f'[INFO] {symbol} gene not in gene_table.')
            chrom, start, end = 'chr1', 0, 0
        return chrom, start, end
