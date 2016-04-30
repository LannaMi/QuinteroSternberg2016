from __future__ import division
import pandas as pd
import numpy as np
import itertools
import scikits.bootstrap as bs
from scipy.optimize import curve_fit

#==============================================================================
# Functions for data processing of gene neighbors
#==============================================================================
def find_neighbors(df, ID='wbid', no_pairs=100):
    """
    Identify and classify by configuration pairs of neighboring genes in
    a list 

    Pairs numbers first and indexes dataframe

    Parameters
    ----------
    df: gene dataFrame with chromosome, start, end, orientation, gene_ID
    ID: name of gene_ID column (according to organism, e.g. wbid, FBgn, ENSG)
    exp_table: path to dataframe with expression data, or pandas df object
    no_pairs: number of neighbors to pair each gene with
    FPKM_to_TPM: normalize each RNAseq exp to TPM as: FPKM_gene/FPKM_sum*1e6

    Returns
    -------
    neighbors : dataframe with neighbors with all initial information plus
    start, end, size, config of intergenic region

    Config is (++,--:parallel,+-:convergent,-+:divergent)
    Start, end and size belong to genomic region from head of the first gene to
    head of the second gene, i.e. size will be shortest for divergent genes,
    longest for convergent genes.

    """

    #load dataframe if not provided yet
    if isinstance(df, basestring): 
        df = pd.read_csv(df)
    #make sure positions are integers and not strings
    df = df.convert_objects(convert_numeric=True)
    #sort by position
    df.sort('start', inplace=True)
    #group by chromosome, returns a list of tuples where the dataframe is the
    #second element (index 1)
    df_chromosomes = df.groupby('chromosome')
    #create new dataframe with labeled columns 
    cols_gene1 = [col + '1' for col in df.columns.tolist()]
    cols_gene2 = [col + '2' for col in df.columns.tolist()]
    neighbors = pd.DataFrame(columns = cols_gene1 + cols_gene2) 
    def pair_genes(chr_genes):
        # function to match gene pairs
        chrom_size = len(chr_genes)
        chrom_size_no_pairs = chrom_size - no_pairs
        chr_index = np.arange(chrom_size_no_pairs)
        pairs = np.arange(no_pairs)
        # generate an array with each element repeated n=no_pairs times to get indices
        # of first gene in neighbor pair
        g1_index = np.repeat(chr_index, no_pairs)
        # generate a list of arrays with length=no_pairs
        # all elements in array increase by 1, 2, 3..no_genes
        # so that each gene1 is paired with no_pairs neighbors
        # this unavoidably leaves last 100 genes without pairing...but sample is
        # probably good enough already
        gene_iter = np.nditer(np.arange(1,chrom_size_no_pairs+1))
        g2_index = np.concatenate([pairs + x for x in gene_iter])
        gene1 = chr_genes.iloc[g1_index]
        gene2 = chr_genes.iloc[g2_index]
        # reindex is necessary to join dataframes
        gene1.index = np.arange(len(gene1))
        gene2.index = np.arange(len(gene1))
        # put neighbors together and add suffixes to columns
        neighbors_chr = gene1.join(gene2, lsuffix='1', rsuffix='2')

        return neighbors_chr

    for chr_name, chr_genes in df_chromosomes: 
        print 'pairing ' + chr_name
        if len(chr_genes)<=no_pairs: 
            print 'too few genes, skipped..'
            continue
        paired_chr = pair_genes(chr_genes)
        neighbors = neighbors.append(paired_chr)

    print 'almost done...'
    neighbors['config'] = neighbors['orientation1'] + neighbors['orientation2']
    neighbors[ID] = neighbors[ID+'1'] + neighbors[ID+'2']
    # get start and end of intergenic region, which depends on gene orientation
    # parallel genes (++) This already have correct start and end
    parallelplus = pd.DataFrame()
    parallelplus[ID] = neighbors[neighbors['config']=='++'][ID]
    parallelplus['start'] = neighbors[neighbors['config']=='++'].start1
    parallelplus['end'] = neighbors[neighbors['config']=='++'].start2

    # parallel genes (--)
    parallelminus = pd.DataFrame()
    parallelminus[ID] = neighbors[neighbors['config']=='--'][ID]
    parallelminus['start'] = neighbors[neighbors['config']=='--'].end1
    parallelminus['end'] = neighbors[neighbors['config']=='--'].end2

    # divergent genes (-+)
    divergent = pd.DataFrame()
    divergent[ID] = neighbors[neighbors['config']=='-+'][ID]
    divergent['start'] = neighbors[neighbors['config']=='-+'].end1
    divergent['end'] = neighbors[neighbors['config']=='-+'].start2

    # convergent genes (+-)
    convergent = pd.DataFrame()
    convergent[ID] = neighbors[neighbors['config']=='+-'][ID]
    convergent['start'] = neighbors[neighbors['config']=='+-'].start1
    convergent['end'] = neighbors[neighbors['config']=='+-'].end2

    # stack correct start and end 
    position = pd.concat([parallelplus, parallelminus, divergent, convergent])

    # some genes will overlap, e.g. in a divergent pair, gene1 ends after gene
    # 2 starts; in a parallel minus pair one gene contains the other.
    start_abs = position[['start','end']].min(axis=1)
    end_abs = position[['start','end']].max(axis=1)
    position['start'] = start_abs
    position['end'] = end_abs
    # add to dataframe
    neighbors = pd.merge(neighbors, position, on=ID)
    # now size can be computed
    neighbors['size'] = neighbors.end - neighbors.start
    return neighbors

def add_exp_pairs(df, exp_table, ID, FPKM_to_TPM): 
    """
    get multiple expression indices for left and right genes of pairs

    Parameters
    -------

    df: dataframe with paired genes, or path to it 
    exp_table: dataframe with expression data in FPKM or TPM, or path to it
    ID: name of gene_ID column (according to organism, e.g. wbid, FBgn, ENSG)
    FPKM_to_TPM: normalize each RNAseq exp to TPM as: FPKM_gene/FPKM_sum*1e6

    Returns
    -------
    gene_1, gene_2 : Dataframes with gene_ID, exp data for l and r gene in pair


    """
    #load dataframe if not provided yet
    if isinstance(df, basestring): 
        df = pd.read_csv(df)
    #load table of gene expression and WBID
    if isinstance(exp_table, basestring): 
        exp_table = pd.read_csv(exp_table)
    exp_table = exp_table.convert_objects(convert_numeric=True)
    exp_numeric = exp_table._get_numeric_data()
    numeric_cols = exp_numeric.columns

    # convert to TPM from FPKM if requested
    if FPKM_to_TPM:
        exp_table[numeric_cols] = exp_numeric/exp_numeric.sum(axis=0)*1e6

    # delete genes that are not detected in any experiment
    # or detected in less than 80% of experiments
    never_expressed = exp_table[(exp_numeric.sum(axis=1)==0)|
        (exp_numeric.count(axis=1)<exp_numeric.shape[1]*0.8)][ID]
    exp_table = exp_table[~(exp_table[ID].isin(never_expressed))]

    # transformed expression values to ranked values for spearman correlation
    exp_table[numeric_cols] = exp_table._get_numeric_data().rank(axis=0)
    gene_1 = pd.DataFrame()
    gene_2 = pd.DataFrame()
    gene_1[ID] = df[ID+'1']
    gene_2[ID] = df[ID+'2']
    gene_1 = pd.merge(gene_1, exp_table, how='left', on=ID)
    gene_2 = pd.merge(gene_2, exp_table, how='left', on=ID)

    return gene_1, gene_2

def compute_corr(self, other):
    """
    Compute pairwise spearman correlation between rows of
    two DataFrame objects.

    Parameters
    ----------
    other : DataFrame

    Returns
    -------
    correls : Series
    """
    # compute pearson correlation 
    # This would be spearmanr if values were ranked as using add_exp_pairs
    this = self._get_numeric_data()
    other = other._get_numeric_data()

    return this.corrwith(other, axis=1) 

def pair_corr(df, ID, exp_table, no_pairs=100, FPKM_to_TPM=False):
    """
    Takes a list of genes, pairs them with nearest neighbors and computes
    correlation in gene expression
    Parameters
    ----------
    df: gene dataFrame with chromosome, start, end, orientation, gene_ID
    ID: name of gene_ID column (according to organism, e.g. wbid, FBgn, ENSG)
    exp_table: path to dataframe with expression data, or pandas df object
    no_pairs: number of neighbors to pair each gene with
    FPKM_to_TPM: normalize each RNAseq exp to TPM as: FPKM_gene/FPKM_sum*1e6

    Returns
    -------
    neighbors : Dataframe with paired genes and spearman correlation for each

    """
    print 'pairing neighbors...'
    neighbors = find_neighbors(df, ID=ID, no_pairs=no_pairs)
    print 'adding expression data...'
    g1_exp, g2_exp = add_exp_pairs(neighbors, exp_table=exp_table, ID=ID, 
            FPKM_to_TPM=FPKM_to_TPM)
    print 'computing correlation...'
    spearmanr = compute_corr(g1_exp, g2_exp)
    neighbors['spearmanr'] = spearmanr
    print 'done'

    return neighbors

def pair_random_corr(df, ID, exp_table, no_pairs=10, FPKM_to_TPM=False):
    """
    Takes a list of genes, pairs them randomly where chr_gene1 != chr_gene2 and
    computes correlation in gene expression

    Parameters
    ----------
    df: gene dataFrame with chromosome, start, end, orientation, gene_ID
    ID: name of gene_ID column (according to organism, e.g. wbid, FBgn, ENSG)
    exp_table: path to dataframe with expression data, or pandas df object
    no_pairs: number of neighbors to pair each gene with
    FPKM_to_TPM: normalize each RNAseq exp to TPM as: FPKM_gene/FPKM_sum*1e6

    Returns
    -------
    neighbors : Dataframe with randomly paired genes and spearman correlation for each

    """

    #load dataframe if not provided yet
    if isinstance(df, basestring): 
        df = pd.read_csv(df)
    no_genes = len(df)
    print 'pairing genes randomly...'
    # create random list of gene pairs of length no_pairs*no_genes and index
    length = no_genes*no_pairs
    random_index = np.random.randint(0,no_genes,size=length)
    random_gene1 = df.iloc[random_index]
    random_index = np.random.randint(0,no_genes,size=length)
    random_gene2 = df.iloc[random_index]
    # reindex is necessary to join dataframes
    random_gene1.index = range(len(random_gene1))
    random_gene2.index = range(len(random_gene1))
    random_neighbors = random_gene1.join(random_gene2, lsuffix='1', rsuffix='2')
    # delete those paired with themselves
    random_neighbors = random_neighbors[~(random_neighbors[ID+'1'] == random_neighbors[ID+'2'])]
    # delete those in the same chromosome
    random_neighbors = random_neighbors[~(random_neighbors['chromosome1'] == random_neighbors['chromosome2'])]

    print 'adding expression data...'
    g1_exp, g2_exp = add_exp_pairs(random_neighbors, exp_table=exp_table, ID=ID, 
            FPKM_to_TPM=FPKM_to_TPM)
    print 'computing correlation...'
    spearmanr = compute_corr(g1_exp, g2_exp)
    random_neighbors['spearmanr'] = spearmanr
    print 'done'

    return random_neighbors

def del_operons(operons_dir, neighbors):
    """
    Removes gene pairs that are in the same operon from neighbors dataframe
    Parameters
    ----------
    operons_dir:  path to list of operons
    neighbors: neighboring gene pairs dataframe

    Returns
    -------
    neighbors_woperons : neighbors stripped of gene pairs in the same operon

    """
    # delete genes in operons
    operons = []
    for line in open(operons_dir):
        # make list of pairs of genes in operons. Remove line break (\n)  and
        # whitespace
        for pair in itertools.permutations(tuple(line.strip('\n').split(',')), 2):
            operons.append(''.join(pair).replace(' ',''))
    # delete operons from dataframe
    neighbors['ggname'] = neighbors['name1'] + neighbors['name2']
    neighbors_woperons = neighbors[(~neighbors['ggname'].isin(operons))]

    return neighbors_woperons

def del_dup_genes(dup_genes, neighbors, ID):
    """
    Removes gene pairs of duplicated genes (sharing a recent common ancestor)
    Parameters
    ----------
    dup_genes:  list of tuples of duplicated genes
    neighbors: neighboring gene pairs dataframe
    ID: name of gene_ID column (according to organism, e.g. wbid, FBgn, ENSG)

    Returns
    -------
    neighbors_wodup : neighbors stripped of gene pairs of duplicated genes

    """
    if isinstance(neighbors, basestring): 
        neighbors = pd.read_csv(neighbors)
    dup_pairs = []
    
    for line in dup_genes:
        # make list of duplicate gene pairs 
        for pair in itertools.permutations(line, 2):
            dup_pairs.append(''.join(pair).replace(' ',''))
    # delete pairs of duplicate genes from dataframe
    neighbors[ID] = neighbors[ID+'1'] + neighbors[ID+'2']
    neighbors_wodup = neighbors[(~neighbors[ID].isin(dup_pairs))]

    return neighbors_wodup

def exp_decay(x, N_0, lmbd, c):
    """
    Exponential decay function
    ----------
    x: data
    N_0: initial value
    lmbd: lambda
    c: baseline

    Returns
    -------
    Exponential decay Function

    """

    return N_0 * np.exp(-lmbd * x) + c

def sliding_median_iqr(neighbors, random=None, compute_random=0, window=1000, p0=None):
    """
    Compute sliding median of spearmanr and size, interquartile range 
    and 95% CI of spearmanr of randomly paired genes

    Parameters
    ----------
    neighbors: neighboring gene pairs dataframe
    window: size of window for sliding median

    Returns
    -------
    rolling_median: sliding median of spearmanr and size with IQR for spearmanr
    median and 95% confidence interval of median from random pairs

    """
    #load dataframe if not provided yet
    if isinstance(neighbors , basestring): 
        neighbors  = pd.read_csv(neighbors)
    if compute_random and isinstance(random , basestring): 
        random  = pd.read_csv(random)

    # sort by size to do sliding window with increasing intergenic distance
    # nans cause error in sliding median
    neighbors  = neighbors.sort('size').dropna()

    print 'computing sliding median...'
    # compute rolling medians. 1000 looks good, less is unnecesserily heavy and noisy.
    rolling_median_spearmanr = pd.rolling_median(neighbors.spearmanr, window)

    print 'computing IQR...'
    # compute interquartile range (IQR). Top 75% and bottom 25%.
    rolling_spearmanr_q1 =  - pd.rolling_quantile(neighbors.spearmanr, window, 0.25) + \
            rolling_median_spearmanr 
    rolling_spearmanr_q3 = pd.rolling_quantile(neighbors.spearmanr, window, 0.75) - \
            rolling_median_spearmanr 
    rolling_median_size = pd.rolling_median(neighbors['size'], window)/1000

    # put it all together
    rolling_median_s = pd.DataFrame({'spearmanr': rolling_median_spearmanr, 
        'size':rolling_median_size, 'q1': rolling_spearmanr_q1, 'q3': rolling_spearmanr_q3})

    # drop all nans from sliding median (first 1000 because of window)
    rolling_median_s = rolling_median_s.dropna()

    # reindex is necessary
    rolling_median_s.index = np.arange(len(rolling_median_s))

    if compute_random:
        print 'computing random pairs median CI'
        # compute 95% confidence interval of median in random pairs
        ci_median = bs.ci(random.spearmanr.dropna().loc[:20000], np.median)
        rolling_median_s['random_lci'] = ci_median[0]
        rolling_median_s['random_hci'] = ci_median[1]

    print 'fitting to exp decay...'
    popt_s, pcov_s = curve_fit(exp_decay, rolling_median_s['size'], 
            rolling_median_s.spearmanr, p0=p0)

    rolling_median_s['popt1'] = popt_s[0]
    rolling_median_s['popt2'] = popt_s[1]
    rolling_median_s['popt3'] = popt_s[2]

    print 'done'
    return rolling_median_s

