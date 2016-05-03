work_dir = 'dir_with_input_files'
import pandas as pd
import sys
sys.path.insert(0, work_dir)
import neighboring_genes_functions as ngf

#==============================================================================
# Run neighboring gene functions for each model organism
#==============================================================================

input_dir = work_dir + 'input/'
output_dir = work_dir + 'output/neighbors/'
# Specify number of neighbors to pair each gene with
no_pairs=100
# Specify number of randomly selected genes to pair each gene with, to set
# correlation baseline
no_pairs_random=100

# Now, for each organism, specify:
# location of expression datasets, location of list of genes, name of the gene
# ID column, and organism name.

#==== S cerevisiae ============================================================
sc = [input_dir + 'yeast/expTPM_sc.txt',
    input_dir + 'yeast/ensgene_sc.txt', 'gname', 'sc']


#==== C elegans ===============================================================
ce = [input_dir + 'worm/RNAseqHillier_expTPM_ce.txt',
    input_dir + 'worm/ensgeneWs220_ce.txt', 'wbid', 'ce']

#==== Drosophila ==============================================================
dm = [ input_dir + 'fly/expTPM_dm.txt',
    input_dir + 'fly/ensgene_dm.txt', 'FBgn', 'dm']

#==== Mouse ===================================================================
mm = [input_dir + 'mouse/expTPM_mm.txt',
    input_dir + 'mouse/ccds_genes_ensembl_mm.txt', 'ENSG', 'mm']

#==== Human ===================================================================
hs = [input_dir + 'human/expTPM_hs.txt',
    input_dir + 'human/ccds_genes_ensembl_hs.txt', 'gname', 'hs']

# Specify the lists of organism data to work on
organisms = [sc, ce, dm, mm, hs]

# Iterate through each list of data
for (exp_dir, gene_list, ID, org) in organisms:

    # pair genes with neighbors and random genes, compute their correlation in
    # expression
    neighbors = ngf.pair_corr(gene_list,ID=ID,exp_table=exp_dir,no_pairs=no_pairs)
    random = ngf.pair_random_corr(gene_list, ID=ID, exp_table=exp_dir , no_pairs=no_pairs_random)

    # delete operons, only necessary for c. elegans
    if org == 'ce': 
        neighbors = ngf.del_operons(input_dir + 'worm/operons.txt', neighbors)
    
    # create location names to save files
    neighbors_dir = output_dir + org + '_' + str(no_pairs) + 'neighbors_TPM_spearmanr_Hillier.txt'
    random_dir = output_dir + org + str(no_pairs_random) + 'random_TPM_spearmanr_Hillier.txt'

    # save correlation datasets
    neighbors.to_csv(neighbors_dir, index=False)
    random.to_csv(random_dir, index=False)

    # compute sliding medians using default window size of 1000 gene pairs
    sliding_median = ngf.sliding_median_iqr(neighbors, random=random, compute_random=1)
    sliding_median.to_csv(output_dir + org + 'slidingMedian_IQR_randomCI.txt', index=False)
