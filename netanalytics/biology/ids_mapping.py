import os 

import numpy as np
import pandas as pd

from pathlib import Path

def check_type(t):
    t = str(t).lower()
    admitted_types = ['entrez_id', 'hgnc_id','hugo_symbol', 'ensembl_gene_id', 'vega_id', 'ucsc_id ',
                'uniprot_ids', 'pubmed_id']
    if not t in admitted_types:
        raise ValueError("You have to insert one between the options for the ids.")
    if t == 'hugo_symbol':
        return 'symbol'
    return t


def get_genes_id(genes, from_='entrez_id', to_='hugo_symbol',
                 return_indices=False):
    """
    genes: list
        List of genes identifier.

    from: string, optional default='entrez_id'
        The type of identifier in the list of genes specified. Note that if the
        type is not the right one the function simply returns NaN as new
        identifiers.
        Options:
            - hgnc_id
            - hugo_symbol
            - ensembl_gene_id
            - vega_id
            - ucsc_id
            - uniprot_ids
            - pubmed_id
    to: string, optional  default ='hugo_symbol'
        The type of output identifier. If there is no match between the input
        and the stored database a NaN is put in its place.
        Options:
            - hgnc_id
            - hugo_symbol
            - ensembl_gene_id
            - vega_id
            - ucsc_id
            - uniprot_ids
            - pubmed_id
    return_indices: boolean, optional default=False
        If True the function returns the indices of the list genes for which
        it was not possible to find a match.
    Returns
    -------
    list:
        The list of new identifiers.
    list: optional
        The list of indices that are not matched with the new identifier.
    """

    from_ = check_type(from_)
    to_ = check_type(to_)
    path = os.path.realpath(__file__)
    path = path.split('/')[:-1]  # if not unix/OX should raise exception 
    path = '/'.join(path)
    gene_mapping = pd.read_table(path+"/databases/HUGO_proteing_encoding_genes.txt",
                                 sep='\t')
    gene_mapping = gene_mapping.set_index(from_)
    gene_mapping = pd.DataFrame(gene_mapping[to_])
    results = []
    not_matched = []
    for i, g in enumerate(genes):
        try:
            _id = gene_mapping.loc[g, to_]
        except KeyError:
            _id = None
            not_matched.append(i)
        results.append(_id)
    if return_indices:
        return results, not_matched
    return results
