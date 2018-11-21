import os
import math

import pandas as pd
import numpy as np

from scipy.stats import hypergeom
from collections import defaultdict


# Benjamini-Hochberg p-value correction for multiple hypothesis testing
def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def _perform_enrichment_bh(clusters, annotations, n_cores=1, verbose=0): #TODO implement parallel execution
    N = len(list(np.unique(annotations.index)))
    data = []
    enrichment = [[] for j in range(len(clusters))]
    cts = 0

    cls = []
    gos = []
    pvals = []

    NE_list = []

    for i, cluster in enumerate(clusters):
        genes_at_least_one = list(set(annotations.index).intersection(
                                  set(list(cluster.ravel()))))
        n = len(genes_at_least_one)
        ann_c = annotations.loc[genes_at_least_one]
        for a in np.unique(ann_c['annotations']):
            genes_annotated = annotations[annotations['annotations']==a]
            K = genes_annotated.shape[0]
            k = ann_c[ann_c['annotations']==a].shape[0]
            pval = hypergeom.sf(k-1, N, K, n)
            cls.append(i)
            gos.append(a)
            pvals.append(pval)
            if pval <= 0.05:
                cts += 1
                enrichment[i].append([a, pval])
        if verbose:
            print("Finished p-value computation cluster "+str(i))

        # Mean normalized entropy:
        d = float(np.unique(ann_c['annotations']).shape[0])  # number of different annotations in cluster c
        if d <= 1.:
            continue
        nec = 0.
        Cl = float(len(cluster))  # number of gene in cluster c
        sum_pi = 0.
        # counting map
        counting_map = dict.fromkeys(list(ann_c['annotations']), 0)
        for a in ann_c['annotations']:
            counting_map[a] += 1.
        counting_map = {k:v for k,v in counting_map.items() if v != 0}
        for a,v in counting_map.items():
            pi = v/Cl  # percentage of genes in c annotated with the considered go term
            sum_pi += pi * math.log(pi)
        nec = (-1. / (math.log(d))) * sum_pi
        NE_list.append(nec)
        if verbose:
            print("Finished mean normalized entropy cluster "+str(i))

    # applying BH correction
    BHpvals = p_adjust_bh(pvals)
    BHcts = 0
    BH_enrichment = [[] for j in range(len(clusters))]
    enriched_genes = []
    for i in np.where(BHpvals<0.05)[0]:
        BHcts += 1
        BH_enrichment[cls[i]].append([gos[i], BHpvals[i]])
    for i, c in enumerate(clusters):
        cluster_set = set()
        enriched_gos = np.array(BH_enrichment[i])
        if enriched_gos.size == 0:
            enriched_genes.append([])
            continue
        genes_at_least_one = list(set(annotations.index).intersection(
                                  set(list(c.ravel()))))
        ann_c = annotations.loc[genes_at_least_one]
        gos = enriched_gos[:,0]
        for g in gos:
           # print(ann_c[ann_c['annotations']==g].index.values)
            cluster_set = cluster_set.union(set(list(ann_c[ann_c['annotations']==g].index.values)))
        enriched_genes.append(list(cluster_set))
        if verbose:
            print("Finished BH correction cluster "+str(i))

    MNE = sum(NE_list) / float(len(NE_list)) if sum(NE_list) != 0 else 0
    enr_cluster = 0
    total_cluster = 0
    for i in range(len(clusters)):
        if len(clusters[i]) > 0:
            total_cluster += 1.
            if len(BH_enrichment[i]) > 0:
                enr_cluster += 1
    perc_enr_cluster = 100. * enr_cluster / total_cluster if total_cluster != 0 else 0

    perc_genes = sum([len(enriched_genes[i]) for i in range(len(clusters))])
    perc_genes = 100.*perc_genes/float(N) if N!=0 else 0
    return [BH_enrichment, MNE, perc_genes, perc_enr_cluster]


def _perform_enrichment_bonferroni(clusters, annotations, n_cores=1): #TODO implement parallel execution
    N = len(np.unique(annotations.index))
    bonferroni_correction = len(list(np.unique(annotations['annotations'])))*len(clusters)
    p_values = []
    genes_to_count = []
    enriched_clusters = 0
    for i, c in enumerate(clusters):
        enriched=False
        genes_at_least_one = list(set(annotations.index).intersection(
                                  set(list(c.ravel()))))
        n = len(genes_at_least_one)
        ann_c = annotations.loc[genes_at_least_one]
        for a in np.unique(ann_c['annotations']):
            genes_annotated = annotations[annotations['annotations']==a]
            K = genes_annotated.shape[0]
            k = ann_c[ann_c['annotations']==a].shape[0]
            pval = hypergeom.sf(k-1, N, K, n)
            if pval < (0.05/bonferroni_correction):
                genes_to_count += \
                        list(np.unique(ann_c[ann_c['annotations']==a].index))
                p_values.append((i, pval, a))
                enriched=True
        if enriched:
            enriched_clusters+=1
    return p_values, len(set(genes_to_count))/N, enriched_clusters/len(clusters)


def _perform_enrichment(clusters, annotations, correction, n_cores=1):
    if  correction.lower() == 'bonferroni':
        return _perform_enrichment_bonferroni(list_clusters, annotations,
                                           n_cores)
    else:
        return  _perform_enrichment_bh(list_clusters,annotations, n_cores)


def _get_database(genes, mode='go'):
    path = os.path.realpath(__file__)
    path = path.split('/')[:-1]  # if not unix/OX should raise exception
    path = '/'.join(path)

    if mode.lower()=='go':
        annotations = pd.read_table(path+"/databases/HSA_GO-BP.LST",
                                    index_col=0, header=None)
    elif mode.lower()=='kegg':
        annotations = pd.read_table(path+"/databases/HSA_Kegg_Pathways.lst",
                                    index_col=0, header=None)
    elif mode.lower()=='reactome':
        annotations = pd.read_table(path+"/databases/HSA_Reactome_Pathways.lst",
                                    sep='\t', skiprows=30, header=None)
    elif mode.lower()=='human_bp':
        annotations = pd.read_table(path+"/databases/Human_BP.lst", sep=';',
                                    header=None,index_col=0)
        annotations.index = [i.lower() for i in annotations.index]
    elif mode.lower()=='human_rr':
        annotations = pd.read_table(path+"/databases/Human_RR.lst", sep=';',
                                    header=None,index_col=0)
        annotations.index = [i.lower() for i in annotations.index]
    elif mode.lower()=='human_rp':
        annotations = pd.read_table(path+"/databases/Human_RP.lst", sep=';',
                                    header=None,index_col=0)
        annotations.index = [i.lower() for i in annotations.index]
    else:
         raise ValueError("The mode you specified is not implemented, please "
                          "try one among go, kegg, reactome, human_bp, "
                          "human_rr, human_rp")
    intersection = list(set(genes).intersection(set(annotations.index)))
    annotations = annotations.loc[intersection]
    to_keep = []
    for i, g in enumerate(annotations[1]):
        if annotations[annotations[1]==g].shape[0] != 1:
            to_keep.append(i)
    annotations = annotations.iloc[to_keep,:]
    to_drop=[]
    for g in annotations.index.values:
        if annotations.loc[g].shape[0] == 0:
            print(g)
            to_drop.append(g)
    annotations.drop(to_drop, inplace=True)
    annotations  = pd.DataFrame(annotations[1])
    annotations.columns = ['annotations']
    return annotations


class Enrichment()
    """
    It performs enrichment on a list of genes using different databases.

    Parameters
    ----------
    genes: list,  length=n
        List of HUGO symbol genes or entrez_id.
    mode: string or list, option default='go'
        The database to use to perform enrichment. Options are 'go', 'kegg',
        'reactome'. 'human_bp', 'human_rr', 'human_rp'.
        If list different databases are used to perform enrichment.
    correction: string, optional default='bh
        The correction to apply to the p-value to perform enrichment.
        Options are 'bh'=Benjamin-Hochberg correction or
        'bonferroni'=Bonferroni correction
    n_cores: int, optional default=1
        Number of cores to use to perform the computation.
    """

    def __init__(genes, mode='go', correction='bh', n_cores=1):
        self.genes = genes
        self.mode = mode
        self.correction = correction
        self.n_cores = n_cores

    def fit(labels):
    """
    Parameters
    ----------
    labels: list, length=n
        List of cluster labels for each gene.

    Returns
    -------
    self
    """

        if (len(self.genes) != len(labels)):
            raise ValueError("The length of genes list and cluster labels must "
                             "be the same, found %d and %d respectively."
                             %(len(self.genes), len(labels)))
        if type(self.genes[0]) == str:
            self.genes = [g.lower() for g in self.genes]
        else:
            self.genes = np.array(self.genes).astype(int)

        list_clusters = []
        for c in np.unique(labels):
            indices_c = np.argwhere(labels==c)
            genes_c = np.array(genes)[indices_c]
            list_clusters.append(genes_c)

        if isinstance(mode, list):
            annotations = []
            self.percentage_enriched_genes_ = []
            self.percentage_enriched_clusters = []
            self.other_results = []
            for m in mode:
                annotations.append(_get_database(self.genes, m))
                res = _perform_enrichment(list_clusters, annotations, self.correction,
                                          self.n_cores)
                self.percentage_enriched_genes_.append(res[-2])
                self.percentage_enriched_clusters_.append(res[-1])
                self.other_results_.append(res[:-2])
        else:
              annotations = _get_database(self.genes, self.mode)
              res = _perform_enrichment(list_clusters, annotations, self.correction,
                                        self.n_cores)
              self.percentage_enriched_genes_ = res[-2]
              self.percentage_enriched_clusters_. = res[-1]
              self.other_results_ = res[:-2]
