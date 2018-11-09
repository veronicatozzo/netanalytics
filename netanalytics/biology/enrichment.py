import pandas as pd
import numpy as np

from scipy.stats import hypergeom


def _perform_enrichment(clusters, dims, annotations):
    N = len(np.unique(annotations.index))
    bonferroni_correction = len(np.unique(annotations['annotations']))*len(dims)

    p_values = []
    genes_to_count = []
    enriched_clusters = 0
    for i, c in enumerate(clusters):
        enriched=False
        genes_at_least_one = list(set(annotations.index).intersection(
                                  set(list(c.ravel()))))
        n = len(genes_at_least_one)
        ann_c = annotations.loc[genes_at_least_one]
        for a in np.unique(ann_c[4]):
            genes_annotated = annotations[annotations['annotations']==a]
            K = genes_annotated.shape[0]
            k = ann_c[ann_c[4]==a].shape[0]
            pval = hypergeom.sf(k-1, N, K, n)
            if pval < (0.05/bonferroni_correction):
                genes_to_count += \
                        list(np.unique(ann_c[ann_c['annotations']==a].index))
                p_values.append((i, pval, a))
                enriched=True
        if enriched:
            enriched_clusters+=1
    return p_values, len(set(genes_to_count))/N, enriched_clusters/len(clusters)


def _get_database(mode='go'):
    if mode.lower()=='go5':
        go_level5 = pd.read_csv("GO_BP_level5.csv", index_col=0)
        intersection = list(set(genes).intersection(set(go_level5.index)))
        annotations = go_level5.loc[intersection]
        annotations = pd.DataFrame(annotations['GOterm'])
        annotations.columns = ['annotations']
    elif mode.lower()=='kegg':
        pathways = pd.read_csv("pathways_kegg.csv", index_col=0)
        annotations  = pd.DataFrame(pathways['pathway'])
        annotations.columns = ['annotations']
    elif mode.lower()=='go':
        go_annotations = pd.read_table("GO_most_specific.gaf", sep='\t',
                                       skiprows=30, header=None)
            go_annotations = go_annotations.set_index(2)
            go_annotations = go_annotations[(go_annotations[6] == 'EXP') |
                                        (go_annotations[6] == 'IDA') |
                                        (go_annotations[6] == 'IPI') |
                                        (go_annotations[6] == 'IMP')]
            go_annotations = go_annotations[go_annotations[8]=='P']
            go_annotations.index = [str(s).lower()
                                    for s in go_annotations.index]
            intersection = list(set(genes).intersection(
                                set(go_annotations.index)))
            annotations = go_annotations.loc[intersection]
            annotations = pd.DataFrame(annotations[4])
            annotations.columns = ['']
     elif:
         raise ValueError("The mode you specified is not implemented, please "
                          "try one between go, go5 or kegg")
    return annotations


def enrichment(genes, labels, type='go'):
    """
    Parameters
    ----------
    genes: list,  length=n
        List of HUGO symbol genes.
    labels: list, length=n
        List of cluster lables for each gene.

    Returns
    -------
    list:
        List of enriched annotations. Each entry is a tuple that has the label
        of the cluster, the p-value and the annotation.
    float:
        The percentage of genes that are enriched at least for one annotation
        in their cluster.
    float: 
        The percentage of clusters that are enriched at least once.
    """

    if (len(genes) != len(labels)):
        raise ValueError("The length of genes list and cluster labels must be "
                         "the same, found %d and %d respectively."
                         %(len(genes), len(labels)))

    annotations = _get_database(mode)

    clusters_dim = []
    list_clusters = []
    for c in np.unique(labels):
        indices_c = np.argwhere(labels==c)
        genes_c = np.array(genes)[indices_c]
        list_clusters.append(genes_c)
        clusters_dim.append(len(genes_c))
    return _perform_enrichment(list_clusters, clusters_dim, annotations)
