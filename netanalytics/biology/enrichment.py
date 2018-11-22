import os
import math
import warnings 
import multiprocessing
import matplotlib

import pandas as pd
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt 

from scipy.stats import hypergeom
from collections import defaultdict
from joblib import Parallel, delayed
from functools import partial

from nmtf.nmtf import SSNMTF
from nmtf.utils import get_clusters

from netanalytics.biology.plot import plot_enrichment_results


# Benjamini-Hochberg p-value correction for multiple hypothesis testing
def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def _shuffle_annotations(annotations):
    genes = np.unique(annotations.index.values)
    indices = np.arange(0, genes.shape[0])
    shuffled = indices.copy()
    np.random.shuffle(shuffled)
    new_annotations = pd.DataFrame(columns=['annotations'])
    for i in range(indices.shape[0]):
        annot = annotations.loc[genes[indices[i]]]
        index = [genes[shuffled[i]] for j in range(annot.shape[0])]
        if annot.shape[0] == 1:
            data = [annot['annotations']]
        else:
            data = [[a] for a in np.array(annot['annotations'])]
        df = pd.DataFrame(data, columns=['annotations'])
        df.index = index
        new_annotations = new_annotations.append(df)
    return new_annotations


def _get_number_of_enriched_annotations(clusters_annotations):
    annots = []
    for c in clusters_annotations:
        for a in c:
            annots.append(a[0])
    return np.unique(annots).shape[0]


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
        genes_at_least_one = list(set(annotations.index.values).intersection(
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
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
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
        return _perform_enrichment_bonferroni(clusters, annotations,
                                           n_cores)
    else:
        return  _perform_enrichment_bh(clusters,annotations, n_cores)


def _get_database(genes, mode='go'):
    path = os.path.realpath(__file__)
    path = path.split('/')[:-1]  # if not unix/OX should raise exception
    path = '/'.join(path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if mode.lower()=='go_bp' or mode.lower() =='go_mf' or mode.lower()== 'go_cc':
            go = pd.read_table(path+"/databases/GO_most_specific.gaf", 
                  skiprows=30, header=None)
            go = go.set_index(2)
            go = go[(go[6] == 'EXP') | (go[6] == 'IDA') | 
                    (go[6] == 'IPI') | (go[6] == 'IMP')]
            if mode.lower()=='go_bp':
                go = go[go[8]=='P']
            elif mode.lower()=='go_mf':
                go = go[go[8]=='F']
            else:
                go = go[go[8]=='C']
            go.index = [str(s).lower() for s in go.index]
            annotations = go[4].to_frame()
            annotations.columns= [1]
        elif mode.lower() == 'kegg_symbol':
            annotations = pd.read_csv(path+"/databases/kegg_pathways.csv", 
                                      index_col=0)
            annotations.index = [i.lower() for i in annotations.index]
        elif mode.lower() == 'kegg_entrez':
            annotations = pd.read_csv(path+"/databases/kegg_entrez_pathways.csv", 
                                      index_col=0)
        elif mode.lower() == 'reactome_symbol':
            annotations = pd.read_csv(path+"/databases/reactome_pathways.csv", 
                                      index_col=0)
            annotations.index = [i.lower() for i in annotations.index]
        elif mode.lower() == 'reactome_entrez':
            annotations = pd.read_csv(path+"/databases/reactome_entrez_pathways.csv", 
                                      index_col=0)
        elif mode.lower()=='go':
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

        index = 1 if not(mode.lower() =='kegg_symbol' or mode.lower() =='reactome_symbol'or
                         mode.lower() =='kegg_entrez' or mode.lower() =='reactome_entrez') else '1'
        for i, g in enumerate(annotations[index]):
            if annotations[annotations[index]==g].shape[0] != 1:
                to_keep.append(i)
        annotations = annotations.iloc[to_keep,:]
        to_drop=[]
        for g in annotations.index.values:
            if annotations.loc[g].shape[0] == 0:
                to_drop.append(g)
        annotations.drop(to_drop, inplace=True)
        annotations  = pd.DataFrame(annotations[index])
        annotations.columns = ['annotations']
    return annotations


def enrichment_significance(clusters, annotations, n_reps = 100, compute_true=False, 
                            no_enriched_ann=None, n_cores=1, verbose=0):
    """
    Params
    ------
    clusters: list, 
        The list of clusters. Each cluster is a list of genes name/id.
    annotations: pandas.DataFrame 
        A dataframe with one columns called 'annotations' and as index the names/ids of the genes. 
    n_reps: int, optional default=100
        The number of times the procedure is repeated. 
    compute_true: boolean, optional default=False
        If True the function computes the enrichment without reshuffling. 
        Note that if False the value of the enrichment (true_enrichment) must be passed. 
    true_enrichment: float, optional default=None
        It is the value of the enrichment obtained without reshuffling. 
        If compute_true is False it must be passed. 
    n_cores: int, optional default=1
        The number of cores to perform the computation of random enrichment in parallel. 
    verbose: boolean, optional default=0
        If True during the computation some progress messages are printed.
        If the process is done in parallel the progress is not shown. 
    
    Returns
    -------
    array-like, shape=(n_reps, )
        The enrichment value obtained at each repetition
    float:
        P-value. The significance of the true enrichment against the reshuffling. 
    """
    if n_reps < 100:
        warnings.warn("It is better to have a number of repetitions higher then 100.")
    if compute_true:
        enriched_ann, _, _, _ = _perform_enrichment_bh(clusters, annotations)
        no_enriched_ann = _get_number_of_enriched_annotations(enriched_ann)
    elif no_enriched_ann is None:
        raise ValueError("You must indicate either to compute the true enrichment or "
                         "provide the true value as input of the function")
    
    def _for_parallel(clusters, annotations):
        annot = _shuffle_annotations(annotations)
        enriched_ann, _, _, _ = _perform_enrichment_bh(clusters, annot)
        return _get_number_of_enriched_annotations(enriched_ann)
    _fp = partial(_for_parallel, clusters, annotations)
    
    if n_cores > 1:
        numbers_annotations = Parallel(n_jobs=n_cores)(delayed(_fp)() for i in range(n_reps))
    else:
        numbers_annotations = np.zeros(n_reps)
        progress = 1
        for i in range(n_reps):
            annot = _shuffle_annotations(annotations)
            enriched_ann, _, _, _ = _perform_enrichment_bh(clusters, annot)
            numbers_annotations[i] = _get_number_of_enriched_annotations(enriched_ann)
            if verbose and i%(n_reps//10) == 0:
                print('Done %d /100'%(progress*10))
                progress +=1
    pval = (np.sum(np.array(numbers_annotations) >= no_enriched_ann) / n_reps)
    return numbers_annotations, pval


def enrichment_analysis_multiple_graphs(graphs, k, genes, results_folder, networks_type, 
                                        _type, 
                                        labels_plot=['GO-BP', 'GO-MF', 'GO-CC', 'KP', 'RP'],
                                        normalize=True,
                                        specific_file_description="",
                                        enrichment_types=['go_bp', 'go_mf', 'go_cc', 
                                                         'kegg_symbol', 'reactome_symbol'],
                                        compute_significance=True, n_cores=1):
    if normalize:
        graphs = graphs.copy()
        graphs= [(g - np.mean(g))/np.std(g) for g in graphs]
    
    def _parallel_SSNMTF(g, i, k, max_iter):
        est = SSNMTF(k=k, max_iter=max_iter)
        est.fit([g])
        return i, est
    pps = partial(_parallel_SSNMTF, k=k, max_iter=1000)
    results = Parallel(n_jobs=n_cores)(delayed(pps)(g, i) for i, g in enumerate(graphs))
          
    ests = [None]*len(graphs)
    for r in results:
        ests[r[0]] = r[1]

    esti = SSNMTF(k=k, max_iter=1000)
    esti.fit(graphs)
    ests.append(esti)
    
    enr = Enrichment(genes, mode=enrichment_types,
                     compute_significance=compute_significance, n_cores=n_cores)
    enr_genes = []
    enr_clusters = []
    if compute_significance:
        significances = []
    
    for est in ests:
        labels = get_clusters(est.G_)
        enr.fit(labels)
        enr_genes.append(enr.percentage_enriched_genes_)
        enr_clusters.append(enr.percentage_enriched_clusters_)
        if compute_significance:
            significances.append(enr.p_value)
    
    with open(results_folder+"/genes_enriched_"+specific_file_description+".pkl", 'wb') as f:
        pkl.dump(enr_genes, f)
    with open(results_folder+"/clusters_enriched_"+specific_file_description+".pkl", 'wb') as f:
        pkl.dump(enr_clusters, f)
    if compute_significance:
        with open(results_folder+"/enrichment_p_value_"+specific_file_description+".pkl", 'wb') as f:
            pkl.dump(significances, f)
        
    matplotlib.rcParams.update({'font.size': 22})
   
    networks_type += ['integrated']
    plot_enrichment_results(enr_genes, results_folder+"/enriched_genes_"+specific_file_description+".pdf",
                            networks_type, labels_plot, "% Enriched genes", _type+' Networks')
    plot_enrichment_results(enr_clusters, results_folder+"/enriched_clusters_"+specific_file_description+".pdf",
                            networks_type, labels_plot, "% Enriched clusters", _type+' Networks')
    

class Enrichment():
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

    def __init__(self, genes, mode='go', correction='bh', 
                 compute_significance=True, n_repetitions=100,
                 n_cores=1, verbose=0):
        self.genes = genes
        self.mode = mode
        self.correction = correction
        self.compute_significance=compute_significance
        self.n_repetitions=n_repetitions
        self.n_cores = n_cores
        self.verbose = verbose

    def fit(self, labels):
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
            genes_c = np.array(self.genes)[indices_c]
            list_clusters.append(genes_c)

        if isinstance(self.mode, list): #TODO if one mode is kegg/reactome fare ids mapping 
            self.percentage_enriched_genes_ = []
            self.percentage_enriched_clusters_ = []
            self.mean_normalized_entropy = []
            self.enriched_annotations_list = []
            def _parallel_enrichment(mode, genes=None, clusters=None, correction=None,
                                     n_cores=1):
                annotations = _get_database(self.genes, mode)
                res = _perform_enrichment(list_clusters, annotations, self.correction,
                                          self.n_cores)
                return mode, res
            
            ppe = partial(_parallel_enrichment, genes=self.genes, clusters=list_clusters, 
                          correction=self.correction, n_cores=self.n_cores)
            results = Parallel(n_jobs=self.n_cores)(delayed(ppe)(m) for m in self.mode)  
            
            for m in self.mode:
                for mode, res in results:
                    if mode==m:
                        self.enriched_annotations_list.append(res[0])
                        self.percentage_enriched_genes_.append(res[-2])
                        self.percentage_enriched_clusters_.append(res[-1])
                        self.mean_normalized_entropy.append(res[1])
                        break
        else:
            annotations = _get_database(self.genes, self.mode)
            res = _perform_enrichment(list_clusters, annotations, self.correction,
                                    self.n_cores)
            self.enriched_annotations_list = res[0]
            self.percentage_enriched_genes_ = res[-2]
            self.percentage_enriched_clusters_ = res[-1]
            self.mean_normalized_entropy = res[1]
        
        if self.compute_significance:
            if isinstance(self.mode, list):
                values = [_get_number_of_enriched_annotations(r) for r in self.enriched_annotations_list]
                self.shuffled_annotations_enriched= []
                self.p_value = []
                for i, m in enumerate(self.mode):
                    annotations = _get_database(self.genes, m)
                    res_s = enrichment_significance(list_clusters, annotations,
                                                  n_reps=self.n_repetitions,
                                                  no_enriched_ann=values[i], 
                                                  n_cores=self.n_cores,
                                                   verbose=np.max(self.verbose-1, 0))
                    self.shuffled_annotations_enriched.append(res_s[0])
                    self.p_value.append(res_s[1])
            else:
                value = _get_number_of_enriched_annotations(self.enriched_annotations_list)
                self.shuffled_annotations_enriched, self.p_value = \
                    enrichment_significance(list_clusters, annotations,
                                            n_reps=self.n_repetitions,
                                            no_enriched_ann=value, n_cores=self.n_cores,
                                            verbose=np.max(self.verbose-1, 0))
