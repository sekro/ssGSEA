"""
GSEA score calculation with multi-processing
as used in:
Deep phenotyping of the prostate tumor microenvironment reveals molecular stratifiers of relapse linked to inflammatory
chemokine expression and aberrant metabolism
Sebastian Krossa, Maria K. Andersen, Elise Midtbust, Maximilian Wess, Antti Kiviaho, Abhibhav Sharma, Trond Viset,
Øystein Størkersen, Guro F. Giskeødegård, Matti Nykter, Alfonso Urbanucci, Morten B. Rye, May-Britt Tessem
bioRxiv 2024.05.13.593822; doi: https://doi.org/10.1101/2024.05.13.593822

Implementation based on:
Rye, M.B., et al., The genes controlling normal function of citrate and spermine secretion are lost in aggressive
prostate cancer and prostate model systems. iScience, 2022. 25(6): p. 104451.

based on:
Markert, E.K., et al., Molecular classification of prostate cancer using curated expression signatures.
Proc Natl Acad Sci U S A, 2011. 108(52): p. 21276–81.

@author
Sebastian Krossa / NTNU-MH-ISB / Nov. 2023 - Nov. 2024
sebastian.krossa@ntnu.no

"""

from typing import List, Tuple, Dict, Union

from multiprocessing import Pool
import time
import pandas as pd
import numpy as np
import os
import json
import argparse
import sys
from scipy.stats import zscore


from typing import List, Tuple, Dict, Union


def calc_one_sided_pval(value, dist):
    if isinstance(dist, list):
        _dist = np.array(dist)
    elif isinstance(dist, np.ndarray):
        _dist = dist
    else:
        print('Warning invalid input for dist')
        return None
    if value >= 0:
        return np.where(_dist > value, 1, 0).sum() / len(_dist)
    else:
        return np.where(_dist < value, 1, 0).sum() / len(_dist)


def gen_permutation_matrix(permutations, n_genes, n_geneset_genes, max_tries=10000):
    _perm_array = None
    rng = np.random.default_rng()
    _base = np.zeros(n_genes).astype(int)
    _base[:n_geneset_genes] = 1
    #print(_base)
    for i in range(max_tries):
        _try = rng.permutation(_base)
        if _perm_array is None:
            _perm_array = _try
        else:
            if not (_perm_array[:, None] == _try).all(-1).any():
                _perm_array = np.vstack((_perm_array, _try))
        if len(_perm_array.flatten()) >= permutations * n_genes:
            break
    if i+1 >= max_tries:
        print('Warning max tries = {} reached, found {} unique permutations'.format(max_tries, len(_perm_array)))
    return _perm_array


def gen_permutation_score_dist(n_permutations, n_genes, n_geneset_genes, posv, negv, alt_mode=True):
    _distribution = []
    _pa = gen_permutation_matrix(permutations=n_permutations, n_genes=n_genes, n_geneset_genes=n_geneset_genes)
    for _p in _pa:
        _tmp_score, _ = calc_GSEA_score(n_genes=n_genes, ovlp=_p, posv=posv, negv=negv, alt_mode=alt_mode)
        _distribution.append(_tmp_score)
    return _distribution


def calc_GSEA_score(n_genes, ovlp, posv, negv, alt_mode=True):
    _result = None
    _scrv = np.zeros(n_genes)  # vector with 0 -> stores the accumulated score for each gene pos?
    _previous = 0  # accmulates the score
    for i in range(len(ovlp)):
        if ovlp[i]:
            _scrv[i] = _previous + posv
        else:
            _scrv[i] = _previous + negv
        _previous = _scrv[i]
    if alt_mode:
        # This is an alternative mode to calculate the GSEA score
        # gives more "finegrading" for gene sets located in in the "middel" of
        # the ranked list then using just the max pos/neg value
        _result = _scrv.max() + _scrv.min()
    else:
        # This is how it is done according to some papers I found
        if abs(_scrv.max()) >= abs(_scrv.min()):
            _result = _scrv.max()
        else:
            _result = _scrv.min()
    return _result, _scrv


def calc_GSEA_score_for_gset(data: Union[pd.DataFrame, np.ndarray],
                             sig_genes: List[str],
                             gene_names: List[str] = None,
                             sample_names: List[str] = None,
                             mean_centering: bool = False,
                             alt_mode: bool = True,
                             normalize: bool = True,
                             n_permutations: Union[int, None] = None) -> Tuple[Dict[str, float], Dict[str, np.ndarray], List]:
    """
    Calculated GSEA scores for each sample (rows) for a gene set found in genes (cols) in the input pandas dataframe
    :param n_permutations: optinal, if not none calculates n permuations of random GSEA scores
    :param sample_names: optional list of sample names - only needed for numpy data
    :param gene_names: list of gene names - only needed for numpy data input
    :param data: input data as pandas DataFrame or numpy array - samples in rows, genes in columns, provide at least
     gene_names when using numpy array as input
    :param sig_genes: Genes in the signature provided as list of str - same gene ids as in input data
    :param mean_centering: Perform mean centering of expression data
    :param alt_mode: return max + min GSEA score if True else return max deviation from 0
    :param normalize: return scores normalized to range(-1, 1)
    :return: dict of scores ({'sample_name': score, ...}), dict of running score vectors ({'sample_name': vector, ...})
    """
    # original GSEA method as defined in Mootha et al.
    # mean centering needed? Rather use zscores as input
    if isinstance(data, np.ndarray) and gene_names is not None:
        _df = pd.DataFrame(data=data, index=sample_names, columns=gene_names)
    elif isinstance(data, pd.DataFrame):
        _df = data
    else:
        print("input data not ok")
        return None, None, None
    if mean_centering:
        # this should be over genes (gene mean centered) not samples - check!!!
        _df = _df.apply(lambda x: x - x.mean())
    N = len(_df.columns)  # number of genes in dataset
    # Number of probes matching gene-set
    SIG = np.isin(_df.columns, sig_genes).sum()
    # negv and posv values for normalization into -1..1 range:
    # max_posv = 1 = no of sig genes * posv
    # max_negv = -1 = (no genes in data - no of sig genes) * negv
    if normalize:
        negv = - 1 / (np.double((N - SIG)))
        posv = 1 / np.double(SIG)
    else:
        negv = -np.sqrt((np.double(SIG)) / (N - SIG))  # scaled or normed factor gene not in geneset
        posv = np.sqrt(np.double((N - SIG)) / (SIG))  # scaled or normed factor gene in geneset
    ncnt = np.invert(np.isin(sig_genes, _df.columns))  # not really used for anything
    results = {}  # results / GSEA score dict
    scrvs = {}  # score vector dict - for plotting / checking...
    for ri, rd in _df.iterrows():
        ovlp = np.isin(rd.sort_values(ascending=False).index, sig_genes)
        results[ri], scrvs[ri] = calc_GSEA_score(n_genes=N, ovlp=ovlp, posv=posv, negv=negv, alt_mode=alt_mode)
    _distribution = None
    if n_permutations is not None and isinstance(n_permutations, int):
        _distribution = gen_permutation_score_dist(n_permutations=n_permutations, n_genes=N, n_geneset_genes=SIG,
                                                   posv=posv, negv=negv, alt_mode=alt_mode)
    return results, scrvs, _distribution


def get_gene_cols(df, substring_to_remove='_goi'):
    """
    Returns a list of column names from df.

    Input:
        df: Dataframe
        appendix_to_remove: substring that will be filtered out of column names.

    Returns list of column names.
    """
    cols = df.columns[np.char.find(df.columns.to_list(), 'ENSG') != -1]
    if substring_to_remove is not None:
        cols = [x.replace(substring_to_remove, '') for x in cols]
    return cols


def apply_z_score(df, alt_input_reads=False):
    if alt_input_reads:
        df = df.apply(zscore)
    else:
        df.loc[:,get_gene_cols(df)] = df.loc[:,get_gene_cols(df)].apply(zscore)
    df.fillna(0, inplace=True) # Genes with 0 expression will result in nan, so we need to replace those with 0.
    return df


def mp_wrapper(wrapper_args):
    arg, kwargs = wrapper_args
    print('spawning process for {}'.format(arg))
    return arg, calc_GSEA_score_for_gset(**kwargs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='CalcGSEA')
    parser.add_argument("input_reads", help="path to input data - expect tab-sep raw RNA reads, samples in columns, genes in rows. Genes need to have ENSG in the name otherwise row will be ignored")
    parser.add_argument("input_sigs_json",
                        help="signature json file")

    parser.add_argument("output", help="folder for output files")
    parser.add_argument("--alt_input_reads", dest='alt_input_reads', action='store_true',
                        help="use this argument if input data is organized as samples in rows, genes in columns of reads or microarray data")
    parser.add_argument("--v", dest='verbose', action='store_true',
                        help="output a lot! info on screen")
    parser.add_argument("--n_procs", help="number of processes to be used during ssGSEA calculation", default=8, type=int)
    parser.set_defaults(alt_input_reads=False)
    parser.set_defaults(verbose=False)
    args = parser.parse_args()
    start_time = time.time()
    if args.verbose:
        print('Verbose mode on!')
        print('Using {} processe(s)'.format(args.n_procs))
    # check and make out folder
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    _sig_d = None
    if os.path.isfile(args.input_sigs_json):
        with open(args.input_sigs_json) as json_file:
            _sig_d = json.load(json_file)
    if os.path.isfile(args.input_reads) and _sig_d is not None:
        if args.alt_input_reads:
            print('Using alternative input data format!')
            if 'xlsx' in args.input_reads:
                df = pd.read_excel(args.input_reads, index_col=0).fillna(0)
                #df_n_z = apply_z_score(df.copy(), alt_input_reads=args.alt_input_reads)
            elif 'csv' in args.input_reads:
                df = pd.read_csv(args.input_reads, index_col=0, header=0, sep=',').fillna(0)
                print('Input is csv - assuming index in column 0, header in row 0, separator = ,')
            else:
                sys.exit("Wrong input format - Expecting input data to be filetype xlsx")
            df_n_z = apply_z_score(df.copy(), alt_input_reads=args.alt_input_reads)
        else:
            df = pd.read_csv(args.input_reads, sep='\t', index_col=0).fillna(0)
            df_n = (df / df.sum()) * df.sum().mean()
            df_n_z = apply_z_score(df_n.T.copy())
        data_genes = df_n_z.columns.to_list()
        sig_d = {}
        merge_d = {}
        for _gsn, _gsgs in _sig_d.items():
            _tgs = []
            for _g in _gsgs:
                if args.alt_input_reads:
                    # use the gene name
                    _g = _g.split(' - ')[0]
                else:
                    # use ENSG...
                    _g = _g.split(' - ')[1]
                if _g in data_genes:
                    _tgs.append(_g)
            sig_d[_gsn] = _tgs
            if '_' in _gsn:
                _pre = _gsn.split('_')[0]
                _suf = _gsn.split('_')[1]
                if _suf in ['down', 'up']:
                    if _pre in merge_d:
                        merge_d[_pre][_suf] = _gsn
                    else:
                        merge_d[_pre] = {}
                        merge_d[_pre][_suf] = _gsn
        if args.verbose:
            print('---head of input data below---')
            print(df_n_z.head())
            print('---head of input data above---')
            print('---description of 5 first columns of input data below---')
            print(df_n_z.iloc[:, 0:5].describe())
            print('---description of input data above---')
            print('input data should show a few rows of samples with genes (ENSG...) or gene names as columns')
            print('description of input data should confirm *zscored* data: mean = 0, std = 1')
            print('if you see no data or wrong format check apply_z_score() and get_gene_cols()')
        concat_l = []
        if args.n_procs > 1:
            mp_args = []
            _df_split = np.array_split(df_n_z, int(args.n_procs / len(sig_d.keys())))
            for gset_n, gset_g in sig_d.items():
                print('---------------------------------------------------------------')
                print('Setting up MP data for {} using {} of {} genes in the signature: {}'.format(gset_n, len(gset_g),
                                                                                       len(_sig_d[gset_n]), gset_g))

                for _df_chunk in _df_split:
                    mp_args.append((gset_n, {'data': _df_chunk, 'sig_genes': gset_g, 'mean_centering': False,
                                           'alt_mode': True, 'n_permutations': 1000}))
            with Pool(args.n_procs) as p:
                _mp_res = p.map(mp_wrapper, mp_args)

            mp_res_data = {}
            for gset_n, _ret_data in _mp_res:
                if gset_n not in mp_res_data:
                    mp_res_data[gset_n] = {
                        'data': {},
                        'scrvs': {},
                        'dist': []
                    }
                _da, _sc, _di = _ret_data
                for k, v in _da.items():
                    mp_res_data[gset_n]['data'][k] = v
                for k, v in _sc.items():
                    mp_res_data[gset_n]['scrvs'][k] = v
                mp_res_data[gset_n]['dist'].extend(_di)
            for gset_n, _d in mp_res_data.items():
                concat_l.append(pd.Series(data=_d['data'], name=gset_n))
        else:
            scrvs_d = {}
            dist_d = {}
            for gset_n, gset_g in sig_d.items():
                print('---------------------------------------------------------------')
                print('Working on {} using {} of {} genes in the signature: {}'.format(gset_n, len(gset_g),
                                                                                       len(_sig_d[gset_n]), gset_g))
                # work with zscores! not reads / normalized reads
                _data, _scrvs, _dist = calc_GSEA_score_for_gset(data=df_n_z, sig_genes=gset_g,
                                                                mean_centering=False, alt_mode=True,
                                                                n_permutations=1000)
                scrvs_d[gset_n] = _scrvs
                dist_d[gset_n] = _dist
                concat_l.append(pd.Series(data=_data, name=gset_n))

        df_GSEA = pd.concat(concat_l, axis=1)
        if args.verbose:
            print('---Results pre-merge---')
            print(df_GSEA.head())
        # WARNING: This is quick and dirty - expecting fixed names in the sig_d dict!!
        for _mn, _parts in merge_d.items():
            df_GSEA[_mn] = 1000 * (df_GSEA[_parts['up']] - df_GSEA[_parts['down']]) / 2
        if args.verbose:
            print('---Results post-merge---')
            print(df_GSEA.head())
        df_GSEA.to_csv(os.path.join(args.output, 'GSEA_scores.csv'))
        print('done - results data file in {} - execution time {} seconds'.format(os.path.join(args.output, 'GSEA_scores.csv'), time.time() - start_time))
