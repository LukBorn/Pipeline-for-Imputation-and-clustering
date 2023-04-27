import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
import tkinter as tk
from tkinter import filedialog
import os
import warnings
class Pipeline_tracker:
    """
    class to keep track of pipeline progress
    for example no GSEA before normalization
     """


trk = Pipeline_tracker


def DIANN_to_df():
    ...

def filt_invalid(df,
                 percent):
    cutoff = (1 - percent)
    df = _add_multiindex(df)
    b = cutoff < (df.isnull().groupby("group", axis=1).sum() / df.groupby("group", axis=1).size())
    filt = b.all(axis=1)
    df.drop(df[filt].index, inplace=True)
    return _remove_multiindex(df)

def log2_transform(df):
    return df.applymap(np.log2)

def custom_impute(df,
                  mnar = 0,
                  width = 0.3,
                  downshift = 1.8):
    """

    :param df:
    :param mnar: max number of non-NAN values to be missing not at random
    :param width:
    :param downshift:
    :return:
    """

    df = _add_multiindex(df)

    # define MNAR values -> impute with 1
    df1 = df.mask(df.groupby('group', axis=1).count() == mnar, 1)

    # impute all values for each column seperately in a copy of df
    df_all_imp = pd.DataFrame(df)
    for i in df_all_imp.columns:
        mean = np.mean(df_all_imp[i]) - (downshift * np.std(df_all_imp[i]))
        std = np.std(df_all_imp[i]) * width
        df_all_imp[i] = df_all_imp[i].isna().apply(lambda v: float(np.random.normal(mean, std, 1)))

    df = df1.where(~df1.isna(), df_all_imp)

    return _remove_multiindex(df)

def normalize():
    #normalize all log 2 foldchanges
    ...

def diff_GSEA(df,
              comparisons):

    alpha = 0.05

    df = _add_multiindex(df)

    for comparison in comparisons:
        cond1, cond2 = tuple(comparison.split("_vs_"))
        if cond1 not in df.columns or cond2 not in df.columns:
            warnings.warn(f"the conditions found in {comparison} cant be found in df")
            print(df.columns)

        group1 = df[df.groupby['group'] == cond1]
        group2 = df[df.groupby['group'] == cond2]


        # Separate the data into the two groups
        # group1 = df.iloc[:, group1_cols].values
        # group2 = df.iloc[:, group2_cols].values

        # Compute the Wilcoxon test statistic and p-values for each row
        tstat, pvalue = wilcoxon(group1, group2, axis=1, alternative='two-sided')

        # Compute the log-fold changes
        logFC = np.log2(np.mean(group1, axis=1) / np.mean(group2, axis=1))

        # Adjust the p-values for multiple testing
        _, pvalue, _, _ = multipletests(pvalue, alpha=alpha, method='fdr_bh')

        return pd.Series(logFC, index=df.index), pd.Series(pvalue, index=df.index)

    ...

def save_clusters(df,
                  cond1 = 'wt',
                  cond2 = 'TET3KO',
                  days = range(0,5),
                  clusters = range(1,9),
                  save = True,
                  dir = 'wt_vs_TET3KO/clusters'):
    """
        saves difference and p-value of all proteins for specific clusters
        into seperate .txt files for IPA
        :param df:
        :param cond1:
        :param cond2:
        :param days: range or list of days to save
        :param clusters: range or list of clusters to save
        :param save: saves if true
        :param dir: directory to save to
        :return:
        """

    if not os.path.exists(dir):
        os.makedirs(dir)

    dict = {}
    for day in days:
        for cluster in clusters:
            indexes = [list(df.columns).index(ID) for ID in
                       [f'{cond1}_d{str(day)}._vs_{cond2}_d{str(day)}._diff',
                        f'{cond1}_d{str(day)}._vs_{cond2}_d{str(day)}._p.val',
                        'Protein.Group',
                        'Protein.Ids',
                        'Protein.Names',
                        'First.Protein.Description',
                        'cluster']
                       ]
            subset = df[df['cluster'] == cluster].iloc[:, indexes]
            name = f'{cond1}_vs_{cond2}_d{day}_clust{cluster}'
            dict[name] = subset
            if save == True:
                dict[name].to_csv(f'{dir}/{name}.txt', sep='\t')
    return dict


def save_all(df,
             cond1='wt',
             cond2='TET3KO',
             days=range(0, 5),
             save=True,
             dir='wt_vs_TET3KO'):
    """
    saves difference and p-value of all proteins regardless of cluster
    into seperate .txt files for IPA
    :param df:
    :param cond1:
    :param cond2:
    :param days: range or list of days to save
    :param save: saves if true
    :param dir: directory to save to
    :return:
    """

    if not os.path.exists(dir):
        os.makedirs(dir)

    dict = {}
    for day in days:
        indexes = [list(df.columns).index(ID) for ID in
                   [f'{cond1}_d{str(day)}._vs_{cond2}_d{str(day)}._diff',
                    f'{cond1}_d{str(day)}._vs_{cond2}_d{str(day)}._p.val',
                    'Protein.Group',
                    'Protein.Ids',
                    'Protein.Names',
                    'First.Protein.Description',
                    'cluster']
                   ]
        subset = df.iloc[:, indexes]
        name = f'{cond1}_vs_{cond2}_d{day}_all'
        dict[name] = subset
        if save == True:
            dict[name].to_csv(f'{dir}/{name}.txt', sep='\t')
    return dict


def _add_multiindex(df):
    # adds a new line of values used for grouping : "group"
    # extracted from name (must be format cond_dx_y with x:day y:replicate)
    new_idx = pd.MultiIndex.from_arrays([
        df.columns,
        df.columns.str.extract("([a-zA-Z0-9]+_d\d+)_\d+", expand=False)
    ], names=["index", "group"])
    df.columns = new_idx
    return df

def _remove_multiindex(df):
    df.columns = df.columns.droplevel(["group"])
    return df