import pandas as pd
import os
import numpy as np

# cluster = pd.read_csv("cluster_means.txt", sep = ",")
#
# Original = pd.read_csv("tet_v_wt_difference_scores.txt", sep = ",")
# Original.set_index('ID', inplace = True)
#
# cluster = cluster[cluster.columns[0:2]]
# cluster.rename(columns ={cluster.columns[0]: "ID"}, inplace = True)
# cluster.set_index('ID', inplace = True)
#
# merged = pd.merge(Original, cluster, left_index=True, right_index=True, how="left")
# merged = pd.read_csv('merged.txt', sep = '\t')

def save_clusters(df,
                  cond1 = 'wt',
                  cond2 = 'TET3KO',
                  days = range(0,5),
                  clusters = range(1,9),
                  save = True,
                  dir = 'wt_vs_TET3KO/clusters'):

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

