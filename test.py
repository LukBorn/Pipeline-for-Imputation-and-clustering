import pandas as pd
import numpy as np

#cluster = pd.read_csv("Cluster_Means_significant_proteins.txt", sep = "\t")
#Original = pd.read_excel("Proteomics_iNGNs_DIA_R_DataOutput_impute_1+Gaussian_with-MeCP2.xlsx")

#cluster = cluster[cluster.columns[0:2]]
#cluster.rename(columns ={cluster.columns[0]: "ID"}, inplace = True)
#cluster.set_index('ID', inplace = True)


#merged = pd.merge(Original, cluster, left_index=True, right_index=True, how="left")
merged = pd.read_csv('merged.txt', sep = '\t')
def extract_cluster(df, cond1, cond2, cluster:int):
    df = df[df['cluster'] == cluster]
    relevant_values = [i for i, s in enumerate(df.columns) if cond1 in s and cond2 in s]
    IDs = [list(df.columns).index(ID) for ID in
           ['Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description']
           ]
    df = df.iloc[:, relevant_values+IDs]
    return df

def generate_dict(df, clusters, save, dir):
    dict = {}
    for i in range(0,5):
        cond1_ = 'wt_d'+str(i)
        cond2_ = 'TET3KO_d' + str(i)
        for j in clusters:
            subset = extract_cluster(df, cond1_, cond2_, j)
            name = [df.columns[i].rstrip('_CI.L') for i, s in enumerate(df.columns) if 'CI.L' in s][0]\
                   + '_clust' + str(j)
            dict[name] = subset
            if save == True:
                dict[name].to_csv(dir +'/'+ name+ '.txt', sep='\t')
    return dict

def save_dict(df, clusters, dir):
    dict = generate_dict(df, clusters)
    for i in list(dict.keys()):
        dict[i].to_csv(dir + '/' + i + '.txt', sep='\t')
    del dict