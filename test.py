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
    indices = [i for i, s in enumerate(df.columns) if cond1 in s and cond2 in s]\
              + [60,61,62,63,64,65]

    df = df.iloc[:, indices]
    return df

def generate_dict(df, clusters, save, dir):
    dict = {}
    for i in range(0,5):
        cond1_ = 'wt_d'+str(i)
        cond2_ = 'TET3KO_d' + str(i)
        for j in clusters:
            name = cond1_ + '_v_' + cond2_ + '_clust' + str(j)
            dict[name] = extract_cluster(df, cond1_, cond2_, j)
            if save == True:
                dict[name].to_csv(dir +'/'+ name, sep='\t')
    return dict

def save_dict(df, clusters, dir):
    dict = generate_dict(df, clusters)
    for i in list(dict.keys()):
        dict[i].to_csv(dir + '/' + i, sep='\t')
    del dict