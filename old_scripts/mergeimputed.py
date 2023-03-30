import numpy as np
import pandas as pd

def main():
    #import NOT imputed, log2 transformed, filtered from file
    
    normal = pd.read_csv(getfile(), sep = "\t", index_col = "T: Protein.Group")
    df1 = cleancopy(normal)

    #import ALL imputed, log2 transformed, filtered from file
    allimputed =  pd.read_csv(getfile(), sep = "\t", index_col = "T: Protein.Group")
    df2 = cleancopy(allimputed)
    
    df = df1.mask(df1.groupby('group', axis=1).count() == 0, 1)
    df = df.where(~df.isna(), df2)

    #remove MultiIndex
    df.columns = df.columns.droplevel(["group"])
    
    #adds the new imputed columns back to original dataframe
    normal[df.columns] = df[df.columns]
    
    #saves dataframe
    normal.to_csv(getfile(), decimal = ",", sep = "\t")   

# returns copy of a, drops unused protein descriptors, adds Multiindex
def cleancopy(a: pd.DataFrame):
    b = pd.DataFrame(a)
    b.drop(columns = ["T: Protein.Ids","T: Protein.Names", "T: Genes", "T: First.Protein.Description"], inplace=True)
    
    #adds a new line of values used for grouping : "group"
    #extracted from name (must be format dx_y with x:day y:replicate)
    new_idx = pd.MultiIndex.from_arrays([
        b.columns,
        b.columns.str.extract("(d\d+)_\d+", expand = False)
    ], names=["index", "group"])
    b.columns = new_idx
    return b   

#opens explorer window and lets you input the file
def getfile():
    import tkinter as tk
    from tkinter import filedialog
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename()
    return file_path

if __name__ == "__main__":
    main()