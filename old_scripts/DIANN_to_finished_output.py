import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog
def main(impute):
    #import DIA-NN output
    og = pd.read_csv(getfile(), sep = ";", index_col = "Protein.Group")

    df1 = pd.DataFrame(og)
    df1.drop(columns=["Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description"], inplace=True)

    # adds a new line of values used for grouping : "group"
    # extracted from name (must be format dx_y with x:day y:replicate)
    new_idx = pd.MultiIndex.from_arrays([
        df1.columns,
        df1.columns.str.extract("(d\d+)_\d+", expand=False)
    ], names=["index", "group"])
    df1.columns = new_idx

    #filters invalid values
    filt = get_invalid(df1, percent=0.7)
    df1.drop(df1[filt].index, inplace = True)

    #log 2 transform every value
    df1 = df1.applymap(np.log2)

    #impute all missing values in copy of df1
    df2 = impute(df1, width=0.3, downshift=1.8)

    #define MNAR values -> impute with 1
    df = df1.mask(df1.groupby('group', axis=1).count() == 0, 1)
    df = df.where(~df.isna(), df2)

    #remove MultiIndex
    df.columns = df.columns.droplevel(["group"])
    
    #adds the new imputed columns back to original dataframe
    og[df.columns] = df[df.columns]
    
    #save dataframe
    try:
        # with block automatically closes file
        with filedialog.asksaveasfile(mode='w', defaultextension=".txt") as file:
            og.to_csv(file.name, sep = "\t")
    except AttributeError:
        # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
        print("The user cancelled save")

#returns a series of Bool -> if none of the values in a group have at least x(70)% valid values
def get_invalid(a: pd.DataFrame, percent: float = 0.7):
    cutoff = (1 - percent)
    b = cutoff < (a.isnull().groupby("group", axis = 1).sum()/a.groupby("group", axis = 1).size())
    c = b.all(axis = 1)
    return c 

def impute_custom(a: pd.DataFrame, width: float = 0.3, downshift: float = 1.8):
    b = pd.DataFrame(a)
    for i in b.columns:
        mean = np.mean(b[i])-(downshift*np.std(b[i]))
        std = np.std(b[i])*width
        b[i] = b[i].isna().apply(lambda v: float(np.random.normal(mean, std, 1)))
        
    return b

#opens explorer window and lets you input the file
def getfile():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename()
    return file_path

if __name__ == "__main__":
    main(impute_custom)