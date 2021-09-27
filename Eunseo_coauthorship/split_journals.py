import pandas as pd
import numpy as np

def main(path):
    df_all = pd.read_csv(path)
    split_dfs = np.array_split(df_all, 10)
    for i in range(len(split_dfs)):
        split_dfs[i].to_csv("2020_split_journals/2020_journals_" + str(i) + ".csv")

if __name__ == "__main__":
    main("2020_all_journal.csv")
