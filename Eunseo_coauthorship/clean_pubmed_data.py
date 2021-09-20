import numpy as np
import pandas as pd
import geograpy
import time
from datetime import timedelta
import multiprocessing as mp
import sys, getopt

'''
Co-Authorship Network Project
Eunseo Sung
Cleaning journal info from pubmed
    - Countries and cities/states
    - Filter authors with <2 papers
'''

# Input: author string (firstname$lastname$affiliation)
# Output: string firstname.lastname.country.city
def get_location(author_lst):
    a_time = time.time()
    a_info_list = author_lst.split("#")
    a_clean_list = []
    for author_str in a_info_list:
        try:
            author_comps = author_str.split("$")
            firstname = author_comps[0]
            lastname = author_comps[1]
            aff_raw = author_comps[2]
            aff_entity = geograpy.get_place_context(text = aff_raw)
            aff_country_cities = aff_entity.country_cities
            country_keys = list(aff_country_cities.keys())
            
            country = country_keys[0]
            city = (aff_country_cities[country])[0]
        except:
            country = "None"
            city = "None"
        clean_str = ".".join([firstname, lastname, country, city])
        a_clean_list.append(clean_str)
    print(timedelta(seconds = time.time() - a_time))
    return ",".join(a_clean_list)

# Input: raw dictionary form df of paper info
# Output: dictionary form of df with location/author info cleaned
def clean_location(df):
    authors = df.loc[:, "author"]
    pool = mp.Pool(mp.cpu_count())
    np_authors = list(authors.to_numpy())
    clean_authors = pool.map(get_location, [author_lst for author_lst in np_authors])
    pool.close()
    df.loc[:, "author"] = clean_authors
    return df

def parse_args(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:",["index="])
    except getopt.GetoptError:
        print("test.py -i <file index>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("test.py -i <file index>")
            sys.exit()
        elif opt in ("-i", "--index"):
            index = arg
    return index

def main(argv):
    file_index = parse_args(argv)
    in_path = "split_journals/2019_journals_" + str(file_index) + ".csv"
    out_path = "split_journals_clean/2019_journals_clean_" + str(file_index) + ".csv"
    print(in_path)
    print(out_path)

    print("read initial data")
    df = pd.read_csv(in_path, index_col = 0)
    print("cleaning locations")
    start_time = time.time()
    clean_df_dict = clean_location(df)
    clean_df = pd.DataFrame.from_dict(clean_df_dict)
    print("--- %s ---" % (timedelta(seconds = time.time() - start_time)))
    clean_df.to_csv(out_path)

if __name__ == "__main__":
    main(sys.argv[1:])