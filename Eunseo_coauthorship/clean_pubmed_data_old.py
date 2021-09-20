import numpy as np
import pandas as pd
import geograpy

'''
Co-Authorship Network Project
Eunseo Sung
Cleaning journal info from pubmed
    - Countries and cities/states
    - Filter authors with <2 papers
'''

# Input: author string (firstname$lastname$affiliation)
# Output: string firstname.lastname.country.city
def get_location(author_str):
    try:
        author_comps = author_str.split("$")
        firstname = author_comps[0]
        lastname = author_comps[1]
        aff_raw = author_comps[2]
        aff_entity = geograpy.get_place_context(text = aff_raw)
        aff_country_cities = aff_entity.country_cities
        country_keys = list(aff_country_cities.keys())
    except:
        return None
   
    try:
        country = country_keys[0]
        city = (aff_country_cities[country])[0]
    except:
        return None 
    clean_str = ".".join([firstname, lastname, country, city])
    return clean_str

# Input: raw dictionary form df of paper info
# Output: dictionary form of df with location/author info cleaned
def clean_location(df_dict):
    for a_key in df_dict["author"]:
        a_info_list = (df_dict["author"][a_key]).split("#")
        a_clean_list = []
        for author_info in a_info_list:
            author_loc = get_location(author_info)
            if author_loc != None:
                a_clean_list.append(author_loc)
        df_dict["author"][a_key] = a_clean_list
    return df_dict

# Input: dictionary form of df with cleaned location/author info
# Output: set of authors with > 2 papers
def get_valid_authors(df_dict, authors_dest):
    author_counts = {}
    for a_key in df_dict["author"]:
        a_info_list = df_dict["author"][a_key]
        for a in a_info_list:
            if a in author_counts:
                author_counts[a] += 1
            else:
                author_counts[a] = 1
    valid_authors = set([a for a, v in author_counts.items() if v >= 2])
    author_str = ("\n".join(list(valid_authors))).strip()
    with open(authors_dest, "wt") as f:
        f.write(author_str)
    return author_counts

def main(path):
    print("read initial data")
    df = pd.read_csv(path, index_col = 0)
    df_dict = df.to_dict()
    # test_df = df.head(10)
    # test_df_dict = test_df.to_dict()
    print("cleaning locations")
    clean_df_dict = clean_location(df_dict)
    clean_df = pd.DataFrame.from_dict(clean_df_dict)
    clean_df.to_csv("2019_all_journal_clean.csv")

    # clean_df_dict = clean_location(test_df_dict)
    # clean_df = pd.DataFrame.from_dict(clean_df_dict)
    # clean_df.to_csv("2019_all_journal_clean.csv")
    print("getting authors")
    get_valid_authors(clean_df_dict, "2019_all_clean_authors.txt")

if __name__ == "__main__":
    main("2019_all_journal.csv")
    # get_location("FangXiao$Xu$Department of Neurobiology and Department of Neurology of First Affiliated Hospital, NHC and CAMS Key Laboratory of Medical Neurobiology, Zhejiang University School of Medicine, Hangzhou, China.")