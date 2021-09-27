import numpy as np
import pandas as pd
import time
import sys, getopt

'''
Co-Authorship Network Project
Eunseo Sung
Filtering authors for network
    - Authors with cleaned affiliation info
    - Filter authors with <2 papers
    - Make author edge list
'''


# Make dictionary with all edges

# Extract dictionary keys (all unique authors)

# Store authors as class

# Combine authors that are "equal"

# Combine keys that are equal in dictionary

# Store to dataframe



class Author:
    def __init__(self, author_str):
        f_name, l_name, country, city = author_str.split(".")
        self.a_str = author_str
        self.f_name = f_name
        self.l_name = l_name
        if country == "None":
            self.country = None
        else:
            self.country = country
        if city == "None":
            self.city = None
        else:
            self.city = city
    
    def __eq__(self, other):
        if (self.f_name == other.f_name) and (self.l_name == other.l_name):
            if (self.country == other.country) and (self.city == other.city):
                return True
            elif (self.country == None) and (self.city == None):
                return True
            elif (self.country == None) and (self.city == other.city):
                return True
            elif (self.country == other.country) and (self.city == None):
                return True
            elif (other.country == None) and (other.city == None):
                return True
            elif (other.country == None) and (other.city == self.city):
                return True
            elif (other.country == self.country) and (other.city == None):
                return True
        return False

# get list of authors
# Input: Author column of cleaned csv file (all authors)
# Output: Dictionary (key = main author, value = equivalent authors)
def list_authors(authors):
    auth_list = []
    for paper in authors:
        a_list = paper.split(",")
        a_obj_list = []
        for a in a_list:
            a_obj_list.append(Author(a))
        auth_list.append(a_obj_list)
    return auth_list

# Input: Author column of cleaned csv file
# Output: Dictionary of author edges
def get_author_dict(auth_list):
    a_dict = dict()
    for a_paper in auth_list:
        for a in a_paper:
            return
            

def main():
    file_index = 0
    in_path = "split_journals_clean/2019_journals_clean_" + str(file_index) + ".csv"
    author_df = pd.read_csv(in_path, index_col = 0)
    author_col = list(author_df.loc[:, "author"])
    get_author_dict(author_col)
    return

if __name__ == "__main__":
    A1 = Author("A", "B", "C", "D")
    A2 = Author("A", "B", "R", "K")
    A3 = Author("J", "E", "W", "S")
    l1 = [A1, A2]
    l2 = [A1, A3]
    print("test1")
    print(A1 in l1)
    print(A2 in l1)
    print(A3 in l1)
    print("test2")
    print(A1 in l2)
    print(A2 in l2)
    print(A3 in l2)
    