import numpy as np
import pandas as pd
import time
import sys, getopt
import itertools

'''
Co-Authorship Network Project
Eunseo Sung
Filtering authors for network
    - Authors with cleaned affiliation info
    - Filter authors with <2 papers
    - Make author edge list
'''




class Author:
    def __init__(self, author_str):
        split_lst = author_str.split(".")
        if len(split_lst) == 2:
            f_name = split_lst[0]
            l_name = split_lst[1]
            country = None
            city = None
        elif len(split_lst) == 3:
            f_name = split_lst[0]
            l_name = split_lst[1]
            country = split_lst[2]
            city = None
        elif len(split_lst) < 2:
            f_name = None
            l_name = None
            country = None
            city = None
        else:
            f_name = split_lst[0]
            l_name = split_lst[1]
            country = split_lst[2]
            city = split_lst[3]

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

    def __str__(self):
        return self.a_str

    def rank_author(self):
        if (self.country == None) and (self.city == None):
            rank = 1
        elif (self.country == None) and (self.city != None):
            rank = 2
        elif (self.country != None) and (self.city == None):
            rank = 3
        else:
            rank = 4
        return rank
        

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


# Returns list of grouped authors and just list of authors
# Input: Column of authors
# Output: 2D list of authors grouped by paper, list of all authors (with overlap),
#         list of length of all authors in order
def get_author_lists(author_col):
    authors_length = []
    authors_all = []
    for group in author_col:
        a_list = group.split(",")
        authors_all += a_list
        authors_length.append(len(a_list))
    return authors_all, authors_length


# Making sublists to group equivalent authors
# Input: list of all authors w/ overlap
# Ouput: Sublists of equivalent author object and author strings (list of np arrays)
def get_author_sublists(authors_all):
    a_obj_list = []
    a_rank_list = []
    for a in authors_all:
        a_obj = Author(a)
        a_obj_list.append(a_obj)
        a_rank_list.append(a_obj.rank_author())
    a_obj_array = np.array(a_obj_list)
    a_rank_array = np.array(a_rank_list)   
    authors_all_array = np.array(authors_all)

    author_str_sublists = []
    while(len(a_obj_list) > 0):
        a_curr = a_obj_array[0]
        eq_index = np.where(a_obj_array == a_curr)
        rank_order = np.argsort(a_rank_array[eq_index])
        author_str_sublists.append((np.array(list(authors_all_array[i] for i in eq_index)))[rank_order])
        a_obj_list = np.delete(a_obj_list, eq_index)
        authors_all_array = np.delete(authors_all_array, eq_index)
        a_rank_array = np.delete(a_rank_array, eq_index)
    return author_str_sublists
    

# Input: author_str_sublists and authors_all
# Output: authors_all where all equivalent authors are replaced by main author
def replace_eq_authors(author_str_sublists, authors_all):
    authors_all_array = np.array(authors_all)
    for sublist in author_str_sublists:
        for i in range(1, len(sublist)):
            authors_all_array = np.where(authors_all_array == sublist[i], sublist[0], authors_all_array)
    return authors_all_array


# Input: authors_all_array where only main author strings are left
#        authors_length with split length for each group
# Output: author strings split by correct paper groupings
def split_author_groups(authors_all_array, authors_length):
    authors_group_fin = []
    cur_ind = 0
    for i in range(len(authors_length)):
        cur_sublist = authors_all_array[cur_ind : cur_ind + authors_length[i]]
        authors_group_fin.append(cur_sublist)
        cur_ind += authors_length[i]
    return authors_group_fin


# Input: authors_group_fin
# Ouput: make edge list for authors
def author_graph_edges(authors_group_fin):
    edges = []
    for group in authors_group_fin:
        edges += list(itertools.combinations(group, 2))
    edges_unique = set(edges)
    edges_final = []
    for e in edges_unique:
        e_weight = edges.count(e)
        e_lst_weight = tuple([e[0], e[1], e_weight])
        edges_final.append(e_lst_weight)
    return edges_final


def main():
    file_index = 0
    in_path = "2019_split_journals_clean/2019_journals_clean_" + str(file_index) + ".csv"
    author_df = pd.read_csv(in_path, index_col = 0)
    author_col = list(author_df.loc[0:1000, "author"])
    authors_all, authors_length = get_author_lists(author_col)
    print("Got all author lists")
    author_str_sublists = get_author_sublists(authors_all)
    print("Got all author sublists")
    authors_all_array = replace_eq_authors(author_str_sublists, authors_all)
    print("Equivalent authors replaced")
    authors_group_fin = split_author_groups(authors_all_array, authors_length)
    print("Authors placed in final groups")
    edges_final = author_graph_edges(authors_group_fin)
    print("Edge list made")
    edge_df = pd.DataFrame(edges_final, columns = ["author1", "author2", "weight"])
    edge_df.to_csv("2020_test_edges.csv", index = False)
    return

if __name__ == "__main__":
    # A1 = Author("A", "B", "C", "D")
    # A2 = Author("A", "B", "R", "K")
    # A3 = Author("J", "E", "W", "S")
    # l1 = [A1, A2]
    # l2 = [A1, A3]
    # print("test1")
    # print(A1 in l1)
    # print(A2 in l1)
    # print(A3 in l1)
    # print("test2")
    # print(A1 in l2)
    # print(A2 in l2)
    # print(A3 in l2)
    main()