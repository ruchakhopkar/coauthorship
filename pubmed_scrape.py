import numpy as np
import sys, getopt
import pandas as pd
import re
from Bio import Entrez
import scipy.sparse as sparse

def search(query, min_date, max_date, ret_start, ret_max):
    Entrez.email = 'persianblue000@gmail.com'
    handle = Entrez.esearch(db='pubmed',
                            sort ='Pub Date',
                            datetype = "PDAT", 
                            mindate = min_date,
                            maxdate = max_date,
                            retstart = str(ret_start), 
                            retmax = str(ret_max),
                            retmode = 'xml',
                            term = query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'persianblue000@gmail.com'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results


def fetch_data(query, min_date, max_date, ret_start, ret_max):
    results = search(query, min_date, max_date, ret_start, ret_max)
    id_list = results['IdList']
    papers = fetch_details(id_list)
    record = Entrez.parse(papers)

    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", \
        "Sep", "Oct", "Nov", "Dec"]

    paper_df = pd.DataFrame(columns = ["pub_year", "pub_month", \
            "journal_name", "paper_title", "author"])

    for i,paper in enumerate(papers['PubmedArticle']):
        try:
            paper_title = (paper['MedlineCitation']['Article']['ArticleTitle'])[:-1]
        except:
            continue
        try:
            journal_name = paper['MedlineCitation']['Article']['Journal']['Title']
            journal_check = bool(re.match(query + "*", journal_name.lower()))
            if journal_check == False:
                continue
        except:
            continue
        try:
            pub_year = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
        except:
            continue
        try:
            pub_month = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Month']
            if pub_month in months:
                pub_month = str(months.index(pub_month))
        except:
            continue

        author_info = []
        if ('AuthorList' in paper['MedlineCitation']['Article']):
            authors = paper['MedlineCitation']['Article']['AuthorList']
            if len(authors) < 2:
                continue
        else:
            continue
        for j in range(len(authors)):
            if (len(authors[j]['AffiliationInfo']) > 0):
                affiliation = authors[j]['AffiliationInfo'][0]['Affiliation']
                try:
                    city = (affiliation.split(",")[-2]).strip()
                    country = ((affiliation.split(",")[-1]).split(".")[0]).strip()
                except:
                    city = None
                    country = None
            else:
                affiliation = None
                city = None
                country = None
            if ('LastName' in authors[j]) and ('ForeName' in authors[j]):
                last_name = authors[j]['LastName']
                first_name = authors[j]['ForeName']
            else:
                last_name = None
                first_name = None
            a_str = ".".join([str(country), str(city), str(first_name), str(last_name)])
            author_info.append(a_str)

        paper_info = {"pub_year":pub_year, "pub_month":pub_month, \
                    "journal_name":journal_name, "paper_title":paper_title, \
                    "author":",".join(author_info)}

        paper_df = paper_df.append(paper_info, ignore_index = True)
    return paper_df

def paper_loop(query, min_date, max_date, est_max, ret_max):
    total_df = pd.DataFrame(columns = ["pub_year", "pub_month", \
                "journal_name", "paper_title", "author"])
    num_loop = int(est_max / ret_max)
    start_i = 0
    for i in range(num_loop):
        print(i)
        try:
            paper_df = fetch_data(query, min_date, max_date, start_i, ret_max)
            total_df = total_df.append(paper_df, ignore_index = True)
            start_i += ret_max
        except:
            break
    return total_df

def combine_df(paths, dest_path):
    combined_df = pd.concat(paths)
    combined_df.to_csv(dest_path)
    return combined_df

def make_graph(df, uw_path, w_path, authors_path):
    authors_all = []
    paper_authors = []
    for i in range(len(df)):
        authors = df['author'][i].split(",")
        authors_all += authors
        paper_authors.append(authors)
    authors_set = set(authors_all)
    authors_index = dict()
    i = 0
    for a in authors_set:
        authors_index[a] = i
        i += 1

    authors_str = "\n".join(list(authors_set))
    authors_str.encode(sys.getdefaultencoding(), errors = 'replace')
    with open(authors_path, "wt", encoding = 'utf-8') as f:
        f.write(authors_str)

    net_w = sparse.lil_matrix((len(authors_set), len(authors_set)))
    for i in range(len(paper_authors)):
        for j in range(len(paper_authors[i])):
            for k in range(j, len(paper_authors[i])):
                x_ind = authors_index[paper_authors[i][j]]
                y_ind = authors_index[paper_authors[i][k]]
                net_w[x_ind, y_ind] += 1
                net_w[y_ind, x_ind] += 1
    
    net_uw = net_w.copy()
    net_uw[net_uw.nonzero()] = 1

    print("All Authors: {}".format(len(authors_all)))
    print("Num Nodes: {}".format(len(authors_set)))
    print("Num Edges: {}".format(net_uw.count_nonzero()))


    
def search_pubmed(queries, min_date, max_date, ret_max, save_sep):
    max_count = 3000000
    journal_dfs = []
    for query in queries:
        journal_dfs.append(paper_loop(query, min_date, max_date, max_count, ret_max))
    if save_sep:
        for i in range(len(queries)):
            filename = str(min_date) + "_" + queries[i] + ".csv"
            journal_dfs[i].to_csv(filename)
    return journal_dfs

def main(search, save_sep, save_comb):
    queries = ['nature', 'science', 'cell', 'plos']
    min_date = 2019
    max_date = 2019
    ret_max = 10000
    if search:
        search_pubmed(queries, min_date, max_date, ret_max, save_sep)
    journal_paths = []
    for q in queries:
        journal_paths.append(str(min_date) + "_" + q + ".csv")

    if save_comb:
        df_comb = combine_df(journal_paths, str(min_date) + "_all.csv")
    else:
        df_comb = pd.read_csv(str(min_date) + "_all.csv")
    
    make_graph(df_comb, str(min_date) + "_all_unweighted_network.csv", \
                        str(min_date) + "_all_weighted_network.csv", \
                        str(min_date) + "_authors.txt")


if __name__ == '__main__':
    main(False, False, False)

"""
keywords: Nature, Cell, Science, PLoS
year: 2019
Num Nodes (authors): 140023
Num Edges (connections): 2778551
Num Authors: 150387
Num Overlapping authors: 10364
"""