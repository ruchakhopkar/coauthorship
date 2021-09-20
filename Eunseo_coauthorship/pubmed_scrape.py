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
        except:
            continue
        try:
            pub_year = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
        except:
            continue
        try:
            pub_month = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Month']
            if pub_month in months:
                pub_month = str(months.index(pub_month) + 1)
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
            else:
                affiliation = None
            if ('LastName' in authors[j]) and ('ForeName' in authors[j]):
                last_name = re.sub('[^a-zA-Z0-9 \n\.]', '', authors[j]['LastName'])
                first_name = re.sub('[^a-zA-Z0-9 \n\.]', '', authors[j]['ForeName'])
            else:
                last_name = None
                first_name = None
            a_str = "$".join([str(first_name), str(last_name), str(affiliation)])
            author_info.append(a_str)

        paper_info = {"pub_year":pub_year, "pub_month":pub_month, \
                    "journal_name":journal_name, "paper_title":paper_title, \
                    "author":"#".join(author_info)}

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

def get_queries(path):
    with open(path, "rt") as f:
        j_list_str = f.read()
    j_list = j_list_str.splitlines()
    j_list = [i for i in j_list if i]
    query_str = ""
    for i in range(len(j_list)):
        if i == (len(j_list) - 1):
            query_str += '\"' + j_list[i] + '\"' + "[journal]"
        else:
            query_str += '\"' + j_list[i] + '\"' + "[journal] OR "
    return query_str

def search_pubmed(queries, min_date, max_date, ret_max):
    max_count = 70000
    journal_dfs = []
    for query in queries:
        journal_dfs.append(paper_loop(query, min_date, max_date, max_count, ret_max))
    return journal_dfs    

def main(query_path, save_path):
    query_str = get_queries(query_path)
    min_date = "2019/01/01"
    max_date = "2019/12/31"
    ret_max = 10000
    journal_dfs = paper_loop(query_str, min_date, max_date, 60000, ret_max)
    journal_dfs.to_csv(save_path)

if __name__ == '__main__':
    main("final_journals_v2.txt", "2019_all_journal.csv")
    # get_queries("final_journals_v2.txt")


# def make_graph(df, uw_path, w_path, authors_path):
#     authors_all = []
#     paper_authors = []
#     for i in range(len(df)):
#         authors = df['author'][i].split(",")
#         authors_all += authors
#         paper_authors.append(authors)
#     authors_set = set(authors_all)
#     authors_index = dict()
#     i = 0
#     for a in authors_set:
#         authors_index[a] = i
#         i += 1

#     authors_str = "\n".join(list(authors_set))
#     authors_str.encode(sys.getdefaultencoding(), errors = 'replace')
#     with open(authors_path, "wt", encoding = 'utf-8') as f:
#         f.write(authors_str)

#     net_w = sparse.lil_matrix((len(authors_set), len(authors_set)))
#     for i in range(len(paper_authors)):
#         for j in range(len(paper_authors[i])):
#             for k in range(j, len(paper_authors[i])):
#                 x_ind = authors_index[paper_authors[i][j]]
#                 y_ind = authors_index[paper_authors[i][k]]
#                 net_w[x_ind, y_ind] += 1
#                 net_w[y_ind, x_ind] += 1
    
#     net_uw = net_w.copy()
#     net_uw[net_uw.nonzero()] = 1

#     print("All Authors: {}".format(len(authors_all)))
#     print("Num Nodes: {}".format(len(authors_set)))
#     print("Num Edges: {}".format(net_uw.count_nonzero()))
