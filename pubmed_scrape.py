import numpy as np
import pandas as pd
from Bio import Entrez

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
            if (query not in journal_name.lower()):
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
        else:
            authors = []
        for j in range(len(authors)):
            if (len(authors[j]['AffiliationInfo']) > 0):
                affiliation = authors[j]['AffiliationInfo'][0]['Affiliation']
                print(affiliation)
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

if __name__ == '__main__':
    query = "Nature*[journal] OR Science*[journal] OR Cell*[journal] OR PLOS*[journal]"
    query1 = 'nature'
    query2 = "Science[journal"
    min_date = 2019
    max_date = 2019
    ret_start = 0
    ret_max = 10000
    nature_df = paper_loop(query1, min_date, max_date, 150000, ret_max)
    nature_df.to_csv("2019_Nature_all.csv")
