import numpy as np
import pandas as pd
from Bio import Entrez

def search(query, min_date, max_date, ret_start):
    Entrez.email = 'persianblue000@gmail.com'
    handle = Entrez.esearch(db='pubmed',
                            sort ='Pub Date',
                            datetype = "PDAT", 
                            mindate = min_date,
                            maxdate = max_date,
                            retstart = str(ret_start), 
                            # retmax = str(ret_max),
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


def fetch_data(query, min_date, max_date, ret_start):
    results = search(query, min_date, max_date, ret_start)
    id_list = results['IdList']
    papers = fetch_details(id_list)
    record = Entrez.parse(papers)

    paper_df = pd.DataFrame(columns = ["pub_year", "pub_month", \
            "journal_name", "paper_title", "country", "city", \
            "first_name", "last_name"])

    for i,paper in enumerate(papers['PubmedArticle']):
        paper_title = (paper['MedlineCitation']['Article']['ArticleTitle'])[:-1]
        journal_name = paper['MedlineCitation']['Article']['Journal']['Title']
        pub_year = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
        pub_month = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Month']

        author_info = []
        if ('AuthorList' in paper['MedlineCitation']['Article']):
            authors = paper['MedlineCitation']['Article']['AuthorList']
        else:
            authors = []
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
            a_info = {"pub_year":pub_year, "pub_month":pub_month, \
                       "journal_name":journal_name, "paper_title":paper_title, \
                       "country":country, "city":city, \
                       "first_name":first_name, "last_name":last_name}
            author_info.append(a_info)

        for a in author_info:
            paper_df = paper_df.append(a, ignore_index = True)
    return paper_df

if __name__ == '__main__':
    query = "Nature[journal] OR Science[journal] OR Cell[journal] OR PLOS[journal]"
    query1 = "Nature[journal]"
    query2 = "Science[journal"
    min_date = 2019
    max_date = 2019
    ret_start = 1
    ret_max = 100
    paper_df = fetch_data(query2, min_date, max_date, ret_start)
    paper_df.to_csv("2019_Nature_test.csv")
