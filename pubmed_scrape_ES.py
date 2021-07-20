# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 11:03:45 2021

@author: rucha
"""

from string import digits
from pymed import PubMed
import pandas as pd

def extract_pubmed_data(pubmed, search_term):
    results = pubmed.query(search_term, max_results=60000)
    articleList = []
    articleInfo = []

    for article in results:
    # Print the type of object we've found (can be either PubMedBookArticle or PubMedArticle).
    # We need to convert it to dictionary with available function
        articleDict = article.toDict()
        articleList.append(articleDict)

    # Generate list of dict records which will hold all article details that could be fetch from PUBMED API
    for article in articleList:
        try:
            #Sometimes article['pubmed_id'] contains list separated with comma - take first pubmedId in that list - thats article pubmedId
            pubmedId = article['pubmed_id'].partition('\n')[0]

            author_info = []
            for a in article['authors']:
                a_dict = {}
                a_dict['lastname'] = a['lastname']
                a_dict['firstname'] = a['firstname']
                a_dict['initials'] = a['initials']
                affiliation_split = (a['affiliation']).split(",")
                a_dict['country'] = (affiliation_split[-1].replace('.', '')).lower()
                try:
                    a_dict['country'] = (((a_dict['country'].strip()).split(" "))[0]).strip()
                except:
                    a_dict['country'] = a_dict['country'].strip()                
                if affiliation_split[-2].isnumeric():
                    remove_digits = str.maketrans('', '', digits)
                    a_dict['city'] = (affiliation_split[-3].translate(remove_digits)).lower()
                    try:
                        a_dict['city'] = ((a_dict['city'].strip().split(" "))[0]).strip()
                    except:
                        a_dict['city'] = a_dict['city'].strip()
                else:
                    remove_digits = str.maketrans('', '', digits)
                    a_dict['city'] = (affiliation_split[-2].translate(remove_digits)).lower()
                    try:
                        a_dict['city'] = ((a_dict['city'].strip().split(" "))[0]).strip()
                    except:
                        a_dict['city'] = a_dict['city'].strip()
                author_info.append(".".join(a_dict.values()))
            if (len(author_info) > 0):
                # Append article info to dictionary 
                articleInfo.append({u'title':article['title'],
                                    u'journal':article['journal'],
                                    u'publication_date':str(article['publication_date']), 
                                    u'authors':author_info})
        except:
            pass
    return articleInfo

def save_to_df(articleInfo):
    # Generate Pandas DataFrame from list of dictionaries
    articlesPD = pd.DataFrame.from_dict(articleInfo)
    #Print first 10 rows of dataframe
    print(articlesPD.head(10))
    articlesPD.to_csv (r'pubmed_info_2019.csv', index = None, header=True)

def scrape_pubmed():
    pubmed = PubMed(tool="PubMedSearcher", email="persianblue000@gmail.com")
    search_term = '(((("nature*"[Journal]) OR ("science*"[Journal])) OR ("cell*"[Journal])) OR ("plos*"[Journal])) AND (("2019/01/01"[Date - Publication] : "2019/12/31"[Date - Publication]))'
    search_term2 = '(((("nature*") OR ("science*")) OR ("cell*")) OR ("plos*")) AND (("2019/01/01"[Date - Publication] : "2019/12/31"[Date - Publication]))'
    articleInfo = extract_pubmed_data(pubmed, search_term2)
    save_to_df(articleInfo) 

def make_graph(path):
    pubmed_df = pd.read_csv("pubmed_info_2019.csv")
    authors = []
    for index, row in pubmed_df.iterrows():
        a = (row['authors'][1:-1]).replace("'", "")
        a = a.replace(" ", "")
        a = a.split(",")
        authors += a
    authors_set = set(authors)
    print(len(authors))
    print(len(authors_set))
    print(authors[:20])

if __name__ == '__main__':
    scrape_pubmed()


