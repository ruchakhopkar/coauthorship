# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 11:03:45 2021
@author: rucha
"""

from pymed import PubMed
import pandas as pd
pubmed = PubMed(tool="PubMedSearcher", email="ruchakhopkar@gmail.com")

## PUT YOUR SEARCH TERM HERE ##

search_term ='(("2020"[Date - Create])) AND (nature OR science OR cell OR plos)'
results = pubmed.query(search_term, max_results=20000)
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
        # Append article info to dictionary 
        articleInfo.append({u'title':article['title'],
                            u'journal':article['journal'],
                            u'publication_date':article['publication_date'], 
                            u'authors':article['authors']})
    except:
        pass

# Generate Pandas DataFrame from list of dictionaries
articlesPD = pd.DataFrame.from_dict(articleInfo)
articlesPD.to_csv (r'E://PTRC//coauthorship//export_dataframe.csv', index = None, header=True) 

#Print first 10 rows of dataframe
print(articlesPD.head(10))