# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 02:36:55 2021

@author: rucha
"""

from Bio import Entrez
import csv

def search(query):
    Entrez.email = 'ruchakhopkar@gmail.com'
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax='20000',
                            retmode='xml', 
                            term=query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'ruchakhopkar@gmail.com'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results

if __name__ == '__main__':
    results = search('cold')
    id_list = results['IdList']
    chunk_size = 50  
    filename='E://PTRC//coauthorship//Book1.csv'
    xyz=[]
    #fields = ['MedlineCitation', 'PubmedData']    
    
    for chunk_i in range(0, len(id_list), chunk_size):
        chunk = id_list[chunk_i:chunk_i + chunk_size]
        try:  
                
                papers = fetch_details(chunk)
                #writer.writeheader() 
                for i, paper in enumerate(papers['PubmedArticle']):
                    print(i)
                    journal=paper['MedlineCitation']['Article']['Journal']['Title']
                    #csvwriter.writerow(journal) 
                    date=paper['MedlineCitation']['DateCompleted']
                    #csvwriter.writerow(date) 
                    title=paper['MedlineCitation']['Article']['ArticleTitle']
                    authorlist=paper['MedlineCitation']['Article']['AuthorList']
                    #csvwriter.writerow(fields) 
                    xyz.append(journal)
                    xyz.append(authorlist)
                    xyz.append(date)
                    xyz.append(title)
                    with open(filename, 'w', newline = '', encoding='UTF-8') as file:
                        csv_writer = csv.writer(file, delimiter=' ',dialect = 'excel')
                        
                        for i in xyz:
                            csv_writer.writerow(str(i))
                        file.close()
        except: 
            pass
        
        
    # writing the fields 
    
        import json
        print(json.dumps(papers['PubmedArticle'][0]))
