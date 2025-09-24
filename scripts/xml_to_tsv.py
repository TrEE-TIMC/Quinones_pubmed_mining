from bs4 import BeautifulSoup
import gzip
import pandas as pd
import os
import argparse

journals = ['Antonie van Leeuwenhoek', 'International journal of systematic and evolutionary microbiology', 'International journal of systematic bacteriology', 'Systematic and Applied Microbiology', 'Archives of Microbiology']


parser = argparse.ArgumentParser()
parser.add_argument('-x', '--xml', required=True, help="xml path")
parser.add_argument('-o', '--output', required=True, help="output path")
args = parser.parse_args()

xml_file = args.xml
print(xml_file)
output = args.output

df = pd.DataFrame(columns=['journal_title', 'year', 'month', 'title', 'abstract', 'PMID', 'xml'])

with gzip.open(xml_file, "rb") as f:
    file = f.read()
    soup = BeautifulSoup(file, "xml")
    all_articles = soup.find_all('PubmedArticle')
    length = len(all_articles)
    xml_full = os.path.split(xml_file)[1]
    xml_short = xml_full.split('.')[0]
    for idx, item in enumerate(all_articles):
        journal_title = item.find('Title').text
        if journal_title in journals:
            print(journal_title)
            pubdate = item.find('PubDate')
            medlinecitation_date = item.find('DateCompleted')
            if pubdate and pubdate.find('Year') is not None:
                year = pubdate.find('Year').text
            else:
                if medlinecitation_date and medlinecitation_date.find('Year') is not None:
                    year = medlinecitation_date.find('Year').text
                else:
                    year = ''
            if pubdate and pubdate.find('Month') is not None:
                month = pubdate.find('Month').text
            else:
                if medlinecitation_date and medlinecitation_date.find('Month') is not None:
                    month = medlinecitation_date.find('Month').text
                else:
                    month = ''
            title = item.find('ArticleTitle').text if item.find('ArticleTitle') is not None else ''
            abstract = item.find('AbstractText').text if item.find('AbstractText') is not None else ''
            pmid = item.find('PMID').text if item.find('PMID') is not None else ''
            row = {
                'journal_title': journal_title,
                'year': year,
                'month': month,
                'title': title,
                'abstract': abstract,
                'PMID' : pmid,
                'xml': xml_short
            }

            df = df.append(row, ignore_index=True)
            print(f'Appending row %s of %s' % (idx+1, length))
        else :
            continue

df.to_csv(output, sep="\t", index=False)
