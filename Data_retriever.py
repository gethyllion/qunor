# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 09:58:35 2022

@author: zhu jun
"""

import pandas as pd
import re
import mygene
import requests
from bs4 import BeautifulSoup
from tqdm import tqdm as tq
import time


class Retriever:
    headers = {'User-Agent': 'qunor onda',
               'Accept-Encoding': 'gzip, deflate',
               'Accept': '*/*',
               'Connection': 'keep-alive'}

    def __init__(self, q):
        self.q = q.split(' ')
        self.pmid = None
        self.meta = []
        self.total_page = None
        self.big_pages = None
        self.batch = None

    def stat(self, a, **kwargu):

        #         if contains:
        #             contains = contains.split(' ')

        condition = '|'.join(self.q)

        result = set(re.findall(pattern=condition, string=a))

        return result

    def get_element(self, path):
        try:
            r = requests.get(path, headers=self.headers)
        except ConnectionError:
            time.sleep(3)

            r = requests.get(path, headers=self.headers)

        soup = BeautifulSoup(r.text, 'lxml')

        abstract = soup.find(id='enc-abstract').text.replace('\n', '')
        title = soup.find(attrs='heading-title').text.replace('\n', '')
        pmid = soup.find(attrs='current-id').text
        meta = {'pmid': pmid, 'title': title,
                'ab': abstract,
                'm_title': self.stat(a=title),
                'm_ab': self.stat(a=abstract),
                'title_match_score':len( self.stat(a=title)),
                'abstract_match_score':len(self.stat(a=abstract))}

        return meta

    def get_total_big_pages(self, ):
        term = '%20'.join(self.q)
        link = f'https://pubmed.ncbi.nlm.nih.gov/?term={term}'
        soup = BeautifulSoup(requests.get(link, self.headers).text, 'lxml')
        the_last_number = soup.find(attrs="of-total-pages").text.split(' ')[-1].replace(',', '')
        the_last_number = int(the_last_number)
        self.total_page = the_last_number

        print(f'total pages = {self.total_page}')
        
        if self.total_page >=5:
            raise AttributeError('too many pages!')

        data = []
        for i in range(1, self.total_page):
            temp = link + f'&pages={i}'
            data.append(temp)
        self.big_pages = data

    def get_pmid(self):
        self.get_total_big_pages()

        data = []
        for i in tq(self.big_pages, desc='GETTIGN BIG PAGES'):
            r = requests.get(i, )
            soup = BeautifulSoup(r.text, 'lxml').find_all(attrs="docsum-pmid")
            if soup is []:
                raise ReferenceError('no pmid')
            data.extend([i.text for i in soup])

        data2 = []
        for i in data:
            link = f"https://pubmed.ncbi.nlm.nih.gov/{i}/"
            data2.append(link)

        self.pmid = data2
        print(f'PMID TOTAL = {len(self.pmid)}')

    def get_meta(self, cotains=None):
        print(self.q)

        self.get_pmid()

        for i in tq(self.pmid, desc='GETTING PMID'):
            res = self.get_element(path=i)
            self.meta.append(res)

        self.meta = pd.DataFrame(self.meta)
        file_name =' and '.join(self.q)
        self.meta.to_csv(f'{file_name}.csv',index_label='pmid')

    def task_q(self,q, pathway):
        pass






def symbol_2_entrez(genes_symbol, scopes='symbol', fields='entrezgene', species='human', returnall=False):
    argu = {'scopes': scopes, 'fields': fields, 'species': species, 'returnall': returnall}

    mg = mygene.MyGeneInfo()
    mg = mg.querymany(genes_symbol, **argu)
    return pd.DataFrame(mg)


def cut_list_into(TheList, space):  # 返回一个迭代器

    if isinstance(TheList, str):
        raise AttributeError(' should not be a string input')

    if not isinstance(TheList, list):
        try:
            TheList = list(TheList)
        except Exception as e:
            print(e, )

    Length_list = len(TheList)

    if space > len(TheList):

        raise AttributeError(f'you virtually cannot cut a list with lenth: {Length_list}',
                             ' every {space}, make a smaller steps.')
    else:
        for i in range(0, len(TheList), space):
            yield TheList[i:i + space]


def get_symbol(gene_q, report=True):
    if isinstance(gene_q, set):  # 去重

        pass
    else:
        gene_q = list(set(gene_q))

    if len(gene_q) <= 999:

        gene_q = [i.split('.')[0] for i in gene_q]
        url = ' http://mygene.info/v3/gene'
        headers = {'content-type': 'application/x-www-form-urlencoded',
                   'species': 'human'}
        q = ','.join(gene_q)

        params = f'ids={q}&fields=name,symbol&species=human'
        res = requests.post(url, data=params, headers=headers)

        df = pd.read_json(res.text)

        return df

    else:

        data_frame_list = []
        gene_q = cut_list_into(gene_q, space=999)
        for i in tq(gene_q, desc='7. 转换成Gene Symbol', ncols=10):
            try:

                data_frame_list.append(get_symbol(i))
            except Exception as e:
                print(e)

        try:
            total_df = pd.concat(data_frame_list)
            return total_df
        except Exception as e:
            print(e, )
            return data_frame_list


def entrez_symbol_dict(dataframe):
    if not isinstance(dataframe, pd.DataFrame):
        raise AttributeError(' input is not a dataframe..')
    else:
        pass
    try:
        return dict(zip(dataframe['query'], dataframe['symbol']))
    except Exception as e:
        print(e)
