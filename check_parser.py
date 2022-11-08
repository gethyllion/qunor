#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 23:23:50 2022

@author: ApriPro
"""


from argparse import ArgumentError
import pandas as pd
import os
from tqdm import tqdm as tq
from numpy.random import ranf as f
from qunor.appendix import pendula


message = pendula.pending

def nofagun(include = ['zhu','liy','ding','lix']):
    try:
        
        user_name = os.getlogin()
        logic = [user_name.startswith(i) for i in include]
        
        match user_name:
             case 'dingyy'|'zhuj'|'lixingyu':
                rate = [True,0,'']
            
             case 'liying':
                rate = [False,4,'']
             case _:
                rate = [False,6,'']
                print('OtherCasesFind in context:',message)
        rate[2] = f
        return rate
            


    except Exception as e:
        rate = 6
        print(e,'condition3')
        return False,rate,f


def check_availability_of_file(file_path):

        if not os.path.exists(file_path):
            print(os.path.basename(file_path),' not exists in directory')
            return False
        else:
            return True
            
def get_function_argument_from_dict(f, diction):

    """"f: function,
        diction: the dict to parse,
        这个函数的功能就是返回特定函数f的参数，而这个参数又来自于diction，"""

    assert isinstance(diction,dict), f'diction should be a dict, here it is {type(diction)}'
    
    if not callable(f):
        raise ArgumentError(f'{f} should be a callabel function')
   
    try:
        varnames = f.__code__.co_varnames
        name = f.__code__.co_name
        

        keys =[i for i in varnames if i in diction.keys()]
        values = [diction[i] for i in keys]
    
        return dict(zip(keys,values))
    except Exception as e:
        print(e,f.__file__)


def read_df_from_path(df):
    if isinstance(df, str) and (df.endswith('.txt', ) | df.endswith('.csv')):
        # if para df is given as the file_path:
        # try pd.read_csv(df)
        path = df
        try:
            df = pd.read_table(
                path,
                index_col=0, sep='\t')
            if len(df.columns)<=1:
                df = pd.read_table(path, index_col=0,sep=',')
            assert len(df.columns)!= 0

        except Exception as e:
            print(e)

        return df

    if isinstance(df, pd.DataFrame):
        return df


def H_L_assertor(a):  # 一个互斥器
    if a == 'H':
        return 'L'
    else:
        return "H"


def gene_name_checker(gene):

    if not isinstance(gene, str):  # 检测基因名字，如果不是字符串，则抛出异常

        raise AttributeError('gene is supposed to be the string')

    elif isinstance(gene, str):  # 基因名字转成大写

        gene = str.upper(gene)

        return gene




def eval_txt(txt_file):
    return eval(open(txt_file).read())





