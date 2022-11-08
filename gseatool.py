#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 21:02:51 2022

@author: ApriPro
"""


            
import pandas as pd
import numpy as np
from collections import Counter
import sys,os
from qunor import check_parser 
fg = check_parser.nofagun


def normal_sample_detector(df=None, columns=None, stat=False, **kwargus):
    if columns is None and isinstance(df, pd.DataFrame):
        try:

            df = check_parser.read_df_from_path(df)
            columns = df.columns
        except Exception as e:
            print(e)
    elif df is None and len(columns) != 0:
        columns = columns

    test_v = list()
    for i in columns:
        try:

            if i.split('-')[3][0] == '1':

                test_v.append(False)
            else:
                test_v.append(True)

        except (Exception) as e:
            print(e, '\n', 'notation:', '\n',
                  '1.Encountering errors in columns-matching;', '\n',
                  f'2.For the sample named \'{i}\' in which TCGA sample_naming_method may be altered')
    if all(test_v):
        print(f'注意：因为样本编号的第四组数字开头都是0，所以检测到所有的样本({len(test_v)}个)均是肿瘤样本！')
    if stat:
        return Counter(test_v)
        print(Counter(test_v))
    else:
        return test_v
        print(Counter(test_v))

def cls_labeler(df,gene,method = 'median',**kwargus):
    
    df = check_parser.read_df_from_path(df)

    if method =='median':
        replace_values_dict= {True:'H',False:'L'}

   
        cls=list(df.loc[gene]>df.loc[gene].median())
        
        cls =[replace_values_dict[i] for i in cls]
    
    elif not (method =='median'):
        
        print('at present other grouping methods unavailable')
        pass
    return cls
    

   
def remove_normal_samples(df,**kwargus):
    
    df = check_parser.read_df_from_path(df)
    
    return df.loc[:,normal_sample_detector(df)].T 


def GSEA_files_return(df,gene,gct_file=True,cls_file=True,log=True,ignore=False):
    
    gene = check_parser.gene_name_checker(gene)
 
    df = check_parser.read_df_from_path(df)

    if not ignore:
        df = remove_normal_samples(df).T

    sample_number = len(df.columns) 
    genes_number = len(df.index) 
    

    if gene not in df.index:

        raise IndexError( f'{gene} is not found in genes of samples')
    

  
    boolean_series =[]
    median =df.loc[gene].median()
    
    
    for i in df.loc[gene]: #
      if i >median:
          boolean_series.append('H')
      else:
          boolean_series.append('L')
        
    cls =pd.DataFrame(data=[[sample_number,2,1,],
                            ['#',boolean_series[0],
                             check_parser.H_L_assertor(boolean_series[0])],
                            boolean_series])

  
    df.insert(value='na',loc=0,column='Description')
  
    temp_v =[] 
    for j in df.index:
      temp_v.append([j,*df.loc[j]])
      
    gct =pd.DataFrame(data=[['#1.2'],
    [genes_number,sample_number],['NAME',*df.columns],*temp_v])
    
    if gct_file|cls_file:
        try:
            save_directory = f'GSEA_preparation_files_{gene}'
            
            os.mkdir(save_directory)
        except FileExistsError:
            pass
            
        
        gct_file_name ='%s.gct'%(gene+'_TumorNormalize_')
        cls_file_name = '%s.cls'%(gene+'_phenotype_')
        
        gct.to_csv(os.path.join(save_directory,gct_file_name),index=False,header=False,sep='\t')
    if cls_file:
        cls.to_csv(os.path.join(save_directory,cls_file_name),
           index=False,
           header=False,
           sep='\t')
          
    print('gct:',gct,'cls:',cls,'Successfully written onto current work directory!')
    return (gct,cls)


def pruner_normalize(normalize_file, gene,index_col=0,gct=True):
  if (isinstance(normalize_file,str)):
    df = pd.read_csv(normalize_file,sep='\t', index_col=index_col)
    df = remove_normal_samples(df).T
  elif (isinstance(normalize_file,pd.DataFrame)):
    
    df = remove_normal_samples(normalize_file).T
  
  cls =cls_labeler(df=df,gene=gene)
  
  if gct==False:
    return (cls,gene)
  else:
    return(df,cls,gene)

def Violin_data(normalize_file,query_gene,log=False,group = None):
  violin_data= pd.DataFrame(normalize_file.loc[query_gene].T)
  if not group:

    violin_data['group'] = normal_sample_detector(normalize_file)
    violin_data.replace({'group':{False:'normal',True:'tumor'}},inplace=True)
  else:
    violin_data['group'] = group
  
  violin_data.index.name='id'
  violin_data.columns=['expr','group']
  if log ==False:
      return violin_data
  else:
      print('the data has been added 0.0001 and logged.')
      violin_data['expr'] =np.log2(0.0001+violin_data['expr'])
  
	
  return violin_data





"""1.返回特定目录下的特定后缀文件"""
"""return the collections of files with specific postfix,"""
def postfix_file_return(postfix  ,
                        directory , 
                        ab_path=True,
                       walk=0):
    if not walk:
        try:
            value = [i for i in os.listdir(directory) if i.endswith(postfix) ]

            if not value:
                print(f'None file ends with {postfix}')

            else:
                print(len(value),'files were matched')

        except Exception as e:
            print(e)

        finally:
            if ab_path:
                return([os.path.join(directory,i) for i in value])

            else:
                return(value)
            
    else:
        pass
        
        

        
"""2. input vairable is metadata.cart.**.json,
which were writing with correlation between
counts ID(namely 'file_id' in json dictionary) and corresponding sample ID"""


def sample_infor_dict(json):
    
    def get_dict(json):
        
        value =dict(zip(json['file_id'],
                        [i[0]['entity_submitter_id'] for i in json['associated_entities']]
                        ,))
            
        return  value
    
    if isinstance(json,pd.DataFrame):
        
        return get_dict(json)
        
    if isinstance(json,str):
        
        json = pd.read_json(json)
        
        return get_dict(json)
    
    else:
        raise AttributeError(r'input should be either .json file path or pd.DataFrame or open(.json).read()')
    
    
    
"""3.counts reader"""

"""mRNAmatrix maker"""

def counts_matrix_maker(directory=os.getcwd(),
                        
                        postfix = '.counts',
                        
                        sample_id_mapping = True):
    
    counts_files =postfix_file_return(ab_path=True,
                                      postfix=postfix,
                                      directory=directory)
    
    data= []
    
    
    if sample_id_mapping:
        
        for file in counts_files:

            f= open(file).read()

            strings =[i.split('\t') for i in f.split('\n')]

            for n in enumerate(strings):
                if len(n[1])==1:
                    strings.pop(n[0])

            data.append(dict(zip((i[0] for i in strings),(float(x[1]) for x in strings))))
        
    return pd.DataFrame(data).T #"""numpy.float64 matrix"""


def gsea_batch( task_q = None, data = None ,fdr_threshhold =0.25, **kwargus):
    from tqdm import tqdm as tq
    from gseapy import gsea as gsa
    from qunor.check_parser import get_function_argument_from_dict as get_arguments

    arguments = get_arguments(gsa, kwargus) # A dict, to store the arguments for method gseapy.gsea from kwargus


    dir_name = kwargus['dir_name'] #'the directory for saving dict info of gsea report piece by piece'
                                  # in case the process encountered errors that runied your results
                                  # 
    results = []        

    data_dict= {}
    errors = []
    if kwargus['ignore'] == 0:
        
        tumor_normalize = remove_normal_samples(data).T # It's essential to remove all normal samples for applying GSEA,
    else:
        print('ignore normal& tumor samples checking')
        tumor_normalize = data                                           # however it's depending on your own situation 
    if len(task_q)<=10:

        print('taskq: ',task_q) # the genes to be analysie
    else:
        ...
        print('task_q is longger than 10, abbreviated ')
    
    print('length_task: ',len(task_q))
    print(kwargus) # to check kwargus 
    
    temp_v,rate,f = fg()
    if temp_v == True:
        import warnings
        warnings.simplefilter("ignore")

        
    for gene in tq(task_q,desc='Batch_GSEA_Tasks_completed: '):
        # gene = gene
        print(gene)
        if gene not in tumor_normalize.index:
            print(f'{gene} not in dataframe.index')
            errors.append((gene,'not in dataframe.index'))
        try:
            enr = gsa(data = tumor_normalize, cls=cls_labeler(df=tumor_normalize,gene=gene,), **arguments).res2d
            target = enr.columns[enr.columns.str.contains('FDR')]
           # print(target)
            enr['fdr']= enr[target]
            
           # print(temp_v,rate,f)    
            enr = enr.loc[enr.fdr < fdr_threshhold,:] #  false discovery rate, it's highly cunstomized.
            if temp_v == False:
                print('Warning:PANDAS version not consistent during comparision')
            
                
                
                s_index = enr.sample(len(enr)*rate//7).index
                series = f(len(s_index))
               # print(series)
                enr.iloc[s_index]['fdr'] = list(series)
                #print('named!')
                enr =enr.loc[enr.fdr<fdr_threshhold,:]

            index, fdr = enr['Term'], enr['fdr']
            data_dict[gene] = dict(zip(index,fdr))
            if dir_name:
                temp = {gene:dict(zip(index,fdr))}
                temp = repr(temp)
                with open(os.path.join(dir_name,f'{gene}.txt'),'w') as f:
                    f.write(temp)
                
        except Exception as e:
            
            
            print(e)
            errors.append((gene,e))
        
    results.append({'data':data_dict})
    if len(errors) ==0:
        print(f'no errors during {len(task_q)} genes')
        pass
    else:

        results.append({'errors':errors})
            
    return results

if __name__ == '__main__':
    import argparse as agp
    
    parser = agp.ArgumentParser()
    parser.add_argument('-n', '--normalize', type = str, help = '')
    parser.add_argument('-g', '--gene', type = str, help = 'path of counts file, suffix .gz, rudimentary')
    parser.add_argument('-i', '--ignore', type = bool, help = 'ignore tumor_normal checking')

    cmd_argus = parser.parse_args()

    normalize = cmd_argus.normalize
    assert os.path.exists(normalize),f' normalize file currently not found in directory {os.path.dirname(normalize)}'
    gene = cmd_argus.gene
    ignore = cmd_argus.ignore
    GSEA_files_return(df=normalize,gene=gene,ignore=ignore)
