
import argparse as agp
from qunor import check_parser
import os
from typing import Iterable
description =f"""Hola, {os.getlogin()}
"""

parser = agp.ArgumentParser(description = description)
parser.add_argument('-n', '--normalize_path', type = str, help = 'normalizeExp.txt')
parser.add_argument('-t', '--task_q_path', default='task.txt', type = str, help = 'the task of gene queque to be applied gsea, with linebreak you shall join your queque in a txt')
parser.add_argument('-r','--res_dict',default='res_dict', type=str, help ='the target directory for results from gsea')
parser.add_argument('-gmt','--gmt',default='1', type=str,help=
""" 1, c2.cp.kegg.v7.5.1.symbols.gmt  
2, c2.all.v7.5.1.symbols.gmt(too tedious to go), 
any other values would be considered as file path, programe then would use a given gmt path""")
parser.add_argument('-ignore','--ignore',default=0, type=int,help=' 如果你的GEO数据已经完全选择了肿瘤样本，请选择1，如果是TCGA数据，请选择0或者不作此参数的调整。ingnore tumor/normal filtering, choose 0 when it is TCGA ')
parser.add_argument('-fdrthresh','--fdrthresh',default=0.25, type=float,help='')

cmd_argus = parser.parse_args()

normalize = cmd_argus.normalize_path # expression dataframe

task = cmd_argus.task_q_path # genes queque
fdr_threshold =cmd_argus.fdrthresh

def read_gene_task_queque_txt(file:str, sep ='\n'):
    """read a file which is a txt,sep='\n'
    """
    assert os.path.exists(file),'task_q file not exists' #check if file path exists
    sep = sep
    with open(file) as f:
        task = f.read()

    try:
        task = eval(task) # if the file is a repr-ed text
    
    except Exception as e:
        print('intended to parse python object with eval() while encountering problems, attempting split the file')
        try:

            task = task.split( sep)

            assert len(task)!=0

            task = [i for i in task if i !=''] # just in case of '' in list, so we remove all '' in task list, blank quotation happens ususally in last line.
            
        except Exception as e:
            raise AttributeError('could not parse task files')

    assert isinstance(task,Iterable),'task should be iterable'
    return task

if __name__ =='__main__':


    from qunor import gseatool as gst
    from qunor import check_parser
    
    import pandas as pd
    import os,time
    
    print(f'start at: {time.ctime()}')
    print('current_wd: ',os.getcwd())
    print('gseatool file in: ', gst.__file__,)
    

    

    
    
    dir_name = cmd_argus.res_dict
    normalize = check_parser.read_df_from_path(normalize)
    gmt =cmd_argus.gmt
    res_dict = cmd_argus.res_dict

    ignore = cmd_argus.ignore

    task = cmd_argus.task_q_path
    task =read_gene_task_queque_txt(task)
    print(
    'total genes number in tumor_normalize: ', len(normalize.index))
    gene_sets1 ="/home/zhuj/data_dicting/c2.cp.kegg.v7.5.1.symbols.gmt"
    gene_sets2 ='/home/zhuj/data_dicting/curated/c2.all.v7.5.1.symbols.gmt'

    
    try:
        dir_name = dir_name
        os.mkdir(dir_name)

    except Exception as e:
        
        print(e)

    if gmt == '1':
        gene_sets = gene_sets1

    elif gmt == '2':
        gene_sets = gene_sets2

    else:
        assert os.path.exists(gmt),f'file {gmt} not exists'
        gene_sets = gmt

    print(f'your are using the following gmt: {gmt}',)
    res = gst.gsea_batch(task_q = task, data = normalize,dir_name = dir_name
    ,permutation_num = 1600,seed = 8, fdr_threshhold = fdr_threshold,permutation_type='phenotype'
    ,method='signal_to_noise', processes = 4
    , no_plot = True, outdir = None, gene_sets = gene_sets, ignore = ignore)
    res = repr(res)
  
  
    with open('results_gsea.python_repr','w') as f:
        
        f.write(res)
        
    print('saved —— results_gsea.txt ',)
    print(time.ctime())
    try:
        from qunor import func
        df = func.make_df(f'{res_dict}')
        df.to_csv(f'{len(task)}genes_gsea_dataframe.txt')
    except Exception as e:
        print(e,)
