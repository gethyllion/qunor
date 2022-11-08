


python_path =['/home/zhuj/miniconda3/envs/nova/bin/python','/home/zhuj/.conda/envs/nova/bin/python']
batch_gsea_path =["/home/zhuj/miniconda3/envs/nova/lib/python3.10/site-packages/qunor/batch_gsea.py",'/home/zhuj/.conda/envs/nova/lib/python3.10/site-packages/qunor/batch_gsea.py']
gmx_path =['/home/hejx/gmx2021.3/bin/gmx','/home/zhuj/md/gmx2021.3/bin/gmx']
vina_path = ['/home/lixingyu/autodock_vina/autodock_vina_1_1_2_linux_x86/bin/vina','/home/lixingyu/autodock_vina/autodock_vina_1_1_2_linux_x86/bin/vina']
obabel_path =['/home/zhuj/miniconda3/bin/obabel','/home/lixingyu/miniconda/miniconda3/bin/obabel']

import os
def get_available_path( path_list ):
    for i in path_list:
        if os.path.exists(i):
            return i
  
        elif path_list[-1] == i and not os.path.exists(i):
            #return None
            raise FileExistsError('no available file in path list')

def decide_which_python(python_path = python_path,):
    
    use_python = get_available_path(python_path)
    return use_python

def decide_which_batch_gsea(batch_gsea_path = batch_gsea_path):

    use_batch = get_available_path(batch_gsea_path)
    return use_batch

def decide_which_obabel(obabel_path = obabel_path):
    use_babel = get_available_path(obabel_path)
    return use_babel

def decide_which_vina(vina_path = vina_path):
    use_vina = get_available_path(vina_path)
    return use_vina

    
from typing import Iterable

class NotIterableError(BaseException):
    pass


def compare_gene_with_path(genes_q, path,log_tumor,thresh=.05, *argus, **kwargus):
    
    genes_satisfied = []

    for i in genes_q:
        try:
            corr = log_tumor.loc[[i]+path,:].T.corr().loc[i]
            if all(corr>thresh):
                genes_satisfied.append(i)
                

            else:
                if kwargus['silence'] == True:
                    pass
                else:
                    print('NONE gene in genes_q ')
        except Exception as e:
            print(i,e)
    print(len(genes_satisfied))
    return genes_satisfied

def gname(a,df=None):
    v = list(df.columns[df.columns.str.contains(a)])
    return v

def iterable_to_txt(iterable_object, file_name = None, **kwargus):
    
    if isinstance(iterable_object,Iterable):

        if 'writetxt' in kwargus.keys():
            if kwargus['writetxt'] == 1:
                
                if hasattr(iterable_object,'name') and not file_name:
                    file_name = iterable_object.name+'.txt'
                    with open(file_name,'w') as f:
                        for i in iterable_object:
                            f.write(i+'\n')
                
                elif not hasattr(iterable_object,'name') and file_name:
                    file_name = file_name +'.txt'
                    with open(file_name, 'w') as f:
                        for i in iterable_object:
                            f.write(i+'\n')
                    pass
            
                else:
                    import time
                    file_name = str(hash(time.ctime())) +'.txt'
                    with open(file_name, 'w') as f:
                        for i in iterable_object:
                            f.write(i+'\n')
                return file_name
        else:
            file_name = str(file_name)
            with open(file_name,'w') as f:
                for i in iterable_object:
                    f.write(i+'\n')
            return file_name
       

    else:
        raise NotIterableError
        
if __name__ =='__main__':
    pass