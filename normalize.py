

from typing import Iterable
from . import DataBuilder, check_parser, gseatool
from qunor.gseatool import remove_normal_samples as remover

from qunor import func
import time,os
import numpy as np, pandas as pd
import pickle as pk
from tqdm import tqdm as tq

from qunor.dict import metabolism
metabolism_path = metabolism.metabolism_path


stemness_data_pickle_path ="/home/zhuj/data_dicting/stem_index/stemness.pickle" 
stemness_pickle = pk.load(open(stemness_data_pickle_path,'rb'))

stemness_note = f"""
stemness_pickle 
是一个通过机器学习预测的癌细胞干性指数（PMID: 29625051）的数据集
这个pickle内包含了2个dict，其中一个dict会map干性指数到每个样本
没有该指数的样本会被过滤,pickle地址：{stemness_data_pickle_path}
"""
samples_cancer_type= stemness_pickle['sample_souce_cancer_type']
samples_stemness = stemness_pickle['sample_stemness']

from qunor.check_parser import check_availability_of_file
class Normalize():

    
# +++++++++++
    """d 
    """
    
    def __init__(self, 
                        normalize_path='./normalizeExp.txt', 
                                deg_path = './edgerOut.xls',         
                                corr = True,
                                **kwargus,):
        

        kwargus = kwargus.copy()
        kwargus['corr'] = corr

        if check_availability_of_file(normalize_path):
      
            self.__normalize_path = normalize_path
            self.__deg_path = deg_path
            self.metabolism_path = metabolism_path
            self(**kwargus)

        else:
            raise FileNotFoundError('normalize_path not exists')


    def __call__(self, *argus, **kwargus) :


        self.__build(*argus,**kwargus)

    def __str__(self):
        
        print('Normalize object in qunor',self.raw)
    

    def __build(self, *argus,**kwargus):
        try:
            self.__log_tumor()
            self.__read_deg()


            match kwargus['corr']:
                case True:

                    self.__correlation_path()
                case False:
                    print('not do corr')
                

    
        except Exception as e:
            print(e)
            ...
    def update(self):
        temp = Normalize(**kwargus)
        return temp

    @property
    def stat(self):
        return gseatool.normal_sample_detector(columns=list(self.raw.columns), stat=True)

    @property
    def raw(self):
        return check_parser.read_df_from_path(self.__normalize_path)

  
    def stemness(self, save = False):
        
        data = self.log_tumor
        data.loc['stemness'] = data.columns.map(samples_stemness)
        original_sample_number = len(data.columns)

        data = data.dropna(axis=1)
        
        after_drop_na_columns = len(data.columns)

        cancer_type = data.columns.map(samples_cancer_type).unique()

        cancer_name  = '_'.join(i for i in cancer_type)
        
        print(f'oiginal_columns: {original_sample_number}, after_drop_na_columns: {after_drop_na_columns} ')

        
        
        if save:
            
            file_name =f'stemness_{cancer_name}.txt'
            data.to_csv(file_name)
            print(f'file saved as {file_name}') 

        self.__stemness = data
        
        return data



    def __log_tumor(self, base_v = 0.01, **kwargus):

        self.log_tumor = np.log2(base_v + self.raw)
        self.log_tumor = remover(self.log_tumor).T
        print('log2 completed,','normal deleted,','access it by self.log_tumor')
    
    def __read_deg(self, fdr_thresh =.05):
        self.deg = pd.read_table(self.__deg_path,index_col=0)
        self.deg = self.deg[self.deg['FDR']<fdr_thresh]
        
        self.up = self.deg[(self.deg.logFC>2)].index
        self.down = self.deg[(self.deg.logFC<-2)].index

    def pull_deg(self,fc):
        data = self.deg.copy()
        deg = data[abs(data.logFC)>fc]
        up,down =deg[deg.logFC>fc].index,deg[deg.logFC<fc].index
        up,down = set(up),set(down)
        print(f'pulling degs with fc =abs({fc})')
        return {'up':up,'down':down}
    def __correlation_path(self, **kwargus):
        



        from qunor import tool4

        deg = self.up

        log_tumor = self.log_tumor

        q_keys = self.metabolism_path.keys()
        q_values = []
        
        for key in tq(q_keys):
            
            path = self.metabolism_path[key]

            value = tool4.compare_gene_with_path(deg, path, log_tumor,silence=True,**kwargus)
            q_values.append(value) 
            print(key,'candidates # :',len(value))
        
        self.candidate = dict(zip(q_keys, q_values))
        print('correlation_metabolism_path_mapped.')
    
    

    
   
    
    
    def violin(self, gene, save=True ):
        from qunor.gseatool import Violin_data
        from qunor.check_parser import gene_name_checker
        
        gene = gene_name_checker(gene)
        file_name = f'{gene}_violin.txt'
        data = Violin_data(self.raw,gene)

        if save==True:
            data.to_csv(file_name)
        else:
            ...
        return data 
  

    def pivot(self,gene):
        metabolism_path = self.metabolism_path.keys()
        data =[]
        for i in metabolism_path:
            temp_v = self.log_tumor.loc[self.metabolism_path[i]+[gene]].T.corr()
            data.append(temp_v.loc[gene])
        return data

    @staticmethod
    def get_gmt(gmt_path = '/home/zhuj/data_dicting/'):
    #    path = path

        gmts = [i for i in os.listdir(gmt_path) if i.endswith('gmt')]
        name = [os.path.basename(i).split('.')[0] for i in gmts]
        file_path = [os.path.join(gmt_path,i) for i in gmts]
        return dict(zip(name,file_path))
    def go_gsea(self, genes, gmt = '1', res ='res_dict',log = 'gsea_nohup',fdr=.25):
        df = self.__normalize_path
        assert os.path.exists(df),'1'
        from qunor import tool4
        file_name = str(hash(time.ctime())).strip('-')

        
        task_q = file_name +'.task_q'
        tool4.iterable_to_txt( genes, file_name =task_q)

        if gmt =='all':
            gmt = "/home/zhuj/data_dicting/c2.all.v7.5.1.symbols.gmt"

        log = log + '_' +file_name +'.log'
        from qunor import tool4
        
        python = tool4.decide_which_python()
        batch_gsea_file = tool4.decide_which_batch_gsea()

        cmd =f'nohup {python} {batch_gsea_file} -n {df} -t {task_q} -r {res} -gmt {gmt} -fdrthresh {fdr} > {log} 2>&1 &'
        print(cmd)
        try:
            os.system(cmd)
            time.sleep(4)
            os.remove(task_q)
        except Exception as e:
            print(e)
    # === 将实例pickle化===
    def go_pickle(self, file_name = None):
        
        match file_name:
            case None:
                """no name provided, we need to see if it is a TCGA datasets
                没有提供存储文件名参数的时候，先判断是否是TCGA的数据"""
                try:
                    samples = self.raw.columns

                    match 'TCGA' in samples[0]:
                        case True:

                            # we record the cancer type(s) into our .pickle, 将癌种名称写入pickle
                            cancer_type = samples.map(samples_cancer_type).unique().dropna()
                            assert len(cancer_type) != 0 ,'after dropna, no cancer types found basing on dictionary, have to manually dessign it'
                            
                            cancer_name = '_'.join(cancer_type)

                            # print(cancer_type,cancer_name)

                            file_name = cancer_name +'_normalize.pickle'
                            self.__cancer_name = cancer_name
                        case False:
                            
                            file_name = str(hash(time.time())) + '_normalize.pickle'

                            
                    with open(file_name,'wb') as f:

                        pk.dump(self,f)
                        print(f'{file_name} saved!')
                            
                except Exception as e:
                    print(e,'normalize.py ')

            case str():

                    with open(str(file_name),'wb') as f:
                                pk.dump(self,f)
                    
        
            case _:
                print('unsupported format')

                
                

if __name__ == '__main__':
    # m=Normalize(normalize_path="/home/zhuj/tcga/coad/RNF183/normalizeExp.txt"
    # ,deg_path="/home/zhuj/tcga/coad/RNF183/edgerOut.xls"
    # # ,only_stemness=True),
    # ,corr=False)
    # m.stemness(save=False)
    # m.go_pickle(file_name=None)
    pass
