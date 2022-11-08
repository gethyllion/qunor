# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 09:49:12 2022

@author: Zhu jun ex Hangzhou
"""
import os,time
import pandas as pd
from qunor import tar, gseatool, check_parser
from qunor.dict_data import gene_biotype_dict as rna_dict


def timmer(f):
    def wrapper(*argus, **kwargus):
        start = time.time()
        res = f(*argus, **kwargus)
        name  = f.__name__
        print(round((time.time()-start),3),f'seconds for process.{name} to be finished')
        return res

    return wrapper

# =============================================================================
# 1. DataBuilder.Matrix类

# =============================================================================
class Matrix:
    """

    DataBuilder.Matrix
            对象的属性
        self.matrix = None  # 矩阵数据的主体，pd.DataFrame
        self.stat = None # 正常和非正常样本的统计，默认False是正常样本， collections.Counter()对象
        self.sample_dict = None  # 字典，记录了 {文件名：样品ID}

        self.__json = None # 是样本数据， pd.DataFrame
        self.__json_path = None # 文件地址
        self.__counts_dir = None  # 解压tar或者总的gz文件后，存放counts压缩文件的目录名

        self.__sample_q = None  # 样品ID
        self.__file_name_q = None  # counts文件的名字，用于字典对照

        self.__json_path = None  # json文件路径
        self.__tar_path = None  # 压缩文件的路径
"""
    __annot = None

    def __init__(self,

                 sample_dict=None,
                 stat=None):
        self.__total_ens_number = None
        self.stat = None
        self.matrix = None  # 矩阵数据的主体，pd.DataFrame
        # =============================================================================
        self.__json = None
        self.__json_path = None
        self.__counts_dir = None  # 解压tar或者总的gz文件后，存放counts压缩文件的目录名
        self.__sample_dict = None  # 字典，记录了 {文件名：样品ID}
        self.__sample_q = None  # 样品ID
        self.__file_name_q = None  # counts文件的名字，用于字典对照
        self.__json_path = None  # json文件路径
        self.__tar_path = None  # 压缩文件的路径
        self.__unmatched_gene_number = None

    # =============================================================================
    # 2.1. 实例方法——加载数据
    #         调用该方法需要传入 2个必需参数，json文件路径 和 压缩文件的路径
    # =============================================================================
    @staticmethod
    def __parse_json(file_path):
        with open(file_path) as f:
            strings = f.read()

        return eval(strings.replace('null', 'None'))

    @staticmethod
    def __return_dict_from_txt(file_path):
        with open(file_path, 'r') as f:
            data = f.read()
        return eval(data)

    @staticmethod
    def __get_ab_path(file_path):
        return os.path.abspath(file_path)


    def __load_data(self, json_path=None, tar_path=None, annot_file=None, **kwargus):

        Matrix.__check_availability_of_file(json_path)

        self.__tar_path = Matrix.__get_ab_path(tar_path)
        self.__json_path = Matrix.__get_ab_path(json_path)

        self.__json = pd.json_normalize(data=Matrix.__parse_json(json_path))

        self.__sample_q = (i for i in pd.json_normalize(data=Matrix.__parse_json(json_path),
                                                        record_path=['associated_entities', ])['entity_submitter_id'])

        self.__file_name_q = (i.replace('.gz', '') for i in pd.read_json(json_path)['file_name'])

        self.__counts_dir = tar.get_counts(tar_path,**kwargus)

        self.__sample_dict = dict(zip(self.__file_name_q, self.__sample_q))

    # =============================================================================
    # 2.2 实例方法——制作矩阵
    # =============================================================================

    def __make_and_map_matrix(self, **kwargs):

        if self.__counts_dir is None or self.__json is None:
            raise ValueError(
                ' malapropos data for constructing the matrix,the instance.json None or instance.__counts_dir is None ')

        else:
            matrix = tar.counts_matrix_maker(directory=self.__counts_dir)

            matrix.columns = matrix.columns.map(self.__sample_dict)

            self.matrix = matrix

            self.stat = gseatool.normal_sample_detector(columns=list(self.__sample_dict.values()), stat=True)

    @classmethod
    def __annotate(cls, annot_file=None, online=False):

        if isinstance(annot_file, str) and annot_file.endswith('.gtf'):
            entrez_dict_class = check_parser.Entrez_dict()

            entrez_dict_class.load_gtf_data(gtf_file=annot_file)

            cls.__annot = entrez_dict_class.GTF_symbol_dict

        if isinstance(annot_file, str) and annot_file.endswith('.txt'):
            D = Matrix.__return_dict_from_txt(annot_file)

            cls.__annot = D

        if online is True:
            print('online querrying would be implemented in near future :)')
            ...
    def rearrange_columns_of_matrix(self):
        temp_col = gseatool.normal_sample_detector(df=self.matrix)
        self.matrix.loc['temp_v'] = temp_col
        new_col_order = self.matrix.T.sort_values(by='temp_v').index
        self.matrix = self.matrix.loc[:, new_col_order].drop('temp_v')
        print('matrix columns rearranged, tumor samples follows the normal ones.')

    def __len__(self):
        return tuple(len(self.matrix.index), len(self.matrix.columns))
    
    
    @timmer
    def self_establish(self, json_path=None, tar_path=None, annot_file=None, **kwargus):
        from collections import Counter
        print(kwargus)
        if annot_file is None:
            try:
                from qunor import check_parser
                annot_file = check_parser.find_dict()
            except Exception as e:
                print(e)
                pass
        Matrix.__check_availability_of_file(annot_file)

        self.__load_data(json_path=json_path, tar_path=tar_path,**kwargus)
        self.__make_and_map_matrix()

        self.matrix.index = [i.split('.')[0] for i in self.matrix.index]
        self.__annotate(annot_file=annot_file)

        self.matrix.index = self.matrix.index.map(self.__annot)
        self.matrix.index.name = 'id'
        self.__unmatched_gene_number = Counter(self.matrix.index.isnull())[True] - 5
        self.__total_ens_number = len(self.matrix.index) - 5  # 这里减去5是因为 有5个报告测序质量的行

        self.matrix = self.matrix.loc[~self.matrix.index.isnull(), :]
        self.__rearrange_columns_of_matrix()

        dir_name = os.getcwd()
        if kwargus['on']:
            dir_name = kwargus['on']

            assert os.path.exists(dir_name)

            
            try:
                os.mkdir(os.path.join(dir_name, 'RNA_data'))
            except FileExistsError as e:
                print(e,r'已经存在了一个文件,尝试覆盖其中的total_RNAmatrix.txt')

        self.matrix.to_csv(os.path.join(dir_name, 'RNA_data/total_RNAmatrix.txt'))
        
        mRNAmatrix = self.matrix
        mRNAmatrix = mRNAmatrix.loc[[i for i in mRNAmatrix.index if rna_dict[i]=='protein_coding'],:]
        self.mRNAmatrix = mRNAmatrix

        mRNAmatrix.to_csv(os.path.join(dir_name, 'mRNAmatrix.txt'))
        
        self.date = time.ctime(time.time())

        self.reporta = pd.DataFrame([{'data':self.date,
                                    'num_tumor_case': (self.stat[True]),
                                     'num_normal_case': (self.stat[False]),
                                     'total_ens_number': self.__total_ens_number,
                                     'unmatched_gene_number': self.__unmatched_gene_number,
                                     'u/t': (self.__unmatched_gene_number / self.__total_ens_number),
                                     'count_dir': self.__counts_dir,
                                     'json_path':self.__json_path,
                                     'tar_path':self.__tar_path,
                                     'mRNA_num':len(self.mRNAmatrix.index)
                                     }]).T

        self.reporta.to_csv(os.path.join(dir_name,'reporta.txt'))


if __name__ =='__main__':
    pass