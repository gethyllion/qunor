from scipy import stats
from seaborn import heatmap,jointplot
import seaborn as sns,pickle
import matplotlib.pyplot as plt
import pandas as pd
from qunor import gseatool as gst
import numpy as np
import os


def make_four_directories():
    directories= ['附件1：数据分析','附件2：原图片','附件4：参考文献','拟题报告','附件3：图片排版']
    for directory in directories:

        try:
            os.mkdir(directory)
    
            print(f'{directory} made!')
        except Exception as e:
            print(f'cannot make {directory}', e, )
            pass

def pearson(x:pd.array or pd.Series or list or tuple
            ,y:pd.array or pd.Series or list or tuple):
    v = stats.pearsonr(x,y)
    return v

def correlation(genes, df):
    if len(genes)!=2:
        raise AttributeError('it should be only 2 two genes in such')
    fig = jointplot(df.loc[genes[0]])

def make_df(dir_name='res_dict'):
    dir_name = dir_name
    files = os.listdir(dir_name)
    files = [eval(open(os.path.join(dir_name,i)).read()) for i in files]
    df = pd.concat(pd.DataFrame(i).T for i in files)
    return df

def heat(matrix):
    fig = heatmap(matrix, annot=True)
    fig.set_xlabel('Gene')
    fig.set_ylabel('Gene')
    plt.savefig(f'{matrix.index.to_list()}_heatmap.pdf', bbox_inches='tight')
    plt.show()

def log_tumor(df):
    tem = gst.remove_normal_samples(df).T
    log = np.log2(tem+0.01)
    return log

def iterable_to_txt(iterable_object,txt_name,sep='\n'):
    with open(txt_name,'w') as f:

        for i in iterable_object:
            f.write(i+sep)

def reg(x, y, file_name = 'corr'):
    
    if hasattr(x,'name') and hasattr(y,'name'):
        
        assert isinstance(x.name, str), 'x.name is not str'
        assert isinstance(y.name, str), 'y.name is not str'
        file_name = x.name + '_' +  y.name +'_'+file_name
    
    corr , p_value = pearson(x,y)

    if p_value > .05:
            print('p>.05!')
    elif p_value < .05:
        
        annot = f'{x.name} corr={round(corr,2)},p<0.05'

        temp = x.copy()
        temp.name = annot
        assert hasattr(x,'name')
        

        f = sns.jointplot(x=temp,y=y,kind='reg')
        
        plt.savefig(f'{file_name}.pdf')
        plt.savefig(f'{file_name}.tiff')
        plt.show()
        return f    
        
def parse_txt(file,sep = '\n'):
    """return a list from sep-txt file"""
    with open(file) as f:
        f = f.read().split(sep)
    
    """try to remove blanks in the list """
    while '' in f:
        f.remove('')
    return f
    
    
def parser_pickle(file,):
    with open(file,'rb') as f:
        return pickle.load(f)

if __name__ =='__main__':
    import pickle
    data =pickle.load(open('x_y.pickle','rb'))
    reg(*data)
    # make_four_directories()
    print('done')
