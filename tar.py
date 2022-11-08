#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 15:39:37 2022

@author: ApriPro
"""


import os,shutil,gzip,tarfile,time
from tqdm import tqdm as tq
import pandas as pd


def unzipp_tar(tar,): 
    
    target_directory = os.path.splitext(tar)[0]
    
    try:
        print(f'1.1. making {os.path.basename(target_directory)}....')
        os.mkdir(target_directory)
    except FileExistsError:
        print(f'1.1. {target_directory} File directory exsists...,Overwriting it')
        pass
        
    tar = tarfile.open(tar,'r')
    
    tar.extractall( target_directory)
    tar.close()
    
    file_number = len(os.listdir(target_directory))
    
    print('\n', "======", file_number, ' files in total ======')
    
    return target_directory

#找到gz文件
def find_postfix_files(target_directory,postfix='.counts.gz'):
    
    postfix_files=[] ##这是 .gz文件的地址集

    print(f'3. copying gz to{target_directory}...')
    
    for r,d,files in tq(list(os.walk(target_directory)),desc=f'查找*{postfix}文件'):
        for file in files:
            if file.endswith(postfix):
                postfix_files.append(os.path.join(r,file))
    
    if len(postfix_files) ==0:
        raise AttributeError(f'No gz found in {target_directory}')
                
    return postfix_files
                    
    
def gunzip_files(filename,replace=True, return_=False):
    
    new_file_name =os.path.basename(filename).replace('.gz','')
    
    dirname =os.path.dirname(filename)
    
    f = gzip.GzipFile(filename=filename)
    
    with open(os.path.join(dirname,new_file_name),'wb') as d:
        d.write(f.read())
        
    f.close()
   
    if replace:
        try:
            os.remove(filename)
        except Exception as e:
            print(e,)
    if return_:
        return new_file_name
    
def get_counts(tar_file,target='Counts_unzipped', **kwargus):
    target = target
    if 'on' in kwargus.keys():
        dir =kwargus['on']
        target = os.path.join(dir,target)
    
    if tar_file.endswith('.tar'):
    # 如果是tar压缩文件：
        
        try:
            os.mkdir(target)
            print('1. target directory build..')
        except FileExistsError:
            pass
        
        directory =unzipp_tar(tar_file)
        print('2. unzipping tar...')
        
        try:
            os.remove(directory)
        except PermissionError as e:
            print(e,'\n',
                  f'注意:貌似权限遇到了挑战，无法自动删除{directory}','\n')
            
        gz_files=find_postfix_files(directory,postfix='.counts.gz')
        
        assert len(gz_files)!=0 ,  f'there are no gz_files in {directory}'
        
        for i in gz_files:
            name_i = os.path.basename(i)
            shutil.copy(i,os.path.join(target,name_i))
        
        
        temp_v= [os.path.join(target, i) for i in os.listdir(target)]
        counts_file_number =len(temp_v)
        
        for i in tq(temp_v, desc='4. 解压counts.gz'):
            gunzip_files(i)
            
            
        print(f'5. {counts_file_number} counts files collected in dir={target}')
        
        return target

    if tar_file.endswith('.gz'):
        print('dealing with .gz, trying to call gunzip in system cmd...')
        # 如果是.gz 文件，调用系统环境的gunzip 解压成tar，
        tar_name = tar_file.replace('.gz','') #tar的名字
        cmd = f'gunzip {tar_file}' #系统命令 gunzip **.gz
        os.system(cmd)
        os.system('exit')
        ## 得到tar，再调用原函数
        return get_counts(tar_name)
        
def postfix_file_return(postfix  ,
                        directory , 
                        ab_path=True,
                       walk=0):
    
    assert isinstance(postfix,str) ,f'postfix shall be given as str type, got:{type(postfix)}'
    
    
    
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
                return(value) #列表，返回directory下特定后缀文件地址
            
    else:
        pass #未来完善更深的遍历 2022年2月28日 —zhuj
        

def counts_matrix_maker(directory:str,
                        save = True,
                        postfix = '.counts',mapping=True,
                        ):
    
#——————————————第一步，调用postfix_file_return，找到所有目录下的counts——————————————————————————————————
#——————————————该函数返回一个list——————————————————————————————————

    counts_files = postfix_file_return(ab_path=True,
                                       
                                      postfix=postfix,directory=directory
                                      
                                      )


    counts_file_number = len(counts_files)
    
    
    
    assert counts_file_number!=0, 'counts_files should not be empty'
    
    
    data = [] #矩阵的主体
    columns = [] # 样本ID
    
#——————————————第二步，遍历该list，合并成一个matrix——————————————————————————————————
    
        
    for file in tq(counts_files,desc=f'6. 合并共计{counts_file_number}个counts文件'): 
        # counts_files 里应该都是counts
        assert file.endswith(postfix), f'{os.path.basename(file)} not ends with {postfix}'
        
        
        file_name= (os.path.basename(file))
        assert file_name!= ''
        
        columns.append(file_name) 
        f= open(file).read()

        strings =[i.split('\t') for i in f.split('\n')]

        for n in enumerate(strings):
            if len(n[1])==1:
                strings.pop(n[0])

        data.append(dict(zip((i[0] for i in strings),(float(x[1]) for x in strings))))
        
    df= pd.DataFrame(data).T
    df.columns=columns
# =============================================================================
    return df #"""numpy.float64 matrix"""
