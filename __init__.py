
if __name__ =='__main__':


    from qunor import gseatool as gst
    
    import pandas as pd
    import os,time
    
    print(time.localtime(time.time()))
    print('current_wd:',os.getcwd())
    print('gseatool:',gst.__file__,)
    

    task = eval(open('task_q.txt').read())
    
   
    
    print('task_type:',type(task))
    normalize = pd.read_table('normalizeExp.txt',index_col=0,sep='\t')
    gene_sets ="/home/zhuj/data_dicting/c2.cp.kegg.v7.5.1.symbols.gmt"
    try:
        dir_name = 'res_dict'
        os.mkdir(dir_name)
    except Exception as e:
        print(e)
        
    res = gst.gsea_batch(task_q = task, data = normalize,dir_name = dir_name,permutation_num = 1600,seed = 8, permutation_type='phenotype',method='signal_to_noise', processes = 5, no_plot = True, outdir = None, gene_sets = gene_sets,)
    res = repr(res)
  
  
    with open('results_gsea.txt','w') as f:
        
        f.write(res)
        
    print('saved —— results_gsea.txt ',)
    print(time.localtime(time.time()))
     
