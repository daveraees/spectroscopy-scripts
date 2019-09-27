def select_files(filelist,common_prefix,common_suffix):
    files_to_load = []
    datalabels = [];
    for filename in filelist:
        datalabel = filename[len(common_prefix):]
        if filename[:len(common_prefix)] == common_prefix:
            if filename[-len(common_suffix):] == common_suffix:
                files_to_load.append(filename)
                datalabels.append(datalabel)
    return (files_to_load,datalabels)  

def multiload_OO(filenames,data_labels=None):
    SSA = pd.DataFrame()
    index = 0
    if data_labels == None:
        data_labels = filenames
    for filename in filenames:
        dataset = pd.read_table(filename,\
                                 skiprows=17,skipfooter=1,header=None,\
                                 names=['wl',data_labels[index]],index_col=0,\
                                 engine='python')
        index += 1
        if sum(SSA.shape) == 0:                
                SSA = dataset;
        else:
                SSA = pd.concat([SSA, dataset], axis=1);
    return SSA
