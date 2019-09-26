# this file 22_dataset_saving.py
# contains functions for saving of transient absorption data
# written by David Rais, 2016, david.rais@atlas.cz

def save_TIMP (dataset, filename, descriptor=''):
    """
    save the dataset in Glotaran format
    """
    datafile = open(filename,mode='xt')
    datafile.writelines([os.path.abspath(os.curdir),'\\', datafile.name,'\n',
                         descriptor,'\n',
                         'Time explicit','\n',
                         'Intervalnr   '+str(len(dataset.columns)),'\n',
                         ])
    dataset.to_csv(datafile,header=True,sep='\t',index_label=False)
    datafile.close()
    return

def sendTAKins2OriginWKS(dataset,wks_title=None,wks_Name=None):
    """
    Sends the data to the opened Origin instance...
    """
    wks = origin.FindWorksheet(wks_title)
    if wks == None: # if the worksheet does not exist, create new one
        pageName = origin.CreatePage(2, wks_title, "origin")
        wks = origin.FindWorksheet(pageName) 
    colNum = 0;
    while True:
        if (wks.Columns.Count) > colNum:
            if len(wks.Columns[colNum].GetData(0,0,-1,-1))>0:
                colNum +=1;
            else:
                break
        else:
            break
    for dataPoint in dataset.index.values:
        origin.PutWorksheet(wks_title,dataPoint,-1,colNum)
        #wks.SetData(dataPoint,-1,colNum)
    wks.columns(colNum).LongName = 'Delay time'
    wks.columns(colNum).Units = 'ps'
    colNum +=1;
    for column in dataset.columns.values:
        for dataPoint in dataset[column]:
            origin.PutWorksheet(wks_title,dataPoint,-1,colNum)
        wks.columns(colNum).LongName = 'delta A'
        wks.columns(colNum).Units = 'OD'
        wks.columns(colNum).Comments = (str(column) + ' nm')
        colNum += 1;
    #wks.LongName = wks_LongName
    if wks.Name != None:
        wks.Name = wks_Name # modify the worksheet label
    
def sendTASpec2OriginWKS(dataset,wks_title=None,wks_Name=None):    
    wks = origin.FindWorksheet(wks_title)
    if wks == None: # if the worksheet does not exist, create new one
        pageName = origin.CreatePage(2, wks_title, "origin")
        wks = origin.FindWorksheet(pageName) 
    colNum = 0;
    while True:
        if (wks.Columns.Count) > colNum:
            if len(wks.Columns[colNum].GetData(0,0,-1,-1))>0:
                colNum +=1;
            else:
                break
        else:
            break
    for dataPoint in dataset.index.values:
        origin.PutWorksheet(wks_title,dataPoint,-1,colNum)
        #wks.SetData(dataPoint,-1,colNum)
    wks.columns(colNum).LongName = 'Wavelength'
    wks.columns(colNum).Units = 'nm'
    colNum +=1;
    for column in dataset.columns.values:
        for dataPoint in dataset[column]:
            origin.PutWorksheet(wks_title,dataPoint,-1,colNum)
        wks.columns(colNum).LongName = 'delta A'
        wks.columns(colNum).Units = 'OD'
        wks.columns(colNum).Comments = (str(column) + ' ps')
        colNum += 1;
    #wks.LongName = wks_LongName
    if wks.Name != None:
        wks.Name = wks_Name
    