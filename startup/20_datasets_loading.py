# file: 20_datasets_loading.py
# written by David Rais, 2016, david.rais@atlas.cz
#

def load_TA (filename,nrows=513):
    '''This is script for loading the csv datafile output from HELIOS transient absorption spectrometer
    returns pandas DataFrame with delay times in column labels and wavelengths in row labels'''
    dataset = pd.read_table(filename, index_col=0, sep=',', decimal='.', header=None, nrows=nrows)                        
    dataset = dataset.T.set_index(0).T # transform the column labels into float64 type
    return dataset

def multiload_TA (filenames):
    """Loads multiple transient absorption dataset into a pandas.Panel
    Alligns the timestamps of each dataset along the first in the list"""
    # load datasets for perpendicular polarization in dict
    all_scans = {n: scan for n, scan in enumerate([load_TA(x) for x in filenames])}
    # align the timestamps of each spectra with first in the list (they differ one or two notches)
    for i in range(len(all_scans)):
        if i > 0:
            all_scans[i].columns = all_scans[0].columns
    # convert the dict dataset into a pd.Panel data type:
    return pd.Panel(all_scans)

def load_TIMP (filename,file_format='Time explicit',mangle_timestamps_add=1e-6):
    '''This is script for loading the transient absorption datafile output from GLOTARAN and TIMP software
    it returns pandas DataFrame with delay times in column labels and wavelengths in row labels
    It only recognizes the time-gated spectra in 'Time explicit' format'''
    # read file header:
    datafile = open(filename,mode='rt')
    headerline = datafile.readline().strip('\n')
    descriptor = datafile.readline().strip('\n')
    #print (headerline.strip('\n'))
    #print (descriptor.strip('\n'))
    formatTIMP = datafile.readline()
    if file_format == 'Time explicit':
        if formatTIMP.strip('\n') == file_format:
            Intervalnr = int(datafile.readline().split()[1])
            #print (Intervalnr)
            timestamps = [np.float(x) for x in datafile.readline().split()]
            timestamps.insert(0,9e99)
            # find and mangle duplicate timestamps by adding some negligible time value:
            for idx in range(len(timestamps)):
                if idx == 0:
                    pass
                else:
                    if timestamps[idx] in timestamps[:idx]:
                        print ('duplicate found:', idx, timestamps[idx])
                        timestamps[idx] = max(timestamps[1:idx]) + mangle_timestamps_add
                        print ('changed to:', idx, timestamps[idx])
            dataset = pd.read_table(datafile, index_col=0, sep='\t', decimal='.', \
                                    header=None, names=timestamps, usecols=range(Intervalnr+1))
        else: 
            dataset = None
            print ('The file seems to display wrong datafile format.')
    else:
        dataset = None
        print ('Unknown file format:', file_format)
    datafile.close()
    return dataset

def multiload_TIMP (filenames):
    """Loads multiple transient absorption dataset into a pandas.Panel
    Alligns the timestamps of each dataset along the first in the list"""
    # load datasets for perpendicular polarization in dict
    all_scans = {n: scan for n, scan in enumerate([load_TIMP(x) for x in filenames])}
    # align the timestamps of each spectra with first in the list (they differ one or two notches)
    for i in range(len(all_scans)):
        if i > 0:
            all_scans[i].columns = all_scans[0].columns
    # convert the dict dataset into a pd.Panel data type:
    return pd.Panel(all_scans)
def load_OO (filename, datalabel=None):
    '''
    This script loads the dataset for steady-state fibre-optic spectrometer from Ocean Optics
    returns DataFrame with row index with wavelengths and absorbances (general spectra) in the column
    '''
    if datalabel == None:
        datalabel = filename;
    dataset = pd.read_table(filename,skiprows=17,skipfooter=1,\
                            header=None,names=['wl',datalabel],index_col=0,engine='python');
    return dataset
    
def multiload_2coll_spec (filenames, datalabels=[], **kwargs):
    '''This script loads the dataset for steady-state fibre-optic spectrometer from Ocean Optics
    returns DataFrame with row index with wavelengths and absorbances (general spectra) in the column
    each column will contain data from single file
    datalabels.
    I accepts list of files to load'''
    SSA = pd.DataFrame()
    count = 0
    if not 'engine' in kwargs.keys():
        kwargs['engine']='python'
    if datalabels == []:
        datalabels = filenames
    for filename in filenames:
        dataset = pd.read_table(filename,skiprows=1,skipfooter=0,usecols=[0,1],\
                              header=None,names=['wl',datalabels[count]],index_col=0,**kwargs);
        count += 1
        if sum(SSA.shape) == 0:
            SSA = dataset
        else:
            SSA = pd.concat([SSA,dataset], axis=1)
    return SSA

def multiload_2coll_spec (filenames, datalabels=[]):
    '''This script loads the dataset for steady-state fibre-optic spectrometer from Ocean Optics
    returns DataFrame with row index with wavelengths and absorbances (general spectra) in the column
    each column will contain data from single file
    datalabels.
    I accepts list of files to load'''
    SSA = pd.DataFrame()
    count = 0
    if datalabels == []:
        datalabels = filenames
    for filename in filenames:
        dataset = pd.read_table(filename,skiprows=1,skipfooter=0,usecols=[0,1],\
                              header=None,names=['wl',datalabels[count]],index_col=0,engine='python');
        count += 1
        if sum(SSA.shape) == 0:
            SSA = dataset
        else:
            SSA = pd.concat([SSA,dataset], axis=1)
    return SSA

def files_select (files, common_prefix,common_suffix,ident_str):
    datalabels = []
    filenames = []
    for filename in files:
        if filename[:len(common_prefix)] == common_prefix:
            if ident_str in filename[len(common_prefix):]:
                filenames.append(filename)
                datalabels.append(filename[len(common_prefix):-len(common_suffix)])
    return (filenames, datalabels)

# calculate rates and uncertainities
def lftm_from_rate_w_err (kVect,kVect_stderr=[]):
    """
    calculate lifetimes from rates with standard errors
    """
    lifetimes=[1/ki for ki in kVect ]
    if len(kVect)==len(kVect_stderr):
        lftmErrs = [ki[1]/(ki[0]**2) for ki in zip(kVect,kVect_stderr) ]
    else:
        lftmErrs = []
    return (lifetimes, lftmErrs)
    
def load_TA_result (filename, lifetimes=[], rates=[], rates_stderrs=[]):
    '''
    Fuction to load single dataset of the results of Target Analysis generated as 
    an ASCII export output of Glotaran software.
    It will throw away multiple X-values collumns, and keep only the first.
    '''
    dataset = pd.read_table(filename, index_col=0)
    dataset = dataset.iloc[:,0::2]
    if len(lifetimes)==0:
        if len(rates)==len(dataset.columns):
            if len(rates_stderrs)==0:
                lifetimes=[1/ki for ki in rates ]
            else:
                lftms,lftmErrs = lftm_from_rate_w_err (rates,rates_stderrs)
                lifetimes=[]
                for lftm in zip(lftms,lftmErrs):
                    orders_lftm = np.log10(lftm)
                    min_order = (min(orders_lftm))
                    #print(min_order)
                    if min_order < -2:        
                        lifetimes.append('%.3f ± %.3f' % lftm)
                    else:
                        if min_order >= -2 and min_order < -1:        
                            lifetimes.append('%.2f ± %.2f' % lftm)
                        else:
                            if min_order >= -1 and min_order < 0:
                                lifetimes.append('%.1f ± %.1f' % lftm)
                            else:
                                lifetimes.append('%d ± %d' % lftm)
        else:
            lifetimes=[('#%d' % l) for l in range(len(dataset.columns)) ]
    dataset.columns=lifetimes
    return dataset

def load_SADS (GTA_fit_params, SADS_filename):
    '''
    Fuction to load single dataset of the results of Target Analysis generated as 
    an ASCII export output of Glotaran software.
    It will throw away multiple X-values collumns, and keep only the first.
    '''
    dataset = pd.read_table(SADS_filename, index_col=0)
    dataset = dataset.iloc[:,0::2]
    rates = GTA_fit_params['kVect']
    rates_stderrs = GTA_fit_params['kVect_stderr']
    SADS = load_TA_result(SADS_filename,rates=rates, rates_stderrs=rates_stderrs)
    dataset.columns=lifetimes
    return dataset
    
def load_TA_result_fits (filename):
    '''
    Fuction to load exported single-trace results of kinetic analysis generated as 
    an ASCII export output of Glotaran software.
    It will separate the multiple curves in 4-column segments "x,data,fit,residuals"
    x could be either time, or wavelength
    '''
    dataset = pd.read_csv(filename, index_col=None) #read the entire table
    subsLength = 4
    TA_result_fits = {}
    for subsetNr in range(dataset.shape[1] // subsLength):
        subset = dataset.iloc[:,subsetNr*subsLength:subsetNr*subsLength+(subsLength)]
        parameter = subset.iloc[:,0].name
        TA_result_fits[float(parameter)] = subset.set_index(parameter) #use the first col and row as indices
    return TA_result_fits

def multiLoad_TA_results (filenames={}, lifetimes=[], rates=[]):
    '''
    Function to create dictionary of the results for single Target analysis.
    '''
    dataset = {}
    if len(lifetimes)==0:
        if len(rates)>0:
            lifetimes=[1/ki for ki in rates ]
    for key in filenames.keys():
        if key == 'SADS':
            SADS = load_TA_result(filenames[key] , lifetimes=lifetimes)
            dataset[key]=SADS
        if key == 'KINS':
            KINS = load_TA_result(filenames[key] , lifetimes=lifetimes)
            dataset[key]=KINS
        if key == 'EADS':
            EADS = load_TA_result(filenames[key] , lifetimes=lifetimes)
            dataset[key]=EADS
        if key == 'DADS':
            DADS = load_TA_result(filenames[key] , lifetimes=lifetimes)
            dataset[key]=DADS
        if key == 'RKINS':
            RKINS = load_TA_result_fits(filenames[key])
            dataset[key]=RKINS
        if key == 'RSPEC':
            RSPEC = load_TA_result_fits(filenames[key])
            dataset[key]=RSPEC
    return dataset
    
def getFitResults (timpresFname):
    """
    this function reads the timpres file produced by Glotaran software, 
    and re-constructs the original dataset, the fitted dataset, SADS, DADS, 
    fraction temporal profiles, chirp parameters, IRF parameters
    """
    def correct_t0_frac(frac,timeLim=2.0,inplace=True):
        """Workaround for correction of time-zero for the temporal profiles
        of concentrations (fractions)
        the funciton finds an inflex point in the early delay time, where the
        time steps are supposed to be in constant linear intervals.
        frac:
            dataset with the index of times and fraction profile in each collumn
        inplace:
            True: makes the correction to the dataFrame in argument
            False: returns the calculated time of the inflex point:
        """
        earlyTotal = gtaFitResult['frac'].loc[:timeLim].sum(axis=1)
        max = earlyTotal.max()
        diff=np.diff(earlyTotal.values,n=2)
        diffTime = earlyTotal.index[1:-1] # without the end time-points
        idxMax = np.argmax(diff)
        tMax = diffTime[idxMax]
        if inplace:
            rVal = None
            frac.index -= tMax
        else:
            rVal = tMax
        return rVal #def correct_t0_frac
    transcript_vectors =\
        {
        'kineticParameters': ('kVect', 'kVect_stderr') ,
        'irfpar': ('IRF' , 'IRF_stderr'),
        'parmu' : ('Parmu','Parmu_stderr'),
        'lamdac': 'centralWL',
        'rms':'residual_stderr',
        'x'  : '_delays',
        'x2' : '_probeWLs'
        } 
    transcript_matrices =\
        {
        'spectra' : ('SADS', 'DADS'),
        'spectraErr' : 'SADS_stderr',
        'concentrations': 'frac',

        }
    transcript_surfaces = \
        {
        'fittedTraces': 'fittedTraces',
        'traces': 'traces',
        }
    transcript_strings = {'datasetName':'fitName'}
    GTA_dict_keys = ['annotations',
         'classdesc',
         'clpequ',
         'coh',
         'concentrations',
         'datasetName',
         'eigenvaluesK',
         'fittedTraces',
         'get_class',
         'intenceIm',
         'irfpar',
         'jvec',
         'kineticParameters',
         'kinscal',
         'lamdac',
         'maxInt',
         'minInt',
         'orheigh',
         'orwidth',
         'oscpar',
         'parmu',
         'partau',
         'prel',
         'residuals',
         'rms',
         'specdisppar',
         'spectra',
         'spectraErr',
         'spectralParameters',
         'traces',
         'type',
         'x',
         'x2']

    with open (timpresFname, 'rb') as resultsFile:
        filestring=resultsFile.read()
    pobj = javaobj.loads(filestring)
    
    gtaFitResult = dict() # objet to store the exracted fit results
    ### STRINGS ###
    for GTAvec in transcript_strings.items():
        vec = getattr(pobj,GTAvec[0])
        if vec != None:
            gtaFitResult[GTAvec[1]] = vec 
    ### VECTORS ###
    for GTAvec in transcript_vectors.items():
        value = getattr(pobj,GTAvec[0])
        if value != None:
            vec = np.array(getattr(pobj,GTAvec[0]))    
            if type(GTAvec[1]) == tuple:
                how_many = len(GTAvec[1])  # in how many vectros the timpres filed should be split
                vec = np.reshape(vec,(-1,how_many))
                for i in range(how_many):
                    gtaFitResult[GTAvec[1][i]] = vec[:,i]
            else:
                gtaFitResult[GTAvec[1]] = vec
        else:
            pass
    Tindex = gtaFitResult['_delays']-gtaFitResult['IRF'][0]
    WLindex = gtaFitResult['_probeWLs']
    ### MARTICES ###
    for GTAvec in transcript_matrices.items():
        javaMatrix = getattr(pobj,GTAvec[0])
        if javaMatrix != None:
            values = np.array(javaMatrix.A)
            shape = np.array(values.shape)
            min_shape = np.min(shape)
            min_idx = np.argmin(shape)
            if min_idx == 0:
                index = WLindex
                values = values.T
            if min_idx == 1:
                index = Tindex
            if type(GTAvec[1]) == tuple:
                how_many = len(GTAvec[1])  # in how many vectros the timpres filed should be split
                grouping = int(min_shape/how_many)
                for i in range(how_many):
                    i_end = i+ grouping
                    colNames = [('%s#%d' % (GTAvec[1][i], n)) for n in range(1,grouping+1)]
                    gtaFitResult[GTAvec[1][i]] = pd.DataFrame(data=values[:,i:i_end],
                                                              index=index, 
                                                              columns=colNames,copy=True)
            else:
                colNames = [('%s#%d' % (GTAvec[1], n)) for n in range(1,min_shape+1)]
                gtaFitResult[GTAvec[1]] = pd.DataFrame(data=values,
                                                       index=index, 
                                                       columns=colNames,
                                                       copy=True)
            if GTAvec[1] == 'frac': # DIRTY workaround corection. 
                                    #Should calculated the fractions from the kinetic constatnt instead
                correct_t0_frac(gtaFitResult[GTAvec[1]])
        else:
            pass
    ### SURFACES ###
    for GTAvec in transcript_surfaces.items():
        javaMatrix = getattr(pobj,GTAvec[0])
        if javaMatrix != None:
            values = np.array(javaMatrix.A)
            shape = np.array(values.shape)
            min_shape = np.min(shape)
            min_idx = np.argmin(shape)
            colNames = [('%s#%d' % (GTAvec[1], n)) for n in range(1,min_shape+1)]
            surface = pd.DataFrame(data=values.T,
                                   columns=gtaFitResult['_delays'], 
                                   index=gtaFitResult['_probeWLs'],
                                   copy=True)
            if 'Parmu' in gtaFitResult.keys():
                gtaFitResult[GTAvec[1]] =   correct_chirp(surface,
                                                          IRF=gtaFitResult['IRF'][1],
                                                          t0=gtaFitResult['IRF'][0],
                                                          dispersion=gtaFitResult['Parmu'],
                                                          centralWL=gtaFitResult['centralWL'])
            else:
                gtaFitResult[GTAvec[1]] =   correct_chirp(surface,
                                                          IRF=gtaFitResult['IRF'][1],
                                                          t0=gtaFitResult['IRF'][0],
                                                          )
        else:
            pass
    return gtaFitResult

def getResultsFromSummary (resultsDirName):
    """
    extracts the information from the glotaran results SUMMARY file in given directory
    returns dictionary with the following fields:
        'iterNum', (number of performed iterations)
        'residual_stderr', ...
        'kVect','kVect_stderr' (vectors of the best fit parameters and estimated standard errors...)
        'IRF', 'IRF_stderr' -ditto-
        'Parmu', 'Parmu_stderr' -ditto-
    """
    for fname in os.listdir(resultsDirName):
        if 'summary' in fname:
            summaryFname = resultsDirName+'/'+fname
    with open(summaryFname, 'r') as sf:
        iters = 0
        rms = 0
        summary = {}
        while True:
            line = sf.readline()
            if line == '': # EOF was hit
                break
            if 'Number of iterations:' in line:
                iters=getScanNum(line)
                summary.update({'iterNum': iters })
            if 'Final residual standard error:' in line:
                rms = get_float(line)
                summary.update({"residual_stderr" :rms[0]})
            if 'Estimated Kinetic parameters: Dataset1:' in line:
                kVect = get_float(line)[1:]
                kVect_stderr = get_float(sf.readline())
                summary.update({'kVect':kVect, 
                                'kVect_stderr':kVect_stderr,})
            if 'Estimated Irf parameters: Dataset1:' in line:
                IRF = get_float(line)[1:]
                IRF_stderr = get_float(sf.readline())
                summary.update({'IRF': IRF,
                                'IRF_stderr':IRF_stderr,})
            if 'Estimated Parmu: Dataset1:' in line:
                Parmu = get_float(line)[1:]
                Parmu_stderr = get_float(sf.readline())
                summary.update({'Parmu':Parmu,
                                'Parmu_stderr':Parmu_stderr})
    return summary
    
def getDatasetResultFnamePart (resultXMLlogFname):
    """
    extracts the information about the filename of the cashed datased 
    from the glotaran results file .xml_overview.xml
    and the corresponding results filename
    """
    with open(resultXMLlogFname, 'r') as sourceFile:
        doc = xmltodict.parse(sourceFile.read())
    #print(doc['GtaResult']['datasets'].keys())
    datasetFnamePart = ''
    resultFname = ''
    try:
        datasetFilePath=doc['GtaResult']['datasets']['datasetFile']['path']
        datasetFnamePart=datasetFilePath.split('/')[1][:-14]
        resultFname = doc['GtaResult']['datasets']['resultFile']['filename']
    except:
        pass
    #print (datasetFilePath)
    return (datasetFnamePart,resultFname)

def average_scans(scans_dict):
    """
    Example data selection code:
    scans_dict = {  key: TA_scans[key] for key in ave_scans }
    """
    mean_TA = None
    count = 0
    for scan in scans_dict.items():
        if type(mean_TA) == type(None):
            mean_TA = scan[1].values
        else:
            mean_TA += scan[1].values
        count += 1
    mean_TA /= count
    mean_TA = pd.DataFrame(index=scan[1].index,
                                    columns=scan[1].columns,
                                    data=mean_TA)
    return mean_TA

def find_summary_load_detail (DatasetName, whichFitName=None, GTArootDir=None):
    # Glotaran analyses folder path
    if GTArootDir == None:
        userprofile = os.getenv('USERPROFILE')
        GTArootDir = userprofile+'/Documents/FS-analyses/DPP_thin_films/'
        
    resultsDirs_spec = find_analyses (DatasetName, GTArootDir)
    if whichFitName == None:
        summary = None
    for summary in resultsDirs_spec:
        datasetFname = summary['datasetFname']
        fitName = summary['fitName']
        timpresFname = GTArootDir + summary['dirName'] + '/' + summary['fitName'] + '.timpres'
        #repreWLS = get_float(repreWLSlist.\
        #                 where(repreWLSlist["datasetFname"]==datasetFname).\
        #                 dropna()['repreWLS'].iloc[0])    
        #if True:
        if fitName==whichFitName:            
            detail = getFitResults(timpresFname)
            summary['detail'] = detail
            break
    return summary

