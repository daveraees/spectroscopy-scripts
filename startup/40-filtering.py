# this file 40-filtering.py
# should contain functions for filtering and manipulation of transient absorption data
# written by David Rais, 2016, david.rais@atlas.cz

def subtract_TA_background (dataset,ave_idx):
    """
    Make the correction of baseline of the dataset by averaging the TA spectra given in the indexes and 
    subtracting the average from the each spectrum in the dataset

    dataset: the input transient absorption spectrotemporal data

    ave_idx: list of positional indices of the TA spectra, which should be used to calculate the baseline TA spectrum

    returns: subtracted dataset
    """
    background = dataset.iloc[:,ave_idx].mean(axis=1);
    TA_data = dataset.copy();
    for i in range (len(TA_data.columns.values)):
        TA_data.iloc[:,i] = TA_data.iloc[:,i] - background;
    return TA_data;

def calculate_anisotropy (parallel, perpendicular):
    """
    Calculate the anisotropy from tr5ansient absorption spectrotemporal data recorded with pump-probe polarization orientated in parallel and perpendicular

    parallel, perpendicular: the input transient absorption spectrotemporal datasets

    returns: dataset with calculated anisotropy

    Warning: The function expects the parallel and perpendicular datasets to have the same time gates, and uses time gates of the parallel also for the perpendicular and final anisotropy datasets
    """
    perpendicular.columns = parallel.columns
    ani = (parallel - perpendicular)/ (parallel + 2 * perpendicular)
    return ani

def calculate_arb(parallel,perpendicular,angle=54.7,align_timestamps=False):
    """
    Calculate the transient absorption signals for arbitrary angle between pump and probe
    polarizations. 
    
    parallel, perpendicular: the input transient absorption spectrotemporal datasets
    
    angle: defauls to 54.7 (magic angle) angle between pump and probe polarizations to be used in calculation, in degrees
    
    align_timestamps: (defaults to False) ?Use time gates of the parallel also for the perpendicular and final anisotropy datasets?

    returns: calculated spectrotemporal datasets
    """
    angle_rad = angle/ 180 * np.pi
    #align the timestamps of the two datasets, if explicitly requested:
    if align_timestamps:
        perpendicular.columns = parallel.columns
    arb = np.cos(angle_rad)**2 * parallel + np.sin(angle_rad)**2 * perpendicular
    return arb

def chirp_correction (WL=[],dispersion=[0],t0=0,centralWL=None):
    """
    This function calculates the pump-probe temporal overlap from the dispersion curve polynomial coefficients,
    available from the Glotaran fitting
    Input values:
    WL
    dispersion
    t0
    centralWL
    returns: value of time shift of the t0
    """
    if centralWL != None:
    #    if len(WL)==1:
    #        t0Shift = 0
    #        polyOrder=1
    #        for coef in dispersion:
    #            t0Shift += dispersion[polyOrder]*(WL-centralWL)^polyOrder
    #            polyOrder+=1
        if len(WL)>0:
            t0Shift=np.zeros(len(WL))
            for idxWL in range (len(WL)):
                polyOrder=1
                for coef in dispersion:
                    t0Shift[idxWL] += dispersion[polyOrder-1]*((WL[idxWL]-centralWL)/100)**polyOrder
                    polyOrder+=1
                t0Shift[idxWL] += t0
    return t0Shift
    
def calcChirpIndex (datasetWLindex,IRF=0,t0=0,dispersion=[],centralWL=0):
    """
    This function creates the chirp correction data. 
    It calculates the t0 position based on given polynomial coefficients
    Paraemters:
    datasetWLindex: array of probe wavelengths
    IRF: gaussian width parameter of teh temporal InstumentResponseFunction
    t0: base time zero correction (added to the chirp data)
    dispersion: polynomial coefficients of the dispersion curve (starting from order 1 upwards)
    centralWL: base wavelength for chirp correction calculation (as used in Glotaran/TIMP fitting software)
    It returns pd.DataFrame.Multiindex instance for embedding it into the original DataFrame dataset object
    """    
    #WLindex = l30ZnTF_PMMA['pump@310nm'].index.values
    t0index = chirp_correction(datasetWLindex, dispersion, t0, centralWL)
    wIRF = np.ones(len(datasetWLindex)) * IRF # the IRF is considered independent on the probe wavelength
    DataIndex = pd.DataFrame ( index=datasetWLindex, data=np.vstack([t0index,wIRF]).T, columns=['t0','wIRF'])
    return DataIndex
    
def correct_chirp(dataset,IRF=0,t0=0,dispersion=[],centralWL=0):
    """
    function to calculate the transformed dataset to compensate for instrument dispersion artifact (chirp corretion)
    using linear interpolation
    inputs:
    dataset: dataframe with TA data
    IRF: width of the instrument response funciton (not used at the moment)
    t0: shift of the origin of delay time at the centralWL
    dispersion: polynomial function describing the relative temporal offset with respect to the central WL
    centralWL: wavelength, at which the temporal overlap occurs at 't0'
    """
    # calculate the chirp correction data from the dispersion curve and basic timeshift (obtained e.g. by fitting)    
    chirpIndex = calcChirpIndex(dataset.index,
                                IRF,
                                t0,
                                dispersion,
                                centralWL)
    # create new dataset with chirp corrected data using simple linear interpolation method
    newCurves = pd.DataFrame(index=dataset.index,columns=(-t0+dataset.columns))
    #newCurves=deepcopy(dataset)
    for wl in dataset.index:
        #print(type(wl))
        times_offset = - chirpIndex.loc[wl,'t0'] + dataset.columns.values 
        ndTrace = dataset.loc[wl,:].astype(type(1.21654e-5), order='K', casting='unsafe')       
        newCurves.loc[wl,:] = np.interp(newCurves.columns.values,times_offset,ndTrace)
    return newCurves

def bg_chirp_correct (dataset,bg_rng=0,fitRes={}):
    """
    Performs background correction with the dataset and calculate chirp-corrected data, 
    using the fit results of chirp dispersion curve and IRF functions (using e.g. Glotaran)
    
    dataset: TA dataset in Pandas DataFrame format
    bg_rng : range of background indices, 
        if integer : number of spectra from the beginning; %%!
        if more int values: integer indices of spectra to be used for calculation
    fitRes : fit results of chirp dispersion curve and IRF functions (using e.g. Glotaran) in dictionary format
    
    Output: 
        bg- and chirp-corrected dataset
    """  
    
    if type(bg_rng) == int:
        dataset = subtract_TA_background(dataset,range(bg_rng))
    else:
            if len(bg_rng) > 1:
                dataset = subtract_TA_background(dataset,bg_rng)
    # calculate the chirp-corrected dataset:
    if len(fitRes) != 0:
        # data from Glotaran global fitting:
        irf1=fitRes['IRF'][0]
        irf2=fitRes['IRF'][1]
        dispersion=fitRes['Parmu']
        centralWL=fitRes['centralWL']
        dataset = correct_chirp(dataset,irf2,irf1,dispersion,centralWL)
    return dataset


    
def fit_chirp_dispersion(fitCoeff_data,dispersion_deg=3,centralWL=0,stats=False,GTA_form=True):
    """
    Calculate the polynomial coefficients approximation of the dispersion curve
    inputs:
    fitCoeff_data : pd.DataFrame/pd.Series (wavelength indexed, spectrally resolved values of t_0 from kinetic fitting (Using e.g. Surface Explorer)
    dispersion_deg : degree of the dispersion polynom (default: 3)
    centralWL : wavelength of the origin of the polynom expansion (in nm)
    stats : output the statistical error estimates for the polynomial coefficients
    GTA_form : output the polynomial expansion in order and units used by Glotaran/TIMP (starting from degree 1 upwards).
        Also the units of the expansion coefficients are in powers of nm/100, strangely. 
        Otherwise it starts from highest degree down to 1, and units are in powers of 1/nm.
    
    outputs:
    t0 : zeroth-order polynomial expansion coefficient - time shift calculated at centralWL
    dispersion : coefficent of polynomial expansion of order 1 and higher (units depend on the parameter GTA_form)
    t0_var : statistical variance of t0 
    dispersion_var : statistical variances of dispersion (in units and order specified by GTA_form parameter
    """
    x = fitCoeff_data.index.values - centralWL
    y = fitCoeff_data['t_0(ps)']
    p,V = np.polyfit(x,y, deg=dispersion_deg, cov=True)
    t0 = p[-1]
    variances = np.diagonal(V)
    t0_var = variances[-1]
       
    if GTA_form: # reorganize the dispersion polynomial coefficients (order and units) in the way Glotaran (TIMP) program computes them
        dispersion = []
        dispersion_var = []
        len_p = len(p[:-1])
        for i in range(len_p):
            j=len_p-i-1
            dispersion.append( p[j]*(100**(i+1)) )
            dispersion_var.append (variances[j]*(100**(i+1)) )
    else:
        dispersion = p[:-1]
        dispersion_var = variances[:-1]
    if stats:
        return (t0,dispersion,t0_var,dispersion_var)
    else:
        return (t0,dispersion)

    
def smOOth (dataset,thrshld=1e-3,rm_span=6,ewma_span=6):
    #from pandas.stats.moments import ewma
    """
    Returns smoothed dataset. performs the smoothing with various procedures.
    Intended to smooth-out the noisy data form OceanOptics fiber optic spectrometer
    Inputs:
    dataset: pandas.DataFrame with wavelenghts in index and spectra in individual columns
    filter paramters: Setting a parameter to 0 disables the filter multiple filters are possible, each applied on the result of the preceding one.
    rm_span > 0, default 6 : Rolling median filter. 
    thrshld > 0 default 1e-3 : Reject datapoints with differentials above some threshold,\
              and remove zeros in the input dataset
    ewma_span > 0, default 6 : Exponentially weighter moving average fitler
    """
    smooth = pd.DataFrame(index=dataset.index, columns=dataset.columns, data=dataset.values,copy=True)
    for column in smooth.columns:
        smooth_col = dataset.loc[:,column]
        # rolling median filter
        if rm_span > 0:
            smooth_col = smooth_col.rolling(rm_span, center=True).median()
        # point-to-point variation filter
        if thrshld > 0:
            smooth_col = smooth_col.where(abs(smooth_col.diff())<thrshld ).dropna()
            smooth_col = smooth_col.where(smooth_col!=0 )
        # Exponentially-weighted moving average
        if ewma_span > 0:
            smooth_col = smooth_col.ewm(span=ewma_span,ignore_na=True).mean()
        #assemble the final
        smooth_col.name = column
        smooth[column]=smooth_col    
    return smooth

def nOrm (dataset,norm_type='min'):
    """
    Returns normalized dataset. performs the normalization with various procedures.
    Inputs:
    dataset: pandas.DataFrame with wavelenghts in index and spectra in individual columns
    paramters: Setting a parameter to 0 disables the filter multiple filters are possible, each applied on the result of the preceding one.
        norm_type: 'min' / 'max'  - in each column separatelly, find the minimum/maximum value and divide column with absolute
        'min'  - in each column separatelly, find the minimum value and divide column with absolute
    """
    norm = pd.DataFrame(index=dataset.index, columns=dataset.columns, data=dataset.values,copy=True)
    for column in norm.columns:
        # rolling median filter
        if norm_type  == 'min':
            #min_value = norm.min()
            norm[column] /= abs(norm[column].min())
        if norm_type  == 'max':
            #min_value = norm.min()
            norm[column] /= abs(norm[column].max())
    return norm
    
def unchirp_fitResult (summary, **kwargs):
    """
    Calculates the interpolated dataset that corrects for chirp parameters found by glotaran analysis
    """
    #if centralWL != None:
    #    summary['centralWL']=centralWL
    # load the detail TA characteristics from the saved result dataset:        
    detail = summary['detail']
    centralWL = float(detail['centralWL'])
    if 'dataset' in kwargs:
        traces = kwargs['dataset']
        # correction for chirp
        traces_chirp = correct_chirp(traces,
                                     IRF=summary['IRF'][1],
                                     t0=summary['IRF'][0],
                                     dispersion=summary['Parmu'],
                                     centralWL=centralWL)
        summary['traces_chirp'] = traces_chirp
    else:    
        pass
    return summary
    
def selectRepSpecTA (dataset, title=None, wl_range=[], wl_ex=[], time_idxs=[],times=[],trngs=[],**kwargs):
    '''
    Should  
        - select representative spectra at selected time delays
        - exclude region of pump light 
        - output the datasets of the representative curves as dataframes
    '''
    repre_spectra = dataset; # in case of no restrictions, plot the whole dataset
    legend_title='delay times (ps)'
    legend=None
    if len (wl_range) == 2:
        repre_spectra = repre_spectra[wl_range[0]:wl_range[1]];
    if len (wl_ex) == 2:
        repre_spectra.loc[wl_ex[0]:wl_ex[1],:] = np.nan
    if len(time_idxs) > 0:
        repre_spectra = repre_spectra.iloc[:,time_idxs];
        legend = [('%.1e' % x) for x in repre_spectra.columns.values]
    if len(times) > 0: # in case the actual times and ranges are specified instead of the indexes, one should search for the times in given ranges
        rspecs = pd.DataFrame(index=repre_spectra.index)
        for time in times:
            rspec = repre_spectra.loc[:,0:time].iloc[:,-1]; # select only the last spectrum from range 0:time!
            if len(rspec.shape) >1:
                if rspec.shape[1]>1:
                    rspec = rspec.mean(axis=1);
            rspecs[time] = rspec;
        repre_spectra = rspecs;
        legend = [('%.1e' % x) for x in repre_spectra.columns.values]
    if len(trngs) > 0: # in case the actual times and ranges are specified instead of the indexes, one should search for the times in given ranges
        rspecs = pd.DataFrame(index=repre_spectra.index)
        for trange in trngs:
            rspec = repre_spectra.loc[:,trange[0]:trange[1]]; # select the given time range
            if len(rspec.shape) >1:
                if rspec.shape[1]>1:
                    rspec = rspec.mean(axis=1); #average the time range
            rspecs[tuple(trange)] = rspec;
        repre_spectra = rspecs;
        legend = [('%s' % str (x)).strip('()').replace(', ','$-$') for x in repre_spectra.columns.values]
        legend_title='delay intrvls.(ps)'
    if 'lgnd_prfx' in kwargs.keys():
        lgnd_prfx=kwargs.pop('lgnd_prfx')
        legend = [('%s; %s' % (lgnd_prfx, leg_text)) for leg_text in legend]        
    repre_spectra.columns = legend
    # do the plotting
    if 'ax' in kwargs.keys():
        ax=kwargs['ax']
        repre_spectra.plot(title=title,**kwargs)
        ax.axhline(linewidth=0.6,color='black')
        ax.set_xlabel('Wavelength (nm)');
        ax.legend(legend,title=legend_title, loc='best',frameon=False);
        ax.set_ylabel('$\Delta$A');
    return repre_spectra;
    
def selectRepKinsTA (dataset, title=None, time_range=[], wl_idxs=[],wlngths=[],bndwidth=1.7,\
                   lgnd_prfx='',**kwargs):
    '''
    TODO: NOT completed yet
    Should select 
        - representative kinetics at selected wavelength channels
    Additionally 
        - output the datasets of the representative curves as dataframes
    '''
    repre_kins = dataset.T # in case of no restrictions, plot the whole dataset
    if len (time_range) == 2:
        repre_kins = repre_kins.loc[time_range[0]:time_range[1]];
    if len(wl_idxs) > 0:
        repre_kins = repre_kins.iloc[:,wl_idxs]
        bndwidth = 0;
    if len(wlngths) > 0: # in case the actual wavelengths are specified instead of the indexes, one should search for the wls in range of bandwidth
        rkins = pd.DataFrame(index=repre_kins.index)
        try:
            if len(bndwidth) == len(wlngths):
                wl_bnds = {wlngths[idx] : bndwidth[idx] for idx in range(len(wlngths))}
        except:
            wl_bnds = {wlngths[idx] : bndwidth for idx in range(len(wlngths))}
        for wl in wl_bnds:
            wl_low = wl - wl_bnds[wl];
            wl_hi = wl + wl_bnds[wl];                
            rkin = repre_kins.loc[:,wl_low:wl_hi];
            if rkin.shape[1]>1:
                rkin = rkin.mean(axis=1);
            rkins[wl] = rkin;
        repre_kins = rkins;
    #generate new legends
    legend_texts = [('%s %d' % (lgnd_prfx, x)) for x in repre_kins.columns ];   
    repre_kins.columns = legend_texts
    # get texts from previous legend, and append them before the current legend texts
    
    # do the plotting, optionally
    if 'ax' in kwargs.keys():        
        ax=kwargs['ax']
        previous_legend = ax.get_legend()
        if previous_legend != None:
            legend_texts = [txt.get_text() for txt in previous_legend.get_texts()] + legend_texts
        if 'logx' in kwargs.keys():
            repre_kins.plot(title=title,**kwargs)
        else:
            repre_kins.plot(title=title,logx=True, **kwargs)
        ax.set_xlabel('Time (ps)');
        ax.set_ylabel('$\Delta$A');
        ax.axhline(linewidth=0.6,color='black')
        # fix the legend
        new_legend_handles, new_labels = ax.get_legend_handles_labels()
        ax.legend(new_legend_handles, legend_texts, title=('nm ± %s' % bndwidth),loc='best',frameon=False); 
        #ax.axis([1e-2,1e+4,-25e-3,2e-2]);
    return repre_kins;
    

 
def average_scans (test_TA,idxs=None) :
    mean_TA = None
    count = 0
    if idxs != None:
        for scan in [('scan%d' % idx) for idx in idxs]:
            if type(mean_TA) == type(None):
                mean_TA = test_TA[scan].values
            else:
                mean_TA += test_TA[scan].values
            count += 1
    else:
        for scan in test_TA.items():
            if type(mean_TA) == type(None):
                mean_TA = test_TA[scan].values
            else:
                mean_TA += test_TA[scan].values
            count += 1
    mean_TA /= count
    mean_TA = pd.DataFrame(index=test_TA[scan].index,
                                    columns=test_TA[scan].columns,
                                    data=mean_TA)
    mean_TA = mean_TA.dropna(axis=0)
    return mean_TA



