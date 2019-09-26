# file: 30-plotting.py
# contains scripts for plottign of various spectral data
# written by David Rais, 2016, david.rais@atlas.cz

def plotcmap (dataset,title=None):
    '''
    Plots the color map out of specified transient absorption dataset
    TODO:
        - active click-able axes for browsing the dataset
    '''
    pp.contourf(dataset.T,30, alpha=.75, logy=True)
    pp.cbar = pp.colorbar()
    pp.cbar.solids.set_edgecolor("face")
    pp.xlabel('wavelength index')
    pp.ylabel('time index')
    if title is not None:
        pp.title(title)
    return

def plot_scans (filenames,datalabels,wlngths,bndwidth,trngs ):
    """
    Plot the mean value of each scan, calculated for given wavelength interval and at given time range,
    in order to verify the stability of the TA signal between the individual scans
    """
    test_TA = dict()

    for (fname,label) in zip(filenames,datalabels):
        if 'Representative' not in fname:
            print (fname,label)
            test_TA[label] = load_TA(fname)   

    for time_rng in trngs:
        fig=pp.figure()
        ax = pp.axes()

        means = pd.DataFrame(index=sorted(test_TA.keys(),key=getScanNum), columns=wlngths)
        names = []

        for scan in test_TA.items():
            RKins=selectRepKinsTA(scan[1], wlngths=wlngths, bndwidth=bndwidth, lgnd_prfx=scan[0])
            names.append(scan[0])
            means.loc[scan[0],:] = RKins.loc[time_rng[0]:time_rng[1],:].mean(axis=0).values
            #print(RKins.loc[time_rng[0]:time_rng[1],:].mean(axis=0).T)
            #means.append(mean)

        means.plot(ax=ax)
        ax.set_title('%s ps'%list(time_rng))
    return test_TA

def plotRepSpecTA (dataset, title=None, wl_range=[], wl_ex=[], time_idxs=[],times=[],trngs=[],**kwargs):
    '''
    TODO: NOT completed yet
    Should plot 
        - representative spectra at selected time delays
        - exclude region of pump light
        - representative kinetics at selected wavelength channels
    Additionally 
        - subtract the baseline by averaging few first spectra (NOT implemented yet)
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
                timeLabel = rspec.columns[0]
                if rspec.shape[1]>1:                    
                    rspec = rspec.mean(axis=1);
            else:
                timeLabel = rspec.name
            rspecs[timeLabel] = rspec;
        repre_spectra = rspecs;
        legend = [('%.1e' % x) for x in repre_spectra.columns.values]
    if len(trngs) > 0: # in case the actual times and ranges are specified instead of the indexes, one should search for the times in given ranges
        rspecs = pd.DataFrame(index=repre_spectra.index)
        for trange in trngs:
            rspec = repre_spectra.loc[:,trange[0]:trange[1]]; # select the given time range
            if len(rspec.shape) >1:
                timeLabel = (rspec.columns[0],rspec.columns[-1])
                if rspec.shape[1]>1:
                    rspec = rspec.mean(axis=1); #average the time range
            else:
                timeLabel = rspec.name
            rspecs[tuple(trange)] = rspec;
        repre_spectra = rspecs;
        legend = [('%s' % str (x)).strip('()').replace(', ','$-$') for x in repre_spectra.columns.values]
        legend_title='delay intrvls.(ps)'
    # do the plotting
    if 'ax' in kwargs.keys():
        pass
    else:
        kwargs['ax'] = pp.axes()
    ax=kwargs['ax']
    repre_spectra.plot(title=title,**kwargs)
    ax.axhline(linewidth=0.6,color='black')
    ax.set_xlabel('Wavelength (nm)');
    ax.legend(legend,title=legend_title, loc='best',frameon=False);
    ax.set_ylabel('$\Delta$A');
    #ax.axis([400,800,-15e-3,5e-3]);
    return repre_spectra;

def plotRepKinsTA (dataset, title=None, time_range=[], wl_idxs=[],wlnghts=[],bndwidth=1.7,\
                   lgnd_prfx='',**kwargs):
    '''
    TODO: NOT completed yet
    Should plot 
        - representative kinetics at selected wavelength channels
    Additionally 
        - subtract the baseline by averaging few first spectra (NOT implemented yet)
        - output the datasets of the representative curves as dataframes
    '''
    repre_kins = dataset.T # in case of no restrictions, plot the whole dataset
    if len (time_range) == 2:
        repre_kins = repre_kins.loc[time_range[0]:time_range[1]];
    if len(wl_idxs) > 0:
        repre_kins = repre_kins.iloc[:,wl_idxs]
        bndwidth = 0;
    if len(wlnghts) > 0: # in case the actual wavelengths are specified instead of the indexes, one should search for the wls in range of bandwidth
        rkins = pd.DataFrame(index=repre_kins.index)
        try:
            if len(bndwidth) == len(wlnghts):
                wl_bnds = {wlnghts[idx] : bndwidth[idx] for idx in range(len(wlnghts))}
        except:
            wl_bnds = {wlnghts[idx] : bndwidth for idx in range(len(wlnghts))}
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
    
    #plot new data
    # do the plotting
    # do the plotting
    if 'ax' in kwargs.keys():
        pass
    else:
        kwargs['ax'] = pp.axes()
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
    
def plot_SADS (datasets, ax=None, styles=None):
    if ax is None:
        ax = pp.axes()
    if styles == None:
        styles = []
        for count in range(20):
            styles.append('.')
    if type(datasets) == list:
        for dataset in datasets:
            dataset.plot(ax=ax,style=styles)
    else:
        datasets.plot(ax=ax,style=styles)
    pp.xlabel ('Wavelength (nm)')
    pp.ylabel ('SADS')
    pp.legend(loc='best',title='lifetime (ps)')
    return

def plot_KINS (datasets, ax=None):
    if ax is None:
        ax = pp.axes()
    if type(datasets) == list:
        for dataset in datasets:
            dataset.plot(logx=True,ax=ax)
    else:
        datasets.plot(logx=True,ax=ax)
    pp.xlabel ('Time (ps)')
    pp.ylabel ('Fraction')
    pp.legend(loc='lower right',title='lifetime (ps)')
    return
    
def TA_scans_evaluate(data_dir='..',common_prefix='',common_suffix='.csv',ident_str='_scan',wlngths=[611],bndwidth=20,time_rngs=[[2.0,3.0],[100,200]],ax=None,nrows=513):
    """
    Loads the individual scans of the delay line of the TA experiment into a dictionary, based on the common part of the filename
    
    returns the dictionary with individual scans, 
    plots the signal from each scan 
    """
    ##########################3 select data
    files = os.listdir(data_dir)
    filenames,datalabels =\
    files_select(files,
                 common_prefix=common_prefix,
                 common_suffix=common_suffix, ident_str=ident_str)
    TA_scans=dict()
    for (fname,label) in zip(filenames,datalabels):
        try:
            TA_scans[label] = load_TA(data_dir+fname,nrows=nrows)
        except: 
            pass
    print('Loaded %d files' % len(TA_scans))

    #########################33 Plotting the individual scan signals
    
    scans_dict = TA_scans
    means = pd.DataFrame(index=sorted(scans_dict.keys()), columns=wlngths)
    scans = scans_dict.items()
    names = []
    for time_rng in time_rngs:
        for scan in scans:
            RKins=selectRepKinsTA(scan[1], wlngths=wlngths, bndwidth=bndwidth, lgnd_prfx=scan[0])
            names.append(scan[0])
            means.loc[scan[0],:] = RKins.loc[time_rng[0]:time_rng[1],:].mean(axis=0).values
            #print(RKins.loc[time_rng[0]:time_rng[1],:].mean(axis=0).T)
            #means.append(mean)
        if type(ax)!= type(pp.axes()):
            ax = pp.axes()
        means.plot(ax=ax)
    return TA_scans