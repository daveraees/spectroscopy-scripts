# this file 50-fitting.py
# should contain functions for fitting of transient absorption data
# written by David Rais, 2016, david.rais@atlas.cz


def fitFixedProfiles (experimental,profiles,wl_rng=[0,-1]):
    """
    Function to decompose the spectro-temporal TA dataset into components with time-independent spectral profiles. It fits weight factors of the compoennts to experimental spectra at each delay time 
    expects DataFrame input datasets
    - experimental TA spectro-temporal array
    - custom spectral profiles for the decomposition; any number of them in the interval (2-N).
    - The wavelengths of the experimental and profiles dataset must match each other.
    wl_rng : range of wavelengths for which to make the fitting.
    Outputs results in form of dictionary with the following keys: value pairs
        - 'coef' : dataframe with the fit coefficients and residuals [a0, a1, .., aN, residuals] for each delay time
        - 'fitCurves' : dataframe with reconstructed spectro-temporal TA dataset using the fitted coefficients
    """
    #create data structure for results
    result_summary=dict()
    indexes = [('a%d' % profile_idx) for profile_idx in range(len(profiles.columns))];
    indexes.append('residuals');
    result_summary['coef'] = pd.DataFrame(columns=experimental.columns, index=indexes);
    result_summary['fitCurves'] = pd.DataFrame(columns=experimental.columns, 
                                               index=experimental.loc[wl_rng[0]:wl_rng[1]].index);    
    # do the fitting using linear algebraic least squares fitting function
    for delayTime_idx in range(len(experimental.columns.values)):
        (x,residuals,rank,s) = np.linalg.lstsq(profiles.loc[wl_rng[0]:wl_rng[1]],\
                                               experimental.loc[wl_rng[0]:wl_rng[1]].iloc[:,delayTime_idx])
        result_summary['coef'].iloc[:,delayTime_idx]=np.concatenate ((x, residuals));
        result_summary['fitCurves'].loc[wl_rng[0]:wl_rng[1]].iloc[:,delayTime_idx]=\
        (profiles.loc[wl_rng[0]:wl_rng[1]] * x).sum(axis=1);
    return result_summary;
    
def decompose_Spec (expData,fit_bases,repreWls=[],scaleWls=[],bndWidth=3,ax=None):
    """
    """
    numBases = len(fit_bases.columns.values)
    namesBases = fit_bases.columns.values
    wl_rng = [fit_bases.index.values[0],fit_bases.index.values[-1]] 
    # do the fitting
    fits = fitFixedProfiles (expData,fit_bases,wl_rng)
    
    # do the plotting
    ExpStyles = [('.%s' % colorLetter) for colorLetter in 'bgrcmyk']
    FitStyles = [('-%s' % colorLetter) for colorLetter in 'bgrcmyk']
    if ax == None:
        ax = pp.axes()
    
    expRepreKins = plotRepKinsTA(expData,
              title=None, wlnghts=repreWls, bndwidth=bndWidth, ax=ax, lgnd_prfx='exp', style=ExpStyles,
              alpha=0.25)
    expRepreKins.index.name = 'Time (ps)'
    
    if len (scaleWls) == 1:
        # scaling of the weight factors with the GSB time profile
        wls_scale = scaleWls
        expScaleKin=plotRepKinsTA(expData,
              title=None, wlnghts=wls_scale, bndwidth=bndWidth, ax=ax, lgnd_prfx='exp', style=ExpStyles[len(repreWls)+1],
              alpha=0.25);
    #print(fits['fitCurves'].shape)    
    fitRepreKins = plotRepKinsTA(fits['fitCurves'],
              title=None, wlnghts=repreWls, bndwidth=bndWidth, ax=ax, lgnd_prfx='fit',
              style=FitStyles)
    fitRepreKins.index.name = 'Time (ps)'
    
    if len (scaleWls) == 1:
        dt_rng=fits['coef'].iloc[0:numBases,[0,-1]].T.index.values
        scaleFitresult = fitFixedProfiles(expScaleKin,fits['coef'].iloc[0:numBases,:].T,wl_rng=dt_rng)
        scaleFitresult['fitCurves'].columns = scaleWls
        #scaleFitresult['fitCurves'].plot(logx=True,ax=ax,style=FitStyles[numBases],)
        scaleFitresult['fitCurves'] = plotRepKinsTA(scaleFitresult['fitCurves'].T,
                                                    title=None,  wl_idxs=[0], ax=ax, 
                                                    lgnd_prfx='fit',
                                                    style=FitStyles[len(repreWls)+1])        
        scaled_fit = -scaleFitresult['coef'].iloc[0:numBases].T.values*fits['coef'].iloc[0:numBases,:].T
        scaled_fit['residuals']=fits['coef'].T['residuals']
        #nameColumns = [('$c_{%s}$' % baseName) for baseName in namesBases]
        #scaled_fit.columns = nameColumns.append('residuals')
        fits = scaled_fit / max(scaled_fit.iloc[:,0:2].sum(axis=1))
        pd.DataFrame.sum
        fitRepreKins = pd.concat([fitRepreKins,scaleFitresult['fitCurves']],axis=1)
        expRepreKins = pd.concat([expRepreKins,expScaleKin],axis=1)
    else:
        fits=fits['coef'].T
    
    # join the representative kinetics in single dataframe, with format similar to Glotaran fit results:
    #repreWls.extend(scaleWls) 
    RepreKinsList = list()
    i = 0
    for colName in expRepreKins.columns:
        RepreKinsList.append(expRepreKins.iloc[:,i])
        RepreKinsList.append(fitRepreKins.iloc[:,i])
        residuals = expRepreKins.iloc[:,i] - fitRepreKins.iloc[:,i]       
        residName = 'residual '+ colName.strip("""abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ,.:'"~!@#$%^&*()_+|`\}{}[] """)
        #[int(s) for s in str.split() if s.isdigit()]
        residuals.name = residName
        RepreKinsList.append(residuals)
        i +=1
    RepreKins = pd.concat(RepreKinsList, axis=1)
    
    # add the description of the axes and some nice legend
    ax.legend(loc='best',frameon=False,title=('nm ± %d' % bndWidth))
    #ax.set_xlabel('Time (ps)')
    #ax.set_ylabel('$\Delta A$ (OD)')
    return (fits,RepreKins)

    