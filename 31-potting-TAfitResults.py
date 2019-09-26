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
                if rspec.shape[1]>1:
                    rspec = rspec.mean(axis=1); #average the time range
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

def plotRepKinsTA (dataset, title=None, time_range=[], wl_idxs=[], wlngths=[], bndwidth=1.7,\
                   lgnd_prfx='',**kwargs):
    '''
    TODO: NOT completed yet
    Should plot 
        - representative kinetics at selected wavelength channels
    Additionally 
        - subtract the baseline by averaging few first spectra (NOT implemented yet)
        - output the datasets of the representative curves as dataframes
    '''
    wlnghts=wlngths
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
    if 'ax' in kwargs.keys():
        pass
    else:
        kwargs['ax'] = pp.axes()
    ax=kwargs['ax']
    previous_legend = ax.get_legend()
    if previous_legend != None:
        legend_texts = [txt.get_text() for txt in previous_legend.get_texts()] + legend_texts
    if 'logx' in kwargs.keys():
        repre_kins.plot(**kwargs)
    else:
        repre_kins.plot(logx=True, **kwargs)
    ax.set_xlabel('Time (ps)');
    ax.set_ylabel('$\Delta$A');
    ax.axhline(linewidth=0.6,color='black')
    # fix the legend
    new_legend_handles, new_labels = ax.get_legend_handles_labels()
    ax.legend(new_legend_handles, legend_texts, title=('nm ± %s' % bndwidth),loc='best',frameon=False); 
    #ax.axis([1e-2,1e+4,-25e-3,2e-2]);
    return repre_kins;
    
def plot_SADS (datasets, styles=None,**kwargs):
    if 'ax' in kwargs.keys():
        pass
    else:
        kwargs['ax'] = pp.axes()
    ax=kwargs['ax']
    if styles == None:
        styles = []
        for count in range(20):
            styles.append('.')
    if type(datasets) == list:
        for dataset in datasets:
            dataset.plot(style=styles, **kwargs)
    else:
        datasets.plot(style=styles, **kwargs)
    pp.xlabel ('Wavelength (nm)')
    pp.ylabel ('SADS')
    pp.legend(loc='best',title='lifetime (ps)')
    return    
    
def plot_RSPEC_fits(traces,fittedTraces,colors='kbgr',fitSty='-',expSty='o',**kwargs ):
    if 'ax' in kwargs.keys():
            pass
    else:
        kwargs['ax'] = pp.axes()
    ax=kwargs['ax']
    fitStyleList = [('-%s' % s) for s in colors]
    expStyleList = [('o%s' % s) for s in colors]
    # experimental traces
    RspecExp = plotRepSpecTA(traces, style=expStyleList, alpha=0.2, **kwargs)
    # fitted curves
    handles,labels = ax.get_legend_handles_labels()
    labels = [('%.2g' % float(x)) for x in labels]
    RspecFit = plotRepSpecTA(fittedTraces, style=fitStyleList, alpha=1.0, **kwargs)
    ax.set_xlim(wl_range)
    ax.legend(handles,labels,title='delay time (ps):')
    #ylim = ax.get_ylim()
    return
    
def plot_RKINS_fits (traces,fittedTraces, time_range=[], 
                     wl_idxs=[], wlngths=[], bndwidth=1.7,
                     colors='kbgr',fitSty='-',expSty='o',**kwargs):
    if 'ax' in kwargs.keys():
        pass
    else:
        kwargs['ax'] = pp.axes()
    ax=kwargs['ax']
    fitStyleList = [('-%s' % s) for s in colors]
    expStyleList = [('o%s' % s) for s in colors]
    
    # experimental traces
    RspecExp = plotRepKinsTA(traces, time_range=time_range, 
                             wl_idxs=wl_idxs, wlngths=wlngths, 
                             bndwidth=bndwidth,style=expStyleList, 
                             alpha=0.2, *kwargs)
    # fitted curves
    handles,labels = ax.get_legend_handles_labels()
    labels = [('%d' % float(x)) for x in labels]
    RspecFit = plotRepKinsTA(fittedTraces, time_range=time_range, 
                             wl_idxs=wl_idxs, wlngths=wlngths, 
                             bndwidth=bndwidth, style=fitStyleList,
                             alpha=1.0, **kwargs)
    #ax.set_xlim(wl_range)
    ax.legend(handles,labels,title='probe wls. (nm):')
    #ax.set_ylim(ylim)
    return

def plot_RKINS_fits_linlog (traces,fittedTraces, time_range=[], 
                             wl_idxs=[], wlngths=[], bndwidth=1.7,
                             colors='kbgr',fitSty='-',expSty='o',**kwargs):
    if 'axs' in kwargs.keys():
        pass
    else:
        n = 1; m = 5;
        gs = gridspec.GridSpec(1,2, width_ratios = [n,m])
        #pp.figure(figsize=(10,8))
        ax = pp.subplot(gs[0,0])
        ax2 = pp.subplot(gs[0,1], sharey = ax)
        pp.setp(ax2.get_yticklabels(), visible=False)
        pp.subplots_adjust(wspace = 0.05)
    
        kwargs['axs'] = (gs,ax,ax2)    
    gs,ax,ax2=kwargs.pop('axs')
    fitStyleList = [('-%s' % s) for s in colors]
    expStyleList = [('o%s' % s) for s in colors]
    #plotting before the break
    RspecExp = plotRepKinsTA(traces.loc[:,:t_break],
                             time_range=time_range, 
                             wl_idxs=wl_idxs, wlngths=wlngths, 
                             bndwidth=bndwidth,style=expStyleList, 
                             alpha=0.2, ax=ax,logx=False,
                             *kwargs)
    handles,labels = ax.get_legend_handles_labels()
    labels = [('%d' % float(x)) for x in labels]
    RspecFit = plotRepKinsTA(fittedTraces.loc[:,:t_break], 
                             time_range=time_range, 
                             wl_idxs=wl_idxs, wlngths=wlngths, 
                             bndwidth=bndwidth,
                             ax=ax,
                             logx=False,
                             style=fitStyleList, 
                             alpha=1.0)
    ax.legend().set_visible(False)
    ax.set_xlim(time_rng[0],t_break)
    #plotting after the break
    RspecExp = plotRepKinsTA(traces.loc[:,t_break:],
                             time_range=time_range, 
                             wl_idxs=wl_idxs, wlngths=wlngths, 
                             bndwidth=bndwidth,
                             ax=ax2, 
                             style=expStyleList, 
                             alpha=0.2)
    RspecFit = plotRepKinsTA(fittedTraces.loc[:,t_break:], 
                             time_range=time_range, 
                             wl_idxs=wl_idxs, wlngths=wlngths, 
                             bndwidth=bndwidth,
                             ax=ax2,
                             style=fitStyleList, 
                             alpha=1.0)
    ax2.set_xlim(t_break+1,time_rng[1])
    
    # hide the spines between ax and ax2
    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax.tick_params(labeltop='off') # don't put tick labels at the top
    ax2.yaxis.tick_right()
    ax2.get_yaxis().set_visible(False)
    d = .015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    n,m = gs.get_geometry()
    #n,m = gs.get_width_ratios()
    on = (n+m)/n; om = (n+m)/m;
    ax.plot((1-d*on,1+d*on),(-d,d), **kwargs) # bottom-left diagonal
    ax.plot((1-d*on,1+d*on),(1-d,1+d), **kwargs) # top-left diagonal
    kwargs.update(transform=ax2.transAxes) # switch to the bottom axes
    ax2.plot((-d*om,d*om),(-d,d), **kwargs) # bottom-right diagonal
    ax2.plot((-d*om,d*om),(1-d,1+d), **kwargs) # top-right diagonal
    # removing unwanted x label
    ax.set_xlabel('')
    # the common legend
    ax2.legend(handles,labels,title='probe wls. (nm):',ncol=RspecExp.shape[1])
    #ax.set_ylim(ylim)
    #ax2.set_ylim(ylim)
    return
    
def plot_SADS_w_err(SADS,rates=[],SADS_err=None,rates_stderrs=[],colors='cmyk',markerSty='-',wl_range=[], CIfactor=1.96, **kwargs):
    if 'ax' in kwargs.keys():
        pass
    else:
        kwargs['ax'] = pp.axes()
    ax=kwargs['ax']
    GAStyleList = [('%s%s' % (markerSty,s)) for s in colors]
    # experimental traces
    plot_SADS(SADS, ax=ax, styles=GAStyleList)
    # fitted curves
    handles,labels = ax.get_legend_handles_labels()
    #labels = [('%d' % float(x)) for x in labels]
    # labels formating
    lftmLabels=[]
    if len(rates) == len(SADS.columns):
        lftms,lftmErrs = lftm_from_rate_w_err (rates,rates_stderrs)
        if len(rates_stderrs) == len(rates):
            for lftm in zip(lftms,lftmErrs):
                orders_lftm = np.log10(lftm)
                min_order = (min(orders_lftm))
                #print(min_order)
                if min_order < -2:        
                    lftmLabels.append('%.3f ± %.3f' % lftm)
                else:
                    if min_order >= -2 and min_order < -1:        
                        lftmLabels.append('%.2f ± %.2f' % lftm)
                    else:
                        if min_order >= -1 and min_order < 0:
                            lftmLabels.append('%.1f ± %.1f' % lftm)
                        else:
                            lftmLabels.append('%d ± %d' % lftm)
        else:
            lftms  = [1/k for k in rates]
            lftmLabels = [('%.2g' % tau) for tau in lftms]
    else:
        lftmLabels = SADS.columns
    if type(SADS_err) != type(None):
        # calculate the confidence interval boundaries at 95% level
        SADS_highLim = np.float64(SADS.values) + np.float64(SADS_err.values)*CIfactor
        SADS_lowLim = np.float64(SADS.values) - np.float64(SADS_err.values)*CIfactor
        # fill between - plot of the error bands
        for comp in range(len(SADS.columns)):
            #count=comp+1
            ax.fill_between(SADS.index.values,                    
                            y1=SADS_lowLim[:,comp], 
                            y2=SADS_highLim[:,comp],
                            facecolor=colors[comp], 
                            #facecolor='b',
                            interpolate=True, alpha=0.5)
    else:
        pass
    ax.legend(handles,lftmLabels,title='lifetime (ps):')
    if len(wl_range) > 0:
        ax.hlines(0,wl_range[0],wl_range[1])
        #ax.set_ylim(ylim)
        ax.set_xlim(wl_range)
    return

def plot_frac_linlog (frac, time_range=[], colors='cmyk',markerSty='-',**kwargs):
    if 'axs' in kwargs.keys():
        pass
    else:
        n = 1; m = 5;
        gs = gridspec.GridSpec(1,2, width_ratios = [n,m])
        #pp.figure(figsize=(10,8))
        ax = pp.subplot(gs[0,0])
        ax2 = pp.subplot(gs[0,1], sharey = ax)
        pp.setp(ax2.get_yticklabels(), visible=False)
        pp.subplots_adjust(wspace = 0.05)
    
        kwargs['axs'] = (gs,ax,ax2)
    
    gs,ax,ax2=kwargs.pop('axs')
    # define plot styls
    GAStyleList = [('%s%s' % (markerSty,s)) for s in colors]

    frac.loc[:1].plot(ax=ax, style=GAStyleList,logx=False)
    ax.legend().set_visible(False)
    frac.loc[1:].plot(ax=ax2, style=GAStyleList,logx=True)
    ax2.legend().set_visible(False)

    ax.set_xlim(time_rng[0],t_break)
    ax2.set_xlim(t_break+1,time_rng[1])

    # hide the spines between ax and ax2
    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax.tick_params(labeltop='off') # don't put tick labels at the top
    ax2.yaxis.tick_right()
    ax2.get_yaxis().set_visible(False)

    d = .015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    n,m=gs.get_geometry()
    on = (n+m)/n; om = (n+m)/m;
    ax.plot((1-d*on,1+d*on),(-d,d), **kwargs) # bottom-left diagonal
    ax.plot((1-d*on,1+d*on),(1-d,1+d), **kwargs) # top-left diagonal
    kwargs.update(transform=ax2.transAxes) # switch to the bottom axes
    ax2.plot((-d*om,d*om),(-d,d), **kwargs) # bottom-right diagonal
    ax2.plot((-d*om,d*om),(1-d,1+d), **kwargs) # top-right diagonal

    ax.set_xlabel('')
    ax2.set_xlabel('Time (ps)')
    return

def plot_GA_summary(summary,repreWLS,wl_range=[],time_range=[],bndwidth=1.7,CIfactor=1.96,colors='kbgrcmy',colorsGA='cmykbgr'):
    # summary parts
    detail = summary['detail']
    bndWidth = 1.7 # bandwidth for spectral averaging ( bndWidth=<1.7 nm usually means no averaging)
    rates = summary['kVect']
    rates_stderrs = summary['kVect_stderr']
    lifetimes = [1/k for k in rates] 
    repreT = lifetimes
    if lifetimes[0] > 1.5:
        repreT.insert(0,0.5)
    traces = detail['traces']
    fittedTraces = detail['fittedTraces']
    SADS = detail['SADS']
    if 'SADS_stderr' in detail.keys():
        SADS_err = detail['SADS_stderr']
    else:
        SADS_err = None
    frac = detail['frac']
    # figure layout definition
    fig = pp.figure(figsize=(12,8))
    gsEXP= gridspec.GridSpec(2, 2) # master grid
    gsEXPkin = gridspec.GridSpecFromSubplotSpec(1, 6, subplot_spec=gsEXP[0,1]) # experimtal
    gsFrac = gridspec.GridSpecFromSubplotSpec(1, 6, subplot_spec=gsEXP[1,1]) # SADS+fracs
    RspecAX = fig.add_subplot(gsEXP[0,0])
    RkinsLinAX = fig.add_subplot(gsEXPkin[:,0], sharey=RspecAX)
    RkinsLogAX = fig.add_subplot(gsEXPkin[:,1:], sharey=RspecAX)
    fracLinAX = fig.add_subplot(gsFrac[:,0], sharex=RkinsLinAX)
    fracLogAX = fig.add_subplot(gsFrac[:,1:], sharex=RkinsLogAX, sharey=fracLinAX)
    sadsAX = fig.add_subplot(gsEXP[1,0], sharex=RspecAX,sharey=RspecAX)
    gsEXP.update(wspace=0.15, hspace=0.1)
    # plotting
    plot_RSPEC_fits(traces,fittedTraces,wl_range=wl_range, times=repreT, ax=RspecAX,colors=colors)
    ylim = RspecAX.get_ylim()
    plot_RKINS_fits_linlog (traces,fittedTraces,colors=colors, 
                            wlngths=repreWLS, bndwidth=bndWidth,
                            fitSty='-',expSty='o',
                            axs=(gsEXPkin,RkinsLinAX,RkinsLogAX))
    plot_SADS_w_err(SADS,rates,SADS_err,rates_stderrs,ax=sadsAX,colors=colorsGA,CIfactor=CIfactor)
    plot_frac_linlog(frac, time_range=time_range, axs=(gsFrac,fracLinAX,fracLogAX),colors=colorsGA)
    #adjust y scale
    RspecAX.set_ylim(ylim)
    RspecAX.set_xlim(wl_range)
    # adjust plot texts
    RkinsLogAX.tick_params(axis='both', left='off', top='off',  labelleft='off', labeltop='off', labelright='off', labelbottom='off')
    RkinsLogAX.set_xlabel('')
    ticks=fracLinAX.get_xticks()
    ticks = [('%.0g' %ticks[i]) for i in range(len(ticks))]
    fracLinAX.set_xticklabels(ticks,visible=True)
    fracLinAX.set_ylabel('Fraction')
    RkinsLinAX.tick_params(axis='both', left='on', top='off',  labelleft='off', labeltop='off', labelright='off', labelbottom='off')
    RkinsLinAX.set_xlabel('')
    return fig