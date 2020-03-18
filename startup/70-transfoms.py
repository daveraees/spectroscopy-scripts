def nm2ev (nm, h=4.135667662e-15,c=299792458):
    """
    transformation of the photon wavelength to energy 
    """
    ev = h*c / ( 1e-9 *nm )
    return ev

def TA_to_eV_scale (datasetNM, h=4.135667662e-15,c=299792458):
    """
    Jacobian transformation of the transient absorption spectra to energy scale
    input : pandas DataFrame of spectra in wavelength
    output : pandas DataFrame of spectra in electronvolts
    """
    nm = datasetNM.index.values # expecting pandas dataframe from TA experiment, scaled in nanometers
    energies = nm2ev (nm, h=h,c=c) 
    #TA_weights = h*c * np.power(energies,-2) # due to ratiometric nature of the TA experiment, no need for transformation of the weights
    datasetEV = pd.DataFrame(index=energies,columns=datasetNM.columns)
    for i in range(len(energies)):
        datasetEV.iloc[i,:] = datasetNM.iloc[i,:]
    return datasetEV

def PL_to_eV_scale (datasetNM, h=4.135667662e-15,c=299792458):
    """
    Jacobian transformation of the *photoluminescence* spectra to energy scale
    input : pandas DataFrame of spectra in wavelength
    output : pandas DataFrame of spectra in electronvolts
    """
    nm = datasetNM.index.values # expecting pandas dataframe from PL experiment, scaled in nanometers
    energies = nm2ev (nm, h=h,c=c) 
    TA_weights = h*c * np.power(energies,-2) # here the Jacobian is important
    datasetEV = pd.DataFrame(index=energies,columns=datasetNM.columns)
    for i in range(len(energies)):
        datasetEV.iloc[i,:] = datasetNM.iloc[i,:]*TA_weights[i]
    return datasetEV

def regularize_time (dataset, start, stop, time_step):
    #dataset = residuals_chirp

    timeExp=dataset.columns
    time_step = 0.05 #
    #start = timeExp[0]
    #stop = 10
    num = (stop - start) / time_step + 1
    tstamps = np.linspace(start, stop, num)
    newCurves = pd.DataFrame(index=dataset.index,columns=tstamps)
    #newCurves=deepcopy(dataset)
    for wl in dataset.index:
        #print(type(wl))
        trace = dataset.loc[wl,:].values
        ndTrace = trace.astype(type(1.21654e-5), order='K', casting='unsafe')
        newCurves.loc[wl,:] = np.interp(tstamps,
                                        timeExp,
                                        ndTrace)
    # attribution, plotting, saving
    return newCurves

def fft_TAevolution (dataset):
    # FFT analysis
    #time_rng = [-1.25, stop]
    #wl_rng = [0,700]

    #dataset = residuals_chirp_reg.loc[wl_rng[0]:wl_rng[1],time_rng[0]:time_rng[1]]



    timeExp=dataset.columns
    time_step = (timeExp[-1] - timeExp[0] ) / (timeExp.shape[-1] -1)

    freq = np.fft.fftfreq(timeExp.shape[-1])/time_step # frequency in THz

    SPreal = pd.DataFrame(index=dataset.index,columns=freq)
    SPimag = pd.DataFrame(index=dataset.index,columns=freq)
    #SPmag = pd.DataFrame(index=dataset.index,columns=freq)
    #newCurves=deepcopy(dataset)
    for wl in dataset.index:
        #print(type(wl))
        trace = dataset.loc[wl,:].values
        spec = np.fft.fft(trace)
        SPreal.loc[wl,:] = spec.real
        SPimag.loc[wl,:] = spec.imag

    # removing the rendundant part
    freqHalfLen = int( SPreal.shape[1] / 2)
    SPreal_red = SPreal.iloc[:,:freqHalfLen]
    SPimag_red = SPimag.iloc[:,:freqHalfLen]
    SPmag = (SPreal_red**2 + SPimag_red**2)**0.5
    SPsinPha = SPimag_red / SPmag
    SPpha = SPsinPha.applymap(lambda x: np.arcsin (x))
    return (SPmag,SPpha)

