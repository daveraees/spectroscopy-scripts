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