# -*- coding: utf-8 -*-
"""
Created on Fri May 6 15:46 2016 

Yield a class that accepts differential equations system
With integrated solver and 
fitting routines

@author: David Rais, david.rais@atlas.cz

"""
import numpy as np
import sympy as syp
from scipy import integrate, optimize
#from sympy.utilities.autowrap import ufuncify
#from strToCode import strToCode

t = syp.symbols ('t') # t, is the symbol for time in differential Eqs.
x1 = syp.Function ('x1')(t)
x2 = syp.Function ('x2')(t)
promenne = (x1,x2)  
dx1 = syp.diff(x1,t)
dx2 = syp.diff(x2,t)
    
class modeleq:
    def __init__(self,symODE,symVar,t,symPar,y0,initGuess):
        """
        Construct the model ODEs with integrator for fitting procedure
        of paramters
        make optimization of parameters (model fitting)
        report fitting result
        """                  
        self.symODE = symODE
#        self.symVar = symVar
#        self.t = t
#        self.symPar = symPar       
        self.y0 = y0
        self.initGuess = initGuess # initial guess of fitting parameter values
        self.ODE_system_lambda = syp.lambdify ((symVar,t,symPar),symODE,'numpy')
        #self.ODE_system_lambda = ufuncify ([symVar,t,symPar],symODE)
        self.n = 0 # iteration counter
        self.fitresult = None
        def ODE_system(y,t,k):
            """
            assign variables, time and optional parameters 
            to differential equation system
            """
            return self.ODE_system_lambda(list(y),t,(list(k)))
        self.f = ODE_system
        return
    def my_ls_func(self,x,y0,teta):
        """definition of function for fit
            x gives evaluation points,
            teta is an array of parameters to be varied for fit"""
        # create an alias to f which passes the optional params    
        f2 = lambda y,t: self.f(y, t, teta)
        # calculate ode solution, return values for each entry of "x"
        r = integrate.odeint(f2,y0,x)
        return r
    def f_resid(self,x_data, y_data,y0,p):
        """ function to pass to optimize.leastsq
            The routine will square and sum the values returned by 
            this function""" 
        solution = self.my_ls_func(x_data,y0,p)
        if y_data.shape==solution.shape:
            residuals = y_data-solution
        else:
            residuals =y_data-np.sum(solution,axis=1)
        return residuals
    def fitScore(self,x_data,y_data,y0,p):
        return np.linalg.norm(self.f_resid(x_data, y_data,y0,p))
    def globalFitScore(self,x_data_global,y_data_global,y0_global,errWeights,p_shared):
        """
        calculate the fitting error sum for the respective datasets with shared parameters        
        """
        fiterr = 0
        #counter = 0
        for idx in range(len(y_data_global)):
            fiterr += self.fitScore(x_data_global[idx],
                                    y_data_global[idx],
                                    y0_global[idx],
                                    p_shared) *errWeights[idx]
            #counter += 1
        return fiterr
    def counter(self, arg):
        self.n +=1
        return
    def minimize (self, x_data_global, y_data_global, y0_global, errWeigths_global, maxIter):
        "wrapper function for the fitting algorithm"
        
        self.fitresult = optimize.minimize (lambda p:self.globalFitScore(x_data_global,
                                                             y_data_global,
                                                             y0_global,
                                                             errWeigths_global,
                                                             p),
                               x0=np.array(self.initGuess), 
                               args=(),
                               method='Nelder-Mead', 
                               jac=None, 
                               hess=None, 
                               hessp=None, 
                               bounds=None, 
                               constraints=(), 
                               tol=None, 
                               #callback=self.counter, 
                               options={'maxiter': maxIter, 
                                        'disp':True,
                                        #'ftol':0.0003,
                                        'maxfev':10000
                                        })
                                        
        return self.fitresult

#(fitresult.x,n,fitresult.message) = optimize.fmin_tnc(lambda p:model.fitScore(x_data,y_data,y0,p),
#                              x0=np.array(list(initialGuess.values())),
#                              fprime=None, 
#                              args=(), 
#                              approx_grad=True, 
#                              bounds=None,
#                              epsilon=1e-08, 
#                              scale=None, 
#                              offset=None, 
#                              messages=15, 
#                              maxCGit=-1, 
#                              maxfun=None, 
#                              eta=-1, 
#                              stepmx=0, 
#                              accuracy=1e-10, 
#                              fmin=0, 
#                              ftol=1.0e-7,
#                              #xtol=1e-7, 
#                              #pgtol=1e-7, 
#                              rescale=-1, 
#                              disp=None, 
#                              #callback=counter,
#                              )


# -*- coding: utf-8 -*-
"""
Created on Fri May 6 15:46 2016
This module should contain functions to output the result of fitting
@author: rais
"""
import pylab as pp
import numpy as np
import pandas as pd

pptitle = ''

def printFitParams(symbolsStringsList,fitresult):
    """ print the best fit parameters 
    and assemble plot title based on the values """   
    global pptitle
    for name in symbolsStringsList:
        fitParamEq = name + '\t=\t' + ('%.4e' % fitresult.x[symbolsStringsList.index(name)])
        print (fitParamEq)
        pptitle += name + (' = %.2e' % fitresult.x[symbolsStringsList.index(name)])
        if symbolsStringsList.index(name)!=0 and not( np.remainder(symbolsStringsList.index(name),3)): pptitle+='\n'
        else: pptitle += '; '
    return 

def printFitMessages(x_data_global, y_data_global, y0_global, errWeigths_global, model, fitresult):
    #print ('Total iterations', model.n)
    print('The fit exit message: ',model.fitresult.message)
    fitScore_global = model.globalFitScore(x_data_global,
                                           y_data_global,
                                           y0_global,
                                           errWeigths_global,
                                           fitresult.x)
    print (('std deviation of the GLOBAL fit is %.5f' % fitScore_global))

def getBestFitCurves(x_data_global, y_data_global, y0_global, model, fitresult):
    """ assemble the  
    fit curves and residuals for each of the datasets"""
    #global pptitle
    dataAndBestFits = list() # list of all the fitted datasets
    fitErrors = list() # registry for the individual fit errors, no weighting applied (standard deviations)
    for idx in range (len(x_data_global)):
        # initialize the output dataset
        if len(y_data_global[idx].shape) > 1:
            bestFitData = pd.DataFrame(index=x_data_global[idx],data=y_data_global[idx]) 
        if len(y_data_global[idx].shape) == 1:
            bestFitData = pd.DataFrame(index=x_data_global[idx])        
            bestFitData['exprnmt' ] = y_data_global[idx] # add experimental values
        bestFitCurves = model.my_ls_func(x_data_global[idx],y0_global[idx],fitresult.x)
        bestFitCurves_sum = np.sum(bestFitCurves,axis=1)
        #bestFitData['mod_sum' ] = bestFitCurves_sum #
        for i in range(bestFitCurves.shape[1]): # add calculated data to the output dataset
            bestFitData[('mod_%i' % i )] = bestFitCurves[:,i] # add modelled values
        fitScore = model.fitScore(x_data_global[idx],y_data_global[idx],y0_global[idx],fitresult.x)
        residues = model.f_resid(x_data_global[idx],y_data_global[idx],y0_global[idx],fitresult.x)
        if len(residues.shape) > 1:
            residues = residues.sum(axis=1)
        bestFitData['resid'] = residues # add residues
        #dataAndBestFits.append(bestFitData.copy(deep=True) )
        dataAndBestFits.append(bestFitData.copy() )
        fitErrors.append(fitScore)     
    return dataAndBestFits

def saveFitCurves(filenameBase, fitCurves):
    """ save the fitcurves in the files"""
    for index in range (len(fitCurves)):
        fitCurves[index].to_csv( ('%s_trace#%i.csv' % (filenameBase,index) ), index_label='time')
    return
        


        