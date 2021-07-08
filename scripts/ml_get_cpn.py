import sys
import socket
import os
import os.path
from optparse import OptionParser
from collections import defaultdict
import tempfile
import time
import subprocess
import gzip

import scipy as scp
import numpy as np
import tables

import scipy
from scipy.stats import norm,poisson

from math import *
import Bio.Statistics.lowess as biostats

#

#def get_loc(msk,wssd,start,end,chr):

def do_lsqfit(x,y,w,pinit,fitfunc):
    assert x.dtype=='float64'
    assert y.dtype=='float64'
    assert w.dtype=='float64'
    werrfunc = lambda p, x, y, w: (y - fitfunc(p, x)) * w
    out = scipy.optimize.leastsq(werrfunc, pinit,args=(x,y,w), full_output=1)
    return out[0]

def do_lsqfit_thruzero(x,y,w,pinit,fitfunc):
    assert x.dtype=='float64'
    assert y.dtype=='float64'
    assert w.dtype=='float64'
    werrfunc = lambda p, x, y, w: (y - fitfunc(p, x)) * w
    out = scipy.optimize.leastsq(werrfunc, pinit,args=(x,y,w), full_output=1)
    return out[0]

def sum_g(dists,i,group_size,n_single):
    sum = np.zeros(dists[i].shape[0])
    i=i*group_size + n_single
    print("\tmade from:")
    for k in range(group_size):
        print("\t",k+i)    
        sum+= dists[k+i]

    sum /= group_size
    return sum

def make_emissions_DCPT(mus,sigmas,fn_out,card):
    gauss_func = lambda p, x: (1/(p[1]*np.sqrt(2*pi))*np.exp(-((x-p[0])**2)/(2*p[1]**2)))    #takes in [mu, sigma]

    max = 10000
    min_p = 1e-100
    #max = mus[card-1]+3*sigmas[card-1]
    #output = "1\n0\nemissions_DCPT\n1\n"
    #now pringting the 1 parent cardinality, and the cardinality of self 
    #output += "%d %d\n"%(card,max)
    output = ""
    
    dists = []
    for i in range(card):
        x_axis = np.arange(max)    
        p = np.array((mus[i],sigmas[i]))
        discrete_g_fit = gauss_func(p,x_axis)+min_p
        discrete_g_fit = discrete_g_fit/discrete_g_fit.sum()
        dists.append(discrete_g_fit)
    
    grouped_dists = []
    n_single = 21
    n_total = 26        
    n_grouped = n_total-n_single #already did 0
    group_size = (card-n_single)/n_grouped # +1 for 0    
    
    #print n_grouped
    #print group_size

    for i in range(n_single):
        grouped_dists.append(dists[i])
        for j in range(max):
            output+=str(dists[i][j]) + " "
        output+="\n"
    
    for i in range(n_grouped):
        print("making group",i)
        sum = sum_g(dists,i,group_size,n_single)
        grouped_dists.append(sum)
        for j in range(max):
            output+=str(sum[j]) + " "
        output+="\n"
    
    #output+="%e "%(discrete_g_fit[j])
    #output+="%f "%(discrete_g_fit[j])
    open(fn_out,'w').write(output)

def allAboutTheSame( A, ep=1e-10 ):
    return np.all( (A>=(A.min()-ep)) & (A<=(A.min()+ep)) )

class LinearRegress:
    def __init__(self,x,y,w):
        self.x=x[np.argsort(x)]
        self.y=y[np.argsort(x)]
        self.w=w[np.argsort(x)]
        
        self.fit_params=None
        
    # xvals_use=array of x-values to use in regression; None means use all
    def fit(self,xvals_use=None,thruZero=False):
        if xvals_use==None:
            liuseX=np.arange(self.x.shape[0])
        else:
            liuseX=np.array([i for i in range(self.x.shape[0]) if abs(self.x[i]-xvals_use).min()<1e-6 ],'int32')
               
            if liuseX.shape[0]==0:
                sys.stderr.write('LinearRegress: warning, none of the specified domain values to use in fit was found: %s in %s\n'%(repr(xvals_use),repr(self.x)))
                liuseX=np.arange(self.x.shape[0])      
     
        if allAboutTheSame( self.x[liuseX] ):
            # implies thruzero.
            # one distinct x value, possibly multiple points
            
            slope0_i = self.y[liuseX]/self.x[liuseX]
            
            wcur = self.w[liuseX]
            if wcur.sum()<=1e-30: 
                wcur[:]=1.
            wcur=wcur/wcur.sum()
                                    
            self.fit_params = [0., (np.dot(slope0_i,wcur)/wcur.sum()) ]
        else:
            pinit = [0., 1.0]
            
            if thruZero:
                fitfunc = lambda p, x: p[1] * x
            else:
                fitfunc = lambda p, x: p[0] + p[1] * x

            wcur = self.w[liuseX]
            if wcur.sum()<=1e-30: 
                wcur[:]=1.
            wcur=wcur/wcur.sum()
                
            self.fit_params = do_lsqfit( self.x[liuseX], self.y[liuseX], wcur, pinit, fitfunc=fitfunc)

    def regress(self,x):
        return self.fit_params[0] + x * self.fit_params[1]
    
    def regrinv(self,y):
        return (y - self.fit_params[0]) / self.fit_params[1]
    
# Do two-piece linear regression
# NB THIS IS NOT SIMULTANEOUS PIECEWISE LINEAR REGRESSION
#  we don't want the higher copy points to influence the early fit.   So fit the first segment first, then the second, 
#  under the constraint that it meet the first
# i.e., fit [0,...2(configurable)] 
#  then fit another segment [3,...] requiring they meet at the same point (but not necessarily imposing zero y-intercept on the second
#  segment were it to be extended back to 0)

# update 10/23/2009 - making the point they meet at be (x_{r+1} - ep), where x_{r+1} is the first point of the right segment
#   rather than how this was previously ( x_{r} )

class TwoPieceNonSimulLinearRegress:
    def __init__(self,x,y,w):
        self.x=x[np.argsort(x)]
        self.y=y[np.argsort(x)]
        self.w=w[np.argsort(x)]
    
        self.xrngDomainLeftInclusive=None  # eg [0.,10.] means [0,10);[10,inf)        
        self.fit_params=None
    
    # xvals_use_0=array of x-values to use in regression first segment; None means use all
    # xvals_use_1=array of x-values to use in regression second segment; None means use all THAT REMAIN AFTER FIRST SEGMENT
    def fit(self,xvals_use_seg0=None,xvals_use_seg1=None,thruZero=False):
        if xvals_use_seg0==None:
            liuseX_0=np.arange(self.x.shape[0])
            liuseX_1=None
            self.xrngDomainLeftInclusive=[0.,1e40]  # second segment starts so high it will never be used
        else:
            liuseX_0=np.array([i for i in range(self.x.shape[0]) if abs(self.x[i]-xvals_use_seg0).min()<1e-6 ],'int32')
            
            if liuseX_0.shape[0]==0:
                sys.stderr.write('LinearRegress: warning, none of the specified domain values to use in fit was found: %s in %s\n'%(repr(xvals_use_seg0),repr(self.x)))
                liuseX_0=np.arange(self.x.shape[0])
                liuseX_1=None
                self.xrngDomainLeftInclusive=[0.,1e40]  # second segment starts so high it will never be used
            else:
                if xvals_use_seg1 == None:
                    # no x points specified for segment 1 -- so use everything larger than the biggest seg0 point
                    liuseX_1 = np.array([ i for i in range(self.x.shape[0]) 
                                         if self.x[i] - self.x[liuseX_0].max() > 1e-10 ], 'int32')
                    
                    if liuseX_1.shape[0]==0:
                        liuseX_1=None
                        self.xrngDomainLeftInclusive=[0.,1e40] # second segment starts so high it will never be used
                    else:
                        self.xrngDomainLeftInclusive=[0., self.x[liuseX_1].min() - 1e-10 ]
                else:
                    liuseX_1=np.array([i for i in range(self.x.shape[0]) if abs(self.x[i]-xvals_use_seg1).min()<1e-6 ],'int32')
                    self.xrngDomainLeftInclusive=[0., self.x[liuseX_1].min() - 1e-10 ]
                
            self.fit_params=[[],[]]
                    
            if allAboutTheSame( self.x[liuseX_0] ):
                # implies thruzero.
                # one distinct x value, possibly multiple points
                
                slope1_i = self.y[liuseX_0]/self.x[liuseX_0]
                
                wcur = self.w[liuseX_0]
                if wcur.sum()<=1e-30: 
                    wcur[:]=1.
                wcur=wcur/wcur.sum()
                                        
                self.fit_params[0] = [0., (np.dot(slope1_i,wcur)/wcur.sum()) ]                
            else:
                pinit=[0., 1.]
                if thruZero: fitfunc=lambda p,x:p[1]*x
                else: fitfunc=lambda p,x:p[0]+p[1]*x
                
                wcur=self.w[liuseX_0]
                if wcur.sum()<=1e-30: 
                    wcur[:]=1.
                wcur=wcur/wcur.sum()
                
                self.fit_params[0]=do_lsqfit( self.x[liuseX_0], self.y[liuseX_0], wcur, pinit, fitfunc=fitfunc)

            if liuseX_1!=None:
                if allAboutTheSame( self.x[liuseX_1] ):
                    
                    # if there's only one point in the right segment, we'll just keep the yintercept and slope from the left seg
                    self.fit_params[1]=self.fit_params[0][:]
                    
                    # Old, pre 10/23/2009
                    # second line segment is defined by 
                    # (x_r,yhat_r); (x_t, y_t)
                    # 
#                    x_r = self.x[liuseX_0].max()
#                    yhat_r = self.fit_params[0][0]+self.fit_params[0][1]*x_r
#                    
#                    yintcSecond = lambda slopeSecond:yhat_r - slopeSecond*x_r
#                    
#                    wcur = self.w[ liuseX_1 ]
#                    if wcur.sum()<=1e-30:
#                        wcur[:]=1.
#                    wcur=wcur/wcur.sum()
#                    
#                    yhat_t = np.dot( self.y[ liuseX_1 ], wcur ) / wcur.sum()
#                    
#                    self.fit_params[1] = [ yintcSecond(yhat_t), (yhat_t-yhat_r) / (self.x[ liuseX_1 ].mean() - x_r) ]
                else:                                  
                    wcur=self.w[liuseX_1]
                    if wcur.sum()<=1e-30: 
                        wcur[:]=1.
                    wcur=wcur/wcur.sum()
                    
                    x_r = self.x[liuseX_1].min() - 1e-6                    
                    
                    # Old, pre 10/23/2009
                    # second line segment goes through
                    # (x_r,yhat_r), where x_r is the maximal x point in the first line segment
#                    x_r = self.x[liuseX_0].max()

                    yhat_r = self.fit_params[0][0]+self.fit_params[0][1]*x_r
                    yintcSecond = lambda slopeSecond:yhat_r - slopeSecond*x_r 
                    fitfuncSecond = lambda slopeSecond,x:( yintcSecond(slopeSecond) + slopeSecond*x )
                    
                    # initialize to slope of first regrn
                    pinit=[ self.fit_params[0][1] ]
                                        
                    sl1 = do_lsqfit( self.x[liuseX_1], self.y[liuseX_1], wcur, pinit, fitfunc=fitfuncSecond)
                    self.fit_params[1] = [ yintcSecond(sl1), sl1 ]
            else:
                self.fit_params[1]=[1e50,0]  # this will force the yrngDomain to be very large for second segment.
             
        print(self.fit_params)
                      
        self.yrngDomainLeftInclusive=[ self.fit_params[0][0]+self.fit_params[0][1]*self.xrngDomainLeftInclusive[0],
                                       self.fit_params[1][0]+self.fit_params[1][1]*self.xrngDomainLeftInclusive[1], ]
                                    
    def regress(self,x):
        yhat_seg0 = self.fit_params[0][0] + x * self.fit_params[0][1]
        yhat_seg1 = self.fit_params[1][0] + x * self.fit_params[1][1]
        
        indic_is_seg0 = np.array( x<self.xrngDomainLeftInclusive[1] ).astype('f')
        indic_is_seg1 = 1.0 - indic_is_seg0 
            
        return (indic_is_seg0*yhat_seg0) + (indic_is_seg1*yhat_seg1) 
    
    def regrinv(self,y):
        xhat_seg0 = (y - self.fit_params[0][0])/self.fit_params[0][1]
        xhat_seg1 = (y - self.fit_params[1][0])/self.fit_params[1][1]
        
        indic_is_seg0 = np.array( y<self.yrngDomainLeftInclusive[1] ).astype('f')
        indic_is_seg1 = 1.0 - indic_is_seg0 
            
        return (indic_is_seg0*xhat_seg0) + (indic_is_seg1*xhat_seg1) 
        
        
class ml_get_cp2_trunc_gaussian:
    def __init__(self,fnfit_params,lUseCp=None,thruZero=False,max_cp=100):
        
        mtxParams = np.array([ [float(x) for x in line.replace('\n','').strip().split()[1:] ] 
                              for line in open(fnfit_params,'r') ],'float64')
        
        cps=mtxParams[:,0]
        mus=mtxParams[:,1]
        if mtxParams.shape[1]>2:
            sigs=mtxParams[:,2]
        else:
            sigs=mus[:]
        if mtxParams.shape[1]>3:
            nbps=mtxParams[:,3]
        else:
            nbps=(0.*cps)+1.
            
        self.mtxParams = mtxParams
            
        self.regressMu = LinearRegress( cps, mus, nbps )
        regressSig = LinearRegress( cps, sigs, nbps )

        self.regressMu.fit(lUseCp, thruZero) 
        regressSig.fit(lUseCp, thruZero)

        self.cns = np.arange(max_cp).astype('float64')
        self.cns[0] = .25 # picked .25 somewhat arbitrarily to reflect the fact that deleted regions
                           # can get hits from repetitive elts
        self.mus = self.regressMu.regress(self.cns)
        self.sigs = regressSig.regress(self.cns)
        
        self.trunc_rescale_ll = np.array( np.log( [ 1.0 - norm.cdf(0.,self.mus[icp],self.sigs[icp]) for icp in range(self.cns.shape[0]) ] ) )
        self.lImpropNorms=[ norm(self.mus[icp],self.sigs[icp]) for icp in range(self.cns.shape[0]) ]
        
        for icp in range(self.cns.shape[0]):
            print('from params: %s: %s -- %d  %.2f  %.2f'%(self.__class__.__name__, fnfit_params, int(self.cns[icp]), self.mus[icp],self.sigs[icp]))
                        
    def get_mle_cp(self,_depth):
        
        depth = _depth.astype('float64')
                
        depthPosReal = depth[ np.where( depth>=0 & ~(np.isnan(depth)|np.isinf(depth)|np.isneginf(depth)) )  ]
        if depthPosReal.shape[0]==0:
            sys.stderr.write('WARNING: no nonneg, real depths in %s.get_mle_cp()\n'%self.__class__ )
            return (0., -np.inf)
        
        #LLimprop = np.array([ (np.log((norm.pdf(depth,self.mus[cp],self.sigs[cp]))+1e-100).sum()) for cp in self.cns ],'float64')
        LL = np.array( [ (np.log( self.lImpropNorms[icp].pdf(depth) + 1e-100 ) - self.trunc_rescale_ll[icp]).sum() for icp in range(self.cns.shape[0]) ], 'float64'  )
                
        LLreal = LL[ np.where( ~(np.isnan(LL)|np.isinf(LL)|np.isneginf(LL)) )  ]
        if LLreal.shape[0] == 0:
            sys.stderr.write('WARNING: no nonzero, real logliks in %s.get_mle_cp()\n'%self.__class__ )
            return (0., -np.inf)
        
        cns = np.where( ~(np.isnan(LL)|np.isinf(LL)|np.isneginf(LL)) )[0]
        return (cns[np.where(LLreal==np.amax(LLreal))[0][0]],np.amax(LLreal))
    
    def get_density_across_cns(self, icp, cns):
        dens = self.lImpropNorms[ icp ].pdf( cns ) / np.exp(self.trunc_rescale_ll[icp])
        return dens 
        
class ml_get_cp2_trunc_gaussian_twoPieceRegress:
    def __init__(self,fnfit_params,lUseCp=[None,None],thruZero=False,max_cp=100):
       
        mtxParams = np.array([ [float(x) for x in line.replace('\n','').strip().split()[1:] ] 
                              for line in open(fnfit_params,'r') ],'float64')
        
        cps=mtxParams[:,0]
        mus=mtxParams[:,1]
        if mtxParams.shape[1]>2:
            sigs=mtxParams[:,2]
        else:
            sigs=mus[:]
        if mtxParams.shape[1]>3:
            nbps=mtxParams[:,3]
        else:
            nbps=(0.*cps)+1.
            
        self.mtxParams = mtxParams
            
        self.regressMu = TwoPieceNonSimulLinearRegress( cps, mus, nbps )
        regressSig = TwoPieceNonSimulLinearRegress( cps, sigs, nbps )
        
        self.regressMu.fit(lUseCp[0], lUseCp[1], thruZero) 
        regressSig.fit(lUseCp[0], lUseCp[1], thruZero)

        self.cns = np.arange(max_cp).astype('float64')
        self.cns[0] = .25 # picked .25 somewhat arbitrarily to reflect the fact that deleted regions
                           # can get hits from repetitive elts
        self.mus = self.regressMu.regress(self.cns)
        self.sigs = regressSig.regress(self.cns)
               
        self.trunc_rescale_ll = np.array( np.log( [ 1.0 - norm.cdf(0.,self.mus[icp],self.sigs[icp]) for icp in range(self.cns.shape[0]) ] ) )
        self.lImpropNorms=[ norm(self.mus[icp],self.sigs[icp]) for icp in range(self.cns.shape[0]) ]
        
        for icp in range(self.cns.shape[0]):
            print('from params: %s: %s -- %d  %.2f  %.2f'%(self.__class__.__name__, fnfit_params, int(self.cns[icp]), self.mus[icp],self.sigs[icp]))
                        
    def get_mle_cp(self,_depth):
        
        depth = _depth.astype('float64')
        
        depthPosReal = depth[ np.where( depth>=0 & ~(np.isnan(depth)|np.isinf(depth)|np.isneginf(depth)) )  ]
        if depthPosReal.shape[0]==0:
            sys.stderr.write('WARNING: no nonneg, real depths in %s.get_mle_cp()\n'%self.__class__ )
            return (0., -np.inf)
        
        #LLimprop = np.array([ (np.log((norm.pdf(depth,self.mus[cp],self.sigs[cp]))+1e-100).sum()) for cp in self.cns ],'float64')
        LL = np.array( [ (np.log( self.lImpropNorms[icp].pdf(depth) + 1e-100 ) - self.trunc_rescale_ll[icp]).sum() for icp in range(self.cns.shape[0]) ], 'float64'  )
        
        LLreal = LL[ np.where( ~(np.isnan(LL)|np.isinf(LL)|np.isneginf(LL)) )  ]
        if LLreal.shape[0] == 0:
            sys.stderr.write('WARNING: no nonzero, real logliks in %s.get_mle_cp()\n'%self.__class__ )
            return (0., -np.inf)
        
        cns = np.where( ~(np.isnan(LL)|np.isinf(LL)|np.isneginf(LL)) )[0]
        return (cns[np.where(LLreal==np.amax(LLreal))[0][0]],np.amax(LLreal))
    
    def get_density_across_cns(self, icp, cns):
        dens = self.lImpropNorms[ icp ].pdf( cns ) / np.exp(self.trunc_rescale_ll[icp])
        return dens 
    
class ml_get_cp2_cont_poisson_twoPieceRegress:
    # lUseCp = the list of copy number BACs .  None means use all.
    # dim0: depths
    # dim1: starts
    # dim2: can-style starts.
    def __init__(self,fnfit_params,lUseCp=[None,None],thruZero=False,max_cp=100):
        
        mtxParams = np.array([ [float(x) for x in line.replace('\n','').strip().split()[1:] ] 
                              for line in open(fnfit_params,'r') ],'float64')
        
        cps=mtxParams[:,0]
        mus=mtxParams[:,1]
        if mtxParams.shape[1]>2:
            nbps=mtxParams[:,2]
        else:
            nbps=(0.*cps)+1.
            
        self.mtxParams = mtxParams
            
        self.regressMu = TwoPieceNonSimulLinearRegress( cps, mus, nbps )
        
        self.regressMu.fit(lUseCp[0], lUseCp[1], thruZero) 

        self.cns = np.arange(max_cp).astype('float64')
        self.cns[0] = .25 # for continuous poisson only, having a zero class implies a log(0) 
                           # picked .25 somewhat arbitrarily to reflect the fact that deleted regions
                           # can get hits from repetitive elts
        self.mus = self.regressMu.regress(self.cns)
        
        for icp in range(self.cns.shape[0]):
            print('from params: %s: %s -- %d  %.2f'%(self.__class__.__name__, fnfit_params, int(self.cns[icp]), self.mus[icp] ))

    def get_mle_cp(self,_depth):

        depth = _depth.astype('float64')

        depthPosReal = depth[ np.where( depth>=0 & ~(np.isnan(depth)|np.isinf(depth)|np.isneginf(depth)) )  ]
        if depthPosReal.shape[0]==0:
            sys.stderr.write('WARNING: no nonneg, real depths in %s.get_mle_cp()\n'%self.__class__ )
            return (0., -np.inf)
        
        N=float(depth.shape[0])
        LL = np.array( [ (depth*np.log(self.mus[icp])).sum() - 
                        np.log( scipy.special.gamma(depth+1.0).clip(0.,1e300) ).sum() -             # gamma can overflow --> inf 
                        N*self.mus[icp] + 1e-100 for icp in range(self.cns.shape[0]) ], 'float64') 
                                                
        LLreal = LL[ np.where( ~(np.isnan(LL)|np.isinf(LL)|np.isneginf(LL)) )  ]
        if LLreal.shape[0] == 0:
            sys.stderr.write('WARNING: no nonzero, real logliks in %s.get_mle_cp()\n'%self.__class__ )
            return (0., -np.inf)
        
        cns = np.where( ~(np.isnan(LL)|np.isinf(LL)|np.isneginf(LL)) )[0]
        return (cns[np.where(LLreal==np.amax(LLreal))[0][0]],np.amax(LLreal))
    
    
    def get_density_across_cns(self, icp, cns):
        dens = ((self.mus[icp])**cns)*np.exp(-self.mus[icp]) / scipy.special.gamma( cns + 1. )
        return dens 
    
class ml_get_cp2_cont_poisson:
    # lUseCp = the list of copy number BACs .  None means use all.
    # dim0: depths
    # dim1: starts
    # dim2: can-style starts.
    def __init__(self,fnfit_params,lUseCp=None,thruZero=False,max_cp=100):
        
        mtxParams = np.array([ [float(x) for x in line.replace('\n','').strip().split()[1:]] 
                              for line in open(fnfit_params,'r') ],'float64')
       
        cps=mtxParams[:,0]
        mus=mtxParams[:,1]
        if mtxParams.shape[1]>2:
            nbps=mtxParams[:,2]
        else:
            nbps=(0.*cps)+1.
            
        self.mtxParams = mtxParams

        self.regressMu = LinearRegress( cps, mus, nbps )
        
        self.regressMu.fit(lUseCp, thruZero) 

        self.cns = np.arange(max_cp).astype('float64')
        self.cns[0] = .25 # for continuous poisson only, having a zero class implies a log(0) 
                           # picked .25 somewhat arbitrarily to reflect the fact that deleted regions
                           # can get hits from repetitive elts
        self.mus = self.regressMu.regress(self.cns)
        
        for icp in range(self.cns.shape[0]):
            print('from params: %s: %s -- %d  %.2f'%(self.__class__.__name__, fnfit_params, int(self.cns[icp]), self.mus[icp] ))
            
    def get_mle_cp(self,_depth):

        depth = _depth.astype('float64')

        depthPosReal = depth[ np.where( depth>=0 & ~(np.isnan(depth)|np.isinf(depth)|np.isneginf(depth)) )  ]
        if depthPosReal.shape[0]==0:
            sys.stderr.write('WARNING: no nonneg, real depths in %s.get_mle_cp()\n'%self.__class__ )
            return (0., -np.inf)
        
        N=float(depth.shape[0])
        LL = np.array( [ (depth*np.log(self.mus[icp])).sum() - 
                        np.log( scipy.special.gamma(depth+1.0).clip(0.,1e300) ).sum() -             # gamma can overflow --> inf 
                        N*self.mus[icp] + 1e-100 for icp in range(self.cns.shape[0]) ], 'float64') 
                        
        LLreal = LL[ np.where( ~(np.isnan(LL)|np.isinf(LL)|np.isneginf(LL)) )  ]
        if LLreal.shape[0] == 0:
            sys.stderr.write('WARNING: no nonzero, real logliks in %s.get_mle_cp()\n'%self.__class__ )
            #from IPython.Shell import IPShellEmbed; ipshell = IPShellEmbed([]); ipshell()
            return (0., -np.inf)
        
        cns = np.where( ~(np.isnan(LL)|np.isinf(LL)|np.isneginf(LL)) )[0]
        return (cns[np.where(LLreal==np.amax(LLreal))[0][0]],np.amax(LLreal))

    def get_density_across_cns(self, icp, cns):
        dens = ((self.mus[icp])**cns)*np.exp(-self.mus[icp]) / scipy.special.gamma( cns + 1. )
        return dens 


def analyze_locations(ml,in_locations,fn_genome_analysis_dir,fn_contigs,fn_mask):
    
    fn_combined_adjusted_wssd = "%s/combined_corrected_wssd/wssd.combined_corrected"%(fn_genome_analysis_dir)

    combined_adjusted_wssd = WssdFile(fn_contigs,
                                      fnWssd=fn_combined_adjusted_wssd,
                                      overwrite=False,
                                      openMode='r')
    
    mask_wssd = DenseTrackSet(fn_contigs,
                              fnWssd=fn_mask,
                              overwrite=False,
                              openMode='r')

    for location in in_locations:
        if(location[0] == "#"):continue
        (name,chr,start,end) = location.split()
        
        corrected_depth = np.nan_to_num(combined_adjusted_wssd.depth["wssd.combined_corrected"][chr][start:end,:,0]).astype(np.float64).sum(1)
        masked_region = mask_wssd["mask"][chr][start:end,:].sum(1)>0
        depth = corrected_depth[np.where(masked_region==0)]
        cp,ll = ml.get_mle_cp(depth)
        print(name,chr,start,end,cp[0],ll)
    
    combined_adjusted_wssd.tbl.close()
    mask_wssd.tbl.close()

def in_completed(fn_completed,indiv):

    for completed_indiv in open(fn_completed,'r').readlines():
        if( completed_indiv.split("/")[9] == indiv): return True

    print("%s does not exist in %s, skipping"%(indiv,fn_completed))
    return False

if __name__=='__main__':
    
    opts = OptionParser()
    opts.add_option('','--in_genomes',dest='fn_in_genomes')
    opts.add_option('','--in_regions',dest='fn_in_locations')
    opts.add_option('','--completed_wssds',dest='fn_completed_wssds')
    opts.add_option('','--contigs',dest='fn_contigs')
    opts.add_option('','--mask_file',dest='fn_mask')

    (o, args) = opts.parse_args()
    
    in_genomes = open(o.fn_in_genomes,"r").readlines()
    in_locations = open(o.fn_in_locations,"r").readlines()
    primary_analysis_dir = "./primary_analysis"    

    for in_genome in in_genomes:
        (genome,fn_wssd_dir,fn_bac_dir) = in_genome.split()
        genome_dir_name = genome.split(".")[0] 
         
        if(not(in_completed(o.fn_completed_wssds,genome))):
            continue

        wssd_lib_dir = "%s/%s"%(fn_wssd_dir,genome)
        bac_analysis_lib_dir = "%s/%s"%(fn_bac_dir,genome_dir_name)
        print(genome, genome_dir_name)
        filtered_libs_dir = "%s/_filtered_libs_analysis"%(bac_analysis_lib_dir)
        if(os.path.exists(filtered_libs_dir)):
            ml = ml_get_cp(bac_analysis_lib_dir)
            analyze_locations(ml,in_locations,"%s/%s"%(primary_analysis_dir,genome),o.fn_contigs,o.fn_mask)    
            

    sys.exit(1)

