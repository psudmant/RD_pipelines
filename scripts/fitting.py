import sys,os
import numpy as np
import scipy as scp
import kg_file_handling as kgf
from optparse import OptionParser
from scipy.stats import norm
from wssd_common import *
from wssd_pw_common import *


class distribution:

	def get_trunc_gauss(self,mu_sig,ret_range):
		mu,sigma = mu_sig
		ndist = norm(mu,sigma)
		x = np.arange(ret_range[0],ret_range[1],ret_range[2])
		correction_factor = 1-norm.cdf(0)	
		print correction_factor
		ret = ndist.pdf(x)/correction_factor
		return ret


class fit:
	def __init__(self,data):
		self.input_data = data

	def fit_by_fmin(self):
		
		x = self.input_data
		N = self.input_data.shape[0]
		trunc_gaussian = lambda p: ( check_p(p)*(-1/(2*p[1]**2))*(np.square(x-p[0]).sum()) - N*np.log(p[1]) - N*np.log(1-scp.stats.norm(p[0],p[1]).cdf(0)) )**2

		check_p = lambda p: (p[0]<=0 or p[1]<=0) and 1e9 or 1
		pinit = [0,1]
		ps = scp.optimize.fmin(trunc_gaussian,pinit)
		ps = (ps[0],ps[1])
		return ps

	def fit_by_simple(self):
		#return(np.median(self.input_data),self.input_data.var())
		return(self.input_data.mean(),self.input_data.var())

	def fit_by_weibull_lsq(self):
		fitfunc = lambda p, x: (p[0]/p[1]*(x/p[1])**(p[0]-1)*np.exp(-(x/p[1])**p[0]))		
		errfunc = lambda p, x, y: (y - fitfunc(p, x))
		
		hist = np.histogram(self.input_data,bins=3000,range=[0,3000])[0].astype(np.float64)

		pdf = hist/hist.sum()
		x=np.arange(3000)
		y=pdf
		pinit = [60, 100]
		out = scp.optimize.leastsq(errfunc, pinit,args=(x,y), full_output=1)
		xs = out[0]
		fit_ys = fitfunc(xs,x)
		xs[1] = (xs[1])
		return xs	

	def fit_by_gaussian_lsq(self,ret_range=None):
		
		fitfunc = lambda p, x: (1/(p[1]*np.sqrt(2*math.pi))*np.exp(-((x-p[0])**2)/(2*p[1]**2)))

		errfunc = lambda p, x, y: (y - fitfunc(p, x))
		hist = np.histogram(self.input_data,bins=3000,range=[0,3000])[0].astype(np.float64)

		pdf = hist/hist.sum()
		x=np.arange(3000)
		y=pdf
		pinit = [60, 100]
		out = scp.optimize.leastsq(errfunc, pinit,args=(x,y), full_output=1)
		xs = out[0]
		fit_ys = fitfunc(xs,x)
		xs[1] = (xs[1])

		if(ret_range!=None):
			cont_poisson_pmf = lambda p: ( (p**x)*np.exp(-p)/scp.special.gamma(x+1) ) 
			x = ret_range
			fitline = cont_poisson_pmf(p[0])
			fitline[np.where(np.isnan(fitline))] = 0
			return xs,fitline
		else:
			return xs

	def fit_by_scipy_stat(self):
		return norm.fit(self.input_data)

	def fit_gamma(self):
		p = scp.stats.gamma.fit(np.r_[self.input_data,1e-10])	
		print p
		return  (p[0],p[2])

	def fit_by_cont_poisson(self,ret_range=None):
			
		x = self.input_data
		#cont_poisson = lambda p:  ( ( x*np.log(p) - p - np.log(scp.special.gamma(x+1))).sum())
		#cont_poisson = lambda p:  ( ( x*np.log(p) - p - np.log(scp.special.gamma(x+1)) ).sum())**2
		#pinit = [1]
		#p = scp.optimize.fmin(cont_poisson,pinit)
		p = [x.mean()]
		hist = np.histogram(self.input_data,bins=3000,range=[0,3000])[0].astype(np.float64)
		print 'param--->', p
		print "max", np.where(hist==np.max(hist))
		if(ret_range!=None):
			cont_poisson_pmf = lambda p: ( (p**x)*np.exp(-p)/scp.special.gamma(x+1) ) 
			x = ret_range
			fitline = cont_poisson_pmf(p[0]) *ret_range[2]
			fitline[np.where(np.isnan(fitline))] = 0
			return p,fitline 
		
		return p
		
	def fit_by_em(self):
		mod_input_data = np.reshape(self.input_data,[self.input_data.shape[0],1])
		lgm = em.GM(1,1,'diag')	
		gmm = em.GMM(lgm,'kmean')
		emaximize = em.EM()
		like = emaximize.train(mod_input_data,gmm,maxiter=30,thresh=1e-8)	
		#print "likelihood", like
		return(gmm.gm.mu,gmm.gm.va)

