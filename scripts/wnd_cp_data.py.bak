import numpy as np
from wssd_pw_common import *
import pysam


class dCGH():
    def __init__(self,fn_dts_test,fn_dts_ref,fn_contigs,wnd_size):
        self.wnd_cp_ref=wnd_cp_indiv(fn_dts_ref,fn_contigs,wnd_size)
        self.wnd_cp_test=wnd_cp_indiv(fn_dts_test,fn_contigs,wnd_size)
        
    def get_cps_by_chr(self,chr,correct=True):
        return self.get_cp_ratio_by_chr(chr,correct)

    def get_cp_ratio_by_chr(self,chr,correct=True):
        ref_cps = self.wnd_cp_ref.get_cps_by_chr(chr,correct=correct)
        sample_cps = self.wnd_cp_test.get_cps_by_chr(chr,correct=correct)

        ratios = np.log(np.cast[np.float32](sample_cps/ref_cps))/np.log(2)
        
        while( not( np.sum(np.isnan(ratios))==0 and np.sum(np.isinf(ratios))==0) ):
            
            nans_locs = arrayToNonZeroContigRanges(np.isnan(ratios))
            infs_locs = arrayToNonZeroContigRanges(np.isinf(ratios))
            
            for i in xrange(nans_locs.shape[0]):
                start,end=nans_locs[i]
                l=(start==0) and end or start-1
                r=end>=(ratios.shape[0]-1) and start-1 or end
                ratios[start:end] = (ratios[l]+ratios[r])/2.0 
            
            for i in xrange(infs_locs.shape[0]):
                start,end=infs_locs[i]
                l=start==0 and end or start-1
                r=end>=(ratios.shape[0]-1) and start-1 or end
                ratios[start:end] = (ratios[l]+ratios[r])/2.0 
        
        assert np.sum(np.isnan(ratios))==0
        assert np.sum(np.isinf(ratios))==0
        
        return ratios    
    
    def get_cp_dup_loci(self, chr):       
        return self.wnd_cp_test.get_cp_dup_loci(chr), self.wnd_cp_ref.get_cp_dup_loci(chr)

    def get_wnds_by_chr(self,chr):
        return self.wnd_cp_ref.get_wnds_by_chr(chr)
    
    def get_overlapping_wnds(self,chr,tbx):
        return self.wnd_cp_ref.get_overlapping_wnds(chr,tbx)

def getNonZeroRanges( A, zeroDefVal=0 ):
    """
    get position of adjactent 1st in a vector 
    """
    dA=np.diff(A)
    N=A.shape[0]
    nzdA,=np.nonzero(dA)
    liRangesAll=np.c_[ np.r_[0,nzdA+1],
    np.r_[nzdA,N-1] +1]
    iiR = A.take(liRangesAll[:,0]) != zeroDefVal
    return liRangesAll[iiR,:]

class wnd_cp_indiv:
    def __init__(self,fn_dts,fn_contigs,wnd_size):
       """
       object for an individual
       """

       wnd = int(fn_dts.split("/")[-1].split("_bp")[0])
       assert int(wnd)==wnd_size
       self.wnd_size = wnd_size
       print fn_dts
       self.wnd_DTS = DenseTrackSet(fn_contigs,
                                   fn_dts,
                                   overwrite=False,
                                   openMode='r')
       
       self.contigs = self.wnd_DTS.mContigNameLen
       self.starts=self.wnd_DTS["starts"]
       self.ends=self.wnd_DTS["ends"]
       self.cps=self.wnd_DTS["copy"]
    
    def get_wnds_by_chr(self,chr):
        return self.starts[chr][:],self.ends[chr][:]
    

    def get_cp_dup_loci(self, chr):
        cps = self.get_cps_by_chr(chr)
        dup_regions = cps>3.5
        ret = getNonZeroRanges(dup_regions)
        return ret
        


    def get_overlapping_wnds(self,chr,tbx):
        wnd_starts, wnd_ends = self.get_wnds_by_chr(chr)
        bnds = np.array([ [int(l[1]),int(l[2])] 
                            for l in tbx.fetch(chr,parser=pysam.asTuple()) ])
         
        start_idxs = np.searchsorted(wnd_starts,bnds[:,0])
        end_idxs = np.searchsorted(wnd_starts,bnds[:,1])
        #print start_idxs
        #print end_idxs
        ret = np.c_[start_idxs,end_idxs]
        
        return ret
    
    def get_gapped_wnds(self,chr,tbx_gaps):
        
        gapped_wnds = []

        for t in tbx_gaps.fetch(chr,parser=pysam.asTuple()):
            _chr,start,end = t
            wnd_start = np.searchsorted(self.starts,start) 
            wnd_end = np.searchsorted(self.starts,end) 
            gapped_wnds.append(tuple([wnd_start,wnd_end]))
        
        return gapped_wnds 

    def get_cps_by_chr(self,chr,correct=False):
        cps=self.cps[chr][:]
        if correct:
            where_not_0 = np.where(cps!=0)
            adj=2-np.median(cps[where_not_0])
            print "correct by:%f median: %f mean: %f"%(adj,np.median(cps),np.mean(cps))
            cps[where_not_0] = cps[where_not_0]+adj
        return cps
              
