from optparse import OptionParser
import json
import glob

from sys import stderr
import os
from wnd_cp_data import wnd_cp_indiv
import time
import numpy as np
import pandas as pd

class gglob:
    """
    a class for loading and saving giant matrices of DTS 
    windowed genome data. A gglob_dir has 2 components:
    1. gglob.index
    2. chr{*}.npz for all chrs 

    TODO:
    add 2 classmethods
    1. basic, load_from_gglob
    2. load_from_DTS (basically just runs the code in the init)
    """
    
    @classmethod
    def get_indivs(cls, DTS_dir, sunk_DTS_dir, DTS_prefix):
        indivs = []
        for f in glob.glob("%s/%s*"%(DTS_dir,DTS_prefix)):
            indiv = f.split("/")[-1].replace(DTS_prefix,"") 
            if os.path.exists("%s/%s%s"%(sunk_DTS_dir, DTS_prefix, indiv )):
                indivs.append(indiv)
            else:
                print >>stderr, "skipping: %s - no associated sunk DTS"%indiv
        return indivs

    @classmethod
    def init_from_DTS(cls, **kwargs):
        """
        requires the below inputs

        gglob.init_from_DTS(DTS_dir = DTS_dir,
                            DTS_prefix = DTS_prefix,
                            sunk_DTS_dir = sunk_DTS_dir,
                            sunk_DTS_prefix = sunk_DTS_prefix,
                            wnd_size = wnd_size,
                            indivs = indivs,
                            contig = contig,
                            fn_contigs = fn_contigs,
                            fn_sunk_contigs = fn_sunk_contigs)

        """
        
        DTS_dir = kwargs["DTS_dir"]
        DTS_prefix = kwargs["DTS_prefix"]
        
        sunk_DTS_dir = kwargs["sunk_DTS_dir"]
        sunk_DTS_prefix = kwargs["sunk_DTS_prefix"]
        
        wnd_size = kwargs['wnd_size']
        wnd_slide = kwargs['wnd_slide']

        indivs = kwargs['indivs']
        contig = kwargs['contig']

        fn_contigs = kwargs['fn_contigs']
        fn_sunk_contigs = kwargs['fn_sunk_contigs']

        DTS_pre="%s/%s"%(DTS_dir, DTS_prefix) 
        sunk_DTS_pre="%s/%s"%(sunk_DTS_dir, DTS_prefix) 
        
        n_indivs = len(indivs)
        
        t = time.time()
        rand_wnd_cp = wnd_cp_indiv("%s%s"%(DTS_pre, indivs[0]), fn_contigs, wnd_size)
        wnd_starts, wnd_ends = rand_wnd_cp.get_wnds_by_chr(contig)
        cp_matrix = np.zeros((n_indivs, wnd_starts.shape[0]))

        rand_sunk_wnd_cp = wnd_cp_indiv("%s%s"%(sunk_DTS_pre, indivs[0]), fn_sunk_contigs, wnd_size)
        sunk_wnd_starts, sunk_wnd_ends = rand_sunk_wnd_cp.get_wnds_by_chr(contig)
        sunk_cp_matrix = np.zeros((n_indivs, sunk_wnd_starts.shape[0]))
        
        correct = not (contig in ["chrY", "chrX"])

        for i, indiv in enumerate(indivs):
            print >> stderr, indiv
            wnd_cp = wnd_cp_indiv("%s%s"%(DTS_pre, indiv),
                                  fn_contigs,
                                  wnd_size)
            
            cp_matrix[i,:] = wnd_cp.get_cps_by_chr(contig, correct=correct) 

            sunk_wnd_cp = wnd_cp_indiv("%s%s"%(sunk_DTS_pre, indiv), 
                                      fn_sunk_contigs,
                                      wnd_size)
            
            sunk_cp_matrix[i,:] = sunk_wnd_cp.get_cps_by_chr(contig, correct=correct) 
        
        return cls(indivs = indivs,
                   wnd_size = wnd_size,
                   wnd_slide = wnd_slide,
                   contig = contig,
                   wnd_starts = wnd_starts,
                   wnd_ends = wnd_ends, 
                   cp_matrix = cp_matrix, 
                   sunk_wnd_starts = sunk_wnd_starts,
                   sunk_wnd_ends = sunk_wnd_ends,
                   sunk_cp_matrix = sunk_cp_matrix)
        
    @classmethod
    def init_from_gglob_dir(cls, gglob_dir, contig, indiv_subset=None):
        
        #open up the index
        idx_data = None
        
        with open("%s/gglob.idx"%gglob_dir) as F:
            idx_data = json.load(F)
        
        indivs = idx_data['indivs'] 
        wnd_size = idx_data['wnd_size'] 
        wnd_slide = idx_data['wnd_slide'] 
        contig = contig
               
        keys = ["wnd_starts","wnd_ends","cp_matrix","sunk_wnd_starts","sunk_wnd_ends","sunk_cp_matrix"]
        fn_in = "%s/%s"%(gglob_dir,contig)

        mats_by_key = {} 
        for k in keys:
            stderr.write("loading %s..."%k)
            stderr.flush()
            t=time.time()
            if not os.path.exists("%s.%s.h5"%(fn_in,k)):
                raise Exception("""path "%s.%s.h5" does not exist!"""%(fn_in,k))

            df = pd.read_hdf("%s.%s.h5"%(fn_in,k),k)
            mats_by_key[k] = df.as_matrix()
            stderr.write("done (%fs)\n"%(time.time()-t))

        wnd_starts = mats_by_key['wnd_starts'][:,0]
        wnd_ends = mats_by_key['wnd_ends'][:,0]
        cp_matrix = mats_by_key['cp_matrix']
        
        sunk_wnd_starts = mats_by_key['sunk_wnd_starts'][:,0]
        sunk_wnd_ends = mats_by_key['sunk_wnd_ends'][:,0]
        sunk_cp_matrix = mats_by_key['sunk_cp_matrix']
        
        if indiv_subset!=None:
            new_indivs = []
            new_indivs_idxs = []
            i=0
            for indiv in indivs:
                if indiv in indiv_subset:
                    new_indivs.append(indiv)
                    new_indivs_idxs.append(i)
                i+=1    
            l = len(new_indivs)
            new_indivs_idxs = np.array(new_indivs_idxs)
            
            cp_matrix = cp_matrix[new_indivs_idxs,:]  
            
            sunk_cp_matrix = sunk_cp_matrix[new_indivs_idxs,:] 
            indivs = new_indivs
        
        #for i,indiv in enumerate(indivs):
        #    print indiv, np.mean(cp_matrix[i]), np.median(cp_matrix[i])

        return cls(indivs = indivs,
                   wnd_size = wnd_size,
                   wnd_slide = wnd_slide,
                   contig = contig,
                   wnd_starts = wnd_starts,
                   wnd_ends = wnd_ends, 
                   cp_matrix = cp_matrix, 
                   sunk_wnd_starts = sunk_wnd_starts,
                   sunk_wnd_ends = sunk_wnd_ends,
                   sunk_cp_matrix = sunk_cp_matrix)
    
    def __init__(self, **kwargs):
        """
        the from_gglob_dir or from_DTS classmethods shoudl be used
        for constructors
        """

        self.indivs = kwargs['indivs'] 
        self.wnd_size = kwargs['wnd_size'] 
        self.wnd_slide = kwargs['wnd_slide'] 
        self.contig = kwargs['contig']
        
        self.wnd_starts = kwargs['wnd_starts']
        self.wnd_ends = kwargs['wnd_ends']
        self.cp_matrix = kwargs['cp_matrix']
        
        self.sunk_wnd_starts = kwargs['sunk_wnd_starts']
        self.sunk_wnd_ends = kwargs['sunk_wnd_ends']
        self.sunk_cp_matrix = kwargs['sunk_cp_matrix']
        
        assert self.sunk_cp_matrix.shape[0] == len(self.indivs) 
        
    def save(self, gglob_dir):
        """
        save contig to a gglob_dir
        """
        fn_out =  "%s/%s"%(gglob_dir, self.contig)

        mats = {"wnd_starts":self.wnd_starts,
                "wnd_ends":self.wnd_ends,
                "cp_matrix":self.cp_matrix,
                "sunk_wnd_starts":self.sunk_wnd_starts,
                "sunk_wnd_ends":self.sunk_wnd_ends,
                "sunk_cp_matrix":self.sunk_cp_matrix
               }
         
        for k, mat in mats.iteritems():
            print >>stderr, "writing out %s..."%k
            t=time.time()
            df = pd.DataFrame(mat)
            df.to_hdf("%s.%s.h5"%(fn_out,k),k,complevel=1,complib='zlib')
            print >>stderr, "done (%fs)"%(time.time()-t)
        
        print >>stderr, "done (%f)"%(time.time()-t)



if __name__=="__main__":
        
    opts = OptionParser()

    usage_txt = "python gglob.py [options] [--init || --setup_chr {chr}]"
    opts = OptionParser(usage=usage_txt)

    opts.add_option('','--in_DTS_dir',dest='DTS_dir', default=None)
    opts.add_option('','--in_sunk_DTS_dir',dest='sunk_DTS_dir', default=None)

    opts.add_option('','--contigs',dest='fn_contigs', default=None)
    opts.add_option('','--sunk_contigs',dest='fn_sunk_contigs', default=None)
    
    opts.add_option('','--wnd_size',dest='wnd_size', type=int, default=None)
    opts.add_option('','--wnd_slide',dest='wnd_slide', type=int, default=None)
    
    opts.add_option('','--DTS_prefix',dest='DTS_prefix', default="500_bp_")
    opts.add_option('','--sunk_DTS_prefix',dest='sunk_DTS_prefix', default="500_bp_")
    
    opts.add_option('','--gglob_dir',dest='gglob_dir')

    opts.add_option('','--setup_chr',dest='setup_chr', default=None)
    opts.add_option('','--init',dest='init', action='store_true', default=False)


    (o, args) = opts.parse_args()
    """
    USAGE:

    first init (it makes an idx file w/ all indivs
    THEN: run the main program

    **NOTE, need to test the old 'build glob dir' again

    """

    if o.init == False and o.setup_chr == None:
        opts.print_help()

    elif o.init:
        indivs = []
        idx_data = {} 
        for f in glob.glob("%s/%s*"%(o.DTS_dir,o.DTS_prefix)):
            indiv = f.split("/")[-1].replace(o.DTS_prefix,"") 
            if os.path.exists("%s/%s%s"%(o.sunk_DTS_dir, o.DTS_prefix, indiv )):
                indivs.append(indiv)
            else:
                print >>stderr, "skiping: %s - no associated sunk DTS"
        
        idx_data = {"indivs":indivs, "wnd_size":o.wnd_size, "wnd_slide":o.wnd_slide}

        with open("%s/gglob.idx"%o.gglob_dir, 'w') as F:
            json.dump(idx_data, F)
    
    else:
        
        idx_data = None
        with open("%s/gglob.idx"%o.gglob_dir) as F:
            idx_data = json.load(F)
        
        DTS_pre="%s/%s"%(o.DTS_dir, o.DTS_prefix) 
        sunk_DTS_pre="%s/%s"%(o.sunk_DTS_dir, o.DTS_prefix) 
        
        wnd_size = idx_data['wnd_size']
        wnd_slide = idx_data['wnd_slide']
        indivs = idx_data['indivs']
        n_indivs = len(indivs)
        contig = o.setup_chr
        
        g = gglob.init_from_DTS(DTS_dir = o.DTS_dir,
                                DTS_prefix = o.DTS_prefix,
                                sunk_DTS_dir = o.sunk_DTS_dir,
                                sunk_DTS_prefix = o.sunk_DTS_prefix,
                                wnd_size = wnd_size,
                                wnd_slide = wnd_slide,
                                indivs = indivs,
                                contig = contig,
                                fn_contigs = o.fn_contigs,
                                fn_sunk_contigs = o.fn_sunk_contigs)
        
        g.save(o.gglob_dir)
        


