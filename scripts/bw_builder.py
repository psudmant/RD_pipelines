from optparse import OptionParser
import json
import glob
import tempfile

from sys import stderr
import sys
import os
import time
import numpy as np
import pandas as pd

from wnd_cp_data import wnd_cp_indiv
from gglob import gglob

class output:
    def __init__(self, contig_prefix, output_contigs):
        self.outdata = []
        self.contig_prefix=contig_prefix
        self.output_contigs=output_contigs
    
    def add(self, contig, start, end, cp):
        self.outdata.append([contig, start, end, rnd_cp])
        
    def output(self, fn_out, indiv, name):            
        color_hash =  { 0 :"229,229,229",
                        1 :"196,196,196",
                        2 :"0,0,0",
                        3 :"0,0,205",
                        4 :"65,105,225",
                        5 :"100,149,237",
                        6 :"180,238,180",
                        7 :"255,255,0",
                        8 :"255,165,0",
                        9 :"139,26,26",
                        10 :"255,0,0"}
        
        s_outdata = [self.outdata[0]]
        for i in xrange(1,len(self.outdata)):
            contig, start, end, cp = self.outdata[i]
            if contig == s_outdata[-1][0] and cp == s_outdata[-1][3]:
                s_outdata[-1][2] = end
            else:
                s_outdata.append(self.outdata[i])
        
        fn_tmp = tempfile.NamedTemporaryFile(mode='w', dir="/tmp").name
        print fn_tmp
        with open(fn_tmp, 'w') as F:
            for l in s_outdata:
                contig, start, end, cp = l
                rnd_cp = min(int(round(cp)), 10)
                print >>F, "\t".join(["%s%s"%(self.contig_prefix,contig),str(start),str(end),indiv,"0","+","0","0",color_hash[rnd_cp],cp]) 
            
        #hg19_contigs = "/net/eichler/vol7/home/psudmant/genomes/contigs/hg19_contigs.txt"
        contigs = self.output_contigs 
        cmd = "/net/eichler/vol7/home/psudmant/local_installations/ucscOLD/ucsc/bin/bedToBigBed %s %s %s"%(fn_tmp,contigs,fn_out)
        track_def = """track type=bigBed name="%s_%s" description="%s_%s" visibility=dense itemRgb="On" dataUrl=%s\n"""%(indiv, 
                                                                                                                         name, 
                                                                                                                         indiv,
                                                                                                                         name,
                                                                                                                         fn_out)
        print cmd 
        ret = os.system(cmd)
        if ret != 0:
            exit(ret)
        print fn_tmp
        os.unlink(fn_tmp)
        fn_out_td = "%s.trackdef"%(fn_out)
        with open(fn_out_td,'w') as F:
            F.write(track_def)


if __name__=="__main__":
        
    opts = OptionParser()
    
    opts.add_option('','--fn_DTS',dest='fn_DTS', default=None)
    opts.add_option('','--contigs',dest='fn_contigs', default=None)
    opts.add_option('','--wnd_size',dest='wnd_size', type=int, default=None)
    #opts.add_option('','--wnd_slide',dest='wnd_slide', type=int, default=None)
    opts.add_option('','--out_dir',dest='out_dir')
    opts.add_option('','--fn_out',dest='fn_out')
    opts.add_option('','--contig_prefix',dest='contig_prefix', default="")
    opts.add_option('','--DTS_prefix',dest='DTS_prefix', default="500_bp_")
    opts.add_option('','--output_contigs',dest='output_contigs', default="/net/eichler/vol7/home/psudmant/genomes/contigs/hg19_contigs.txt")

    
    (o, args) = opts.parse_args()
    #usage, init, then run
    
    indiv = o.fn_DTS.split("/")[-1].replace("500_bp_","")
    wnd_cp = wnd_cp_indiv(o.fn_DTS, o.fn_contigs, o.wnd_size)
    """
    outstr:
    chr start end indiv 0 0 0 color
    """
    
    c_out = output(o.contig_prefix, o.output_contigs) 
    for contig in wnd_cp.contigs:   
        print  >>stderr, contig
        
        cps = wnd_cp.get_cps_by_chr(contig)
        wnd_starts, wnd_ends = wnd_cp.get_wnds_by_chr(contig)
         
        prev_start = 0
        for i in xrange(0, cps.shape[0]-1):
            s, e = wnd_starts[i], wnd_ends[i]
            mid = (s+e)/2

            n_s, n_e = wnd_starts[i+1], wnd_ends[i+1]
            n_mid = (n_s+n_e)/2
                   
            end = (n_mid+mid)/2
            c_out.add(contig, prev_start, end, cps[i])
            prev_start = end
         
        e = wnd_ends[-1]
        c_out.add(contig, prev_start, e, cps[-1])
    
    fn_out = "%s/%s_%s.bb"%(o.out_dir, indiv, o.fn_out) 
    c_out.output(fn_out, indiv, o.fn_out)
                
     
