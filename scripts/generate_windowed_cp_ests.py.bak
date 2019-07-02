import os,sys
import os.path
import cPickle
import pickle 
from optparse import OptionParser
from wssd_common_v2 import *
from wssd_pw_common import *
from  ml_get_cpn import *
from sys import stderr

if __name__=='__main__':
    opts = OptionParser()
    opts.add_option('','--contig_file',dest='fn_contig_file')
    opts.add_option('','--mask_file',dest='fn_mask')
    opts.add_option('','--genome',dest='in_genome') #INDEX FILE
    opts.add_option('','--wnd_width',dest='wnd_width',default=None,type=int)
    opts.add_option('','--wnd_pickle',dest='wnd_pickle',default=None)
    opts.add_option('','--wnd_contig_file',dest='fn_wnd_contig_file')
    opts.add_option('','--param_file',dest='fn_param',default='_filtered_libs_analysis.GC3v2-NOWM-pad36/fit_params-CALKAN-NOWM-pad36-simple-dim0-ed-all.per_bac_params')
    opts.add_option('','--combined_corrected_name',dest='fn_comb_corr',default='wssd.combined_corrected.GC3.v2')
    opts.add_option('','--out_prefix',dest='fn_out_prefix',default=None)
    opts.add_option('','--sunk_based',dest='sunk_based',action='store_true',default=False)
    
    (o, args) = opts.parse_args()
    
    in_genome=open(o.in_genome,"r").readline() 
    (in_genome,in_wssd_dir,in_bac_dir,in_chunk_dir,in_primary_analysis_dir) = in_genome.split()
    
    mask = DenseTrackSet (o.fn_contig_file,
                                                o.fn_mask,
                                                overwrite=False,
                                                openMode='r' )
    print>>stderr, "%s/%s"%(o.fn_out_prefix,in_genome)    
    out_wnd_DTS = DenseTrackSet(o.fn_wnd_contig_file,
                                                            "%s_%s"%(o.fn_out_prefix,in_genome),
                                                            overwrite=True,
                                                            openMode='w')
    
    out_wnd_DTS.addGroup("copy")
    out_wnd_DTS.addGroup("starts")
    out_wnd_DTS.addGroup("ends")
    out_wnd_DTS['copy'].addArray(tables.Float32Atom(),[])
    out_wnd_DTS['starts'].addArray(tables.UInt32Atom(),[])
    out_wnd_DTS['ends'].addArray(tables.UInt32Atom(),[])
    ###WE ONLY NEED THE STARTS because start[k],start[k+1] == start[k], end[k]
    
    print >>stderr, "input genome: %s"%in_genome
    print >>stderr,"loading regions..." 
    F_region_pickle=open(o.wnd_pickle,"rb")
    regions_chrms,regions_coords,regions_wnds = cPickle.load(F_region_pickle) 
    print >>stderr,"done "

#read in base genome data
    print >>stderr,"getting...copies"

    g_data=genome_data(in_genome,"%s/%s/combined_corrected_wssd/%s"%(in_primary_analysis_dir,
                                                                        in_genome,o.fn_comb_corr),
                                "wssd.combined_corrected",
                                "%s/%s/%s"%(in_bac_dir,in_genome,o.fn_param),
                                o.fn_contig_file,
                                sunk_based = o.sunk_based)
    base_regressions=[]
    
    curr_chr=None
    curr_wnd_bin=0

    for xi in xrange(len(regions_chrms)):
        chr=regions_chrms[xi]
        start,end= regions_coords[xi]
        wnds = regions_wnds[xi]
        
        #if chr!="chr20": continue

        if curr_chr!=chr:
            curr_chr=chr
            curr_wnd_bin=0

        print "REGION",chr,start,end,wnds.shape, o.wnd_width
        regressions=g_data.get_regression(chr,start,end,wnds,mask,int(o.wnd_width))
        l_regressions=regressions.shape[0]
        #print regressions
        out_wnd_DTS['starts'][chr][curr_wnd_bin:curr_wnd_bin+l_regressions] = wnds[:,0]+start
        out_wnd_DTS['ends'][chr][curr_wnd_bin:curr_wnd_bin+l_regressions] = wnds[:,1]+start
        out_wnd_DTS['copy'][chr][curr_wnd_bin:curr_wnd_bin+l_regressions] = regressions
        #out_wnd_DTS.tbl.flush()
        curr_wnd_bin+=l_regressions



