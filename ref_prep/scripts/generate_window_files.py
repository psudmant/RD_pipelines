import os,sys
import os.path
import cPickle
import pickle 
import numpy as np
from optparse import OptionParser
import wssd_common
import wssd_pw_common
from sys import stderr
from collections import defaultdict


def get_sliding_regions(mask,wnd_width,slide_by, sunk_based=False):
    
    assert wnd_width%slide_by == 0
    
    slides_in_wnd = wnd_width/slide_by
    regions_chrms,regions_coords,regions_wnds = [],[],[]
    
    for chr,l in mask.mContigNameLen.iteritems():
        if "random" in chr: continue
        print chr 
        if sunk_based:
            curr_mask = mask['isntSunkOrIsMasked'][chr][:]
        else:
            curr_mask = mask['mask'][chr][:,:].sum(1)>0
        chr_wnd_bnds,max_wnd = wssd_pw_common.getBoundsForEqualWeightedWindows(curr_mask,0,l,slide_by)

        regions_chrms.append(chr)
        regions_coords.append(tuple([0,l]))
        #because you are sliding - you now need to offset the ENDs
        #note, the last few wnds are increasingly smaller..
        if chr_wnd_bnds is not None:
            print chr_wnd_bnds[0:200]
            chr_wnd_bnds[0:-(slides_in_wnd-1),1] = chr_wnd_bnds[(slides_in_wnd-1):,1]
            chr_wnd_bnds[-(slides_in_wnd-1):,1] = chr_wnd_bnds[-1:,1]

            # BN 2018: Fix off by one error with single window contigs
            if len(chr_wnd_bnds) == 1 and chr_wnd_bnds[0][1] == l:
                chr_wnd_bnds[0][1] = l - 1
            regions_wnds.append(chr_wnd_bnds)
            print chr_wnd_bnds[0:200] 
            del curr_mask
        else:
            regions_wnds.append(np.array([[0, l-1]]))
        
    return regions_chrms,regions_coords,regions_wnds

if __name__=='__main__':
    opts = OptionParser()
    opts.add_option('','--contig_file',dest='fn_contig_file')
    opts.add_option('','--mask_file',dest='fn_mask')
    opts.add_option('','--base_genome',dest='base_genome') #INDEX FILE
    opts.add_option('','--wnd_width',dest='wnd_width',default=1000,type=int)
    opts.add_option('','--slide_by',dest='slide_by',default=0,type=int)
    opts.add_option('','--in_regions',dest='fn_in_regions',default=None) #NONE
    opts.add_option('','--region_output',dest='region_out',default=None)
    opts.add_option('','--sunk_based',dest='sunk_based',action='store_true',default=False)
    
    (o, args) = opts.parse_args()
    mask = wssd_common.DenseTrackSet (o.fn_contig_file,
                                      o.fn_mask,
                                      overwrite=False,
                                      openMode='r' )

    print "Determining Regions..." 
    if o.slide_by!=0:
        regions_chrms,regions_coords,regions_wnds = get_sliding_regions(mask,o.wnd_width,o.slide_by,sunk_based=o.sunk_based)
    else:
        regions_chrms,regions_coords,regions_wnds = get_regions(None,mask,wnd_width=int(o.wnd_width),sunk_based=o.sunk_based)
    print "Regions Determined..."
    print "Number of region chunks  %d..."%len(regions_chrms)
    
    n_winds_by_chr=defaultdict(int)
    for chunk_i in xrange(len(regions_chrms)):
        print >>stderr,regions_chrms[chunk_i],regions_wnds[chunk_i].shape[0]
        n_winds_by_chr[regions_chrms[chunk_i]]+=regions_wnds[chunk_i].shape[0]

    #read in base genome data
        
    output_regions_contigs=open("%s.contigs"%o.region_out,"w")
    for contig, len in n_winds_by_chr.iteritems():
        output_regions_contigs.write("%s\t%d\n"%(contig,len))
    output_regions_contigs.close()

    output_regions=open(o.region_out,"wb")
    cPickle.dump([regions_chrms,regions_coords,regions_wnds], output_regions)
    output_regions.close()
