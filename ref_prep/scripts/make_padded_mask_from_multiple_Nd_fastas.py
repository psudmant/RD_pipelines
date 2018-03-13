from optparse import OptionParser
from pygr import *

from pygr import worldbase
import pygr
from pygr.seqdb import SequenceFileDB

import numpy as np
from kitz_wssd.wssd_common import *
from wssd_pw_common import *


def pad(masked,pad):
    masked_locs = arrayToNonZeroContigRanges(masked)
    contig_len = masked.shape[0]
    starts = masked_locs[:,0]-pad
    ends = masked_locs[:,1]+pad

    starts[np.where(starts<0)] = 0
    ends[np.where(ends>contig_len)] = contig_len
    if starts.shape[0]>0:
        flip_to_true = np.unique(np.concatenate([np.arange(aa,bb) for (aa,bb) in zip(starts,ends)]))
        masked[flip_to_true] = 1

    return masked

if __name__=="__main__":
    opts = OptionParser()
    opts.add_option('','--input_fa',dest='fn_fa',default=None)
    opts.add_option('','--input_RM_fa',dest='fn_RM_fa',default=None)
    opts.add_option('','--input_TRF_fa',dest='fn_TRF_fa',default=None)
    opts.add_option('','--fn_mask_out',dest='fn_mask_out',default=None)
    opts.add_option('','--contigs',dest='fn_contigs',default=None)
    opts.add_option('','--pad',dest='pad',default=30,type=int)
    opts.add_option('','--old_to_new_contigs',dest='old_to_new_contigs')
    opts.add_option('','--limit_to_contigs',dest='fn_limit_to_contigs',default=None)
    opts.add_option('','--exclude_TRF_from',dest='fn_exclude_TRF_from',default=None)

    (o,args) = opts.parse_args()
    #####RIGHT NOW THIS ONLY ACCEPTS FASTAS ->BUT< could take coords too

    #######NOTE< kindof a cop-out, i'm just setting the one track (out of 3) to be the whole mask... ignoring the different
    #######levels for RM, TRF, etc. could change later 
 #old_to_new_contigs = o.old_to_new_contigs!= None and dict([tuple([l.rstrip().split()[0],l.rstrip().split()[1]]) for l in open(o.old_to_new_contigs,'r').readlines()]) or None

    limit_to_contigs = None; 
    if o.fn_limit_to_contigs != None:    
         limit_to_contigs = [line.rstrip() for line in open(o.fn_limit_to_contigs).readlines() ]
     
    exclude_TRF_from = None; 
    if o.fn_exclude_TRF_from != None:    
         exclude_TRF_from = [line.rstrip() for line in open(o.fn_exclude_TRF_from).readlines() ]

    input_fa = SequenceFileDB(o.fn_fa)
    input_fa_RM = SequenceFileDB(o.fn_RM_fa)
    input_fa_TRF = SequenceFileDB(o.fn_TRF_fa)

    mask_track = DenseTrackSet(
                         o.fn_contigs,
                         o.fn_mask_out,
                         True,
                         'w',
                         compression=True )


    ###only one group, called, mask
    grp = "mask"
    mask_track.addGroup( grp ) 
    mask_track[grp].addArray( tables.BoolAtom(), [3] )    

    for contig in input_fa:
        print contig
        contig = contig.replace(".","_")
        #fa_contig=old_to_new_contigs!=None and old_to_new_contigs[contig] or contig
        char_fa_RM = np.array(str(input_fa_RM[contig]),'c')
        char_fa_TRF = np.array(str(input_fa_TRF[contig]),'c')
        print char_fa_RM[0:100]
        print char_fa_TRF[0:100]

        if limit_to_contigs == None or contig in limit_to_contigs:
            if  exclude_TRF_from !=None and contig in exclude_TRF_from:
                masked = (char_fa_RM=="N")
            else:
                masked = (char_fa_RM=="N")|(char_fa_TRF=="N")

            masked = pad(masked,o.pad)
        else:
            print "IGNORING CURRENT CONTIG", contig
            masked = (char_fa_RM=="?")|(char_fa_TRF=="?") ###SHOULD EVALUATE TO ALL FALSES!
            
        print masked[0:100]

        ###here, really should be putting the RM in 0, and the TRF in 1... or, other way around, can't remember :)
        mask_track['mask'][contig][:,0] = masked    

