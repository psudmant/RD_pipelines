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

import numpy as np
import tables
import pygr.Data

from kitz_wssd.wssd_common import *
from wssd_pw_common import *


#def get_loc(msk,wssd,start,end,chr):


def create_DenseTrackSet(outTableFile,contigLengths,grp):

    print "creating  Dense Track Set..."
    track = DenseTrackSet(
                         contigLengths,
                         outTableFile,
                         overwrite=True, ###############################DANGERDANGERDANGER
                         openMode='w',
                         compression=True )


    ###only one group, called, mask
    track.addGroup( grp ) 
    track[grp].addArray( tables.BoolAtom(), [] )    
    
    print "done"
    return track

def set_regions(annot_locs,dts_annotation,grp):
    
    for contig in dts_annotation.mContigNameLen:
        
        print "setting contig:",contig,"..."
        load_contig = dts_annotation[grp][contig][:]

        starts = annot_locs[np.where(annot_locs['chr']==contig)]['start']
        ends = annot_locs[np.where(annot_locs['chr']==contig)]['end']
    
        if starts.shape[0] == 0: continue
    
        flip_to_true = np.unique(np.concatenate([np.arange(aa,bb) for (aa,bb) in zip(starts,ends)]))
        load_contig[flip_to_true] = 1    

        dts_annotation[grp][contig][:] = load_contig

############
#this funciton creates an n*3 array of the input annotations
#
###########
def load_positions(fn_input,chr_col,start_col,end_col,header,delim=None):
    
    start_col = int(start_col)
    end_col = int(end_col)
    chr_col = int(chr_col)
    
    nskiprows = (header and 1) or 0

    dt = np.dtype([('chr', np.str_, 8), ('start', np.uint32),('end',np.uint32)])

    print "getting columns %d:%d:%d"%(chr_col,start_col,end_col)
    positions = np.loadtxt(fn_input, dtype=dt, delimiter=delim, skiprows=nskiprows, usecols=(chr_col,start_col,end_col))    
    positions = np.sort(positions,axis=0,order=['chr'])
    return positions


if __name__=='__main__':
    
    opts = OptionParser()
    opts.add_option('','--outTable',dest='fnoutTable')
    opts.add_option('','--contigLengths',dest='fnContigLengths')
    #opts.add_option('','--overwrite',default=False,action='store_true',dest='overwrite')
    opts.add_option('','--pad',dest='pad',default=0)
    opts.add_option('','--grp',dest='grp',default="annotation")
    opts.add_option('','--input_annotations',dest='fn_input_annotations')
    opts.add_option('','--chr:start:end_cols',dest='chr_start_end', default="0:1:2")
    opts.add_option('','--header',default=False,action='store_true',dest='header')

    print "annotation"    

    (o, args) = opts.parse_args()

    chr_col,start_col,end_col = o.chr_start_end.split(":")    

    annot_locs = load_positions(o.fn_input_annotations,chr_col,start_col,end_col,o.header)

    annotation=create_DenseTrackSet(o.fnoutTable,o.fnContigLengths,o.grp)
    set_regions(annot_locs,annotation,o.grp)    

    


        
    #print mask_fns.values()

    



