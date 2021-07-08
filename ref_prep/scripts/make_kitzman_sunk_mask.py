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
    windows = [np.arange(aa,bb) for (aa,bb) in zip(starts,ends)]
    if len(windows) > 0:
        flip_to_true = np.unique(np.concatenate(windows))
        masked[flip_to_true] = 1
    return masked

if __name__=="__main__":
    opts = OptionParser()
    opts.add_option('','--input_fa',dest='fn_fa',default=None)
    opts.add_option('','--input_RM_fa',dest='fn_RM_fa',default=None)
    opts.add_option('','--input_TRF_fa',dest='fn_TRF_fa',default=None)
    opts.add_option('','--input_sunk_locs',dest='fn_sunk_locs',default=None)
    opts.add_option('','--fn_sunkmask_out',dest='fn_sunkmask_out',default=None)
    opts.add_option('','--contigs',dest='fn_contigs',default=None)
    opts.add_option('','--pad',dest='pad',default=36,type=int)
    opts.add_option('','--override_all_sunk',dest='all_sunk',default=False,action='store_true')
    
    opts.add_option('','--alt_assembly_to_curr_contig',dest='fn_alt_assembly_to_curr',default=None)
    opts.add_option('','--alt_assembly_sunk_track',dest='fn_alt_sunk_track',default=None)
    opts.add_option('','--alt_assembly_contigs',dest='fn_alt_contigs',default=None)
    opts.add_option('', '--starts_only', dest='starts_only', action='store_true', default=False)

    (o,args) = opts.parse_args()
    #####RIGHT NOW THIS ONLY ACCEPTS FASTAS ->BUT< could take coords too

    input_fa = SequenceFileDB(o.fn_fa)
    input_fa_RM = SequenceFileDB(o.fn_RM_fa)
    input_fa_TRF = SequenceFileDB(o.fn_TRF_fa)
    
    #GET SUNK LOCS
    dt = np.dtype([('contig', np.str_, 128), ('start', np.int64),('end',np.int64)])
    sunk_positions = np.loadtxt(o.fn_sunk_locs, dtype=dt, delimiter="\t")

    print sunk_positions
    #positions = np.sort(positions,axis=0,order=['chr'])

    #DONE

    sunkmask_track = DenseTrackSet(
                         o.fn_contigs,
                         o.fn_sunkmask_out,
                         True,
                         'w',
                         compression=True )
    
    alt_assembly_to_curr = None    
    alt_assembly_sunkmask = None    

    if o.fn_alt_assembly_to_curr != None:
        print "BROKEN! DON'T KNOW WHY!"
        exit()
        alt_assembly_sunkmask = DenseTrackSet(
                                 o.fn_alt_contigs,
                                 o.fn_alt_sunk_track,
                                 False,
                                 'r',
                                 compression=True )
        alt_assembly_to_curr = {}
        for line in open(o.fn_alt_assembly_to_curr,'r'):
            alt_chr,alt_start,alt_end,curr_chr,curr_start,curr_end = line.rstrip().split()
            alt_start,alt_end,curr_start,curr_end = int(alt_start),int(alt_end),int(curr_start),int(curr_end)
            alt_assembly_to_curr[tuple([alt_chr,alt_start,alt_end])] =  tuple([curr_chr,curr_start,curr_end])

    #######groups
    #/isntSunkOrIsMasked/
    #/isUnmasked/
    #/isSunk/

    grps = ['isntSunkOrIsMasked','isUnmasked','isSunk']

    for i in xrange(3):
        sunkmask_track.addGroup(grps[i])
        sunkmask_track[grps[i]].addArray(tables.BoolAtom(), [])


    for contig in input_fa:
        print contig

        char_fa_RM = np.array(str(input_fa_RM[contig]),'c')
        char_fa_TRF = np.array(str(input_fa_TRF[contig]),'c')
        masked = (char_fa_RM=="N")|(char_fa_TRF=="N")
        masked = pad(masked,o.pad)


        sunk_starts = np.sort(sunk_positions[np.where(sunk_positions["contig"]==contig)]['start'])
        sunk_ends = np.sort(sunk_positions[np.where(sunk_positions["contig"]==contig)]['end'])
        if o.all_sunk:
            sunk_vect = np.ones(masked.shape[0])
        elif alt_assembly_to_curr != None:
            sunk_vect = np.zeros(masked.shape[0])
            for alt,curr in     alt_assembly_to_curr.iteritems():
                alt_chr,alt_start,alt_end = alt            
                curr_chr,curr_start,curr_end = curr
                if curr_chr != contig: continue
                
                sunk_vect[curr_start:curr_end] = alt_assembly_sunkmask['isSunk'][alt_chr][alt_start:alt_end]
        else:
            sunk_vect = np.zeros(masked.shape[0])
        if o.starts_only:
            sunk_vect[w_sunks] = 1
        else:
            for s, e in zip(sunk_starts, sunk_ends):
                sunk_vect[s:e] = 1
        sunk_vect = sunk_vect>0
    
        sunkmask_track['isntSunkOrIsMasked'][contig][:] = ~sunk_vect|masked    
        sunkmask_track['isUnmasked'][contig][:] = masked    
        sunkmask_track['isSunk'][contig][:] = sunk_vect    

