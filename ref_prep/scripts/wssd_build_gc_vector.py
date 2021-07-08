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
from pygr import seqdb

from kitz_wssd.wssd_common import *
import scipy.signal as sp_sig


#def get_loc(msk,wssd,start,end,chr):

def get_chr_correction(DB,chr,GC_width):	
	
		print "generating %s GC contents len:%d"%(chr,len(DB[chr]))
		seqData = DB[chr]
		seqStr = str(seqData).upper()		
		char_seqStr = np.array(seqStr,'c')
		gc_barray	 = (char_seqStr=='C')|(char_seqStr=='G')
		GC_conv = np.ones(2*GC_width+1)

		if gc_barray.shape[0] < GC_conv.shape[0]:
			p_conv = float(np.sum(gc_barray))/gc_barray.shape[0]
			i_conv = round(p_conv*(2*GC_width+1))
			conv = np.ones(gc_barray.shape)*i_conv
			print 'n gcs:', np.sum(gc_barray)
			print 'frac gcs:', p_conv
			print 'expected in 401:', i_conv
			print 'expected frac:', float(i_conv)/401
			return conv
		else:
			conv = np.convolve(gc_barray,GC_conv,'same') 
			GC_count = conv
			conv[:GC_width] = conv[GC_width+1]
			conv[(conv.shape[0]-GC_width-1):] = conv[conv.shape[0]-GC_width-1]
			#GC_percent = conv/(2*GC_width+1)
			print "done"
			return conv

def create_gc_DenseTrackSet(outTableFile,contigLengths,grp,overwrite_b):

	
	print "creating GC Dense Track Set..."
	GC = DenseTrackSet(
						 contigLengths,
						 outTableFile,
						 overwrite_b,
						 'w',
						 compression=True )

	GC.addGroup( grp )
	GC[grp].addArray( tables.UInt16Atom(), [] )
	print "done"
	return GC

if __name__=='__main__':
	
	opts = OptionParser()
	opts.add_option('','--outTable',dest='fnoutTable')
	opts.add_option('','--contigLengths',dest='fnContigLengths')
	opts.add_option('','--gc_width',dest='gc_width')
	opts.add_option('','--overwrite',default=True,action='store_false',dest='overwrite')
	opts.add_option('','--pygr_sequence',dest='pygr_seq')
	opts.add_option('','--fastq_sequence',dest='fn_fastq_seq',default=None)
	
	(o, args) = opts.parse_args()
	
	GC_width = int(o.gc_width)

	if(o.pygr_seq != None):
		if(o.pygr_seq.upper() == "HG18"):
			seqs = pygr.Data.Bio.Seq.Genome.Human.hg18()
		elif(o.pygr_seq.upper() == "HG19"):
			seqs = pygr.Data.Bio.Seq.Genome.Human.hg19()
		elif(o.pygr_seq.upper() == "CHIMPY"):
			seqs = pygr.Data.Bio.Seq.Genome.chimp.chrY()
		elif(o.pygr_seq.upper() == "BACS" or o.pygr_seq.upper() == "CONTROL_BACS"):
			seqs = pygr.Data.Bio.Seq.Genome.Human.control_bacs()
	elif(o.fn_fastq_seq!=None):
		seqs = seqdb.SequenceFileDB(o.fn_fastq_seq)
	else:
		print "no sequence file input... exiting"
		sys.exit(1)

	GC_content = {}

	grp = "GC_content"
	GC_DT = create_gc_DenseTrackSet(o.fnoutTable,o.fnContigLengths,grp,o.overwrite)


	for contig in seqs:
		if(contig in GC_DT[grp]):
			print "loading %s..."%(contig)
			#if(chr=="chr21"):
			#GC_content["chr21"] = get_chr_correction(hg18,"chr21",GC_width)
			#print GC_content["chr21"][15e6:15.01e6]
			GC_depth = get_chr_correction(seqs,contig,GC_width)
			print len(GC_depth)
			print GC_depth
			print GC_DT[grp][contig][:].shape
			print GC_DT[grp][contig][:]
			GC_DT[grp][contig][:] = GC_depth
			print len(GC_DT[grp][contig][:])
		else:
			print contig

	del GC_DT

######################
######################
######################
######################
######################
######################
######################
if __name__=='__NULL__':    
	opts = OptionParser()
	opts.add_option('','--inWssd',dest='fngrpWssd')
	opts.add_option('','--inMaskVec',dest='fngrpMaskVec')
	opts.add_option('','--contigLengths',dest='fnContigLengths')


	fnMsk,grpMask=o.fngrpMaskVec.split(':')[0], o.fngrpMaskVec.split(':')[1]
	fnWssd,grpWssd=o.fngrpWssd.split(':')[0], o.fngrpWssd.split(':')[1]
	msk = DenseTrackSet(
                  o.fnContigLengths,
                  fnMsk,
                  overwrite=False,
                  openMode='r',
                  compression=False )
    
	wssd = WssdFile( o.fnContigLengths,
                     fnWssd,
                     overwrite = False,
                     openMode = 'r' )        
		
