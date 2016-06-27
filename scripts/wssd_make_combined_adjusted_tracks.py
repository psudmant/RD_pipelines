import sys
import os
import os.path
from optparse import OptionParser
from collections import defaultdict
import tempfile
import time

import numpy as np
import tables
import pygr.Data

from wssd_common import *
from wssd_pw_common import *
import kg_file_handling as kgf
#def get_loc(msk,wssd,start,end,chr):

def mkdir(dir,file):

  ls_dir = os.listdir(dir)
  if(not(file_exists(ls_dir,file))):
    command = "mkdir %s/%s"%(dir,file)
    print command
    os.system(command)

def file_exists(ls,file):
  for f in ls:
    if(f==file):
      return 1
  return 0


def in_corrupted(genome):
	corrupted_genomes = open("/net/gs/vol1/home/psudmant/EEE_Lab/1000G/analysis_files/input_genomes/corrupted").readlines()
	for g in corrupted_genomes:
		g = g.rstrip()
		if(g==genome):
			return True
	
	return False
	

def get_passing_libs(bac_analysis_lib_dir,max_correction,overide_thresh):
	
	passing_libs = {}

	libls = os.listdir(bac_analysis_lib_dir)
	for lib in libls:
		if(lib[0:1]=="_" or lib[0] == "."):
			continue
		lib_qc_dir = "%s/%s"%(bac_analysis_lib_dir,lib)
		if(kgf.lib_pass_qc(lib_qc_dir,overide_thresh)):
			print "lib passed: %s"%(lib)
			gc_correction_vects = []	
			#gc_correction_vects.append(load_gc_correction_vect(lib_qc_dir,0))
			#gc_correction_vects.append(load_gc_correction_vect(lib_qc_dir,1))
			gc_correction_vects.append(kgf.get_GC_depth_correction_from_dir(lib_qc_dir,max_correction))
			#print gc_correction_vects[0]

			passing_libs[lib] = gc_correction_vects

	return passing_libs

def load_wssds(passing_libs,lib_wssd_dir,fnContigLengths,wssd_file_name,alt_input_contigs):

	if alt_input_contigs!=None:
		fnContigLengths = alt_input_contigs
	
	lib_wssds = {}	
	for lib,lib_correction in passing_libs.iteritems():
		print lib	
		fn_lib_wssd = "%s/%s/%s"%(lib_wssd_dir,lib,wssd_file_name)	
		#fn_lib_wssd = "%s/%s/hg18rmsk.wssd"%(lib_wssd_dir,lib)	
		print fn_lib_wssd
		lib_wssd  = WssdFile(fnContigLengths,
												fn_lib_wssd,
												overwrite=False,
												openMode='r')

		lib_wssds[lib] = lib_wssd
	return lib_wssds

def close_wssds(lib_wssds):
	for lib,lib_wssd in lib_wssds.iteritems():
		lib_wssd.tbl.close()						


if __name__=='__main__':
	
	opts = OptionParser()
	opts.add_option('','--contigLengths',dest='fnContigLengths')
	opts.add_option('','--gc_width',dest='gc_width')
	opts.add_option('','--inGC',dest='fngrpGC')
	opts.add_option('','--in_genomes',dest='fn_in_genomes')
	opts.add_option('','--no_correction',default=False,action='store_true',dest='no_correction')
	opts.add_option('','--overide_thresh',default=None,dest='overide_thresh')
	opts.add_option('','--sex_pop_index',dest='fn_sex_pop_index')
	opts.add_option('','--append_to_name',dest='i_append_to_name',default="")
	opts.add_option('','--max_correction',dest='max_correction',default=3)
	opts.add_option('','--input_wssd_file_names',dest='wssd_file_name',default="hg18rmsk")
	opts.add_option('','--alt_wssd_lambda',dest='alt_wssd_lambda',default=None)
	opts.add_option('','--alt_primary_analysis_lambda',dest='alt_primary_analysis_lambda',default=None)
	
	opts.add_option('','--alt_input_contigs',dest='alt_input_contigs',default=None)
	opts.add_option('','--map_alt_input_to_output_contigs',dest='map_alt_input_to_output_contigs',default=None)

	(o, args) = opts.parse_args()

	max_correction = int(o.max_correction)
	genome_info = kgf.genome_info(o.fn_in_genomes,o.fn_sex_pop_index)


	map_alt_input_to_output_contigs = None
	map_alt_output_to_input_contigs = None
	
	if(o.alt_input_contigs!=None):
		map_alt_input_to_output_contigs = dict([[line.split()[0],line.rstrip().split()[1]] for line in open(o.map_alt_input_to_output_contigs,'r').readlines()])	
		map_alt_output_to_input_contigs = dict([[line.rstrip().split()[1],line.rstrip().split()[0]] for line in open(o.map_alt_input_to_output_contigs,'r').readlines()])	

	if(o.no_correction):
		append_to_name = o.i_append_to_name + ".no_correction"		
		print ".............not applying correction.............."
	else:
		append_to_name = o.i_append_to_name

	####fnMsk,grpMask=o.fngrpMaskVec.split(':')[0], o.fngrpMaskVec.split(':')[1]

	fn_gc_table,gc_grp=o.fngrpGC.split(':')[0],o.fngrpGC.split(':')[1]
	GC_width = int(o.gc_width)

	in_genomes = open(o.fn_in_genomes,"r").readlines()

	GC = DenseTrackSet (o.fnContigLengths,
									fn_gc_table,
									overwrite=False,
									openMode='r' )

	for in_genome in in_genomes:
		(genome,fn_wssd_dir,fn_bac_dir,chunk_dir,genome_out) = in_genome.split()
		if(not(genome_info.genomes[genome].passed_qc)):
			print "genome failed"
			continue

		if o.alt_wssd_lambda !=None:
			fn_wssd_mod = eval(o.alt_wssd_lambda)
			fn_wssd_dir = fn_wssd_mod(fn_wssd_dir)

		if o.alt_primary_analysis_lambda:
			fn_primary_analysis_mod = eval(o.alt_primary_analysis_lambda)
			genome_out = fn_primary_analysis_mod(genome_out)	

		genome_orig_fn = genome.split(".")[0]
		wssd_lib_dir = "%s/%s"%(fn_wssd_dir,genome_orig_fn)
		bac_analysis_lib_dir = "%s/%s"%(fn_bac_dir,genome_orig_fn)
	
		print genome_out
		mkdir(genome_out,genome)
		fn_indiv_out = "%s/%s/"%(genome_out,genome)
		mkdir(fn_indiv_out,"combined_corrected_wssd")

		fn_corrected_depths_out = "%s/combined_corrected_wssd/wssd.combined_corrected%s"%(fn_indiv_out,append_to_name)	

		mkdir(fn_indiv_out,"output")

		if o.overide_thresh == None:
			overide_thresh = None
		else:
			overide_thresh = float(o.overide_thresh)

		passing_libs = get_passing_libs(bac_analysis_lib_dir,max_correction,overide_thresh)
		lib_wssds = load_wssds(passing_libs,wssd_lib_dir,o.fnContigLengths,o.wssd_file_name,o.alt_input_contigs)
	

		#1 get libraries
		#2 make new wssd track file
		print fn_corrected_depths_out
		adjusted_wssd_track = WssdFile(o.fnContigLengths,
                                       fnWssd=fn_corrected_depths_out,
                                       overwrite=True,
                                       openMode='w',
                                       compression=True,
                                       datatype=tables.Atom.from_dtype(np.dtype(np.float16))
                                      )
		
		grpName = "wssd.combined_corrected"

		#!#!#!
		adjusted_wssd_track.addTrackSet(grpName)
		

		#50M*3*2
		for chr in GC.mContigNameLen:
			chunk_iter = 5e7
			currChromLen = adjusted_wssd_track.mContigNameLen[chr]
			chunkStart = 0
			chunkEnd = min(chunkStart+chunk_iter,currChromLen)

			while chunkStart < currChromLen-1:
				print "working on contig:%s:%d:%d..."%(chr,chunkStart,chunkEnd)
				curr_working_contig = np.zeros((chunkEnd-chunkStart,3,2))		
				for lib,gc_correction_vects in passing_libs.iteritems():
					print "getting depths (and starts!)..."

					lib_wssd_chr=chr
					if o.alt_input_contigs != None:
						lib_wssd_chr=map_alt_output_to_input_contigs[chr]

					currlib_depths_and_starts = lib_wssds[lib].depth["wssd"][lib_wssd_chr][chunkStart:chunkEnd,:,:]
					for k in range(2): #for dim 0 and 1
						print "getting gc correction lib %s ..."%( lib)
						#HARD CODED gc_correction_vects[0] because we use the gc correction from the 
						#
						adjustment = gc_correction_vects[0].take(GC[gc_grp][chr][chunkStart:chunkEnd])
						sum = currlib_depths_and_starts[:,:,k].sum(0).sum()	
						sum_sanity_check = 0
						for i in range(3):
							#w = float(currlib_depths_and_starts[:,i,k].sum())/sum
							#sum_sanity_check+=w
							if(o.no_correction):
								curr_working_contig[:,i,k]+=currlib_depths_and_starts[:,i,k]
							else: #DEFAULT
								curr_working_contig[:,i,k]+=currlib_depths_and_starts[:,i,k]*adjustment
								#curr_working_contig[:,i,k]+=currlib_depths_and_starts[:,i,k]*adjustment*w
						print "sum sanity check: ", sum_sanity_check					
				
				print "curr_working_contig shape " ,curr_working_contig.shape
				#np.clip(curr_working_contig,0,1e30,curr_working_contig)
				print "curr_working_contig shape " ,curr_working_contig.shape
				adjusted_wssd_track.depth[grpName][chr][chunkStart:chunkEnd,:,:] = curr_working_contig
				
				chunkStart+=chunk_iter
				chunkEnd = min(chunkStart+chunk_iter,currChromLen)


		#3 make adjustment track
		#4 add adjustment track to each input accordingly
		close_wssds(lib_wssds)

