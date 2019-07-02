import sys,os
import numpy as np
import scipy as sp
#import scipy.interpolate.fitpack2 as sp_i_fp2
import kg_file_handling as kgf
#from kg_file_handling import *
#from scipy import *
from optparse import OptionParser

from wssd_common import *
from wssd_pw_common import *

from calc_gc_depth import calc_gc_depth

import Bio.Statistics.lowess as biostats

def die(str):
    print str
    assert(0)

class GC_correction:

    def __init__(self,GC_width):
        #self.sum_read_depth = {}
        #self.n_locations = {}
        self.width_range = 2*GC_width+2

        self.total_sum_depth_at_gc = np.zeros(self.width_range)
        self.total_count_at_gc = np.zeros(self.width_range)

        #for i in range(2*GC_width+2):
        #    self.sum_read_depth[i] = 0
        #    self.n_locations[i] = 0

    #def create_equiv_map(self,mask_chunk,superDup_chunk,DGV_chunk,GC_count_chunk):

    def add_wssd_chunk(self,depth_chunk,mask_chunk,superDup_chunk,DGV_chunk,GC_count_chunk,genomic_gaps_chunk):

        cp2_unmasked_bool = np.cast[np.uint8](((mask_chunk) | (superDup_chunk) | (DGV_chunk) | genomic_gaps_chunk) ) #make sure NOT masked, not DGV or GC_count
        print "this many bases have a GC_count of 0...", ((GC_count_chunk==0)[cp2_unmasked_bool]).sum()
        sum,sum_depth_at_gc,count_at_gc = calc_gc_depth(depth_chunk,GC_count_chunk,cp2_unmasked_bool,self.width_range)

        #print sum_depth_at_gc
        #print count_at_gc

        self.total_sum_depth_at_gc += sum_depth_at_gc
        self.total_count_at_gc += count_at_gc

    def get_and_output_udepth_gc(self,fn):

        f = open(fn,'w')

        correction_tuple=kgf.get_GC_depth_correction_vect(self.total_sum_depth_at_gc,self.total_count_at_gc,GC_width)

        (GCp,ave_depths,GCp2_8,lowess_depth,correction,mu) = correction_tuple
        f.write("%d\n"%(GC_width))

        for i in range(len(GCp)):
            f.write("%f %f %d %d\n"%(GCp[i],correction[i],self.total_sum_depth_at_gc[i],self.total_count_at_gc[i]))

        return correction_tuple

    def output_udepth_gc(self,fn):
        f = open(fn,'w')
        f.write("%d\n"%(GC_width))

        for i in range(self.width_range):
            f.write("%f %d %d\n"%(float(i)/self.width_range,self.total_sum_depth_at_gc[i],self.total_count_at_gc[i]))

        f.close()

def calc_gc_correction(fn_wssd_file,lib_bac_analysis_dir,GC_width,lib_wssd,mask_track,gc_track,DGV_flat_track,superDup_flat_track,genomic_gaps_flat_track,dim):
    #GC correctino is calculated for the 2 copy bacs
    GC_corr = GC_correction(GC_width)

    for curr_contig in lib_wssd.mContigNameLen:
        if(curr_contig=="chrX" or curr_contig=="chrY" or curr_contig.find("random")!=-1 or curr_contig.find("alt")!=-1 or curr_contig.find("chrUn")!=-1): continue
        if(curr_contig=="X" or curr_contig=="Y" or curr_contig.find("random")!=-1): continue
        #if curr_contig != "1": continue
        #chunk_iter = 5e6
        chunk_iter = 5e7
        start  = 0
        contig_len = lib_wssd.mContigNameLen[curr_contig]

        while(start < contig_len):
            end = min(start+chunk_iter,contig_len)
            print "loading and analyzing: %s %d %d"%(curr_contig,start,end)
            mask_chunk = mask_track["mask"][curr_contig][start:end,:].sum(1)>0
            depth_chunk = np.cast['uint64'](lib_wssd.depth["wssd"][curr_contig][start:end,dim,:].sum(1))
            superDup_flat_chunk = superDup_flat_track["annotation"][curr_contig][start:end]
            DGV_flat_chunk = DGV_flat_track["annotation"][curr_contig][start:end]
            GC_count_chunk = gc_track["GC_content"][curr_contig][start:end]
            genomic_gaps_chunk = genomic_gaps_flat_track["annotation"][curr_contig][start:end]

            GC_corr.add_wssd_chunk(depth_chunk,mask_chunk,superDup_flat_chunk,DGV_flat_chunk,GC_count_chunk,genomic_gaps_chunk)
            start = end

    fn = "%s/cn2-gc-depth-WG-W%d-dim%d.txt"%(lib_bac_analysis_dir,GC_width,dim)
    GC_corr.output_udepth_gc(fn)

def get_lib_bac_wssd_pairs(genome,fn_wssd_dir,fn_bac_dir,wssd_file):

    fn_wssd_dir = "%s/%s"%(fn_wssd_dir,genome)
    fn_bac_dir = mkdir(fn_bac_dir,genome)
    lib_names  = os.listdir(fn_wssd_dir)
    #print fn_wssd_dir

    dirbac_fnwssd_pairs = {}
    for lib_name in lib_names:
        if(lib_name[0] == "_"): continue
        mkdir(fn_bac_dir,lib_name)
        newtuple = tuple(("%s/%s"%(fn_bac_dir,lib_name),"%s/%s/%s"%(fn_wssd_dir,lib_name,wssd_file)))
        dirbac_fnwssd_pairs[lib_name] = newtuple

    return dirbac_fnwssd_pairs

if __name__=='__main__':

    opts = OptionParser()
    opts.add_option('','--sex_pop_index',dest='fn_sex_pop_index')

    opts.add_option('','--hg_mask_file',dest='fn_hg_mask')
    opts.add_option('','--hg_contig_file',dest='fn_hg_contigs')
    opts.add_option('','--hg_gc_content_file',dest='fn_hg_gc_content_table')
    opts.add_option('','--wssd_file',dest='wssd_file',default="hg18rmsk.wssd")

    opts.add_option('','--GC_width',dest='GC_width',default=200)
    opts.add_option('','--fn_input_genomes',dest='fn_input_genomes')
    opts.add_option('','--fn_DGV_flatfile',dest='fn_DGV_flatfile')
    opts.add_option('','--fn_superdup_flatfile',dest='fn_superdup_flatfile')
    opts.add_option('','--fn_genomic_gaps_flatfile',dest='fn_genomic_gaps_flatfile')

    #opts.add_option('','--input_genome',dest='input_genome')

    (o, args) = opts.parse_args()

    input_genome_lines = open(o.fn_input_genomes,'r').readlines()
    print "getting genome_info..."
    genome_info = kgf.genome_info(o.fn_input_genomes,o.fn_sex_pop_index)
    print "done"

    mask_track = DenseTrackSet(o.fn_hg_contigs,
                                            fnWssd=o.fn_hg_mask,
                                            overwrite=False,
                                            openMode='r')

    gc_track = DenseTrackSet(o.fn_hg_contigs,
                                            fnWssd=o.fn_hg_gc_content_table,
                                            overwrite=False,
                                            openMode='r')

    DGV_flat_track = DenseTrackSet(o.fn_hg_contigs,
                                            fnWssd=o.fn_DGV_flatfile,
                                            overwrite=False,
                                            openMode='r')

    superDup_flat_track = DenseTrackSet(o.fn_hg_contigs,
                                            fnWssd=o.fn_superdup_flatfile,
                                            overwrite=False,
                                            openMode='r')

    genomic_gaps_flat_track = DenseTrackSet(o.fn_hg_contigs,
                                            fnWssd=o.fn_genomic_gaps_flatfile,
                                            overwrite=False,
                                            openMode='r')

    GC_width = int(o.GC_width)


    in_genomes = open(o.fn_input_genomes,'r').readlines()

    for in_genome in in_genomes:
        #(genome,fn_wssd_dir,fn_bac_dir,chunk_dir) = in_genome.split()
        (genome,fn_wssd_dir,fn_bac_dir,chunk_dir,genome_out) = in_genome.split()
        genome_id = genome.split(".")[0]

        dirbac_fnwssd_pairs = get_lib_bac_wssd_pairs(genome_id,fn_wssd_dir,fn_bac_dir,o.wssd_file)
        print genome
        print "*",dirbac_fnwssd_pairs

        for lib,bacdir_fnwssd_pair in dirbac_fnwssd_pairs.iteritems():

            print "\t",lib, bacdir_fnwssd_pair

            lib_wssd = WssdFile(o.fn_hg_contigs,
                                                    fnWssd=bacdir_fnwssd_pair[1],
                                                    overwrite=False,
                                                    openMode='r')
            dims = [0,1]

            for dim in dims:

                fn_cn2_gc_depth = "%s/cn2-gc-depth-WG-W%d-dim%d.txt"%(bacdir_fnwssd_pair[0],GC_width,dim)
                #if(os.path.exists(fn_cn2_gc_depth)):continue

                calc_gc_correction(bacdir_fnwssd_pair[1],bacdir_fnwssd_pair[0],GC_width,lib_wssd,mask_track,gc_track,DGV_flat_track,superDup_flat_track,genomic_gaps_flat_track,dim)

            lib_wssd.tbl.close()

