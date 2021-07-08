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
import math

from sys import stderr

import Bio.Statistics.lowess as biostats
import kg_file_handling as kgf

from wssd_common import *
from  ml_get_cpn import *

#def get_loc(msk,wssd,start,end,chr):

def isnan(x):
    return math.isnan(x)

def file_exists(ls,file):
    for f in ls:
        if(f==file):
            return 1
    return 0

def mkdir(dir,file,test=False):
    
    if test:
        return "%s/%s"%(dir,file)

    ls_dir = os.listdir(dir)
    if(not(file_exists(ls_dir,file))):
        command = "mkdir %s/%s"%(dir,file)
        os.system(command)
    return "%s/%s"%(dir,file)

class genome_data:

    def __init__(self,genome,fn_wssd,grp,fn_fit_params,fnContigLengths,sunk_based=False,lUseCp=[1,2,4,6]):
    
        self.wssd  = WssdFile(fnContigLengths,
                                                fn_wssd,
                                                overwrite=False,
                                                openMode='r')
        self.grp = grp 

        self.sunk_based = sunk_based

        print(fn_fit_params)
        print(fn_wssd)
        self.fn_wssd = fn_wssd
        print("**********")
        with open(fn_fit_params,'r') as f:
            mtxParams = np.array([ [float(x) for x in line.replace('\n','').strip().split()[1:]]
                              for line in f],'float64')
        cps=mtxParams[:,0]
        mus=mtxParams[:,1]
        sigs=mtxParams[:,2]
        nbps=mtxParams[:,3]

        print(cps)
        print(mus)
        print(sigs)
        print(nbps)
        
        self.mtxParams = mtxParams
        thruZero=True
        #lUseCp = [1,2,4,6]
        #lUseCp = [1,2]
        self.regressMu = LinearRegress( cps, mus, nbps )
        self.regressMu.fit(lUseCp, thruZero)

        self.cns = np.arange(300).astype('float64')
        self.cns[0] = 0.05 # for continuous poisson only, having a zero class implies a log(0) 
                                             # picked 0.05 somewhat arbitrarily to reflect the fact that deleted regions
                                             # can get hits from repetitive elts
        #self.regress_cn = regressMu.regress(self.cns)

    def get_regression(self,chr,start,end,wnds,mask,nbpWidth,mask_chr=None):
        
        
        if self.sunk_based:
            sdepths_working_chunk = np.nan_to_num(self.wssd.depth[self.grp][chr][start:end,0,1]).astype(np.float64)
        else:
            sdepths_working_chunk = np.nan_to_num(self.wssd.depth[self.grp][chr][start:end,:,0]).astype(np.float64).sum(1)

        
        mask_chr = mask_chr==None and chr or mask_chr
        #IF IT's A SUNK MASK -> different 
        if self.sunk_based:
            masked_contig = mask['isntSunkOrIsMasked'][mask_chr]
            if masked_contig.ndim == 2:
                masked_region = mask['isntSunkOrIsMasked'][mask_chr][start:end,:]
            elif masked_contig.ndim == 1:
                masked_region = mask['isntSunkOrIsMasked'][mask_chr][start:end]
        else:
            masked_region = mask['mask'][mask_chr][start:end,:].sum(1)>0

        #if sum(~masked_region) < nbpWidth:
        #    return -1

        sdepths_working_chunk[np.where(masked_region==1)] = 0
        
        csum_sdepths = np.cumsum(sdepths_working_chunk)[wnds[:,1]]
        csum_sdepths_shift = np.cumsum(sdepths_working_chunk)[wnds[:,0]]

        """ 
        ORIGINAL CODE
        THIS DOESN'T THINK ABOUT THE STARTS AT ALL! 
        assumes windows abut each other
        """
        #csum_sdepths = np.cumsum(sdepths_working_chunk)[wnds[:,1]-1]
        #csum_sdepths_shift = np.r_[0,csum_sdepths[0:csum_sdepths.shape[0]-1]]
        depths = csum_sdepths-csum_sdepths_shift
        
        regression = self.regressMu.regrinv(depths.astype('float32')/nbpWidth)
        return regression

    def get_cp_over_locus(self,chr,start,end,wnds,mask,nbpWidth):
        
        regr=self.get_regression(chr,start,end,wnds,mask,nbpWidth)
        #print chr,start,end
        #print regr

        return np.median(regr)

    def get_single_bp_regression(self,chr,start,end):
        if self.sunk_based:
            return self.regressMu.regrinv(np.nan_to_num(self.wssd.depth[self.grp][chr][start:end,0,1]).astype(np.float64))
        else:
            return self.regressMu.regrinv(np.nan_to_num(self.wssd.depth[self.grp][chr][start:end,:,0]).astype(np.float64).sum(1))



def get_wnds_mean(wssd,chr,start,end,wnds,mask,wnd_width,sunk_based=False):
    
        if sunk_based:
            sdepths_working_chunk = np.nan_to_num(wssd.depth["wssd.combined_corrected"][chr][start:end,0,1]).astype(np.float64)
        else:
            sdepths_working_chunk = np.nan_to_num(wssd.depth["wssd.combined_corrected"][chr][start:end,:,0]).astype(np.float64).sum(1)

        #IF IT's A SUNK MASK -> different 
        if sunk_based:
            masked_region = mask['isntSunkOrIsMasked'][chr][start:end]
        else:
            masked_region = mask['mask'][chr][start:end,:].sum(1)>0

        #    return -1
        print(sdepths_working_chunk)
        sdepths_working_chunk[np.where(masked_region==1)] = 0
        #print sdepths_working_chunk.shape
        #print sdepths_working_chunk
        #print wnds
        csum_sdepths = np.cumsum(sdepths_working_chunk)[wnds[:,1]-1]

        n_unmasked_bp = np.sum(~masked_region)

        csum_sdepths_shift = np.r_[0,csum_sdepths[0:csum_sdepths.shape[0]-1]]
        depths = csum_sdepths-csum_sdepths_shift
        depths = depths.astype('float32')/wnd_width
        return depths

def hash_bed_val_by_loc(input_bed,val_loc=3,dict_vals=False,alt_key=None):

    hash_ret = {}
    header = None
    F_in = open(input_bed,'r')

    if dict_vals:
        header = F_in.readline()
        
        for line in F_in.readlines():
            sheader = header.split()
            sline = line.split()
            chr,start,end = sline[0:3]
            c_dict = {}
            for i in range(len(sheader)):
                c_dict[sheader[i]] = sline[i] 
            
            if alt_key!=None:
                key = sline[alt_key]
            else:
                key = "%s:%s-%s"(chr,start,end)
            hash_ret[key] = c_dict

        return hash_ret,header
    else:
        for line in F_in.readlines():
            chr,start,end = line.split()[0:3]
            val = line.split()[val_loc]
            
            if alt_key!=None:
                key = sline[alt_key]
            else:
                key = "%s:%s-%s"(chr,start,end)
            hash_ret[key] = val
        return hash_ret


#####
##This functino creates tuples from fields while reading in a file
#optionally you can specify a filter with a lambda function
def read_tuples(fn_tuples,t_field1_s,t_field1_e,t_field2_s,t_field2_e,weight_field=None,splitter=None,test_func='lambda line: True',joiner=":"):

    #print fn_tuples
    f = eval(test_func)
    if weight_field!=None:
        tuples = [(joiner.join(line.rstrip().split(splitter)[t_field1_s:t_field1_e]),joiner.join(line.rstrip().split(splitter)[t_field2_s:t_field2_e]),{'weight':line.rstrip().split(splitter)[weight_field]}) for line in open(fn_tuples).readlines() if f(line)]
        
    else:
        tuples = [(joiner.join(line.rstrip().split(splitter)[t_field1_s:t_field1_e]),joiner.join(line.rstrip().split(splitter)[t_field2_s:t_field2_e])) for line in open(fn_tuples).readlines() if f(line)]

    return tuples


def read_table_no_header(fn_input,header,hash=False,hash_off_s=0,hash_off_e=1,splitter=None):
    return read_table_in(fn_input,header,hash,hash_off_s,hash_off_e,splitter)

#Header should be a lis
#THIS FUNCTION SHOULD TAKE CARE OF ALL READ TABLE (hash or list)
#that you ever want, the read_table_no_header,and read_table are just there 
#for reverse compatibliilty
def read_table_in(F,header_list=None,hash=False,hash_off_s=0,hash_off_e=1,splitter=None):

    if not isinstance(F,file):
        F = open(F,'r')

    if header_list==None:
        header = F.readline().rstrip().split(splitter)
    else:
        header = header_list

    if hash:
        print("WARNING hashing values on key %d-%d %s, collisions will not be resolved"%(hash_off_s, hash_off_e,":".join(header[hash_off_s:hash_off_e])), file=stderr)
        table = dict([(":".join(line.split(splitter)[hash_off_s:hash_off_e]),dict([(header[i],item) for i,item in enumerate(line.rstrip().split(splitter))])) for line in F.readlines()])
        F.seek(0)
        if header_list==None: F.readline()
        for line in F.readlines():    
            table[":".join(line.split(splitter)[hash_off_s:hash_off_e])]["_0"] = line
            i+=1
    else:
        table = [dict([(header[i],item) for i,item in enumerate(line.rstrip().split(splitter))]) for line in F.readlines()]
        i=0    
        F.seek(0)
        if header_list==None: F.readline()
        for line in F.readlines():    
            table[i]["_0"] = line
            i+=1
    
    return table


def read_table(input_table,splitter = "\t",alt_header=None):
    F_table = open(input_table,'r')
    
    if alt_header==None:
        header = F_table.readline().rstrip().split(splitter)
    else:
        header = alt_header.rstrip().split(splitter)

    hashed_lines = [dict([(header[i],item.rstrip()) for i,item in enumerate(line.rstrip().split(splitter))]) for line in F_table.readlines() ]

    F_table.close()
    F_table = open(input_table,'r')
    
    if alt_header==None:
        header = F_table.readline().rstrip().split(splitter)
    else:
        header = alt_header.rstrip().split(splitter)

    i=0

    
    for line in F_table.readlines():    
        hashed_lines[i]["_0"] = line
        i+=1

    return hashed_lines,header


#######
#Takes a bed file and turns it into fast accessible numpy tables

class regions_from_bed:
    
    def __init__(self, fn, vals=False,names=False,delim=None):

        self.vals = vals
        self.names = names
        if vals:
            dt = np.dtype([('chr', np.str_, 8), ('start', np.uint32),('end',np.uint32),('val',np.float32)])
        elif names:
            dt = np.dtype([('chr', np.str_, 8), ('start', np.uint32),('end',np.uint32),('name',np.str_,32)])
        else:
            dt = np.dtype([('chr', np.str_, 8), ('start', np.uint32),('end',np.uint32)])
            
        self.regions = np.loadtxt(fn, dtype=dt, delimiter=delim)
        self.regions = np.sort(self.regions,axis=0,order=['chr','start'])
    
    def get_locations_over_interval(self,chr,start,end,truncate_to_query=False):
        #####
        #given a region this function returns all the intervals in the bed file that extend within that region
        #useful for getting all the say genes in a region and their start and end coords
        #truncate to query means 
        #####
        #print self.regions
        chr_spec_regions = self.regions[np.where(self.regions['chr']==chr)]
        region_intervals = chr_spec_regions[np.where(np.logical_and(chr_spec_regions['start']<end,chr_spec_regions['end']>start))]

        if self.vals:
            ret_vals = region_intervals['val']
        elif self.names:
            ret_vals = region_intervals['name']
        else:
            ret_vals = None

        if truncate_to_query:
            for i in range(region_intervals['start'].shape[0]):
                region_intervals['start'][i] = max(start, region_intervals['start'][i])    
                region_intervals['end'][i] = min(end, region_intervals['end'][i])    

        ret_intervals = np.c_[region_intervals['start'],region_intervals['end']]
        return ret_intervals,ret_vals
    
    def region_overlap(self,chr,start,end):
        #######
        ##Takes in  region and returns whether a region in the bed file completely or partially overlaps that interval
        ## used for checking if genes are totally or partially duped. 
        #return values
        # 0: no overlap
        # 1: partial overlap
        # 2: total overlap
        #######
        
        overlapping_regions,ret_vals = self.get_locations_over_interval(chr,start,end)
        #print overlapping_regions
        if overlapping_regions.shape[0] == 0:
            return 0
        elif overlapping_regions.shape[0] > 1:
            return 1
        else:
            if overlapping_regions[0,0] <= start and     overlapping_regions[0,1] >= end:
                return 2
            else:
                return 1

    def frac_represented(self,chr,start,end):
        #get the fraction of the region represented
        #by regions in the bed file
        overlapping_regions,ret_vals = self.get_locations_over_interval(chr,start,end)
        total = np.sum(overlapping_regions[:,1] - overlapping_regions[:,0])
        #return 0 
        return total/float(end-start)
    
    def frac_n_bp_overlap(self,chr,start,end):
        #returns fraction overlapped by items,n_items,nbp overlapped     
        overlapping_regions,ret_vals = self.get_locations_over_interval(chr,start,end)
        n_regions = overlapping_regions.shape[0]
    
        overlapping_regions[:,0] = np.amax(np.c_[overlapping_regions[:,0],np.zeros(n_regions)+start],axis=1)
        overlapping_regions[:,1] = np.amin(np.c_[overlapping_regions[:,1],np.zeros(n_regions)+end],axis=1)
        
        covered_array = np.zeros(end-start)
        for i in range(n_regions):
            covered_array[overlapping_regions[i,0]-start:overlapping_regions[i,1]-start]+=1
        
        covered = np.sum(covered_array>0)
        return covered/float(end-start), n_regions, covered
    
    
class hash_from_file:

  def __init__(self,fn,hash_on=0):

    input = open(fn,'r')
    hash_on = 0
    columns_names = input.readline().split()
    self.hashed_data = {}
    for line in input:
      sline = line.split()
      curr_data  = {}
      for i,column in enumerate(columns_names):
        if i==hash_on: continue
        curr_data[column] = sline[i]
      self.hashed_data[sline[hash_on]] = curr_data


  def __getitem__(self,i):
    return self.hashed_data[i]
    

def get_regions(fn,mask,chunk_iter=3e7,wnd_width=1000,in_chr=None,in_start=None,in_end=None,sunk_based=False,fillout=False):

    ###Fill in end means: if the window is 1000 and the region % 1000 ~= 0, mak the last window 1000+ regions%1000, 
    ###example region si 1200, one window of size 1200

    regions_chrs = []
    regions_coords = []    
    wnds = []
    
    if fn!=None and in_chr==None:
        lines = open(fn,'r').readlines()
        for line in lines:
            if line[0] == "#": continue
            chr,start,end = line.split()
        
            start = int(start)
            end = int(end)
            
            currChromLen = mask.mContigNameLen[chr]

            if sunk_based:
                curr_mask = mask['isntSunkOrIsMasked'][chr][start:end]
            else:
                curr_mask = mask['mask'][chr][start:end,:].sum(1)>0
            #curr_mask = mask['mask'][chr][start:end,:].sum(1)>0

            #curr_mask = mask['mask'][chr][0:currChromLen,:].sum(1)>0
            chr_wnd_bnds,max_wnd = getBoundsForEqualWeightedWindows(curr_mask,int(start)-1,int(end),wnd_width,fillout)
    
            if chr_wnd_bnds != None:
                wnds.append(chr_wnd_bnds)
                regions_chrs.append(chr)
                regions_coords.append(tuple([int(start),int(end)]))
        
        return regions_chrs,regions_coords,wnds

    elif(in_chr!=None and in_start!=None and in_end!=None):
            chr = in_chr
            start = in_start
            end = in_end

            currChromLen = mask.mContigNameLen[chr]
            
            #curr_mask = mask['mask'][chr][start:end,:].sum(1)>0
            if sunk_based:
                curr_mask = mask['isntSunkOrIsMasked'][chr][start:end]
            else:
                curr_mask = mask['mask'][chr][start:end,:].sum(1)>0


            #curr_mask = mask['mask'][chr][0:currChromLen,:].sum(1)>0
            chr_wnd_bnds,max_wnd = getBoundsForEqualWeightedWindows(curr_mask,int(start)-1,int(end),wnd_width,fillout)
    
            if chr_wnd_bnds != None:
                wnds.append(chr_wnd_bnds)
                regions_chrs.append(chr)
                regions_coords.append(tuple([int(start),int(end)]))

    else:
    
        for chr in mask.mContigNameLen:
            if chr.find("random") != -1:
                print("skipping ", chr)
                continue
            
            if in_chr!=None and chr != in_chr: continue

            print("working on ", chr)
            #chunk_iter = 3e7
            currChromLen = mask.mContigNameLen[chr]
            chunkStart = 0
            #curr_mask = mask['mask'][chr][0:currChromLen,:].sum(1)>0
            
            if sunk_based:
                curr_mask = mask['isntSunkOrIsMasked'][chr][0:currChromLen]
            else:
                curr_mask = mask['mask'][chr][0:currChromLen,:].sum(1)>0
        
            chr_wnd_bnds,max_wnd = getBoundsForEqualWeightedWindows(curr_mask,0,currChromLen,wnd_width,fillout)
            print(chr_wnd_bnds)
            print(max_wnd)

            chr_wnd_bnds_by_chunk_idx = np.r_[np.where(np.diff(np.floor(chr_wnd_bnds[:,1]/chunk_iter))==1)[0],chr_wnd_bnds.shape[0]-1]
            

            print(chr_wnd_bnds_by_chunk_idx)

            wnd_bnds_chunks_idx = np.c_[np.r_[0,chr_wnd_bnds_by_chunk_idx[0:chr_wnd_bnds_by_chunk_idx.shape[0]-1]],chr_wnd_bnds_by_chunk_idx[:]]  
            print("chr_wnd_bnds_by_chunk_idx", chr_wnd_bnds[chr_wnd_bnds_by_chunk_idx])
            print("wnd_bnds_chunks_idx", wnd_bnds_chunks_idx)
            for xi in range(wnd_bnds_chunks_idx.shape[0]):
                
                idx_start,idx_end =  wnd_bnds_chunks_idx[xi,:]
                chunk_start = chr_wnd_bnds[wnd_bnds_chunks_idx[xi,0]][0] 
                chunk_end = chr_wnd_bnds[wnd_bnds_chunks_idx[xi,1]][1] 
                
                associated_windows = chr_wnd_bnds[idx_start:idx_end]-chunk_start
                print(associated_windows)

                print(chunk_start)
                print(chunk_end)
                regions_chrs.append(chr)
                regions_coords.append(tuple([chunk_start,chunk_end]))
                wnds.append(associated_windows)
            
        return regions_chrs,regions_coords,wnds    

            #while chunkStart < currChromLen-1:
                #fnout = "%s/%s_%d_%d"%(fn_gbin_out,chr,chunkStart,chunkEnd)
            #    regions_chrs.append(chr)
            #    regions_coords.append(tuple([chunkStart,chunkEnd]))
    
            #    wnd_starts.append(chr_wnd_bnds[curr_chunk_wnd_bnds_idx]-chunkStart)

            #    curr_chunk_wnd_bnds_idx = np.where(chr_wnd_bnds[:,1]<chunk_iter)
            #    chunkEnd = min(chr_wnd_bnds[curr_chunk_wnd_bnds_idx[0,-1],1],currChromLen)

            #    chunkStart+=chunk_iter
            #    chunkEnd = min(chunkStart+chunk_iter,currChromLen)

    return regions_chrs,regions_coords,wnds



def getBoundsForEqualWeightedWindows(   maskVec,
                                        cooStart,
                                        cooStop,
                                        nbpUnmaskedPerWnd,
                                        fillout=False  ):
    maskReg=~(maskVec)
    #maskReg=~(maskVec[cooStart:cooStop])
    csMaskReg=maskReg.astype('int32').cumsum()

    nTotalUnmasked = csMaskReg[-1]
    if nTotalUnmasked==0:
        return None, None

    assert nTotalUnmasked>0, nTotalUnmasked

    igapdWndStarts = np.arange(0,nTotalUnmasked,nbpUnmaskedPerWnd)

    irngungWnd = csMaskReg.searchsorted(igapdWndStarts,'left')

    irngungWnd = np.c_[ irngungWnd[0:irngungWnd.shape[0]-1], irngungWnd[1:] ]

    if fillout:
        if irngungWnd.shape[0] == 0:
            irngungWnd = np.c_[np.array([0]),np.array([maskVec.shape[0]])]
        elif irngungWnd[-1,1]<maskVec.shape[0]:
            irngungWnd[-1,1] = maskVec.shape[0]
        
    #corngWnd = cooStart + irngungWnd
    corngWnd = irngungWnd
    
    if corngWnd.shape[0] == 0:
        ########ADDED Oct 26 2011 for dudes shorter than the length of the window, ie, rep elements
        return np.array([[0, maskVec.shape[0] ]]), maskVec.shape[0]
        #return None, None
    #from IPython.Shell import IPShellEmbed; ipshell = IPShellEmbed([]); ipshell()
    
    
    
    
    return corngWnd,corngWnd[-1][1]




def get_genome_pop_hash(index_files,sex_index):

    YORUBA = ["YORUBA","YRI"]
    HAN_CHINESE = ["HAN"]
    CEPH = ["CEPH","CEU","EUROPEAN"]
    JAPANESE = ["JAPANESE"]
    KOREAN = ["KOREAN"]

    ASIAN =  JAPANESE + KOREAN + HAN_CHINESE
    
    #pops = {"YORUBA":YORUBA,"HAN_CHINESE":HAN_CHINESE,"CEPH":CEPH,"JAPANESE":JAPANESE,"KOREAN":KOREAN}
    pops = {"YORUBA":YORUBA,"ASIAN":ASIAN,"CEPH":CEPH}



    def get_pop(indiv_pop):
        indiv_pop = indiv_pop.split()[0]
        indiv_pop = indiv_pop.split("_")[0]
        indiv_pop = indiv_pop.split("-")[0]
        indiv_pop = indiv_pop.upper()
        for pop,ids in pops.items():
            if(indiv_pop in ids):
                return pop.upper()

        print("not found",indiv_pop)


    genome_pop_hash = {}

    for i in range(len(index_files)):
        kgf.init(index_files[i],0,sex_index[i])
        for genome in kgf.genomes:
            #print genome, kgf.genomes[genome].population
            pop = (get_pop(kgf.genomes[genome].population)).upper()
            genome_pop_hash[genome] = pop
            #print pop
                    
    return genome_pop_hash



def get_adjusted_depth(wssd,grpWssd,GC_content_table,grpGC,chr,start,end,gc_correction_vect,wssd_dim,v=False): 
    #print grpWssd
    #print msk[grpMsk][chr][start:end]
    #print wssd.depth[grpWssd][chr][start:end,:,0].sum(1).mean(0)
    #print GC_content_table[grpGC][chr][start:end]

    #print chr
    if(v):print("getting gc correction...")
    adjustment = gc_correction_vect.take(GC_content_table[grpGC][chr][start:end])    
    if(v):print("getting depths...")
    depths = np.nan_to_num(wssd.depth[grpWssd][chr][start:end,:,wssd_dim]).astype(np.float64).sum(1)
    if(v):print("performing adjustment...")
    adjusted_depth = depths*adjustment
    #adjusted_depth[np.where(adjusted_depth<0)]=0

    if(v):print("original depth track len %d"%(adjusted_depth.shape[0]))
    return adjusted_depth



def arrayContigRangeAndVal( A ):
    dA=np.diff(A)
    N=A.shape[0]
    nzdA,=np.nonzero(dA)
    print(np.nonzero(dA))
    print(nzdA)
    liRangesAll=np.c_[ np.r_[0,nzdA+1],np.r_[nzdA,N-1],A[np.r_[0,nzdA+1]] ]
    return liRangesAll


###########
##This is a confusing name, what it means is
##get the entire contig (in this case it's usually for bacs) with spliced out masked regions
##BUT -> adjusted for read GC 
###########
def get_total_fuguized_adjusted_depth(wssd,grpWssd,GC_content_table,grpGC,contig_name,gc_correction_vect,mask_vect,chr,start,end,wssd_dim):

    #print contig_name
    #print wssd.depth[grpWssd][contig_name].shape[0]
    non_fugu_depth = get_adjusted_depth(wssd,grpWssd,GC_content_table,grpGC,chr,start,end,gc_correction_vect,wssd_dim)
    fuguized_depth = non_fugu_depth[np.where(mask_vect==False)]
    return fuguized_depth    

def make_adjusted_depth_binary(indiv,wssd_lib_dir,bac_analysis_lib_dir,chr,start,end,fn_gc_table,fn_contig,mask,fnout):

    adjusted_depth = get_region_from_indiv(indiv,wssd_lib_dir,bac_analysis_lib_dir,chr,start,end,fn_gc_table,fn_contig)
    mask_chunk = np.cast['i4'](mask["mask"][chr][np.arange(start,end),:].sum(1)>0)

    adjusted_depth = np.cast['i4'](np.around(adjusted_depth))
    output_matrix = np.c_[adjusted_depth,mask_chunk]
    print(output_matrix)
    print(output_matrix.dtype)
    output_matrix.tofile(fnout)        

def load_gc_correction_vect(lib_dir,wssd_dim):

    fn = "%s/cn2-gc-depth-dim%d.txt"%(lib_dir,wssd_dim)
    f = open(fn,'r')
    len_line = f.readline()
    len = int(len_line.rstrip("\n"))
    line = f.readline()
    gc_correct_vect = []
    n_corrections = 0
    while(line):
        n_corrections+=1
        gc_correct_vect.append(float(line.split()[1]))
        line = f.readline()

    assert(n_corrections == 2*len+2)
    gc_correct_vect = np.array(gc_correct_vect)
    return gc_correct_vect

def get_region_from_indiv(indiv,wssd_lib_dir,bac_analysis_lib_dir,chr,start,end,fn_gc_table,fn_contig):

    #mask_table = DenseTrackSet(fn_contig,
#                                                        fn_mask_table,
#                                                        overwrite=False,
#                                                        openMode='r') 
    gc_content_table =   DenseTrackSet(fn_contig,
                                                        fn_gc_table,
                                                        overwrite=False,
                                                        openMode='r') 

    read_depth = np.zeros(end-start)
    print(read_depth)
    print(end-start)

    libls = os.listdir(bac_analysis_lib_dir)
    for lib in libls:
        if(lib[0:1]=="_"):
            continue
        lib_qc_dir = "%s/%s"%(bac_analysis_lib_dir,lib)
        if(kgf.lib_pass_qc(lib_qc_dir)):
            print("library %s passed qc: loading..."%(lib))
            fn_wssd = "%s/%s/hg18rmsk.wssd"%(wssd_lib_dir,lib)
            gc_correction_vect = load_gc_correction_vect(lib_qc_dir)

            print("opening wssd %s"%(fn_wssd))

            wssd = WssdFile(fn_contig,
                                            fn_wssd,
                                            overwrite=False,
                                            openMode='r') 

            print("wssd opened..")
            grpWssd = "wssd"
            grpGC = "GC_content"
            read_depth+=get_adjusted_depth(wssd,grpWssd,gc_content_table,grpGC,chr,start,end,gc_correction_vect,0)
        
            wssd.tbl.close()
            del wssd

    gc_content_table.tbl.close()
    del gc_content_table

    return read_depth

def start_new_track_in_bed(fn_bed_out,name,viewLimits=[0,26],yLine='on'):
    f = open(fn_bed_out,'a')
    f.write("""track type=bedGraph name="%s" visibility=full autoScale=off color=0,200,0 altColor=0,0,0 viewLimits=%d:%d yLineMark=2 yLineOnOff=%s\n"""%(name,viewLimits[0],viewLimits[1],yLine))
    f.close()

def init_bed(fn_bed_out):
    f = open(fn_bed_out,'w')
    f.close()

def output_wiggle(outStarts,values,chr,name,fnOut,step):

    print("outputting %s..."%(fnOut))
    fout = open(fnOut,"w")
    fout.write("track type=wiggle_0 name=%s_depths\nvariableStep span=%d chrom=%s\n"%(name,step,chr))

    for i in range(len(outStarts)):
        fout.write("%d %d\n"%(outStarts[i],values[i]))
    #output = np.c_(outStarts,values)
    print("done")


class contig:
  def __init__(self,key,name,length):
    self.key,self.length,self.name = (key,name,length)

def get_contig_objects(fn_contig):

    contigs_bykey = {}
    contigs_byname = {}
    contigs_list = open(fn_contig,"r").readlines()
    for contig_line in contigs_list:
        key,len = contig_line.rstrip("\n").split("\t")
        name = key
        new_contig = contig(key,int(len),name)
        contigs_bykey[key] = new_contig
        contigs_byname[name] = new_contig
  
    return (contigs_byname,contigs_bykey)





def calkan_window_avg_masked(depth,mask,start,end,width,slide,type='mean'): 

    print(mask.shape[0])
    print(depth.shape[0])
    assert(mask.shape[0]==depth.shape[0])

    cat_unmasked = depth[np.where(mask==0)]    
    kernel = np.ones(width,np.uint32)    
    convolution = np.convolve(cat_unmasked.astype(np.uint16),kernel,'full')

    n_windows = ((cat_unmasked.shape[0]/width)-1)*(width/slide)+1
    print(cat_unmasked.shape[0])
    print(width)
    print(slide)
    print("**")


    if(type=='mean'):
        convolution = convolution.astype(np.float64)/width

    convStarts = np.arange(width-1,(width-1)+n_windows*slide,slide)
    catStarts = np.arange(0,n_windows*slide,slide)
    catStartsArray = np.zeros(cat_unmasked.shape[0],np.bool)
    catStartsArray[catStarts] = 1
    
    #catStartsArray is an array of length cut_unmasked with a 1 every "width" position, 
    #ie, indicating where a Start should be
    
    depth = convolution[convStarts]
    print(n_windows, start, end)
    print(convolution[convStarts].shape[0])
    assert(n_windows==convolution[convStarts].shape[0])

    startsArray = np.zeros(mask.shape[0],np.bool)
    
    startsArray[np.where(mask==0)] = catStartsArray
    outStarts,= np.where(startsArray==1)
    outStarts += start
    outEnds = np.hstack((outStarts[1:outStarts.shape[0]],outStarts[-1]+slide))

    return(depth,outStarts,outEnds)

# for an array, return inclusive ranges which are continuously
# equal to a single value (excluding ranges equal to zeroDefVal)    
##########WARNING, this reports in bed form with -1 from the end coord
##########Need to add 1 to the second coordinatez!!!!!
#array([1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1])
#>>> arrayToNonZeroContigRanges(k)
#array([[ 0,  0],
#       [ 4,  6],
#       [ 8,  8],
#       [12, 12]])
#######WARNING, the last value is inclusive! It should have a +1!
def arrayToNonZeroContigRanges( A, zeroDefVal=0):
    ####April 12 2012 I added a +1 to this: np.r_[nzdA,N-1]+1] so it was zero based open coords
    dA=np.diff(A)
    N=A.shape[0]
    nzdA,=np.nonzero(dA)
    liRangesAll=np.c_[ np.r_[0,nzdA+1],
    np.r_[nzdA,N-1] +1]
    iiR = A.take(liRangesAll[:,0]) != zeroDefVal
    return liRangesAll[iiR,:]

#output bed file in form 1-N inclusive of end_coord
#ie give genomic coordinates
def output_bed_of_mask(mask,chr,start_coord,end_coord,fnout=None):

    #outf = open(fnout,"w")
    #print mask
    #mask = mask['mask'][chr][start_coord-1:end_coord,:].sum(1)>0
    mask = mask['mask'][chr][start_coord:end_coord,:].sum(1)>0
    #print mask
    ranges = arrayToNonZeroContigRanges(mask) 
    ranges = ranges+start_coord
    N = ranges.shape[0] 
    
    for i in range(N):
        #outf.write("%s %d %d\n"%(chr,ranges[i][0],ranges[i][1]))
        print("%s\t%d\t%d"%(chr,ranges[i][0],ranges[i][1]))
    #outf.close()

def make_bed_start_end_val(chr,starts,ends,val,fnout,name,append=False):

    if(not(append)):
        outf = open(fnout,"w")
        outf.write("""track type=bedGraph name="%s"\n"""%(name))
    else:
        outf = open(fnout,"a")

    N = starts.shape[0]
    #print starts.shape[0]
    #print ends.shape[0]
    #print val.shape[0]
    for i in range(N):
        outf.write("%s %d %d %f\n"%(chr,starts[i],ends[i],val[i]))

def make_big_wig(fn_temp_bed,fnOutBigWig,fnOutTrackDef,fn_contig_file):
    
    #fnOutBigWig = "%s/%s_%s.bw"%(browser_track_dir,genome,o.name)
    #fnOutTrackDef = "%s/%s_%s.trackdef"%(browser_track_dir,genome,o.name)

    print(fnOutBigWig)
    print(fnOutTrackDef)

    command = "~/local_installations/ucscOLD/ucsc/bin/bedGraphToBigWig %s %s %s"%(fn_temp_bed,fn_contig_file,fnOutBigWig)
    os.system(command)
    print(command)

    track_def = """track type=bigWig name="%s_%s_%s" description="%s_%s" color=%s visibility=full autoScale=off alwaysZero=on maxHeightPixels=50:50:20 viewLimits=0.0:10 windowingFunction=mean smoothingWindow=3 yLineMark=2.0  yLineOnOff=on dataUrl=%s\n"""%(o.name,genome,sex,o.name,genome,color_str,fnOutBigWig)
    open(fnOutTrackDef,'w').write(track_def)
    print(fnOutTempBedgr)
    
    if o.out_raw_bed:
        command = "mv %s %s/%s_%s.bed"%(fnOutTempBedgr, browser_track_BED_dir,genome,o.name)
        os.system(command)
    else:
        os.unlink(fnOutTempBedgr)
    os.unlink(fnOutTempBedgrHM)

