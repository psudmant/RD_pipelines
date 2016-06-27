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

import scipy as scp
import numpy as np
import tables

from fitting import *

#from math import *
import Bio.Statistics.lowess as biostats

from wssd_common import *
from wssd_pw_common import *

from scipy.stats import norm


#def get_loc(msk,wssd,start,end,chr):

def mkdir(dir,file):

    ls_dir = os.listdir(dir)
    if(not(file_exists(ls_dir,file))):
        command = "mkdir %s/%s"%(dir,file)
        print command
        os.system(command)
    return "%s/%s"%(dir,file)

def file_exists(ls,file):
  for f in ls:
    if(f==file):
            return 1
  return 0

def do_lsqfit_w(x,y,w,pinit,fitfunc):

    werrfunc = lambda p, x, y, w: (y - fitfunc(p, x)) * w
    out = scp.optimize.leastsq(werrfunc, pinit,args=(x,y,w), full_output=1)
    return out[0]

def do_lsqfit(x,y,pinit,fitfunc):

    werrfunc = lambda p, x, y: (y - fitfunc(p, x))
    out = scp.optimize.leastsq(werrfunc, pinit,args=(x,y), full_output=1)
    return out[0]

def fit_line(mus,cps,force_zero=0):

    if(force_zero):
        fitfunc = lambda p, x: p[1] * x
    else:
        fitfunc = lambda p, x: p[0] + p[1] * x


    if(mus.shape[0] == 1):
        fit_params_u = [0,mus[0]/cps[0]]
        return fit_params_u    

    #print cps, mus, sigs
    pinit = [1.0, -1.0]
 
    fit_params_u = do_lsqfit(cps,mus,pinit,fitfunc)

    #full_range_cns = np.arange(card)
    #fit_mus = fitfunc(fit_params_u,full_range_cns)
    #fit_sigs = fitfunc(fit_params_sig,full_range_cns)
    #fit_sigs = np.abs(fit_sigs)
    if(force_zero):
        fit_params_u[0] = 0
    return (fit_params_u)


def sum_g(dists,i,group_size,n_single):
    sum = np.zeros(dists[i].shape[0])
    i=i*group_size + n_single
    print "\tmade from:"
    for k in range(group_size):
        print "\t",k+i    
        sum+= dists[k+i]

    sum /= group_size
    return sum

class sexIndex:
    def __init__(self,sex_index_file):
        self.sex_id = {}
        for sex_line in open(sex_index_file,'r').readlines():
            self.sex_id[sex_line.split()[0]] = sex_line.split()[1]

    def __getitem__(self,key):
        return self.sex_id[key]    
    
class region:
    
    def __init__(self,chr,start,end,cp,name=None):
        self.chr =chr;
        self.start = int(start);
        self.end = int(end);
        self.cp = cp
        #print "\t", cp
        if(name == None):
            self.name = "%s-%d-%d"%(self.chr,self.start,self.end)    
        else:
            self.name = name
    
    def get_cp(self,indiv):
        global genome_info
    
        #genome_info = kgf.genome_info(o.fn_in_genomes,o.fn_sex_pop_index)

        #q_indiv = indiv
        #if indiv.find('sampled') != -1:
        #    q_indiv = indiv.split("-")[0]

        if(self.chr == "chrX"):
            if(genome_info.get_sex(indiv) == "M"):
                return 1
            else:
                return 2
        else:
            return int(self.cp)

def import_bac_regions(fn_in_bacs,regions,only_use_cps):

    

    lines = open(fn_in_bacs,'r').readlines()
    for line in lines:
        if(line[0] == "#"): continue
        (name,chr,start,end) = line.split()[0:4]
        cp=int(name.split(",")[1])
    
        name = "%d_%s_%s_%s_%s"%(cp, name, chr,start,end)
        print name
        newregion = region(chr,start,end,cp,name)
        if only_use_cps == None or cp in only_use_cps:
            regions.append(newregion)
            print "using region:", cp, chr, start, end

def import_tomas_regions(fn_in_tomas,regions):
    lines = open(fn_in_tomas,'r').readlines()
    for line in lines:
        if(line[0] == "#"): continue
        (chr,start,end) = line.split()[0:3]
    
        if(line.split()[3]=="*"):
            continue
        elif(line.split()[3]=="#" or line.split()[3]=="%"):
            name = "%s_%s_%s_%s"%(line.split()[4], chr,start,end)
            cp = line.split()[4]
        else:
            name = "%s_%s_%s_%s"%(line.split()[7], chr,start,end)
            cp = line.split()[7]

        newregion = region(chr,start,end,cp,name)
        regions.append(newregion)

def output_regions_bed(regions):
    for region in regions:
        print region.chr,region.start,region.end


def add_depth_to_fig(fignum,depth,sp,name,n):

    #subplot(n,1,sp)
    #plot(depth,"-",linewidth=.25)

    plot(depth,'-',linewidth=.25)
    legend([name],loc='lower right')
    xlabel(name)

def analyze_regions(regions,input_wssd,mask_wssd,output_file,group):

    n = len(regions)
    #subplot(n,1,1)
    #subplots_adjust(hspace=1)

    k=0
    for region in regions:
        f=figure(1)
        f.set_figwidth(20)
        f.set_figheight(5)
        #f.set_figheight(100)
        depth = np.nan_to_num(input_wssd.depth[group][region.chr][region.start:region.end,:,0]).astype(np.float64).sum(1)
        masked = mask_wssd["mask"][region.chr][region.start:region.end,:].sum(1)>0
        depth[np.where(masked>0)] = 0
        add_depth_to_fig(1,depth,k,region.name,n)    
  
        figure_name = "%s%d.control_region_plot.jpg"%(output_file,k)
        ylabel("depth")
        print figure_name
        savefig(figure_name,format='png')
        close(1)
        
        k+=1

def calkan_smooth(depth,window_size,slide_by,masked):
    window = np.ones(window_size)
    masked = (masked>0)

    
    #mask_conv = np.convolve(masked,window,'some')

    collapsed_depths = depth[np.where(masked==0)]
    conv = np.convolve(collapsed_depths,window,'some')
    new_depths = np.zeros(depth.shape[0])
    
    #new_depths[np.where(masked==0)] = conv        
    new_depths = conv
    new_depths = new_depths/1000

    new_depths = new_depths[np.arange(0,new_depths.shape[0]-1,slide_by)]
    return new_depths
    
def simple_fit_params(regions,region_depths,cat_regions_by_cp,input_wssd,mask_wssd,window_size,slide_by,out_file_name,indiv,dim):

    depth_cp = {}

    fn_out_file_per_bac = "%s.per_bac_params"%(out_file_name)
    fout_per_bac = open(fn_out_file_per_bac,'w')
    
    fn_out_file_per_1kb = "%s.per_1kb_params"%(out_file_name)
    fout_per_1kb = open(fn_out_file_per_1kb,'w')
    
    #fn_out_file_per_1bp = "%s.per_1bp_params"%(out_file_name)
    #fout_per_1bp = open(fn_out_file_per_1bp,'w')

    for region in regions:
        p = "%f %f"%(region_depths[region.name].mean(),sqrt(region_depths[region.name].var()))
        reg_len = region_depths[region.name].shape[0]
        fout_per_bac.write("%s %d %s %d\n"%(region.name,region.get_cp(indiv),p,reg_len))    


        start = 0
        while(start<reg_len):
            end = min(reg_len,start+1000)
            p = "%f %f"%(region_depths[region.name][start:end].mean(),sqrt(region_depths[region.name][start:end].var()))
            #p = str(region_depths[region.name][start:end].mean())
            fout_per_1kb.write("%s.%d %d %s %d\n"%(region.name,start,region.get_cp(indiv),p,end-start))    
            start+=1000
        
        #start = 0
        #while(start<reg_len):
        #    end = min(reg_len,start+1)
        #    p = "%f %f"%(region_depths[region.name][start:end].mean(),sqrt(region_depths[region.name][start:end].var()))
        #    fout_per_1bp.write("%s.%d %d %s %d\n"%(region.name,start,region.get_cp(indiv),p,end-start))    
        #    start+=1
            
    fout_per_1kb.close()
    #fout_per_1bp.close()
    fout_per_bac.close()

    params = {}
    fn_out_file = "%s.cat_bac_params"%(out_file_name)
    fout = open(fn_out_file,'w')
    
    for cp,depth in cat_regions_by_cp.iteritems():
        print "fitting simple"
        p = str(depth.mean())
        p = "%f %f"%(depth.mean(),sqrt(depth.var()))

        params[cp] = p
        print p
        reg_len = depth.shape[0]    
        fout.write("%d %s %d\n"%(cp,p,reg_len))    

    fout.close()

def fit_params(regions,region_depths,cat_regions_by_cp,input_wssd,mask_wssd,window_size,slide_by,out_file_name,indiv,dim,type):

    depth_cp = {}

    fn_out_file_per_bac = "%s.per_bac_params"%(out_file_name)
    fout_per_bac = open(fn_out_file_per_bac,'w')

    #fn_out_file_per_1kb = "%s.per_bac1kb_params"%(out_file_name)
    #fout_per_bac = open(fn_out_file,'w')

    for region in regions:
        p = fit_dist_to_data(region_depths[region.name],type)
        p=[str(ap) for ap in p]
        reg_len = region_depths[region.name].shape[0]
        fout_per_bac.write("%s %d %s %d\n"%(region.name,region.get_cp(indiv)," ".join(p),reg_len))    
    
        #if(window_size>1):    
        #    depth = calkan_smooth(depth,window_size,slide_by,masked)    

    fout_per_bac.close()

    #hists = {}
    #f=figure(1)
    #f.set_figwidth(40)
    #f.set_figheight(40)
    #leg = []

    params = {}
    
    #print "create_histogram - fit and plot"

    fn_out_file = "%s.cat_bac_params"%(out_file_name)
    fout = open(fn_out_file,'w')
    
    for cp,depth in cat_regions_by_cp.iteritems():
        #hists[cp] = np.histogram(depth,bins=3000,range=[0,3000])[0].astype(float64)    
        #print cp , hists[cp]
        #hists[cp] = hists[cp]/hists[cp].sum()
        #plot(hists[cp][0:150])
        #leg.append(str(cp))
        print "fitting %s"%(type)
        #depth = depth[0:depth.shape[0]:10]
        p = fit_dist_to_data(depth,type)
        p=[str(ap) for ap in p]

        params[cp] = p
        print p
        reg_len = depth.shape[0]    
        fout.write("%d %s %d\n"%(cp," ".join(p),reg_len))    

    fout.close()
        #fit = fit_dist(np.arange(0,3000),p,type)
        #plot(fit[0:150],'--')        
        #leg.append(str(cp)+"-fit")
    
    #legend(leg)
    #print out_file_name
    #savefig("%s.png"%(out_file_name),format='png')
    #close(1)
    #sys.exit(1)

    #return params

def copy(bac,indiv):
    if(kgf.genomes[indiv].sex == "F" and bac.is_x):
        return 2*bac.cp
    return bac.cp

def fit_dist(x,p,type='gaussian'):
    if(type=='gaussian' or type == 'gaussian_lsq'):
        return norm.pdf(x,p[0],p[1])
    if(type == 'cont_poisson'):
        return (exp(-p)*p**x)/scp.special.gamma(x+1)    
    if(type == 'gamma'):
        return scp.stats.gamma.pdf(x,p[0],p[1])
    if(type=="simple"):
        return np.zeros(x.shape[0])

def fit_dist_to_data(input_data,type='lsq'):

    new_fit = fit(input_data)
    #return fit_by_em(input_data)
    if(type=='simple'):
        return new_fit.fit_by_simple()
        #return(input_data.mean(),input_data.var())
    if(type=='gaussian'):
        return new_fit.fit_by_fmin()
    if(type=='cont_poisson'):
        return new_fit.fit_by_cont_poisson()
    if(type=='gamma'):
        return new_fit.fit_gamma()
    if(type=='gaussian_lsq'):
        return new_fit.fit_by_gaussian_lsq()    
    if(type=='weibull_lsq'):
        return new_fit.fit_by_weibull_lsq()
    else:
        print "ERROR, no type %s"%(type)
        sys.exit(1)

def graph_line(ys,xs,ps,outfile):
    fitfunc = lambda p, x: p[0] + p[1] * x
    
    figure(1)
    plot(xs,ys,".")
    xrange = np.array((0,np.amax(xs)))
    plot(xrange,fitfunc(ps,xrange),'-')
    savefig("%s.png"%(outfile),format='png')
    close(1)

def output_parameter_file(ps,cps,bac_loc,indiv,type,dim,additional_info = ''):

    fn_param_file = "%s/%s/_filtered_libs_analysis/fit_params%s_%s_dim%d"%(bac_loc,indiv,additional_info,type,dim)    
    fout = open(fn_param_file,'w')
    for i in range(len(cps)):
        strout = "%d"%(cps[i])
        for p in ps[i]:
            strout += " %f"%(p)
            
        fout.write("%s 1\n"%(strout))


def load_region_depths(regions,wssd_input,mask_wssd,dim,input_group,indiv,edit):
    
    region_depths = {}
    cat_regions_by_cp = {}

    print "getting regions %d..."%(dim)

    for region in regions:
        print "getting %s"%(region.name)


        if(edit == -1):
            depth = np.nan_to_num(wssd_input.depth[input_group][region.chr][region.start:region.end,:,dim]).astype(np.float64).sum(1)
        else:
            depth = np.nan_to_num(wssd_input.depth[input_group][region.chr][region.start:region.end,edit,dim]).astype(np.float64)
        
        masked = mask_wssd["mask"][region.chr][region.start:region.end,:].sum(1)>0
        depth = depth[np.where(masked==0)] #JUST THE DEPTHS FOR THE HIST 
    
        if(region.name in region_depths):
            print "ERROR - duplicate region name %s"%(region.name)
            sys.exit(1)

        region_depths[region.name] = depth

        if(region.get_cp(indiv) in cat_regions_by_cp):
            cat_regions_by_cp[region.get_cp(indiv)] = np.r_[cat_regions_by_cp[region.get_cp(indiv)],depth]
        else:    
            cat_regions_by_cp[region.get_cp(indiv)] = depth

    print "done"
    return region_depths,cat_regions_by_cp

def fit_params_on_input_genomes(fn_input_genomes,regions,mask_wssd,fn_contigs,wssd_file_name,bac_type,output_dir_name,types,edits):
          
    input_genomes = open(fn_input_genomes,'r').readlines()

    #edits = [-1,0] #-1 stands for :, 0, 1 2 are edit 0, 1, 2
    edit_str = {-1:"all",0:"0",1:"1",2:"2"}

    for input_genome in input_genomes:
        (genome,wssd_libs_loc,bac_loc,chunked_reads_loc,primary_analysis_loc) = input_genome.split()

        fn_in_wssd = "%s/%s/combined_corrected_wssd/%s"%(primary_analysis_loc,genome,wssd_file_name)
        indiv = genome.split(".")[0]
        
        ####output_dir
        indiv_bac_dir = "%s/%s"%(bac_loc,indiv)
        
        #fn_output_dir = "%s/%s/_filtered_libs_analysis"%(bac_loc,indiv)
        fn_output_dir = mkdir(indiv_bac_dir,output_dir_name) 
    
        print fn_in_wssd
        if(not(os.path.exists(fn_output_dir)) or not(os.path.exists(fn_in_wssd))): continue

        input_group = "wssd.combined_corrected" 

        print fn_in_wssd
        adjusted_wssd = WssdFile(fn_contigs,
                                                fnWssd=fn_in_wssd,
                                                overwrite=False,
                                                openMode='r')

        #mkdir("./",genome)
        for dim in [0,1]:

            for edit in edits:
        
                region_depths,cat_regions_by_cp = load_region_depths(regions,adjusted_wssd,mask_wssd,dim,input_group,indiv,edit)
                
                                        
                
                simple_fit_params(regions,region_depths,cat_regions_by_cp,adjusted_wssd,mask_wssd,1,1,"%s/fit_params-%s-%s-dim%d-ed-%s"%(fn_output_dir,bac_type,'simple',dim,edit_str[edit]),indiv,dim)
                
                if(types[0]=="simple"): continue
                for type in types:        
                    #analyze_regions(regions,adjusted_wssd,mask_wssd,fn_output_dir,"wssd.combined_corrected")
                    fit_params(regions,region_depths,cat_regions_by_cp,adjusted_wssd,mask_wssd,1,1,"%s/fit_params-%s-%s-dim%d-ed-%s"%(fn_output_dir,bac_type,type,dim,edit_str[edit]),indiv,dim,type)
        
        adjusted_wssd.tbl.close()

def choose_inputs(good_cps,ps_hash):
    ps = []
    cps = []    
    for cp in good_cps:
        if(cp in ps_hash):
            cps.append(cp)
            ps.append(ps_hash[cp])        
    
    cps = np.array(cps)
    ps = np.array(ps)

    return cps,ps

def analyze_yoruba_trio_uncorrected(combined_libs_dir,fn_input_genomes,regions,mask_wssd,fn_contigs,window,slide):
    
    input_genomes = open(fn_input_genomes,'r').readlines()

    for input_genome in input_genomes:
        (genome,wssd_libs_loc,bac_loc,chunked_reads_loc) = input_genome.split()


        if(genome[0:6] =="NA1850"):
            genome = genome.split(".")[0]
            fn_in_wssd = "%s/%s/combined/hg18rmsk.wssd"%(wssd_libs_loc,genome)
            print fn_in_wssd
            adjusted_wssd = WssdFile(fn_contigs,
                                                fnWssd=fn_in_wssd,
                                                overwrite=False,
                                                openMode='r')

            genome = "%s-unadjusted"%(genome)
            mkdir("./",genome)
            fn_output_dir = "./%s/"%(genome)
            #analyze_regions(regions,adjusted_wssd,mask_wssd,fn_output_dir,"wssd.combined_corrected")
            create_histograms(regions,adjusted_wssd,mask_wssd,window,slide,"%s/histogram1_1"%(fn_output_dir),"wssd")
            adjusted_wssd.tbl.close()


    #ylim((0,2*depth.mean()))
if __name__=='__main__':

    #this program takes in "control region" locations and then creates 
    #parameter fits for different parameters

    opts = OptionParser()
    #opts.add_option('','--in_tomas_regions',dest='fn_in_tomas',default=None)
    opts.add_option('','--in_bac_regions',dest='fn_in_bacs')
    opts.add_option('','--in_mask',dest='fn_in_mask')
    opts.add_option('','--in_contigs',dest='fn_contigs')
    opts.add_option('','--input_genomes',dest='fn_input_genomes')
    #opts.add_option('','--primary_analysis_dir',dest='fn_primary_analysis')
    #opts.add_option('','--sex_file',dest='fn_sex_index')
    opts.add_option('','--sex_pop_index',dest='fn_sex_pop_index')
    opts.add_option('','--bac_type',dest='bac_type')
    opts.add_option('','--wssd_file_name',dest='wssd_file_name',default="wssd.combined_corrected")
    opts.add_option('','--output_directory',dest='output_directory',default="_filtered_libs_analysis")
    opts.add_option('','--type',dest='type',default=None)
    opts.add_option('','--edits',dest='edits',default="-1")
    opts.add_option('','--only_use_cps',dest='only_use_cps',default=None)
    
    #opts.add_option('','--out_location',dest='fn_out_location')

    (o, args) = opts.parse_args()
    #sexIdx = sexIndex(o.fn_sex_index)
    genome_info = kgf.genome_info(o.fn_input_genomes,o.fn_sex_pop_index)

    regions = []

    #if(o.fn_in_tomas != None):
    #    print "importing PTR-HG18 overlap (Tomas regions)"
    #    import_tomas_regions(o.fn_in_tomas,regions)


    edits = [int(e) for e in o.edits.split(":")]
    
    print "edits: ",edits

    if o.only_use_cps != None:
        only_use_cps = [int(e) for e in o.only_use_cps.split(":")]
    else:
        only_use_cps = None
    
    import_bac_regions(o.fn_in_bacs,regions,only_use_cps)

    print o.bac_type

    types = ['gaussian_lsq','cont_poisson','gaussian']

    if(o.type!=None):
        if o.type in types or o.type == 'simple':
            types = [o.type]
        else:
            print o.type    
            sys.exit(1)

    mask_wssd = DenseTrackSet(o.fn_contigs,
                                                        fnWssd=o.fn_in_mask,
                           overwrite=False,
                            openMode='r')

    #output_regions_bed(regions)

    #analyze_regions(regions,input_wssd,mask_wssd,o.fn_out_file)
    #create_histograms(regions,input_wssd,mask_wssd,1,1,"hg18_self1_1","wssd")
    #analyze_yoruba_trio_uncorrected(o.fn_primary_analysis,o.fn_input_genomes,regions,mask_wssd,o.fn_contigs,1,1)
    fit_params_on_input_genomes(o.fn_input_genomes,regions,mask_wssd,o.fn_contigs,o.wssd_file_name,o.bac_type,o.output_directory,types,edits)

