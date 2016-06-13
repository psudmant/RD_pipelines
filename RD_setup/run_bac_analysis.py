import sys,os
import numpy as np
import scipy as sp
#import scipy.interpolate.fitpack2 as sp_i_fp2
import kg_file_handling as kgf
#from kg_file_handling import *
######import pygr.Data
from pygr import worldbase
from pygr import seqdb
#from scipy import *

#from pylab import *
#import matplotlib
import scipy.optimize
import scipy.stats
from optparse import OptionParser

from wssd_common import *
from wssd_pw_common import *

import Bio.Statistics.lowess as biostats

class bac():

	def get_GC(self):
		global GC_width

		seqstr = str(self.sequence).upper()

		char_seqstr = np.array(seqstr,'c')
		gORc_array = (char_seqstr=='C')|(char_seqstr=='G')
		GC_conv = np.ones(2*GC_width+1)
		conv = np.convolve(gORc_array,GC_conv,'same')
		self.GC_count = conv
		self.GC_percent = conv/(2*GC_width+1)

	def update_cp(self,cp):
		self.cp = int(cp)

	def get_unmasked_depth(self,depth):
		return depth[np.where(self.mask == 0)]

	def get_unmasked_median(self,depth):
		return np.median(self.get_unmasked_depth(depth))

	def get_unmasked_mean(self,depth):
		return self.get_unmasked_depth(depth).mean()

	def make_mask(self,mask_table,chr,start,end):
		self.mask = mask_table["mask"][chr][start:end].sum(1).astype(np.bool)
		self.unmasked_len = np.where(self.mask==False)[0].shape[0]

	def set_is_x(self,is_x):
		self.is_x = is_x

	def __init__(self,name,seq,chr,start,end,mask_table,gc_content_table,source,cpn=None,is_x = False):
		self.is_x = is_x
		self.source = source
		self.mask_table = mask_table
		self.gc_content_table = gc_content_table
		self.chr = chr
		self.start = start
		self.end = end
		self.sequence = seq
		self.name = name
		self.get_GC()
		self.make_mask(mask_table,chr,start,end)
		if(cpn!=None):
			self.cp = int(cpn)
		#self.gc_content_vect = gc_content_vect

def file_exists(ls,file):
	for f in ls:
		if(f==file):
			return 1
	return 0


def load_control_bac_objects(fn_bac_list,mask_table,gc_content_table,genome):

	control_bacs = {}
	total_unmasked_len = 0
	if genome.upper() == "HG18":
		#HG_sequence = pygr.Data.Bio.Seq.Genome.Human.hg18()
		#HG_sequence = worldbase.Bio.Seq.Genome.HUMAN.hg18() #chopped out march 9th
		HG_sequence = seqdb.SequenceFileDB("/net/eichler/vol7/home/psudmant/genomes/fastas/hg18/hg18.fa") #added march 9th
	elif genome.upper() == "HG19":
		#HG_sequence = pygr.Data.Bio.Seq.Genome.Human.hg19()
		#HG_sequence = worldbase.Bio.Seq.Genome.HUMAN.hg19() #commented march 9
		HG_sequence = seqdb.SequenceFileDB("/net/eichler/vol7/home/psudmant/genomes/fastas/hg19/hg19.fa") #added march 9th
	else:
		# If the given genome isn't in the predefined list, assume it is actually the FASTA file for the genome itself.
		HG_sequence = seqdb.SequenceFileDB(genome)

	for line in open(fn_bac_list,'r').readlines():
		if(line[0]=="#"): continue
		#print "loading %s"%(line.rstrip())
		(name,chr,start,end) = line.split()[0:4]

		start = int(start)
		end = int(end)
		#seq = HG_sequence[chr][start:end]
		seq = "ATAT"
		#print seq
		id = "%s%s:%i-%i"%(name,chr,start,end)
		if(chr.upper()=="CHRX"):
			cpn=1
		else:
			cpn = int(name.split(",")[1])

		#gc_content_vect = gc_content_table["GC_content"][chr][start:end]
		new_bac = bac(id,seq,chr,start,end,mask_table,gc_content_table,genome,cpn)
		if(chr=="chrx" or chr=="chrX"):
			new_bac.set_is_x(True)
		control_bacs[id]=new_bac

		total_unmasked_len+=new_bac.unmasked_len
		print new_bac.unmasked_len , "\t" , len(new_bac.sequence) , "\t", new_bac.cp

	return (total_unmasked_len, control_bacs)


#######
#DEPRICATED

#def load_additional_control_regions(control_bacs,fn_hg18_ptr_inv,mask_table,gc_content_table,genome):
#
#	total_unmasked_len = 0
#	if genome.upper() == "HG18"
#		HG_sequence = pygr.Data.Bio.Seq.Genome.Human.hg18()
#	if genome.upper() == "HG18"
#		HG_sequence = pygr.Data.Bio.Seq.Genome.Human.hg19()
#
#	for line in open(fn_hg18_ptr_inv,'r').readlines():
#		if(line[0]=="#"): continue
#		#print "loading %s"%(line.rstrip())
#		(chr,start,end,cpstate) = line.split()[0:4]
#		if(cpstate == "*"):
#			continue
#		elif(cpstate == "**"):
#			(chr,start,end,cpn) = line.split()[4:8]
#		elif(cpstate == "#" or cpstate=="%"):
#			cpn = line.split()[4]
#		else:
#			print "something is wrong with additional_control regions file"
##			print line
#			sys.exit(1)
#		start = int(start)
#		end = int(end)
#		seq = HG_sequence[chr][start:end]
#		#print seq
#		id = "ptr_hg18_const_%s:%i-%i"%(chr,start,end)
#		#gc_content_vect = gc_content_table["GC_content"][chr][start:end]
#		new_bac = bac(id,seq,chr,start,end,mask_table,gc_content_table,'hg18',cpn)
#		if(chr=="chrx" or chr=="chrX"):
#			new_bac.set_is_x(True)
#		control_bacs[id]=new_bac
#
#		total_unmasked_len+=new_bac.unmasked_len
#		print new_bac.unmasked_len , "\t" , len(new_bac.sequence) , "\t", new_bac.cp
#
#	return total_unmasked_len

def mkdir(dir,file):

	ls_dir = os.listdir(dir)
	if(not(file_exists(ls_dir,file))):
		command = "mkdir %s/%s"%(dir,file)
		print command
		os.system(command)

def die(str):
	print str
	assert(0)


def copy(bac,indiv):
	global genome_info
	if(genome_info.get_sex(indiv) == "F" and bac.is_x):
		return 2*bac.cp
	return bac.cp

def get_depth_vs_cn(cp_num,bac_wssd,hg_wssd,gc_correction_tup,control_bacs,indiv,wssd_dim):
	#filename = "%s.%d"%(prefix,cp_num)
	#print "loading %s..."%(filename)
	#depth_vs_cn	= np.fromfile(filename, 'float', -1, ' ')

	#should now load the GC corrected reads

	(GCp,ave_depths,GCp2_8,lowess_depth,correction,mu) = gc_correction_tup
	bac_fuguized_depths = []
	depth_vs_cn = np.zeros(0)

	print "copy number:%d"%(cp_num)
	for bac in control_bacs:
		curr_bac = control_bacs[bac]
		curr_cp = curr_bac.cp
		if(copy(curr_bac,indiv) == cp_num):
			print "\t%s %d *%s*"%(bac,curr_bac.unmasked_len,curr_bac.source)
			if(curr_bac.source != 'bac'):
				bac_fuguized_depth = get_total_fuguized_adjusted_depth(hg_wssd,"wssd",curr_bac.gc_content_table,"GC_content",bac,correction,curr_bac.mask,curr_bac.chr,curr_bac.start,curr_bac.end,wssd_dim)
			elif(curr_bac.source == 'bac'):
				bac_fuguized_depth = get_total_fuguized_adjusted_depth(bac_wssd,"wssd",curr_bac.gc_content_table,"GC_content",bac,correction,curr_bac.mask,curr_bac.chr,curr_bac.start,curr_bac.end,wssd_dim)
			else:
				print "ERROR undefined source", control_bacs[bac].source
				sys.exit(1)

			print bac_fuguized_depth.shape[0]
			depth_vs_cn = np.r_[depth_vs_cn,bac_fuguized_depth]

	#depth_vs_cn = np.r_[bac_fuguized_depths]
	cn = np.ones(len(depth_vs_cn))*cp_num
	return (depth_vs_cn,cn)


def weighted_fit_line_get_corr(x,y,w,f,color,line_type='-'):

	fitfunc = lambda p, x: p[0] + p[1] * x
	print x
	print y
	print w
	errfunc = lambda p, x, y, w: (y - fitfunc(p, x)) * w

	pinit = [1.0, -1.0]
	out = scipy.optimize.leastsq(errfunc, pinit,args=(x,y,w), full_output=1)
	xs = out[0]

	fit_line = y*xs[1]+xs[0]
	pcorr = scipy.stats.pearsonr(x,fit_line)

	simplex = np.array((0,40))
	simpley = (simplex*xs[1])+xs[0]

	do_plot = False
	if do_plot:
		figure(f)
		plot(simplex,simpley,line_type,color=color)

	return pcorr


def fit_line_get_corr(x,y,f,color,line_type='-'):


	return [1.0]
	do_plot = False
	fitfunc = lambda p, x: p[0] + p[1] * x
	errfunc = lambda p, x, y: (y - fitfunc(p, x))

	pinit = [1.0, -1.0]
	print x.shape[0]
	print y.shape[0]
	out = scipy.optimize.leastsq(errfunc, pinit,args=(x,y), full_output=1)
	xs = out[0]

	fit_line = y*xs[1]+xs[0]
	pcorr = scipy.stats.pearsonr(x,fit_line)

	simplex = np.array((0,np.amax(x)+2))
	simpley = (simplex*xs[1])+xs[0]

	if do_plot:
		figure(f)
		plot(simplex,simpley,line_type,color=color)

	return pcorr

def add_gc_plot_fit_interpolation(gc_correction_tup,color,f1):
	#global GC_width
	#fnGC_hist = "%s/%sGC_vs_depth_CN2"%(indiv_lib_output_dir,prefix)
	#GCp,ave_depths,GCp2_8,lowess_depth,correction,mu = kgf.get_GC_depth_correction_vect(fnGC_hist,GC_width)

	GCp,ave_depths,GCp2_8,lowess_depth,correction,mu = gc_correction_tup
	do_plot = False


	if do_plot:
		figure(f1)

		#lowy = biostats.lowess(GCp,ave_depth,f=.25)
		l=plot(GCp,lowess_depth,"r--")
		l=plot(GCp,ave_depths, "%s."%(color))
		l=plot(np.array([0,1]),np.array([mu,mu]),"%s-"%(color))
		l=plot(GCp,correction,"c-.+")

def get_sdevs(npobs):
	sdevs = []
	for npob in npobs:
		sdevs.append(npob.std())
	return sdevs

def get_avs(npobs):
	avs = []
	for npob in npobs:
		avs.append(npob.mean())
	return avs

def get_cns(control_bacs,indiv):
	cns = []

	for bac in control_bacs:
		curr_bac = control_bacs[bac]
		curr_cp = copy(curr_bac,indiv)

		if(not(curr_cp in cns) and curr_cp < 100):
			cns.append(curr_cp)

	cns.sort()
	print cns
	return cns

#def analyze_data(indiv,lib,indiv_lib_output_dir,prefix,ed,color,f1,f2,f3,ed_struct,offset):
def  analyze_bac_depths(indiv,indiv_lib_output_dir,ed,color,f1,f2,f3,ed_struct,bac_wssd,hg_wssd,gc_correction,control_bacs,wssd_dim):
	#file = "%s.txt.GC_vs_depth_CN2"
	#file = "%s/%sbp_depths"%(indiv_lib_output_dir,prefix)

	offset = 0

	cns = get_cns(control_bacs,indiv)

	cn_arrays = []
	depth_vs_cn_arrays = []
	sum_cn2_read_depth = 0
	for cn in cns:
		(depth_vs_cn_array,cn_array) = get_depth_vs_cn(cn,bac_wssd,hg_wssd,gc_correction,control_bacs,indiv,wssd_dim)
		cn_arrays.append(cn_array)
		depth_vs_cn_arrays.append(depth_vs_cn_array)

		if(cn==2):
			sum_cn2_read_depth = depth_vs_cn_array.sum()

	x_with_x=np.hstack(cn_arrays)
	y_with_x=np.hstack(depth_vs_cn_arrays)

	pcorr_with_x = fit_line_get_corr(x_with_x,y_with_x,f2,color)

	cn_simple = np.array(cns) + offset

	avs = get_avs(depth_vs_cn_arrays)

	sdevs = get_sdevs(depth_vs_cn_arrays)

	simplex = np.array(cns)

	print x_with_x
	print y_with_x
	print x_with_x.shape[0]
	print y_with_x.shape[0]
	pcorr_avs_with_x = fit_line_get_corr(simplex,np.array(avs),f3,color)

	fit_line_get_corr(x_with_x,y_with_x,f3,color,"--")

	ed_struct['corr_all'] = pcorr_with_x[0]
	ed_struct['corr_all_no_x'] = pcorr_with_x[0]
	ed_struct['corr_avs'] = pcorr_avs_with_x[0]
	ed_struct['corr_avs_no_x'] = pcorr_avs_with_x[0]
	ed_struct['us'] = avs
	ed_struct['sigs'] = sdevs


	print x_with_x, len(x_with_x)
	print "y values", y_with_x, len(y_with_x)
	print cn_simple
	print avs

	do_plot = False


	if do_plot:
		figure(f2)
		l=plot(x_with_x+offset,y_with_x, '.',markersize=1)
		#legend(['test'])
		setp(l, 'color',color)

		figure(f3)
		errorbar(cn_simple, avs, sdevs, fmt='o', color=color,ecolor=color)

	return sum_cn2_read_depth


def get_means(list1,list2):
	ret = []
	for i in range(len(list1)):
		if(float(list1[i])==0 or float(list2[i]==0)):
			ret.append(0)
		else:
			ret.append(float(list1[i])/float(list2[i]))
	return ret

#def run_bac_analysis(output_dir,indiv,lib,lib_dir,total_unmasked_len,mask_file):
def run_bac_analysis(output_dir,indiv,lib,total_unmasked_len,gc_correction_tup,bac_wssd,hg_wssd,control_bacs,wssd_dim):

	do_plot = False

	indiv_output_dir = "%s/%s"%(output_dir,indiv)
	mkdir(indiv_output_dir,lib)
	indiv_lib_output_dir = "%s/%s"%(indiv_output_dir,lib)

	eds = ["012"]
	#figure(1)
	colors = ['g']
	i=0
	#FIRST WE MAKE 3 figures, f1, f2, and f3, fill the figures, get the stats, then, finally, make the web page

	if do_plot:
		f1 = figure(1)
		f2 = figure(2)
		f3 = figure(3)
	ed_structs = []
	#offsets = [0,.3,.6]
	summary_file = "%s/summary_stats_dim%d.txt"%(indiv_lib_output_dir,wssd_dim)
	sum_stat_file = open(summary_file,"w");
	output_header = "edit_distance,effective_coverage,corr_all,corr_all_no_x,corr_avs,corr_avs_no_x,cv1,cv2,cv6,cv14,cv16,cv36\n"
	sum_stat_file.write(output_header)

	for ed in eds:

		ed_struct = dict(corr_all=-1,corr_all_no_x=-1,corr_avs=-1,corr_avs_no_x=-1,us=[],sigs=[])
		ed_structs.append(ed_struct)

		(GCp,ave_depths,GCp2_8,lowess_depth,correction,mu) = gc_correction_tup
		add_gc_plot_fit_interpolation(gc_correction_tup,colors[i],1)
		sum_cn2_read_depth = analyze_bac_depths(indiv,indiv_lib_output_dir,ed,colors[i],1,2,3,ed_struct,bac_wssd,hg_wssd,gc_correction_tup,control_bacs,wssd_dim)

		print ed_struct
		print sum_cn2_read_depth, total_unmasked_len
		print "effective coverage: ", float(sum_cn2_read_depth)/total_unmasked_len

		ec = float(sum_cn2_read_depth)/total_unmasked_len
		print ec

		cv = get_means(ed_struct["us"],ed_struct["sigs"])

		#print ed_struct
		output_result = "%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n"%(ed,ec,ed_struct["corr_all"],ed_struct["corr_all_no_x"],ed_struct["corr_avs"],ed_struct["corr_avs_no_x"],0,0,0,0,0)
		sum_stat_file.write(output_result)
		i+=1
		print ed

	print "finished LIB--------------\n-----------\n",lib

	if do_plot:
		#NOW
		figure(1)
		legend(['e012','e012-fit','mean depth','correction factor'])
		figure_name = "%s/GCplot_dim%d.jpg"%(indiv_lib_output_dir,wssd_dim)
		xlim((0,1))
		xlabel("GC%")
		ylabel("average depth")
		savefig(figure_name,format='png')
		print figure_name
		close(1)

		figure(2)
		legend(['e012','e012'],loc=2)
		figure_name = indiv_lib_output_dir+"/bp_read_depth_plot_dim%d.jpg"%(wssd_dim)
		xlabel("cp number")
		ylabel("gc adjusted read depth")
		savefig(figure_name,format='png')
		close(2)

		figure(3)
		legend(['fit on avg depth','fit on all points (weighted by data)'])
		figure_name = indiv_lib_output_dir+"/avg_bp_read_depth_plot_dim%d.jpg"%(wssd_dim)
		xlabel("cp number")
		ylabel("gc adjusted read depth")
		savefig(figure_name,format='png')
		close(3)


class GC_correction:

	def __init__(self,GC_width):
		self.sum_read_depth = {}
		self.n_locations = {}

		for i in range(2*GC_width+2):
			self.sum_read_depth[i] = 0
			self.n_locations[i] = 0

	def add_wssd_file(self,wssd,bac,wssd_dim):
		#depth = np.loadtxt(file,'float',usecols=[2])
		depth = wssd.depth["wssd"][bac.chr][bac.start:bac.end,wssd_dim].sum(1)

		for n_gc,sum_depth in self.sum_read_depth.iteritems():

			#print "****"
			#print "n_locations",self.n_locations[n_gc]
			#print "sum_read_depth",self.sum_read_depth[n_gc]

			equiv = ((bac.GC_count==n_gc) & (bac.mask==0)) #make sure NOT masked
			#print "SUM SUM" , equiv.sum()
			#print depth[equiv]
			#print depth[equiv].sum()

			self.sum_read_depth[n_gc]+=(depth[equiv]).sum()
			self.n_locations[n_gc]+=equiv.sum()
			#print "n_gc",n_gc
			#print "GC_count",bac.GC_count[0:100]
			#print "depth",(depth[0:100])[(equiv[0:100])]
			#print "n_locations",self.n_locations[n_gc]
			#print "sum_read_depth",self.sum_read_depth[n_gc]

	def get_and_output_udepth_gc(self,fn):
		f = open(fn,'w')

		#print self.sum_read_depth
		#print self.n_locations
		#print GC_width
	 	correction_tuple=kgf.get_GC_depth_correction_vect(self.sum_read_depth.values(),self.n_locations.values(),GC_width)

		(GCp,ave_depths,GCp2_8,lowess_depth,correction,mu) = correction_tuple
		#print GCp
		#print ave_depths
		#print GCp2_8
		#print lowess_depth
		#print correction
		#print mu
		#print "*********"
		f.write("%d\n"%(GC_width))

		for i in range(len(GCp)):
			f.write("%f %f\n"%(GCp[i],correction[i]))

	#	f.write("%d %f %d %d\n"%(n_gc,udepth[-1],self.sum_read_depth[n_gc],self.n_locations[n_gc]))
		return correction_tuple


def calc_gc_correction(indiv_lib_output_dir,indiv,bac_wssd,hg_wssd,control_bacs,GC_width,wssd_dim):
	#GC correctino is calculated for the 2 copy bacs
	GC_corr = GC_correction(GC_width)
	for bac_name,bac in control_bacs.iteritems():
		if(bac.cp==2):
			if(bac.source == 'hg18'):
				GC_corr.add_wssd_file(hg_wssd,bac,wssd_dim)
			elif(bac.source == 'bac'):
				GC_corr.add_wssd_file(bac_wssd,bac,wssd_dim)
			else:
				print "error, undefinced ousrce ",bac.source
				sys.exit(1)

	fn = "%s/cn2-gc-depth-dim%d.txt"%(indiv_lib_output_dir,wssd_dim)
	correction_tuple =  GC_corr.get_and_output_udepth_gc(fn)
	return correction_tuple
	#print totalu
	#print udepth

class hist:

	def __init__(self,cn):
		self.cn = cn
		self.bins = zeros(3000)

	def add_file(self,file):
		print "loading..." , file
		data = np.loadtxt(file)
		self.bins+=data

	def get_pdf(self):
		self.pdf = self.bins/np.sum(self.bins)

	def add_to_plot(self,plot_handle):
		do_plot = False

		if do_plot:
			figure(plot_handle)
			self.get_pdf()
			self.x_axis = arange(3000);
			plot(self.x_axis,self.pdf)

	def build_from_data(self,data):
		self.bins = np.histogram(data,bins=3000,range=[0,3000])[0].astype(float64)

class fuguized_bac_collector:

	###MAKE SURE ADJUSTED~
	def __init__(self,bac_wssd,hg_wssd,control_bacs,cn,indiv,gc_correction_vect,wssd_dim):

		self.data = np.zeros(0)
		self.cn = cn
		for bac in control_bacs:
			curr_bac = control_bacs[bac]
			#print bac, curr_bac.chr, curr_bac.start, curr_bac.end, curr_bac.source
			if(copy(curr_bac,indiv) == cn):
				if(curr_bac.source=="bac"):
					bac_fuguized_depth = get_total_fuguized_adjusted_depth(bac_wssd,"wssd",curr_bac.gc_content_table,"GC_content",bac,gc_correction_vect,curr_bac.mask,curr_bac.chr,curr_bac.start,curr_bac.end,wssd_dim)
				elif(curr_bac.source=="hg18"):
					bac_fuguized_depth = get_total_fuguized_adjusted_depth(hg_wssd,"wssd",curr_bac.gc_content_table,"GC_content",bac,gc_correction_vect,curr_bac.mask,curr_bac.chr,curr_bac.start,curr_bac.end,wssd_dim)
				else:
					print "source error " , curr_bac.source
					sys.exit(1)

				self.data = r_[self.data,bac_fuguized_depth]
				#self.data = r_[self.data,control_bacs[bac].get_unmasked_depth(bac_wssd.depth["wssd"][bac][:,:,0].sum(1))]

	def add_file(self,bac_wssd,hg_wssd,control_bacs,gc_correction_vect,indiv,wssd_dim):

		added_wssd = np.zeros(0)
		for bac in control_bacs:
			curr_bac = control_bacs[bac]
			if(copy(curr_bac,indiv) == self.cn):
				#bac_fuguized_depth = get_total_fuguized_adjusted_depth(bac_wssd,"wssd",gc_content_table,"GC_content",bac,gc_correction_vect,control_bacs[bac].mask)
				if(curr_bac.source=="bac"):
					bac_fuguized_depth = get_total_fuguized_adjusted_depth(bac_wssd,"wssd",curr_bac.gc_content_table,"GC_content",bac,gc_correction_vect,curr_bac.mask,curr_bac.chr,curr_bac.start,curr_bac.end,wssd_dim)
				elif(curr_bac.source=="hg18"):
					bac_fuguized_depth = get_total_fuguized_adjusted_depth(hg_wssd,"wssd",curr_bac.gc_content_table,"GC_content",bac,gc_correction_vect,curr_bac.mask,curr_bac.chr,curr_bac.start,curr_bac.end,wssd_dim)

				added_wssd = r_[added_wssd,bac_fuguized_depth]
				#added_wssd = r_[added_wssd,control_bacs[bac].get_unmasked_depth(bac_wssd.depth["wssd"][bac][:,:,0].sum(1))]

		self.data += added_wssd

	def get_hist(self):
		new_hist = hist(self.cn)
		new_hist.build_from_data(self.data)
		return new_hist


def do_fit_indiv(cn,fuguized_bac,fig_num):

	#WEIBULL = k/LAMBDA*((x/lambda)^k-1)*exp(-x/lambda)^k)
	#k lambda
	#fitfunc = lambda p, x: (p[0]/p[1])*((x/p[1])**(p[0]-1))*exp(-(x/p[1])**p[0])

	#mu sigma
	fitfunc = lambda p, x: (1/(p[1]*sqrt(2*pi))*exp(-((x-p[0])**2)/(2*p[1]**2)))
	errfunc = lambda p, x, y: (y - fitfunc(p, x))

	h=fuguized_bac.get_hist()
	h.get_pdf()

	x=np.arange(3000)
	y=h.pdf
	pinit = [60, 100]
	print "doing fit"
	out = scipy.optimize.leastsq(errfunc, pinit,args=(x,y), full_output=1)
	print out
	xs = out[0]

	print xs
	fit_ys = fitfunc(xs,x)

	figure(fig_num)
	plot(x,fit_ys,"--")


	#fit_line = y*xs[1]+xs[0]
	#pcorr = scipy.stats.pearsonr(x,fit_line)
	return xs


def do_fit_curves(hist_labels,cns,fuguized_bac_collection,fig_num,indiv_directory,wssd_dim):

	fit_params_u = []
	total_bp = []
	fit_params_sig = []
	strhists = []
	for i in range(len(hist_labels)):
		cn = cns[i]
		label = hist_labels[i]
		strhists.append(str(cn))
		print "fitting parameters for ", cn , "..."
		fp=do_fit_indiv(cn,fuguized_bac_collection[label],fig_num);
		fit_params_u.append(fp[0])
		fit_params_sig.append(fp[1])
		total_bp.append(fuguized_bac_collection[label].data.shape[0])

	legend(strhists)
	xlim((0,1000))
	mkdir(indiv_directory,"_filtered_libs_analysis")
	figurename = "%s/_filtered_libs_analysis/combined_libs_hist_fits-dim%d.png"%(indiv_directory,wssd_dim)
	savefig(figurename,format='png')
	print figurename
	close(fig_num)

	np_hist_cns = np.array(cns)
	np_fit_params_u =np.array(fit_params_u)
	np_fit_params_sig =np.array(fit_params_sig)
	np_total_bp = np.array(total_bp)
	w = np_total_bp.astype(float64)/np_total_bp.sum()

	figure(1)
	figure(2)
	print w
	pcorr_u=weighted_fit_line_get_corr(np_hist_cns,np_fit_params_u,w,1,"g")
	pcorr_sig=weighted_fit_line_get_corr(np_hist_cns,np_fit_params_sig,w,2,"g")
	#pcorr_u=fit_line_get_corr(np_hist_cns,np_fit_params_u,1,"g")
	#pcorr_sig=fit_line_get_corr(np_hist_cns,np_fit_params_sig,2,"g")

	figure(1)
	figurename = "%s/_filtered_libs_analysis/fit_u_dim%d.png"%(indiv_directory,wssd_dim)
	plot(np_hist_cns,np_fit_params_u,"r.")
	t = "u vs cn, pearson correlation of:",pcorr_u[0]
	title(t)
	savefig(figurename,format='png')
	print figurename
	close(1)

	figure(2)
	figurename = "%s/_filtered_libs_analysis/fit_sig_dim%d.png"%(indiv_directory,wssd_dim)
	plot(np_hist_cns,np_fit_params_sig,"r.")
	t = "sig vs cn, pearson correlation of:",pcorr_sig[0]
	title(t)
	savefig(figurename,format='png')
	print figurename
	close(2)

	print "u corr:",pcorr_u
	print "sig corr:", pcorr_sig

	fn_fit_u_sigs = "%s/_filtered_libs_analysis/fit_params_dim%d"%(indiv_directory,wssd_dim)
	f_fit_u_sigs = open(fn_fit_u_sigs,"w")

	for k in range(len(cns)):
		f_fit_u_sigs.write("%d %f %f %d\n"%(cns[k],fit_params_u[k],fit_params_sig[k],np_total_bp[k]))
	f_fit_u_sigs.close()








#########
#######DEPRECATED
#def combine_libs_analyze_hist(wssd_bac_dir,wssd_hg_dir,indiv_directory,indiv,control_bacs,fnbac_contigs,fn_hg18_contigs):
#
#	accepted_libs_list = []
#	ls_libs = os.listdir(indiv_directory)
#
#
#	cns = get_cns(control_bacs,indiv)
#	hist_labels = cns
#
#	fuguized_bac_collection = {}
#	fuguized_bac_collection_starts = {}
#	#fuguized_bac_collection[i] = fugized_bac_collector(
#
#	n_accepted_libs = 0
#	for lib in ls_libs:
#		if(lib[0:1] == "_" or lib[0:3]=="do_"):continue
#
#		analysis_lib_dir = "%s/%s"%(indiv_directory,lib)
#		wssd_bac_lib_dir = "%s/%s"%(wssd_bac_dir,lib)
#		fn_bac_wssd = "%s/control_bac_wssd"%(wssd_bac_lib_dir)
#
#		wssd_hg_lib_dir = "%s/%s"%(wssd_hg_dir,lib)
#		fn_hg_wssd = "%s/hg18rmsk.wssd"%(wssd_hg_lib_dir)
#
#		if(simple_check_failed_lib(analysis_lib_dir)):
#			continue
#
#		gc_correction_vect = load_gc_correction_vect(analysis_lib_dir,0)
#		gc_correction_vect_starts = load_gc_correction_vect(analysis_lib_dir,1)
#
#		bac_wssd = WssdFile(fnbac_contigs,
#										fn_bac_wssd,
#										overwrite = False,
#										openMode = 'r' )
#
#		hg_wssd = WssdFile(fn_hg18_contigs,
#										fn_hg_wssd,
#										overwrite = False,
#										openMode = 'r' )
#
#		if(kgf.lib_pass_qc(analysis_lib_dir)):
#			print "LIB PASSED"
#
#			n_accepted_libs +=1
#			accepted_libs_list.append(lib)
#
#			for i in range(len(hist_labels)):
#				cn = cns[i]
#				label = hist_labels[i]
#				print "making fuguized_bac_collector %d"%(label)
#				if(label in fuguized_bac_collection):
#					fuguized_bac_collection[label].add_file(bac_wssd,hg_wssd,control_bacs,gc_correction_vect,indiv,0)
#					fuguized_bac_collection_starts[label].add_file(bac_wssd,hg_wssd,control_bacs,gc_correction_vect_starts,indiv,1)
#				else:
#					fuguized_bac_collection[label] = fuguized_bac_collector(bac_wssd,hg_wssd,control_bacs,cn,indiv,gc_correction_vect,0)
#					fuguized_bac_collection_starts[label] = fuguized_bac_collector(bac_wssd,hg_wssd,control_bacs,cn,indiv,gc_correction_vect_starts,1)
#
#		else:
#			print "LIB BAD!"
#		print "**********"
#
#		####################
#
#		hists = {}
#		for i in range(len(hist_labels)):
#			print i
#			cn = cns[i]
#			label = hist_labels[i]
#			single_fbc = fuguized_bac_collector(bac_wssd,hg_wssd,control_bacs,cn,indiv,gc_correction_vect,0)
#			newhist = single_fbc.get_hist()
#			hists[label] = newhist
#
#		strhists = []
#
#		do_plot = False
#		if do_plot:
#			figure(1)
#
#		for i in hist_labels:
#			hists[i].add_to_plot(1)
#			strhists.append(str(i))
#
#		legend(strhists)
#		xlim((0,1000))
#		figurename = "%s/hist_reads.png"%(analysis_lib_dir)
#		if do_plot:
#			savefig(figurename,format='png')
#			print figurename
#			close(1)
#
#		bac_wssd.tbl.close()
#		hg_wssd.tbl.close()
#		del hg_wssd
#		del bac_wssd
	#####################
#
#	if(n_accepted_libs==0):
##		print "no accepted libs"
#	return
#
#
#	if do_plot:
#		f=figure(1)
#		f.set_figwidth(50)
#		f.set_figheight(30)
#
#	strhists=[]
#	for i in hist_labels:
##		print "accessing fuguized_bac_collection %d"%(i)
#		hist_inst = fuguized_bac_collection[i].get_hist()
#		hist_inst.add_to_plot(1)
#		strhists.append(str(i))
#
#	do_fit_curves(hist_labels,cns,fuguized_bac_collection,1,indiv_directory,0)
#
#	if do_plot:
#		f=figure(2)
#		f.set_figwidth(10)
#		f.set_figheight(5)
#	strhists=[]
#	for i in hist_labels:
#		print "accessing fuguized_bac_collection %d"%(i)
#		hist_inst = fuguized_bac_collection_starts[i].get_hist()
#		hist_inst.add_to_plot(2)
#		strhists.append(str(i))
#
#	do_fit_curves(hist_labels,cns,fuguized_bac_collection_starts,2,indiv_directory,1)


def get_lib_corr_vs_cov(dir,ed=3):
	file = open("%s/summary_stats_dim0.txt"%(dir),"r")
	line = file.readline()
	for i in range(ed):
		line = file.readline()

	sline = line.split(",")
	return (float(sline[2],sline[1]))



def get_indiv_corr_vs_cov(indiv,indiv_dir):

	lib_ls = os.listdir(indiv_dir)
	for lib in lib_ls:
		lib_dir = "%s/%s"%(indiv_dir,lib_dir)
		(lib_corr,eff_cov) = get_lib_corr_vs_cov(dir)


def get_batch(batches_dir,indiv):
	ls_batches = os.listdir(batches_dir)
	for batch in ls_batches:
		ls_indivs = os.listdir("%s/%s"%(batches_dir,batch))
		for i in ls_indivs:
			if(i==indiv):
				return batch

	return None



##########DEPRECATED

#def select_libs_build_hists(input_bac_dir,input_hg_dir,analysis_dir,control_bacs,fnbac_contigs,fn_hg18_contigs,subset):
#		#select_libs_build_hists(o.input_bac_dir,o.input_hg18_dir,analysis_dir,control_bacs,o.fn_bac_contigs,o.fn_hg18_contigs)
#
#	#FIRST GO THROUGH ALL THE LIBS AND build a plot of correlation vs effective coverage ignoring the x chrom
#
#	ls_indivs = os.listdir(analysis_dir)
#	for indiv in ls_indivs:
#		if(subset != None):
#			if(not(indiv in subset)): continue
#		indiv_dir = "%s/%s"%(analysis_dir,indiv)
#		os.system("rm -r  %s/_filtered_libs_analysis"%(indiv_dir))
#		#wssd_bac_dir = "%s/%s/%s"%(input_bac_dir,get_batch(batches_dir,indiv),indiv)
#		#wssd_hg18_dir = "%s/%s/%s"%(input_hg18_dir,get_batch(batches_dir,indiv),indiv)
#		wssd_bac_dir = "%s//%s"%(input_bac_dir,indiv)
#		wssd_hg18_dir = "%s/%s"%(input_hg_dir,indiv)
#		#(lib_corrs,effectiv_covs) = get_indiv_corr_vs_cov(indiv,indiv_dir)
#		combine_libs_analyze_hist(wssd_bac_dir,wssd_hg18_dir,indiv_dir,indiv,control_bacs,fnbac_contigs,fn_hg18_contigs)

def simple_check_failed_lib(indiv_lib_output_dir):
	summary_file = "%s/summary_stats_dim0.txt"%(indiv_lib_output_dir)
	sum_stat_file = open(summary_file,"r");
	rd = float(sum_stat_file.readlines()[1].split(",")[1])
	if(rd==0):
		return True
	return False


def check_failed_lib(bac_wssd,hg_wssd,control_bacs):

	total_sum = 0
	for bac in control_bacs:
		curr_bac = control_bacs[bac]
		if(curr_bac.source!="bac"):
			total_sum+=hg_wssd.depth["wssd"][curr_bac.chr][curr_bac.start:curr_bac.end,:,0].sum(1).sum()
		elif(curr_bac.source=="bac"):
			total_sum+=bac_wssd.depth["wssd"][curr_bac.chr][curr_bac.start:curr_bac.end,:,0].sum(1).sum()
		else:
			die("terrible horrendous error")

		if(total_sum>10000):
			return False
	return True

def make_failed_output(indiv_lib_output_dir):
	summary_file = "%s/summary_stats_dim0.txt"%(indiv_lib_output_dir)
	sum_stat_file = open(summary_file,"w");
	output_header = "edit_distance,effective_coverage,corr_all,corr_all_no_x,corr_avs,corr_avs_no_x,cv1,cv2,cv6,cv14,cv16,cv36\n"
	sum_stat_file.write(output_header)
	output_header = "0,0,0,0,0,0,0,0,0,0,0,0\n"
	sum_stat_file.write(output_header)


def	run_batch_analysis(genome_info,output_dir,total_unmasked_len,input_hg_dir,control_bacs,GC_width,fn_hg18_contigs,wssd_dim,input_wssd_name,curr_genome):

	print input_hg_dir

	ls = os.listdir("%s"%(input_hg_dir))
	for dir in ls:
		if not dir==curr_genome: continue
		indiv = dir

		indiv_directory="%s/%s"%(input_hg_dir,indiv)
		mkdir(output_dir,indiv)
		ls_libs = os.listdir(indiv_directory)
		for lib in ls_libs:

			lib_bac_dir = "%s/%s"%(indiv_directory,lib)
			if(lib[0:1] == "_" or lib[0:3]=="do_" or not(os.path.isdir(lib_bac_dir))):continue

			indiv_output_dir = "%s/%s"%(output_dir,indiv)
			indiv_lib_output_dir = "%s/%s"%(indiv_output_dir,lib)
			#os.system("rm -r %s"%(indiv_lib_output_dir))
			mkdir(indiv_output_dir,lib)

			lib_hg18_dir = "%s/%s/%s"%(input_hg_dir,indiv,lib)
			#wssd_dir = "%s/wssd_out/"%(lib_dir)
			wssd_hg18_dir = "%s/wssd_out/"%(lib_hg18_dir)
			eds = ["0","01","012"]
			print "analyzing %s \noutputing to %s"%(lib_hg18_dir,indiv_lib_output_dir)
			#OPEN THE WSSD FILE

			fnhg_wssd = "%s/%s"%(lib_hg18_dir,input_wssd_name)
			#fnhg_wssd = "%s/hg18rmsk.wssd"%(lib_hg18_dir)

			print fnhg_wssd
			#bac_wssd = WssdFile(fn_bac_contigs,
			#								fnbac_wssd,
			#								overwrite = False,
			#								openMode = 'r' )

			bac_wssd = None


			hg_wssd = WssdFile(fn_hg18_contigs,
											fnhg_wssd,
											overwrite = False,
											openMode = 'r' )

			if(check_failed_lib(bac_wssd,hg_wssd,control_bacs)):
				make_failed_output(indiv_lib_output_dir)
				continue

			#make_bac_plots(control_bacs,bac_wssd,indiv_lib_output_dir,mask,mask_grps)
			print "output_dir %s"%(output_dir)
			gc_depths_file = "%s/%s/%s/cn2-gc-depth-WG-W%d-dim%d.txt"%(output_dir,indiv,lib,GC_width,wssd_dim)
			print gc_depths_file
			#gc_correction_tup = #calc_gc_correction(indiv_lib_output_dir,indiv,bac_wssd,hg_wssd,control_bacs,GC_width,wssd_dim)
			gc_correction_tup = kgf.get_GC_depth_correction_from_file(gc_depths_file,3)
			run_bac_analysis(output_dir,indiv,lib,total_unmasked_len,gc_correction_tup,bac_wssd,hg_wssd,control_bacs,wssd_dim)

			hg_wssd.tbl.close()
			del hg_wssd

def get_lib_cov_corr(lib_dir):

	f=open("%s/summary_stats_dim0.txt"%(lib_dir),"r")
	line = f.readline()
	line = f.readline()

	if(not(line)):
		return None

	sline = line.split(",")
	ec = sline[1]
	corr = sline[2] #INCLUDES X!!!!

	return(ec,corr)




def analyze_coverage_to_correlation(analysis_dir):

	#anal_by_alg	= {}
	print "************************"
	print "analyze_coverage_to_correlation"
	print "************************"

	s_analysis_dir =  analysis_dir.split("/")
	pilot = s_analysis_dir[len(s_analysis_dir)-2]

	print s_analysis_dir
	pass_thresh = 0;
	total = 0
	ls_algs = os.listdir(analysis_dir)
	for alg in ls_algs:
		if(alg != "mrsfast"):
			continue
		print analysis_dir
		alg_dir = "%s/%s"%(analysis_dir,alg)
		print alg_dir
		ls_indivs = os.listdir(alg_dir)

		anal_by_indivs = {}

		for indiv in ls_indivs:
			ecs = []
			corrs = []
			ecs_corrs = {}
			ecs_corrs["ecs"] = ecs
			ecs_corrs["corrs"] = corrs

			anal_by_indivs[indiv] = ecs_corrs

			indiv_dir = "%s/%s"%(alg_dir,indiv)
			ls_libs = os.listdir(indiv_dir)
			for lib in ls_libs:
				#if(lib[0:1] == "_"): continue
				if(lib[0:1] == "_" or lib[0:3]=="do_" or lib[0]=="."):continue
				total+=1
				lib_dir = "%s/%s"%(indiv_dir,lib)
				if(get_lib_cov_corr(lib_dir) != None):
					if(kgf.lib_pass_qc(lib_dir)):
						pass_thresh+=1
					(ec,corr) = get_lib_cov_corr(lib_dir)
					ecs.append(ec)
					corrs.append(corr)

	do_plot  = False

	if do_plot:
		figure(1)
	indivs = []
	for indiv in anal_by_indivs:
		indivs.append(indiv)
		ecs_corrs = anal_by_indivs[indiv]
		ecs = np.array(ecs_corrs["ecs"])
		corrs = np.array(ecs_corrs["corrs"])
		marker = '.'
		if(indiv[2:4] == "19" and pilot=="Pilot2"):
			marker = '^'
		if do_plot:
			plot(ecs,corrs,marker)


	thresh = kgf.get_lib_pass_thresh()
	plot(np.array((0,12)),np.array((thresh,thresh)),"r--")
	indivs.append("threshold")
	if(pilot=="Pilot2"):
		legend(indivs,loc=4)
	xlabel("effective coverage bp")
	ylabel("correlation with bacs (X inc)")
	title("%s - correlation vs ec - %d/%d meet threshold"%(pilot,pass_thresh,total))
	figurename = "%s/%s"%(analysis_dir,"/_all_alg_analysis/ec_vs_corr_dim0.png")


	if do_plot:
		savefig(figurename,format='png')
		print figurename
		close(1)


#def run_dist_analysis(index_file,total_unmasked_len,control_bacs,batches,input_dir,output_dir,sex_index,fnbac_contigLengths,mask,mask_grps,GC_width,gc_content_table):
def run_dist_analysis(genome_info,total_unmasked_len,control_bacs,output_dir,input_hg_dir,GC_width,fn_hg_contigs,wssd_dim,input_wssd_name,curr_genome):
	#kgf.init(index_file,2,sex_index)
	#for batch in batches:
	run_batch_analysis(genome_info,output_dir,total_unmasked_len,input_hg_dir,control_bacs,GC_width,fn_hg_contigs,wssd_dim,input_wssd_name,curr_genome)


def	load_table(fn,fn_contig):
	msk = DenseTrackSet(
                  fn_contig,
                  fn,
                  overwrite=False,
                  openMode='r')
	return msk


def add_depth_to_fig(fignum,depth,sp,name,n,mean,median,mask,mask_grps):

	rm = mask["mask"][name][:,0]
	wm = mask["mask"][name][:,1]
	trf = mask["mask"][name][:,2]
	#print rm
	#print wm
	#print trf

	#rm = individual_masks["repeat_mask"].contig[name]
	#wm = individual_masks["window_mask"].contig[name]
	#trf = individual_masks["trf_mask"].contig[name]

	#figure(fignum)
	subplot(n,1,sp)
	plot(depth,"-",linewidth=.25)
	plot(np.array([0,depth.shape[0]]),np.array([mean,mean]),'r--')
	plot(np.array([0,depth.shape[0]]),np.array([median,median]),'g-*')

	plot(rm+wm+trf,'-r',linewidth=.25)
	plot(rm*-1,'g-',linewidth=.25)
	plot(wm*-2,'m-',linewidth=.25)
	plot(trf*-1.5,'c-',linewidth=.25)

	legend([name])
	#ylim((0,2*depth.mean()))
	xlabel(name)


def make_bac_plots(bacs,wssd,analysis_dir_out,mask,mask_grps):

	do_plot = False
	##############UNDER CONSTRUCTION
	#####INPUT -> bacs wssd, plot the wssd

	#masking_study_dir = "/net/eichler/ebod/eichler_testing/kitz/wssd/misc/aug2009_masking_study/"
	#genome_collection = "non1kg"

	#mkdir(nonkg_output_dir,masking_type)
	#nonkg_out = "%s/%s"%(nonkg_output_dir,masking_type)

	#mkdir(kg_output_dir,masking_type)
	#Ekg_out = "%s/%s"%(kg_output_dir,masking_type)

	n_bacs = len(bacs)

	if do_plot:
		f=figure(1)


		f.set_figwidth(100)
		f.set_figheight(100)
		subplot(n_bacs,1,1)
		subplots_adjust(hspace=.02)


	sp = 1
	for bac_name,bac in bacs.iteritems():
		print "plotting %s..."%(bac_name)
		assert(bac_name==bac.name)
		#bac = "control_bac_perbase_depth_ed012.txt.%s.txt"%(bac.name)
		depth = wssd.depth["wssd"][bac.name][:,:,0].sum(1)

		if do_plot:
			add_depth_to_fig(1,depth,sp,bac.name,n_bacs,bac.get_unmasked_mean(depth),bac.get_unmasked_median(depth),mask,mask_grps)
		#print "max at,", np.where(depth==depth.max())
		#print np.where(depth==depth.max())[0][0]
		#print bac.sequence[np.where(depth==depth.max())[0].min()-200:np.where(depth==depth.max())[0].max()]
		sp+=1

	#legend(leg)
	if do_plot:
		figure_name = "%s/lib_specific_depths.jpg"%(analysis_dir_out)
	#subplots_adjust(wspace=90000)
	#xlabel("pos")
	#f.set_figheight(100)

	#fset_size_inches((100,10))
		ylabel("depth")
		print figure_name
		savefig(figure_name,format='png')
		close(1)


if __name__=='__main__':

	control_bacs = {}

	#bac_list_file = "./baccopies"
	#mask_file = "/net/gs/vol1/home/psudmant/EEE_Lab/control_bacs/bacs_aggressivemask.ntable"
	#mask_file = "/net/gs/vol1/home/psudmant/genomes/masking/output/control_bacs-rm-100.0wm.mask_table"
	fnbac_contigLengths = "/net/gs/vol1/home/psudmant/EEE_Lab/1000G/1000genomesScripts/kitz_wssd/wssd_contigLengths_BACS.txt"

	opts = OptionParser()
	opts.add_option('','--sex_index',dest='fn_sex_pop_index')
	#opts.add_option('','--index_file',dest='index_file')

	#opts.add_option('','--input_bac_dir',dest='input_bac_dir')
	#opts.add_option('','--input_hg18_dir',dest='input_hg18_dir')
	#opts.add_option('','--analysis_root_dir',dest='analysis_root_dir')
	#opts.add_option('','--algorithm',dest='algorithm')
	opts.add_option('','--do_dist_analysis',default=False,action='store_true',dest='do_dist_analysis')
	opts.add_option('','--analyze_coverage_to_correlation',default=False,action='store_true',dest='analyze_cov_to_corr')
	opts.add_option('','--build_hists',default=False,action='store_true',dest='build_hists')
	#opts.add_option('','--batches',dest='batches')
	#opts.add_option('','--additional_control_regions',dest='fn_additional_control_regions')

	opts.add_option('','--hg_mask_file',dest='fn_hg_mask')
	opts.add_option('','--hg_contig',dest='fn_hg_contigs')
	opts.add_option('','--hg_gc_content',dest='fn_hg_gc_content_table')
	opts.add_option('','--genome',dest='genome')

	#opts.add_option('','--bac_mask_file',dest='fn_bac_mask')
	#opts.add_option('','--bac_contig',dest='fn_bac_contigs')
	#opts.add_option('','--bac_gc_content',dest='fn_bac_gc_content_table')
	opts.add_option('','--bac_list_file',dest='fn_bac_list_file')
	opts.add_option('','--input_genomes',dest='fn_input_genomes',default=None)
	opts.add_option('','--gc_width',dest='GC_width',default=200)
	opts.add_option('','--input_wssd_file_names',dest='input_wssd_file_names')


	(o, args) = opts.parse_args()
	GC_width = o.GC_width
	#analysis_dir = "%s/%s"%(o.analysis_root_dir,o.algorithm)
	#output_dir = analysis_dir

	#bac_mask = load_table(o.fn_bac_mask,o.fn_bac_contigs)
	#bac_gc_content_table = load_table(o.fn_bac_gc_content_table,o.fn_bac_contigs)

	hg_mask = load_table(o.fn_hg_mask,o.fn_hg_contigs)
	hg_gc_content_table = load_table(o.fn_hg_gc_content_table,o.fn_hg_contigs)

	#(total_unmasked_len,control_bacs) = load_control_bac_objects(bac_list_file,bac_mask,bac_gc_content_table)
	(total_unmasked_len,control_bacs) = load_control_bac_objects(o.fn_bac_list_file,hg_mask,hg_gc_content_table,o.genome)

	print len(control_bacs)
	print total_unmasked_len
	#(additional_unmasked_len) = load_additional_control_regions(control_bacs,o.fn_additional_control_regions,hg18_mask,hg18_gc_content_table)
	#print "additional", additional_unmasked_len
	#total_unmasked_len+=additional_unmasked_len

	print len(control_bacs)

	mask_grps = ["RM","WM","TRF"]

	genome_info = kgf.genome_info(o.fn_input_genomes,o.fn_sex_pop_index)

	for in_genome_line in open(o.fn_input_genomes,'r').readlines():
		(genome_name,fn_wssd_dir,fn_bac_dir,chunk_dir,genome_out) = in_genome_line.split()

		print genome_info.genomes[genome_name].sex
		print genome_name

		indiv = genome_name.split(".")[0]

		if(o.do_dist_analysis):
			#batches = o.batches.rstrip("\n").split(":")
			#run_dist_analysis(o.index_file,total_unmasked_len,control_bacs,batches,o.input_dir,output_dir,o.sex_index,fnbac_contigLengths,mask,mask_grps,GC_width,gc_content_table)
			print "running dist analysis..."
			run_dist_analysis(genome_info,total_unmasked_len,control_bacs,fn_bac_dir,fn_wssd_dir,o.GC_width,o.fn_hg_contigs,0,o.input_wssd_file_names,indiv)
			run_dist_analysis(genome_info,total_unmasked_len,control_bacs,fn_bac_dir,fn_wssd_dir,o.GC_width,o.fn_hg_contigs,1,o.input_wssd_file_names,indiv)
		if(o.analyze_cov_to_corr):
			analyze_coverage_to_correlation(o.analysis_root_dir)
		if(o.build_hists):
			print "ARE YOU SURE YOU WANT TO BUILD HISTS!"
			print "CURRENTLY the FIT PARAMETERS SCRIPT is building hists"
			#select_libs_build_hists(o.input_bac_dir,o.input_hg18_dir,analysis_dir,control_bacs,o.fn_bac_contigs,o.fn_hg18_contigs,subset)


		exit(0)

		#analyze_mask_study()
		#exit(1)

