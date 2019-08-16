import sys,os
import subprocess
from datetime import datetime
import time
import subprocess
import numpy as np
import scipy as sp


import scipy.optimize
#import scipy.stats

#import scipy.interpolate.fitpack2 as fp2 

import Bio.Statistics.lowess as biostats

def mkdir_file_exists(ls,file):
    for f in ls:
        if(f==file):
            return 1
    return 0

def mkdir(dir,file):
    ls_dir = os.listdir(dir)
    if(not(mkdir_file_exists(ls_dir,file))):
        command = "mkdir %s/%s"%(dir,file)
        os.system(command)
    return "%s/%s"%(dir,file)

class library:
    def __init__(self,name,paired,insert_size,center):
        self.center = center
        self.name = name
        self.files = []
        self.paired = paired
        if(self.paired == "PAIRED" and insert_size != ""):
            self.insert_size = int(insert_size)
        else:
            self.insert_size = -1

    def addfile(self,file):
        self.files.append(file)
    def update_paired(self,paired,insert_size):
        if(paired != self.paired):
            self.paired = "BOTH"
        if(self.paired == "BOTH" or self.paired == "PAIRED"):
            if(insert_size != "" and insert_size != "0"):
                if(self.insert_size!=-1 and (self.insert_size!=int(insert_size)) and assert_lib_insert==1):
                    print("#",str(self.insert_size) , "UPDATING TO " , insert_size)
                    print(self.name)
                self.insert_size = int(insert_size)
                #self.insert_size = 0

class file:
    
    def get_read_len(self,output_dir):
        command  = "zcat -q  %s/%s/%s/%s | head -q 2>/dev/null" %(output_dir,self.pilot_dir,self.sample,self.name)
        stdout_handle = os.popen(command, "r")
        text = stdout_handle.read()
        stext = text.split("\n")
        return len(stext[1])
    def __init__(self,path,md5sum,mate_pair,paired,pilot_dir,withdrawn=0):
        self.pilot_dir = pilot_dir
        spath = path.split("/")
        self.path = path 
        self.name = spath[len(spath)-1] 
        self.withdrawn=withdrawn
        self.md5sum = md5sum
        self.sample = spath[1]
        self.mate_pair = mate_pair
        self.paired = paired
        #now get the read length

class genome:
    def __init__(self,name,population,sex):
            self.name = name    
            self.libraries = []
            self.population = population
            self.sex = sex
    def addlibrary(self,library):
        self.libraries.append(library)


def run_file_check(md5sum=0):
    
    #check that all files exist
    if(check_files_exist()):
        print("\n\n\n------------\nFILES MISSING\n------------\n\n\n")
    #check that there are no files to be removed
    if(check_no_withdrawn()):
        print("\n\n\n------------\nWITHDRAWN_FILES_EXIST\n------------\n\n\n")
    #check_sums    
    if(md5sum):
        if(check_sums()):
            print("\n\n\n------------\nCHECKSUMS FAILED\n------------\n\n\n")
            


def get_line_params(x,y):
    
    line_func = lambda p, x: (p[0]*x+p[1])
    errfunc = lambda p, x, y: (y - line_func(p, x))
    pinit = [0,0]
    out = scipy.optimize.leastsq(errfunc, pinit,args=(x,y), full_output=1)
    xs = out[0]

    #    m = (y[1]-y[0])/(x[1]-x[0])
    #b = y[1]-m*x[1]    
    #return m,b
    return xs

def get_GC_depth_correction_from_dir(fn_bac_dir,max_scale_factor=0,correct_range=False):
    fn_bac_gc_depths = "%s/cn2-gc-depth-WG-W200-dim0.txt"%(fn_bac_dir)

    gc_tuple = get_GC_depth_correction_from_file(fn_bac_gc_depths,max_scale_factor,correct_range)
    GCp,ave_depths,GCp2_8,lowess_depth,correction,mu = gc_tuple
    return correction

def get_GC_depth_correction_from_file(fn_bac_gc_depths,max_scale_factor,correct_range=False):

    dt = np.dtype([('gc',np.float32),('depth',np.uint32),('n',np.uint32)])

    GC_width = int(open(fn_bac_gc_depths,'r').readline())
    bac_gc_depths = np.loadtxt(fn_bac_gc_depths,dtype=dt,skiprows=1)
    
    #print correct_range
    gc_tuple = get_GC_depth_correction_vect(bac_gc_depths['depth'],bac_gc_depths['n'],GC_width,max_scale_factor,correct_range=correct_range)
    
    GCp,ave_depths,GCp2_8,lowess_depth,correction,mu = gc_tuple
    return gc_tuple

def get_GC_depth_correction_vect(sum_depths,n_bases,GC_width,max_correction_factor,calc_correction_factor=True,correct_range=False):

    #max_correction_factor = (max_correction_factor==-1) and 999999 or max_correction_factor
    #print "getting GC correction factor: max scale factor %d"%(max_correction_factor)


    assert(sum_depths.shape[0]==2*GC_width+1+1) #have to count 0 as well~ so, 41+1

    frac_bases = n_bases.astype(np.float64)/n_bases.sum()

    ave_depths = sum_depths/n_bases.astype(np.float64)
    GCp = np.arange(0,2*GC_width+1+1)/float((2*GC_width+1))

    ave_depths[np.where(np.isnan(ave_depths))] = 0

    #now, chop off the nans, make the array only as wide as the non-nans    
    ave_depths2_8 = ave_depths[np.where(np.logical_and(GCp>=0.25,GCp<=.75))]
    GCp2_8 = GCp[np.where(np.logical_and(GCp>=0.25,GCp<=.75))]
    GCp0_2 = GCp[np.where(GCp<0.25)]
    GCp8_10 = GCp[np.where(GCp>.75)]

    #now apply lowess on this
    lowess_depth = biostats.lowess(GCp2_8,ave_depths2_8,f=.15)
    #lowess_depth = biostats.lowess(GCp2_8,ave_depths2_8,f=.1) #USING .1... not working... wy?
    #lowess_depth = biostats.lowess(GCp2_8,ave_depths2_8,f=.25)


    line_func = lambda p, x: (p[0]*x+p[1])

    k=5
    y1 = lowess_depth[0:k+1]
    x1 = GCp2_8[0:k+1]

    l = lowess_depth.shape[0]
    y2 = lowess_depth[l-k:l+1]
    x2 = GCp2_8[l-k:l+1]

    print(lowess_depth)
    print("fit lowess on",ave_depths2_8)
    print(GCp0_2)
    print(GCp8_10)
    print(x1,y2)
    print(x2,y2)

    p1 = get_line_params(x1,y1)
    p2 = get_line_params(x2,y2)



    left_line = line_func(p1,GCp0_2)
    right_line = line_func(p2,GCp8_10) 

    lowess_depth = np.r_[left_line,lowess_depth,right_line]                    

    #lowess_depth[np.where(lowess_depth<=0)=np.min(lowess_depth[np.where(lowess_depth>0)])
    mu = sum_depths.astype(np.float64).sum()/n_bases.astype(np.float64).sum()
    #correction = lowess_depth - mu
    lowess_depth = np.clip(lowess_depth,1e-10,1e30)
    correction = mu/lowess_depth

    if(correct_range):
        #print "correcting in range .75"
        GC_max = 0.75
        GC_min = 0.2
        iGC_max = (np.where(GCp>GC_max))[0][0]
        iGC_min = (np.where(GCp<GC_min))[0][0]
        max_correction_factor = correction[iGC_max]
        correction=np.clip(correction,1.0/max_correction_factor,max_correction_factor)
    elif(max_correction_factor<=0):
        correction=np.ones(correction.shape[0])
    else:
        correction=np.clip(correction,1.0/max_correction_factor,max_correction_factor)

    return GCp,ave_depths,GCp2_8,lowess_depth,correction,mu
    #return GCp,ave_depths,y_interpolate,correction,mu





def get_lib_pass_thresh(thresh):
    #return 0.85
    if thresh!= None:
        return thresh
    return 0.85

def get_pop(input_pop,verbose):
    if(input_pop.upper() in ["YORUBA","YRI_1"]): return "Yoruba"    
    if(input_pop.upper() in ["KOREAN","JAPANESE","HAN-CHINESE","HAN"]): return "Asian"    
    if(input_pop.upper() in ["CEPH","CEU_1"]): return "European"    

    if verbose:    
        print("ethnicity %s ?????"%(input_pop))
    return input_pop
    sys.exit(1)

    return "NONE"

class genome_info:

    def get_sex(self,name):

        name = name.split('.')[0]

        if name.find("sampled") != -1:
            name = name.split("-")[0]        

        return self.sex_tags[name]    
            
    
    class genome:

        def __init__(self,name,sex,pop,mapping_dir,bac_dir,seq_dir,primary_analysis_dir):
            self.genome_name = name
            self.indiv_name = name.split(".")[0]
            self.sex = sex
            self.pop = pop
            self.passed_qc = False
            self.coverage = -1
            self.mu_corr=-1

            self.mapping_dir = mapping_dir
            self.bac_dir = bac_dir
            self.seq_dir = seq_dir
            self.primary_analysis_dir = primary_analysis_dir

    def __init__(self,in_genomes,sex_ethnicity_index,QC_check=False,verbose=False,ignore_1kg=False,alternate_pop_column=None):
        self.genomes_list = open(in_genomes,'r').readlines()
        self.sex_tags = {} 
        self.pop_tags = {} 
        self.pops = []
        self.genomes_by_pop = {}
        
        for sex_pop_line in open(sex_ethnicity_index,'r').readlines():
            (indiv,sex,pop) = sex_pop_line.split()[0:3] 
            if alternate_pop_column!=None:
                pop = sex_pop_line.split()[alternate_pop_column] 
            
            pop = (get_pop(pop,verbose)).upper()
            self.sex_tags[indiv.split(".")[0]] = sex
            self.pop_tags[indiv.split(".")[0]] = pop
            if not pop in self.pops:
                self.pops.append(pop)
                self.genomes_by_pop[pop] = {}
        
        
        self.genomes = {}
    
        if(QC_check):
            if verbose:
                print("loading QC data...")
        for genome_info in self.genomes_list:
            (genome_name,mapping_dir,bac_dir,seq_dir,primary_analysis_dir) = genome_info.split()
            if genome_name[0] == "#": continue
            indiv_name = genome_name.split(".")[0]
        
            ################THIS IS A CHECKER TO SEE IF IT"S A 'special' subsampled genome
            pop_sex_tag_name = indiv_name
            if indiv_name.find("sampled") != -1:
                #indiv_name = indiv_name.split("-")[0]        
                pop_sex_tag_name = indiv_name.split("-")[0]        
            

            new_genome = self.genome(genome_name,self.sex_tags[pop_sex_tag_name],self.pop_tags[pop_sex_tag_name],mapping_dir,bac_dir,seq_dir,primary_analysis_dir)
            self.genomes[genome_name] = new_genome
            self.genomes_by_pop[self.pop_tags[genome_name]][genome_name] = new_genome
            
            if(QC_check):
                total_cvg = 0
                mu_corr = 0
                ls_libs = os.listdir("%s/%s"%(bac_dir,indiv_name))
                
                lib_count=0
                for lib in ls_libs:
                    if(lib[0] == "_" or lib[0]=="."):continue
                    libdir = "%s/%s/%s"%(bac_dir,indiv_name,lib)
                    
                    if(lib_pass_qc(libdir)):
                        new_genome.passed_qc = True
                        total_cvg+=get_cvg(libdir)
                        lib_count+=1
                        mu_corr+=get_corr(libdir)

                new_genome.coverage=total_cvg
                
                if lib_count !=0:
                    new_genome.mu_corr = mu_corr/lib_count
            else:
                new_genome.passed_qc = True
        if(QC_check):
            if verbose:
                print("done")
            
def get_cvg(dir):
    
    sumF = open("%s/%s"%(dir,"summary_stats_dim0.txt"),"r")
    line = sumF.readline()            
    line = sumF.readline()            
    if(not(line)):
        return 0
    cvg = float(line.split(",")[1])
    return cvg


def get_corr(dir):
    sumF = open("%s/%s"%(dir,"summary_stats_dim0.txt"),"r")
    line = sumF.readline()            
    line = sumF.readline()            
    if(not(line)):
        return 0

    sline = line.split(",")
    ed = sline[0]
    ec = float(sline[1])
    corr = float(sline[2])
    return corr    

def lib_pass_qc(dir,thresh=None):
    sumF = open("%s/%s"%(dir,"summary_stats_dim0.txt"),"r")
    line = sumF.readline()            
    line = sumF.readline()            
    if(not(line)):
        return 0

    sline = line.split(",")
    ed = sline[0]
    ec = float(sline[1])
    corr = float(sline[2])
    corr_no_x = float(sline[3])
    
    if(corr > get_lib_pass_thresh(thresh)):
        return 1
    else:
        return 0

def file_exists(ls,file):
    for f in ls:
        if(f==file.name):
            return 1
    return 0


def get_sexes(sex_file):
    sexes = {}
    fsex = open(sex_file,"r")
    line=fsex.readline()

    while(line):
        (sample_id,sex) = line.split()[0:2]
        if(sample_id in sexes):
            if sexes[sample_id] != sex:
                print("ERROR - sample in sex file twice")
                print(sample_id)
                assert(0)
        else:
            sexes[sample_id] = sex

        line=fsex.readline()
    return sexes

class kg_file_handler:

    def __init__(self,index_file,assert_lib_insert_in,sex_file,analyze_pilot,output_dir,analyze_platform = "ILLUMINA"):



        self.output_dir = output_dir 
        self.genomes = {}
        self.libraries= {}
        self.pilot_dir = ""
        
        sex_hash = get_sexes(sex_file)

        #global assert_lib_insert
        assert_lib_insert = assert_lib_insert_in
        seqfile = open(index_file,"r")

        self.pilot_dir = "Pilot"+analyze_pilot.split()[3]
        line = seqfile.readline()
        line = seqfile.readline()



        for line in seqfile.readlines():
            sline = line.rstrip().split("\t")
            #(FASTQ_FILE,MD5,    RUN_ID,    STUDY_ID,    STUDY_NAME,    CENTER_NAME,SUBMISSION_ID,    SUBMISSION_DATE,SAMPLE_ID,SAMPLE_NAME,POPULATION,    EXPERIMENT_ID,    INSTRUMENT_PLATFORM,    INSTRUMENT_MODEL,    LIBRARY_NAME,    RUN_NAME,    RUN_BLOCK_NAME,    INSERT_SIZE,    LIBRARY_LAYOUT,    PAIRED_FASTQ,    WITHDRAWN,    WITHDRAWN_DATE,    COMMENT,    READ_COUNT,    BASE_COUNT,ANALYSIS_GROUP) = sline;
            (FASTQ_FILE,MD5,    RUN_ID,    STUDY_ID,    STUDY_NAME,    CENTER_NAME,SUBMISSION_ID,    SUBMISSION_DATE,SAMPLE_ID,SAMPLE_NAME,POPULATION,    EXPERIMENT_ID,    INSTRUMENT_PLATFORM,    INSTRUMENT_MODEL,    LIBRARY_NAME,    RUN_NAME,    RUN_BLOCK_NAME,    INSERT_SIZE,    LIBRARY_LAYOUT,    PAIRED_FASTQ,    WITHDRAWN,    WITHDRAWN_DATE,    COMMENT,    READ_COUNT,    BASE_COUNT) = sline;


            ###ONLYU LOOK AT ILLUMINA!
            if not STUDY_NAME == analyze_pilot:
                continue
        
            if(not(INSTRUMENT_PLATFORM==analyze_platform)):
                line = seqfile.readline()
                continue

            
            #MAKE A NEW FILE OBJECT
            if(WITHDRAWN == "1"):
                newFile = file(FASTQ_FILE,MD5,PAIRED_FASTQ,LIBRARY_LAYOUT,self.pilot_dir,1)
            else:
                newFile = file(FASTQ_FILE,MD5,PAIRED_FASTQ,LIBRARY_LAYOUT,self.pilot_dir,0)    

            #See if this genome exists...    
            if(not(SAMPLE_NAME in self.genomes)):
                if(not(SAMPLE_NAME in sex_hash)):
                    print("error no sex for %s" %(SAMPLE_NAME))
                    print(sline) 
                    assert(0)
                sex = sex_hash[SAMPLE_NAME]
                new_genome = genome(SAMPLE_NAME,POPULATION,sex)
                self.genomes[SAMPLE_NAME] = new_genome;    
            
            #ASSIGN THAT FILE TO A LIBRARY -> IF THAT LIBRARY does not exist, 
            #then, make sure you add it to the correct genome
            LIBRARY_NAME = LIBRARY_NAME.split(".")[0]
            if(not(LIBRARY_NAME in self.libraries)):
                new_lib = library(LIBRARY_NAME,LIBRARY_LAYOUT,INSERT_SIZE,CENTER_NAME)
                self.libraries[LIBRARY_NAME] = new_lib;
                self.genomes[SAMPLE_NAME].addlibrary(new_lib);    
                self.libraries[LIBRARY_NAME].addfile(newFile)
            else:
                self.libraries[LIBRARY_NAME].update_paired(LIBRARY_LAYOUT,INSERT_SIZE)
                self.libraries[LIBRARY_NAME].addfile(newFile)

            #if(FASTQ_FILE == "data/NA18547/sequence_read/ERR000285_1.recal.fastq.gz"):
            #    print libraries[LIBRARY_NAME].insert_size
            #    print libraries[LIBRARY_NAME].paired
            #    for file in libraries[LIBRARY_NAME].files:
            #        print file
            #    sys.exit(1)

    def check_file_sum(self,file,sample):    
        args = "%s//%s/%s/%s" %(self.output_dir,self.pilot_dir,sample,file.name)
        print(args)
        sum = subprocess.Popen(["md5sum", args],stdout=subprocess.PIPE).communicate()[0]
        sum = sum.split()[0]
        print("calculated md5sum:" + sum)
        print("index file md5sum:" + file.md5sum)
        if(sum == file.md5sum):
            print("---md5sum match---")
            return 1
        else:
            print("---md5sum FAIL!---")
            print(args)
            return 0 

    def generate_dir_struct(self):
        
        for key in self.genomes:
            command = "mkdir  %s/%s/%s" % (self.output_dir,self.pilot_dir,key)
            print(command)
            os.system(command)


    def check_sums(self):
        print("analyzing checksums...")
        ret = 0
        print(self.genomes)
        for key in self.genomes:
            ls = os.listdir("%s/%s/%s" %(self.output_dir,self.pilot_dir,key))    
            for lib in self.genomes[key].libraries:
                for file in lib.files:
                    if(not(file_exists(ls,file)) and file.withdrawn==0):
                        print("file missing! " + file.path)
                    elif(file.withdrawn==0):
                        if(not(self.check_file_sum(file,key))):
                            ret +=1
                            print("CHECKSUM FAILED! " + file.path)

        return ret

    def clean_feeze(self,name):

        global genomes
        for key in self.genomes:
            command = "rm -r %s/%s/%s/%s" %(self.output_dir,self.pilot_dir,key, name)
            print(command)

    def make_sim_link_freeze(self):
        
        if(check_files_exist()):
            print("\n\n\n------------\nWARNING! FILES MISSING\n------------\n\n\ncontinuing\n")
        #check that there are no files to be removed
        if(check_no_withdrawn()):
            print("\n\n\n------------\nWARNING! WITHDRAWN_FILES_EXIST\n------------\n\n\ncontinuing\n")
        
        today = datetime.today()
        ts = time.time()
        freezeName = "%d-%d-%d_lf_%f" %(today.year,today.month,today.day,ts)

        print("generating local freeze:" , freezeName)
        count = 0
        for key in self.genomes:
            ls = os.listdir("%s/%s/%s"%(self.output_dir,self.pilot_dir,key))    
            command = "mkdir %s/%s/%s/%s" %(self.output_dir,self.pilot_dir,key, freezeName)
            print(command)
            os.system(command)
            currdir = "%s/%s/%s" %(self.output_dir,self.pilot_dir,key)
            for lib in self.genomes[key].libraries:
                seCounter = 0
                peCounter = 0
                peHash = {}  #keep track of paired end indices, so both paired end files have a isngle index
                for file in lib.files:
                    if(file.withdrawn == 0):
                    #    print file.name
                        filelen = file.get_read_len()
                        if(file.paired == "PAIRED" and file.mate_pair != ""):
                            lr = int(file.name.split("_")[1][0]) 
                            if(lib.insert_size==-1): 
                                ins_size = "UNK"
                            else:
                                ins_size = str(lib.insert_size) 

                            if(file.mate_pair in peHash):
                                peHash[file.name] = pindex
                            else:
                                peHash[file.path] = peCounter
                                pindex = peCounter
                                peCounter+=1

                            filename = "pe%d-ins%s_%s_%s.%d.fastq.gz"%(lr,ins_size,filelen,lib.name,pindex)
                        else:
                            filename = "se_%s_%s.%d.fastq.gz"%(filelen,lib.name,seCounter)
                            seCounter+=1
                            
                        command = "ln -s %s/%s %s/%s/%s"%(currdir,file.name,currdir,freezeName,filename)
                        print(command)
                        os.system(command)
                            #count +=1
                            #if(count%100 == 0):
                            #print count


    def rm_withdrawn_files(self):
        for key in self.genomes:
            ls = os.listdir("%s/%s/%s"%(self.output_dir,self.pilot_dir,key))    
            for lib in self.genomes[key].libraries:
                for file in lib.files:
                    if((file_exists(ls,file)) and file.withdrawn==1):
                        command = "rm self.output_dir/%s/%s/%s" %(self.pilot_dir,key, file.name) 
                        #print "MISSING %s" % (file.name)
                        print(command)    

    def check_no_withdrawn(self):
        ret = 0
        for key in self.genomes:
            ls = os.listdir("%s/%s/%s" %(self.output_dir,self.pilot_dir, key))    
            for lib in self.genomes[key].libraries:
                for file in lib.files:
                    if((file_exists(ls,file)) and file.withdrawn==1):
                        print("withdrawn file exists! " + file.path)
                        ret +=1
        return ret;

    def check_files_exist(self):
        ret = 0
        for key in self.genomes:
            ls = os.listdir("%s/%s/%s" %(self.output_dir,self.pilot_dir, key))    
            for lib in self.genomes[key].libraries:
                for file in lib.files:
                    if(not(file_exists(ls,file)) and file.withdrawn==0):
                        print("file missing! " + file.path)
                        ret +=1
        return ret;

    def generate_download_script(self,aspera_speed=90):
        for key in self.genomes:
            pilot_dir = mkdir("%s"%(self.output_dir),self.pilot_dir)
            genome_dir = mkdir(pilot_dir,"%s"%(key))
            #ls = os.listdir("%s/%s/%s" %(self.output_dir,self.pilot_dir, key))    
            ls = os.listdir(genome_dir)    
            for lib in self.genomes[key].libraries:
                for file in lib.files:
                    if(not(file_exists(ls,file)) and file.withdrawn==0):
                        command = "ascp -i ~/EEE_Lab/asperaetc/asperaweb_id_dsa.putty -l %dM anonftp@ftp-trace.ncbi.nlm.nih.gov:/1000genomes/ftp/pilot_data/%s %s/%s/%s" % (aspera_speed,file.path,self.output_dir,self.pilot_dir,key)
                        #print "MISSING %s" % (file.name)
                        print(command)    

                    
    def check_single_end_reads(self):
        
        for key in self.genomes:
            for lib in self.genomes[key].libraries:
                if(lib.paired == "SINGLE" or lib.paired == "BOTH"):
                    ls = os.listdir("%s/%s/%s" %(self.output_dir,self.pilot_dir, key))    
                    #for file in lib.files:
                        #if(not(file_exists(ls,file)) and file.withdrawn==0):
                        #    print "missing -> ",key, lib.name, file.name, lib.paired
                            #print "\t\tEXISTS " , str(file.withdrawn)
                        #if((file_exists(ls,file)) and file.withdrawn==1):
                        #    print "delete me -> ",key, lib.name, file.name, lib.paired
                        #if((file_exists(ls,file)) and file.withdrawn==0):
                        #    print "ok"
                        #else:
                        #    print "deleted"
                            

    #def get_pilot(self,index_file):
        
    #    seqfile = open(index_file,"r")
    #    line = seqfile.readline()
    #    line = seqfile.readline()
    #    sline = line.split("\t")
    #    study_id = sline[4]
    #    sstudy_id = study_id.split()
    #    pilot = sstudy_id[2] + sstudy_id[3]
    #    return pilot

genomes = {}
libraries = {}
pilot_dir = ""
assert_lib_insert=1


