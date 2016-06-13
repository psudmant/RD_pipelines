import sys
import socket
import os
import os.path
from optparse import OptionParser
#from collections import defaultdict
import time
import gzip
import socket
import re

import math

import numpy as np
import numpy.ma as ma

import tables

try:
    _cur_hostname = socket.gethostname()
except Error(e):
    sys.stderr.write( 'Error occurred getting hostname' )
    
def debug_output( msg, indent=6 ):
    print '%s{WSSD %s} [%s] %s'%( ''.join(' '*indent), _cur_hostname, time.ctime(), msg )

class DenseTrackSet(object):
    
    def __init__( self,
                  fnContigLengths,
                  fnWssd,
                  overwrite,
                  openMode,
                  groupsToCheck=[],
                  compression=False ):
        
        self.compress = compression
        
        assert os.path.exists(fnContigLengths)
        if openMode=='r':
            assert not overwrite
            assert os.path.exists(fnWssd), fnWssd
        
        debug_output('WssdBase: reading contig lengths from file %s'%fnContigLengths)        
        
        self.mContigNameLen = {}
        for l in open(fnContigLengths,'r'):
            l=l.replace('\n','').split('\t')
            self.mContigNameLen[l[0]]=int(l[1])
        
        debug_output('WSSD space: %d contigs totaling %d bp'%( len(self.mContigNameLen), sum(self.mContigNameLen.values()) ))
        
        if overwrite or not os.path.exists(fnWssd): 
            self.tbl = tables.openFile( fnWssd, 'w' )
        else:
            if openMode=='r':
                self.tbl = tables.openFile( fnWssd, 'r' )
            else:
                self.tbl = tables.openFile( fnWssd, 'a' )
    
    class DenseTrackSetGroup(object):
        def __init__(self, dst, h5grpname):
            self.dst=dst
            self.h5grpname = h5grpname
            
        def __getitem__( self, contigName ):
            return self.dst.tbl.getNode( self.dst.tbl.getNode('/%s'%self.h5grpname) , contigName )
        
        def __contains__(self, contigName):
            return contigName in self.dst.tbl.getNode('/%s'%self.h5grpname)
        
        def addArray( self, dtype, lncols ):            
            for contigName in self.dst.mContigNameLen:
                debug_output('setting up contig %s'%contigName)
#                self.dst.tbl.createArray(
#                                 self.dst.tbl.getNode('/%s'%self.h5grpname), 
#                                 contigName,
#                                 np.zeros( tuple( [self.dst.mContigNameLen[contigName]] + lncols ), 'uint16' )) 
                
                ar=self.dst.tbl.createCArray(
                                 self.dst.tbl.getNode('/%s'%self.h5grpname), 
                                 contigName,
                                 dtype, 
                                 tuple( [self.dst.mContigNameLen[contigName]] + lncols ), 
                                 '', 
                                 filters=self.dst.compress and tables.Filters(complevel=1,complib='lzo') or None)        

                chunkStart, chunkEnd = 0, 0
                while chunkStart < self.dst.mContigNameLen[ contigName ]:
                    chunkEnd = min( chunkStart+int(50e6), self.dst.mContigNameLen[ contigName ] - 1 )
                        
                    if len(lncols)==3:
                        ar[chunkStart:(chunkEnd+1),:,:,:]=0
                    elif len(lncols)==2:
                        ar[chunkStart:(chunkEnd+1),:,:]=0
                    elif len(lncols)==1:
                        ar[chunkStart:(chunkEnd+1),:]=0
                    elif len(lncols)==0:
                        ar[chunkStart:(chunkEnd+1)]=0
                        
                    chunkStart = chunkEnd + 1
            
    def addGroup( self, grpName ):
        assert not grpName in self, grpName
        self.tbl.createGroup( self.tbl.root, grpName )
        return self[grpName]
        
    def __contains__(self, grpName):
        return '/%s'%grpName in self.tbl.root
    
    def __getitem__(self, grpName):
        return DenseTrackSet.DenseTrackSetGroup(self, grpName)
    
    def __del__(self):
        self.tbl.close()
    

def sam_get_chrom(l):
    return l[2]

repatCigarMatchOrIns = re.compile('[0-9]+[MI]+')
# assumes presence of an NM tag.    
# this is probably screaming to be cython'd.
def sam_line_to_start_stop_ed(l):
    lTags=l[11].split(':')
    i=0
    while lTags[i*3]!='NM':
        i+=1
    
    refLen = sum([int(x[:-1]) for x in repatCigarMatchOrIns.findall( l[5] )])
    refCozStart=int(l[3])-1
        
    return (l[2],
            refCozStart,
            refCozStart + refLen - 1,
            int(lTags[i*3+2]) )
        
#        
#class WssdTest:
#    def __init__(self, fnContigLengths, fnWssd, overwrite, maxEditDist, openMode):
#        
#        debug_output('WssdBase: reading contig lengths from file %s'%fnContigLengths)        
#        
#        self.mContigNameLen = {}
#        for l in open(fnContigLengths,'r'):
#            l=l.replace('\n','').split('\t')
#            self.mContigNameLen[l[0]]=int(l[1])
#        
#        self.lContigs=sorted(self.mContigNameLen.keys())
#        self.lofsContigStarts=np.array( [0]+[self.mContigNameLen[contigName] for contigName in self.lContigs][:-1], 'int32' ).cumsum()
#        self.mContigNameOfsStart = dict(zip(self.lContigs,self.lofsContigStarts))
#        self.totalLen=sum(self.mContigNameLen.values())
#        
#        debug_output('WSSD space: %d contigs totaling %d bp'%( len(self.mContigNameLen), sum(self.mContigNameLen.values()) ))
#        
#        self.tbl = np.memmap( fnWssd, 'int16', openMode, shape=( self.totalLen, maxEditDist+1, 2 ), order='C')
#    
#    
#    def populateFromSAMFiles(self, 
#                             lFns,
#                             maxEditDist=2,
#                             areCompressed=False,
#                             ignoreHitsToUnknownContigs=True,
#                            ):
#        
#        nLinesProcd, nHitsCounted = 0, 0
#        
#        isLastChromInMem=False
#        lastChrom, curChrom = None, None
#        curChromBaseOfs = None
#
#        tbl=self.tbl
#        
#        for fnPlacement in lFns:            
#            if areCompressed:            
#                filIn = gzip.open(fnPlacement,'r')
#            else:   
#                filIn = open(fnPlacement,'r')
#                
#            debug_output('working on %s'%fnPlacement)
#                
#            l=filIn.readline().replace('\n','').rstrip()
#            if len(l)==0: continue
#            
#            while len(l)>0:
#                l=l.split('\t')
#                
#                curChrom, curCoStart, curCoStop, curEditDist = sam_line_to_start_stop_ed(l)
#                
#                edDist=0
#            
#                # if we read out the entire last chromosome
#                # (and that chromosome is changing), write it back
#                if lastChrom != curChrom:
#                    curChromBaseOfs = self.mContigNameOfsStart[curChrom]
#                    
#                nLinesProcd += 1
#                if nLinesProcd % 10000 == 0 :
#                    debug_output('processed %d lines'%nLinesProcd)
#
#                if curEditDist>maxEditDist:  
#                    l=filIn.readline().replace('\n','').rstrip()
#                    continue
#                
#                tbl[ (curChromBaseOfs+curCoStart-1):(curChromBaseOfs+curCoStop), curEditDist, 0 ] += 1  # record read depth
#                tbl[ (curChromBaseOfs+curCoStart-1), curEditDist, 1 ] += 1              # record read start
#                nHitsCounted += 1
#                
#                l=filIn.readline().replace('\n','').rstrip()
#            
#                lastChrom = curChrom        
#            
#        return nHitsCounted, nLinesProcd
        

class WssdFile(DenseTrackSet):
    
    class DepthGetter:
        def __init__(self, wssd):
            self.wssd=wssd
        def __contains__(self, tracksetBaseName):
            return self.wssd.__contains__('depthAndStarts_%s'%tracksetBaseName)
        def __getitem__(self, tracksetBaseName):
            return self.wssd['depthAndStarts_%s'%tracksetBaseName]
    
    def __init__( self,
                  fnContigLengths,
                  fnWssd, 
                  overwrite,
                  openMode,
                  tracksetsToCheck=[],
                  compression=False,
                  datatype=tables.UInt16Atom() ):

        DenseTrackSet.__init__( 
                           self,    
                           fnContigLengths,
                           fnWssd,
                           overwrite,
                           openMode,
                           groupsToCheck=tracksetsToCheck,
                           compression=compression )
        
        self.datatype = datatype
        self.depthgetter = WssdFile.DepthGetter(self)
        
    def addTrackSet( self, tracksetBaseName, maxEditDist=2 ):
        self.addGroup( 'depthAndStarts_%s'%tracksetBaseName )
        # N x (E+1) x 2 
        # N=length of chrom, E=max edit dist, (0:depth, 1: starts)
        self[ 'depthAndStarts_%s'%tracksetBaseName ].addArray( self.datatype , [maxEditDist+1,2] )  
        
    def __getattribute__(self,k):
        if k=='depth':
            return self.depthgetter
        else:
            return DenseTrackSet.__getattribute__(self,k)
        
    def __contains__(self,k):
        return DenseTrackSet.__contains__(self, 'depthAndStarts_%s'%k) or \
               DenseTrackSet.__contains__(self, k)
        
    ########################################################################
    
    # requirement: 
    # each of the SAM files in the list lFns must be sorted with respect to chromosome
    #  (though not necessarily coordinate)

    # this is quite a bit more memory-efficient
    # requires reads all short and ~same length
    def populateFromSAMFilesChunkify(self, 
                             lFns,
                             tracksetBaseName,
                             maxEditDist,
                             areCompressed=False,
                             runFullyInMem=False,
                             ignoreHitsToUnknownContigs=True,
                             placementSizeThresh=1e9,
                            ):
        
        # TODO: impl runFullyInMem
        assert not runFullyInMem
        
        nLinesProcd, nHitsCounted = 0, 0
        
        isLastChromInMem=False
        lastChrom, curChrom = None, None
        curArray = None
        
        assert len(lFns)<=256,'time to fix the file handle bug jacob!'
        
        lFiles = areCompressed and [gzip.open(fn,'r') for fn in lFns ] or [open(fn,'r') for fn in lFns]
        lNextLine = [ f.readline().replace('\n','').rstrip() for f in lFiles ]
        lFilesDone = [ len(l)==0 for l in lNextLine ]
        lChromNextLine = [ (not lFilesDone[i]) and sam_get_chrom(lNextLine[i].split('\t')) or None for i in xrange(len(lFiles)) ]
                
        while not all(lFilesDone): 
            # find the next chromosome to process
            chromToProcess = min([c for c in lChromNextLine if c!=None])
            debug_output('%d files not done, %s is next chromosome to process'%( len([k for k in lFilesDone if k]), chromToProcess ))
            liFilesToProcess = [ i for i in xrange(len(lFiles)) if lChromNextLine[i]==chromToProcess ]
            # if we encounter a chrom/contig name that we don't know about then skip it in all files.
            while not all(lFilesDone) and not chromToProcess in self.depth[tracksetBaseName]: 
                assert ignoreHitsToUnknownContigs, chromToProcess
                for ifile in liFilesToProcess:
                    filCur=lFiles[ifile]
                    l=lNextLine[ifile]
                    lNextLine[ifile]=None
                    curChrom = chromToProcess
                    while len(l)>0 and curChrom==chromToProcess:
                        l=filCur.readline().replace('\n','').rstrip()
                        curChrom = len(l)>0 and sam_get_chrom(l.split('\t')) or None
                    lNextLine[ifile]=l
                    lChromNextLine[ifile]=curChrom
                    lFilesDone[ifile]=(curChrom==None)
                    
                if len([c for c in lChromNextLine if c!=None])>0:
                    chromToProcess = min([c for c in lChromNextLine if c!=None])
                    liFilesToProcess = [ i for i in xrange(len(lFiles)) if lChromNextLine[i]==chromToProcess ]
                    debug_output('%d files not done, %s is next chromosome to process'%( len([k for k in lFilesDone if k]), chromToProcess ))
                else:
                    liFilesToProcess = []
                 
            if all(lFilesDone):break
                                
            inextRec=0    
            lCurStarts,lLens,lCurEdDist=np.zeros((100000,),'uint32'), np.zeros((100000,),'uint16'), np.zeros((100000,),'uint8')
                
            for ifile in liFilesToProcess:
                debug_output('working on %s chrom %s'%(lFns[ifile],chromToProcess))
                filCur=lFiles[ifile]
                l=lNextLine[ifile]
                lNextLine[ifile]=None
                curChrom=chromToProcess
                while len(l)>0 and curChrom==chromToProcess:
                    l=l.split('\t')                    
                    curChrom, curCozStart, curCozStop, curEditDist = sam_line_to_start_stop_ed(l)
                                        
                    nLinesProcd += 1
                    if nLinesProcd % 10000 == 0 :
                        debug_output('processed %d lines (inextrec: %d)'%(nLinesProcd,inextRec))
        
                    if curEditDist>maxEditDist:  
                        l=filCur.readline().replace('\n','').rstrip()
                        curChrom = len(l)>0 and sam_get_chrom(l.split('\t')) or None
                        continue
                 
                    if inextRec%100000==0:
                        if inextRec==lCurStarts.shape[0]:
                            lCurStarts.resize( (inextRec+100000,) )
                            lCurStarts[inextRec:]=0
                            lLens.resize( (inextRec+100000,) )
                            lLens[inextRec:]=0
                            lCurEdDist.resize( (inextRec+100000,) )
                            lCurEdDist[inextRec:]=0
                        
                    lCurStarts[inextRec]=curCozStart
                    lLens[inextRec]=curCozStop-curCozStart+1
                    lCurEdDist[inextRec]=curEditDist
                    
                    inextRec+=1
                            
                    nHitsCounted += 1
                    
                    l=filCur.readline().replace('\n','').rstrip()
                    curChrom = len(l)>0 and sam_get_chrom(l.split('\t')) or None
                    
                lNextLine[ifile]=l
                lChromNextLine[ifile]=curChrom
                lFilesDone[ifile]=(curChrom==None)
                print lFilesDone
                
            if inextRec > 0:
                # TODO put these back but for some reason ipython debug does not like them
                #lCurStarts.resize((inextRec,))
                #lLens.resize((inextRec,))
                #lCurEdDist.resize((inextRec,))
                lCurStarts=lCurStarts[:inextRec]
                lLens=lLens[:inextRec]
                lCurEdDist=lCurEdDist[:inextRec]
                
                ipermCurHits=np.argsort(lCurStarts[:inextRec])
                lCurStarts = lCurStarts[ipermCurHits]
                lLens = lLens[ipermCurHits]
                lCurEdDist = lCurEdDist[ipermCurHits]
                                
                maxLen = lLens.max()
                                
                # populate the depth matrix a chunk at a time to avoid filling RAM
                curChromLen = self.depth[tracksetBaseName][ chromToProcess ].shape[0]
                nextChunkStart = 0
                iChunk = 0
                while nextChunkStart < self.depth[tracksetBaseName][ chromToProcess ].shape[0]:
                    nextChunkEnd = min( curChromLen - 1, 
                                        int(nextChunkStart + 10e6) - 1 )
                    debug_output('reading in chrom %s (len %d) chunk %d (%d-%d)'%(
                        chromToProcess, curChromLen, iChunk, nextChunkStart, nextChunkEnd))
                    curArray = self.depth[tracksetBaseName][ chromToProcess ][nextChunkStart:(nextChunkEnd+1),:,:]
                    curArray = curArray.astype('uint32')
                    
                    # 3 cases; bulk of the action happens in (1)
                    # 1. Reads starting AND ending within this window, given tightly-bounded max read length
                    ics_l = lCurStarts.searchsorted( nextChunkStart, side='left' )
                    ics_r = lCurStarts.searchsorted( max(nextChunkStart,nextChunkEnd + 1 - maxLen), side='right' )
                    # if there are any reads _starting_ within this window:
                    # [nextChunkStart,nextChunkEnd+1-maxLen]
                    if (ics_r - 1) - (ics_l) + 1  > 0 :
                        for ics in xrange(ics_l, ics_r):
                            curCozStart = lCurStarts[ ics ]
                            curEditDist = lCurEdDist[ ics ]
                            
                            rngroOfs = (curCozStart - nextChunkStart,
                                        curCozStart+lLens[ics] - nextChunkStart )
                            
                            assert curCozStart >= nextChunkStart, (lCurStarts[ics],lCurEdDist[ics],lLens[ics])
                            assert curCozStart+lLens[ics]-nextChunkStart-1 <= nextChunkEnd,  (lCurStarts[ics],lCurEdDist[ics],lLens[ics])
                            
                            curArray[ rngroOfs[0]:rngroOfs[1], curEditDist, 0 ] += 1  # record read depth
                            curArray[ rngroOfs[0], curEditDist, 1 ] += 1              # record read start
                            
                    # 2. Reads starting before and potentially ending within this window given tightly-buonded max read length
                    # ie starts in 
                    # [nextChunkStart-maxLen+1,nextChunkStart-1]
                    if nextChunkStart>0:
                        ics_l = lCurStarts.searchsorted( max(0,nextChunkStart - maxLen + 1), side='left' )
                        ics_r = lCurStarts.searchsorted( nextChunkStart - 1, side='right' )
                        if (ics_r - 1) - (ics_l) + 1  > 0 :
                            for ics in xrange(ics_l, ics_r):
                                curCozStart = lCurStarts[ ics ]
                                curEditDist = lCurEdDist[ ics ]
                                
                                rngroOfs = (max(0,curCozStart - nextChunkStart),
                                            min(curCozStart + lLens[ics], nextChunkEnd) - nextChunkStart)

                                if rngroOfs>=0: # if this read does not overlap with this window, don't count it
                                    curArray[ rngroOfs[0]:rngroOfs[1], curEditDist, 0 ] += 1
                                    # don't need to record the start, b/c before the window.
                    
                    # 3. Reads within and potentially ending outside this window given tightly-buonded max read length
                    # ie starts in 
                    # [nextChunkEnd-maxLen, nextChunkEnd]
                    ics_l = lCurStarts.searchsorted( max(nextChunkStart,nextChunkEnd + 2 - maxLen), side='left' )
                    ics_r = lCurStarts.searchsorted( nextChunkEnd, side='right' )
                    if (ics_r - 1) - (ics_l) + 1 > 0:
                        for ics in xrange(ics_l, ics_r):
                            curCozStart = lCurStarts[ ics ]
                            curEditDist = lCurEdDist[ ics ]
                                                        
                            rngroOfs = (curCozStart - nextChunkStart,
                                        min(curCozStart + lLens[ics], nextChunkEnd) - nextChunkStart )
                                    
                            assert min(rngroOfs)>0
                                                        
                            curArray[ rngroOfs[0]:rngroOfs[1] , curEditDist, 0 ] += 1
                            curArray[ rngroOfs[0], curEditDist, 1 ] += 1
                                                
                    debug_output('writing out chrom %s (len %d) chunk %d (%d-%d)'%(
                        chromToProcess, curChromLen, iChunk, nextChunkStart, nextChunkEnd))
                        
                    np.clip(curArray,0,2**16-1,curArray)
                    curArray=curArray.astype('uint16')                        
                    self.depth[tracksetBaseName][ chromToProcess ][nextChunkStart:(nextChunkEnd+1),:,:] = curArray
                    
                    nextChunkStart = nextChunkEnd + 1
                
                    del curArray
                
            del lCurStarts
            del lLens
            del lCurEdDist
        
        for f in lFiles: f.close()        
        
        return nHitsCounted, nLinesProcd
    
    # needs more memory
    # works with an entire chrom at a time.
    def populateFromSAMFiles(self, 
                             lFns,
                             tracksetBaseName,
                             maxEditDist,
                             areCompressed=False,
                             runFullyInMem=False,
                             ignoreHitsToUnknownContigs=True,
                             placementSizeThresh=1e9,
                            ):
        
        # TODO: impl runFullyInMem
        assert not runFullyInMem
        
        nLinesProcd, nHitsCounted = 0, 0
        
        isLastChromInMem=False
        lastChrom, curChrom = None, None
        curArray = None
        
        assert len( lFns ) <= 256 , 'time to fix the file handle bug jacob!' 
        
        lFiles = areCompressed and [gzip.open(fn,'r') for fn in lFns ] or [open(fn,'r') for fn in lFns]
        lNextLine = [ f.readline().replace('\n','').rstrip() for f in lFiles ]
        lFilesDone = [ len(l)==0 for l in lNextLine ]
        lChromNextLine = [ (not lFilesDone[i]) and sam_get_chrom(lNextLine[i].split('\t')) or None for i in xrange(len(lFiles)) ]
                
        while not all(lFilesDone): 
            # find the next chromosome to process
            chromToProcess = min([c for c in lChromNextLine if c!=None])
            debug_output('%d files not done, %s is next chromosome to process'%( len([k for k in lFilesDone if k]), chromToProcess ))
            liFilesToProcess = [ i for i in xrange(len(lFiles)) if lChromNextLine[i]==chromToProcess ]
            # if we encounter a chrom/contig name that we don't know about then skip it in all files.
            while not all(lFilesDone) and not chromToProcess in self.depth[tracksetBaseName]: 
                assert ignoreHitsToUnknownContigs, chromToProcess
                for ifile in liFilesToProcess:
                    filCur=lFiles[ifile]
                    l=lNextLine[ifile]
                    lNextLine[ifile]=None
                    curChrom = chromToProcess
                    while len(l)>0 and curChrom==chromToProcess:
                        l=filCur.readline().replace('\n','').rstrip()
                        curChrom = len(l)>0 and sam_get_chrom(l.split('\t')) or None
                    lNextLine[ifile]=l
                    lChromNextLine[ifile]=curChrom
                    lFilesDone[ifile]=(curChrom==None)

                if len([c for c in lChromNextLine if c!=None])>0:
                    chromToProcess = min([c for c in lChromNextLine if c!=None])
                    liFilesToProcess = [ i for i in xrange(len(lFiles)) if lChromNextLine[i]==chromToProcess ]
                    debug_output('%d files not done, %s is next chromosome to process'%( len([k for k in lFilesDone if k]), chromToProcess ))
                else:
                    liFilesToProcess = []
                
            if all(lFilesDone):break
                
            debug_output('reading in full chrom %s'% chromToProcess)
            curArray = self.depth[tracksetBaseName][ chromToProcess ][:,:,:].astype('int32')
                
            for ifile in liFilesToProcess:
                debug_output('working on %s chrom %s'%(lFns[ifile],chromToProcess))
                filCur=lFiles[ifile]
                l=lNextLine[ifile]
                lNextLine[ifile]=None
                curChrom=chromToProcess
                while len(l)>0 and curChrom==chromToProcess:
                    l=l.split('\t')                    
                    curChrom, curCozStart, curCozStop, curEditDist = sam_line_to_start_stop_ed(l)
                    
                    nLinesProcd += 1
                    if nLinesProcd % 10000 == 0 :
                        debug_output('processed %d lines'%nLinesProcd)
        
                    if curEditDist>maxEditDist:  
                        l=filCur.readline().replace('\n','').rstrip()
                        curChrom = len(l)>0 and sam_get_chrom(l.split('\t')) or None
                        continue
                    
                    curArray[ (curCozStart):(curCozStop+1), curEditDist, 0 ] += 1  # record read depth
                    curArray[ curCozStart, curEditDist, 1 ] += 1              # record read start
        
                    nHitsCounted += 1
                    
                    l=filCur.readline().replace('\n','').rstrip()
                    curChrom = len(l)>0 and sam_get_chrom(l.split('\t')) or None
                    
                lNextLine[ifile]=l
                lChromNextLine[ifile]=curChrom
                lFilesDone[ifile]=(curChrom==None)
                
            debug_output('writing out full chrom %s'% chromToProcess)
            np.clip(curArray,0,2**16-1,curArray)
            curArray=curArray.astype('uint16')
            self.depth[tracksetBaseName][ chromToProcess ][:,:,:] = curArray
                
        for f in lFiles: f.close()        
        
        return nHitsCounted, nLinesProcd    

    def populateFromSAMFilesOneAtATime(self, 
                             lFns,
                             tracksetBaseName,
                             maxEditDist,
                             areCompressed=False,
                             runFullyInMem=False,
                             ignoreHitsToUnknownContigs=True,
                             placementSizeThresh=1e9,
                            ):
        
        # TODO: impl runFullyInMem
        assert not runFullyInMem
        
        nLinesProcd, nHitsCounted = 0, 0
        
        isLastChromInMem=False
        lastChrom, curChrom = None, None
        curArray = None
                
        for fnPlacement in lFns:            
            if areCompressed:            
                filIn = gzip.open(fnPlacement,'r')
            else:   
                filIn = open(fnPlacement,'r')
                
            debug_output('working on %s'%fnPlacement)
                
            l=filIn.readline().replace('\n','').rstrip()
            if len(l)==0: continue
            
            while len(l)>0:
                l=l.split('\t')
                
                curChrom, curCozStart, curCozStop, curEditDist = sam_line_to_start_stop_ed(l)
                
                edDist=0
            
                # if we read out the entire last chromosome
                # (and that chromosome is changing), write it back
                if lastChrom != curChrom:
                    if isLastChromInMem:
                        debug_output('writeback full chrom %s'% curChrom)
                        self.depth[tracksetBaseName][ lastChrom ][:,:,:] = curArray
                        curArray = None
                        isLastChromInMem = False
                        
                    if not curChrom in self.depth[tracksetBaseName]:
                        assert ignoreHitsToUnknownContigs, curChrom
                        curArray = None
                        l=filIn.readline().replace('\n','').rstrip()
                        lastChrom = curChrom
                        continue

                    if os.stat( fnPlacement ).st_size > placementSizeThresh:
                        isLastChromInMem = True
                        curArray = self.depth[tracksetBaseName][ curChrom ][:,:,:]
                        debug_output('read in full chrom %s'% curChrom)
                    else:
                        isLastChromInMem = False
                        curArray = self.depth[tracksetBaseName][ curChrom ]
                    
                nLinesProcd += 1
                if nLinesProcd % 10000 == 0 :
                    debug_output('processed %d lines'%nLinesProcd)

                if curEditDist>maxEditDist:  
                    l=filIn.readline().replace('\n','').rstrip()
                    continue
                
                curArray[ curCozStart:(curCozStop+1), curEditDist, 0 ] += 1  # record read depth
                curArray[ curCozStart, curEditDist, 1 ] += 1              # record read start

                nHitsCounted += 1
                
                l=filIn.readline().replace('\n','').rstrip()
            
                lastChrom = curChrom
    
        if isLastChromInMem:
            debug_output('writeback full chrom %s'% lastChrom)
            self.depth[tracksetBaseName][lastChrom][:,:,:] = curArray
            curArray = None
            isLastChromInMem = False        
            
            
        return nHitsCounted, nLinesProcd
    

    ########################################################################
    
    # coordinates are 1-based and inclusive
    # provide a list of edit distances which should be summed together lEdsToSum
    # returns an array of depths, and an array of starts 
    def getDepthsAndStartsPerBase(self, 
                                  tracksetBaseName, 
                                  chrom, 
                                  cooStart, 
                                  cooStop,
                                  lEdsToSum=[0]):
        assert cooStart>=1 and cooStart<=self.mContigNameLen[chrom], (chrom, cooStart, cooStop)
        assert cooStop>=1 and cooStop<=self.mContigNameLen[chrom], (chrom, cooStart, cooStop)
        
        ard=self.depth[tracksetBaseName][chrom][(cooStart-1):cooStop,:,:]
        ard=ard[:,lEdsToSum,:]
        ard=ard.sum(1)
        
        retDepths = ard[:,0]
        retStarts = ard[:,1]
        return retDepths, retStarts
 
# Average base-by-base vector (eg coverage)
# using equally-sized/spaced windows.   
#
# IN v: full base-by-base vector
# IN cooStart: [1, start coordinate
# IN cooStop: ,1] stop coordinate
# IN wndWidth: int, >1 of averaging window width
# IN kernel: 'sum' or 'mean'
# OUT (pos,      #  array of N positions ([1,)for window starts
#      v)        #  array of N values for windows
def perbaseToWindowAverage( V, 
                            cooStart,
                            cooStop,
                            wndWidth,
                            kernel='sum',  # can be sum, or mean
                           ):
    if kernel=='sum':
        krn = np.ones( (wndWidth,), 'int32' )
    elif kernel=='mean':
        krn = (1.0/float(wndWidth))*np.ones((wndWidth,), 'float64')
        
    Vconv = np.convolve( V, krn, 'full' )
    Vconv = Vconv[ (wndWidth-1):V.shape[0]:wndWidth ]
    
    outStarts = cooStart+np.arange(Vconv.shape[0])*wndWidth
    
    return outStarts,Vconv

def perbaseToWindowAverageMasked( 
                            V, 
                            arMasked,
                            cooStart,
                            cooStop,
                            wndWidth,
                            kernel='sum',  # can be sum, or mean
                           ):
    assert arMasked.shape==V.shape, (arMasked.shape,V.shape)
    
    
    krn = np.ones( (wndWidth,), 'int32' )
        
    Vconv = np.convolve( V*((1-arMasked).astype('uint64')), krn, 'full' )
    Vconv = Vconv[ (wndWidth-1):Vconv.shape[0]:wndWidth ]
        
    assert Vconv.shape[0]==(V.shape[0]/wndWidth + min(V.shape[0]%wndWidth,1))
        
    nBpMasked = np.convolve( arMasked, np.ones(wndWidth,'int32'), 'full' )
    nBpMasked = nBpMasked[ (wndWidth-1):nBpMasked.shape[0]:wndWidth ]
    
    assert nBpMasked.shape[0]==Vconv.shape[0]
    
    outStarts = cooStart+np.arange(Vconv.shape[0])*wndWidth
    
    nBpUnmasked = wndWidth - nBpMasked

    if kernel=='mean':
        nBpUnmasked[-1] = cooStop - outStarts[-1] + 1
        
        Vconv = (ma.masked_array( Vconv, nBpUnmasked == 0 ) / nBpUnmasked.astype('float64')).filled(0.)
        
    return outStarts,Vconv,nBpUnmasked


# find boundaries for tiling windows with a constant number of NONmasked bases inside
# which therefore can potentially differ in width

# -> [cooStart_i, cooStop_i], i.e., Nx2.  Coords are [1,1].

def findBoundariesForEqualWeightedWindows( maskVec,
                                           cooStart,
                                           cooStop,
                                           nbpUnmaskedPerWnd  ):
    maskReg=~(maskVec[cooStart-1:cooStop])
    csMaskReg=maskReg.astype('int32').cumsum()
    
    nTotalUnmasked = csMaskReg[-1]
    
    assert nTotalUnmasked>0, nTotalUnmasked
    
    igapdWndStarts = np.arange(0,nTotalUnmasked,nbpUnmaskedPerWnd)
    
    irngungWnd = csMaskReg.searchsorted(igapdWndStarts,'left')
    irngungWnd = np.c_[ irngungWnd, np.r_[ irngungWnd[1:]-1, csMaskReg.shape[0]-1 ] ]
    
    corngWnd = cooStart + irngungWnd
    
    return corngWnd
    
def perbaseToWindowAverageMaskedVariableBoundaries( 
                            V, 
                            arMasked,
                            cooStart,  # overall bounds start, inclusive, [11]
                            cooStop,   # overall bounds end, inclusive, [11]
                            corngWnd,  # window boundaries, inclusive, [11]
                            kernel='sum',  # can be sum, or mean
                           ):
    
    assert arMasked.shape==V.shape, (arMasked.shape,V.shape)
    
    Vwndavg=np.zeros( (corngWnd.shape[0],), 'float64' )
    nBpUnmasked=np.zeros( (corngWnd.shape[0],), 'int32' )
    
    ofsrngWnd = np.c_[ corngWnd[:,0] - 1 - (cooStart - 1), corngWnd[:,1] - (cooStart - 1) ]
    
    for iwnd in xrange(corngWnd.shape[0]):
        nBpUnmasked[iwnd] = (ofsrngWnd[iwnd,1]-ofsrngWnd[iwnd,0]+1) - arMasked[ ofsrngWnd[iwnd,0]:ofsrngWnd[iwnd,1]+1 ].sum()
        Vwndavg[iwnd] = ( V[ ofsrngWnd[iwnd,0]:ofsrngWnd[iwnd,1]+1 ] * (~arMasked[ofsrngWnd[iwnd,0]:ofsrngWnd[iwnd,1]+1]).astype('f') ).sum()

    if kernel=='mean':
        Vwndavg = (ma.masked_array(Vwndavg,nBpUnmasked==0)/nBpUnmasked.astype('f')).filled(0.)
        
    return Vwndavg, nBpUnmasked