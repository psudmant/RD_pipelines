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
except Exception, e:
    sys.stderr.write( 'Error occurred getting hostname' )
    
try:
    import pysam
except Exception, e:
    sys.stderr.write('pysam import error')
    
def debug_output( msg, indent=6 ):
    print >>sys.stderr, '%s{WSSD %s} [%s] %s'%( ''.join(' '*indent), _cur_hostname, time.ctime(), msg )

###########################################################################

def resizeIfNotNone(A, growBy, inextRec):
    if A!=None:
        newSize=A.shape[0]+growBy
        A.resize( newSize, refcheck=False )
        A[inextRec:] = 0
        
###########################################################################
        
def sam_get_chrom(l):
    return l[2]

repatCigarMatchOrIns = re.compile('[0-9]+[MI]+')
# assumes presence of an NM tag.    
# this is probably screaming to be cython'd.

# commenting out the version of this designed to parse busted-ass
# SAM files from mr/mrs fast.  instead changed to meet SAM spec as far as 
# addl tags being each in a tab-delim'd field
#def sam_line_to_start_stop_ed(l):
#    lTags=l[11].split(':')
#    i=0
#    while lTags[i*3]!='NM':
#        i+=1
#    
#    refLen = sum([int(x[:-1]) for x in repatCigarMatchOrIns.findall( l[5] )])
#    refCozStart=int(l[3])-1
#        
#    return (l[2],
#            refCozStart,
#            refCozStart + refLen - 1,
#            int(lTags[i*3+2]) )
    
def sam_line_to_start_stop_ed(l):
    i=11
    lTags=l[i].split(':')
    while lTags[0]!='NM':
        i+=1
        lTags=l[i].split(':')
    
    refLen = sum([int(x[:-1]) for x in repatCigarMatchOrIns.findall( l[5] )])
    refCozStart=int(l[3])-1
        
    return (l[2],
            refCozStart,
            refCozStart + refLen - 1,
            int(lTags[2]) )

###########################################################################

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
        
        filContigLengths = open(fnContigLengths,'r')
        l=filContigLengths.readline()
        if l[0]=='@':
            while len(l)>0:
                l=l.replace('\n','').split('\t')
                if l[0]=='@SQ':
                    m={}
                    for kv in l[1:]: m[kv.split(':')[0]]=kv.split(':')[1]
                    if 'SN' in m and 'LN' in m:
                        self.mContigNameLen[m['SN']]=int(m['LN'])
                l=filContigLengths.readline()            
        else:
            while len(l)>0:
                l=l.replace('\n','').split('\t')
                self.mContigNameLen[l[0]]=int(l[1])
                l=filContigLengths.readline()
        
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
        #if alt_grpName!=None: grpName = alt_grpName[grpName]
        return DenseTrackSet.DenseTrackSetGroup(self, grpName)
    
    def __del__(self):
        self.tbl.close()
    
    ########### populate-related methods
        
    # arrays inside mFieldLVals must be sorted coming in.  Fields are added (clipped, etc)
    def handleCurrentHitsInChunks(self,
                                   mFieldLVals,  # starts and ends are [00] - 0-based, inclusive.  
                                   chromChunk,
                                   corngChunk,  # [00] - 0-based, inclusive
                                   tracksetBaseName,
                                   **kwargs):

        # get the subset of hits that might fall inside the target region.
        
        nAllHits = mFieldLVals['starts'].shape[0]
        maxHitLen = mFieldLVals['lens'].max()
        
        # now find hits that overlap this window
        
                    
        # look at lens - what's the longest sequence?
        mFieldLValsCurChunk = {}
                    
        lStartsClipped=np.clip(mFieldLVals['starts'],corngChunk[0],corngChunk[1])
        lEndsClipped=np.clip(mFieldLVals['starts']+mFieldLVals['lens']-1,corngChunk[0],corngChunk[1])
        lStartsInside = (mFieldLVals['starts'] == lStartsClipped)
        lEndsInside = ((mFieldLVals['starts']+mFieldLVals['lens']-1) == lEndsClipped)

#        if corngChunk[0]==60000000 and chromChunk=='chr19':
#            from IPython.Shell import IPShellEmbed; ipshell = IPShellEmbed([]); ipshell()

        lAreInChunk = np.where( lStartsInside | lEndsInside )[0]

        if lAreInChunk.shape[0]==0: return
        
        for k in mFieldLVals:
            mFieldLValsCurChunk[k] = mFieldLVals[k][ lAreInChunk ]
        
        mFieldLValsCurChunk['startsClipped'] = lStartsClipped[ lAreInChunk ]
        mFieldLValsCurChunk['endsClipped'] = lEndsClipped[ lAreInChunk ]
        mFieldLValsCurChunk['startsInside'] = lStartsInside[ lAreInChunk ]
        mFieldLValsCurChunk['endsInside'] = lEndsInside[ lAreInChunk ]
                
        debug_output('handling %d hits in %s:%d-%d'%( lAreInChunk.shape[0], chromChunk, corngChunk[0], corngChunk[1] ))

        # only read in/write out as little as necessary
        corngChunkClipped = ( mFieldLValsCurChunk['startsClipped'].min(), mFieldLValsCurChunk['endsClipped'].max() )

        #from IPython.Shell import IPShellEmbed; ipshell = IPShellEmbed([]); ipshell()

        debug_output('enterring')  #jdbg

        self.populateCurChunkHits( mFieldLValsCurChunk,
                                   chromChunk,
                                   corngChunkClipped,
                                   tracksetBaseName,
                                   **kwargs )
        
        debug_output('returning2')  #jdbg
        

    def fieldsRequiredToPopulate(self):
        return set( ['starts','lens','eds'] )

    # expect: starts, lengths, startsInside, endsInside, startsClipped, endsClipped,
    #   possibly other fields as requested by self.fieldsRequiredToPopulate
    #
    # override this method in subclasses to populate differently
    #  (i.e., based on mapping uniqueness, edit distance, etc)
    def populateCurChunkHits(self,
                             mFieldLVals,
                             chromChunk,
                             corngChunk, 
                             tracksetBaseName,
                             **kwargs ): 

        curChromLen = self.mContigNameLen[chromChunk]

        debug_output('reading in chrom %s (len %d): %d-%d'%(
            chromChunk, curChromLen, corngChunk[0], corngChunk[1]))

        curArray = self[tracksetBaseName][ chromChunk ][corngChunk[0]:(corngChunk[1]+1)]
        curArray = curArray.astype('uint32')

        debug_output('tallying coverage in %s:%d-%d'%(chromChunk, corngChunk[0], corngChunk[1]))

        lStartsClipped,lEndsClipped=mFieldLVals['startsClipped'],mFieldLVals['endsClipped']
        for ir in xrange(lStartsClipped.shape[0]):
            curArray[ lStartsClipped[ir]:lEndsClipped[ir]+1 ] += 1
        
        np.clip(curArray,0,2**16-1,curArray)
        curArray=curArray.astype('uint16')                        
        self[tracksetBaseName][ chromChunk ][ corngChunk[0] :(corngChunk[1]+1)] = curArray

        debug_output('writing out chrom %s (len %d) %d-%d'%(
            chromChunk, curChromLen, corngChunk[0], corngChunk[1]))
            
        debug_output('returning')  #jdbg
            
        del curArray
    
    def populateFromSAMFilesChunkify(self, 
                                     lFns,
                                     tracksetBaseName,
                                     maxEditDist,
                                     areCompressed=False,
                                     ignoreHitsToUnknownContigs=True,
                                     placementSizeThresh=1e9,
                                     **kwargs
                                    ):
        
        if isinstance(self, WssdFile):
            fxnGetGroup=lambda gn:self.depth[gn]
        else:
            fxnGetGroup=lambda gn:self[gn]
        
        nLinesProcd, nHitsCounted = 0, 0
        
        isLastChromInMem=False
        lastChrom, curChrom = None, None
        curArray = None
        
        assert len(lFns)<=256,'time to fix the file handle bug jacob!'
        
        lFiles = areCompressed and [gzip.open(fn,'r') for fn in lFns ] or [open(fn,'r') for fn in lFns]        
        lNextLine = [ f.readline().rstrip() for f in lFiles ]
        lFilesDone = [ len(l)==0 for l in lNextLine ]
        
        # skip past SAM headers
        while any( [len(lNextLine[ifile])>0 and lNextLine[ifile][0]=='@' for ifile in xrange(len(lFns))] ):
            for ifile in xrange(len(lFns)):
                if len(lNextLine[ifile])==0:
                    lFilesDone[ifile]=True
                elif lNextLine[ifile][0]=='@': 
                    lNextLine[ifile]=lFiles[ifile].readline().rstrip() 
        
        lChromNextLine = [ (not lFilesDone[i]) and sam_get_chrom(lNextLine[i].split('\t')) or None for i in xrange(len(lFiles)) ]
        
        while not all(lFilesDone): 
            # find the next chromosome to process
            chromToProcess = min([c for c in lChromNextLine if c!=None])
            debug_output('%d files not done, %s is next chromosome to process'%( len([k for k in lFilesDone if k]), chromToProcess ))
            liFilesToProcess = [ i for i in xrange(len(lFiles)) if lChromNextLine[i]==chromToProcess ]
            # if we encounter a chrom/contig name that we don't know about then skip it in all files.
            while not all(lFilesDone) and not chromToProcess in fxnGetGroup(tracksetBaseName): 
                assert ignoreHitsToUnknownContigs, chromToProcess
                for ifile in liFilesToProcess:
                    filCur=lFiles[ifile]
                    l=lNextLine[ifile]
                    lNextLine[ifile]=None
                    curChrom = chromToProcess
                    while len(l)>0 and curChrom==chromToProcess:
                        l=filCur.readline().rstrip()
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
                 
            if all(lFilesDone):break # 
                                
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
                        l=filCur.readline().rstrip()
                        curChrom = len(l)>0 and sam_get_chrom(l.split('\t')) or None
                        continue
                 
                    if inextRec%100000==0:
                        if inextRec==lCurStarts.shape[0]:
                            resizeIfNotNone(lCurStarts, 100000, inextRec)
                            resizeIfNotNone(lLens, 100000, inextRec)
                            resizeIfNotNone(lCurEdDist, 100000, inextRec)
                        
                    lCurStarts[inextRec]=curCozStart
                    lLens[inextRec]=curCozStop-curCozStart+1
                    lCurEdDist[inextRec]=curEditDist
                    
                    inextRec+=1
                            
                    nHitsCounted += 1
                    
                    l=filCur.readline().rstrip()
                    curChrom = len(l)>0 and sam_get_chrom(l.split('\t')) or None
                    
                lNextLine[ifile]=l
                lChromNextLine[ifile]=curChrom
                lFilesDone[ifile]=(curChrom==None)
                
            if inextRec > 0:
                lCurStarts=lCurStarts[:inextRec]
                lLens=lLens[:inextRec]
                lCurEdDist=lCurEdDist[:inextRec]
                
                ipermCurStarts = np.argsort( lCurStarts )
                lCurStarts = lCurStarts[ipermCurStarts]
                lLens = lLens[ipermCurStarts]
                lCurEdDist = lCurEdDist[ipermCurStarts]
                
                curChromLen = fxnGetGroup(tracksetBaseName)[ chromToProcess ].shape[0]
                cozNextChunkStart = 0
                iChunk = 0
                while cozNextChunkStart < self.mContigNameLen[ chromToProcess ]:
                    cozNextChunkEnd = min( curChromLen - 1, int(cozNextChunkStart + 10e6) - 1 )
                                                                        
                    self.handleCurrentHitsInChunks({'starts': lCurStarts,
                                                     'lens': lLens,
                                                     'eds': lCurEdDist},
                                                     chromToProcess, 
                                                     (cozNextChunkStart,cozNextChunkEnd),
                                                     tracksetBaseName,
                                                     **kwargs)
                    
                    cozNextChunkStart = cozNextChunkEnd + 1
                    
        return nHitsCounted, nLinesProcd
                    
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

    def populateCurChunkHits(self,
                             mFieldLVals,
                             chromChunk,
                             corngChunk, 
                             tracksetBaseName,
                             **kwargs ): 

        curChromLen = self.mContigNameLen[chromChunk]

        debug_output('reading in chrom %s (len %d): %d-%d'%(
            chromChunk, curChromLen, corngChunk[0], corngChunk[1]))

        curArray = np.nan_to_num(self.depth[tracksetBaseName][ chromChunk ][corngChunk[0]:(corngChunk[1]+1),:,:])
        curArray = curArray.astype('uint32')
        
        debug_output('tallying coverage in %s:%d-%d'%(chromChunk, corngChunk[0], corngChunk[1]))

        lofsStarts,lofsEnds=mFieldLVals['startsClipped']-corngChunk[0],mFieldLVals['endsClipped']-corngChunk[0]
        lStartsInChunk=mFieldLVals['startsInside']
        lEds = mFieldLVals['eds']
        for ir in xrange(lStartsInChunk.shape[0]):
            if ir%100000==0: debug_output('popul %d'%ir)  #jdbg
            curArray[ lofsStarts[ir]:lofsEnds[ir]+1, lEds[ir], 0 ] += 1
            if lStartsInChunk[ir]: 
                curArray[ lofsStarts[ir], lEds[ir], 1 ] += 1
        
        
        debug_output('tally1')  #jdbg
                
        np.clip(curArray,0,2**16-1,curArray)
        debug_output('tally2')  #jdbg
        curArray=curArray.astype('uint16')
        debug_output('tally3')  #jdbg                        
        self.depth[tracksetBaseName][ chromChunk ][ corngChunk[0] :(corngChunk[1]+1),:,:] = curArray
        debug_output('tally4')  #jdbg

        debug_output('writing out chrom %s (len %d) %d-%d'%(
            chromChunk, curChromLen, corngChunk[0], corngChunk[1]))
            
        del curArray
        
        debug_output('deld')  #jdbg

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
        
        ard=np.nan_to_num(self.depth[tracksetBaseName][chrom][(cooStart-1):cooStop,:,:]).astype(np.float64)
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
    
    igapdWndEnds = np.r_[ np.arange(0,nTotalUnmasked,nbpUnmaskedPerWnd)[1:], nTotalUnmasked ]
    
    irngungWnd = csMaskReg.searchsorted(igapdWndEnds,'right')
        
    irngungWnd = np.c_[ np.r_[0,irngungWnd[:-1]], irngungWnd-1 ]
    
    corngWnd = cooStart + irngungWnd
    
    #from IPython.Shell import IPShellEmbed; ipshell = IPShellEmbed([]); ipshell()
    
    return corngWnd

# find boundaries for tiling windows; try for a constant number of NONmasked bases 
# but don't allow the windows to exceed a given width

# --> [cooStart_i, cooStop_i], i.e., Nx2, coords are [1,1]
def findBoundariesForAboutEqualWeightedWindowsMaxWidth( 
                                           maskVec,
                                           cooStart,
                                           cooStop,
                                           nbpUnmaskedPerWnd,
                                           maxWndLen,
                                           supressEmptyWindows=True  ):
    maskReg=~(maskVec[cooStart-1:cooStop])
    csMaskReg=maskReg.astype('int32').cumsum()
    
    nTotalUnmasked = csMaskReg[-1]
    
    assert nTotalUnmasked>0, nTotalUnmasked
    
    lirngungWnd = []
    
    igapdNextWndStart = 0
    iungNextWndStart = 0
    
    while igapdNextWndStart < nTotalUnmasked:
        # look, up to maxWndLen bases ahead
        nUnmaskedCurCandWnd = csMaskReg[ min(iungNextWndStart + maxWndLen - 1, csMaskReg.shape[0]-1) ] - ( iungNextWndStart>0 and csMaskReg[ iungNextWndStart - 1 ] or 0 )
        
        # are there enough unmasked bases in this maximum searchable region
        if nUnmaskedCurCandWnd < nbpUnmaskedPerWnd :
            # no: just go to the end
            iungCurWndEnd = min(iungNextWndStart + maxWndLen - 1, csMaskReg.shape[0]-1)
            
            if nUnmaskedCurCandWnd > 0 or not supressEmptyWindows:
                # is the window completely empty?
                    lirngungWnd.append( [iungNextWndStart,iungCurWndEnd]  )
        
            iungNextWndStart = iungCurWndEnd + 1      
            
            igapdNextWndStart = csMaskReg[ iungNextWndStart - 1]      
        else:
            # yes: look for the break
            iungCurWndEnd = iungNextWndStart+(csMaskReg[ iungNextWndStart:min(iungNextWndStart + maxWndLen, csMaskReg.shape[0]) ] -
                                              csMaskReg[ iungNextWndStart ] ).searchsorted( nbpUnmaskedPerWnd, 'left' )
            
            lirngungWnd.append( [iungNextWndStart, iungCurWndEnd] )
            
            iungNextWndStart = iungCurWndEnd + 1
            
            igapdNextWndStart = csMaskReg[ iungNextWndStart - 1 ]
    
    lirngungWnd = np.array( lirngungWnd, 'int32' )
    
    corngWnd = cooStart + lirngungWnd
    
    #from IPython.Shell import IPShellEmbed; ipshell = IPShellEmbed([]); ipshell()
    
    return corngWnd

# find boundaries for tiling windows; try for a constant number of NONmasked bases 
# but don't allow the windows to exceed a given width
# 
# also enforce window boundaries at certain points (i.e., windows may not extend thru these points)

# --> [cooStart_i, cooStop_i], i.e., Nx2, coords are [1,1]
def findBoundariesForAboutEqualWeightedWindowsMaxWidthEnforcedBreaks( 
                                           maskVec,
                                           cooStart,
                                           cooStop,
                                           lcoBreaks,
                                           nbpUnmaskedPerWnd,
                                           maxWndLen,
                                           supressEmptyWindows=True  ):
    maskReg=~(maskVec[cooStart-1:cooStop])
    csMaskReg=maskReg.astype('int32').cumsum()
    
    nTotalUnmasked = csMaskReg[-1]
    
    if nTotalUnmasked==0:
        return np.zeros((0,2),'int32')
    
    lirngungWnd = []
    
    igapdNextWndStart = 0
    iungNextWndStart = 0
    
    liungBreaks = lcoBreaks - cooStart
    while len(liungBreaks)>0 and liungBreaks[0]==0:
        liungBreaks=liungBreaks[1:]
    while len(liungBreaks)>0 and liungBreaks[-1]==cooStop-cooStart+1:
        liungBreaks=liungBreaks[:-1]
    
    while igapdNextWndStart < nTotalUnmasked:
        
        # look up to the nearest of 
        #  (1). maxWndLen bp ahead,  (2). to the next of the lcoBreaks, or (3). to the end
        #
        # -> how far ahead _would_ we look if no lcoBreaks
        iungCurMaxCandWnd = min(iungNextWndStart + maxWndLen - 1, csMaskReg.shape[0]-1)
        
        if liungBreaks.shape[0]>0:
            # is there a break between the start of this window (iungNextWndStart) and there?
            iinextBreak_l=liungBreaks.searchsorted(iungNextWndStart,'left')
            iinextBreak_r=liungBreaks.searchsorted(iungCurMaxCandWnd,'right')
            if iinextBreak_l!=iinextBreak_r:
                # yes, so take the leftmost one
                iungCurMaxCandWnd =  max( 0, liungBreaks[ iinextBreak_l ] - 1  )
        
        # look, up to maxWndLen bases ahead
        nUnmaskedCurCandWnd = csMaskReg[ iungCurMaxCandWnd ] - ( iungNextWndStart>0 and csMaskReg[ iungNextWndStart - 1 ] or 0 )
        
        # are there enough unmasked bases in this maximum searchable region?

        iungCurWndEnd = iungNextWndStart+(csMaskReg[ iungNextWndStart:iungCurMaxCandWnd+1 ] -
                                          ( iungNextWndStart>0 and csMaskReg[ iungNextWndStart - 1 ] or 0)
                                           ).searchsorted( min(nUnmaskedCurCandWnd,nbpUnmaskedPerWnd), 'left' )

        if nUnmaskedCurCandWnd > 0 or not supressEmptyWindows:
            # is the window completely empty?
                lirngungWnd.append( [iungNextWndStart,iungCurWndEnd]  )
        
        iungNextWndStart = iungCurWndEnd + 1
        
        igapdNextWndStart = csMaskReg[ iungNextWndStart - 1 ]
    
    lirngungWnd = np.array( lirngungWnd, 'int32' )
    
    corngWnd = cooStart + lirngungWnd
    
    #from IPython.Shell import IPShellEmbed; ipshell = IPShellEmbed([]); ipshell()
    
    return corngWnd
    
# find boundaries for tiling windows; try for a constant number of NONmasked bases 
# but don't allow the windows to exceed a given width
# 
# joinLastWndsWithUnderNbpUnmsk: if the last window overall or the last window
#   before a break has fewer than this number of unmasked positions then join it with
#   the preceeding window (as long as there is one before the preceeding break)
# 
# also enforce window boundaries at certain points (i.e., windows may not extend thru these points)

# --> [cooStart_i, cooStop_i], i.e., Nx2, coords are [1,1]
def findBoundariesForAboutEqualWeightedWindowsMaxAndMinWidthEnforcedBreaks( 
                                           maskVec,
                                           cooStart,
                                           cooStop,
                                           lcoBreaks,
                                           minNbpUnmaskedPerWnd,
                                           maxWndLen,
                                           joinLastWndsWithUnderNbpUnmsk=100,
                                           supressEmptyWindowsFewerThanNbpUnmsk=1  ):
    maskReg=~(maskVec[cooStart-1:cooStop])
    csMaskReg=maskReg.astype('int32').cumsum()
    
    nTotalUnmasked = csMaskReg[-1]
    
    if nTotalUnmasked==0:
        return np.zeros((0,2),'int32')
    
    lirngungWnd = []
    
    igapdNextWndStart = 0
    iungNextWndStart = 0
    
    liungBreaks = lcoBreaks - cooStart
    while len(liungBreaks)>0 and liungBreaks[0]==0:
        liungBreaks=liungBreaks[1:]
    while len(liungBreaks)>0 and liungBreaks[-1]==cooStop-cooStart+1:
        liungBreaks=liungBreaks[:-1]
    
    while igapdNextWndStart < nTotalUnmasked:
        
        # look up to the nearest of 
        #  (1). maxWndLen bp ahead,  (2). to the next of the lcoBreaks, or (3). to the end
        #
        # -> how far ahead _would_ we look if no lcoBreaks
        iungCurMaxCandWnd = min(iungNextWndStart + maxWndLen - 1, csMaskReg.shape[0]-1)
        # is there a break between the start of this window (iungNextWndStart) and there?
        iinextBreak_l=liungBreaks.searchsorted(iungNextWndStart,'left')
        iinextBreak_r=liungBreaks.searchsorted(iungCurMaxCandWnd,'right')
        if iinextBreak_l!=iinextBreak_r:
            # yes, so take the leftmost one
            iungCurMaxCandWnd =  max( 0, liungBreaks[ iinextBreak_l ] - 1  )
        
        # look, up to maxWndLen bases ahead
        nUnmaskedCurCandWnd = csMaskReg[ iungCurMaxCandWnd ] - ( iungNextWndStart>0 and csMaskReg[ iungNextWndStart - 1 ] or 0 )
        
        # are there enough unmasked bases in this maximum searchable region?

        iungCurWndEnd = iungNextWndStart+(csMaskReg[ iungNextWndStart:iungCurMaxCandWnd+1 ] -
                                          ( iungNextWndStart>0 and csMaskReg[ iungNextWndStart - 1 ] or 0)
                                           ).searchsorted( min(nUnmaskedCurCandWnd,minNbpUnmaskedPerWnd), 'left' )

        if nUnmaskedCurCandWnd >= supressEmptyWindowsFewerThanNbpUnmsk:
            lirngungWnd.append( [iungNextWndStart,iungCurWndEnd]  )
        
        iungNextWndStart = iungCurWndEnd + 1
        
        igapdNextWndStart = csMaskReg[ iungNextWndStart - 1 ]
    
    lirngungWnd = np.array( lirngungWnd, 'int32' )
    
    lcoBreaksCheck=np.array( sorted(list(set( lcoBreaks.tolist()+[cooStop] ))) , 'int32' ) 
    _corngWnd=cooStart+lirngungWnd
    for jbrkCheck in xrange(lcoBreaksCheck.shape[0]):
        # find the first window starting after this break
        # while the next window is before the next break:   
        #    no - only one window here, continue
        # find the last window starting after this break and before (<=) the next break
        #   does it have <= joinLastWndsWithUnderNbpUnmsk bp unmasked?
        #   if so merge it 
        iwndFirstAfterCurBrk=np.searchsorted( _corngWnd[:,0], lcoBreaksCheck[jbrkCheck], 'left' )-1
        if jbrkCheck==lcoBreaksCheck.shape[0]-1:
            iwndLastBeforeNextBrk=_corngWnd.shape[0]-1
        else:
            iwndLastBeforeNextBrk=np.searchsorted( _corngWnd[:,0], lcoBreaksCheck[jbrkCheck+1], 'left' )-1
        assert iwndFirstAfterCurBrk > 0
        assert iwndLastBeforeNextBrk > 0 
        while iwndLastBeforeNextBrk>iwndFirstAfterCurBrk:
            sys.stderr.write('   jbrkCheck:%d (%d)  wnd0:%d-%d  wnd%d:%d-%d  iwndLastBeforeNextBrk:%d ,  iwndLastBeforeNextBrk:%d  / %d\n'%(
                        jbrkCheck, lcoBreaksCheck[jbrkCheck], _corngWnd[0,0], _corngWnd[0,1],
                        _corngWnd.shape[0]-1, _corngWnd[-1,0], _corngWnd[-1,1],
                        iwndLastBeforeNextBrk,iwndFirstAfterCurBrk,_corngWnd.shape[0]))
            
            corngWndConsideringMerging=_corngWnd[iwndLastBeforeNextBrk,:]
            nUnmBpInLastWnd=csMaskReg[ corngWndConsideringMerging[1]-cooStart ] - ( corngWndConsideringMerging[0]>cooStart and csMaskReg[ corngWndConsideringMerging[0]-cooStart ] or 0)
            if nUnmBpInLastWnd >= joinLastWndsWithUnderNbpUnmsk:
                break

            sys.stderr.write('   %d-%d : break %d %s (next: %s)  firstwnd %d/%d %d-%d %dbp (%d) lastwnd %d %d-%d %dbp\n'%(
                cooStart, cooStop, jbrkCheck, '%d'%lcoBreaksCheck[jbrkCheck],
                jbrkCheck==lcoBreaksCheck.shape[0]-1 and 'NA' or '%d'%lcoBreaksCheck[jbrkCheck+1],
                iwndFirstAfterCurBrk, _corngWnd.shape[0], 
                _corngWnd[iwndFirstAfterCurBrk,0], _corngWnd[iwndFirstAfterCurBrk,1],
                ((maskReg[_corngWnd[iwndFirstAfterCurBrk,0]-cooStart:_corngWnd[iwndFirstAfterCurBrk,1]-cooStart]).astype('int32')).sum(),
                nUnmBpInLastWnd,
                iwndLastBeforeNextBrk, _corngWnd[iwndLastBeforeNextBrk,0], _corngWnd[iwndLastBeforeNextBrk,1],
                ((maskReg[_corngWnd[iwndLastBeforeNextBrk,0]-cooStart:_corngWnd[iwndLastBeforeNextBrk,1]-cooStart+1]).astype('int32')).sum(),
                 ))

            #from IPython.Shell import IPShellEmbed; ipshell = IPShellEmbed([]); ipshell()

            # combine the last window with the previous one
            assert iwndLastBeforeNextBrk>1
            corngWndTemp=_corngWnd[:,:]

            _corngWnd=corngWndTemp[:iwndLastBeforeNextBrk,:]
            _corngWnd[-1,1]=corngWndTemp[iwndLastBeforeNextBrk,1]
            _corngWnd=np.r_[ _corngWnd, corngWndTemp[iwndLastBeforeNextBrk+1:,:] ]

            iwndFirstAfterCurBrk=np.searchsorted( _corngWnd[:,0], lcoBreaksCheck[jbrkCheck], 'left' ) - 1
            
            if jbrkCheck==lcoBreaksCheck.shape[0]-1:
                iwndLastBeforeNextBrk=_corngWnd.shape[0]-1
            else:
                iwndLastBeforeNextBrk=np.searchsorted( _corngWnd[:,0], lcoBreaksCheck[jbrkCheck+1], 'left' ) -1
            
            assert iwndFirstAfterCurBrk > 0
            assert iwndLastBeforeNextBrk > 0 

    
    corngWnd=_corngWnd
    
    #from IPython.Shell import IPShellEmbed; ipshell = IPShellEmbed([]); ipshell()
    
    return corngWnd
    
#def combineLastWindowIfBelowThresh( lcorngWnd,lcoBreaks, thresh ):
#    if abs(lcorngWnd[-1,1]-lcorngWnd[-1,0]+1)<thresh and lcorngWnd.shape[0]>1:
#        lcorngWnd[-2,1]=lcorngWnd[-1,1]
#        return lcorngWnd[:-1,:]
#    else:
#        return lcorngWnd

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

def extractPerBase(
    trkgrp,               # x[chrom][coord:coord].   (i.e., x=wssd.depths;  can use x.wssd from fxnDataAccessor to get back "up" to wssd object; or x.dst if densetrackset)...
    lchromQuery, 
    lcorngQuery,            # inclusive, [11]
    fxnDataAccessor = lambda t,chrom,cor:t[chrom][cor[0]-1:cor[1]], # get the data: f( trackset, chrom, corng[11] ) --> array size [N,]  , where N=corng[1]-corng[0]+1
    fxnProcessData = None,  # reproc the data in some way: f( trackset, chrom, corng[11], dataVector( size[m,] ) ) --> array size [m,]  : m<=N (b/c of masking)
    trkMask = None,
    maxChunkSize = 5000000,
    arTypeStr='float32' ):
    
    arVals=np.zeros((0,),arTypeStr)
    lcoVals=np.zeros((0,),'int32')
    lIdxRegionPerWnd = np.zeros((0,),'int32')
    
    for iq in xrange(len(lchromQuery)):
        chrom = lchromQuery[iq]
        
        corzWndCur = ( lcorngQuery[iq][0] - 1, 
                       min(lcorngQuery[iq][1]-1, lcorngQuery[iq][0] - 1 + maxChunkSize - 1 ) )
        
        while corzWndCur[0]<=lcorngQuery[iq][1]-1:
            
            if trkMask!=None:
                mask = trkMask[ chrom ][ corzWndCur[0]:corzWndCur[1]+1 ]

            #sys.stderr.write('getting %s:%d-%d\n'%(chrom,corzWndCur[0],corzWndCur[1]))
                            
            arWndCurAll = fxnDataAccessor( trkgrp, chrom, (corzWndCur[0]+1,corzWndCur[1]+1) )
            
            if trkMask!=None:
                liMask = np.where(~mask)[0]
                arWndCurUnmasked = arWndCurAll[ liMask ]
                
                lcoVals = np.r_[ lcoVals, corzWndCur[0]+liMask+1 ]
            else:
                arWndCurUnmasked = arWndCurAll
                
                lcoVals = np.r_[ lcoVals, np.arange(corzWndCur[0]+1,corzWndCur[1]+1+1) ]
                
            if fxnProcessData!=None:
                arWndProcd=fxnProcessData( trkgrp, chrom, (corzWndCur[0]+1,corzWndCur[1]+1), arWndCurUnmasked )
            else:
                arWndProcd = arWndCurUnmasked
        
            arVals = np.r_[ arVals, arWndProcd ]

            lIdxRegionPerWnd = np.r_[ lIdxRegionPerWnd, iq*np.ones(arWndProcd.shape[0],'int32') ]
            
            corzWndCur=( corzWndCur[1]+1,
                         min(corzWndCur[1]+1+maxChunkSize-1, lcorngQuery[iq][1] - 1 ) )
            

            
    return lcoVals, lIdxRegionPerWnd, arVals


#  lcorngWndBnds:  Lx2, where L is the number of windows; [11] coords
#  lIdxRegionPerWnd: L: is the index of the query region corresponding to each window
#  arVals: L, the per-window values.
def extractPerWnd( 
    trkgrp,      # x[chrom][coord:coord].   (i.e., x=wssd.depths;  can use x.wssd from fxnDataAccessor to get back "up" to wssd object; or x.dst if densetrackset)...
    lchromQuery,
    lcorngQuery,
    wndSize,
    wndVariable=True,
    fxnDataAccessor=lambda t,chrom,cor:t[chrom][cor[0]-1:cor[1]], # get the data: f( trackset, chrom, corng[11] ) --> array size [N,]  , where N=corng[1]-corng[0]+1,  
    fxnProcessDataPerWnd=lambda t,chrom,cor,V:V.mean(),  # reproc the data in some way: f( trackset, chrom, corng[11], dataVector( size[m,] ) ) --> one value  (m<=N bc of masking)
    trkMask=None,
    arTypeStr='float32' ):

    if wndVariable: assert trkMask!=None, 'must provide a mask if using variable-sized wnds...'
        
    arVals=np.zeros((0,),arTypeStr)
    lcorngWndBnds = np.zeros((0,2),'int32')
    lIdxRegionPerWnd = np.zeros((0,),'int32')
    
    for iq in xrange( len(lchromQuery) ):
        chrom = lchromQuery[iq]
        
        corngCur = lcorngQuery[iq]

        if trkMask!=None:
            mask = trkMask[ chrom ][ corngCur[0]-1:corngCur[1] ]
        else:
            mask = np.zeros( (corngCur[1]-corngCur[0]+1), 'bool' ) 
        
        if wndVariable:
            lcorngWndBndsCur = \
                findBoundariesForEqualWeightedWindows(trkMask[chrom], 
                                                  corngCur[0], 
                                                  corngCur[1], 
                                                  nbpUnmaskedPerWnd=wndSize )
            lcorngWndBndsCur -= 1  # [11] --> [00]
        else:
            lcorngWndBndsCur = np.zeros((0,2),'int32')
            
            corzWndCur = ( lcorngQuery[iq][0] - 1,
                           min(lcorngQuery[iq][1] - 1, lcorngQuery[iq][0] - 1 + wndSize - 1 ) )
            
            while corzWndCur[0]<lcorngQuery[iq][1]-1:
                lcorngWndBndsCur = np.r_[ lcorngWndBndsCur, np.array( [corzWndCur], 'int32' ) ]
                
                corzWndCur=( corzWndCur[1]+1,
                             min(corzWndCur[1]+1+wndSize-1, lcorngQuery[iq][1] - 1) )
                
        
        lcorngWndBnds = np.r_[ lcorngWndBnds, lcorngWndBndsCur ]
        lIdxRegionPerWnd = np.r_[ lIdxRegionPerWnd, iq*np.ones(lcorngWndBndsCur.shape[0],'int32') ]

        arValsWnds = np.zeros( (lcorngWndBnds.shape[0],), arTypeStr )

        # go thru windows fetching
        for iwnd in xrange( lcorngWndBnds.shape[0] ):
            
            corzWndCur = lcorngWndBnds[iwnd,:]
            
            arWndCurAll = fxnDataAccessor( trkgrp, chrom, (corzWndCur[0]+1,corzWndCur[1]+1) )

            if trkMask!=None:
                arWndCurUnmasked = arWndCurAll[ np.where(~(mask[corzWndCur[0]-(corngCur[0]-1):(corzWndCur[1]+1-(corngCur[0]-1))]) )[0] ]
            else:
                arWndCurUnmasked = arWndCurAll

            valWndProcd = fxnProcessDataPerWnd( trkgrp, chrom, (corzWndCur[0]+1,corzWndCur[1]+1), arWndCurUnmasked )
        
            arValsWnds[iwnd] = valWndProcd
                        
        arVals = np.r_[ arVals, arValsWnds ]
        
    return lcorngWndBnds+1, lIdxRegionPerWnd, arVals

#  lcorngWndBnds:  Lx2, where L is the number of windows; [11] coords
#  lIdxRegionPerWnd: L: is the index of the query region corresponding to each window
#  arVals: L, the per-window values.

# Here fxnProcessDataPerWnd is vectorized and takes a masked array
#   so it's Nxm, N = num windows, m = width each window, and mask
#   indicates if masking applies (and/or last window truncated)
def extractPerVarWndQuick( 
    trkgrp,      # x[chrom][coord:coord].   (i.e., x=wssd.depths;  can use x.wssd from fxnDataAccessor to get back "up" to wssd object; or x.dst if densetrackset)...
    lchromQuery,
    lcorngQuery,
    wndSize,
    wndVariable=False,
    fxnDataAccessor=lambda t,chrom,cor:t[chrom][cor[0]-1:cor[1]], # get the data: f( trackset, chrom, corng[11] ) --> array size [N,]  , where N=corng[1]-corng[0]+1,  
    fxnProcessWnds=lambda t,chrom,cor,V:V.mean(1),  # reproc the data in some way: f( trackset, chrom, corng[11], size[N,m] ) --> N values
    trkMask=None,
    arTypeStr='float32' ):

    assert wndVariable
        
    arVals=np.zeros((0,),arTypeStr)
    lcorngWndBnds = np.zeros((0,2),'int32')
    lIdxRegionPerWnd = np.zeros((0,),'int32')
    
    for iq in xrange( len(lchromQuery) ):
        chrom = lchromQuery[iq]
        
        corngCur = lcorngQuery[iq]

        if trkMask!=None:
            mask = trkMask[ chrom ][ corngCur[0]-1:corngCur[1] ]
            
            if (~mask).sum() == 0:
                continue
        else:
            mask = np.zeros( (corngCur[1]-corngCur[0]+1), 'bool' ) 
        
        lcorngWndBndsCur = \
            findBoundariesForEqualWeightedWindows(trkMask[chrom], 
                                              corngCur[0], 
                                              corngCur[1], 
                                              nbpUnmaskedPerWnd=wndSize )
        lcorngWndBndsCur -= 1  # [11] --> [00]
        
            
        lcorngWndBnds = np.r_[ lcorngWndBnds, lcorngWndBndsCur ]
        lIdxRegionPerWnd = np.r_[ lIdxRegionPerWnd, iq*np.ones(lcorngWndBndsCur.shape[0],'int32') ]

        arValsWnds = np.zeros( (lcorngWndBnds.shape[0],), arTypeStr )
        
        arValsAll = fxnDataAccessor( trkgrp, chrom, (lcorngWndBndsCur[0,0]+1, lcorngWndBndsCur[-1,1]+1) )
        liUM = np.where(~mask)[0]
        liiUMWndStart = np.arange( 0, liUM.shape[0], wndSize )
        liiUMWndEnd = np.r_[ liiUMWndStart[:-1] - 1, liUM.shape[0]-1 ]
        #np.clip( liiUMWndEnd, 0, liUM.shape[0]-1, liiUMWndEnd  )
        
        assert np.all(liUM[liiUMWndStart][1:] == (lcorngWndBndsCur[:,0]-(corngCur[0]-1))[1:])
        
        arValsAllUM = arValsAll[ liUM ]
        nValsCur=arValsAllUM.shape[0]
        arValsAll = np.resize( arValsAllUM, (nValsCur/wndSize+min(1,nValsCur%wndSize), wndSize) )

        arValsAllButLastWnd = arValsAll[:-1,:]
        arValsAllButLastWnd = ma.masked_array( arValsAllButLastWnd, mask=False )
        arValsAllLastWnd = arValsAll[-1,:]
        
        arValsAllLWM = np.zeros( (wndSize,), 'bool' )
        if nValsCur<arValsAll.shape[0]*arValsAll.shape[1]:
            arValsAllLWM[-(arValsAll.shape[0]*arValsAll.shape[1] - nValsCur):] = 1
        
        arValsAllLastWnd = ma.masked_array( arValsAllLastWnd, mask=arValsAllLWM )
        
        if len(arValsAllButLastWnd.shape)==1:
            arValsAllButLastWnd.reshape( (1,arValsAllButLastWnd.shape[0]) )
            
        arValsAllLastWnd=ma.reshape( arValsAllLastWnd, (1,arValsAllLastWnd.shape[0]) )

        if arValsAllButLastWnd.shape[0]==0:
            #arValsWnds = fxnProcessWnds( trkgrp, chrom, (lcorngWndBndsCur[0,0]+1,lcorngWndBndsCur[0,1]+1), arValsAll )
            arValsWnds = fxnProcessWnds( trkgrp, chrom, (lcorngWndBndsCur[0,0]+1,lcorngWndBndsCur[0,1]+1), arValsAllLastWnd )
        else:
            arValsWnds = np.r_[ fxnProcessWnds( trkgrp, chrom, (lcorngWndBndsCur[0,0]+1,lcorngWndBndsCur[-2,1]+1), arValsAllButLastWnd ),
                                fxnProcessWnds( trkgrp, chrom, (lcorngWndBndsCur[-1,0]+1,lcorngWndBndsCur[-1,1]+1), arValsAllLastWnd ) ]
               
        arVals = np.r_[ arVals, arValsWnds ]
        
    return lcorngWndBnds+1, lIdxRegionPerWnd, arVals

###

# returns array(summary val per wnd), array(#unmasked bases per wnd)

def extractPerVarWndQuickDefinedWindows( 
    trkgrp,      # x[chrom][coord:coord].   (i.e., x=wssd.depths;  can use x.wssd from fxnDataAccessor to get back "up" to wssd object; or x.dst if densetrackset)...
    chromQuery,
    lcorngWnds,
    fxnDataAccessor=lambda t,chrom,cor:t[chrom][cor[0]-1:cor[1]], # get the data: f( trackset, chrom, corng[11] ) --> array size [N,]  , where N=corng[1]-corng[0]+1,  
    fxnProcessWnds=lambda t,chrom,cor,V:V.mean(1),  # reproc the data in some way: f( trackset, chrom, corng[11], size[N,m] ) --> N values
    trkMask=None,
    arTypeStr='float32' ):
            
    corngWholeQuery=(lcorngWnds[0,0], lcorngWnds[-1,1])
    lofsWnds = lcorngWnds - corngWholeQuery[0]
    arSummaryVals,arNumUnmasked=np.zeros((lcorngWnds.shape[0],),arTypeStr),np.zeros((lcorngWnds.shape[0],),'int32')

    if trkMask!=None:
        mask = trkMask[ chromQuery ][ corngWholeQuery[0]-1:corngWholeQuery[1] ]
    else:
        mask = np.zeros( (corngWholeQuery[1]-corngWholeQuery[0]+1), 'bool' ) 
        
    arValsAll = fxnDataAccessor( trkgrp, chromQuery, (corngWholeQuery[0]+1, corngWholeQuery[1]+1) )
    
    arValsAllM = ma.masked_array( arValsAll, mask=mask )

    for iw in xrange(lofsWnds.shape[0]):
        arSummaryVals[iw] = fxnProcessWnds( 
                trkgrp, 
                chromQuery, 
                (lcorngWnds[iw,0], lcorngWnds[iw,1]), 
                arValsAllM[ lofsWnds[iw,0]:lofsWnds[iw,1]+1 ] )
        arNumUnmasked[iw] = (~arValsAllM.mask[ lofsWnds[iw,0]:lofsWnds[iw,1]+1 ]).sum()
    
    return arSummaryVals, arNumUnmasked
