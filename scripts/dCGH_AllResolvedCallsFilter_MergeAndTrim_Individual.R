#############################################################################
# This Rscript is for converting the dCGH all_resolved calls into a         #
# filtered and annotated CNV file.                                          # 
# Only calls with unique (nonSD non Gap) bases are kept                     # 
#############################################################################

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Parameter Files
#segDupFile <- "/net/eichler/vol3/home/bcoe/hg38/segDups_Gaps_merged_hg38.bed"
#repeatFile <- "/net/eichler/vol3/home/bcoe/hg38/hg38_noWM_pad36.bed"
#CNVCalls <- read.table("all_resolved_calls.gz",header=TRUE,as.is=TRUE)
#DTS_FilePrefix <- "/net/eichler/vol26/7200/SAGE_genomes/nobackups/SVResults/dCGH/DTS_files/500_bp_slide/500_bp_"

if (length(args)!=6) {
  stop("Usage: Rscript script.R sampleID segDupFile repeatFile CNVCalls DTS_FilePrefix Outfile", call.=FALSE)
} 

sampleID <- as.character(args[1])
segDupFile <- as.character(args[2])
repeatFile <- as.character(args[3])
CNVCalls <- read.table(as.character(args[4]), header=T, as.is=T)
DTS_FilePrefix <- as.character(args[5])
Outfile <- as.character(args[6])
library(rhdf5)

# helper function for output
writeCall <- function(chr,start,end,sample,copies,SDSize,rptSize){
  #this fucntion had inheritance checking in the per family version.
  #here its just a cat.
  cat(paste(chr,start,end,sample,copies,end-start,SDSize,rptSize,"\n"))
}

#helper function for SD bases
fractionSD <- function(chr,start,end,segDups){
  #as a safety returns 1 for 0 length regions
  dupBases <- 0
  #find containing
  #SD s------e
  #     s--e
  if (dim(segDups[segDups$chr == chr & segDups$start <= start & segDups$end >= end,])[1]>0){
    dupBases = end-start
  }else{
    #find overlap left
    #SD  s------e
    #        s---e
    leftOverlaps <- segDups[segDups$chr == chr 
                            & segDups$start <= start 
                            & start <= segDups$end
                            & segDups$end <= end,]
    if (dim(leftOverlaps)[1]>0){
      dupBases <- dupBases + (max(leftOverlaps$end) - start)
    }
    #find overlap right
    #SD   s------e
    #   s---e
    rightOverlaps <- segDups[segDups$chr == chr 
                             & segDups$start >= start
                             & segDups$start <= end 
                             & segDups$end >= end,]
    if (dim(rightOverlaps)[1]>0){
      dupBases <- dupBases + (end - min(rightOverlaps$start))
    }
    #find internal
    #SD   s---e
    #   s-------e
    internalOverlaps <- segDups[segDups$chr == chr 
                                & segDups$start >= start 
                                & segDups$end <= end,]
    if (dim(internalOverlaps)[1]>0){
      dupBases <- dupBases + sum(internalOverlaps$end-internalOverlaps$start)
    }
  }
  if (end <= start){
    return(c(1,1))
  }else{
    return(c(dupBases / (end - start),dupBases))
  }
}

# helper function to return unique segments of a call
generateUniqueRegions <- function(chr,bps,bpe,SDEdges){
  subRegions <- data.frame(chr=character(), start=integer(), end=integer(),stringsAsFactors=FALSE)
  intersectedEdges <- SDEdges[SDEdges$chr == chr & SDEdges$pos >= bps & SDEdges$pos <= bpe,]
  addedsubRegions <- 0
  if (dim(intersectedEdges)[1]>0){
    for (i in 1:dim(intersectedEdges)[1]){
      # is first edge a start?
      subRegionsAdd <- data.frame(chr=character(), start=integer(), end=integer(),stringsAsFactors=FALSE)
      if (addedsubRegions == 0){
        if (intersectedEdges$Type[i]=="start"){
          subRegionsAdd <- data.frame(chr=as.character(chr), start=as.integer(bps), end=intersectedEdges$pos[i],stringsAsFactors=FALSE)
        }else {
          #two options end of segment or start of next
          if (i < dim(intersectedEdges)[1]){
            subRegionsAdd <- data.frame(chr=as.character(chr), start=intersectedEdges$pos[i], end=intersectedEdges$pos[i+1],stringsAsFactors=FALSE)
          }else{
            subRegionsAdd <- data.frame(chr=as.character(chr), start=intersectedEdges$pos[i], end=as.integer(bpe),stringsAsFactors=FALSE)
          }
        }
        # next edge onwards continues in same manner but only if we hit a end
        #if this is an end add to next point
        #if this is a start move on.  
      }else if (intersectedEdges$Type[i]=="end"){
        if (i < dim(intersectedEdges)[1]){
          subRegionsAdd <- data.frame(chr=as.character(chr), start=intersectedEdges$pos[i], end=intersectedEdges$pos[i+1],stringsAsFactors=FALSE)
        }else{
          subRegionsAdd <- data.frame(chr=as.character(chr), start=intersectedEdges$pos[i], end=as.integer(bpe),stringsAsFactors=FALSE)
        }
      }
      addedsubRegions <- 1
      subRegions<-rbind(subRegions,subRegionsAdd)
    }
  }else{
    subRegionsAdd <- data.frame(chr=character(), start=integer(), end=integer(),stringsAsFactors=FALSE)
    subRegionsAdd <- data.frame(chr=as.character(chr), start=as.integer(bps), end=as.integer(bpe),stringsAsFactors=FALSE)
    subRegions<-rbind(subRegions,subRegionsAdd)
  }
  return(subRegions)
}


medianCN <- function(testRegions,sample){
  dtscp_all <- numeric(0)
  for (i in 1:dim(testRegions)[1]){
    dtscp <- h5read(paste(DTS_FilePrefix,sample,sep="")
                    ,paste("copy/",testRegions$chr[i],sep=""))
    dtscp_starts <- h5read(paste(DTS_FilePrefix,sample,sep="")
                           ,paste("starts/",testRegions$chr[i],sep=""))
    dtscp_stops <- h5read(paste(DTS_FilePrefix,sample,sep="")
                          ,paste("ends/",testRegions$chr[i],sep=""))
    #dtscp_region <- dtscp[dtscp_starts >= testRegions$start[i] & dtscp_stops <= testRegions$end[i]]
    
    dtscp_r <- dtscp[dtscp_starts <= testRegions$end[i] & dtscp_stops >= testRegions$start[i]]
    dtscp_starts_r <- dtscp_starts[dtscp_starts <= testRegions$end[i] & dtscp_stops >= testRegions$start[i]]
    dtscp_stops_r <- dtscp_stops[dtscp_starts <= testRegions$end[i] & dtscp_stops >= testRegions$start[i]]
    dtscp_pos_r <- round((dtscp_stops_r-dtscp_starts_r)/2,0)+dtscp_starts_r
    dtscp_region <- dtscp_r[dtscp_pos_r >= testRegions$start[i] & dtscp_pos_r <= testRegions$end[i]] 
    dtscp_all <- c(dtscp_all,dtscp_region)
  }
  if (length(dtscp_all)>0){
    return(median(dtscp_all))
  }else{
    return(NaN)
  }
}

innerEdges <- function(chr,bps,bpe,sample,copies,segDups){
  # new version will trim edges to first threshold passing DTS window
  # if start or end are in SD leave as is
  leftOverlaps <- segDups[segDups$chr == chr 
                          & segDups$start <= bps 
                          & bps <= segDups$end
                          & segDups$end <= bpe,]

  rightOverlaps <- segDups[segDups$chr == chr 
                           & segDups$start >= bps
                           & segDups$start <= bpe 
                           & segDups$end >= bpe,]
  dtscp <- h5read(paste(DTS_FilePrefix,sample,sep="")
                  ,paste("copy/",chr,sep=""))
  dtscp_starts <- h5read(paste(DTS_FilePrefix,sample,sep="")
                         ,paste("starts/",chr,sep=""))
  dtscp_stops <- h5read(paste(DTS_FilePrefix,sample,sep="")
                        ,paste("ends/",chr,sep=""))
  #change start and end if they are not SD
  nbps <- bps
  nbpe <- bpe
  if (!(dim(rightOverlaps)[1]>0)){
    if (copies < 2){
      nbpe <- dtscp_stops[dtscp_starts > bps - 1000 & dtscp_stops < bpe + 1000 & dtscp <= 1.5]
    }else{
      nbpe <- dtscp_stops[dtscp_starts > bps - 1000 & dtscp_stops < bpe + 1000 & dtscp >= 2.5]
    }
  }
  if (!(dim(leftOverlaps)[1]>0)){
    if (copies < 2){
      nbps <- dtscp_starts[dtscp_starts > bps - 1000 & dtscp_stops < bpe + 1000 & dtscp <= 1.5]
    }else{
      nbps <- dtscp_starts[dtscp_starts > bps - 1000 & dtscp_stops < bpe + 1000 & dtscp >= 2.5]
    } 
  }
  if (length(nbps)>1 & length(nbpe)>1){
    pos<-round((nbpe-nbps)/2,0)+nbps
    return(c(min(pos),max(pos)))
  }else{
    return(c(min(nbps),max(nbpe)))
  }
}

#helper function that takes the merged call to output and runs all genotyping
masterWriter <- function(chr,start,end,indiv,CN){
  if(CN <1.5 | CN > 2.5){
    #rewrite start and end to DTS clipped ends
    newedges<-innerEdges(chr,start,end,indiv,CN,segDups)
    start <- newedges[1]
    end <- newedges[2]
    #test SD
    SD <- fractionSD(chr,start,end,segDups)
    rpt <- fractionSD(chr,start,end,repeats)
    #generate some temp variables
    regionCN <- numeric(0)
    if (SD[1] < 0.8){
      #generate non SD base windows
      testRegions <- generateUniqueRegions(chr,start,end,SDEdges)
      if (dim(testRegions)[1]>0){        
        regionCN <- medianCN(testRegions,indiv)
        if (!is.nan(regionCN)&!is.na(regionCN)){
          if ((regionCN < 1.5 | regionCN > 2.5) & regionCN < 4.5){
            writeCall(chr,start,end,indiv,round(regionCN,3),SD[2],rpt[2])
          }
        }
      }
    }
  }
}


# prepare segdups and gaps for intersectionTesting
segDups <- read.table(segDupFile)
#segDups <- read.table("SegDups_hg19_Merged.txt")
colnames(segDups) <- c("chr","start","end")
SDStarts <- data.frame(segDups$chr,segDups$start, rep("start",dim(segDups)[1]))
colnames(SDStarts) <- c("chr","pos","Type")
SDEnds <- data.frame(segDups$chr,segDups$end, rep("end",dim(segDups)[1]))
colnames(SDEnds) <- c("chr","pos","Type")
SDEdges<-rbind(SDStarts,SDEnds)
SDEdges <- SDEdges[order(SDEdges$chr,SDEdges$pos),]
rm(SDStarts,SDEnds)

#prepare repeat masking file for intersecting
repeats <- read.table(repeatFile)
colnames(repeats) <- c("chr","start","end")

#filter to only sample
CNVCalls<-CNVCalls[CNVCalls$indiv %in% sampleID,]

#First Pass add CN back to calls
CNVCalls<-CNVCalls[order(CNVCalls$indiv,CNVCalls$chr,CNVCalls$start),]
CN <- rep(2,dim(CNVCalls)[1])
for (i in 1:dim(CNVCalls)[1]){
#for (i in 1:200){
  testRegions <- generateUniqueRegions(CNVCalls$chr[i],CNVCalls$start[i],CNVCalls$end[i],SDEdges)
  if (dim(testRegions)[1]>0){
    if(sum(testRegions$end-testRegions$start) > 500){
      CN[i] <- medianCN(testRegions,CNVCalls$indiv[i])
    }
  }
}

CNVCalls$CN <- CN
# delete CN = 2 aka not tested lines also pre-remove over 5 copies (later we remove 4.5 copies...)
# also drop median CN between 1.4 and 2.4 as these will not pass later.
CNVCalls <- CNVCalls[CNVCalls$CN != 2 & (CNVCalls$CN < 1.4 | CNVCalls$CN > 2.4) & CNVCalls$CN <= 5 & !is.nan(CNVCalls$CN),]

# Test Calls and Write output
sink(Outfile)
cat("chr\tstart\tend\tsample\tcopies\tsize\tSDSize\tRepeatSize\n")

#start at call #1 check next merger, repeat until not merge candidate

i <- 1

doneFlag <- 0

#initiate currentCall
currentCall <- CNVCalls[i,]
z<- i + 1
j <- dim(CNVCalls)[1]

while (doneFlag == 0){
  # am I the last call? - output and stop
  if (i >= j) {
    #write out
    masterWriter(currentCall$chr, currentCall$start, currentCall$end
                 , currentCall$indiv, currentCall$CN) 
    sink()
    doneFlag <- 1
    #is the next call a merge option (both <2 or both >2)?
  }else if(currentCall$chr == CNVCalls$chr[z] &
             ((currentCall$CN < 2 & CNVCalls$CN[z] < 2)
              | (currentCall$CN > 2 & CNVCalls$CN[z] > 2))
           & (CNVCalls$start[z] - currentCall$end) < 1000000){
    #define percent SD of Gap
    gapSD <- fractionSD(currentCall$chr,currentCall$end,CNVCalls$start[z],segDups)
    #define percent repeat of Gap only for small gaps
    if ((CNVCalls$start[z] - currentCall$end) < 5000){
      gaprpt <- fractionSD(currentCall$chr,currentCall$end,CNVCalls$start[z],repeats)
    }else{
      gaprpt <-rep(0,2)
    }
    # define median CN of Gap
    dtscp <- h5read(paste(DTS_FilePrefix,currentCall$indiv,sep="")
                    ,paste("copy/",currentCall$chr,sep=""))
    dtscp_starts <- h5read(paste(DTS_FilePrefix,currentCall$indiv,sep="")
                           ,paste("starts/",currentCall$chr,sep=""))
    dtscp_stops <- h5read(paste(DTS_FilePrefix,currentCall$indiv,sep="")
                          ,paste("ends/",currentCall$chr,sep=""))
    dtscp_region <- dtscp[dtscp_starts>currentCall$end & dtscp_stops<CNVCalls$start[z]]
    #check if gap is large enough to have a CN
    if (length(dtscp_region)>0){
      gapCN <- median(dtscp_region)
    }else{
      if(currentCall$CN<2){
        gapCN = 1.4
      }else{
        gapCN = 2.6
      }
    }
    #test merging citeria
    # test call segdup (defined or CN > 3.5 and <50kbp for del?) in between OR CN = expected median 1.5
    if (((currentCall$CN < 2 & gapCN <= 1.5)
         |(currentCall$CN < 2 & gapCN >= 3.5 
           & (CNVCalls$start[z]-currentCall$end < 50000))
         |(currentCall$CN > 2 & gapCN >= 2.5)
         | (gapSD[1] + gaprpt[1]) > 0.75
         | (CNVCalls$start[z] - currentCall$end) < 1000
         | min((CNVCalls$start[z] - currentCall$end)/(CNVCalls$end[z]-CNVCalls$start[z]),(CNVCalls$start[z] - currentCall$end)/(CNVCalls$end[z]-CNVCalls$start[z])) < 0.01)
        & (CNVCalls$start[z] - currentCall$end) < 1000000){
      currentCall$end <- CNVCalls$end[z] 
      #did we just add the last call?
      if (z == j){
        masterWriter(currentCall$chr, currentCall$start, currentCall$end
                     , currentCall$indiv, currentCall$CN) 
        sink()
        doneFlag <- 1
      }else{
        z <- z + 1
      }
      #write current and move to next
    }else{
      masterWriter(currentCall$chr, currentCall$start, currentCall$end
                   , currentCall$indiv, currentCall$CN) 

      i<-z
      z<-i+1
      currentCall <- CNVCalls[i,]
    }
    # if not a merge candidate output and move to next
  }else{ 
    masterWriter(currentCall$chr, currentCall$start, currentCall$end
                 , currentCall$indiv, currentCall$CN) 
    i<-z
    z<-i+1
    currentCall <- CNVCalls[i,]
  }
}




