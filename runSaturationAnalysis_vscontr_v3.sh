#!/bin/bash 

set -ex

# functions
source $UPSCb/src/bash/functions.sh

# test 
#isEnvVarSet $UPSCb

### FOR SE DATA ONLY!!!  ######

# usage 
USAGETXT=\
"
Usage: runSaturationAnalysis.sh [options] <output directory> <bwa index> <input FastQ file> <control file>

Options:
  -c clean temporary data (FastQ, BAM and MACS2)
  -g genome size in bp (defaults to 400 Mbp (Aspen))
  -t number of threads to use
  
Note:
  The fastq files are expected to be gzipped
"

# defaults
CLEAN=1
CPUs=12
GSIZE=180000000
# options
while getopts cg:t: option
do
  case "$option" in
	    c) CLEAN=1;;
	    g) GSIZE=$OPTARG;;
	    t) CPUS=$OPTARG;;
		  \?) ## unknown flag
		  usage;;
  esac
done
shift `expr $OPTIND - 1`

# input file / output directory
if [ $# -lt 4 ]; then
  abort "This function expects four arguments at least."
fi

if [ ! -d $1 ]; then
  abort "The output directory does not exist."
fi

if [ ! -f $2 ]; then
  abort "The index file does not exist."
fi

if [ ! -f $3 ]; then
  abort "The input file does not exist."
fi

if [ ! -f $4 ]; then
  abort "The control file does not exist."
  else
    cont=$4
fi


# length of the input file divided by 10
fsize=$(expr $(gunzip -c $3 | wc -l) "/" 40)

# basename of the input file #removed _[1,2] for SE
fnam=$(basename ${3/.f*q.gz/})

# create the out dir
out=$1/$fnam

# remove the out
shift
inx=$1
shift

# loop
for i in {1..9}; do

  pct=$(expr $i "*" 10)-percent
  
  outp=$out/$pct
  if [ ! -d $outp ]; then
    mkdir -p $outp
  fi
  
     # sampleN
    if [ ! -f $outp/$pct ]; then 
      sampleN -n $(expr $fsize "*" $i) -o $outp/$pct $1
    fi
    
    # BWA
    if [ ! -f $outp/${pct}.sorted.bam ]; then 
      $UPSCb/pipeline/runBwamem.sh -t $CPUs -s -f $outp/${pct}.fq.gz -i $inx -o $outp 
    fi
  
    # MACS2
    # no need to subset the control, MACS2 scales it down
    $UPSCb/pipeline/runMacs2.sh -g $GSIZE $outp/${pct}.sorted.bam $cont $outp 

    # clean
    if [ $CLEAN -eq 1 ]; then
      rm -f $outp/${pct}.fq.gz $outp/${pct}.sorted.bam $outp/${pct}.sorted_treat_pileup.bdg
    fi
 
done
