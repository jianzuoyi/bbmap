#!/bin/bash

function usage(){
echo "
Written by Brian Bushnell
Last modified July 27, 2015

Description:  Maps reads, then filters an assembly by contig coverage.
Intended to reduce misassembly rate of SPAdes by removing suspicious contigs.

Usage:  postfilter.sh in=<reads> ref=<contigs> out=<filtered contigs>


Standard Parameters:
in=<file>           File containing input reads.
in2=<file>          Optional file containing read mates.
ref=<file>          File containing input assembly.
cov=covstats.txt    File to write coverage stats generated by pileup.
out=filtered.fa     Destination of clean output assembly.
outdirty=<file>     (outd) Destination of removed contigs; optional.
ow=f                (overwrite) Overwrites files that already exist.
app=f               (append) Append to files that already exist.
zl=4                (ziplevel) Set compression level, 1 (low) to 9 (max).
int=f               (interleaved) Determines whether input reads are considered interleaved.

Filtering Parameters:
minc=2              (mincov) Discard contigs with lower average coverage.
minp=95             (minpercent) Discard contigs with a lower percent covered bases.
minr=6              (minreads) Discard contigs with fewer mapped reads.
minl=400            (minlength) Discard shorter contigs.
trim=0              (trimends) Trim the first and last X bases of each sequence.

Mapping Parameters (unlisted params will use BBMap defaults)
minhits=2
maxindel=0
tipsearch=0
bw=20
rescue=f


Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

Other parameters will be passed directly to BBMap.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx800m"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
	if [[ $set == 1 ]]; then
	return
	fi
	freeRam 800m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

function postfilter() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
		module load samtools/1.4
	fi
	local CMD="java $EA $z -cp $CP assemble.Postfilter $@"
	echo $CMD >&2
	eval $CMD
}

postfilter "$@"
