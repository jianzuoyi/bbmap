#!/bin/bash
#mergesam in=<infile> out=<outfile>

function usage(){
echo "
Written by Brian Bushnell
Last modified March 8, 2017

Description:  Concatenates sam files, keeping only the header from the first.

Usage:  mergesam.sh <files> out=<file>

Java Parameters:
-Xmx    This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
        -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

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

z="-Xmx200m"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
}
calcXmx "$@"

function mergesam() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
		module load samtools/1.4
	fi
	local CMD="java $EA $z -cp $CP jgi.MergeSam $@"
	echo $CMD >&2
	eval $CMD
}

mergesam "$@"
