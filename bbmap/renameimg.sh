#!/bin/bash
#rename in=<infile> out=<outfile>

function usage(){
echo "
Written by Brian Bushnell
Last modified August 23, 2017

Description:  Renames img records to be prefixed by their id.

Usage:  renameimg.sh in=auto out=renamed.fa.gz


Parameters:
in=         3-column tsv with imgID, taxID, and file path.
            These files will have their sequences renamed and concatenated.
img=        Optional, if a different (presumably bigger) file will be used for taxonomic assignment.
            For example, in could be a subset of img, potentially with incorrect taxIDs.


Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
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

z="-Xmx1g"
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

function rename() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP tax.RenameIMG $@"
	echo $CMD >&2
	eval $CMD
}

rename "$@"
