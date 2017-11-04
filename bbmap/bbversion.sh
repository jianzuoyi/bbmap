#!/bin/bash
#bbversion

usage(){
echo "
Written by Brian Bushnell
Last modified May 4, 2017

Description:  Prints the BBTools version.

Usage:  bbversion.sh
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

CP="$DIR""current/"

bbversion() {
	#if [[ $NERSC_HOST == genepool ]]; then
	#	module load oracle-jdk/1.8_64bit
	#fi
	local CMD="java -Xmx80m -cp $CP driver.BBVersion $@"
	eval $CMD
}

bbversion "$@"
