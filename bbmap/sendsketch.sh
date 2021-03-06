#!/bin/bash
#sendsketch in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified  August 28, 2017

Description:  Compares a Sketches to others, and prints their kmer identity.
The files can be sketches made by sketch.sh, or fasta files.
This is similar to comparesketch.sh, but connects to a remote server
which is holding the reference sketeches in memory.

Usage:  sendsketch.sh in=file


Standard parameters:
in=<file>           Sketch or fasta file to compare.
out=stdout          Can be set a file instead.
local=f             For local files, have the server load the sketches.
                    Allows use of whitelists; recommended for Silva.
                    Local can only be used when the client and server access 
                    the same filesystem - e.g., Genepool and Cori.
address=            Address of remote server.  Default address:
                    https://nt-sketch.jgi-psf.org/sketch
                    You can also specify these abbreviations:
                       nt:      nt server
                       refseq:  Refseq server
                       silva:   Silva server
                       img:     IMG server (not yet available)
                    Using an abbreviation automatically sets the address, 
                    the blacklist, and k.
blacklist=<file>    Ignore keys in this sketch file.  Additionaly, there are
                    built-in blacklists that can be specified:
                       nt:      Blacklist for nt
                       refseq:  Blacklist for Refseq
                       silva:   Blacklist for Silva
                       img:     Blacklist for IMG
amino=f             Use amino acid mode.

Sketch-making parameters:
mode=single         Possible modes, for fasta input:
                       single: Generate one sketch per file.
                       sequence: Generate one sketch per sequence.
size=10000          Size of sketches to generate, if autosize=f.
                    For raw PacBio data, 100k is suggested.
maxfraction=0.01    (mgf) Max fraction of genomic kmers to use.
autosize=t          Use flexible sizing instead of fixed-length.
autosizefactor=1    Multiply the default size of sketches by this factor.
k=31                Kmer length, 1-32.  To maximize sensitivity and 
                    specificity, dual kmer lengths may be used:  k=31,24
                    Dual kmers are fastest if the shorter is a multiple 
                    of 4.  Query and reference k must match.
                    Currently, JGI-hosted nt and RefSeq use k=31,24
                    while JGI-hosted Silva uses k=31, but this is set for you.
keyfraction=0.2     Only consider this upper fraction of keyspace.
samplerate=1        Set to a lower value to sample a fraction of input reads.
                    For raw reads (rather than an assembly), 1-3x coverage
                    gives best results, by reducing error kmers.  Somewhat
                    higher is better for high-error-rate data like PacBio.
minkeycount=1       Ignore kmers that occur fewer times than this.  Values
                    over 1 can be used with raw reads to avoid error kmers.
sketchheapfactor=4  If minkeycount>1, temporarily track this many kmers until
                    counts are known and low-count kmers are discarded.

Taxonomy-related parameters:
tree=<file>         Specify a TaxTree file.  On Genepool, use tree=auto.
                    Only necessary for use with printtaxa and level.
                    Assumes comparisons are done against reference sketches
                    with known taxonomy information.
level=2             Only report the best record per taxa at this level.
                    Either level names or numbers may be used.
                       -1: disabled
                        1: subspecies
                        2: species
                        3: genus
                       ...etc

Output columns:
printall=f          Enable all output columns.
printani=t          (ani) Print average nucleotide identity estimate.
completeness=t      Genome completeness estimate.
score=f             Score (used for sorting the output).
printmatches=t      Number of kmer matches to reference.
printlength=f       Number of kmers compared.
printtaxid=t        NCBI taxID.
printimg=f          IMG identifier (only for IMG data).
printgbases=f       Number of genomic bases.
printgkmers=f       Number of genomic kmers.
printgsize=t        Estimated number of unique genomic kmers.
printgseqs=t        Number of sequences (scaffolds/reads).
printtaxname=t      Name associated with this taxID.
printname0=t        (pn0) Original seqeuence name.
printfname=t        Query filename.
printtaxa=f         Full taxonomy of each record.
printcontam=t       Print contamination estimate, and factor contaminant kmers
                    into calculations.  Kmers are considered contaminant if
                    present in some ref sketch but not the current one.
printunique=t       Number of matches unique to this reference.
printnohit=t        Number of kmers that don't hit anything.
printucontam=f      Contam hits that hit exactly one reference sketch.

Other output parameters:
minhits=3           (hits) Only report records with at least this many hits.
minani=0            (ani) Only report records with at least this ANI (0-1).
minwkid=0.0001      (wkid) Only report records with at least this WKID (0-1).
records=20          Report at most this many best-matching records.
color=family        Color records at the family level.  color=f will disable.
                    Colors work in most terminals but may cause odd characters
                    to appear in text editors.  So, color defaults to f if 
                    writing to a file and t if writing to stdout.
format=2            Available formats:
                      1: Original format; deprecated.
                      2: Customizable with different columns.
                      3: 3-column format for alltoall: query, ref, ANI.

Metadata flags (optional, for the query sketch header):
taxid=-1            Set the NCBI taxid.
imgid=-1            Set the IMG id.
spid=-1             Set the sequencing project id (JGI-specific).
name=               Set the name (taxname).
name0=              Set name0 (normally the first sequence header).
fname=              Set fname (normally the file name).


Java Parameters:
-Xmx            This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

For more detailed information, please read /bbmap/docs/guides/BBSketchGuide.txt.
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

z="-Xmx4g"
z2="-Xms4g"
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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

sendsketch() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP sketch.SendSketch $@"
#	echo $CMD >&2
	eval $CMD
}

sendsketch "$@"
