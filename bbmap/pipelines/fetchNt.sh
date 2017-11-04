#Fetches and sketches nt.
#Be sure taxonomy is updated first!
#To use this script outside of NERSC,
#add "taxpath=/path/to/taxonomy_directory/"
#or, alternatively, replace "auto" with the appropriate paths

module load pigz
wget -nv ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz
mv nt.gz nt.fa.gz
time gi2taxid.sh -Xmx63g in=nt.fa.gz out=renamed.fa.gz pigz=32 zl=8 tree=auto accession=auto table=auto ow shrinknames
time sortbyname.sh -Xmx63g in=renamed.fa.gz out=sorted.fa.gz ow taxa tree=auto gi=ignore fastawrap=255 zl=8 pigz=32 minlen=60

time sketchblacklist.sh -Xmx63g in=sorted.fa.gz prepasses=1 tree=auto taxa taxlevel=species ow out=blacklist_nt_species_1000.sketch mincount=1000 k=31,24
time bbsketch.sh -Xmx63g in=sorted.fa.gz out=taxa#.sketch mode=taxa tree=auto accession=null gi=null files=31 ow unpigz minsize=300 prefilter autosize blacklist=blacklist_nt_species_1000.sketch k=31,24
