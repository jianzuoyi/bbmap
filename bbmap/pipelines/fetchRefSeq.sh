#Fetches and sketches RefSeq.
#Be sure taxonomy is updated first!
#To use this script outside of NERSC,
#add "taxpath=/path/to/taxonomy_directory/" to each command
#or, alternatively, replace "auto" with the appropriate paths

module load pigz
time wget -nv ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*genomic.fna.gz
time cat *genomic.fna.gz > all.fa.gz

time gi2taxid.sh -Xmx63g in=all.fa.gz out=renamed.fa.gz tree=auto table=auto accession=auto zl=6
time sortbyname.sh -Xmx63g in=renamed.fa.gz out=sorted.fa.gz zl=6 pigz=32 taxa tree=auto gi=ignore fastawrap=255 minlen=60
time sketchblacklist.sh -Xmx63g in=sorted.fa.gz prepasses=1 tree=auto taxa taxlevel=species ow out=blacklist_refseq_species_300.sketch mincount=300 k=31,24
time bbsketch.sh -Xmx63g in=sorted.fa.gz out=taxa#.sketch mode=taxa tree=auto accession=null gi=null files=31 ow unpigz minsize=400 prefilter autosize blacklist=blacklist_refseq_species_300.sketch k=31,24
