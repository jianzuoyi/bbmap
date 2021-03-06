#This script is designed to preprocess and map data for variation calling of 2x150bp reads from Illumina HiSeq 2500.
#Some numbers and steps may need adjustment for different data types or file paths.
#For large genomes, tadpole and bbmerge may need the flag "prefilter=2" to avoid running out of memory 


#Remove duplicates
clumpify.sh in=reads.fq.gz out=clumped.fq.gz dedupe optical

#Remove low-quality regions
filterbytile.sh in=clumped.fq.gz out=filtered_by_tile.fq.gz

#Trim adapters
bbduk.sh in=filtered_by_tile.fq.gz out=trimmed.fq.gz ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=100 ref=bbmap/resources/adapters.fa ftm=5 ordered

#Remove synthetic artifacts and spike-ins.  Add "qtrim=r trimq=8" to also perform quality-trimming at this point, but not if quality recalibration will be done later.
bbduk.sh in=trimmed.fq.gz out=filtered.fq.gz k=27 ref=bbmap/resources/sequencing_artifacts.fa.gz,bbmap/resources/phix174_ill.ref.fa.gz ordered

#Map to reference
bbmap.sh in=filtered.fq.gz out=mapped.sam.gz bs=bs.sh pigz unpigz ref=reference.fa

#Call variants
callvariants.sh in=mapped.sam.gz out=vars.txt vcf=vars.vcf.gz ref=reference.fa ploidy=1 prefilter



#Optional error-correction and recalibration for better quality:

#Generate recalibration matrices
calctruequality.sh in=mapped.sam.gz vcf=vars.vcf.gz

#Recalibrate.  This can be done on the mapped reads instead of remapping, but if error-correction is desired it needs to be done on the unmapped reads.
bbduk.sh in=filtered.fq.gz out=recal.fq.gz recalibrate ordered

#Error-correct by overlap
bbmerge.sh in=recal.fq.gz out=ecco.fq.gz ecco strict mix ordered

#Quality-trim, if not already done earlier
bbduk.sh in=ecco.fq.gz out=qtrimmed.fq.gz qtrim=r trimq=8 ordered

#Re-map to the reference
bbmap.sh in=ecco.fq.gz out=mapped2.sam.gz pigz unpigz

#Re-call variants
callvariants.sh in=mapped2.sam.gz out=vars2.txt vcf=vars2.vcf.gz ref=reference.fa ploidy=1 prefilter
