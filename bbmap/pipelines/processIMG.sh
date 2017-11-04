#!/bin/bash

#Combines IMG into a single flat file, renamed by taxonomy, and sketches it.
#This script only works on Genepool or other systems connected to projectb.

#time renameimg.sh in=auto out=renamed.fa.gz zl=6
time renameimg.sh in=auto imghq out=renamed.fa.gz fastawrap=255 zl=6
time sketchblacklist.sh -Xmx31g in=renamed.fa.gz prepasses=1 tree=auto taxa taxlevel=species ow out=blacklist_img_species_300.sketch mincount=300 k=31,24
time sketch.sh -Xmx31g in=renamed.fa.gz out=img#.sketch files=31 mode=img tree=auto img=auto gi=null ow blacklist=blacklist_img_species_300.sketch k=31,24
