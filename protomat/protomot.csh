#!/bin/csh
#		protomot.run 
#      Extract & modify .lis files produced by PROTOMOT
#      Edited for Blocks10/Prosite14
# NOTE:  Set $prosite before running & create liss/ subdirectory

set prosite=/fish/prosite/prosite.dat

unalias rm
#  0. Run protomot
protomot $prosite all none none
#
#  1. Remove known duplicates (same PDOC entry, don't end in _2 or _3)
#     NOTE: These are skipped by protomot if they end in _2 or _3
rm PS00080.lis
rm PS00108.lis
rm PS00109.lis
rm PS00119.lis
rm PS00135.lis
rm PS00137.lis
rm PS00138.lis
#rm PS00182.lis  same PDOC as PS00180, seqs are subset of PS00180, see below
rm PS00193.lis
rm PS00308.lis
rm PS00320.lis
#rm PS00336.lis PS00337.lis   different seqs from PS00146, keep
rm PS00382.lis
#rm PS00472.lis  different seqs from PS00471, keep
rm PS00560.lis
rm PS00589.lis
rm PS00637.lis
rm PS00638.lis
rm PS00639.lis
rm PS00640.lis
rm PS00673.lis
rm PS00687.lis
#rm PS00690.lis   different seqs from PS00039, keep it
rm PS00805.lis
rm PS00950.lis
rm PS00969.lis
rm PS01080.lis PS01258.lis PS01259.lis PS01260.lis
#rm PS01045.lis  _2 so protomot won't extract
rm PS01056.lis
#rm PS01076.lis  _2 so protomot won't extract
rm PS01122.lis
#rm PS01124.lis  _2 so protomot won't extract
rm PS01132.lis
rm PS01163.lis
rm PS01174.lis
rm PS50008.lis
rm PS50011.lis
rm PS50045.lis
rm PS50054.lis PS50055.lis
rm PS50065.lis
rm PS50067.lis
rm PS50075.lis
#
#   1a. Remove subsets contained in other groups
#	PS00033 is a subset of PS00027
rm PS00033.lis
#	PS00181 & PS00182 are subsets of PS00180
rm PS00181.lis
rm PS00182.lis
#	PS00228 is a subset of PS00227
rm PS00228.lis
#	PS00238 is a subset of PS00237
rm PS00238.lis
#
#  3. Remove groups for sites
#rm PS0000[1-9].lis rm PS0001[0-7].lis 
rm PS00010.lis PS00011.lis PS00012.lis PS00014.lis
#rm PS00294.lis
rm PS00342.lis PS00409.lis
#
#  4. Add groups skipped because of PROTOMOT's "_2" or "_3" criterion
#     NOTE:  PROTOMOT will extract "_4", "_5", etc.
protomot $prosite PS00023 PS00023 none
protomot $prosite PS00353 PS00353 none
#protomot $prosite PS00372 PS00372 none
protomot $prosite PS00424 PS00424 none
protomot $prosite PS00486 PS00486 none
protomot $prosite PS00573 PS00573 none
#	PS00574 is _2 after _1, but get it anyway
protomot $prosite PS00574 PS00574 none
protomot $prosite PS00599 PS00599 none
protomot $prosite PS00600 PS00600 none
protomot $prosite PS01039 PS01039 none
protomot $prosite PS01213 PS01213 none
protomot $prosite PS01240 PS01240 none
#
#The following pairs have the same PDOC, distinct names, but different seqs
#PS00039 & PS00690
#PS00059 & PS01162
#PS00336 & PS00337
#PS00471 & PS00472
#PS00685 & PS00686
#
mv PS*.lis liss
#
#	Really big families (> 400 seqs+dups) have to be processed separately
mv liss/PS00018.lis liss/PS00018.big
mv liss/PS00022.lis liss/PS00022.big
mv liss/PS00027.lis liss/PS00027.big
mv liss/PS00028.lis liss/PS00028.big
mv liss/PS00107.lis liss/PS00107.big
mv liss/PS00211.lis liss/PS00211.big
mv liss/PS00237.lis liss/PS00237.big
mv liss/PS01033.lis liss/PS01033.big
mv liss/PS01187.lis liss/PS01187.big
#
exit(0)

