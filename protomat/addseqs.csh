#!/bin/csh
#		addseqs.csh
#	Try to add complete sequences to existing blocks
#	Must be in order & all must have calibrated alignment score
#		>= 800
# NOTE: Set $swiss (fasta format) and create addblks/
#
unalias rm
unalias mv

set swiss = /fish/temp/swiss35.uni

set liss = (liss/PS*.lis)
foreach lis ($liss)
   set temp = $lis:t; set ps = $temp:r
   set pros = $ps.pros

#	Extract all complete sequences, not just those for .lst file
   uextract $lis $swiss -o$pros -n >& /dev/null
#
   cat blks/$ps*.blk > $ps.blks
   addseqs $ps.blks $pros $ps.newblks >& /dev/null

   blweight $ps.newblks addblks/$ps.blks P M >& /dev/null

   rm $pros $ps.blks $ps.newblks 
end
exit(0)
