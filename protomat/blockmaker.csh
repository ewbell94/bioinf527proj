#!/bin/csh 
#      blockmaker.csh   Run PROTOMAT on user-provided sequences
#   Type 		blockmaker.csh <name> 
#   and name.pros is the file of sequences in fasta format
#	EG  blockmaker.csh example/lipo

#   Or visit	http://www.blocks.fhcrc.org/blockmkr/make_blocks.html

unalias rm
unalias mv

#    Split up the input file name
#  If $1 = ~ptest/tmp/AA/AA, then $r = ~ptest/tmp/AA/AA and $t = AA
#  motomat dumps all its files into the current directory

set r = $1:r
set t = $r:t
rm -f $r.out $r.warn $t*.blk >& /dev/null
touch $r.out
touch $r.warn

#	Find out how many sequences there are
set nseq = (`grep -c ">" $r.pros`)

#	
if ( $nseq > 400 ) then
   echo -n "ERROR: Your input has $nseq sequences which" >> $r.warn
   echo " exceeds Block Maker's limit of 400 sequences" >> $r.warn
   exit(-1)
endif


#============================================================================
#	Make blocks from the proteins (motifj executes motomat)
echo "              **BLOCKS from MOTIF**" >> $r.out
echo " " >> $r.out
#	For dups:  motifj 4 -$r.pros 0 dups 17 (0 seqs => n/2 to start)
motifj 4 -$r.pros >& /dev/null
#	motifj writes $t.motifj.pros with MINIMUM seq len marked for gibbsj
rm $t.motifj.pros

#  Run motomat again so sequences aren't clumped/re-ordered
motomat $t.mot 1 1 -10 >& /dev/null
rm motomat.err $t.mot $t.plt $t*.old

if ( -e $t.blk || -e {$t}A.blk ) then
   cat $t*.blk > $r.blks
#  Add sequences to blocks if more than one block and add sequence weights
   if ( -e {$t}A.blk ) then
      addseqs $r.blks $r.pros $r.addblks 
      blweight $r.addblks $r.mblks P M 
      rm $r.addblks
   else
      blweight $r.blks $r.mblks P M 
   endif
   rm $t*.blk $r.blks

#	Display blocks in a "multiple alignment" format
   blalign $r.mblks >> $r.out

#	Make cobbler sequence = $r.mcob
   echo "TY	2"	     > $r.cf
   echo "BL	$r.mblks"   >> $r.cf
   echo "DB	$r.pros"    >> $r.cf
   echo "OU	$r.mcob"    >> $r.cf
   echo "SU	$BLIMPS_DIR/docs/default.iij" >> $r.cf
   cobbler $r.cf >& /dev/null
   rm $r.cf
#		Show the cobbler sequence now
   echo "" >> $r.out
   echo "" >> $r.out
   echo "      **COBBLER sequence from MOTIF**" >> $r.out
   cat $r.mcob | tr '\015' ' ' >> $r.out
else
   echo "ERROR: No blocks produced by MOTIF" >> $r.out
endif

echo "" >> $r.out
echo "" >> $r.out
echo "     **BLOCKS in searchable format**" >> $r.out
cat $r.mblks >> $r.out

#	Output
cat $r.warn $r.out

exit(0)
