#!/bin/csh
#    This is uextract.csh
#    Runs uextract on all .lis files 
#    Uextract runs motifj which runs motomat.
#    Removes extracted proteins when done.
#    Requires programs uextract, motifj, motomat, blweight in $path

#	NOTE: set $swiss (fasta format) before running and create
#		lsts/, blks/ and mots/
#
limit coredumpsize 1k
unalias rm
unalias cp
unalias mv

#WARNING: Be sure $swiss has full titles so uextract picks up FRAGMENT
set swiss = /fish/temp/swiss35.uni

#   Remove old statistics files
mv uextract.dat uextract.dat.bak
mv motomat.err motomat.err.bak

set list = (liss/PS*.lis)
foreach lis ($list)
   set psl = $lis:t
   set ps = $psl:r
#	uextract executes motifj which executes motomat
   uextract $lis $swiss > /dev/null
#	$ps.lst lists the sequences actually used to make the blocks,
#	differs from $ps.lis because it should not include fragments,
#	or sequences with very similar SWISS-PROT names
   mv $ps.lst lsts

#	Add sequence weights to the blocks produced by motomat
   if (-e $ps.blk || -e $ps{A}.blk) then
      set blks = ($ps*.blk)
      foreach blk ($blks)
         blweight $blk blks/$blk P M >& /dev/null
      end

#	$ps.mot is the binary file of motifs produces by motifj and
#	read by motomat
      mv $ps.mot mots
      compress mots/$ps.mot
#	This files were produced by motomat
      rm $ps.plt $ps*.blk
   endif

#	Remove the extracted sequences
   rm $ps.motifj.pros pros
   rm -r $ps
end
exit(0)
