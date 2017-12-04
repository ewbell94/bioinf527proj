#!/bin/csh
#    This is PS00094.csh
#    Runs  PROTOMAT on lsts/PS00094.lst file
#    Assmes the proteins are in PS00094/
#
unalias rm
unalias cp
unalias mv

#   Remove old statistics files
mv motomat.err motomat.err.bak >& /dev/null

set list = (lsts/PS00094.lst)
foreach lis ($list)
   set psl = $lis:t
   set ps = $psl:r
   motifj 4 $lis >& PS00094.out
   mv $ps*.blk blks
   rm $ps.plt $ps.mot $ps.motifj.pros
end
exit(0)
