/* (C) Copyright 1991 by Fred Hutchinson Cancer Research Center         */
/*    motomat2.c is second part of motomat.c; see comments there        */
/*======================================================================*/
/*   prune_blocks sets TopScore so there are <= MAXBLK surviving blocks.*/
/*   NOTE:  Not currently used. */
/*======================================================================*/
int prune_blocks(blocks)
struct merged_motif *blocks;
{
   struct temp *temp;
   int b, ntemp;

   temp = (struct temp *) malloc(Total_Motifs*sizeof(struct temp));
   if (temp == NULL)
   {
      fprintf(stderr,"\n%.8s", Title);
      fprintf(stderr, " prune_blocks: Unable to allocate temp structure!\n");
      restart();
   }
/*---------  Sort blocks by trimmed score  */
   ntemp = 0;
   for (b=0; b<Total_Motifs; b++)
      if (blocks[b].nmotif > 0)
      {
	 temp[ntemp].value = blocks[b].t_score;
	 temp[ntemp].index = b;
	 temp[ntemp++].flag = 0;
      }
   qsort(temp, ntemp, sizeof(struct temp), tempcmp);
/*    A descending sort would make it easier... */
   if (ntemp > MAXBLK)
   {
      b = ntemp - MAXBLK;  ntemp = MAXBLK;
   }
   else            b = 0;            /*  take min score if <= 25 blocks */
   TopScore = temp[b].value;
   printf("\nRevised TopScore=%d", TopScore);
   free(temp);
   return(ntemp);
}  /*  end of prune_blocks */
/*=====================================================================*/
/*   order_blocks orders all blocks within each sequence in ascending
      order by the left-most amino acid in each block.  Blocks may
      overlap.  Orders only unmerged blocks that have not been dropped
      because of low score.                                            */
/*=====================================================================*/
void order_blocks(blocks)
struct merged_motif *blocks;
{
   struct temp *temp;
   int s, b, ntemp, i;

   temp = (struct temp *) malloc(Total_Motifs*sizeof(struct temp));
   if (temp == 0)
   {
      fprintf(stderr,"\n%.8s", Title);
      fprintf(stderr," order_blocks: Unable to allocate temp structure!\n");
      restart();
   }
   for (s=0; s<NumSeqs; s++)
   {
      ntemp = 0;
      for (b=0; b<Total_Motifs; b++)
	 if (!blocks[b].dropped)
	 {
	    temp[ntemp].value = blocks[b].leftpos[s] - blocks[b].t_loffset;
	    temp[ntemp].index = b;
	    temp[ntemp++].flag = 0;	/* field not used */
	 }
      qsort(temp, ntemp, sizeof(struct temp), tempcmp);
      if (Debug) printf("\nseq=%d: ", s);
      for (i=0; i<ntemp; i++)
      {
	 b = temp[i].index;
	 blocks[b].position[s] = i;
	 if (Debug) printf("%d,%d ", temp[i].index, temp[i].value);
      }
   }
   for (b=0; b<Total_Motifs; b++)
      if (!blocks[b].dropped)
      {
	 blocks[b].minpos = Total_Motifs; blocks[b].maxpos = 0-Total_Motifs;
	 for(s=0; s<NumSeqs; s++)
	 {
	    if (blocks[b].position[s] > blocks[b].maxpos)
	       blocks[b].maxpos = blocks[b].position[s];
	    if (blocks[b].position[s] < blocks[b].minpos)
	       blocks[b].minpos = blocks[b].position[s];
	 }
      }

   free(temp);
}  /*  end of order_blocks */

/*=====================================================================*/
/*   tempcmp sorts a temp structure in ascending order by value, index
       and flag when it is used by qsort.                              */
/*=====================================================================*/
int tempcmp(t1, t2)
struct temp *t1, *t2;
{
   if (t1->value != t2->value)  return(t1->value - t2->value);
   if (t1->index != t2->index)  return(t1->index - t2->index);
   return(t1->flag - t2->flag);

}  /*  end of tempcmp */

/*=====================================================================*/
/*   build_dag creates a directed acyclic graph where the nodes are the
       motif blocks and the arcs represent all orderings of the blocks
       among all sequences.  The graph is represented by an adjacency
       matrix.  Only those cells of the matrix for which the arc from
       the row block to the column block is positive are allocated.  An
       arc is positive if the column block comes after the row block in
       at least RSignif sequences.
       NOTE: In order to insure the graph contains no cycles, RSignif
       must be no less than NumSeqs/2, since different sequences may be
       in different arcs. */
/*  3/18/91 Changed so that b1 must start before b2 in at least RSignif
    seqs AND b1 cannot overlap b2 in those seqs in order for there to
    be an arc from b1 to b2 */
/*=====================================================================*/
void build_dag(blocks)
struct merged_motif *blocks;
{
   int row, col, s, diff, nposrc, nposcr, maxdiffrc, maxdiffcr;
   int distrc, distcr;

   for (row=0; row<Total_Motifs; row++)
      if (!blocks[row].dropped )
	      blocks[row].in_degree = blocks[row].out_degree = 0;

   for (row=0; row<Total_Motifs; row++)
      if (!blocks[row].dropped)
      {
	 for (col=row+1; col<Total_Motifs; col++)
	    if (!blocks[col].dropped)
	    {
	       nposrc = nposcr = 0;  maxdiffrc = maxdiffcr = -999;
	       for (s=0; s<NumSeqs; s++)
	       {  /* diff is block position, dist is overlap */
		  diff = blocks[col].position[s]-blocks[row].position[s];
		  distrc = left(s, col, blocks)-right(s, row, blocks)-1;
		  distcr = left(s, row, blocks)-right(s, col, blocks)-1;
		     /* positive is from row block to col block */
		  if (diff >= 0 && distrc >= 0)
		  {
		     nposrc += 1;
		     if (diff > maxdiffrc) maxdiffrc = diff;
		  }
		  else if (diff < 0 && distcr >= 0)
		  {
		     diff = 0-diff;
		     nposcr += 1;
		     if (diff > maxdiffcr) maxdiffcr = diff;
		  }
	       }  /*  end of sequence s */
	       if (nposrc >= RSignif)
	       {
		  Dag[row][col] = makematrix();
		  Dag[row][col]->npos = nposrc;
		  Dag[row][col]->maxdiff = maxdiffrc;
		  blocks[row].out_degree = blocks[row].out_degree + 1;
		  blocks[col].in_degree = blocks[col].in_degree + 1;
/*		  for (s=0; s<NumSeqs; s++)
		    Dag[row][col]->dist[s] =
		      left(s, col, blocks) - right(s, row, blocks) - 1;
*/
	       }
	       if (nposcr >= RSignif)
	       {
		  Dag[col][row] = makematrix();
		  Dag[col][row]->npos = nposcr;
		  Dag[col][row]->maxdiff = maxdiffcr;
		  blocks[col].out_degree = blocks[col].out_degree + 1;
		  blocks[row].in_degree = blocks[row].in_degree + 1;
/*		  for (s=0; s<NumSeqs; s++)
		    Dag[col][row]->dist[s] =
		      left(s, row, blocks) - right(s, col, blocks) - 1;
*/
	       }
	 }  /*  end of col */
      }  /*  end of row */
      if (Debug)
	 for (row=0; row<Total_Motifs; row++)
	 {
	    printf("\nrow %.3d: ", row);
	    for (col=0; col<Total_Motifs; col++)
	       if (Dag[row][col] != NULL)
		  printf(" %.3d", Dag[row][col]->npos);
	       else
		  printf("   x");
	 }

}  /*  end of build_dag */
/*======================================================================*/
/*   Make and return a new matrix structure. */
/*======================================================================*/
struct matrix *makematrix()
{
   struct matrix *new;

   new = (struct matrix *) malloc(sizeof(struct matrix));
   if (new == NULL)
   {
      fprintf(stderr,"\n%.8s", Title);
      fprintf(stderr," makematrix: Unable to allocate matrix structure!\n");
      restart();
   }
   new->npos = 0;
   new->maxdiff = 0;
   new->mark = -1;
   return(new);
}  /*  end of makematrix */
/*=======================================================================*/
/*   makepath allocates and intializes a new path.                       */
/*=======================================================================*/
struct path *makepath()
{
  struct path *new;
  int s, row, col;

  new  = (struct path *) malloc(sizeof(struct path));
  if (new == NULL)
  {
     fprintf(stderr,"\n%.8s", Title);
     fprintf(stderr," makepath: Unable to allocate path structure!\n");
/*---  Have to free enough memory to reload motomat ----------------------*/
     for (row=0; row<Total_Motifs; row++)
	 for (col=0; col<Total_Motifs; col++)
	    if (Dag[row][col] != NULL)  free(Dag[row][col]);
     restart();
  }
  new->nblocks = new->nbest = new->naas = new->totmotif = new->totident = 0;
  new->totscore = (unsigned long) 0;
  new->nseqs = 0;
  for (s=0; s<MAXSEQS; s++)
     new->seqs[s] = YES;
  new->first_block = new->first_best = NULL;
  new->next_path = NULL;
  return(new);
}  /* end of makepath */
/*=======================================================================*/
/*   ins_path inserts a new path at the end of a list of paths.          */
/*=======================================================================*/
void ins_path(list, new)
struct path *list, *new;
{
   struct path *cur;

   cur = list;
   while (cur->next_path != NULL)
     cur = cur->next_path;

   cur->next_path = new;
   list->nblocks = list->nblocks + 1;  /* This is the # of paths */
   if (list->nblocks > 1500)
   {
      fprintf(stderr,"\n%.8s", Title);
      fprintf(stderr," ***MOTOMAT STOPPING BECAUSE >1500 PATHS***");
      restart();
   }

}  /* end of ins_path */

/*=======================================================================*/
/*   copypath copies an existing path into a new path.                   */
/*=======================================================================*/
struct path *copypath(old)
struct path *old;
{
   struct path *new;
   struct block_list *bold, *bnew, *blas;
   int s;

   new = makepath();
   new->nblocks = old->nblocks; new->nbest = old->nbest;
   new->naas = old->naas; new->totscore = old->totscore;
   new->nseqs = old->nseqs;
   new->totmotif = old->totmotif;
   new->totident = old->totident;

   for (s=0; s<MAXSEQS; s++)
      new->seqs[s] = old->seqs[s];

   bold = old->first_block;
   if (bold != NULL)
   {                            /*  copy first block */
      bnew = makebllist();
      bnew->b = bold->b;
      new->first_block = bnew;
      bold = bold->next_block;
      while (bold != NULL)  /*  copy other blocks */
      {
	 blas = bnew;
	 bnew = makebllist();
	 bnew->b = bold->b;
	 blas->next_block = bnew;
	 bold = bold->next_block;
      }
   }
/*------Have to copy best path pointers: old->first_best, bold->next_best,
	bold->prev_best-----*/

   return(new);
}  /* end of copypath */
/*=======================================================================*/
/*   free_path frees all memory allocated to a path  .                   */
/*=======================================================================*/
void free_path(old)
struct path *old;
{
   struct block_list *bold, *bfree;

   if (old != NULL)
   {
      bold = old->first_block;
      while (bold != NULL)            /*  free blocks */
      {
	 bfree = bold;
	 bold = bold->next_block;
	 free(bfree);
      }
      free(old);	 	   /*  free path */
   }

}  /* end of free_path */
/*========================================================================*/
/*    find_paths finds all positive paths through the DAG.
       Does not consider whether or not blocks in the
       path overlap.  Stops building a path when fewer than RSignif
       sequences are represented.              */
/*=======================================================================*/
void find_paths(paths, blocks)
struct path *paths;
struct merged_motif *blocks;
{
   int row, col, more, ntemp, i;
   struct path *newpath;
   struct temp *temp;

/*---- First add all single blocks as possible paths of one block each --*/
   for (i=0; i<Total_Motifs; i++)
   {
      if (!blocks[i].dropped)
      {
	 newpath = makepath();
	 ins_bllist(i, newpath);
         if (Debug) printf("\nfind_paths: single %d", i);
	 best_path(newpath, paths, blocks);   /* was ins_path(paths,newpath) */
      }
   }

   temp = (struct temp *) malloc(Total_Motifs*sizeof(struct temp));
   if (temp == NULL)
   {
      fprintf(stderr,"\n%.8s", Title);
      fprintf(stderr," find_paths: Unable to allocate temp structure!\n");
      restart();
   }
   ntemp = 0;
/*--Process rows in order of lowest in-degree (topological sort of graph
    nodes) If two rows have same in-degree, the one with the larger
    out-degree sorts first ---*/
   for (row=0; row <Total_Motifs; row++)
      if (!blocks[row].dropped)
      {
	  temp[ntemp].index = row;
	  temp[ntemp].value = blocks[row].in_degree;  /* ascending */
	  temp[ntemp++].flag = Total_Motifs - blocks[row].out_degree;
						     /* descending */
      }

   qsort(temp, ntemp, sizeof(struct temp), tempcmp);
   for (i=0; i<ntemp; i++)
   {
      row = temp[i].index;
      more = NO;
      for (col=0; col<Total_Motifs; col++)
      {
	 if (Dag[row][col] != NULL &&
	     Dag[row][col]->mark < 0 &&
	     Dag[row][col]->npos >= RSignif)     /*  was NumSeqs */
		more = YES;
      }
/* ---------  See if there are more starting paths in this row ----------*/
      if (more)
      {
	 newpath = makepath();
	 ins_bllist(row, newpath);
         if (Debug) printf("\n\nfind_paths: more-----------");
         /*  don't need to call best_path() here because would
             have already been checked with singles
	 best_path(newpath, paths, blocks);  */
	 follow_arcs(paths, newpath, row, blocks);
      }
      else if (ntemp == 1)          /*  Only one block */
      {
	 newpath = makepath();
	 ins_bllist(row, newpath);
	 best_path(newpath, paths, blocks);
      }
   }

   free(temp);
}  /*  end of find_paths */

/*======================================================================*/
/*   follow_arcs does a recursive depth first search of the DAG,
     finding all postive paths.                                         */
/*======================================================================*/
void follow_arcs(paths, curpath, row, blocks)
struct path *paths, *curpath;
int row;
struct merged_motif *blocks;
{
   int col1, ntemp, ntemp1, nbranch, i, savnpath, s, nseq;
   struct path *newpath, *savpath;
   struct block_list *tblock;
   struct temp *temp;

   savpath = newpath = NULL;
   savnpath = paths->nblocks;  /* This is really the # of paths */
/*------  See if there is an unmarked fully positive column, and take
	  the one with the closest next position if there is ------------*/
/*-------  Make a list of unmarked columns for this row, sorted in
	     increasing order of maxdiff -----------*/
   temp = (struct temp *) malloc(Total_Motifs*sizeof(struct temp));
   if (temp == NULL)
   {
      fprintf(stderr,"\n%.8s", Title);
      fprintf(stderr," follow_arcs: Unable to allocate temp structure!\n");
      restart();
   }
/*--- Continue topological sort of graph nodes (lowest in-degree next)---*/
   ntemp = 0;
   for (col1=0; col1<Total_Motifs; col1++)
      if (Dag[row][col1] != NULL &&
	  Dag[row][col1]->mark < savnpath &&
	  Dag[row][col1]->npos >= RSignif)
	  {
	     temp[ntemp].index = col1;
	     temp[ntemp].value = blocks[col1].in_degree;
	     temp[ntemp++].flag = Total_Motifs - blocks[col1].out_degree;
	  }
    qsort(temp, ntemp, sizeof(struct temp), tempcmp);
    if (Debug)
    {
	printf("\n\nfollow:  row=%d, ntemp=%d", row, ntemp);
	printf("\n curpath blocks=");
	tblock = curpath->first_block;
	while (tblock != NULL)
	{
	    printf("%d ", tblock->b);
	    tblock = tblock->next_block;
	}
	printf("\n curpath_seqs:");
	for (s=0; s<NumSeqs; s++)
	    printf("%d=%d ", s, curpath->seqs[s]);
    }

/*--- There are ntemp possible ways to continue this path, but not
   all of them may have enough sequences remaining.  ntemp1 is the
   number of ways to continue that still have enough sequences -----*/
    ntemp1 = ntemp;
    for (i=0; i<ntemp; i++)
    {
       nseq = NumSeqs;
       for (s=0; s<NumSeqs; s++)
	  if (curpath->seqs[s] == NO ||
	     blocks[temp[i].index].position[s] < blocks[row].position[s])
		   nseq--;
       if (Debug)
	  printf("\n row=%d temp[i].index=%d nseq=%d",
		     row, temp[i].index, nseq);
       if (nseq >= RSignif)
	  temp[i].flag = YES;   /* NOTE: re-using flag (was out_degree)*/
       else
       {
	  temp[i].flag = NO;
	  if (Dag[row][temp[i].index] != NULL)
	      Dag[row][temp[i].index]->mark = savnpath;
	  ntemp1--;
       }
    }

/*-- If there is no way to continue the path, insert it in list of paths --*/
    if (ntemp1 <= 0 && curpath->nblocks > 0)
    {
       /*   Remove this if Paths are evaluated each time a block is added*/
       best_path(curpath, paths, blocks); 		/*NOTENOTE*/
       free(temp);
       if (savpath != NULL) free(savpath);
       if (Debug) printf("\nEnd of follow_arcs for %d", row);
    }
    else   /* continue along */
    {
       nbranch = 0;           /* number of branches at this point */
       for (i=0; i<ntemp; i++)
	  if (temp[i].flag == YES &&
	      Dag[row][temp[i].index] != NULL &&
	      Dag[row][temp[i].index]->mark < savnpath &&
	      Dag[row][temp[i].index]->npos >= RSignif)
	  {
	      nbranch += 1;
	      if (nbranch > 1)             /* start a new path */
		 newpath = copypath(savpath);
	      else   /* continue along the current path, and... */
	      {      /*  save it for future branches*/
		 newpath = curpath;
		 if (ntemp > 1) savpath = copypath(curpath);
	      }
/*------------  Add new block to path & update sequences in path --- */
	      for (s=0; s<NumSeqs; s++)
		if (blocks[temp[i].index].position[s] < blocks[row].position[s])
		    newpath->seqs[s] = NO;
	      ins_bllist(temp[i].index, newpath);
/*---Add this if evaluate the path so far 
              best_path(newpath, paths, blocks);        NOTENOTE*/
/*-----------  Mark blocks already in this path and continue */
	      tblock = newpath->first_block;
	      while (tblock != NULL)
	      {
		 if (Dag[tblock->b][temp[i].index] != NULL)
		     Dag[tblock->b][temp[i].index]->mark = savnpath;
		 tblock = tblock->next_block;
	       }
	       follow_arcs(paths,newpath,temp[i].index,blocks);
	  }
       }   /*  end of else ntemp1 > 0 */

}  /*  end of follow_arcs */
/*======================================================================*/
/*  ins_bllist inserts a block at the end of the block list for a path. */
/*======================================================================*/
void ins_bllist(b, path)
int b;
struct path *path;
{
   struct block_list *bcur, *bnew;

   bnew = makebllist();
   bnew->b = b;

   bcur = path->first_block;
   if (bcur == NULL)
      path->first_block = bnew;
   else
   {
      while(bcur->next_block != NULL)
	 bcur = bcur->next_block;
      bcur->next_block = bnew;
   }
   path->nblocks += 1;

}  /*  end of ins_bllist */
/*======================================================================*/
/*   Allocate and return a new block_list structure */
/*======================================================================*/

struct block_list *makebllist()
{
   struct block_list *new;

   new = (struct block_list *) malloc(sizeof(struct block_list));
   if (new == NULL)
   {
      fprintf(stderr,"\n%.8s", Title);
      fprintf(stderr," makebllist: Unable to allocate block_list structure!\n");
      restart();
   }
   new->b = -1;
   new->minprev = 9999;  new->maxprev = -9999;
   new->next_block = new->next_best = new->prev_best = NULL;
   return(new);
}  /*  end of makebllist */
/*======================================================================*/
/*  best_paths finds the best sub-path in each path.
     Paths are "best" if the blocks are in the same order in RSignif
     sequences, and if adjacent blocks don't overlap for any sequence.
     If there are multiple adjacent overlapping blocks, the "best"
     block is the one with the highest score and with the most motifs.
     When best_paths gets the paths, the blocks are already in the
     same order in RSignif sequences, however they may overlap in some or
     all of the sequences in different ways (if two blocks overlap the
     same way in all sequences within the motif region they would have
     been merged earlier).    */
/*======================================================================*/
struct path *best_paths(paths, blocks)
struct path *paths;
struct merged_motif *blocks;
{
   struct path *cpath, *allpaths;

/*-----  Looks at all paths in  "paths" & keeps best one in "allpaths" ---*/
   allpaths = makepath();
   cpath = paths->next_path;
   while (cpath != NULL)
   {
      best_path(cpath, allpaths, blocks);
      cpath = cpath->next_path;
   }
   return(allpaths->next_path);
}  /*  end of best_paths */
/*======================================================================*/
/*   Given a new path, compute the best non-overlapping sub-path and
     compare it with the current best path, possibly replacing the best
     path with it. "paths" points to one path here, the current best path */
/*========================================================================*/
void best_path(newpath, paths, blocks)
struct path *newpath, *paths;
struct merged_motif *blocks;
{
   struct path *bestpath, *temppath;
   struct block_list *cblock, *pblock, *lblock;
   int maxdist, mindist, s, b, lastb, pb;
   int allb[MAXSEQS], thisb[MAXSEQS], saveb[MAXSEQS];
   int tots, temp, thisseq, saveseq, prop;
   double dtemp;

/*---- Note on variables:  There are three sets
      1. Current block in path =   {cblock, b, thisb[s], thisseq}
      2. Last best block in path = {lblock, lastb, saveb[s], saveseq}
	 This block is next in line to be added, once all blocks that
	 overlap with it are checked.
      3. Last block in best path = {pblock, pb, allb[s], tots}
	 This is the current state of the path     --------*/
/*-----Note on path: path->seqs[s] flags all variables in the path.
    However, since some overlapping blocks will be dropped when the best
    path is selected, more sequences than those flagged may end up in
    the best path */

   bestpath = paths->next_path;
   temppath = copypath(newpath);
   paths->nblocks += 1;		/* number of paths tested so far */

   if (Debug)
   {
      printf("\nbest_path: paths->nblocks=%d", paths->nblocks);
      tots=0;
      for (s=0; s<NumSeqs; s++)
	if (temppath->seqs[s] == YES) tots++;
	 printf("\n New path: nblocks=%d, tots=%d",
		 temppath->nblocks, tots);
   }

   cblock = temppath->first_block;             /*  Get first block */
   b = cblock->b;
   tots = saveseq = NumSeqs;
   for (s=0; s<NumSeqs; s++) 			/* all seqs to start */
	  allb[s] = saveb[s] = YES;

   if (Debug) printf(" B%d %.1s%.1s%.1s", b,
		num_to_aachar(blocks[b].aa[0]),
		num_to_aachar(blocks[b].aa[1]),
		num_to_aachar(blocks[b].aa[2]));

   lastb = b; 			/* last block in whole list */
   lblock = cblock;
   pb = -1;				/* last block in best list */
   pblock = NULL;
   cblock = cblock->next_block;
   while (cblock != NULL)
   {
       b = cblock->b;
/*----- For all sequences still in the path in which b follows lastb,
   compute the distance from the end of lastb to the beginning of b
   to see if they overlap, and see if there are enough non-overlapping
   sequences left. allb[s] are seqs in current path, which includes
   neither lastb nor b. saveb[s] are seqs with current path plus lastb,
   which hasn't been added yet.-------*/
       maxdist = -9999; mindist = 9999; thisseq = 0;
       for (s=0; s<NumSeqs; s++)
       {
	    thisb[s] = NO;
	    if (saveb[s] == YES &&
		Dag[lastb][b] != NULL &&
		blocks[b].position[s]-blocks[lastb].position[s] >= 0)
	    {
	       temp = left(s, b, blocks) - right(s, lastb, blocks) - 1;
	       if (temp < mindist) mindist = temp;
	       if (temp > maxdist) maxdist = temp;
	       if (temp >= 0) {thisb[s] = YES; thisseq++;} /*no overlap*/
	    }
	    else
	       thisb[s] = NO;
       }
       if (Debug)
       {   printf("\n (%d, %d) B%d %.1s%.1s%.1s", mindist, maxdist, b,
		num_to_aachar(blocks[b].aa[0]),
		num_to_aachar(blocks[b].aa[1]),
		num_to_aachar(blocks[b].aa[2]));
	   printf("\n   thisb: ");
	   for (s=0; s<NumSeqs; s++)
	      printf("%d=%d ", s, thisb[s]);
       }

/*---- If b never overlaps lastb, or if b overlaps lastb in few enough
   sequences, then use lastb & replace lastb with b--*/
	 if (thisseq >= RSignif)  /* Enough non-overlapping sequences */
	 {
	    temppath->nbest += 1;
	    temppath->naas += blocks[lastb].t_domain + 1;
            dtemp = sqrt(blocks[lastb].nmotif);
	    dtemp *= (double) blocks[lastb].t_score;
	    temppath->totscore += (unsigned long) dtemp;
	    temppath->totmotif += blocks[lastb].nmotif;
	    temppath->totident += blocks[lastb].nident;
	    lblock->prev_best = pblock;
	    if (pblock == NULL)
	       temppath->first_best = lblock;
	    else
	       pblock->next_best = lblock;

	    pblock = lblock; pb = lastb;
	    lblock = cblock; lastb = b;
	    tots = saveseq;  saveseq = thisseq;
	    for (s=0; s<NumSeqs; s++)
	    {
	       if (saveb[s] == NO)  allb[s] = NO;
	       saveb[s] = thisb[s];
	    }
	   if (Debug)
	   {
	      printf("\n    allb: ");
	      for (s=0; s<NumSeqs; s++)
		 printf("%d=%d ", s, allb[s]);
	   }
	 }  /* end of if still enough sequences */
/*-- If b overlaps lastb in one or more sequences currently in path,
     & has a higher score, then replace lastb with b --*/
	 else if (mindist < 0 &&
		  ( (double) (blocks[b].nmotif*blocks[b].t_score -
		       blocks[lastb].nmotif*blocks[lastb].t_score) ) > 0)
	 {
	    lblock = cblock; lastb = b;
/*  Now thisb[] is incorrect since lastb has been eliminated, so figure
    out saveb[] directly, unless this is still 1st block */
	    if (pb >= 0) 		/* Any blocks in path yet ? */
	    {
	       saveseq = 0;
	       for (s=0; s<NumSeqs; s++)
	       {
		  saveb[s] = NO;
		  if (allb[s] == YES)
		  {
		      temp = left(s,b,blocks) - right(s,pb,blocks) - 1;
		      if (Dag[pb][b] != NULL &&
			  blocks[b].position[s]-blocks[pb].position[s] >= 0 &&
			  temp >= 0)
		      { saveb[s] = YES; saveseq++;}
		   }
	       }  /* end of for s */
	    }  /* end of if pb */
	 }   /*  end of if overlap  and better score */
	 cblock = cblock->next_block;
      }
/*------------  Final block in the path ------------------------------- */
      temppath->nbest += 1;
      temppath->naas += blocks[lastb].t_domain + 1;
      dtemp = sqrt(blocks[lastb].nmotif);
      dtemp *= (double) blocks[lastb].t_score;
      temppath->totscore += (unsigned long) dtemp;
      temppath->totmotif += blocks[lastb].nmotif;
      temppath->totident += blocks[lastb].nident;
      lblock->prev_best = pblock;
      if (pblock == NULL)
	 temppath->first_best = lblock;
      else
	 pblock->next_best = lblock;
      tots = saveseq;
      for (s=0; s<NumSeqs; s++)
	 allb[s] = saveb[s];
      if (Debug)
      {
	   printf("\n   saveb, allb: ");
	   for (s=0; s<NumSeqs; s++)
	     printf("%d=%d,%d ", s, saveb[s], allb[s]);
      }

/*---------Modify path score by fraction of sequences included in it ---*/
    dtemp = (double) tots*100/NumSeqs;  /* fraction of sequences in path */
    prop = round(dtemp);
    dtemp *= (double) temppath->totscore;
    temppath->totscore = (unsigned long) dtemp;

/*---------------------------------------------------------------------*/
   if (Debug)
   {
      printf("\n ");
      printf("Best sub-path:%d blocks, %d AAs, total score=%ld",
	       temppath->nbest, temppath->naas, temppath->totscore);
      printf(",\n total motifs=%d, total idents=%d",
	       temppath->totmotif, temppath->totident);
      printf(", total sequences=%d, prop=%ld", tots, prop);
   }

/*---------------------------Update best path ---------------------------*/
      if (bestpath == NULL ||
	  temppath->totscore > bestpath->totscore ||
	  (temppath->totscore == bestpath->totscore &&
	   temppath->totmotif > bestpath->totmotif) )
      {
	  paths->next_path = temppath;      /* New best path */
	  free_path(bestpath);  /*  Throw away old best path */
	  if (Debug) printf("\n >>>NEW BESTPATH");
	  /*---- update info for best path ----*/
	  temppath->nseqs = tots;
	  for (s=0; s<NumSeqs; s++)
	      temppath->seqs[s] = allb[s];
	  /*---- first block in best path, distance from start of seq ---*/
	  lblock=temppath->first_best;
	  mindist=9999; maxdist=-9999;
	  for (s=0; s<NumSeqs; s++)
	     if (temppath->seqs[s] == YES)
	     {
		   temp=left(s,lblock->b,blocks);
		   if (temp <mindist) mindist=temp;
		   if (temp >maxdist) maxdist=temp;
	     }
	  lblock->minprev=mindist; lblock->maxprev=maxdist;
	  /*---- rest of the blocks in best path ----*/
	  cblock=lblock->next_best;
	  while (cblock != NULL)
	  {
	     mindist=9999; maxdist=-9999;
	     for (s=0; s<NumSeqs; s++)
		if (temppath->seqs[s] == YES)
		{
		   temp=left(s,cblock->b,blocks)-right(s,lblock->b,blocks)-1;
		   if (temp <mindist) mindist=temp;
		   if (temp >maxdist) maxdist=temp;
		}
	     cblock->minprev=mindist; cblock->maxprev=maxdist;
	     lblock = cblock;
	     cblock=cblock->next_best;
	  }
      }
      else
	  free_path(temppath);  /* Throw away new path */

}  /*  end of best_path */
/*=======================================================================
   If there are DUPS, take a look at the sequences left out of the best
   path and see if they can be added. They can be if they have all the
   blocks in some non-overlapping order, not necessarily the same order
   as in the best path sequences.
======================================================================*/
void check_seqs(path, blocks)
struct path *path;
struct merged_motif *blocks;
{
   int s, b, lastb, ntemp, i, overlap;
   struct block_list *cblock;
   struct temp *temp;

   temp = (struct temp *) malloc(Total_Motifs*sizeof(struct temp));
   if (temp == NULL)
   {
      fprintf(stderr,"\n%.8s", Title);
      fprintf(stderr," check_seqs: Unable to allocate temp structure!\n");
      restart();
   }
   for (s=0; s<NumSeqs; s++)
   {
      if (path->seqs[s] == NO)
      {
	ntemp = 0;
	cblock = path->first_best;
	while (cblock != NULL)
	{
	   temp[ntemp].index = cblock->b;
	   temp[ntemp].value = blocks[cblock->b].position[s];
	   temp[ntemp++].flag = 0;
	   cblock = cblock->next_best;
	}
/*------- Sort blocks in best path by their order in the sequence, which
     is probably different than their order in the best path.
     Then check to see whether they overlap in the sequence. If they
     don't, then add the sequence to the best path ----*/
	qsort(temp, ntemp, sizeof(struct temp), tempcmp);
	lastb = temp[0].index; overlap = NO;
	for (i=1; i<ntemp; i++)
	{
	   b = temp[i].index;
	   if ( (left(s, b, blocks) - right(s, lastb, blocks) - 1) < 0 )
		overlap = YES;
	   lastb = b;
	}
	if (overlap == NO)
	{ path->seqs[s] = YES; path->nseqs += 1; }
      }  /*  end of s == NO */
   }

   free(temp);
}  /* end of check_seqs */
/*====================================================================*/
/*    print_best prints the best path.                                */
/*====================================================================*/
void print_best(path, blocks, seqs)
struct path *path;
struct merged_motif *blocks;
struct sequences *seqs;
{
  struct block_list *cblock;
  int i, s, b, nseq, nb, allb[MAXSEQS], plt_cursor;
  char ac[MAXLINE], id[MAXLINE], de[MAXLINE], temp[MAXLINE], *ptr;
  FILE *plt;

/*------ Assume Title has form ">ACxxxxx ;ID;Description;" ----------*/
  if (strlen(Title) < 10)                   /* Be sure! */
  {
     strcpy(temp, ">ACxxxxx ;none;");              /*  AC and ID */
     strcat(temp, Title); strcat(temp, ";");       /* DE */
  }
  else  strcpy(temp, Title);
  Title[0] = id[0] = ac[0] = de[0] = '\0';
  if (temp[0] != '>') temp[0] = '>';
  ptr = strtok(temp, " ;"); 
  if (ptr != NULL)
  {
      strcpy(ac, ptr);                           /*  AC  */
      ptr = strtok(NULL, ";");
      if (ptr != NULL)
      {
          strcpy(id, ptr);                       /*  ID  */
          ptr = strtok(NULL, ";");               /*  DE  */
          if (ptr != NULL) strcpy(de, ptr);
      }
  }
  /*-------  stretch out or shorten the AC to make it exactly long enough */
  for (i=strlen(ac); i<9; i++) ac[i] = 'x';
  ac[9] = '\0';
  if (!strlen(id)) strcpy(id, "none"); 
  if (strlen(id) > 21) id[21] = '\0';
  if (!strlen(de)) strcpy(de, "none"); 
  if (strlen(de) > 61) de[61] = '\0';

/*------------Write out blocks------------------------------------------*/
/*--   Open the plotting file ---*/
  temp[0] = '\0';
/*  s = strcspn(Mot_Filename, ".");
  strncat(temp, Mot_Filename, s)
*/
  /* 5/11/92  changed .plt file name */
  strcpy(temp, Blk_Filename); strcat(temp, ".plt");
  if ( (plt=fopen(temp, "w+t")) == NULL)
     printf("\nCannot open file %s, no plotting output.\n", temp);
  plt_cursor = 0;			/*  x-axis plotting value */
/*--   Problem with passing path->seqs to save_block.... */
  for (s=0; s<NumSeqs; s++)
     allb[s] = path->seqs[s];

  printf("\nBest path has %d sequences out of %d:", path->nseqs, NumSeqs);
  nb = 0;
  cblock = path->first_best;
  while (cblock != NULL)
  {
     b = cblock->b;
     printf("\n (%d, %d) B%d %.1s%.1s%.1s", cblock->minprev, cblock->maxprev,
		cblock->b,
		num_to_aachar(blocks[cblock->b].aa[0]),
		num_to_aachar(blocks[cblock->b].aa[1]),
		num_to_aachar(blocks[cblock->b].aa[2]));
     nb += 1;
				     /*---  Add path info to title ----*/
     if (nb == 1 && path->nbest == 1) ac[8] = '\0';  /* no A if only 1 blk */
     else if (nb == 1 && path->nbest > 1) ac[8] = 'A';
     else if (nb == 2)  ac[8] = 'B';
     else if (nb == 3)  ac[8] = 'C';
     else if (nb == 4)  ac[8] = 'D';
     else if (nb == 5)  ac[8] = 'E';
     else if (nb == 6)  ac[8] = 'F';
     else if (nb == 7)  ac[8] = 'G';
     else if (nb == 8)  ac[8] = 'H';
     else if (nb == 9)  ac[8] = 'I';
     else if (nb == 10) ac[8] = 'J';
     else if (nb == 11) ac[8] = 'K';
     else if (nb == 12) ac[8] = 'L';
     else if (nb == 13) ac[8] = 'M';
     else if (nb == 14) ac[8] = 'N';
     else if (nb == 15) ac[8] = 'O';
     else if (nb == 16) ac[8] = 'P';
     else if (nb == 17) ac[8] = 'Q';
     else if (nb == 18) ac[8] = 'R';
     else if (nb == 19) ac[8] = 'S';
     else if (nb == 20) ac[8] = 'T';
     else if (nb == 21) ac[8] = 'U';
     else if (nb == 22) ac[8] = 'V';
     else if (nb == 23) ac[8] = 'W';
     else if (nb == 24) ac[8] = 'X';
     else if (nb == 25) ac[8] = 'Y';
     else if (nb > 25)  ac[8] = 'Z';
     sprintf(Title, "%s (%d,%d);%s;%s",
	     ac, cblock->minprev, cblock->maxprev, id, de);
     save_block(&blocks[b], seqs, allb);
     if (plt != NULL)
  save_plot(&plt_cursor, cblock->minprev, cblock->maxprev, plt, &blocks[b]);
     cblock = cblock->next_best;
  }
  if (plt != NULL) fclose(plt);

  /*----------Summary of the best path --------------------------------*/
  printf("\n%d blocks, %d AAs, total score=%ld, total motifs=%d",
	       path->nbest, path->naas, path->totscore, path->totmotif);
  printf(", total conserved=%d", path->totident);
  printf("\n total sequences=%d out of %d", path->nseqs, NumSeqs);

  /*------Print seqs in path----------------------------------------------*/
  printf("\nSequences in the best path:");
  nseq = 0;
  for (s=0; s<NumSeqs; s++)
  {
    if (path->seqs[s] == YES)
    {
	if (nseq%5 == 0) printf("\n");
	nseq++;
	printf("%.3d:%10s ", s, seqs->name+SNAMELEN*s);
    }
  }

  /*------Print seqs NOT in path----------------------------------------*/
  printf("\nSequences not in the best path:");
  nseq = 0;
  for (s=0; s<NumSeqs; s++)
  {
    if (path->seqs[s] == NO)
    {
	if (nseq%5 == 0) printf("\n");
	nseq++;
	printf("%.3d:%10s ", s, seqs->name+SNAMELEN*s);
    }
  }


}  /*  end of print_best */
/*==================================================================*/
/*      print a path  */
/*==================================================================*/
void print_path(path)
struct path *path;
{
  struct block_list *cblock;

  printf("\nPath:");
  cblock = path->first_block;
  while (cblock != NULL)
  {
     printf(" %d", cblock->b);
     cblock = cblock->next_block;
  }
}  /*  end of print_path */
/*=======================================================================
    save_plot adds plotting information for a block to the plot file
     and updates the plotting cursor, which is the x-axis value.
=========================================================================*/
void save_plot(cursor, minprev, maxprev, plt, block)
int *cursor, minprev, maxprev;
FILE *plt;
struct merged_motif *block;
{
   int c, temp;

   /*-----  Plot zeros or ones between blocks ---------*/
   for (c=0; c<maxprev; c++)
   {
      if (c < minprev) temp = 0;
      else             temp = 1;
      if (plt != NULL) fprintf(plt, "%d  %d\n", *cursor, temp);
      *cursor = *cursor + 1;
   }
   /*------  Plot the block scores -------------------*/
   temp = block->loffset - block->t_loffset;		/* offset */
   for (c=0; c <= block->t_domain; c++)
   {
      if (plt != NULL)
	 fprintf(plt, "%d  %d\n", *cursor, block->scores[c+temp]);
      *cursor = *cursor + 1;
   }
}  /* end of save_plot */
/*====================================================================*/
/*   Return the left-hand block position */
/*======================================================================*/

int left(s, b, blocks)
int s, b;
struct merged_motif *blocks;
{
   return(blocks[b].leftpos[s] - blocks[b].t_loffset);
}  /*  end of left */
/*====================================================================*/
/*  Return the right-hand block position */
/*======================================================================*/

int right(s, b, blocks)
int s, b;
struct merged_motif *blocks;
{
   return(left(s, b, blocks) + blocks[b].t_domain);
}  /*  end of right */
