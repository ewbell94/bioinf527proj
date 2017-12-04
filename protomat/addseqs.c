/*    COPYRIGHT 1998 Fred Hutchinson Cancer Research Center
	Requires BLIMPS 3.2 libraries to compile
        Requires environment variable BLIMPS_DIR be set               */
/*    addseqs.c   Reads a file of blocks and a file of sequences.
		  Assumes all of the input blocks were made from
		  the input sequences, so that the offsets in the
		  blocks are correct wrt the sequence file.
		  Adds sequences not in the blocks to the blocks.
		  Only adds a sequence if the blocks are in order
		  in all the blocks (ignores "dups") and if the
		  alignment score of all of the blocks is >= CUTOFF.

           addseqs <in_blocks_file> <seqs_file> <out_blocks_file>

--------------------------------------------------------------------
 6/19/97 J. Henikoff
 6/23/97 Added theoretical calibration (pssmdist())
 6/26/97 Change segments scoring below CUTOFF to Xs 
         Added order check, but it's just left to right, should be
	 centered around the highest scoring block, like blksort.
 1/26/98 Check best alignment for all blocks before adding a sequence,
	 and only add it if all blocks score above CUTOFF and are in
	 the correct order. Changed CUTOFF from 1100 to 800.
 1/28/98 Update the distance from previous blocks.
 1/29/98 Puts added sequences in an additional cluster.
         Writes out addseqs.dat file.
 2/ 7/98 Fixed double AC from fix_ac()
====================================================================*/

#define ADDSEQS_C_
#define EXTERN
#define MAXNAME 80	/* Maximum file name length */
#define MAXWIDTH 60     /* Max block width */
#define MAXAA 26	/* alphabet size */
#define AASALL 26	/* alphabet size */
#define AAS 21	
#define MAXWIDTH 60
#define CUTOFF 800	/* Min. calibrated cutoff score */

/*		Stuff for pssmdist()   */
#define TPS 100			/* number of known TPs */
#define Search 185371		/* number of seqs in database */
#define SearchAA 58639837	/* number of aas in database */
#define MAXSCORE 6000

#include "blocksprogs.h"

struct block_list {
   int nblock;		/* number of blocks in family */
   int nseq;		/* number of seqs originally in blocks */
   int nadd;		/* number of seqs added to blocks */
   int mindist;		/* min dist from last block */
   int maxdist;		/* max dist from last block */
   int bestpos;		/* for current sequence */
   int bestscore;	/* for current sequence */
   Block *block;
   Matrix *pssm, *pssm_frq;
   struct block_list *next;
};

void addseq();
void best_pos();
void scale_weights();
void fix_ac();
void add_cluster();
/*		List of blocks routines  */
struct block_list *make_blist();
void insert_blist();
void free_blist();
/*		Block calibration routine */
struct score_struct {
	double ways, prob;
};
int pssmdist();

int NumSeqs, TN995, TPabove;

/*=======================================================================*/
/*
 * main
 */

int main(argc, argv)
     int argc;
     char *argv[];
{
   FILE *fseq, *fblk, *fout, *fqij, *fdat;
   Block *block, *newblock;
   Sequence *sequence;
   struct block_list *blist, *bcur, *bprev;
   char frqname[MAXNAME], qijname[MAXNAME];
   char infile[MAXNAME], outfile[MAXNAME], seqsfile[MAXNAME];
   char *ptr, *blimps_dir;
   int i, j, not_in, in_order, above_cut, width, pos, endpos;

   ErrorLevelReport = 2;	/* 5 for No blimps errors */

   if (argc < 3)
   {
      printf("ADDSEQS: Copyright 1997 by the Fred Hutchinson Cancer");
      printf(" Research Center\n");
      printf("Adds sequences to blocks.\n");
      printf("USAGE:  addseqs <in_block_file> <seqs_file> <out_blocks_file>\n");
      printf("                <in_block_file> = input blocks\n");
      printf("		      <seqs_file>     = fasta format sequences\n");
      printf("                <out_blocks_file> = output blocks\n");
   }
/*----------------Get the input blocks------------------------------*/
   if (argc > 1)
      strcpy(infile, argv[1]);
   else
   {
      printf("\nEnter name of file containing blocks: ");
      gets(infile);
   }
   if ( (fblk=fopen(infile, "r")) == NULL)
   {
      printf("\nCannot open file %s\n", infile);
      exit(-1);
   }

/* ------------2nd arg = sequence file -----------------------*/
   if (argc > 2)
      strcpy(seqsfile, argv[2]);
   else
   {
      printf("\nEnter name of sequence file: ");
      gets(seqsfile);
   }
   if ( (fseq=fopen(seqsfile, "r")) == NULL)
   {
      printf("\nCannot open file %s\n", seqsfile);
      exit(-1);
   }
   /* ------------3rd arg = output blocks -----------------------*/
   if (argc > 3)
      strcpy(outfile, argv[3]);
   else
   {
      printf("\nEnter name of output blocks file: ");
      gets(outfile);
   }
   if ( (fout=fopen(outfile, "w")) == NULL)
   {
      printf("\nCannot open file %s\n", outfile);
      exit(-1);
   }

   /*--------- Get the information to make PSSMs ---------------------*/
   blimps_dir = getenv("BLIMPS_DIR");
   frqname[0] = '\0';
   if (blimps_dir != NULL) sprintf(frqname, "%s/docs/", blimps_dir);
   strcat(frqname, "default.amino.frq");
   /*  Creates global array frequency[]  */
   load_frequencies(frqname);
   qijname[0] = '\0';
   if (blimps_dir != NULL) sprintf(qijname, "%s/docs/", blimps_dir);
   strcat(qijname, "default.qij");
   Qij = NULL;
   if ( (fqij=fopen(qijname, "r")) != NULL) Qij = load_qij(fqij);
   fclose(fqij);
   RTot = LOCAL_QIJ_RTOT;

   /*=======================================================================*/

   /*  Read all of the blocks into memory  */
   blist = make_blist();
   while ((block = read_a_block(fblk)) != NULL)
   {
       pb_weights(block);		  /* add position-based weights */
       scale_weights(block, 1);	           /* weights add to # of seqs */
       insert_blist(blist, block);	 /* also computes PSSMs & calibrates */
   }
   fclose(fblk);
   if (blist->nblock == 0)
   {
      printf("No blocks found in %s .\n", infile);
      exit(-1);
   }

   /*======================================================================*/
   /*---- Process each block for every sequence ---------------------------*/
   NumSeqs = 0;
   while (!feof(fseq) &&
      (sequence = read_a_sequence(fseq, FASTA, AA_SEQ)) != NULL)
   {
     NumSeqs += 1;
     /*---------------Process each block---------------------------------*/
     /*   Assumes all the blocks are from the same family!  */
     bcur = blist->next;
     endpos = -99;	/* first block is never out of order */
     not_in = in_order = above_cut = TRUE;
     while (bcur != NULL && bcur->block != NULL)
     {
         for (i = 0; i < bcur->block->num_sequences; i++)
         {
            j = (int) strlen(bcur->block->sequences[i].name);
            if ((int) strlen(sequence->name) < j) j = strlen(sequence->name);
          if ((strncmp(sequence->name, bcur->block->sequences[i].name, j)) == 0)
                not_in = FALSE;
         }
         if (not_in)
         {
           /* Get best position & score for this block in this sequence */
	   bcur->bestpos = bcur->bestscore = -999; 
           best_pos(bcur, sequence);
           if (bcur->bestscore < CUTOFF) above_cut = FALSE;
           if (bcur->bestpos < endpos)
           {
               printf(" out of order %d, previous end %d\n", 
		       bcur->bestpos, endpos);
               in_order = FALSE;
           }
           else
           {  endpos = bcur->bestpos + bcur->block->sequence_length; }
         }
         bcur = bcur->next;
     }  /* end of block  */
     if (not_in && in_order && above_cut)
     {
        bcur = blist->next; bprev = NULL;
        while (bcur != NULL && bcur->block != NULL)
        {
            addseq(bprev, bcur, sequence);
            if (!blist->nadd) add_cluster(bcur->block);
	    bprev = bcur;
            bcur = bcur->next;
        }
        blist->nadd += 1;
     }
     free_sequence(sequence);
   }  /* end of sequence */
   if (NumSeqs == 0)
   {
      printf("No sequences found in %s .\n", seqsfile);
      exit(-1);
   }
   fclose(fseq);
   printf("\n");

   /* Add sequence weights & output the final blocks */
   bcur = blist->next;
   while (bcur != NULL && bcur->block != NULL)
   {
      pb_weights(bcur->block);		/* add position-based weights */
      scale_weights(bcur->block, 0);    /* max weight = 100 */
      fix_ac(bcur);
      output_block(bcur->block, fout);
      bcur = bcur->next;
   }  /* end of block  */
   fclose(fout);

   /*----------Append statistics to addseqs.dat----------------------*/
   if ( (fdat=fopen("addseqs.dat", "a")) != NULL)
   {
      fprintf(fdat, "%s %d %d %d %d\n",
              blist->next->block->number, blist->nblock, blist->nseq, 
              blist->next->block->num_sequences, NumSeqs);
      fclose(fdat);
   }

   printf("%d sequences added\n", blist->nadd);
   exit(0);

}  /* end of main */

/*=======================================================================
	Add sequence to the end of the blist block
========================================================================*/
void addseq(bprev, blist, seq)
struct block_list *bprev, *blist;
Sequence *seq;
{
   int j, pos, posj, prevdist, clmax, newseq;

      printf("%s: ", blist->block->number);
      pos = blist->bestpos;

      printf("Adding %s at position %4d\n", seq->name, pos);

      if (blist->block->num_sequences + 1 > blist->block->max_sequences)
      {
         resize_block_sequences(blist->block);
         /*  resize() doesn't fix cluster pointers?  */
/*
         blist->block->clusters[0].sequences = &(blist->block->sequences[0]);
*/
      }

      newseq = blist->block->num_sequences;
      /*  Block points to first position */
      blist->block->sequences[newseq].position = pos + 1;
      strcpy(blist->block->sequences[newseq].name, seq->name);
      strcpy(blist->block->sequences[newseq].info, seq->info);
      blist->block->sequences[newseq].max_length =
             blist->block->sequence_length;
/*    blist->block->sequences[newseq].length =
             blist->block->sequence_length; 
*/
      blist->block->num_sequences += 1;

      clmax = blist->block->num_clusters - 1;
      blist->block->clusters[clmax].num_sequences += 1; 

      /*  Force all sequences into a single cluster */
/*
      blist->block->num_clusters = 1;
      blist->block->clusters[0].num_sequences = blist->block->num_sequences;
      blist->block->clusters[0].sequences = &(blist->block->sequences[0]);
*/

      for (j=0; j< blist->block->sequence_length; j++)
      {
         if (blist->bestscore >= CUTOFF)
         {
            posj = pos + j;
            if (posj >= 0 && posj < seq->length)
               blist->block->sequences[newseq].sequence[j] =
                      seq->sequence[posj];
            else
               blist->block->sequences[newseq].sequence[j] = 23;    /*   X */
         }
         else
         {
            blist->block->sequences[newseq].sequence[j] = 23;    /*   X */
         }
      }

      /*   Update distance from previous block if necessary  */
      prevdist = pos;
      if (bprev != NULL)
      {
         prevdist -= (bprev->block->sequences[newseq].position + bprev->block->sequence_length);
      }
      if (prevdist < blist->mindist) blist->mindist = prevdist;
      if (prevdist > blist->maxdist) blist->maxdist = prevdist;
    
}  /* end of addseq */

/*======================================================================
      Scale sequence weights so maximum value is 100
=======================================================================*/
void scale_weights(block, stype)
Block *block;
int stype;
{
   double maxweight, minweight, sumweight, factor;
   int seq;

   sumweight = maxweight = 0.0; minweight = 999999.9;
   for (seq = 0; seq < block->num_sequences; seq++)
   {
      sumweight += block->sequences[seq].weight;
      if (block->sequences[seq].weight > maxweight) 
           maxweight = block->sequences[seq].weight;
      if (block->sequences[seq].weight < minweight) 
           minweight = block->sequences[seq].weight;
   }

   factor = 1.0;
   /*    Force maximum weight to be 100  */
   if (stype == 0)  factor = 100. / maxweight;
   /*    Force sum of weights to be number of sequences */
   else if (stype == 1) factor = (double) block->num_sequences / sumweight;
   /*    Force sum of weights to be stype */
        else if (stype > 1) factor = (double) stype / sumweight;

   for (seq = 0; seq < block->num_sequences; seq++)
   {
      block->sequences[seq].weight *= factor;
   }
}  /* end of scale_weights */

/*=======================================================================
            Score all alignments of matrix with sequence & return
	     position of best alignment
============================================================================*/
void best_pos(blist, seq)
struct block_list *blist;
Sequence *seq;
{
  int scan_pos, seq_pos, mat_pos, best_position;
  double seq_score, max_seq_score;

  max_seq_score = 0.0;
  
  /* for every alignment of the matrix and sequence score the alignment. */
  /* there are (seq_length + matrix_length  - 1) different alignments. */
  /* Note: seq_pos is the relative position of the 1st matrix column to the */
  /*       1st sequence column */
  /* Note: the indexing is shifted to make the calculation of scan_pos */
  /*       (see below) easier/faster */
  for (seq_pos= -blist->pssm->length+1; seq_pos < seq->length; seq_pos++) {


    /* 
     * score this alignment
     */
    seq_score = 0.0;

    /* for each position in the matrix add the weight to the total score */
    /* Note: mat_pos is the current column of the matrix */
    for (mat_pos=0; mat_pos < blist->pssm->length; mat_pos++) {

      /* make sure the sequence and the matrix overlap at this point */
      /* Note: scan_pos is where the current matrix column is in the */
      /*       sequence */
      scan_pos = seq_pos + mat_pos;
      if ((scan_pos >= 0) && (scan_pos < seq->length)) { /* if in the seq */
	seq_score += blist->pssm->weights[seq->sequence[scan_pos]][mat_pos];
      }
      else {			/* not in the sequence */
	seq_score += blist->pssm->weights[aa_atob['-']][mat_pos]; 
	/* score as if it was a gap */
      }

    } /* end score this alignment */

    
      if (seq_score > max_seq_score) {
	max_seq_score = seq_score;
	best_position = seq_pos;
      }
    
  } /* end of pairwise lineup */
  
  /*  Normalize score */
  if (blist->pssm->percentile > 0)
  {
     printf("%s %s calibrated ", seq->name, blist->block->number);
     max_seq_score *= 1000.0;
     max_seq_score /= (double) blist->pssm->percentile;
  }
  blist->bestpos = best_position;
  blist->bestscore = round(max_seq_score);
  printf("best score = %5d at %d\n", blist->bestscore, blist->bestpos);

}  /* end of best_pos */

/*=======================================================================
     routines for a list of blocks
========================================================================*/
struct block_list *make_blist()
{
   struct block_list *new;
   
   new = (struct block_list *) malloc (sizeof(struct block_list));
   new->nblock = new->nseq = new->nadd = 0;
   new->bestscore = new->bestpos = new->maxdist = -9999;
   new->mindist = 9999;
   new->block = NULL;
   new->next = NULL;

   return(new);
}  /* end of make_blist */

/*=======================================================================
  PSSMs are computed & block calibrated here. The 99.5
  and strength values are not separate fields in the block record, but
  are extra fields stored in the PSSM record.
=======================================================================*/
void insert_blist(blist, block)
struct block_list *blist;
Block *block;
{
   struct block_list *cur;
   char *ptr, ctemp[80];

   cur = blist;
   while (cur->next != NULL) cur = cur->next;
   cur->next = make_blist();
   cur->next->block = block;

   cur->next->pssm = block_to_matrix(block, 3);		/* log odds  PSSM */
   cur->next->pssm_frq = block_to_matrix(block, 2);	/* frequency PSSM */

   /*   Get min & max distances:
	block->ac = "PS00094A; distance from previous block=(0,1141)"  */
   strcpy(ctemp, block->ac);
   ptr = strtok(ctemp, "(");
   if (ptr != NULL)
   { 
      ptr = strtok(NULL, ",");
      if (ptr != NULL)
      {
          cur->next->mindist = atoi(ptr);
          ptr = strtok(NULL, ")");
          if (ptr != NULL)
          { cur->next->maxdist = atoi(ptr);}
      }
   }

   /*  Calibrate the block if it isn't already & add to block->bl line */
   if (cur->next->pssm->percentile <= 0)
   {
      cur->next->pssm->percentile = pssmdist(cur->next->pssm, 0, frequency, NULL);
      block->percentile = cur->next->pssm->percentile;
/*
      ptr = strtok(block->bl, "\n\r");
      sprintf(block->bl, "%s; 99.5%%=%d", ptr, cur->next->pssm->percentile);
*/
   }
   if (cur->next->pssm->strength <= 0)
   {
      cur->next->pssm->strength = pssmdist(cur->next->pssm, 1, frequency, cur->next->pssm_frq);
      block->strength = cur->next->pssm->strength;
/*
      ptr = strtok(block->bl, "\n\r");
      sprintf(block->bl, "%s; strength=%d", ptr, cur->next->pssm->strength);
*/
   }

   /*  Update the list totals */
   blist->nblock += 1;
   blist->nseq = block->num_sequences;
 
}  /* end of insert_blist */

/*=======================================================================*/
void free_blist(blist)
struct block_list *blist;
{
   struct block_list *cur, *last;
   cur = last = blist;
   while (cur->next != NULL)
   {
      last = cur;  cur = cur->next;
   }
   if (cur != blist)
   {
      free(cur);
      last->next = NULL;
      free_blist(last);
   }
   else free(blist);

}  /* end of free_blist */

/*=======================================================================
    ftype == 0 => use freqs (TN distribution)
    ftype == 1 => use obs_freqs (TP distribution)
========================================================================*/
int pssmdist(matrix, ftype, freqs, obs_freqs)
Matrix *matrix, *obs_freqs;
int ftype;
double *freqs;
{
  struct score_struct *last, *this, ends[MAXSCORE];
  struct score_struct middle[MAXSCORE], scores[2][MAXSCORE];
  int col, aa, minvalue, maxvalue, minscore, maxscore, mincol, maxcol;
  int x, score, minfirst, minlast, pflag;
  double cum, aligns, report, ftemp, dtemp, probwt[MAXWIDTH];

  /*   compute which score to report */
  if (ftype == 1) report = (double) 0.5 * TPS;		/* median TP */
  else            report = (double) 0.005 * Search;	/* 99.5% TN  */
  minscore = maxscore = 0;
  maxvalue = -1;		/* assumes no negative scores */
  minvalue = 9999;
  probwt[0] = 20.0;
  for (col = 0; col < matrix->length; col++)
  {
     /*  Intialize probability weights with powers of 20 */
     if (col > 0) probwt[col] = probwt[col - 1] * 20.0;
     mincol = 9999; maxcol = -1;
/*   for (aa = 0; aa < AASALL; aa++) */
     for (aa = 1; aa <= 20; aa++)
     {
        if (matrix->weights[aa][col] > maxvalue)
            maxvalue = matrix->weights[aa][col];
	if (matrix->weights[aa][col] < minvalue)
	    minvalue = matrix->weights[aa][col];
        if (matrix->weights[aa][col] > maxcol)
            maxcol = matrix->weights[aa][col];
	if (matrix->weights[aa][col] < mincol)
	    mincol = matrix->weights[aa][col];
     }
     maxscore += maxcol;
     minscore += mincol;
     if (col == 0) minfirst = mincol;
     if (col == (matrix->length - 1) ) minlast = mincol;
   }
   /*    Probability weights: weights applied to probs depending
         on how many columns of the block are aligned, the
         full alignment width gets about 19/21 = .905, alignments
         off either end of the sequence get the remainder */
   dtemp = probwt[matrix->length - 1] * 21.0 - 40.0;
   cum = 0.0;
   for (col = 0; col < matrix->length; col++)
   {
      probwt[col] = probwt[col] * 19.0 / dtemp;
      cum += probwt[col];
      if (col < matrix->length - 1) cum += probwt[col];
   }
   /* ---- Min. score could be just first or last column ----*/
   if (minfirst < minscore) minscore = minfirst;
   if (minlast < minscore)  minscore = minlast;
/*
   if (ftype == 1) printf("\nTrue positive distribution\n");
   else            printf("\nRandom distribution\n");
   printf("minscore without ends=%d\n", minscore);
   printf("minvalue=%d maxvalue=%d minscore=%d maxscore=%d\n",
	   minvalue, maxvalue, minscore, maxscore);
*/
   /*-----------------------------------------------------------------*/
   if (maxscore > MAXSCORE) maxscore = MAXSCORE - 1;
   last = scores[0]; this = scores[1];
   for (x = minvalue; x <= maxscore; x++)
   {  last[x].ways = last[x].prob = this[x].ways = this[x].prob = 0.0; 
      ends[x].ways = ends[x].prob = middle[x].ways = middle[x].prob = 0.0; }

   /*---------Initialize from first column -------------------------------*/
   /*   use obs_freqs->weights[aa][col] instead of freqs[aa] for TP probs */
   col = 0;
/* for (aa=0; aa < AASALL; aa++) */
   for (aa=1; aa <= 20; aa++)
   {
      x = matrix->weights[aa][col];
      last[x].ways += 1.0;
      if (ftype == 1) ftemp = obs_freqs->weights[aa][col] / 100.;
      else            ftemp = freqs[aa];
      last[x].prob += ftemp;   
   }

   /*---- Now enumerate all possible scores, expanding one column
          at a time ----------*/
   for (col=1; col < matrix->length; col++)
   {
      /*--------- Save the alignments hanging off the left end ------*/
      /*    There are currently col+1 columns of the block aligned  */
      for (x=minvalue; x <= maxscore; x++)
      {
            ends[x].ways += last[x].ways;
            ends[x].prob += last[x].prob * probwt[col - 1];
      }
/*    for (aa=0; aa < AASALL; aa++) */
      for (aa=1; aa <= 20; aa++)
      {
         for (x=minvalue; x <= maxscore; x++)
         {
            if (last[x].ways > 0)
            {
               score = x + matrix->weights[aa][col];
               this[score].ways += last[x].ways;
               if (ftype == 1) ftemp = obs_freqs->weights[aa][col] / 100.;
               else            ftemp = freqs[aa];
               this[score].prob += last[x].prob * ftemp;
            }
         } /* end of score x */
      }  /* end of aa */
      /*---------   Switch the arrays ------------------------------*/
      if (this == scores[1])
      {  last = scores[1]; this = scores[0]; }
      else
      {  last = scores[0]; this = scores[1]; }
      for (x = minvalue; x <= maxscore; x++)
      {  this[x].ways = this[x].prob = 0.0;  }
   }  /* end of col */

   /*-------- last now has the final counts of all combinations
      of column scores for full blocks, and ends has the counts
      for alignments off the left end. Still have to get counts
      for alignments off the right end and need two arrays to do
      it. So have to keep last results in another array--------*/
   for (x=minvalue; x <= maxscore; x++)
   {
	middle[x].ways = last[x].ways;
	middle[x].prob = last[x].prob * probwt[matrix->length - 1];
	last[x].ways = last[x].prob = 0.0;
   }

   /*-------- Get the alignments hanging off the right end -------*/
   /*    Initialize with last column   */
   col = matrix->length - 1;
   for (aa=1; aa <= 20; aa++)
   {
      x = matrix->weights[aa][col];
      last[x].ways += 1.0;
      if (ftype == 1) ftemp = obs_freqs->weights[aa][col] / 100.;
      else            ftemp = freqs[aa];
      last[x].prob += ftemp;   
   }
   for (col = matrix->length - 2; col >= 1; col--)
   {
      /*--------- Save the gapped alignments off the right end ------*/
      /*  There are currently length-col columns of the block aligned */
      for (x=minvalue; x <= maxscore; x++)
      {
            ends[x].ways += last[x].ways;
            ends[x].prob += last[x].prob * probwt[matrix->length - col - 2];
      }
      for (aa=1; aa <= 20; aa++)
      {
         for (x=minvalue; x <= maxscore; x++)
         {
            if (last[x].ways > 0)
            {
               score = x + matrix->weights[aa][col];
               this[score].ways += last[x].ways;
               if (ftype == 1) ftemp = obs_freqs->weights[aa][col] / 100.;
               else            ftemp = freqs[aa];
               this[score].prob += last[x].prob * ftemp;
            }
         } /* end of score x */
      }  /* end of aa */
      /*---------   Switch the arrays ------------------------------*/
      if (this == scores[1])
      {  last = scores[1]; this = scores[0]; }
      else
      {  last = scores[0]; this = scores[1]; }
      for (x = minvalue; x <= maxscore; x++)
      {  this[x].ways = this[x].prob = 0.0;  }
   }  /* end of column */
   /*--------- Save the gapped alignments off the right end ------*/
   /*  Need to get the length - 1 counts from the right   */
   for (x=minvalue; x <= maxscore; x++)
   {
         ends[x].ways += last[x].ways;
         ends[x].prob += last[x].prob * probwt[matrix->length - 2];
   } 

   /*   number of alignments done in a hypothetical search   */
   /*   Always assume one alignment per TP sequence          */
   if (ftype == 1) aligns = TPS;
   else aligns = (double) SearchAA + Search * (matrix->length - 1); 

   /* Count the TP scores above 99.5% of TN scores, x is the score */
   if (ftype == 1 && TN995 > 0)
   {
      x = maxscore; TPabove = 0.0;
      while (x > TN995)
      {
         if (middle[x].prob > 0.0 || ends[x].prob > 0.0)
         {
            TPabove += (middle[x].prob + ends[x].prob) * aligns;
         }
         x--;
      }
   }

   /* Print the top scores, x is the score */
   x = maxscore; cum = 0.0;
/*
   printf("\nCum      Score   Expected\n");
*/
   while (cum < report && x >= 0)
   {
       if (middle[x].prob > 0.0 || ends[x].prob > 0.0)
       {
/*
          if (cum > 100.0)
             printf("%8.4f %4d %8.4f\n", 
		cum, x, (middle[x].prob + ends[x].prob) * aligns);
*/
          cum += (middle[x].prob + ends[x].prob) * aligns;
       }
       x--;
   }
   if (ftype == 1) printf("Median TP score = %d, TPabove = %lf\n",
                           x, TPabove);
   else            printf("99.5 TN score = %d\n", x);
   return(x);
}  /* end of pssmdist */
/*================================================================
    fix_ac() updates distances on the block->ac line
	block->ac = "PS00094A; distance from previous block=(0,1141)"
================================================================*/
void fix_ac(blist)
struct block_list *blist;
{
   sprintf(blist->block->ac, "%s; distance from previous block=(%d,%d)",
           blist->block->number, blist->mindist, blist->maxdist);
} /* end of fix_ac */
/*=====================================================================
    Add a cluster pointing to the last sequence in the block
    Blimps assumes the sequences are in cluster order within the block
=======================================================================*/
void add_cluster(block)
Block *block;
{
   /*  Open up an extra cluster for the additional sequences */
   if (block->num_clusters + 1 > block->max_clusters)
   {
      resize_block_clusters(block);  /* updates block->max_clusters */
   }
   /*  First sequence will have been added to last cluster */
   block->clusters[ block->num_clusters - 1 ].num_sequences -= 1;
   block->clusters[ block->num_clusters ].num_sequences = 1;
   block->clusters[ block->num_clusters ].sequences = 
          &(block->sequences[ block->num_sequences - 1 ]);
   block->num_clusters += 1;
}  /* add_cluster */
