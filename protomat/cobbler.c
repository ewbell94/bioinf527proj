/* COPYRIGHT 1998 Fred Hutchinson Cancer Research Center                   */

/*    cobbler.c   Make consensus sequences from blocks 
	Always requires:
		BL, OU
	Embedded sequence from block requires:
		DB
	Embedded external sequence requires (FR & QM for alignment):
		SQ, FR, QM
        Consensus based on substitution scores requires:
		SU, FR, QM
	Consensus based on PSSM requires:
		FR, QM (not SU)
	Consensus based on percentage requires:
		PC (not SU or QM)  See NOTE above.
		Using 0% will select the aa with the max weight in each column.
		Using negative % will generate random consensus.

	cobbler <config file>
		Required: BL, OU, FR, QM and (PC or SU)
		Optional: (DB or SQ) and MA 
			  
		TY Consensus type:
			1 Percent of weighted counts ( percent_consensus() )
			  Requires PC
			2 Scores with weighted counts ( score_consensus() )
			  Requires SU
			3 Scores with weighted counts + pseudo-counts
				( score_pseudos() )
			  Requires SU, FR & QM
			4 Maximum PSSM residue ( pssm_consensus() )
			  Requires FR & QM; if MA is provided will also
			  print a profile.
		BL input file of blocks to make consensus from; must
		   contain sequence weights.
		OU output file of consensus blocks, or embedded
		   sequences if DB or SQ is specified.
		PC Percent identity for consensus residue in a position
		   of a block. PC = 0 => use most common residue; PC < 0
		   => use a residue at random. Or use SU instead.
		SU Name of substitution matrix to use to determine
		   consensus residues, instead of using PC. Takes residue
		   with highest pairwise score with other residues in position.
		MA Name of PSSM to use when making profile. Only used
		   when embedding sequence from block right now (DB).
		   This triggers profile output = .prf file in swat format.
		DB file containing sequences in blocks in BL
		   Consensus blocks will be embedded in the
		   first sequence most like the consensus & written
		   to file OU.
		SQ file containing query sequences NOT in blocks in BL.
		   First family of blocks in BL will be aligned with
		   sequence using a PSSM; then consensus blocks will
		   be embedded at the alignments if all blocks are
		   aligned in the correct order without overlap and
	           written to file OU. Requires FR and QM.
		FR Name of frequency matrix to use to make PSSM for SQ.
		QM Name of qij matrix to make PSSM for SQ.
		TR Trimmed length: length of cobbler sequence before &
		   after start of blocks to include in output sequence
		   (entire sequence is output by default)

	NOTE: Problem - can't embed external sequence with percentage
	      consensus because QM without SU now means compute
	      PSSM consensus ...

--------------------------------------------------------------------
10/16/93 J. Henikoff
====================================================================*/

#define COBBLER_C_
#define EXTERN
#define MAXNAME 80	/* Maximum file name length */
#define MAXLINE 240	/* Maximum length of input file line */
#define MAXAA 25	/* Dimension of aa_ arrays */
#define MAXWIDTH 80	/* Maximum block width */
#define CPROP .25	/* Default Minimum proportion of consensus weight */
#define PRF_SWAT 0	/* 0 = swat, 1 = profilesearch */
#define PRF_PFS 1
#define RAND_MAX 2147483647	/* 2^31 - 1 */
#define YES 1
#define NO 0
#define SHORT_TRIM 10	/* Default trimmed length */

#include "blocksprogs.h"
#include <sys/time.h>

/*   Blimps routines and data structures   */
extern int ErrorLevelReport;
extern struct float_qij *Qij;
extern double RTot;
extern double frequency[MATRIX_AA_WIDTH];
extern struct float_qij *load_qij();
extern int load_frequencies();
extern Matrix *block_to_matrix();

/*
 * Global variables and data structures
 */

/*  List of blocks structure:
    First entry has no block, just nblock, nseq, totwidth & minseq,
    other entries in list have pointers to the blocks, minprev & seqprev
*/
struct block_list {		/* list of blocks for a family */
   int nblock;				/* number of blocks in family */
   int nseq;				/* number of sequences in blocks */
   int totwidth;			/* total width of blocks in list  */
   int minseq;				/* sequence most like consensus */
   int maxseq;				/* sequence least like consensus */
   int minprev;				/* min distance from previous block*/
   int maxprev;
   int minseqprev;			/* previous distance for minseq */
   int maxseqprev;			/* previous distance for maxseq */
   int query_pos;
   int *consensus;			/* consensus for this block */
   Block *block;
   Matrix *matrix;
   struct block_list *next;		/* next block in ABC order */
   struct block_list *next_min;		/* next block in minseq order */
};
struct seqseq { int seq; };

struct sorttemp {
   int position;			/* sequence position */
   struct block_list *blist;		/* block list */
};

int read_cf();
void process_family();
void make_consensus();
void percent_consensus();
void score_consensus();
void score_pseudos();
void cons_diff();
void closest();
void print_consensus();
void embed_consensus();
void extract_seq();
void query_consensus();
int best_pos();
void pssm_consensus();

void make_random();
int ran0();

struct block_list *make_blist();
void insert_blist();
void free_blist();
void order_seq();
void load_sij();
int sortcmp();

int SeqType;			/*  why??? blimps stuff */
int Type;			/* Consensus type */
int PrfType;			/* Profile type */
int ShortTrim, Short;		/* Sequence trimming info */
double CProp;			/* Min weight % for consensus residue */
int Scores[MAXAA][MAXAA];	/*  Scoring matrix */
Matrix *Profile;		/* Baseline PSSM for making a profile */
FILE *Fbl, *Fou, *Fpr, *Fsq, *Fdb;	/* Config. files */
/*======================================================================
 * main
 */

int main(argc, argv)
     int argc;
     char *argv[];
{
  FILE *cfp;
  Block *block;
  Matrix *matrix;
  struct block_list *blist;
  char cfname[MAXNAME], save_number[10];
  struct timeval tv;

  ErrorLevelReport = 5;			/* don't want to see them */
  /* initialize random number generator */
  gettimeofday(&tv, NULL);
  srandom(tv.tv_sec^tv.tv_usec);

   if (argc <= 1)
   {
       printf("cobbler: Make a consensus sequence from a family of blocks\n");
       printf("USAGE:  cobbler <config_file>\n");
   }
/*------------1st arg = configuration file -------------------------------*/
   if (argc > 1)
      strcpy(cfname, argv[1]);
   else
   {
      printf("\nEnter name of configuration file: ");
      gets(cfname);
   }
   if ( (cfp=fopen(cfname, "r")) == NULL)
   {
      printf("\nCannot open file %s\n", cfname);
      exit(-1);
   }
   if ( (Type = read_cf(cfp)) < 0)
   {
      printf("Configuration file errors\n");
      exit(-1);
   }

  /*-------------------Read in each family of blocks -------------------*/
  save_number[0] = '\0';
  blist = NULL;
  while ((block = read_a_block(Fbl)) != NULL)
  {
     if (strncmp(save_number, block->number, 7) != 0)
     {  /*  new family */
         process_family(blist, Fou, Fpr); /* process previous family */
         /*   initialize new family */
         strncpy(save_number, block->number, 7); save_number[7] = '\0';
         if (Fdb == NULL && Fsq == NULL) 
              fprintf(Fou, ">%s %s\n", save_number, block->de); 
         blist = make_blist();			/* header for list */
         insert_blist(blist, block);		/* 1st block in list */
     }
     else      /*  same family */
     {
         insert_blist(blist,block);
     }
  }     /*  end of block */

  /*  process previous family */
  process_family(blist, Fou, Fpr);
   
  fclose(cfp); fclose(Fbl); fclose(Fou); fclose(Fpr);
  if (Fdb != NULL) fclose(Fdb);
  if (Fsq != NULL) fclose(Fsq);
  exit(0);

}  /* end of main */
/*=======================================================================*/
int read_cf(cfp)
FILE *cfp;
{
   char line[MAXLINE], keyword[20], *ptr;
   char blname[MAXNAME], ouname[MAXNAME], frname[MAXNAME];
   char suname[MAXNAME], qmname[MAXNAME], sqname[MAXNAME];
   char dbname[MAXNAME], maname[MAXNAME], ctemp[MAXNAME];
   int cscores, freqs;
   FILE *fqm, *fsu, *fma;

   CProp = CPROP;
   PrfType = PRF_SWAT;
   Type = -1;
   RTot = 5.0;
   Short = NO;
   cscores = freqs = NO;
   
   Fbl = Fou = Fsq = Fdb = fqm = fma = NULL;
   while (!feof(cfp) && fgets(line, MAXLINE, cfp) != NULL)
   {
      ptr = strtok(line, " \t\r\n");
      if (ptr != NULL && ptr[0] != ';')
      {
         strcpy(keyword,ptr);
         ptr = strtok(NULL, " \t\r\n;");
         if (ptr != NULL)
         {
            if (strncasecmp(keyword, "TY", 2) == 0)
            {
               Type = atoi(ptr);
            }
            else if (strncasecmp(keyword, "BL", 2) == 0)
            {
               strcpy(blname, ptr);
               if ( (Fbl = fopen(blname, "r")) == NULL)
                  printf("\nCannot open BL file %s\n", blname);
            }
            else if (strncasecmp(keyword, "OU", 2) == 0)
            {
               strcpy(ouname, ptr);
               if ( (Fou = fopen(ouname, "w")) == NULL)
                  printf("\nCannot open OU file %s\n", ouname);
               else printf("Consensus written to %s\n", ouname);
            }
            else if (strncasecmp(keyword, "SQ", 2) == 0)
            {
               strcpy(sqname, ptr);
               if ( (Fsq = fopen(sqname, "r")) == NULL)
                  printf("\nCannot open SQ file %s\n", sqname);
            }
            else if (strncasecmp(keyword, "DB", 2) == 0)
            {
               strcpy(dbname, ptr);
               if ( (Fdb = fopen(dbname, "r")) == NULL)
                  printf("\nCannot open DB file %s\n", dbname);
            }
            else if (strncasecmp(keyword, "FR", 2) == 0)
            {
               strcpy(frname, ptr);
               if (load_frequencies(frname) < 0)
                  printf("\nCannot open FR file %s\n", frname);
               else freqs = YES;
            }
            else if (strncasecmp(keyword, "PC", 2) == 0)
            {
               CProp = (double) atoi(ptr) / 100.;
               if (CProp > 1.0 || CProp < 0.0) CProp = -1.0;
            }
            else if (strncasecmp(keyword, "SU", 2) == 0)
            {
               strcpy(suname, ptr);
               if ( (fsu = fopen(suname, "r")) == NULL)
                  printf("\nCannot open SU file %s\n", suname);
               else
               {  cscores=YES;  load_sij(fsu, Scores);  fclose(fsu);  }
            }
            else if (strncasecmp(keyword, "MA", 2) == 0)
            {
               strcpy(maname, ptr);
               if ( (fma = fopen(maname, "r")) == NULL)
                  printf("\nCannot open MA file %s\n", maname);
               else
               {   Profile = read_a_matrix(fma);  fclose(fma);  }
            }
            else if (strncasecmp(keyword, "QM", 2) == 0)
            {
               strcpy(qmname, ptr);
               if ( (fqm = fopen(qmname, "r")) == NULL)
                  printf("\nCannot open QM file %s\n", qmname);
               else
               { Qij = load_qij(fqm); RTot = 5.0; fclose(fqm);   }
            }
            else if (strncasecmp(keyword, "TR", 2) == 0)
            {
                 Short = YES;
                 ShortTrim = atoi(ptr);
                 if (ShortTrim < 0) ShortTrim = SHORT_TRIM;
            }
         }   /* end of if ptr */
      }  /* end of if ptr */
   }  /* end of cfp */

   if (Profile != NULL)
   {
      strcpy(ctemp, ouname); strcat(ctemp, ".prf");
      if ( (Fpr = fopen(ctemp, "w")) == NULL)
      {
            printf("\nCannot open OU file %s\n", ctemp);
            free_matrix(Profile);
      }
      else printf("Profile written to %s\n", ctemp);
   }

   /*----------------------------------------------------------------*/
   /*    Figure out the type now   */
   /*  1  Percent
       2  Substitution matrix + counts
       3  Substitution matrix + counts + pseudo-counts
       4  PSSM
   */
   if (Type < 0)		/* TY wasn't specified explicitly */
   {
      if (cscores) Type = 3;
      else if (Qij != NULL && freqs) Type = 4;
      else Type = 1;
   }
   else				/* TY was specified explicitly */
   {
      if (Type == 1)
      {  cscores = NO;  }
      else if (Type == 2)
      {
         if (!cscores)
         {  printf("Missing SU input file name\n"); Type = -1; }
      }
      else if (Type == 3)
      {
         if (!cscores)
         {  printf("Missing SU input file name\n"); Type = -1; }
         if (Qij == NULL || !freqs)
         {  printf("FR and QM required with Type 3\n");  Type = -1; }
      }
      else if (Type == 4)
      {
         if (Qij == NULL || !freqs)
         {  printf("FR and QM required with Type 4\n");  Type = -1; }
      }
   }
   if (Type > 4)
   {  printf("Invalid consensus type: %d\n", Type); Type = -1; }
   if (Fbl == NULL)
   {  printf("Missing BL input file name\n"); Type = -1; }
   if (Fou == NULL)
   {  printf("Missing OU output file name\n"); Type = -1; }
   if (Fdb != NULL && Fsq != NULL)
   {  printf("Cannot specify both DB and SQ\n"); Type = -1; }
   if (Fsq != NULL && (Qij == NULL || !freqs))
   {  printf("FR and QM required with SQ\n");  Type = -1; }
   if (Fpr != NULL && (Qij == NULL || !freqs))
   {  printf("FR and QM required with MA\n");  Type = -1; }
   

   return(Type);
}  /* end of read_cf */
/*=======================================================================
========================================================================*/
void process_family(blist, ofp, pfp)
struct block_list *blist;
FILE *ofp, *pfp;
{

  if (blist != NULL && blist->nblock > 0)
  {
      switch (Type)
      {
         case 1:			/* percentage */
		printf("Type 1: Consensus based on percentage: %d\n",
			(int) (100. * CProp) );
		break;
	 case 2:			/* subst. with counts */
		printf("Type 2: Consensus based on substitution matrix");
                printf(" with counts (no pseudo-counts)\n");
		break;
	 case 3:			/*subst with counts + pseudos */
		printf("Type 3: Consensus based on substitution matrix");
		printf(" with counts + pseudo-counts\n");
		break;
	 case 4:			/* PSSM */
		printf("Type 4: Consensus based on PSSM\n");
		break;
	 default:
		break;
      }
     if (Type == 1 && CProp < 0.0)
         make_random(blist, ofp);
     else
     {
        make_consensus(blist); 
        if (Fdb != NULL)      embed_consensus(blist, ofp, pfp);
        else if (Fsq != NULL) query_consensus(blist, ofp, pfp);
             else             print_consensus(blist, ofp);
     }

     free_blist(blist);
   }

}   /* end of process_family  */
/*=======================================================================
      Turn a list of blocks into a consensus sequence
      Count the number of positions each sequence in the block differs
      from the consensus. The master sequence order is that in the 
      first block of the family.
========================================================================*/
void make_consensus(blist)
struct block_list *blist;
{
   int seq, pos, aa, width, dist;
   struct block_list *bcur;
   Block *block, *first_block;
   struct seqseq *sseq;

   sseq = (struct seqseq *) malloc (blist->nseq * sizeof(struct seqseq));
   first_block = NULL;
   bcur = blist->next;
   while(bcur != NULL && bcur->block != NULL)
   {
      /*  process each block separately  */
      block = bcur->block;

      /*-----  Figure out the consensus segment for this block ------------*/
      /*  (Puts result in bcur->consensus )   */
      switch (Type)
      {
         case 1:			/* percentage */
		percent_consensus(bcur);
		break;
	 case 2:			/* subst. with counts */
		score_consensus(bcur);
		break;
	 case 3:			/*subst with counts + pseudos */
		score_pseudos(bcur);
		break;
	 case 4:			/* PSSM */
		pssm_consensus(bcur);
		break;
	 default:
		break;
      }

      /*---------Count the number of differences b/w each seq & consensus--*/
      /*  (Puts result in first_block->sequences[seq].undefined)   */
      if (first_block == NULL) first_block = block;
      cons_diff(bcur, first_block, sseq);

      bcur = bcur->next;
   }  /* end of block */

   /*----------Determine the closest sequence ---------*/
   /*   (Puts results in block list records)           */
   /*  Extract blist->next->block->sequences[blist->minseq].name  */
   closest(blist, first_block, sseq);

   free(sseq);
}  /* end of make_consensus */
/*========================================================================
        Figure out the consensus for this block based on percentage
        & store it in the block list record
=========================================================================*/
void percent_consensus(bcur)
struct block_list *bcur;
{
   double sum[MAXAA+1], tot, maxw;
   int seq, pos, aa, width, maxaa, dupmax;
   Block *block;

   block = bcur->block;
   width = block->sequence_length;

   /*-------------------------------------------------------------------*/
   for (pos = 0; pos < width; pos++)
   { 
      tot = 0.0;
      maxw = maxaa = dupmax = -1;
      for (aa = 0; aa < MAXAA+1; aa++) sum[aa] = 0.0;
  
      for (seq = 0; seq < block->num_sequences; seq++)
      {
         if (block->residues[seq][pos] >= 0 &&
                block->residues[seq][pos] < MAXAA)
         {
            sum[block->residues[seq][pos]] += block->sequences[seq].weight;
            tot += block->sequences[seq].weight;
         }
         else printf("Unknown residue: %d\n", block->residues[seq][pos]);
      }

      for (aa = 0; aa < MAXAA+1; aa++)
      {
         if (sum[aa] > maxw)
         {
            maxw = sum[aa]; maxaa = aa;
         }
         else if (sum[aa] > 0.0 && sum[aa] == maxw)
                  dupmax = aa;
      }
      /*-------  Print "x" if maximum weight for any aa is < CProp% ------*/
      if ((maxw / tot) < CProp) maxaa = MAXAA;

      /*--------Save the consensus sequence------------------------*/
      if (maxaa >= 0 && maxaa < MAXAA )
           bcur->consensus[pos] = maxaa;
      else
           bcur->consensus[pos] = MAXAA;	/* X==23 */
   }  /*  end of pos */

}   /*  end of percent_consensus  */
/*========================================================================
        Figure out the consensus for this block & store it in the
        block list record - use highest total pairwise score
        Prints x if the best total score is negative
=========================================================================*/
void score_consensus(bcur)
struct block_list *bcur;
{
   double sum[MAXAA+1], maxw, weight;
   int seq1, seq2, pos, width, nseq, maxaa, aamark[MAXAA];
   int aa1, aa2;
   Block *block;

   block = bcur->block;
   width = block->sequence_length;
   nseq = block->num_sequences;

   /*-------------------------------------------------------------------*/
   for (pos = 0; pos < width; pos++)
   { 
      /*    good sums can be negative */
      maxw = -999999.99;
      maxaa = -1;
      for (aa1 = 0; aa1 < MAXAA+1; aa1++)
      {
         sum[aa1] = 0.0;
         aamark[aa1] = NO;	/* does this aa appear in the pos ? */
      }
  
      /*   Mark the aas that appear in this position */
      for (seq1 = 0; seq1 < nseq; seq1++)
      {
            aa1 = block->residues[seq1][pos];
            aamark[aa1] = YES;
      }
/*		Pair counts code
      for (seq1 = 0; seq1 < nseq; seq1++)
         for (seq2 = seq1 + 1; seq2 < nseq; seq2++)
         {
            aa1 = block->residues[seq1][pos];
            aa2 = block->residues[seq2][pos];
            weight =  block->sequences[seq1].weight;
            weight *= block->sequences[seq2].weight;
            sum[aa1] += weight * Scores[aa1][aa2];
            if (aa2 != aa1) sum[aa2] += weight * Scores[aa1][aa2];
            aamark[aa1] = aamark[aa2] = YES;
         }
*/

      /*    Only choose among the aas that appear in the pos  */
      for (aa1 = 0; aa1 < MAXAA; aa1++)
         if ( aamark[aa1] )
         {
            sum[aa1] = 0.0;
            for (seq1 = 0; seq1 < nseq; seq1++)
            {
               aa2 = block->residues[seq1][pos];
               weight =  block->sequences[seq1].weight;
               sum[aa1] += weight * Scores[aa1][aa2];
            }
         }

      /*    Only choose among the aas that appear in the pos  */
      for (aa1 = 0; aa1 < MAXAA; aa1++)
      {
         if (aamark[aa1] && sum[aa1] > maxw)
         {
            maxw = sum[aa1]; maxaa = aa1;
         }
      }
      if (maxw <= 0.0) maxaa = MAXAA;	/*  Use X if non-positive sum */

      /*--------Save the consensus sequence------------------------*/
      if (maxaa > 0 && maxaa < MAXAA )
           bcur->consensus[pos] = maxaa;
      else
           bcur->consensus[pos] = MAXAA;	/* X==23 */
   }  /* end of pos */

}  /*  end of score_consensus */

/*========================================================================
        Figure out the consensus for this block & store it in the
        block list record - use highest total pairwise score
        Prints x if the best total score is negative
	Variation using pseudo-counts
=========================================================================*/
void score_pseudos(bcur)
struct block_list *bcur;
{
   double sum[MAXAA+1], maxw, weight;
   int seq1, seq2, pos, width, nseq, maxaa, aamark[MAXAA];
   int aa1, aa2;
   Block *block;
   Matrix *matrix;

   block = bcur->block;
   width = block->sequence_length;
   nseq = block->num_sequences;
   matrix = block_to_matrix(block, 20);		/* counts+pseudos */

   /*-------------------------------------------------------------------*/
   for (pos = 0; pos < width; pos++)
   { 
      /*    good sums can be negative */
      maxw = -999999.99;
      maxaa = -1;
      for (aa1 = 0; aa1 < MAXAA+1; aa1++)
      {
         sum[aa1] = 0.0;
         aamark[aa1] = NO;	/* does this aa appear in the pos ? */
      }
  
      /*   Mark the aas that appear in this position */
      for (seq1 = 0; seq1 < nseq; seq1++)
      {
            aa1 = block->residues[seq1][pos];
            aamark[aa1] = YES;
      }

      /*    Only score the aas that appear in the pos  */
      for (aa1 = 0; aa1 < MAXAA; aa1++)
      {
         sum[aa1] = 0.0;
         if ( aamark[aa1] )
         {
            /*  Just the major 20 aas */
            for (aa2=1; aa2 <= 20; aa2++)
               sum[aa1] += matrix->weights[aa2][pos] * Scores[aa1][aa2];
         }
      }

      /*    Only choose among the aas that appear in the pos  */
      for (aa1 = 0; aa1 < MAXAA; aa1++)
      {
         if (aamark[aa1] && sum[aa1] > maxw)
         {
            maxw = sum[aa1]; maxaa = aa1;
         }
      }
      if (maxw <= 0.0) maxaa = MAXAA;	/*  Use X if non-positive sum */

      /*--------Save the consensus sequence------------------------*/
      if (maxaa > 0 && maxaa < MAXAA )
           bcur->consensus[pos] = maxaa;
      else
           bcur->consensus[pos] = MAXAA;	/* X==23 */
   }  /* end of pos */

   free_matrix(matrix);
}  /*  end of score_pseudos */
/*=======================================================================
      Count the number of differences b/w each seq & consensus
      (Puts result in first_block->sequences[seq].undefined) 
========================================================================*/
void cons_diff(bcur, first_block, sseq)
struct block_list *bcur;
Block *first_block;
struct seqseq *sseq;
{
   int width, seq, pos;
   Block *block;

      block = bcur->block;
      width = block->sequence_length;
      /*----- Count the number of differences in this block ------*/
      for (seq = 0; seq < block->num_sequences; seq++)
         block->sequences[seq].undefined = 0;
      for (seq = 0; seq < block->num_sequences; seq++)
         for (pos = 0; pos < width; pos++)
            if (bcur->consensus[pos] >= 0 && bcur->consensus[pos] < MAXAA &&
                (int) block->residues[seq][pos] != bcur->consensus[pos])
                   block->sequences[seq].undefined += 1;
      /*------Total diffs for all blocks accumulated in first block -----*/
      if (first_block != NULL && first_block != block)
      {
         /*--- put block sequences in same order as first_block -----*/
         order_seq(sseq, first_block, block);
         for (seq=0; seq < block->num_sequences; seq++)
             first_block->sequences[seq].undefined += 
                   block->sequences[ sseq[seq].seq ].undefined;
      }

}  /*   end of cons_diff */
/*=====================================================================
   Determine which sequence is most like the consensus
   and update seqprev fields with its distances between blocks
=======================================================================*/
void closest(blist, first_block, sseq)
struct block_list *blist;
Block *first_block;
struct seqseq *sseq;
{
   Block *block, *save_block;
   struct block_list *bcur;
   int i, mindiff, maxdiff, width, seq, save_seq;
   struct sorttemp stemp[26];

   /*--------Determine the closest sequence -------------------------*/
   /*  (cumulative differences are stored in first_block)            */
   mindiff = 999; maxdiff = save_seq = 0;
   for (seq=0; seq < blist->nseq; seq++)
      if (first_block->sequences[seq].undefined < mindiff)
      {
         mindiff = first_block->sequences[seq].undefined;
         blist->minseq = seq;
      }
      else if (first_block->sequences[seq].undefined > maxdiff)
      {
         maxdiff = first_block->sequences[seq].undefined;
         blist->maxseq = seq;
      }
   printf("%s: First sequence most  like consensus is %d %s (%d)\n",
            first_block->number, blist->minseq,
            first_block->sequences[blist->minseq].name, mindiff);
   printf("%s: First sequence least like consensus is %d %s (%d)\n",
            first_block->number, blist->maxseq,
            first_block->sequences[blist->maxseq].name, maxdiff);

   /*  The order of the blocks in minseq may be different than the
       ABC order, if repeats were allowed by PROTOMAT; get the
       correct order in minseq */
   i = 0;
   bcur = blist->next;
   while(bcur != NULL && bcur->block != NULL)
   {
      block = bcur->block;
      if (block == first_block)
      {   seq = blist->minseq;  }
      else
      {
         order_seq(sseq, first_block, block);
         seq = sseq[blist->minseq].seq;
      }
      stemp[i].blist = bcur;
      stemp[i].position = bcur->block->sequences[seq].position;
      i++;
      bcur = bcur->next;
   }
   qsort(stemp, blist->nblock, sizeof(struct sorttemp), sortcmp);
  
   bcur = blist;
   for (i=0; i< blist->nblock; i++)
   {
      bcur->next_min = stemp[i].blist;
      bcur = bcur->next_min;
   }
   bcur->next_min = NULL;

   /*-------------------------------------------------------------------*/
   /*------ Determine the distance between blocks in the closest seq ---*/
   save_seq = 0;
   bcur = blist->next_min;
   while(bcur != NULL && bcur->block != NULL)
   {
      block = bcur->block;
      width = block->sequence_length;
/*      if (block == first_block)  */
      if (bcur == blist->next_min)
      {
         bcur->minseqprev = block->sequences[blist->minseq].position - 1;
         save_seq = blist->minseq;
      }
      else
      {
         order_seq(sseq, first_block, block);
         seq = sseq[blist->minseq].seq;
         /*  from end to last block to start of this one ... */
         bcur->minseqprev = block->sequences[seq].position - 
    (save_block->sequences[save_seq].position + save_block->sequence_length);
         save_seq = seq;
      }
/*printf("%d\n", bcur->minseqprev); */
/*    bcur->minseqprev = bcur->minprev; */

      save_block = block;
      bcur = bcur->next_min;
   }  /* end of block */

}  /* end of closest */

/*=======================================================================
  Imbed the consensus sequence for a list of blocks
    Extract blist->next->block->sequences[blist->minseq].name 
    If Profile exists, output a profile as well in PrfType format
========================================================================*/
void embed_consensus(blist, ofp, pfp)
struct block_list *blist;
FILE *ofp, *pfp;
{
   int sngle, aa, pos, savpos, totpos, width, aacount[25];
   int firstpos, firstaa, lastaa;
   struct block_list *bcur;
   Block *block, *first_block;
   char *sequence, ctemp[12];

   /*   Figure out about how long the sequence is */
   totpos = 0;
   bcur = blist->next;
   while(bcur != NULL && bcur->block != NULL)
   {
      block = bcur->block;
      width = block->sequence_length;
      totpos += width + bcur->minseqprev;
      bcur = bcur->next;
   }

   /* Get the sequence, which is at least totpos long */
   /* Problem with this if there is a lot of sequence after the last block */
   totpos *= 2;
   if (totpos < 1000) totpos = 1000;
   sequence = (char *) malloc( totpos * sizeof(char) );
   extract_seq(Fdb, blist->next->block->sequences[blist->minseq].name,
               sequence);

   for (aa=0; aa < 25; aa++) aacount[aa] = 0;

   /*----------------------------------------------------------------*/
   /*     Now print out the sequence in lower case with the consensus
          blocks embedded in upper case    */

   /*  Set the first and last aas to print */
   firstaa = 0; lastaa = (int) strlen(sequence);
   if (Short)
   {
      bcur = blist->next_min;
      /*  First block */
      if (bcur->block != NULL)
      {
         first_block = bcur->block;
         firstaa = bcur->minseqprev - ShortTrim;
         if (firstaa < 0) firstaa = 0;
      }
      lastaa = 0;
      while(bcur != NULL && bcur->block != NULL)
      {
         lastaa += bcur->minseqprev;
         lastaa += bcur->block->sequence_length;
         bcur = bcur->next_min;
      }
      lastaa += ShortTrim;
   }
   if (firstaa < 0) firstaa = 0;
   if (lastaa > (int) strlen(sequence)) lastaa = (int) strlen(sequence);
   if (lastaa < firstaa)
   {   firstaa = 0; lastaa = (int) strlen(sequence);   }
   
   /*     Headings for both the COBBLER sequence & PROFILE  */
   totpos = firstaa;
   ctemp[1] = '\0';
   strcpy(ctemp, blist->next->block->number);
   if (ctemp[0] == 'P' && ctemp[1] == 'S')
   {  ctemp[0] = 'B'; ctemp[1] = 'L'; }
   ctemp[7] = '\0';		/* family name */
   fprintf(ofp, ">%s %s from %d to %d with embedded consensus blocks\n",
           ctemp, blist->next->block->sequences[blist->minseq].name,
           firstaa+1, lastaa);
   if (pfp != NULL && strlen(sequence) )
   {
      if (PrfType == PRF_PFS)
      {
         fprintf(pfp, "(Peptide) Length: %d\n", (int) strlen(sequence) );
         fprintf(pfp, "COBBLER Profile for %s %s from %d to %d\n",
           ctemp, blist->next->block->sequences[blist->minseq].name,
           firstaa+1, lastaa);
         fprintf(pfp, "Cons  ");
         for (aa=1; aa < 23; aa++) fprintf(pfp, "  %c  ", aa_btoa[aa]);
         fprintf(pfp,"Gap Len ..\n");
       }
       else
       {
          fprintf(pfp, "#COBBLER Profile for %s %s from %d to %d\n#\n      ",
           ctemp, blist->next->block->sequences[blist->minseq].name,
           firstaa+1, lastaa);
          for (aa=1; aa < 21; aa++) fprintf(pfp, "  %c  ", aa_btoa[aa]);
          fprintf(pfp, "\n");
       }
   }

   /* -------------------Now the body-------------------------------- */
   /*  This ASSUMES the blocks are in order in the chosen sequence,
	but they may not be if PROTOMAT allowed repeats, etc
	Need to fix this - eg BL50006 in Blocks 9.0   */
   bcur = blist->next_min;
   while(bcur != NULL && bcur->block != NULL)
   {
      block = bcur->block;
      width = block->sequence_length;
      if (block == first_block) firstpos = firstaa;
      else firstpos = 0;
      /*   Before the block  */
      for (pos = firstpos; pos < bcur->minseqprev; pos++)
      {	
          if (lastaa > totpos)
          {
             ctemp[0] = sequence[totpos];
             fprintf(ofp, "%c", tolower(ctemp[0]) );
          }
          else
          {
             ctemp[0] = 'X';
             fprintf(ofp, "X");
          }
	  sngle = aa_atob[ ctemp[0] ];
          aacount[sngle]++;
          if (pfp != NULL && lastaa )
          {
             fprintf(pfp, "%c    ", tolower(ctemp[0]) );
             if (PrfType == PRF_PFS)
             {
                for (aa=1; aa < 23; aa++)
                   fprintf(pfp, "%4d ",
                     (int) round(10.0 * Profile->weights[aa][sngle]));
                fprintf(pfp, " 100 100\n");
             }
             else
             {
                for (aa=1; aa < 21; aa++)
                   fprintf(pfp, "%4d ",
                     (int) round(Profile->weights[aa][sngle]));
                fprintf(pfp, "\n");
             }
          }

          totpos++;
          if ( totpos%MAXWIDTH == 0) fprintf(ofp, "\n");
      }
      /*   Within the block region */
      for (pos = 0; pos < width; pos++)
      {
         if (bcur->consensus[pos] >= 0 && bcur->consensus[pos] < MAXAA )
         {
              ctemp[0] = aa_btoa[ bcur->consensus[pos] ];
              fprintf(ofp, "%c", ctemp[0]);
              aacount[ bcur->consensus[pos] ]++;
         }
         else
         {		/* no consensus residue selected */
/*
             if (lastaa >= totpos)
             {
                ctemp[0] = sequence[totpos];
                fprintf(ofp, "%c", tolower(ctemp[0]) );
             }
             else fprintf(ofp, "x");
*/
             ctemp[0] = 'x';
             fprintf(ofp, "x");
             aacount[ aa_atob['X'] ]++;
         }

         if (pfp != NULL && lastaa )
         {
/*			Uses the 1/3 bit PSSM made by insert_blist()    */
            fprintf(pfp, "%c    ", ctemp[0]);
            if (PrfType == PRF_PFS)
            {
               for (aa=1; aa < 23; aa++)
                fprintf(pfp, "%4d ",
                  (int) round(10.0 * bcur->matrix->weights[aa][pos]));
               fprintf(pfp, " 100 100\n");
            }
            else		/* swat format */
            {
               for (aa=1; aa < 21; aa++)
                fprintf(pfp, "%4d ",
                  (int) round(bcur->matrix->weights[aa][pos]));
               fprintf(pfp, "\n");
            }
         }

         totpos++;
         if ( totpos%MAXWIDTH == 0) fprintf(ofp, "\n");
      }
      bcur = bcur->next_min;
   }  /* end of a block */

   /* ---Print the rest of the sequence now after the blocks------------- */
   if (lastaa > totpos)
   {
      savpos = totpos;
      for (pos = savpos; pos < lastaa; pos++)
      {
          ctemp[0] = sequence[pos]; 
          fprintf(ofp, "%c", tolower(ctemp[0]) );
	  sngle = aa_atob[ ctemp[0] ];
          aacount[sngle]++;
          if (pfp != NULL && lastaa )
          {
             fprintf(pfp, "%c    ", tolower(ctemp[0]));
             if (PrfType == PRF_PFS)
             {
                for (aa=1; aa < 23; aa++)
                   fprintf(pfp, "%4d ",
                      (int) round(10.0 * Profile->weights[aa][sngle]));
                fprintf(pfp, " 100 100\n");
             }
             else
             {
                for (aa=1; aa < 21; aa++)
                   fprintf(pfp, "%4d ",
                      (int) round(Profile->weights[aa][sngle]));
                fprintf(pfp, "\n");
             }
          }

          totpos++;
          if ( totpos%MAXWIDTH == 0) fprintf(ofp, "\n");
      }
      fprintf(ofp, "\n");
   }  /* end of if sequence was found */
   else fprintf(ofp, "\n");

   if (pfp != NULL && lastaa && PrfType == PRF_PFS)
   {
      fprintf(pfp, "*    ");
      for (aa=1; aa<23; aa++)
         fprintf(pfp, "%4d ", aacount[aa]);
      fprintf(pfp, "\n");
   }

   free(sequence);
}   /* end of embed_consensus  */
/*=======================================================================
	Extract the named sequence from a fasta format file of sequences
	aa_atob[] values:  0- 1A 2R ...20V 21B 22Z 23X,J,O,U 24* 25%,$,...
                              27space
        If the name is null, extract the first sequence
========================================================================*/
void extract_seq(pfp, name, sequence)
FILE *pfp;
char *name, *sequence;
{
   int spos, found, nseq, itemp, i;
   char line[MAXLINE], *ptr;
 
   found = NO;
   nseq = spos = 0;
   sequence[0] = '\0';
   rewind(pfp);
   while (!found & (fgets(line, MAXLINE, pfp)) != NULL)
   {
      if (line[0] == '>')
      {
         nseq++;
/*
         if (!strlen(name) || (strncmp(line+1, name, strlen(name)) == 0) )
*/
         /*  Look for sequence name anywhere on title line */
         if (!strlen(name) || (strstr(line, name) != NULL) )
         {
            found = YES;
            while ( fgets(line, MAXLINE, pfp) != NULL &&
                    line[0] != '>')
            {
	 	ptr = strtok(line, "\r\n");
                /*  Have to get rid of spaces, etc */
                for (i=0; i < (int) strlen(line); i++)
                {
                    itemp = aa_atob[ line[i] ];
                    if (itemp >= 1 && itemp <= 23)
                        sequence[spos++] = line[i];
                }
            }
         }
      }
   }
   if (!found) printf("Sequence %s not found\n", name);
   sequence[spos] = '\0';
}  /* end of extract_seq  */

/*=======================================================================
  Print out the consensus sequence for a list of blocks
========================================================================*/
void print_consensus(blist, ofp)
struct block_list *blist;
FILE *ofp;
{
   int pos, itemp, width;
   struct block_list *bcur;
   Block *block;

   bcur = blist->next_min;
   while(bcur != NULL && bcur->block != NULL)
   {
      block = bcur->block;
      width = block->sequence_length;
      /*------- Print leading X's for minimum preceding distance ------*/
      for (pos = 0; pos < bcur->minseqprev; pos++)
      {
          fprintf(ofp, "X");
          if ( (pos+1)%MAXWIDTH == 0) fprintf(ofp, "\n");
      }
      itemp = (int) bcur->minseqprev/MAXWIDTH;
      itemp = bcur->minseqprev - itemp * MAXWIDTH;
      if ((itemp + width) > MAXWIDTH) fprintf(ofp, "\n");
  
      for (pos = 0; pos < width; pos++)
      {
         if (bcur->consensus[pos] >= 0 && bcur->consensus[pos] < MAXAA )
              fprintf(ofp, "%c", aa_btoa[ bcur->consensus[pos] ]);
         else
              fprintf(ofp, "x");
      }
      fprintf(ofp, "\n");
      bcur = bcur->next_min;
   }  /* end of block */
   fprintf(ofp, "\n");
}   /*  end of print_consensus  */

/*=======================================================================
      Turn a block into a random consensus sequence
========================================================================*/
void make_random(blist, ofp)
struct block_list *blist;
FILE *ofp;
{
   double rrand;
   int seq, pos, aa, width, nseq, randaa[MAXWIDTH], itemp;
   char ctemp[MAXNAME], *ptr;
   struct block_list *bcur;
   Block *block;

   bcur = blist->next;
/*
   if (bcur->block != NULL)
       fprintf(ofp, ">%s %s\n", bcur->block->number, bcur->block->de);
*/
   while (bcur !=NULL && bcur->block != NULL)
   {
      block = bcur->block;
      width = block->sequence_length;
      nseq = block->num_sequences;
      for (pos = 0; pos < width; pos++)
      { 
         rrand = (double) ran0() / (RAND_MAX + 1.0);  /* bw 0.0 and 1.0 */
         seq = (int) (nseq - 1) * rrand;
         randaa[pos] = block->residues[seq][pos];
      }

      /*------- Print leading X's for minimum preceding distance ------*/
      for (pos = 0; pos < bcur->minprev; pos++)
      {
          fprintf(ofp, "X");
          if ( (pos+1)%MAXWIDTH == 0) fprintf(ofp, "\n");
      }
      itemp = (int) bcur->minprev/MAXWIDTH;
      itemp = bcur->minprev - itemp * MAXWIDTH;
      if ((itemp + width) > MAXWIDTH) fprintf(ofp, "\n");
  
      /*--------Print the consensus block now ---------------------------*/
      for (pos = 0; pos < width; pos++)
      {
         if (randaa[pos] >= 0 && randaa[pos] < MAXAA )
            fprintf(ofp, "%c", aa_btoa[randaa[pos]]);
         else fprintf(ofp, "x");
      }
      fprintf(ofp, "\n");
      bcur = bcur->next;
   }  /*  end of block */
}  /* end of make_random */
/*=================================================================
   From Press, et al, "Numerical Recipes in C", pp.207-208
   Improved random number generator (breaks up sequential correlations)

   Assumes random number generator has been seeded elsewhere (does
   not call srand).
   Assumes caller will divide results by RAND_MAX+1 (so it can be
   used the same as rand()).
========================================================================*/
int ran0()
{
   static int y, v[98];
   static int iff=0;			/* first call flag */
   int j;

   if (iff == 0)	/* fill array on first call */
   {
      iff=1;
      for (j=0; j<=97; j++) v[j] = random();
      y=random();
   }
   j = 1 + 97.0 * y / (RAND_MAX+1.0);
   if (j >=0 && j <= 97)
   {
      y=v[j];
      v[j] = random();
   }
   else y=random();
   return (y);
}  /* end of ran0 */
/*=======================================================================
     routines for a list of blocks
========================================================================*/
struct block_list *make_blist()
{
   struct block_list *new;
   
   new = (struct block_list *) malloc (sizeof(struct block_list));
   new->nblock = new->nseq = new->totwidth = new->minprev = new->maxprev = 0;
   new->minseq = new->minseqprev = new->query_pos = 0;
   new->maxseq = new->maxseqprev = 0;
   new->consensus = NULL;
   new->block = NULL;
   new->matrix = NULL;
   new->next = new->next_min = NULL;

   return(new);
}  /* end of make_blist */

/*---------------------------------------------------------------------
       Inserts the block in a list of blocks
       Computes the 1/3 bit PSSM for the block
-----------------------------------------------------------------------*/
void insert_blist(blist, block)
struct block_list *blist;
Block *block;
{
   struct block_list *cur;
   char *ptr, ctemp[MAXNAME];

   /*------ Accumulate totals in header record ------*/
   blist->nblock += 1;
   blist->nseq = block->num_sequences;
   blist->totwidth += block->sequence_length;

   /*--- Insert a new record for the current block at the end of the list ---*/
   cur = blist;
   while (cur->next != NULL)
      cur = cur->next;

   cur->next = make_blist();
   cur->next_min = cur->next;
   cur->next->consensus = (int *) malloc(block->sequence_length * sizeof(int));
   cur->next->block = block;
   /*  Following requires global vars Qij and frequency[]  */
   /*   NOTE:  For Gribskov-type PSSM change 6 to 30 in block_to_matrix()
               and Qij should be a substitution matrix    */
   if (Qij != NULL)
      cur->next->matrix = block_to_matrix(block, 6);	/* third bit PSSM */
   /*------- get minimum preceding distance ------*/
   cur->next->minprev = cur->next->maxprev = 0;
   strcpy(ctemp, block->ac);
   ptr = strtok(ctemp, "(");
   if (ptr != NULL)
   {
      ptr = strtok(NULL, ",");
      if (ptr != NULL) cur->next->minprev = atoi(ptr);
   }
   if (cur->next->minprev < 0) cur->next->minprev = 0;
   
}  /* end of insert_blist */

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
      if (cur->consensus != NULL) free(cur->consensus);
      if (cur->block!= NULL) free_block(cur->block);
      if (cur->matrix != NULL) free_matrix(cur->matrix);
      free(cur);
      last->next = last->next_min = NULL;
      free_blist(last);
   }
   else free(blist);

}  /* end of free_blist */
/*=======================================================================*/
/*  Sequences may not be in the same order both blocks. If not, then
    set sseq[s1] = s2  where
    b1->sequences[s1].name == b2->sequences[s2].name
=========================================================================*/
void order_seq(sseq, b1, b2)
struct seqseq *sseq;
Block *b1, *b2;
{
   int nseq, i1, i2;

   nseq = b1->num_sequences;
   if (b2->num_sequences < nseq) nseq = b2->num_sequences;
   for (i1 = 0; i1 < nseq; i1++)
   {
      if (b1 == b2) sseq[i1].seq = i1;
      else
      {
         sseq[i1].seq = -1;
         i2 = 0;
         while (sseq[i1].seq < 0 && i2 < nseq)
         {
            if (strcmp(b1->sequences[i1].name, b2->sequences[i2].name) == 0)
               sseq[i1].seq = i2;
            i2++;
         }
      }
   }
}  /*  end of order_seq */
/*======================================================================
	Read an integer-valued substitution matrix
	MAXAA == 25 - should use 24 here?
========================================================================*/
void load_sij(fin, scores)
FILE *fin;
int scores[MAXAA][MAXAA];
{
   char line[132], *ptr;
   int alpha[MAXAA], nrows, ncols, row, col, i;

/*----------Read file until first non-blank line --------------*/
/* Skip comments at beginning of file - 1st char = #, > or ;   */
   line[0] = '\0';
   while (((int) strlen(line) < 1 || line[0]=='#' || line[0]=='>' || 
                                     line[0]==';')
          && fgets(line, sizeof(line), fin) != NULL)
	    ;
/*------See if the first line has characters on it ------------*/
   for (col=0; col < MAXAA; col++) alpha[col] = -1;
   if (strstr(line, "A") != NULL)	/* This line has characters */
   {
      row = 0;	/* # of alphabetic characters on the line */
      for (i=0; i< (int) strlen(line); i++)
      {
	 col = aa_atob[ line[i] ];
	 if (col >= 0 && col < MAXAA)
	 {
	    alpha[row] = col;
	    row++;
	 }
	 else if (isalpha(line[i])) row++;
      }
   }
/*-------Get the data values now ------------*/
   for (row=0; row<MAXAA; row++)
     for (col=0; col<MAXAA; col++)
	scores[row][col] = -999;		/* Null value */
   nrows = 0;
   line[0] = '\0';
   while (fgets(line, sizeof(line), fin) != NULL)
   {
      if ((int) strlen(line) > 1 && nrows < MAXAA)
      {
	 if (alpha[nrows] >= 0 && alpha[nrows] < MAXAA)
	 {
	    row = alpha[nrows]; ncols = 0;
	    ptr = strtok(line, " ,\n");
	    while (ptr != NULL)
	    {
	       if (strspn(ptr, "+-0123456789") == strlen(ptr))
	       {
		  col = alpha[ncols];
		  if (col >= 0 && col < MAXAA)
		     scores[row][col] = atoi(ptr);
		  ncols++;
	       }
	       ptr = strtok(NULL, " ,\n");
	    }
	 }
	 nrows++;
      }
   }

/*-------If some entries are still missing, assume symmetry ---------*/
   for (row=0; row<MAXAA; row++)
   {
     for (col=0; col<MAXAA; col++)
     {
	if (scores[row][col] == -999) scores[row][col] = scores[col][row];
	if (row==20 && scores[row][col]==-999)	/*  B -> D */
	   scores[row][col] = scores[3][col];
	if (row==21 && scores[row][col]==-999)	/* Z -> E */
	   scores[row][col] = scores[6][col];
	if (col==20 && scores[row][col]==-999)	/*  B -> D */
	   scores[row][col] = scores[row][3];
	if (col==21 && scores[row][col]==-999)	/* Z -> E */
	   scores[row][col] = scores[row][6];
     }
   }

}  /* end of load_sij */
/*=======================================================================
   Read in query sequence, align it with blocks, embed consensus
   at alignment points if all blocks align in the right order without
   overlapping.
========================================================================*/
void query_consensus(blist, ofp, pfp)
struct block_list *blist;
FILE *ofp, *pfp;
{
   int pos, lastpos, totpos, savpos, width;
   char name[15], ctemp[2];
   Sequence *seq;
   Block *block;
   struct block_list *bcur;

   seq = read_a_sequence(Fsq, UNI, AA_SEQ);
   /*   If query is shorter than the sum of the block widths
        just print a message  */
   bcur = blist->next;
   width = 0;
   while (bcur != NULL && bcur->block != NULL)
   {
      width += bcur->minprev;
      width += bcur->block->sequence_length;
      bcur = bcur->next;
   }
   if (seq->length < width)
   {
      printf("Query sequence %s is shorter than minimum path width\n",
              seq->name);
      free_sequence(seq);
      seq = NULL;
   }
 
   /*   Align the sequence with each block & get the value of
        minseqprev  */
   lastpos = 0;
   if (seq != NULL)
   {
      bcur = blist->next;
      while (bcur != NULL && bcur->block != NULL)
      {
         width = bcur->block->sequence_length;
/*       matrix = block_to_matrix(bcur->block, 3);   scaled 0-99  */
/*		Uses 1/3 bit PSSM made by insert_blist()          */
         bcur->query_pos = best_pos(bcur->matrix, seq); 
         if (bcur->query_pos < lastpos - 1)
         {
            printf("WARNING: Aligned blocks overlap!!! ");
            printf(" block=%s, pos=%d, last pos = %d\n",
             bcur->block->number, bcur->query_pos, lastpos);
         }
	 else
	 {
            printf(" block=%s, pos=%d\n",
             bcur->block->number, bcur->query_pos);
         }
         bcur->minseqprev = bcur->query_pos - lastpos;
         lastpos = bcur->query_pos + width;
         bcur = bcur->next;
      }
      /*     Now print out the sequence in lower case with the consensus
          blocks embedded in upper case    */
      totpos = 0;
      ctemp[1] = '\0';
      fprintf(ofp, ">%s with embedded consensus blocks\n", seq->name);
      bcur = blist->next;
      while(bcur != NULL && bcur->block != NULL)
      {
         block = bcur->block;
         width = block->sequence_length;
         for (pos = 0; pos < bcur->minseqprev; pos++)
         {

             if (seq->length >= totpos)
             {
                ctemp[0] = aa_btoa[ seq->sequence[totpos] ];
                fprintf(ofp, "%c", tolower(ctemp[0]) );
             }
             else
                fprintf(ofp, "X");
             totpos++;
             if ( totpos%MAXWIDTH == 0) fprintf(ofp, "\n");
         }
         for (pos = 0; pos < width; pos++)
         {
            if (bcur->consensus[pos] >= 0 && bcur->consensus[pos] < MAXAA )
              fprintf(ofp, "%c", aa_btoa[ bcur->consensus[pos] ]);
            else
            {
/*
             if (seq->length >= totpos)
             {
                ctemp[0] = aa_btoa[ seq->sequence[totpos] ];
                fprintf(ofp, "%c", tolower(ctemp[0]) );
             }
             else fprintf(ofp, "x");
*/
                fprintf(ofp, "x");
            }
            totpos++;
            if ( totpos%MAXWIDTH == 0) fprintf(ofp, "\n");
         }
         bcur = bcur->next;
      }
      /*  print the rest of the sequence now */
      if (seq->length >= totpos)
      {
         savpos = totpos;
         for (pos = savpos; pos < seq->length; pos++)
         {
             ctemp[0] = aa_btoa[ seq->sequence[pos] ];
             fprintf(ofp, "%c", tolower(ctemp[0]) );
             totpos++;
             if ( totpos%MAXWIDTH == 0) fprintf(ofp, "\n");
         }
         fprintf(ofp, "\n");
      }
      free_sequence(seq);
   } /*  end of if sequence was found */
  
}   /*  end of query_consensus */
/*=======================================================================
            Score all alignments of matrix with sequence & return
	     position of best alignment
============================================================================*/
int best_pos(matrix, seq)
Matrix *matrix;
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
  for (seq_pos= -matrix->length+1; seq_pos < seq->length; seq_pos++) {


    /* 
     * score this alignment
     */
    seq_score = 0.0;

    /* for each position in the matrix add the weight to the total score */
    /* Note: mat_pos is the current column of the matrix */
    for (mat_pos=0; mat_pos < matrix->length; mat_pos++) {

      /* make sure the sequence and the matrix overlap at this point */
      /* Note: scan_pos is where the current matrix column is in the */
      /*       sequence */
      scan_pos = seq_pos + mat_pos;
      if ((scan_pos >= 0) && (scan_pos < seq->length)) { /* if in the seq */
	seq_score += matrix->weights[ seq->sequence[scan_pos] ][mat_pos];
      }
      else {			/* not in the sequence */
	seq_score += matrix->weights[aa_atob['-']][mat_pos]; 
	/* score as if it was a gap */
      }

    } /* end score this alignment */

/*    alignments_done++; */
    
      if (seq_score > max_seq_score) {
	max_seq_score = seq_score;
	best_position = seq_pos;
      }
    
  } /* end of pairwise lineup */
  
  return(best_position);

}  /* end of best_pos */
/*========================================================================
        Figure out the consensus for this block & store it in the
        block list record - use highest PSSM score
        Prints x if the best total score is negative
=========================================================================*/
void pssm_consensus(bcur)
struct block_list *bcur;
{
   double maxw, weight;
   int aa1, seq1, pos, width, nseq, maxaa, aamark[MAXAA];
   Block *block;

   block = bcur->block;
   width = block->sequence_length;
   nseq = block->num_sequences;

   /*-------------------------------------------------------------------*/
   for (pos = 0; pos < width; pos++)
   { 
      /*    good sums can be negative */
      maxw = -999999.99;
      maxaa = -1;
      for (aa1 = 0; aa1 < MAXAA+1; aa1++)
         aamark[aa1] = NO;	/* does this aa appear in the pos ? */
  
      /*   Mark the aas that appear in this position */
      for (seq1 = 0; seq1 < nseq; seq1++)
      {
            aa1 = block->residues[seq1][pos];
            aamark[aa1] = YES;
      }

      /*    Only choose among the aas that appear in the pos  */
      for (aa1 = 0; aa1 < MAXAA; aa1++)
         if ( aamark[aa1] )
         {
            if (bcur->matrix->weights[aa1][pos] > maxw)
            {   maxw = bcur->matrix->weights[aa1][pos]; maxaa = aa1;  }
         }

      if (maxw <= 0.0) maxaa = MAXAA;	/*  Use X if non-positive sum */

      /*--------Save the consensus sequence------------------------*/
      if (maxaa > 0 && maxaa < MAXAA )
           bcur->consensus[pos] = maxaa;
      else
           bcur->consensus[pos] = MAXAA;	/* X==23 */
   }  /* end of pos */

}  /*  end of pssm_consensus */
/*=======================================================================*/
int sortcmp(t1, t2)
struct sorttemp *t1, *t2;
{
   return(t1->position - t2->position);
}  /* end of sortcmp */
