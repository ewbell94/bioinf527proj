/*======================================================================*/
/* (C) Copyright 1991 by Fred Hutchinson Cancer Research Center         */
/*     motomat.c   Block extension and assembly program                 */
/*   USE: motomat <input file> <cutoff score> <cluster %> <drop sd*10|asis>  */
/*                 EG:  motomat PS00094.mot 1 80 -10                    */
/*        input file = .mot file from motifj.c (this is a binary file)  */
/*        cutoff score = average PAM score for all possible pairs of    */
/*         amino acids in a column of the block (0-2500).               */
/*         If the cutoff score = 1, motomat will determine the best     */
/*         block extension. This is the recommended choice. If the      */
/*         cutoff score =2, lots of intermediate output is produced.    */
/*        cluster % = Percent identical AAs to cluster seqs in .blk file*/
/*        drop sd = sd of blocks score * 10, blocks scoring lower will  */
/*         be dropped and not assembled. EG -25 => -2.5 sds             */
/*         If drop sd == asis, blocks will not be merged, dropped or    */
/*           extended; duplicate blocks will be dropped                 */
/*        outputs blocks with the extension ".blk" and file name made   */
/*         up of the 1st 7 chars of the input file title plus one       */
/*         character for each assembled block; EG PS00094A.blk          */
/*   Reference:  S. Henikoff & J. Henikoff, NAR 19:23, 6565-6572 (1991)
      See PROTOMAT documentation file for description of algorithms.    */
/*----------------------------------------------------------------------
 6/20/90  J. Henikoff
>>>>>>>>>>>>>>>>>>>>>  Blocks 9.0  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 6/29/97 Changed execlp() to use argv[0]
10/21/97 Round proportion of sequences instead of truncate in best_path()
10/28/97 Added pb_weights() to score_block() & score_cols()
11/ 8/97 Make min drop score == highpass value
11/10/97 Weight blocks in best_path() by nident instead of nmotif
         Fixed bug in check_pos()
	 Added list of seqs not in best path to print_best()
11/11/97 Weight blocks in best_path() by nmotif from motifj
12/30/97 Added reduce_overlap()
 1/19/98 Weight blocks in best_path() by sqrt(nmotif) instead of nmotif.
 1/21/98 Fixed argv[0] error in restart(). Removed pb_weights() error message.
 1/25/98 Moved reduce_overlap() to after blocks are dropped for low score.
 1/28/98 Fixed reduce_overlap() bug, wrong block trimmed
         New option (-asis) => don't drop any blocks, don't extend,
		don't reduce overlap, etc.
 2/12/98.1 Added version date (VerDate[]). Drop blocks < 4 wide.
           Don't shorten block in reduce_overlap() if resulting block
	   will be < MIN_DOMAIN_WIDTH=10 wide.
           Fixed overcounting of block->nident in score_block()
	   Added Overlap = #times reduce_overlap() shortened a block to
		motomat.err
 2/16/98.1 Drop blocks < 5 wide.
 2/28/98.1 Modified trim_block() to never make block narrower.
==========================================================================*/
#include <math.h>
/*#include <process.h>*/
#include <unistd.h>
#include "motifj.h"

/*--------------------Function prototypes--------------------------------*/
void dump_data();
struct merged_motif *init_mm();
void merge_blocks();
int check_align();
int check_pos();
int check_overlap();
void dup_blocks();
void reduce_overlap();
void trim_block();
void extend_block();
void print_block();
void print_head();
void save_block();
void cluster_seqs();
void print_seq();
void write_seq();
void score_cols();
void score_block();
void pb_weights();
void restart();
/*    Shotgun assembly stuff;   code is in motomat2.c */
int prune_blocks();
void order_blocks();
int tempcmp();
void build_dag();
void find_paths();
void next_arc();
struct path *makepath();
struct block_list *makebllist();
void ins_path();
struct path *copypath();
void free_path();
struct matrix *makematrix();
void follow_arcs();
void ins_bllist();
struct path *best_paths();
void best_path();
void check_seqs();
void print_best();
void print_path();
void save_plot();
int left();
int right();
/*---------------Routines from motmisc.c -----------------------------*/
char *num_to_aachar();
int aachar_to_num();
void getscore();
struct split_name *split_names();

/*------------------Global vars------------------------------------------*/
char VerDate[12] = " 2/28/98.1";
int Debug = NO;		/* Flag for debugging print-out */
int AsIs = NO;		/* if argv[5] exists  */
int MinScore;		/* Minimum column score for block extension */
int ClThres;		/* Clustering threshold (%) */
int DropScore;		/* # of sds from mean before block is dropped */
long TopScore;		/* Actual drop score */
char Mot_Filename[FNAMELEN];    /* Name of the .mot input file */
char Blk_Filename[FNAMELEN];	/* Name of the .blk output file */
int Version, RunType, Signif, Dups, Distance;	/*  motifj parameters */
int RSignif;					/* reduces significance */
int NumSeqs;					/* number of sequences */
int Total_Motifs;				/* # motifs from motifj */
int Overlap;					/* # reduce_overlap()  */
char Title[MAXLINE], ctemp[MAXLINE];
int Gibbs = NO;		/* Gibbs input flag */
struct matrix *Dag[RELEVANT_MOTIFS][RELEVANT_MOTIFS];	/* Graph */
struct score *SMatrix;				/* Scoring matrix */
char Argv0[MAXLINE];		/* Save argv[0] for restart */

/*=======================================================================*/
void main(argc, argv)
int argc;
char *argv[];
{
   FILE *mot, *err;
   char temp[12], *ptr;
   int i, j, b, c, totb, motauto, maxdups;
   int lastscore, lastthres, lastdrop;
   int pseqs[MAXSEQS];
   double ftemp, fmean, fsdev;
   struct sequences *seqs;
   struct motif_struct *motif;
   struct merged_motif *blocks;
   struct path *paths;
   struct split_name *motsplit;

   printf("\nMOTOMAT %s: ", VerDate);
   printf("(C) Copyright 1991 by Fred Hutchinson Cancer Research Center");
   printf("\nPlease cite: S.Henikoff & J.Henikoff, NAR 19:6565-6572 (1991).");
/*-------Get the scoring matrix and high pass filter values --------*/
   SMatrix = (struct score *) malloc(sizeof(struct score));
   if (SMatrix != NULL)     getscore(SMatrix);
   else {printf("\nCANNOT ALLOCATE SMatrix\n"); exit(-1);}

   strcpy(Argv0, argv[0]);

/*----------------Get the input file -------------------------------*/
   if (argc > 1)
      strcpy(Mot_Filename, argv[1]);
   else
   {
      printf("\n");
      system("ls -al *.mot");
      printf("\nEnter name of file containing motifs: ");
      gets(Mot_Filename);
   }
   ptr = strrchr(Mot_Filename, '.');   /* check for extension */
   if (ptr == NULL) strcat(Mot_Filename, ".mot");
   if ( (mot=fopen(Mot_Filename, "rb")) == NULL)
   {
      printf("\nCannot open file %s\n", Mot_Filename);
      exit(-1);
   }
   motsplit=split_names(Mot_Filename);   /* split up name */
   Blk_Filename[0] = '\0';
   strncat(Blk_Filename, Mot_Filename+motsplit->dir_len, motsplit->name_len);
   Blk_Filename[motsplit->name_len] = '\0';

/*----------------Read the data file ------------------------------------*/
   Title[0] = '\0';
   NumSeqs = Total_Motifs = -1;
   fread(&Version, sizeof(int), 1, mot);
   if (Version != VERSION)
   {
      printf("\n\nInput file was created using wrong motifj version!\n");
      fclose(mot);
      exit(-1);
   }
   fread(&RunType, sizeof(int), 1, mot);
   fread(&Signif, sizeof(int), 1, mot);
   fread(&Dups, sizeof(int), 1, mot);
   fread(&Distance, sizeof(int), 1, mot);
   fread(&NumSeqs, sizeof(int), 1, mot);		/*  # sequences */
   fread(&Total_Motifs, sizeof(int), 1, mot);		/* # motifs */

   seqs = (struct sequences *) malloc(sizeof (struct sequences));
   seqs->num = NumSeqs;

   seqs->len = (int *) malloc(NumSeqs * sizeof (int));
   fread(seqs->len, sizeof(int), NumSeqs, mot);	/*  sequence lengths */

   seqs->totlen = 0;
   seqs->offlen = (int *) malloc(NumSeqs * sizeof (int));
   for (i=0; i<NumSeqs; i++)
   {
     *(seqs->offlen+i) = seqs->totlen;	/*  offset to ith sequence */
     seqs->totlen += *(seqs->len+i);
   }
					/* sequence names */
   seqs->name = (char *) malloc(SNAMELEN * NumSeqs * sizeof (char));
   fread(seqs->name, SNAMELEN*sizeof(char), NumSeqs, mot);

   seqs->seq = (char *) malloc(seqs->totlen * sizeof(char));
   fread(seqs->seq, sizeof(char), seqs->totlen, mot);	/* sequences */

   if (Total_Motifs > 0)
   {
      motif = (struct motif_struct *) 
          malloc(Total_Motifs * sizeof(struct motif_struct));
      /*---  For Gibbs, Total_Motifs is the maximum number, may be fewer */
      for (i=0; i<Total_Motifs; i++) motif[i].domain = -1;
      if (RunType == 9)				/* output from gibbs */
      {
         Gibbs = YES;
         if (!feof(mot))  fread(Title, 100*sizeof(char), 1, mot);
         printf("\nAssuming input from Gibbs\n");
         i = 0;
         while (!feof(mot) && i<Total_Motifs)
         {
            fread(&motif[i].domain, sizeof(int), 1, mot);
            fread(&motif[i].freq, sizeof(int), 1, mot);
            for (j=0; j< motif[i].freq; j++)
            {
               fread(&motif[i].seq_no[j], sizeof(int), 1, mot);
               motif[i].seq_no[j] -= 1;
               fread(&motif[i].pos[j], sizeof(int), 1, mot);
               motif[i].pos[j] -= 1;
            }
            motif[i].dups = motif[i].score = 0;
            motif[i].mots = 1;
            motif[i].distance1 = 0;
            motif[i].aa1 = motif[i].aa2 = motif[i].aa3 = 'X';
            motif[i].distance2 = motif[i].domain - 1;
            i++;
         }
         if (i < Total_Motifs) Total_Motifs = i;
      }
      else
      {
         fread(motif, sizeof(struct motif_struct), Total_Motifs, mot);
         if (!feof(mot))  fread(Title, MAXLINE*sizeof(char), 1, mot);
      }
   }
   else motif = NULL;
   fclose(mot);

   /*--------  Beware of incomplete records ---------------------*/
   i = 0;
   while (i < Total_Motifs && motif[i].domain > 0) i++;
   if (i < Total_Motifs) Total_Motifs = i;

/*---Fix title for use by print_best, save_block, etc.  -----------------*/
/*--   General form is >ACxxxxx ;ID;Description;$ -------*/
   if (strlen(Title) < 10 || Title[0] != '>')
   {
      strcpy(ctemp, Title);
      if (!strlen(ctemp)) strcpy(ctemp, "none");
      strcpy(Title, ">");
      strcat(Title, Blk_Filename);	/* use the filename as AC */
      if ( (strlen(Title)) < 9)
	 for (i=strlen(Title); i<9; i++) strcat(Title, "x");
      else if ( (strlen(Title)) > 9)     Title[9] = '\0';
      strcat(Title, ";none;"); strcat(Title, ctemp); strcat(Title, ";$");
   }
   /*------  There have to be some semicolons or other routines mess up ---*/
   else
   {
      if (strstr(Title,";") == NULL)
      {
	 strcpy(ctemp, Title+9); Title[9]='\0';
	 strcat(Title, ";none;"); strcat(Title, ctemp); strcat(Title, ";");
      }
      for (i=1; i<9; i++) if (Title[i] == ' ') Title[i] = 'x';
   }
   /*---  Strip off the MOTIFJ stuff from uextract, if present ------*/
   ptr = strstr(Title, "MOTIFJ=[");
   if (ptr != NULL) Title[strlen(Title)-strlen(ptr)] = '\0';
   printf("\n%s\n", Title);
   printf("\n%d blocks read", Total_Motifs);
   printf(" using parameters [%d, %d, %d].", Signif, Dups, Distance);
   printf(" %d sequences.", NumSeqs);

/*------------------Get the parameters ----------------------------------*/
   if (argc <= 3) RunType = 0;			/* ask for the args */
   else RunType = 3;				/* don't ask for args  */
   motauto = NO;  /* MinScore==1 is a flag for auto parameter selection  */
   if (argc > 2)
   {
      strcpy(temp, argv[2]);
      MinScore = atoi(temp);
   }
   if (MinScore == 2)
   {
      Debug = YES; MinScore = 1;
   }
   if (MinScore == 1)		/* Automatic determination of parameters */
   {
      motauto = YES;  RunType = 3;
   }
   else if (MinScore < 100 || MinScore > 2500) MinScore = MINSCORE;
   lastscore = MinScore;

   if (argc > 3)		/* Clustering threshold percentage */
   {
      strcpy(temp, argv[3]);
      ClThres = atoi(temp);
   }
   if (ClThres <= 0 || ClThres > 100) ClThres = CLTHRES;
   lastthres = ClThres;

   if (argc > 4)		/* Block drop score stand. dev. */
   {
      if (strcasecmp(argv[4], "asis") == 0) AsIs = YES;
      else DropScore = atoi(argv[4]);
   }
   else if (motauto)
   {
      DropScore = -10;		/* default = -1sd */
   }
   if (DropScore < -30 || DropScore > 30) DropScore = DROPSCORE;
   lastdrop = DropScore;

   printf("\nProcessing blocks in %s using a cut-off score of %d ",
		      Mot_Filename, MinScore);
   printf("\n   and a clustering threshold of %d percent.", ClThres);


/*----------------------Quit if nothing to do -------------------------*/
   if (Debug)
      dump_data(motif, seqs);
   if (Total_Motifs < 1)
   {
      printf("\nQuitting becasuse no motifs to work on.\n");
      exit(0);
   }

/*-------------------Redirect standard error output --------------------*/
   err = freopen("motomat.err", "a", stderr);

/*------------------Initialize sequences to save to all-----------------*/
   for (i=0; i<NumSeqs; i++)
      pseqs[i] = YES;

/*----------------------Create and merge the blocks --------------------*/
   blocks = init_mm(motif);
   free(motif);			/* Done with motifs now */

/*---------------------Merge motifs if not from Gibbs ------------------*/
   if (!Gibbs && !AsIs)
   {
      printf("\n\nMerging raw blocks...");
      merge_blocks(blocks, 0);
   }

/*-------Adjust the block width and recalculate the column scores -------*/
/*-------Do this over MAX_MERGE_WIDTH-----------------------------------*/
   printf("\nRecalculating block column scores");
   for (b=0; b<Total_Motifs; b++)
      if (blocks[b].nmotif > 0)
      {
	 if (!AsIs && blocks[b].domain+1 < MAX_MERGE_WIDTH)
	 {
	    blocks[b].loffset = blocks[b].t_loffset =
	     (int) (blocks[b].loffset +
		    (MAX_MERGE_WIDTH-blocks[b].domain-1)/2);
	    blocks[b].domain = blocks[b].t_domain = MAX_MERGE_WIDTH-1;
	 }
	 printf(".");
	 score_cols(&blocks[b], seqs);
      }


/*--------- Extend, cluster & print out the top merged blocks -----------*/
      printf("\nExtending blocks...");
      totb=0;				/* total surviving blocks */
      RSignif=Signif;			/* reduced significance */
      for (b=0; b<Total_Motifs; b++)
      {
	 printf("\nBlock %d: %.1s%.1s%.1s",
		 b, num_to_aachar(blocks[b].aa[0]),
		    num_to_aachar(blocks[b].aa[1]),
		    num_to_aachar(blocks[b].aa[2]));
	 if (blocks[b].nmotif <= 0) { printf(" has been merged.");}
	 else 
         {
            if (!AsIs)
	    {
	    do {
	       if (RunType == 0 || RunType == 2)
	       {
		  printf("\nEnter minimum block extension column score");
		  printf(" or 1 for best block [100-2500 or 1; %d]: ",
		    lastscore);
		  gets(temp);
		  MinScore = atoi(temp);
		  if ((MinScore > 1 && MinScore < 100) ||
		       MinScore <= 0 || MinScore > 2500)
			      MinScore = lastscore;
		  lastscore = MinScore;
		  printf("Enter clustering identity threshold [1-100; %d]: ",
		    lastthres);
		  gets(temp);
		  ClThres = atoi(temp);
		  if (ClThres <= 0 || ClThres > 100) ClThres = lastthres;
		  lastthres = ClThres;
		  printf("\nProcessing block %d in %s using a cut-off score of %d ",
		     b, Mot_Filename, MinScore);
		  printf("\n   and a clustering threshold of %d percent.", ClThres);
	       } /* end of runtype = 0 or 2 */
	       totb++;
	       if (MinScore == 1)
	       {
		  extend_block(&blocks[b]);
		  if (blocks[b].t_domain+1 < MIN_DOMAIN_WIDTH)
		  {
		     MinScore = SMatrix->highpass;
		     while (blocks[b].t_domain+1 < MIN_DOMAIN_WIDTH &&
			    MinScore >= (SMatrix->highpass - 50))
		     {
			printf("\nExtended width was %d",
			     blocks[b].t_domain+1);
			printf(", re-extending with cutoff score of %d.",
			     MinScore);
			trim_block(&blocks[b]);
			MinScore -= 50;
		     }
		     MinScore = 1;
		  }
		  if (blocks[b].t_domain+1 < 5)
                  { 
	             blocks[b].dropped = YES;
	             printf("\nBlock %d has been dropped for low width:", b);
	             print_head(&blocks[b]);
                  }
	       }  /* end of if MinScore==1 */
	       else
	       {
		   trim_block(&blocks[b]);
	       }

	       if (RunType==0 || RunType==2)
	       {
		   score_block(&blocks[b], seqs);
		   cluster_seqs(&blocks[b], seqs);
		   print_block(&blocks[b], seqs);
		   printf("\n\nType S[ENTER] to save this block,");
		   printf(" or N[ENTER] to skip to next block,");
		   printf("\nor anything else[ENTER] to re-extend this block...");
		   c=getchar(); getchar();    /*  2nd get is for ENTER */
		   if (c == 's' || c == 'S')
		   {
		      save_block(&blocks[b], seqs, pseqs);
		      c = 'N';
		   }
	       }
	    }  while ((RunType==0 || RunType==2) && c != 'N' && c != 'n');

	  }  /*  end of if !AsIs  */
        }  /* end of else not merged */
     }  /*  end of for b */
/*-----------Merge extended blocks if not from Gibbs -----------*/
   if (!Gibbs && !AsIs)
   {
      printf("\n\nMerging extended blocks....");
      merge_blocks(blocks, 1);
   }
   else
   {
      printf("\n\nRemoving duplicate blocks....");
      dup_blocks(blocks);
   }

/*---------------------------------------------------------------------*/
/*  Shotgun assembly of all unmerged, high-scoring blocks.             */
/*---------------------------------------------------------------------*/
   if (RunType == 0 || RunType == 2)
   {
      printf("\nEnter drop standard deviations [-30 to 30; %d]: ", lastdrop);
      gets(temp);
      DropScore = atoi(temp);
      if (DropScore < -30 || DropScore > 30) DropScore = lastdrop;
      lastdrop = DropScore;
   }
   printf("\nUsing a drop score of %d standard deviations or %d.",
	  (int) DropScore/10, SMatrix->highpass);
/*--------- Re-score blocks and compute mean and sd -----------*/
   totb = maxdups = 0;
   fmean = fsdev = 0.0;
   for (b=0; b<Total_Motifs; b++)
   {
      if (blocks[b].nmotif <= 0)
	  blocks[b].dropped = YES;
      else
      {
	 score_block(&blocks[b], seqs);
	 fmean += (double) blocks[b].t_score;
	 fsdev += (double) blocks[b].t_score*blocks[b].t_score;
	 totb++;
      }
   }  /*  end of for b */
/*-------- Compute mean and sd  ----------------*/
   if (totb > 1)
   {
      fsdev = fsdev - (double) fmean*fmean/totb;
      fsdev = (double) fsdev/(totb-1);
      if (fsdev > 0.0)
         fsdev = sqrt(fsdev);			/* standard deviation */
      else fsdev = 0.0;
      fmean = (double) fmean/totb;			/* mean */
   }
   else
      fsdev = 0.0;
   printf("\nMean score=%.2f, standard deviation=%.2f", fmean, fsdev);
/*-----If sdev==0.0, then either all blocks have same score, or there
       is only one block; don't want to drop them all ------------------*/
   if (fsdev > 0.01) ftemp = fmean + fsdev * DropScore / 10;
   else ftemp = 1.0;
   TopScore = (long) ftemp - 1.0;	            /*  Round down */
   /*    Possibly raise the drop score, but be careful not to drop
	 all the blocks!  */
   if (fmean >= SMatrix->highpass &&
       TopScore < SMatrix->highpass) TopScore = SMatrix->highpass;
   if (!AsIs)
   {
      printf("\nDropping blocks with score below %ld...", TopScore);
      printf("\n unless they have more than one merged motif");
   }
   totb = 0;
   for (b=0; b<Total_Motifs; b++)
      if (!AsIs &&
          blocks[b].dropped == NO &&
	  blocks[b].nmotif == 1 &&
	  blocks[b].t_score < (int) TopScore)
      {
	  blocks[b].dropped = YES;
	  printf("\nBlock %d has been dropped for low score:", b);
	  print_head(&blocks[b]);
      }
      else if (blocks[b].dropped == NO)
      {
	  totb++;
	  if (blocks[b].dups >maxdups) maxdups = blocks[b].dups;
	  cluster_seqs(&blocks[b], seqs);
	  printf("\nBlock %d:", b);
	  print_block(&blocks[b], seqs);
      }

/*   if (totb > MAXBLK)
      totb = prune_blocks(blocks);
*/
   printf("\nSurviving blocks=%d, max dups=%d", totb, maxdups);
/*------------------------------------------------------------------*/
/*         Check each pair of surviving
           blocks and look at those that overlap. By this point,
           they don't overlap consistently in all the sequences or
           they would have been merged. However, if they overlap
           consistently in most of the sequences and if that overlap is
           small, then trim one or both blocks so they don't overlap
	   any more. (Addresses the PS00031 problem without adding
	   gaps). See motomat2.c lines 140-157.
*/
   Overlap = 0;
   if (!AsIs)
   {
      printf("\nReducing block overlaps ...\n");
      reduce_overlap(blocks, seqs);
   }

/*---- Reduce significance level if there are dups. However, be
  sure that RSignif > (NumSeqs-RSignif) to prevent cycles in the DAG ---*/
   RSignif = Signif - ((int) maxdups/4);
   if (RSignif <= (NumSeqs-RSignif) )
      RSignif = ((int) NumSeqs/2) + 1;  /* R must be more than half N */
   if (RSignif < 2) RSignif = Signif;	/* not sure what to do here */
   if (RSignif < Signif)
      printf(", reduced significance=%d", RSignif);

   if (totb > 0)			/* Don't try if no blocks left */
   {
      printf("\n\nLooking for best path...");
      order_blocks(blocks);
      build_dag(blocks);
      paths = makepath();
      find_paths(paths, blocks);
      printf(", npath=%d", paths->nblocks);
      if (paths->next_path != NULL)
      {
	  printf(", nseqs=%d\n", paths->next_path->nseqs);
	  if (maxdups > 0) check_seqs(paths->next_path, blocks);
	  print_path(paths->next_path);
	  print_best(paths->next_path, blocks, seqs);
      }
   }
   else
      printf("\nNo surviving blocks");

   printf("\n");
/*----------Print a summary of results in stderr--------------------*/
   fprintf(stderr, "%.8s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
      Title, NumSeqs, Signif, Dups, Distance, Total_Motifs, MinScore,
      ClThres, DropScore, totb, paths->nblocks, paths->next_path->nbest,
      paths->next_path->totmotif, paths->next_path->totscore,
      paths->next_path->naas, paths->next_path->nseqs, Overlap);
   fclose(err);
   exit(0);
}    /*  end of main */

/*======================================================================
     Print the motifs and sequences
========================================================================*/
void dump_data(motif, seqs)
struct motif_struct *motif;
struct sequences *seqs;
{
   int i, j, k, offset, jj, kk;

   printf("\nNumSeqs=%d, Total_Motifs=%d", NumSeqs, Total_Motifs);
   k=0;
   for (i=0; i<NumSeqs; i++)
   {
      printf("\nSequence #%d: %10s %d",
	      i, seqs->name+SNAMELEN*i, *(seqs->len+i) );
      if (i>0) k += *(seqs->len+i-1);
      printf(" %.5s", seqs->seq+k);
   }

/*------NOTE:  motif[].domain is the ACTUAL domain width, while in
      blocks[].domain it is the displacement, or ACTUAL-1 -----*/
   for (i=0; i<Total_Motifs; i++)
   {
      printf("\n\nMOTIF %d: %d(%d)%d(%d)%d", i, motif[i].aa1,
	motif[i].distance1, motif[i].aa2, motif[i].distance2, motif[i].aa3);
      printf("  domain=%d, score=%d, freq=%d, dups=%d, mots=%d",
	motif[i].domain, motif[i].score, motif[i].freq, motif[i].dups,
        motif[i].mots);
      offset = (int) (motif[i].domain-motif[i].distance1-motif[i].distance2-1)/2;
/*>>>>>> Confused about this!
      offset++;
>>>>>>>>*/
      for (j=0; j<NumSeqs+motif[i].dups; j++)
      {
	 jj = motif[i].seq_no[j];
	 kk = motif[i].pos[j] - offset;
	 printf("\n%10s (% 5d) ", seqs->name+SNAMELEN*jj, kk);
	 for (k=0; k<motif[i].domain; k++)
	 {
	    kk = motif[i].pos[j] + k - offset;
	    if (kk >= 0 && kk < *(seqs->len+jj) )
	       printf("%.1s", seqs->seq+seqs->offlen[jj]+kk);
	    else  printf(".");
	 }
	 printf("\n                   ");
	 for (k=0; k<motif[i].domain; k++)
	    if (k == offset)  printf("[");
	    else if (k == offset+motif[i].distance1) printf("|");
	    else if (k == offset+motif[i].distance1+motif[i].distance2)
		 printf("]");
	    else printf("_");
      }
   }
}   /*  end of dump_data */

/*=======================================================================*/
/*  init_mm creates an array of merged motifs, each consisting of one    */
/*   motif, from an array of motifs created by H. Smith's program.       */
/*=======================================================================*/
struct merged_motif *init_mm(motif)
struct motif_struct *motif;
{
   struct merged_motif *new;
   int m, s, ms;
   double ftemp;

   new = (struct merged_motif *) malloc(Total_Motifs * sizeof(struct merged_motif));
   if (new == NULL)
   {
      printf("\nmerged_motif: Unable to allocate merged_motif structure!\n");
      exit(-1);
   }
   for (m=0; m<Total_Motifs; m++)
   {
      new[m].dropped = NO;
      new[m].aa[0] =(motif[m].aa1);
      new[m].aa[1] =(motif[m].aa2);
      new[m].aa[2] =(motif[m].aa3);
      /* award a bonus if some motifs were merged in MOTIFJ */
/*      if (motif[m].mots > 1) new[m].nmotif = 2; 
*/
      new[m].nmotif = 1;
      new[m].nmotif = motif[m].mots;
      new[m].nident = 0;
      ftemp = (double) motif[m].domain;
      new[m].max_score = (int) (motif[m].score/sqrt(ftemp));
      /*--NOTE: block domain is displacement, not actual with  ---*/
      new[m].domain = new[m].t_domain = motif[m].domain - 1;
      new[m].distance = motif[m].distance1 + motif[m].distance2;
      new[m].dups = motif[m].dups;
      /*--- loffset is not really offset, but 1st position ???--*/
      /* NOTE: If the block width is the motif width, loffset is zero */
      new[m].loffset = new[m].t_loffset =
       (int) (motif[m].domain-motif[m].distance1-motif[m].distance2-1)/2;
      new[m].t_score = -1;       /* -1 => block has not been trimmed */
      for (s=0; s<NumSeqs; s++)
      {
	 new[m].leftpos[s] = -1;
	 new[m].cluster[s] = -1;
	 new[m].position[s] = -1;
      }
      new[m].maxpos = new[m].minpos = -1;
      new[m].in_degree = new[m].out_degree = -1;
      for (ms=0; ms<NumSeqs+motif[m].dups; ms++)
      {
	s = motif[m].seq_no[ms];
	if (new[m].leftpos[s] < 0)   /* Just take the 1st occurence */
	   new[m].leftpos[s] = motif[m].pos[ms];
      }
      for (ms=0; ms<motif[m].domain; ms++)
	 new[m].scores[ms] = motif[m].scores[ms];
   };

   return (new);
}  /* end of init_mm */

/*=====================================================================
   merge_blocks merges blocks if 1) ALL sequences are aligned the
     same in all blocks and 2) the blocks overlap.  Flag=0 means
     check overlap in motif region, flag=1 means check overlap anywhere
     in extended region.
     NOTE:  .t_loffset & .t_domain are used instead of .loffset and
	    .domain so the code will work for both raw (flag==0) and
	    extended (flag=1) blocks - since .t_loffset is initialized
	    to .loffset and .t_domain to .domain.
     If the total width of the merged blocks would exceed
     MAX_MERGE_WIDTH they are not merged. The amount of their overlap
     is divided between them to make two contiguous, non-overlapping
     blocks.
========================================================================*/
void merge_blocks(blocks, flag)
struct merged_motif *blocks;
int flag;
{
   int b1, b2, s, right, left, ldiff, rdiff, lb1, lb2, rb1, rb2;
   int overlap, ov1, ov2, bleft, bright, temp;

   for (b1=0; b1<Total_Motifs-1; b1++)
      if (blocks[b1].nmotif > 0)
	 for (b2=b1+1; b2 <Total_Motifs; b2++)
	 {
	    if (blocks[b2].nmotif > 0 &&
		check_align(&blocks[b1], &blocks[b2]) &&
		check_overlap(&blocks[b1], &blocks[b2], flag) )
	    {
/*--- lb1,2 and rb1,2 are the end positions of blocks b1 and b2
      respectively for sequence 0.  They are used to compute the relative
      positions of the two blocks.  This technique only works if 1) All
      sequences are aligned the same way in both blocks, and 2) sequence
      0 extends all the way across both blocks (?) ----*/
	       lb1 = blocks[b1].leftpos[0]-blocks[b1].t_loffset;
	       lb2 = blocks[b2].leftpos[0]-blocks[b2].t_loffset;
	       rb1 = lb1 + blocks[b1].t_domain;
	       rb2 = lb2 + blocks[b2].t_domain;
	       ldiff = lb1 - lb2;	      /*-- difference at left end*/
	       rdiff = rb1 - rb2;  	     /*-- difference at right end*/
	       left = UMIN(lb1, lb2);
	       right = UMAX(rb1, rb2);
/*--------  Merge only if the total width is not too big ----------*/
	       if ( (right-left) < MAX_MERGE_WIDTH )
	       {

/*---- Update the block information ---*/
		  printf("\nMerging motif %d into motif %d", b2, b1);
		  blocks[b1].nmotif = blocks[b1].nmotif + 1;
		  blocks[b2].nmotif = 0;
		  if (blocks[b2].dups > blocks[b1].dups)
		     blocks[b1].dups = blocks[b2].dups;
		  if (!flag) blocks[b1].domain = right - left;
		  blocks[b1].t_domain = right - left;

		  left = UMIN(blocks[b1].leftpos[0], blocks[b2].leftpos[0]);
		  right = UMAX(blocks[b1].leftpos[0]+blocks[b1].distance,
			   blocks[b2].leftpos[0]+blocks[b2].distance);
		  blocks[b1].distance = right - left;
		  blocks[b1].t_loffset = left - UMIN(lb1, lb2);
		  if (!flag) blocks[b1].loffset = left - UMIN(lb1, lb2);
		  for (s=0; s<NumSeqs; s++)
		     blocks[b1].leftpos[s] =
		       UMIN(blocks[b1].leftpos[s], blocks[b2].leftpos[s]);
	       }
/*----------------------------------------------------------------------*/
/*--- Instead of merging, divide the blocks into non-overlapping blocks */
	       else if (flag)
	       {
		  bleft = bright = -1;
		  overlap = temp = 0;
		  if (ldiff < 0 && rdiff < 0)  /*  b1 is leftmost */
		  {
		     overlap = 1 + rb1 - lb2;
		     bleft = b1; bright = b2;
		  }
		  else if (ldiff > 0 && rdiff > 0) /* b2 is leftmost */
		  {
		     overlap = 1 + rb2 - lb1;
		     bleft = b2; bright = b1;
		  }
		  if (overlap > 0 && bleft >= 0 && bright >= 0)
		  {
		     ov1 = (int) (overlap/2); ov2 = overlap - ov1;
		     blocks[bleft].domain = blocks[bleft].domain - ov1;
		     if ( (blocks[bleft].loffset + blocks[bleft].distance) >
			  blocks[bleft].domain)
			 blocks[bleft].distance =
			    blocks[bleft].domain - blocks[bleft].loffset;
		     blocks[bleft].t_domain = blocks[bleft].t_domain - ov1;

		     blocks[bright].domain = blocks[bright].domain - ov2;
		     blocks[bright].t_domain = blocks[bright].t_domain - ov2;
		     blocks[bright].loffset = blocks[bright].loffset - ov2;
		     if (blocks[bright].loffset < 0)
		     {
			 temp = blocks[bright].loffset;
			 blocks[bright].distance = blocks[bright].distance
			      + temp;
			 blocks[bright].loffset = 0;
		     }
		     blocks[bright].t_loffset = blocks[bright].t_loffset - ov2;
		     if (blocks[bright].t_loffset < 0)
			blocks[bright].t_loffset = 0;
		     for (s=0; s<NumSeqs; s++)
			blocks[bright].leftpos[s] =
			   blocks[bright].leftpos[s] - temp;
		     printf("\nBlocks %d and %d have been divided",
			     bleft, bright);
		  }
	       }
	    }
	 }
}  /*  end of merge_blocks */
/*=====================================================================
   check_align returns YES if all the sequences in two blocks have the
     same alignment.
========================================================================*/
int check_align(b1, b2)
struct merged_motif *b1, *b2;
{
   int s, diff;

   s=0;
   diff = b1->leftpos[s]-b2->leftpos[s];
   do  {
      s++;
   }  while (s<NumSeqs &&
	     b1->leftpos[s]-b2->leftpos[s] == diff);
   if (s == NumSeqs) return(YES);
   else return(NO);
}  /*  end of check_align */

/*=====================================================================
   check_pos returns YES if all the sequences in two blocks have the
     same left position.
========================================================================*/
int check_pos(b1, b2)
struct merged_motif *b1, *b2;
{
   int s;

   s=0;
/*	Used to be
	     b1->leftpos[s] == b2->leftpos[s])
*/
   while (s<NumSeqs &&
	  b1->leftpos[s] - b1->t_loffset == b2->leftpos[s] - b2->t_loffset)
        s++;
   if (s == NumSeqs) return(YES);
   else return(NO);
}  /*  end of check_pos */

/*=====================================================================
   check_overlap returns YES if the first sequence in each of two
     blocks overlap anywhere within the blocks' motif area.
   Flag=0 means check for overlap in motif region; flag=1 means check
   anywhere in the extended region.
   May not catch a block that overlaps in only some of the sequences,
	but flag=2 forces it to check all sequences.
========================================================================*/
int check_overlap(b1, b2, flag)
struct merged_motif *b1, *b2;
int flag;
{
   int left1, left2, s, over;

   if (flag == 0)      /* check motif regions */
   {
/*   Find start of 1st sequence in both blocks  */
      left1 = b1->leftpos[0];
      left2 = b2->leftpos[0];
      if ( (left2 >= left1 && left2 <= left1 + b1->distance) ||
	   (left1 >= left2 && left1 <= left2 + b2->distance) )
	      return(YES);
      else    return(NO);
   }
   else if (flag == 1)        /*  check full regions  */
   {
      left1 = b1->leftpos[0] - b1->t_loffset;
      left2 = b2->leftpos[0] - b2->t_loffset;
      if ( (left2 >= left1 && left2 <= left1 + b1->t_domain) ||
	   (left1 >= left2 && left1 <= left2 + b2->t_domain) )
	      return(YES);
      else    return(NO);
   }
   else 	/* check ALL seqs in full regions */
   {
      s = 0; over = NO;
      while (s<NumSeqs && !over)
      {
         left1 = b1->leftpos[s] - b1->t_loffset;
         left2 = b2->leftpos[s] - b2->t_loffset;
         if ( (left2 >= left1 && left2 <= left1 + b1->t_domain) ||
	      (left1 >= left2 && left1 <= left2 + b2->t_domain) )
	      over = YES;
         s++;
      }
      return(over);
   }
}  /*  end of check_overlap  */

/*========================================================================
      dup_blocks eliminates blocks that are exact duplicates by setting
       the merged flag. Blocks are duplicates if they have the same
       width & same offsets for all sequences
=========================================================================*/
void dup_blocks(blocks)
struct merged_motif *blocks;
{
   int b1, b2;

   for (b1=0; b1<Total_Motifs-1; b1++)
      if (blocks[b1].nmotif > 0)
	 for (b2=b1+1; b2 <Total_Motifs; b2++)
	 {
	    if (blocks[b2].nmotif > 0 &&
                blocks[b1].domain == blocks[b2].domain &&
                   check_pos(&blocks[b1], &blocks[b2]))
                      blocks[b2].nmotif = 0;
          }
}  /* end of dup_blocks */

/*========================================================================
   reduce_overlap trims blocks that overlap only slightly but consistently
     in most sequences
=========================================================================*/
void reduce_overlap(blocks, seqs)
struct merged_motif *blocks;
struct sequences *seqs;
{
   int b1, b2, s, diff[MAXSEQS], bleft, bright, left, right, over[MAXSEQS];
   int samediff, minover;


   for (b1=0; b1<Total_Motifs-1; b1++)
   {
      if (blocks[b1].nmotif > 0)
      {
	 for (b2=b1+1; b2 <Total_Motifs; b2++)
	 {
//      printf("b2 %d/%d \n",b2,Total_Motifs);

	    if ( blocks[b2].nmotif > 0 &&
                 check_overlap(&blocks[b1], &blocks[b2], 2) &&
		 !check_align(&blocks[b1], &blocks[b2])        )
            {
            /*>>>  Find out how many overlap & by how much
                if it's most & they all overlap by the
                same amount, then trim the left one by reducing
                ->t_domain & re-score it  */
//               printf("b2 %d/%d -overlap \n",b2,Total_Motifs);

               samediff = YES;
               minover = 0;
               for (s=0; s<NumSeqs; s++)
               {
                  diff[s] = blocks[b2].leftpos[s] - blocks[b1].leftpos[s];
                  if (s > 0)
                  {
                     if ( (diff[0] <= 0 && diff[s] > 0) ||
                          (diff[0] >= 0 && diff[s] < 0)   )  samediff = NO;
                  }
                  if (diff[s] < 0) {  bleft = b2; bright = b1; }
                  else             {  bleft = b1; bright = b2; }
                  left =  blocks[bleft].leftpos[s] - blocks[bleft].t_loffset + blocks[bleft].t_domain;
                  right = blocks[bright].leftpos[s] - blocks[bright].t_loffset;
                  over[s] = right - left - 1;
                  if (over[s] < minover) minover = over[s];
               }  /* end of sequence s */
               /* Now check:  1) diff[s] same sign for all s
                  2) over[s] has a limited number of values, will be < 0
                  if there is an overlap: EG, if over = 0 for 2/10 seqs
                  and = -2 for 8/10, then trim 2 off the right end of bleft.
                  Count #different negative values in over[], if the
                  largest negative value is not too small, trim left block
                  by that amount. EG -10 is probably too much overlap,
                  but -5 may not be if left block is wide enough */
               if (samediff && minover < 0 && minover >= -5)
               {
                  /* Don't shorten if resulting block < 4 wide  */
                  if ( (blocks[bleft].t_domain + minover) > MIN_DOMAIN_WIDTH )
                  {
                     blocks[bleft].t_domain += minover;
		               score_block(&blocks[bleft], seqs);
		               Overlap++;
                     printf("Block %d has been trimmed by %d to eliminate overlap with Block %d.\n", 
                             bleft, minover, bright);
                  }
               }
            } /* end of if b1 and b2 overlap inconsistently */
          }  /* end of b2 */
      } /* end of if b1 is still active */
   }   /* end of b1 */
}  /*  end of reduce_overlap */

/*=======================================================================
    trim_block  trims a block around the motifs.                      
    It is mis-named, since what it actually does is EXTEND the block out
    from the outer edges of the merged motifs up to the full width of
    the block received from motifj. So it is the motifj block that is
    trimmed.							
=======================================================================*/
void trim_block(block)
struct merged_motif *block;
{
   int left, right;

   left = block->loffset;
	 do
	 {
	    left -= 1;
	 }  while (left >= 0 && block->scores[left] >= MinScore);

   right = block->loffset + block->distance;
	 do
	 {
	    right += 1;
	 } while (right <= block->domain &&
		  block->scores[right] >= MinScore);

   /* update block only if it is wider than before  */
   if ( (right - left - 2) > block->t_domain)
   {
      block->t_loffset = block->loffset - left - 1;
      block->t_domain = right - left - 2;
   }
}  /*  end of trim_block */
/*=======================================================================*/
/*    extend_block  extends a block out from the motifs.                 */
/*  Computes all possible blocks using all possible left and right
    extensions. Takes block with maximum block score = Smith's score
    divided by sqrt(width).                                              */
/*=======================================================================*/
void extend_block(block)
struct merged_motif *block;
{
   int left, middle, right, i, rn, ln, r, l, max_l, max_r;
   double ftemp, motif_sum, max_score, l_score, r_score, l_sum, r_sum;

   left = block->loffset;		   /* left edge of motif */
   middle = block->distance + 1;          /* width of motif */
   right = left + block->distance;         /* right edge  */
   motif_sum = 0;
   for (i=left; i<=right; i++)
      if (block->scores[i] >= SMatrix->highpass)
	 motif_sum += (block->scores[i]-SMatrix->highpass);
   ftemp = (double) middle;
   max_score = (double) ( (motif_sum/sqrt(ftemp)) );
   max_l = left;  max_r = right;

/*-----------  Extend right only -----------------------------*/
   r_sum = motif_sum; rn = middle;
   for (r=right+1; r<=block->domain; r++)
   {
     rn++;
     if (block->scores[r] >= SMatrix->highpass)
	   r_sum += (block->scores[r]-SMatrix->highpass);
     ftemp = (double) rn;
     r_score = (double) ( (r_sum/sqrt(ftemp)) );
     if (r_score > max_score)
     {
	max_score = r_score; max_r = r; max_l = left;
     }
   }

/*--Extend left & try all possible right extensions for each left -----*/
   l_sum = motif_sum; ln = middle;
   for (l=left-1; l>=0; l--)
   {
      ln++;
      if (block->scores[l] >= SMatrix->highpass)
		l_sum += block->scores[l]-SMatrix->highpass;
      ftemp = (double) ln;
      l_score = (double) ( (l_sum/sqrt(ftemp)) );
      if (l_score > max_score)
      {
	 max_score = l_score; max_l = l; max_r = right;
      }
      r_sum = l_sum; rn = ln;
      for (r=right+1; r<=block->domain; r++) /*--- all possible right ---*/
      {
	 rn++;
	 if (block->scores[r] >= SMatrix->highpass)
		 r_sum += (block->scores[r]-SMatrix->highpass);
	 ftemp = (double) rn;
	 r_score = (double) ( (r_sum/sqrt(ftemp)) );
	 if (r_score > max_score)
	 {
	    max_score = r_score; max_r = r; max_l = l;
	 }
      }
   }

/*------Now max_r and max_l are the best right & left extensions ------- */
   block->t_loffset = block->loffset - max_l;
   block->t_domain = max_r - max_l;
   block->t_score = max_score;
}  /*  end of extend_block */

/*=======================================================================
    Print out a trimmed block.
=========================================================================*/
void print_block(block, seqs)
struct merged_motif *block;
struct sequences *seqs;
{
   int s, i, j, x, sout;
   char tempname[30];

   print_head(block);
   if (RunType==0 || RunType == 2 || Debug)
   {
	 printf("\nScore              ");
	 for (i= block->loffset-block->t_loffset;
	      i<=block->loffset-block->t_loffset+block->t_domain;
	      i++)
	 {
	    x = (int) (block->scores[i]/100);
	    if (x >= 10) strcpy (tempname, "*");
	    else  kr_itoa(x, tempname, 10);
/*   Print sequence characters at the motif extremes
      Just takes 1st sequence value... */
	    if (i == block->loffset ||
		i == block->loffset+block->distance)
	    {
	       j = block->leftpos[0] - block->loffset + i;
	       strncpy(tempname, seqs->seq+seqs->offlen[0]+j, 1);
	    }
	    printf("%1s", tempname);
	 }

/*---- Print out in cluster order -------------*/
	 sout = 0;			/* # seqs printed for block */
	 for (s=0; s<NumSeqs; s++)
	    if (block->cluster[s] < 0)    /* unclustered first */
	    {
	       print_seq(block, seqs, s);
	       sout++;
	    }
	 i = 0;				/* cluster # */
	 while (sout < NumSeqs)
	 {
	    for (s=0; s<NumSeqs; s++)
	       if (block->cluster[s] == i)
	       {
		  print_seq(block, seqs, s);
		  sout++;
	       }
	    i++;
	 }
   }  /* end of RunType check */
}   /*  end of print_block */
/*======================================================================*/
void print_head(block)
struct merged_motif *block;
{
   printf("\n %.1s%.1s%.1s (extended width=%d, score=%d",
		num_to_aachar(block->aa[0]),
		num_to_aachar(block->aa[1]),
		num_to_aachar(block->aa[2]),
		block->t_domain+1, block->t_score );
   printf("  motifs=%d, motif width=%d, conserved=%d)",
		block->nmotif, block->distance+1, block->nident);
}  /* end of print_head */

/*=======================================================================
    Save a block.  Include the sequences flagged in pseqs.
=========================================================================*/
void save_block(block, seqs, pseqs)
struct merged_motif *block;
struct sequences *seqs;
int pseqs[MAXSEQS];			/* seqs to print */
{
    int s, i, sout, prevsout, pumseqs, best;
    FILE *blk;
    char blk_filename[FNAMELEN], partname[FNAMELEN], tempname[FNAMELEN];
    char *ptr, mem[MAXLINE];

/*---------------  Build the output file name --------------------------*/
/*  5/11/92 Changed to take name of .blk file from .mot file & to leave
	    .blk in current directory */
/*    blk_filename[0] = '\0';
    strncpy(tempname, Title+1, 9);       Use 1st 8 chars of title...
    ptr = strtok(tempname, " (;");	 Up to one of these chars
    strcpy(blk_filename, ptr);  strcpy(partname, ptr);
*/
    strcpy(blk_filename, Blk_Filename);
    /*  check for block letter */
    if (Title[8] != ' ' && Title[8] != '(' && Title[8] != ';')
    {
       i = strlen(blk_filename);
       strncat(blk_filename, Title+8, 1);	/* append block letter */
       blk_filename[i+1] = '\0';
    }
    strcpy(partname, blk_filename);	/* without extension */
    strcat(blk_filename, ".blk");
    printf("\nSaving block %.1s%.1s%.1s in file %s",
		num_to_aachar(block->aa[0]),
		num_to_aachar(block->aa[1]),
		num_to_aachar(block->aa[2]),
		blk_filename);
/*----------  Check to see if file already exists ---------------------*/
    if ((blk=fopen(blk_filename, "r")) != NULL)
    {
       printf("\nFile %s already exists,", blk_filename);
       if (RunType == 0 || RunType == 2)
       {
	  printf(" to avoid overwriting it");
	  printf("\n enter alternative filename [%s]: ", blk_filename);
	  gets(tempname);
	  if (strlen(tempname) > 1)
	  {
	      strcpy(blk_filename, tempname);
	      printf("\nSaving block %.1s%.1s%.1s in file %s",
		num_to_aachar(block->aa[0]),
		num_to_aachar(block->aa[1]),
		num_to_aachar(block->aa[2]),
		blk_filename);
	   }
       }
       else		/* In automatic modes, save the old file */
       {
	  strcat(partname, ".old");
	  printf(" copying it to %s\n", partname);
	  sprintf(mem, "cp %s %s", blk_filename, partname);
	  if (system(mem) != 0)
	     printf("\n    File copy was unsuccessful!");
       }
       fclose(blk);
    }
/*-------------- Open and write the output file -----------------------*/
    if ((blk=fopen(blk_filename, "w+t")) == NULL)
    {
	printf("\nCannot open %s", blk_filename);
	exit(-1);
    }
/*--------- Pseudo-EMBL; Parse Title -------------------------*/
/*  If title went through print_best, format is:  >ACa (x,y);ID;DE;$
				otherwise it is:  >AC ;ID;DE;$       */
    if ( strstr(Title, "(") != 0) best=YES; else best=NO;
    tempname[0] = partname[0] = '\0';
    strcpy(mem, &Title[1]);                   /* Skip Title[0] = '>' */
    ptr = strtok(mem, "; ");  if (ptr != NULL) strcpy(tempname, ptr);
    if (best)
    {
       ptr = strtok(NULL, ")");  if (ptr != NULL) strcpy(partname, ptr);
    }
    ptr = strtok(NULL,";");
    if (ptr != NULL) fprintf(blk, "ID   %s; BLOCK\n", ptr);
    else             fprintf(blk, "ID   none\n");
    fprintf(blk, "AC   %s;", tempname);
    if (best)  fprintf(blk, " distance from previous block=%s)\n", partname);
    else       fprintf(blk, "\n");
    ptr = strtok(NULL, ";");
    if (ptr != NULL) fprintf(blk, "DE   %s\n", ptr);
    else             fprintf(blk, "DE   none\n");
/*    fprintf(blk, "DO   \n");  */
    if (Gibbs)
    fprintf(blk, "BL   %.1s%.1s%.1s gibbs=[%d,%d,%d] motomat=[%d,%d,%d]",
		num_to_aachar(block->aa[0]),
		num_to_aachar(block->aa[1]),
		num_to_aachar(block->aa[2]),
	        Signif, Dups, Distance,
	        MinScore, ClThres, DropScore);
    else
    fprintf(blk, "BL   %.1s%.1s%.1s motif=[%d,%d,%d] motomat=[%d,%d,%d]",
		num_to_aachar(block->aa[0]),
		num_to_aachar(block->aa[1]),
		num_to_aachar(block->aa[2]),
	        Signif, Dups, Distance,
	        MinScore, ClThres, DropScore);
    fprintf(blk, " width=%d", block->t_domain+1);
/*-------------- Write sequences in cluster order --------------------*/
    pumseqs = 0;		/* # seqs to be printed for block */
    for (s=0; s<NumSeqs; s++)
       if (pseqs[s] == YES) pumseqs += 1;
    fprintf(blk, " seqs=%d\n", pumseqs);
    sout = 0;			/* # seqs printed for block */
    for (s=0; s<NumSeqs; s++)
      if (pseqs[s] == YES && block->cluster[s] < 0)   /* unclustered first */
      {
	  write_seq(blk, block, seqs, s);
	  fprintf(blk, "\n\n");
	  sout++;
      }
    i = 0;				/* cluster # */
    while (sout < pumseqs)
    {
       prevsout = sout;
       for (s=0; s<NumSeqs; s++)
	 if (pseqs[s] == YES && block->cluster[s] == i)
	 {
	    write_seq(blk, block, seqs, s);
	    fprintf(blk, "\n");
	    sout++;
	  }
	 if (sout != pumseqs && sout > prevsout)
	    fprintf(blk, "\n"); /* Extra line bw clusters*/
	 i++;
     }
     fprintf(blk, "//\n");
     fclose(blk);
}   /*  end of save_block */

/*======================================================================*/
/*    Cluster sequences in a trimmed block based on the number of       */
/*    identities within the block.                                      */
/*      1. Compute number of identities for each possible pair of seqs. */
/*         Results stored in lower half of matrix (pairs).              */
/*      2. Use clustering threshold % of # of AAs in trimmed block.     */
/*      3. Cluster recursively by traversing cols, rows of matrix.      */
/*======================================================================*/
void cluster_seqs(block, seqs)
struct merged_motif *block;
struct sequences *seqs;
{
   int nclus, npair, threshold, s1, s2, l1, l2, px, i, i1, i2;
   int oldclus, minclus;
   struct pair *pairs;

   npair = (NumSeqs*(NumSeqs-1))/2;
   pairs = (struct pair *) malloc (npair * sizeof(struct pair));
   if (pairs == NULL)
   {
      printf("\ncluster_seqs: Unable to allocate pair structure!\n");
      exit(-1);
   }
   threshold = (int) (ClThres*(block->t_domain+1))/100;

/*    Compute scores for all possible pairs of sequences            */
   for (s1=0; s1<NumSeqs-1; s1++)   		/* col = 0, n-2     */
   {
      l1 = block->leftpos[s1] - block->t_loffset;
      for (s2=s1+1; s2<NumSeqs; s2++)		/* row = col+1, n-1 */
      {
	 l2 = block->leftpos[s2] - block->t_loffset;
	 px = INDEX(NumSeqs, s1, s2);
	 pairs[px].score = 0;
	 pairs[px].cluster = -1;
	 for (i=0; i<=block->t_domain; i++)
	 {
	    i1 = l1+i;  i2 = l2+i;
	    if (i1 >= 0 && i1 < *(seqs->len+s1) &&
		i2 >= 0 && i2 < *(seqs->len+s2) &&
		strncmp(seqs->seq + seqs->offlen[s1] + i1,
			seqs->seq + seqs->offlen[s2] + i2, 1) == 0 )
		   pairs[px].score += 1;
	 }
      }  /* end of s2 */
   }  /* end of s1 */

/*  Print scores */
/*   printf("\nThreshold=%d", threshold);
   for (s2=1; s2<NumSeqs; s2++)
   {
      printf ("\n");
      for (s1=0; s1<s2; s1++)
      {
	 px = INDEX(NumSeqs, s1, s2);
	 printf(" %.3d", pairs[px].score);
      }
    }
*/

/*-------Cluster if score exceeds threshold by scanning cols (s1) */
   for (s1=0; s1<NumSeqs; s1++)
      block->cluster[s1] = -1;			/* clear out old values */
   nclus = 0;
   for (s1=0; s1<NumSeqs-1; s1++)   		/* col = 0, n-2     */
      for (s2=s1+1; s2<NumSeqs; s2++)		/* row = col+1, n-1 */
      {
	 px = INDEX(NumSeqs, s1, s2);
	 if (pairs[px].score >= threshold)
	 {
	    if (block->cluster[s1] < 0)	     /*  s1 has no cluster yet */
	    {
	       if (block->cluster[s2] < 0)    /* new cluster */
	       {
		  block->cluster[s1] = nclus++;
		  block->cluster[s2] = block->cluster[s1];
	       }
	       else  			/* s2 has a cluster, use it */
		  block->cluster[s1] =  block->cluster[s2];
	    }
            /* use s1's cluster if it has one and s2 doesn't */
            else if (block->cluster[s1] >= 0 && block->cluster[s2] < 0)
	       block->cluster[s2] = block->cluster[s1];
            /* merge the two clusters into the lower number */
            else if (block->cluster[s1] >= 0 && block->cluster[s2] >= 0)
            {
               minclus = block->cluster[s1]; oldclus = block->cluster[s2];
               if (block->cluster[s2] < block->cluster[s1])
               { minclus = block->cluster[s2]; oldclus = block->cluster[s1]; }
               for (i1=0; i1<NumSeqs; i1++)
                  if (block->cluster[i1] == oldclus)
                      block->cluster[i1] = minclus;
            }
	 }
      }   /* end of s2 */

   free(pairs);
}  /*  end of cluster_seqs */
/*======================================================================*/
/*  The blocks score is the sum of the column scores at or over HighPass
    divided by the square root of the block width.
    Scores only trimmed blocks. */
/*======================================================================*/
void score_block(block, seqs)
struct merged_motif *block;
struct sequences *seqs;
{
   int a, s, col, aa[20], pos, offset, flag;
   double ftemp, block_sum, ident_sum, weights[MAXSEQS];
   char c1, c[2];

   pb_weights(block, seqs, weights);		/* add sequence weights */

   block->nident = 0;				/* could be a re-count */
   block->t_score = 0;
   block_sum = ident_sum = 0.0;
   offset = block->loffset - block->t_loffset;
   for (col=0; col<=block->t_domain; col++)
   {
      flag = NO;			/* flag this col if an ident */
      for (a=0; a<20; a++)  aa[a] = 0;
      for (s=0; s<NumSeqs; s++)
      {
	  pos = block->leftpos[s] - block->t_loffset + col;
	  if (pos >= 0 && pos < *(seqs->len+s) )
	  {
	     strncpy(c, seqs->seq + seqs->offlen[s] + pos, 1);
	     c1 = c[0];
	     a = aachar_to_num(c1);
	     aa[a]++;
	  }
      }
      /*  Note: may be more than one aa in a col at least Signif times,
	  but are only counting such a column once */
      for (a=0; a<20; a++) if (aa[a] >= Signif) flag = YES;
      if (flag)
      {
          block->nident += 1;
	  ident_sum += (block->scores[col+offset] - SMatrix->highpass);
      }
      if (block->scores[col+offset] >= SMatrix->highpass)
	  block_sum += (block->scores[col+offset] - SMatrix->highpass);
   }  /*  end of col */
/*printf(" block_sum=%d, ident_sum=%d", (int) block_sum, (int) ident_sum); */
/*   ftemp = (double) block->t_domain + 1 + block->nident; */
   ftemp = (double) block->t_domain + 1;
   if (ftemp > 0)
      block->t_score = (int) ( block_sum/sqrt(ftemp) );
}  /* end of score_block */

/*======================================================================
	Put position-based sequence weights in the weights array
	col = block column, pos = sequence position
======================================================================*/
void pb_weights(block, seqs, weights)
struct merged_motif *block;
struct sequences *seqs;
double weights[MAXSEQS];
{
   struct pb_counts *pb;
   double factor, dtemp;
   int seq, pos, col, aa, width;
   char c1, c[2];

   width = block->t_domain;
   pb = (struct pb_counts *) malloc(width*sizeof(struct pb_counts));

   for (col = 0; col < width; col++)
   { 
      pb[col].diffaas = 0.0;
      for (aa = 0; aa < 20; aa++)
         pb[col].naas[aa] = (double) 0.0;
   }
  
   for (col = 0; col < width; col++)
      for (seq = 0; seq < NumSeqs; seq++)
      {
	 pos = block->leftpos[seq] - block->t_loffset + col;
	 if (pos >= 0 && pos < *(seqs->len+seq) )
	 {
	     strncpy(c, seqs->seq + seqs->offlen[seq] + pos, 1);
	     c1 = c[0];
	     aa = aachar_to_num(c1);
             if (aa >= 0 && aa <= 19)
             {
                pb[col].naas[aa] += 1;
             }
/*
             else
             {
                printf("pb_weights:%d ignored\n", aa);
             }
*/
         }
      }
   factor = 1.0;
   for (col = 0; col < width; col++)
   { 
      for (aa = 0; aa < 20; aa++)
      {
         if (pb[col].naas[aa] > 0.0)
         {
            pb[col].diffaas += 1;	/* # of different types of aas in col */
         }
      }
   }

   for (seq = 0; seq < NumSeqs; seq++)
   {
      weights[seq] = 0.0;
      for (col = 0; col < width; col++)
      {
	 pos = block->leftpos[seq] - block->t_loffset + col;
	 if (pos >= 0 && pos < *(seqs->len+seq) )
	 {
	     strncpy(c, seqs->seq + seqs->offlen[seq] + pos, 1);
	     c1 = c[0];
	     aa = aachar_to_num(c1);
             dtemp = pb[col].diffaas * pb[col].naas[aa];
             if (dtemp > 0.0) weights[seq] += 1.0 / dtemp;
         }
      } 
   }
 
   /*  Normalize weights to add to NumSeqs  */
   dtemp = 0.0;
   for (seq = 0; seq < NumSeqs; seq++) dtemp += weights[seq];
   for (seq = 0; seq < NumSeqs; seq++) 
	weights[seq] *= ((double) NumSeqs/dtemp);

   free(pb);
}  /* end of pb_weights */

/*======================================================================*/
/*    Print a sequence                                                  */
/*======================================================================*/
void print_seq(block, seqs, s)
struct merged_motif *block;
struct sequences *seqs;
int s;
{
   int left, right, i;

   left = block->leftpos[s] - block->t_loffset;
   right = left + block->t_domain;
   printf("\n%10s (% 5d) ", seqs->name+SNAMELEN*s, left+1);
   for (i=left; i<=right; i++)
   {
      if (i >= 0 && i < *(seqs->len+s) )
	  printf("%.1s", seqs->seq+seqs->offlen[s]+i );
      else
	  printf(".");
   }
   printf(" % 3d", block->cluster[s]);
}  /*  end of print_seq */

/*======================================================================*/
/*   Write  a sequence to an output file                                */
/*======================================================================*/
void write_seq(out, block, seqs, s)
FILE *out;
struct merged_motif *block;
struct sequences *seqs;
int s;
{
   int left, right, i;

   left = block->leftpos[s] - block->t_loffset;
   right = left + block->t_domain;
   /*  Print out the actual starting position of the sequence in
       the block, not the offset to it  */
   fprintf(out, "%10s (% 5d) ", seqs->name+SNAMELEN*s, left+1);
   for (i=left; i<=right; i++)
   {
     if (i >= 0 && i < *(seqs->len+s) )
	fprintf(out, "%.1s", seqs->seq+seqs->offlen[s]+i );
     else
	fprintf(out, "X");
   }
}  /*  end of write_seq */

/*======================================================================*/
/*   Calculate column scores for a block.
     Looks at all possible pairs of sequence AAs in a block column
     and takes average score. If one or both sequences has no AA in
     a column, because the block goes off one end or the other of the
     sequence, uses the score for 'J' */
/*======================================================================*/
void score_cols(block, seqs)
struct merged_motif *block;
struct sequences *seqs;
{
   int col, s1, s2, sum, npair, col1, col2;
   double weights[MAXSEQS];
   char c1, c2, c[2];

   pb_weights(block, seqs, weights);		/* add sequence weights */

   npair = (int) NumSeqs * (NumSeqs-1) / 2;  /* all possible pairs */
   for (col=0; col<=block->domain; col++)
   {
      sum = 0;				/* sum of all AA pairs */
      for (s1=0; s1<NumSeqs; s1++)
      {
	 col1 = block->leftpos[s1] - block->loffset + col;
	 if (col1 >= 0 && col1 < *(seqs->len+s1) )
	 {
	    strncpy(c, seqs->seq+seqs->offlen[s1]+col1, 1);
	    c1 = c[0];
	 }
	 else { c1 = 'X';}		/* off end of seq 1 */
	 for (s2=s1+1; s2<NumSeqs; s2++)
	 {
	    col2 = block->leftpos[s2] - block->loffset + col;
	    if (col2 >= 0 && col2 < *(seqs->len+s2) )
	    {
	       strncpy(c, seqs->seq+seqs->offlen[s2]+col2, 1);
	       c2 = c[0];
	       sum += SMatrix->scores[aachar_to_num(c1)][aachar_to_num(c2)];
	    }
	    else { c2 = 'X';}		/* off end of seq 2 */
	 }  /* end of s2 */
      }  /* end of s1 */
      block->scores[col] = (int) ((long int) sum * 100 / npair);
   }  /* end of col */
}  /*  end of score_cols */
/*======================================================================*/
/*   Increase drop score and restart motomat.                           */
/*======================================================================*/
void restart()
{
   char arg2[5], arg3[5], arg4[5];
   int row, col;

   printf("\n>>>MOTOMAT out of memory!");
   if (DropScore <= 28)
   {
      DropScore += 2;
      printf("  Restarting with DropScore = %d\n", DropScore);
      kr_itoa(MinScore, arg2, 10);
      kr_itoa(ClThres, arg3, 10);
      kr_itoa(DropScore, arg4, 10);
/*---  Have to free enough memory to reload motomat ----------------------*/
      for (row=0; row<Total_Motifs; row++)
	 for (col=0; col<Total_Motifs; col++)
	    if (Dag[row][col] != NULL)  free(Dag[row][col]);
      execlp(Argv0, Argv0, Mot_Filename, arg2, arg3, arg4, NULL);
      exit(0);
   }
   else exit(0);
}  /* end of restart  */
/*======================================================================*/
#include "motomat2.c"
