/*------------------------------------------------------------------------*/
/*                      Motif identifier, Version 1.12.1                  */
/*									  */
/*                       Thomas Annau, Hamilton Smith                     */
/*                             August 19, 1989                            */
/*  Reference:  Hamilton O. Smith, Thomas M. Annau, Srinivasan
	Chandrasegaran, "Finding sequence motifs in groups of
	functionally related proteins", PNAS, Vol 87 (January 1990),
	pp. 826-830.                                                      */
/*             Modified 2/9/90 by Ken Fasman, BME Systems Inc. for        */
/*                     Microsoft & draft ANSI standard C                  */
/*-----------------------------------------------------------------------
    NOTES for Silicon Graphics compilers:  Change EOF to '\377'           */
/*------------------------------------------------------------------------
       motifj RunType Batch_Filename Signif Dups Distance Drop           
	Defaults: Signif = (int) ((NumSeqs+1)/2) + 1
		  Dups = 0
		  Distance = 17
		  Drop = 0
*/
/*------------------------------------------------------------------------*/
/* 3/21/90 Modified for PROTOMAT by J. Henikoff, Fred Hutchinson Cancer
           Research Center (JGH)
   USE: Interactive:  motifj (prompts for responses)
        Command line: motifj RunType lisname [Signif Dups Distance Drop]
         RunType=0,1,2,3 or 4. Recommended types are:
           0 interactive, non-iterative (similar to original motif program)
           1 non-interactive, non-iterative,
             all parameters specified on command line
	   2 interactive, iterative (will decrease Signif until >= MOTAUTO3
             motifs are found)
           3 non-interactive, iterative (non-interactive version of 2)
           4 shuffled, non-interactive, iterative (will increase Signif
             until <= MOTAUTO4 motifs are found, then executes RunType 3),
             all parameters determined by motifj
           5 shuffled, non-interactive, iterative (like RunType 2/3, not 4)
           See PROTOMAT documentation for more information.
         lisname=name of file containing list of proteins to analyze.
           First line can be a title starting with ">"
           Second line can be a directory, EG /home/proteins/PS00094/,
	    if missing the protein files are assumed to be in pros/.
           Assumes each subsequent line lists the name of a file which
            has extension ".pro". So if lisname contains "MTB1_BREEP",
            the file name must be "MTB1_BREEP.pro".
           All protein files must be in fasta format with a title line
	    beginning ">".
           Proteins must be in upper or lower case letters, the usual 20 amino
            acid abbreviations are recognized, plus B (changed to D) and Z
            (changed to E). X, O and J are changed to X and treated as
            place holders. All other characters are ignored.
           If lisname filename is preceded by "-", it is assumed that
	    the file contains all of the proteins to be analyzed in
            fasta format.
	   PLEASE VERIFY THAT YOUR PROTEINS ARE BEING READ CORRECTLY!!!
         Signif=Motif significance level, required number of proteins
           containing a motif.
         Dups=Maximum number of duplications of a motif among all sequences.
         Distance=Maximum motif distance. 
         Drop=drop score, motif score below which motifs will be dropped.

 >>>>>>>>>>>>>>>>>>>>>>>>Blocks 9.0 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 10/21/95 Write out sequences used in .motifj.pros file, write_pros()
 10/23/95 Flag min/max sequence lengths in .motifj.pros file
 11/14/95 Leave Title intact (don't remove MOTIFJ=[]) for mot2blk
  7/11/96 Changed to re-execute program in argv[0] instead of "motifj"
  6/29/97 Execute write_pros() right away.
 10/21/97 Added pb_weights()
 10/28/97 Modified cut_motifs() to take actual overlap into account;
          lowered min #mots for runtype 3 (MOTAUTO3) from 12 to 6 because
          now expect fewer motifs. Raised max #mots for runtype 4 (MOTAUTO4)
          from 2 to 3 to reduce running time.
 11/ 4/97 Modified cut_motifs() to group two motifs only if they match
	  in at least half of their common sequences 
 11/ 7/97 Modified cut_motifs() to drop excess motifs in the same group.
          Expect to now produce fewer motifs, but spread across more
          of the sequence length (groups are roughly equivalent to a
	  geographical region of the sequences).
  1/21/98 Removed pb_weights() error message.
*/
/*--------------------------------------------------------------------JGH*/

#include <time.h>
#include <math.h>
#include "motifj.h"

struct realign_struct {
	int sequence[MAXSEQS];
	int motif_pos[MAXSEQS];
	long int max_score[MAXSEQS], max2_score[MAXSEQS];
	int max_pos[MAXSEQS], max2_pos[MAXSEQS];
};
/*
Function Prototypes
*/
void main();
void print_motif();
int pamidx();
int wpamidx();
int cut_motifs();
void map_seqs();
int compare_scores();
int compare_groups();
int compare_motgrps(); 
int compare_subgrps(); 
FILE *getfile();
int getseqs();
int shuffle();
void print_stats();
int get_count();
int check_dup();
void score_mot();
void pb_weights();
void realign();
void write_motifs();
void write_pros();

/*----------These routines are in motmisc.c ----------------------------*/
int aachar_to_num();
char *num_to_aachar();
void pr_num_to_aa();
void pr_num_to_aa_space();
void getscore();
struct split_name *split_names();

/* --------   Global variables ---------------------------------------------*/
char Version[12] = "11/ 7/97.2";
int Short_Form, First;	
int NumSeqs, RunType, Signif, Dups, Distance;
unsigned int Drop;
int Len[MAXSEQS] = {0};
char Seqname[MAXSEQS][SNAMELEN];
int *Seq[MAXSEQS];
struct score *SMatrix;
char Batch_Filename[FNAMELEN], Mot_Filename[FNAMELEN];
char Title[MAXLINE];
struct split_name *Batsplit;

/*---------------------------------------------------------------------------*/
/*                             MAIN PROCEDURE                                */
/*---------------------------------------------------------------------------*/

void main(argc, argv)
int argc;
char *argv[];
{
/*   NOTE: Be careful of integers defined as char to save space:
	    char is        8bits, valid range is -128 to +127
	    unsigned char  8                        0 to +255
	    int            16                  -32768 to +32767
	    unsigned int   16                       0 to +65535
	    long           32
*/
int i, j, k, n;
int batch, shuffled, prevsig, prevdup, save_motifs, save_dups;
int total_motifs, real_total_motifs, found, duplicate;
int count, minseqs, d1, d2, group, max_group, domain_width;
unsigned long total_score = 0;
long int CPU_time;
char *ptr, intemp[MAXLINE], seq_name[FNAMELEN], filepath[FNAMELEN]; 
char aa1, aa2, aa3;
char chsig[5], chdup[5], chdis[5], chdrop[6]; 
FILE *fp, *bf;
struct motif_struct *motif;
aa_type *aa3_obs;

   printf("\nMOTIFJ Version %s", Version);
   printf("\nPlease cite Hamilton O. Smith, et al,");
   printf(" PNAS 87 (1990), pp. 826-830.\n");

   /*---------Get the scoring matrix and high pass filter values ------------*/
   SMatrix = (struct score *) malloc (sizeof(struct score));
   if (SMatrix != NULL)    getscore(SMatrix);
   else {printf("\nCANNOT ALLOCATE SMatrix\n"); exit(-1);}

   /*--------Determine the run type -----------------------*/
   RunType = 0;	
   if (argc > 1) RunType = atoi(argv[1]);
   if (RunType < 0 || RunType > 5)  RunType = 0;	
   printf("\nUsing RunType=%d", RunType);

   /*------------------------------------------------------------------------*/
   /*              Get amino acid sequences from files.                      */
   /*------------------------------------------------------------------------*/
   Batch_Filename[0] = intemp[0] = '\0';
   if (argc > 2) strcpy(intemp, argv[2]);
   else 
   {
     printf("\n\nEnter input file name: ");
     gets(intemp);
   }
   if (strlen(intemp) && intemp[0] == '-')
   {
      batch = NO;
      strcpy(Batch_Filename, intemp+1);
   }
   else
   {
      batch = YES;
      strcpy(Batch_Filename, intemp);
   }

   /* Open input file */
   Batsplit = split_names(Batch_Filename);
   if ((bf = fopen(Batch_Filename, READ)) == NULL)
   {
      printf("Cannot open file %s not found\n", Batch_Filename);
      exit(-1);
   }

   if (batch)
     printf("\nReading batch file %s", Batch_Filename); 
   else
     printf("\nReading sequence file %s", Batch_Filename);

   NumSeqs = 0;
   Title[0] = '\0';
   if (batch)
   {
     /*    Process the batch file title, if there is one */
     fgets(Title, MAXLINE, bf);
     if (strlen(Title) && Title[0] != '>')
     {
	Title[0] = '\0';
	rewind(bf);
     }
     else
     {
	printf("\n%s", Title);
     }
     /* Read in all sequences (go through batch file)  */
     //printf("\nFilepath: %s", filepath);

     
     strcpy(filepath, PROTEIN_SUBDIRECTORY);      /* Initialize path */
     while (!feof(bf) && fgets(intemp, MAXLINE, bf) != NULL)
     {
        fp = getfile(intemp, filepath, Seqname[NumSeqs]);
        //printf("(%s)\n", );
        if (fp != NULL)
        {
	   n = getseqs(fp);		/* updates NumSeqs, etc.*/
    	   fclose(fp);
        }  /* end of if fp is not NULL */
     }  /* end of bf file */
   }  /* end of if batch */
   else
   {
     /*   Try to find a Title */
     intemp[0] = '\0';
     strncat(intemp, Batch_Filename,
             Batsplit->dir_len + Batsplit->name_len);
     intemp[Batsplit->dir_len + Batsplit->name_len] = '\0';
     strcpy(filepath, intemp); strcat(filepath, ".lst");
     if ((fp = fopen(filepath, READ)) == NULL)
     {
        strcpy(filepath, intemp); strcat(filepath, ".lis");
        fp = fopen(filepath, READ);
     }
     if (fp != NULL)
     {
        fgets(Title, MAXLINE, fp);
        fclose(fp);
     }
     n = getseqs(bf);				/* updates NumSeqs */
   }
   
   fclose(bf);
   if (NumSeqs < MINSEQS)
   {
      printf("\nToo few sequences to process(need at least 2).\n");
      exit(0);
   }
   else
   {
      printf("\n");
      for (n=0; n< NumSeqs; n++)
      {
	 printf("Sequence #%d:\t%s\t\t", n + 1, Seqname[n]);
    	 printf("Length of sequence #%d =  %d\n", n + 1, Len[n]);
      }
      printf("\n");
   }
   /*    Write out the sequences for Gibbs */
   write_pros();

   /*------------------------------------------------------------------------*/
   /*  If shuffling was selected, shuffle the sequences before analysis.     */
   /*------------------------------------------------------------------------*/
   shuffled = shuffle();

   /*------------------------------------------------------------------------*/
   /*        Perform brief statistical analysis on amino acid data.          */
   /*------------------------------------------------------------------------*/
   print_stats();

   /*---------------------------------------------------------------------*/
   /*     Get the parameters                                              */
   /*---------------------------------------------------------------------*/
   /*      Initialize the parameters        */
   Signif = Dups = Distance = Drop = 0;
   if (strlen(Title))
   {
	ptr = strstr(Title, "MOTIFJ=[");
	if (ptr != NULL)
        {
            strcpy(intemp, ptr);
/*          Title[strlen(Title)-strlen(ptr)] = '\0';
*/
            /*  Also look at title line for RunType 3 & 5 ? */
            if (RunType == 1 || RunType == 4)
            {
            /*  Read parameters off the title line */
	    /*  Assumes rigid format from uextract is:
                  MOTIFJ=[RunType,Signif,Dups,Distance]*/
            ptr = strtok(intemp, "["); ptr = strtok(NULL, ","); /* RunType */
	    ptr = strtok(NULL, ","); Signif = atoi(ptr);
	    ptr = strtok(NULL, ","); Dups = atoi(ptr);
	    ptr = strtok(NULL, "]"); 	/* don't use Distance */
	    Distance = 17;
            }
        }
   }
   //minseqs = (int) ((NumSeqs+1)/2) + 1;    /* want more than half of seqs */

   minseqs = (int) ((NumSeqs+1)/2) + 1;    /* want more than 1/3 of seqs */
   printf("minseqs: %d\n",minseqs);
   if (argc > 3) Signif = atoi(argv[3]);
   else if (RunType == 1)
   {
      Signif = (int) ((NumSeqs+1)/2) + 1;
   }
   else if (RunType == 0 || RunType == 2)
   {
     printf("\nEnter significance level [2-%d; %d]: ", NumSeqs, NumSeqs); 
     gets(intemp);   Signif = atoi(intemp); 
   }
   else if (RunType == 4 || Signif == 0)
     Signif = minseqs;	
   if ((RunType == 3 || RunType == 4 || RunType == 5)
       && Signif < minseqs) Signif = minseqs;
   if (Signif == 0 || Signif > NumSeqs) Signif = NumSeqs;
   if (Signif < MINSEQS)                Signif = MINSEQS;

   if (argc > 4) Dups = atoi(argv[4]);	
   else if (RunType == 0 || RunType == 2) 
   {
     Dups = MAXFREQ - NumSeqs;
     printf("Enter number of allowable internal duplicates [0-%d; 0]: ",
	     Dups);
     gets(intemp);  Dups = atoi(intemp); 
   }
   if (Dups + NumSeqs > MAXFREQ) Dups = MAXFREQ-NumSeqs;
   if (Dups < 0) Dups = 0;
   save_dups = Dups;
   if (shuffled) Dups = 0;

   if (argc > 5) Distance = atoi(argv[5]);
   else if (RunType == 0 || RunType == 2)
   {
     Distance = -1;
     printf("Enter search width [1-24; 17]: "); /* Distance=0 causes abort*/
     gets(intemp);
     if (strlen(intemp)) Distance = atoi(intemp);
     else                Distance = 17;
   }
   else Distance = 17;
   if (Distance <= 0 || Distance > MAX_DISTANCE) Distance = MAX_DISTANCE;

   if (argc > 6) Drop = atoi(argv[6]);
   else if (RunType == 0 || RunType == 2)
   {
     Drop = 0;			/* 18540 = 2500*55/sqrt(55)  */
     printf("Enter drop score [0-18540; %d]: ", Drop, Drop);
     gets(intemp);
     if (strlen(intemp)) Drop = atoi(intemp);
   }
   else Drop = 0;
   if (Drop <=0 || Drop > 18540) Drop = 0;

/*-------------------------------------------------------------------------*/
/*		 Allocate memory to dynamic structures.                    */
/*-------------------------------------------------------------------------*/
  printf("Allocating and initializing arrays...\n");

  aa3_obs = (aa_type *) malloc(sizeof(aa_type));
  if (aa3_obs == NULL) {
	printf("Couldn't allocate AA3_OBS pointer.\n");
	exit(-1);
	}
  for (i = 0; i < 20; i++) {
    for (j = 0; j < 20; j++) {
      for (n = 0; n < Distance; n++) {          /* was MAX_DISTANCE   JGH*/
	(*aa3_obs)[i][j][n] = (unsigned char *) malloc(Distance);   /*JGH*/
	  if ((*aa3_obs)[i][j][n] == NULL) {
		printf("Couldn't allocate AA3_OBS array.\n");
		exit(-1);
		}
	for (k = 0; k < Distance; k++)		/* was MAX_DISTANCE   JGH*/
	  (*aa3_obs)[i][j][n][k] = 0;
	}
      }
    }

  motif = (struct motif_struct *)
	      malloc(sizeof(struct motif_struct) * MAX_MOTIFS);
  if (motif == NULL) {
	printf("Couldn't allocate MOTIF.\n");
	exit(-1);
	}

/*---------------------------------------------------------------------------*/
/*   Find motifs by looking for three amino acids separated by identical     */
/*               distances in Signif number of proteins.                     */
/*---------------------------------------------------------------------------*/

REDO: CPU_time = time(0);			/* Record starting time */
    printf("\nRunType = %d, Normalized drop score = %d", RunType, Drop);
    printf("\nMotif parameters = [%d, %d, %d]", Signif, Dups, Distance);
    if (shuffled) printf(" with shuffled sequences.\n");
    else          printf("\n");
    total_motifs = real_total_motifs = 0;

  /* Step through the 20 amino acids as the first amino acid in a motif.
     Use a matrix to tabulate frequencies of all possible motifs built on
     that first amino acid.  The maximum spacing between amino acids of
     the motif is limited by the user's entry of the Distance parameter.

     aa1	First amino acid in motif.
     d1		Distance separating first aa from second.
     aa2	Second amino acid in motif.
     d2		Distance separating second aa from third.
     aa3	Third amino acid.

     (*aa3_obs)[aa2][aa3][d1][d2]  is the frequency of the motif
				   (matrix is redone for each aa1).
     */

  /* Start loop for first amino acid. */
  for (aa1 = 0; aa1 < 20; aa1++)               /* First amino acid */
  {
    printf(".");                                /* Look busy */

    /*>>>> RunType=4:  Only continue if haven't found too many yet */
    if (RunType != 4 ||
       (RunType == 4 && total_motifs <= MOTAUTO4))
    {
/*----------------------------------------------------------------------------*/
    /* Tabulate frequency matrix for all possible motifs within limits
       defined by user (d1 and d2 range from 0 to Distance): */
      for (n = 0; n < NumSeqs; n++)               /* Look through all seqs */
        for (i = 0; i < Len[n]; i++)              /* Step through each seq */
	  if (Seq[n][i] == aa1)                /* Search for aa1 */
          {
	  /* Found aa1, now assimilate Distance*Distance number of motifs
	     into array: */
	    for (j = i + 1; j < Len[n] && j < i + Distance + 1; j++)
	      for (k = j + 1; k < Len[n] && k < j + Distance + 1; k++)
                if (Seq[n][j] < 20 && Seq[n][k] < 20)  /* count only AAs */
	        /* (*aa3_obs) [  aa2    ] [  aa3    ] [ d1  ] [ d2  ]++ */
	           (*aa3_obs) [Seq[n][j]] [Seq[n][k]] [j-i-1] [k-j-1]++;
	  }
/*----------------------------------------------------------------------------*/

      /* Scan array for motifs that occurred at least Signif number of times */
      for (aa2 = 0; aa2 < 20; aa2++)
        for (aa3 = 0; aa3 < 20; aa3++)
	  for (d1 = Distance - 1; d1 >= 0; d1--)
	    for (d2 = Distance - 1; d2 >= 0; d2--)   
            {
/*----------------------------------------------------------------------------*/
              if ( (RunType != 4 ||
                    (RunType == 4 && total_motifs <= MOTAUTO4) ) &&
	          (*aa3_obs)[aa2][aa3][d1][d2] >= (unsigned char) Signif )
              {
                /* Record all information about motif in motif array: */
	        motif[total_motifs].aa1 = aa1;
	        motif[total_motifs].aa2 = aa2;
	        motif[total_motifs].aa3 = aa3;
	        motif[total_motifs].distance1 = d1 + 1;
	        motif[total_motifs].distance2 = d2 + 1;
		motif[total_motifs].freq = motif[total_motifs].dups = 0;
	        motif[total_motifs].group = motif[total_motifs].sub_group = 0;
                motif[total_motifs].mots = 1;
		count = get_count(&motif[total_motifs]);

      /* Now freq is the total frequency of the motif and count is the number
	 of different sequences the motif occurred in.  Therefore the number
	 of duplicates is (freq - count).  The next step is to make sure
	 motif occurs Dups or less in each sequence and that freq is over the
	 significance level.  If so, store the motif in the motif array: */

              if (motif[total_motifs].dups <= Dups && count >= Signif)   
              {
	        duplicate = NO;
		if (total_motifs)
		   duplicate = check_dup(&motif[total_motifs],
                                         &motif[total_motifs-1]);

              /* If it's not an obvious duplicate, add it to the motif list */
	      if (!duplicate)   
              {
	        domain_width = d1 + d2 + 3;		/* width of motif */
	        if (domain_width < MIN_DOMAIN_WIDTH)
	          domain_width = MIN_DOMAIN_WIDTH;
	        if (domain_width > MAX_DOMAIN_WIDTH)
	          domain_width = MAX_DOMAIN_WIDTH;
	        motif[total_motifs].domain = domain_width;

		score_mot(&motif[total_motifs]);

	        /*  Display motif */
	        pr_num_to_aa(aa1); pr_num_to_aa(aa2);
	        pr_num_to_aa(aa3); printf(",");
	        total_motifs++;

	        /* Flooded with too many motifs, sort and cut back to 
	                 RELEVANT_MOTIFS */
	        if (total_motifs == MAX_MOTIFS)
                {
	             save_motifs = total_motifs;
	             total_motifs = cut_motifs(motif, save_motifs);
                }
/*-------------------------------------------------------------------------*/

	  }  /* end of if !duplicate */
	}  /* end of if freq - count <= Dups  */
      }  /* end of if >= Signif */
      (*aa3_obs)[aa2][aa3][d1][d2] = 0;       /* Reset scoring matrix to 0 */

    } /* end of scan array */
   } /* >>>>  end of RunType check */
  }  /* End of aa1 loop */
/*------------------------------------------------------------------------*/
/*               Sort data and weed out redundant information.            */
/*------------------------------------------------------------------------*/
  printf("\n\n");

  /* Cut back to RELEVANT_MOTIFS, sorted by frequency and score */
  save_motifs = total_motifs;                                 /*JGH*/
  if (total_motifs > 0)
  {
      if (RunType != 4) total_motifs = cut_motifs(motif, save_motifs);
      else if (total_motifs > RELEVANT_MOTIFS)
		total_motifs = RELEVANT_MOTIFS;
  }
  real_total_motifs = total_motifs;

  printf("\nCPU time to find and sort motifs: %ld seconds.",
						time(0) - CPU_time);
/*------------------------------------------------------------------------*/
/*                      Consider results & decide what to do next         */
/*------------------------------------------------------------------------*/
  found = 0;
  if (real_total_motifs > 0)
	found = Signif;               /* signif level with motifs found JGH*/

  /* User has option of saving only top N scores:  */
  printf("\nFound %d motifs.", real_total_motifs);  
  total_motifs = real_total_motifs;
  if (RunType==0 || (RunType==2 && real_total_motifs >= MOTAUTO3))     /*4JGH*/
  {
     printf(" Enter number to save or -1 to redo [1-%d; %d]: ",
	      real_total_motifs, real_total_motifs);
     gets(chsig); total_motifs = atoi(chsig);
     if (total_motifs == 0 || total_motifs > real_total_motifs)
	 total_motifs = real_total_motifs;
     else if (total_motifs < 0) 
     {
       printf("\nPrevious significance level was %d.", Signif);
       printf(" Enter new level [2-%d; %d]: ", NumSeqs, Signif); 
       prevsig = Signif; gets(intemp);  Signif = atoi(intemp);  
       if (Signif < 2 || Signif > NumSeqs) Signif = prevsig;  
       printf("Previous number of internal duplications was %d.", Dups); 
       printf(" Enter new number [0-;%d]: ", Dups);  
       prevdup = Dups; gets(intemp);  Dups = atoi(intemp);  
       if (Dups < 0 || !strlen(intemp) ) Dups = prevdup;  
       printf(" Using previous distance  %d\n", Distance);
       goto REDO;
     }
  }

/*------ If nothing was found, try again automatically if specified     JGH*/
    if ( (RunType == 2 || RunType == 3 || RunType == 5) &&
	  real_total_motifs < MOTAUTO3 && Signif > MINSEQS)          /*6JGH*/
    {
       Signif--;
       if (Signif >= MINSEQS)  goto REDO;
    }
    if (RunType == 4)
    {
       prevsig = Signif;
       if (real_total_motifs > MOTAUTO4)      /* found too many */
       {
	  Signif++;
	  if (Signif < NumSeqs) goto REDO;
	  else                  goto EXIT;  /*too many motifs at all signif*/
       }
       else if (real_total_motifs <= MOTAUTO4 && found <= 0)
       {                              /* haven't found enough motifs yet  */
	  Signif = (int) ( (prevsig+1)/2 );
	  if (Signif == prevsig) Signif--;
	  if (Signif >= minseqs) goto REDO;
	  else                        /* no motifs at any signif */
          {  Signif = prevsig;
             goto EXIT;
          }
       }
       else if (real_total_motifs <= MOTAUTO4 && found > 0)
       {                   /*  some found, but not too many at this level */
	  goto EXIT;
       }
    }
    if (total_motifs < 1) goto EXIT;
/*NOTE:  Should be able to free aa3_obs now, but causes problems!   */

/*------------------------------------------------------------------------*/
/*                          Print out motifs.                             */
/*   NOTE:  print_motif() modifies the motif[] structures!             JGH*/
/*------------------------------------------------------------------------*/
/*   Default output is always short form.                              JGH*/
   Short_Form = YES; First = NO;
   if (RunType==0 || RunType==2)
   {
     printf("Enter 'g' to print first motif in each group: ");
     gets(chsig);
     if (chsig[0] == 'g' || chsig[0] == 'G') First = YES;
     printf("Enter 'l' for long form:  "); 
     gets(chsig);
     if (chsig[0] == 'l' || chsig[0] == 'L') Short_Form = NO;
  }
/*-------------------------------------------------------------------------*/
  /* Separate groups of related proteins with horizontal line: */
  /* Print out all motifs: */
  /*  Sort motifs by group first for printing.   JGH*/
  qsort(motif, total_motifs, sizeof(struct motif_struct), compare_motgrps);
  max_group = 0;
  for (j = 0; j < total_motifs; j++)
  {
    if (motif[j].freq > 0 && motif[j].group > max_group) {
      group = motif[j].group;
      max_group = group;
      printf("\n_________________________________Group #%d_________________________________\n\n", group);
      print_motif(&motif[j]);

      total_score += (unsigned long) motif[j].score;

      for (k = j + 1; k < total_motifs && motif[j].freq > 0; k++)
	if (motif[k].freq > 0 && motif[k].group == group)
	  {
	      if (!Short_Form && !First) printf("\n");
	      if (!First) print_motif(&motif[k]);
	      total_score += (unsigned long) motif[k].score;
	    }
    }
  }  /* end of for j */

  printf("\n___________________________________________________________________________\n\n");
  printf("Printed %d motifs out of %d total.\nMean score: %ld \n",
	total_motifs, real_total_motifs,
	(unsigned long) total_score/total_motifs);

/*-----------------------------------------------------------------------*/
/*                   Write the motifs out to a file                   JGH*/
/*-----------------------------------------------------------------------*/
  write_motifs(total_motifs, motif);

/*------------------------------------------------------------------------*/
/*                       Print out motif maps.                            */
/*------------------------------------------------------------------------*/
  if (max_group > 9) max_group = 9;
  n = 0;
  if (RunType == 0 || RunType == 2)
  {
     printf("Enter number of motif groups to map (maximum %d)", max_group);
     printf(" or 0 to assemble blocks: ");
     gets(chsig); n = atoi(chsig);
  }
  else  n = 0; 
  if (n <= 0) goto EXIT;
  if (n > max_group) n = max_group;
  max_group = n;
  map_seqs(max_group, total_motifs, motif);

/*------------------------------------------------------------------------*/
/*     Finish up                                                          */
/*------------------------------------------------------------------------*/
EXIT:
   /* Free up all memory and quit. */
   free(motif);
   for(i = 0; i < NumSeqs; free(Seq[i++]));                     /*JGH*/
/*   for (i = 0; i < 20; i++)       **causes program to hang**
    for (j = 0; j < 20; j++)
      for (n = 0; n < Distance; n++)
	free ( (*aa3_obs)[i][j][n] );
*/
   free(aa3_obs);

/*---- Get maximum normalized block score for blocks with the motif
       in Signif sequences -------*/
   if (RunType == 4)
   {
      n = 0;      /* Get maximum normalized block score */
      for (j=0; j<real_total_motifs; j++)
      {
	 if ((motif[j].freq-motif[j].dups)==Signif)
	 {
	    k = (int) motif[j].score;
	    if ( k > n)  n = k;
	 }
      }
      printf("\nStarting motifj...\n");
      kr_itoa(Signif, chsig, 10);
      kr_itoa(save_dups, chdup, 10);
      kr_itoa(Distance, chdis, 10);
      kr_itoa(n, chdrop, 10);
      fflush(stdout);
      if (batch == NO)
      {  strcpy(intemp, "-"); strcat(intemp, Batch_Filename); }
      else strcpy(intemp, Batch_Filename);
      execlp(argv[0], argv[0], "3", intemp,
	     chsig, chdup, chdis, chdrop, NULL);
   }
   else if (total_motifs > 0)
   {
      printf("\nStarting motomat...\n");
      fflush(stdout);
      execlp("motomat", "motomat", Mot_Filename, "1", NULL);
   }
   else  printf("\n");
   exit(0);
}  /*  end of main */
/*------------------------------------------------------------------------*/
/*                       END OF MAIN PROCEDURE                            */
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/*                   Function to print out motifs.                        */
/*   Also calls re-alignment routine, this is because some of the         */
/*   information from the re-alignment is printed, but not saved          */
/*------------------------------------------------------------------------*/
void print_motif(pat)
struct motif_struct *pat;
{
int i, x, aa_freq[20], aa, firstaa;
int offset, p_index, conserved;
int motseq, seq, pos, count;
char prototype[MAX_DOMAIN_WIDTH];
struct realign_struct realign_data;

  /* Print out the three aa sequence and the score: */
  printf("\nMotif: ");
  pr_num_to_aa(pat->aa1);
  pr_num_to_aa(pat->aa2);
  pr_num_to_aa(pat->aa3);
  printf(" in %d proteins (%d duplicates), Merged %d, Score %d\n", 
         pat->freq - pat->dups, pat->dups, pat->mots, pat->score);
  if (!Short_Form) printf("\n");

  count = pat->freq - pat->dups;
  offset = (int) (pat->domain - pat->distance1 - pat->distance2 - 1) / 2;

  /* Calculate frequency of each aa in each column of the block
     and print the conservation line heading */
  for (i = 0; i < 20; i++) aa_freq[i] = 0;
  for (i = 0; i < pat->domain; i++)
  {
    /* See if column is composed of a single amino acid (if so, c = 1) or
       is composed of more than 50% of a single amino acid (if so, record
       in prototype array for 50% degenerate pattern): */
    conserved = YES;
    firstaa = -1;
    for (x = 0; x < pat->freq; x++)
    {
      pos = pat->pos[x] + i - offset;
      if (pos >= 0 && pos < Len[ pat->seq_no[x] ])
      {
         aa = Seq[ pat->seq_no[x] ][pos];
         aa_freq[aa]++;
         if (x == 0) firstaa = aa;
         else if (aa != firstaa) conserved = NO;
      }
      else conserved = NO;	/* Any null in column => not conserved */
    }
    prototype[i] = -1;
    for (x = 0; x < 20; x++)
    {
      if (aa_freq[x] > count / 2) prototype[i] = x;
      aa_freq[x] = 0;
    }

    /* Print out conservation evaluation line (line over block composed
       of +, * and amino acid symbols): */
    p_index = pat->scores[i];
    if (conserved)
       pr_num_to_aa_space(Seq[pat->seq_no[0]][pat->pos[0] + i - offset]);
      else if(p_index < SMatrix->highpass) printf("  ");
      else if(p_index < 1000) printf("+ ");
      else printf("* ");
  }


  /* Realign every sequence against block and print it out: */
  realign(pat, &realign_data);

  if (!Short_Form) printf("  Score (Position) for best two fits.\n");
  for (x = 0; x < NumSeqs; x++) 
  {
      seq = realign_data.sequence[x];
      if (!Short_Form && x == count) 
      {
	for (i = 0; i < pat->domain - 1; i++) printf("--");
	printf("\n");
      }

      /* Print out sequence, using '.' when off one end: */
      for (i = 0; i < pat->domain; i++)
      {
        pos = realign_data.max_pos[x] + i - offset;
	if (!Short_Form && pos >= 0 && pos < Len[seq] )
	  pr_num_to_aa_space(Seq[seq][pos]);
	else if (!Short_Form) printf(". ");
      }

      /* Report scores and positions */
      if (!Short_Form)
      {
	 printf("  %3ld (%3d)  %3ld (%3d)  %s",
		realign_data.max_score[x]/10000,
                realign_data.max_pos[x] + 1,
		realign_data.max2_score[x]/10000,
                realign_data.max2_pos[x] + 1, 
                Seqname[seq]);
         /*----  Indicate seqs with motif that were realigned ------*/
         /*  Find first corresponding pat-> sequences, could be > 1 if
             mot->dups is > 0 */
         i = 0;
         motseq = -1;
         while (motseq < 0 && i < pat->freq)
         {
            if (pat->seq_no[i] == seq)
            {
                motseq = seq;
                if (pat->pos[i] != realign_data.motif_pos[x])
	           printf(" (Motif at %d)", realign_data.motif_pos[x] + 1); 
            }
            i++;
         }
         printf("\n");
      }
  }  /*  end of sequence  x */

  /* Print prototype (50% degenerate pattern): */
  if (!Short_Form)
  {
       printf("\n");
       for (i = 0; i < pat->domain; i++)
	 pr_num_to_aa_space(prototype[i]);
       printf("  Most common amino acids\n\n");
  }
}     /*  end of print_motif */

/*------------------------------------------------------------------------*/
/*                 Function to calculate un-weighted pam index            */
/*------------------------------------------------------------------------*/

/* Performs pairwise comparison of all amino acids in a column (col) of
   length len.  Reports average score multiplied by 100 to retain two
   place accuracy while keeping scores to integers for speed: */
int pamidx(len, col)
int len, col[MAXFREQ];
{
int *col_i, *col_j;
int val = 0, real_len = len, *end_of_col;

  end_of_col = col + len;

  /* Go through all two letter combinations.
     Note: -1 marks an empty column (which is ignored in calculation). */
  for (col_i = col; col_i < end_of_col - 1; col_i++)
    if (*col_i != -1) {
      for (col_j = col_i + 1; col_j < end_of_col; col_j++)
	if (*col_j != -1) val += SMatrix->scores[*col_i][*col_j];
      }
     else real_len--;

  if (*(end_of_col - 1) == -1) real_len--;
  if (real_len > 1)
		return((int) ((long int) val * 100 / (real_len * (real_len - 1) / 2)));
  else return(0);
} /* end of pamidx */
/*------------------------------------------------------------------------*/
/*                 Function to calculate weighted pam index               */
/*------------------------------------------------------------------------*/

/* Performs pairwise comparison of all amino acids in a column (col) of
   length len.  Reports average score multiplied by 100 to retain two
   place accuracy while keeping scores to integers for speed: */
int wpamidx(len, col, wts)
int len, col[MAXFREQ];
double wts[MAXFREQ];
{
  int seq_i, seq_j, nval;
  double val;

  nval = 0;
  val = 0.0;
  /* Go through all two letter combinations.
     Note: -1 marks an empty column (which is ignored in calculation). */
  for (seq_i = 0; seq_i < len; seq_i++)
    if (col[seq_i] != -1)
    {
      for (seq_j = seq_i + 1; seq_j < len; seq_j++)
	if (col[seq_j] != -1)
        {
           val += (wts[seq_i]*wts[seq_j]*
                (double)SMatrix->scores[col[seq_i]][col[seq_j]]);
           nval++;
        }
     }

   if (nval > 0)
   {
      val *= 100.0;
      val /= (double) nval;
      return( round(val) );
   }
   else return(0);

}  /* end of wpamidx */
/*======================================================================
   Remove motifs that are slightly offset from a more significant
     motif or a higher scoring motif, in that order of priority.
     (More significant means more sequences have the motif).
   Lower motif is compared with higher & dropped if the starting positions
     of the two motifs in all sequences in the lower are the same 
     distance apart & that distance is <= 10.
     If this is true for at least 4 sequences, then the lower is put in
     the higher's group.
   Modified alg to eliminate more low-scoring motifs:
     Lower motif is discarded if it the same distance apart as indicated
     above indicates any overlap between the two motifs.
=======================================================================*/
int cut_motifs(motif, total_motifs)
struct motif_struct *motif;
int total_motifs;
{
  int i, j, k, x, real_total_motifs, dist, matches, group, sub_group;
  int seq_list[MAXSEQS], countk, widthk, widthj, mincount, maxcount;
  int firstk, lastk, maxkeep;

/*---- First, descending sort of motifs by freq-dups & score -------------*/
  printf(":\ncutting %d", total_motifs);
  qsort(motif, total_motifs, sizeof(struct motif_struct), compare_scores);

/*-----Clear out groups & get ready to re-group ---------*/
  for (j=0; j<total_motifs; j++)
  {
      motif[j].group = motif[j].sub_group = 0;
  }

  real_total_motifs = total_motifs;
  group = sub_group = 1;
  for (j = 0; j < total_motifs; j++)
  {
/*-----Second, if normalized score is < Drop, drop motif------------------*/
    if (motif[j].freq > 0 && (int) motif[j].score  < Drop)
    {  motif[j].freq = 0; real_total_motifs--; }

/*-----Third, if two motifs overlap, drop the
    lower scoring one by setting its freq to 0 ---------------------------*/
    if (motif[j].freq > 0)
    {
       /*   Set the group number  */
       if (motif[j].group < 1)
       {
          motif[j].group = group++;
          motif[j].sub_group = 1;
          sub_group = 2;
       }

       /*  Compare motif[j] with all motifs below it   */
       /* motif[k] is dropped if all of its sequences are in [j],
          and all overlap. The motifs don't have to be the
          SAME distance apart in all sequences ... */
       widthj = (int) motif[j].distance1 + (int) motif[j].distance2 + 3;
       for (k = j + 1; k < total_motifs && motif[j].freq > 0; k++)
         if (motif[k].freq > 0)   /*  may have already been dropped */
	 {
	     matches = 0;
             /*  seq_list[i] == # of times sequence i is in motif[k] */
             for (x = 0; x < NumSeqs; x++) seq_list[x] = 0;
             widthk = (int) motif[k].distance1 + (int) motif[k].distance2 + 3;
             countk = motif[k].freq - motif[k].dups;
	     for (x = 0; x < motif[k].freq; x++)
	     {
	        /* Make sure sequence numbers match before comparing
		   positions, if a sequence is represented by dups
                   will only use the first position */
                seq_list[ motif[k].seq_no[x] ]++;
                if (seq_list[ motif[k].seq_no[x] ] == 1)
                {
                   i = 0;
                   while (i < motif[j].freq &&
		          motif[k].seq_no[x] != motif[j].seq_no[i]) 
                              i++;
                   if (i < motif[j].freq )
                   {
	              /* Do we overlap?  */
                      dist = motif[j].pos[i] - motif[k].pos[x];
                      /*  dist >= 0:  kkkkkkkkkkkk
                                                 jjjjjjjjjjjjjjjj
                          dist <= 0:  jjjjjjjjjjjj
                                                 kkkkkkkkkkkkkk         */
                      if ((dist >= 0 && dist <= widthk)  ||
                          (dist <= 0 && (0-dist) <= widthj) )
	                        matches++;
                   }
                 }
	      }
	      /* Mark redundant motifs with a zero freq so that during
	       the next sort, those motifs will drop to bottom. */
	      if (matches == countk )
              {
	         motif[k].freq = 0;
	         motif[j].mots = motif[j].mots + 1;
	         /* Dropping motif out, so keep track of how many are left: */
	         real_total_motifs--;
	      }
              /*  Same group if at least half the common sequences match */
	      else if (matches >= round( (float) countk / 2.0) &&
                       motif[k].group < 1)
              {
	         motif[k].group = motif[j].group;
	         motif[k].sub_group = sub_group++;
	      }
	  }  /* end of if motif[k].freq > 0 */
    }  /* end of if motif[j].freq > 0 */
  }  /*  end of for j */
  printf(" %d", real_total_motifs);

/*-------- Fourth, if there are more than RELEVANT_MOTIFS left, then
           save no more than #mots/#groups in each group  ------------*/
  /*   sort by group, sub_group and freq-dups , but group 0 are first */
  if (real_total_motifs > RELEVANT_MOTIFS)
  {
/*NOTE>>>  Might be better to use compare_motgrps() instead ...  */
  qsort(motif, total_motifs, sizeof(struct motif_struct), compare_subgrps);
  maxkeep = round( (double) real_total_motifs / group );
  if (maxkeep < 1 ) maxkeep = 1;
  firstk = 0;
  while (firstk < total_motifs && motif[firstk].group == 0) firstk++;
  for (j=1; j < group; j++)
  {
      maxcount = 0;
      mincount = NumSeqs;
      k = firstk;
      while (k < total_motifs && motif[k].group == j)
      {
         if (motif[k].freq > 0)
         {
            if (motif[k].sub_group > maxkeep)
            {
                motif[k].freq = 0;
	        real_total_motifs--;
            }
/*
            countk = motif[k].freq - motif[k].dups;
            if (countk < mincount) mincount = countk;
            if (countk > maxcount) maxcount = countk;
*/
         }
         k++;
      }
      lastk = k - 1;
/*
      if (mincount < maxcount)
      {
         for (k = lastk; k > firstk; k--)
         {
            if ((motif[k].freq - motif[k].dups) == mincount)
            {
                motif[k].freq = 0;
	        real_total_motifs--;
            }
         }
      }
*/
      firstk = lastk + 1;
  } /* end of group j */
  }   /* end of RELEVANT_MOTIF check */
  printf(" %d", real_total_motifs);

/*-------- Fifth, sort again by freq-dups and score so 
           zero frequency motifs drop out ----------------------------*/
  /*---NOTE: This sort can scramble group order! ------*/
  qsort(motif, total_motifs, sizeof(struct motif_struct), compare_scores);

  if (real_total_motifs > RELEVANT_MOTIFS)
      real_total_motifs = RELEVANT_MOTIFS;
  if (real_total_motifs < 0) real_total_motifs = 0;
  printf(" %d:", real_total_motifs);
  printf(" last score = %d, groups = %d\n",
          motif[real_total_motifs - 1].score, group);

   return(real_total_motifs);
}  /* end of cut_motifs */
/*=======================================================================*/
/*    Map motifs in all sequences                                        */
/*=======================================================================*/
void map_seqs(max_group, total_motifs, motif)
int max_group, total_motifs;
struct motif_struct *motif;
{
  struct group_struct *groups;
  int i, j, x, group_no;

  groups = (struct group_struct *)
	      malloc(sizeof(struct group_struct) * total_motifs);
  if (groups == NULL) {
	printf("Couldn't allocate GROUPS.\n");
	exit(-1);
	}

  for (x = 0; x < NumSeqs; x++)  {
    group_no = 0;
    for (i = 0; i < total_motifs; i++)
      for (j = 0; j < motif[i].freq; j++)
	if (motif[i].seq_no[j] == x) {
	  groups[group_no].position = motif[i].pos[j];
	  groups[group_no].group_no = motif[i].group;
	  groups[group_no++].sub_no = motif[i].sub_group;
	  }
    if (group_no > 0) {
      qsort(groups, group_no, sizeof(struct group_struct), compare_groups);
      i = group_no;
      for (group_no = 0, j = 1; j < i; j++)
	if (groups[group_no].group_no == groups[j].group_no &&
	    groups[j].position - groups[group_no].position < 12)
	  groups[group_no].sub_no = 0;
	else { groups[++group_no].position = groups[j].position;
	       groups[group_no].group_no = groups[j].group_no;
	       groups[group_no].sub_no = groups[j].sub_no;}
      group_no++;
      }

    printf("                    Protein %s:\n", Seqname[x]);

    j = -1;
    for (i = 0; i < group_no; i++)  {
      for (j++; j < groups[i].position / 6; j++) printf(".");
      if (groups[i].group_no <= max_group) printf("%d",groups[i].group_no);
	else printf(".");
      if (groups[i].sub_no > 0) printf("%c",(char) groups[i].sub_no + 64);
      }
    for (j++; j < Len[x] / 6; j++) printf(".");
    printf(" %d\n\n",Len[x]);
    }
   free(groups);
}  /*  end of map_seqs */

/*------------------------------------------------------------------------*/
/*  Routines to convert amino acid symbols to numbers and vice versa.     */
/*    These routines can be found in "motmisc.c"                          */
/*------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*              Routines used by the sorting procedure                   */
/*-----------------------------------------------------------------------*/

/* Function to compare two motifs by frequency and score, in that
   order of priority.  Used to sort motif array in descending order
   by frequency, number of "merged" motifs and score: */
int compare_scores(motif_1, motif_2)
struct motif_struct *motif_1, *motif_2;
{
int sequences;

  /* Frequencies without the internal duplicates */
  sequences = (motif_2->freq - motif_2->dups) -
	      (motif_1->freq - motif_1->dups);

  /* If frequencies are identical, sort by #of motifs: */
  if (!sequences) sequences = (motif_2->mots - motif_1->mots);
   else return(sequences);
  /* If frequencies and motifs identical, sort by score */
  if (!sequences) return(motif_2->score - motif_1->score);
   else return(sequences);
}
/*====================================================================*/
/* Compare two motifs by group, count = freq-dups, #merged motifs,
    and score */
int compare_motgrps(motif_1, motif_2)
struct motif_struct *motif_1, *motif_2;
{
   int diff;

   diff = (motif_1->group - motif_2->group);
   if (diff) return(diff);
   else
   {
      diff = (motif_2->freq - motif_2->dups) -
             (motif_1->freq - motif_1->dups);
      if (diff) return(diff);
      else
      {
         diff = (motif_2->mots - motif_1->mots);
         if (diff) return(diff);
         else return(motif_2->score - motif_1->score);
      }
   }
}  /*  end of compare_motgrps()  */
/*====================================================================*/
/* Compare two motifs by group and sub_group   */
int compare_subgrps(motif_1, motif_2)
struct motif_struct *motif_1, *motif_2;
{
   int diff;

   diff = (motif_1->group - motif_2->group);
   if (diff) return(diff);
   else
   {
      diff = (motif_1->sub_group - motif_2->sub_group);
      if (diff) return(diff);
      else
      {
         diff = (motif_2->freq - motif_2->dups) -
                (motif_1->freq - motif_1->dups);
         return(diff);
      }
   }
}  /*  end of compare_subgrps()  */
 
/*======================================================================*/
/* Function to compare two groups by position.  Used to rank groups
   in ascending order by position for map_seqs()                        */
int compare_groups(group_1, group_2)
struct group_struct *group_1, *group_2;
{
  return(group_1->position - group_2->position);
}  /* end of compare_groups */

/*=======================================================================
     Opens a sequence file, reads line[], updates filepath[] and
     seqname[]
=======================================================================*/
FILE *getfile(line, filepath, seqname)
char *line, *filepath, *seqname;
{
   FILE *fp;
   char *ptr, ctemp[FNAMELEN];
   
   fp = NULL;
   seqname[0] = '\0';
   /*     Check for new path  */
   if (strstr(line, "/") != NULL || strstr(line, ":") != NULL)
   {
	ptr = strtok(line, " \n\t\r");
	strcpy(filepath, ptr);
   }
   /*    If not a title, assume it's a file name */
   else if (line[0] != '>')
   {
	ptr = strtok(line, " \n\t\r");
	if (ptr != NULL)
	{
		if ((int) strlen(ptr) > SNAMELEN - 1) ptr[SNAMELEN - 1] = '\0';
		strcpy(seqname, ptr);
		strcpy(ctemp, filepath);
		strcat(ctemp, seqname);
		strcat(ctemp, PROTEIN_EXTENSION);
		if ( (fp = fopen(ctemp, READ) ) == NULL)
			printf("***%s not found\n", ctemp);
	}
   }

   return(fp);
}  /*  end of getfile  */

/*======================================================================
       Read all sequences in file, fasta format assumed
	Updates global variables NumSeqs, Seqname[], Seq[], Len[]
========================================================================*/
int getseqs(fp)
FILE *fp;
{
   int n, ns, itemp, i, j;
   char line[MAXLINE], *ptr;

   n = 0;
   ns = NumSeqs - 1;
   
   while (!feof(fp) && fgets(line, MAXLINE, fp) != NULL)
   {
	if (line[0] == '>')
	{
	   n++; ns++;
	   ptr = strtok(line+1, " \t\r\n");
	   if ((int) strlen(ptr) > SNAMELEN - 1) ptr[SNAMELEN - 1] = '\0';
	   strcpy(Seqname[ns], ptr);
	   Seq[ns] = (int *) malloc(MAX_LENGTH * sizeof(int));
	   if (Seq[ns] == NULL)
	   {
		printf("Couldn't allocate Seq[%d]\n", ns);
		exit(-1);
	   }
	   Len[ns] = j = 0;
        }
        else if (Seq[ns] != NULL && Len[ns] < MAX_LENGTH + (int) strlen(line) )
        {
	   for (i=0; i < (int) strlen(line); i++)
           {
                itemp = aachar_to_num(line[i]);
		if (itemp >= 0 && itemp <= 20)
                {
  	  	   Seq[ns][j++] = itemp;
		   Len[ns] += 1;
                }
           }
	}
   }

   NumSeqs = ns + 1;
   return(n);
}  /*  end of getseqs */
/*=====================================================================
	See if sequences should be shuffled, and do it
========================================================================*/
int shuffle()
{
   int shuff, i, j, k, n;
   char temp, shuffle[3];

   shuff = NO;
   if (RunType == 0 || RunType == 2)
   {
     printf("\n");
     printf("Use shuffled sequences [y/n; n]? ");
     gets(shuffle);
   }
   else if (RunType == 4 || RunType == 5)
     strcpy(shuffle, "y");
   else  
     strcpy(shuffle,"n");

   /*---- For RunTypes 4 or 5, sequences are always shuffled the same way.
       This is assured by reseting the generator with the same seed ---- */
   if (shuffle[0] == 'y' || shuffle[0] == 'Y')
   {
      shuff = YES;
      randomize();
      for (n = 0; n < NumSeqs; n++)
      {
          j = Len[n];
          for(i = 0; i < j; i++)
          {
             if (RunType == 4 || RunType == 5)
                 srand(i);                      /* i is the amino acid index */
             k = rand() % j;
             temp = Seq[n][i];
             Seq[n][i] = Seq[n][k];
             Seq[n][k] = temp;
          }
      }
   }

   return(shuff);
}  /* end of shuffle */
/*=====================================================================
	Print some basic statistics
	Expects global variables NumSeqs & Seq[][]
========================================================================*/
void print_stats()
{
   unsigned int aa_obs[20], total_aa;
   int n, i;

   for (i=0; i < 20; i++) aa_obs[i] = 0;
   total_aa = 0;

  for (n = 0; n < NumSeqs; n++)  {
    total_aa += Len[n];			/* Compute total number of aa's in */
    for (i = 0; i < Len[n]; i++)	/* all sequences and individual    */
      aa_obs[Seq[n][i]]++;		/* frequencies of the aa's.	   */
    }

  for (i = 0; i < 20; i++) {
    printf("  ");
    pr_num_to_aa(i);			/* Print amino acid's letter code */
    printf("   ");
    }
  printf("\n");

  /* Print aa frequency table: */
  for (i = 0; i < 20; i++) printf("%5d ", aa_obs[i]);
  printf("\n");

  /* Print frequencies as fraction of total number of aa's: */
  for (i = 0; i < 20; i++)  printf("%5.3f ", (float) aa_obs[i]/total_aa);
  printf("\n");

  printf("Total amino acids in %d proteins: %d\n\n", NumSeqs, total_aa);
}  /* end of print_stats */
/*======================================================================*/
   /* Found a three amino acid motif that occurs at least Signif
      number of times.  Scan seqs for the motif 
      and see how many sequences motif occurred in (variable count): */
int get_count(motif)
struct motif_struct *motif;
{
   int d1, d2, count, n, i, flag;

   count = 0;
   motif->freq = 0;
   for (i=0; i < MAXFREQ; i++)
   {  motif->seq_no[i] = motif->pos[i] = 0; }
   d1 = motif->distance1 - 1;
   d2 = motif->distance2 - 1;
   for (n = 0; n < NumSeqs; n++)   
   {
      flag = 0;
      for (i = 0; i < Len[n]; i++)
            if (Seq[n][i]==motif->aa1 && Seq[n][i+1+d1]==motif->aa2 &&
                                      Seq[n][i+2+d1+d2]==motif->aa3 &&
                                       motif->freq < MAXFREQ)
            {
                motif->seq_no[motif->freq] = n;
                motif->pos[motif->freq++] = i;
                flag++;
            }

            if (flag) count++;
    }
    motif->dups = motif->freq - count;
    return(count);
}  /*  end of get_count */
/*==========================================================================*/
	/* Weed out sequential duplicates that start with the same first
	   amino acid in the same position (eliminates alternative
	   combinations of three aa's from motif of four aa's or higher): */
int check_dup(motif1, motif2)
struct motif_struct *motif1, *motif2;
{
   int duplicate, x;

   duplicate = YES;
   for (x = 0; x < motif1->freq && duplicate; x++)
      duplicate = (motif1->pos[x] == motif2->pos[x] &&
	           motif1->seq_no[x] == motif2->seq_no[x]);
   return(duplicate);
}  /* end of check_dup */
/*=======================================================================
	Makes an aa column vector for pamidx
========================================================================*/
void score_mot(motif)
struct motif_struct *motif;
{
   int offset, col, pos, seq, iseq, p_index, aa[MAXFREQ];
   double weights[MAXFREQ];

   /*   compute sequence weights and nident */
   pb_weights(motif, weights);

   offset = (int) (motif->domain - motif->distance1 - motif->distance2 - 1) / 2;
   motif->score = 0;
   for (col = 0; col < motif->domain; col++) 
   {
      for (iseq = 0; iseq < motif->freq; iseq++)
      {
         pos = motif->pos[iseq] + col - offset;
         seq = motif->seq_no[iseq];
         if (pos >= 0 && pos < Len[seq]) aa[iseq] = Seq[seq][pos]; 
         else                            aa[iseq] = -1;
      }
      motif->scores[col] = (p_index=wpamidx(motif->freq, aa, weights));
/*
      motif->scores[col] = (p_index=pamidx(motif->freq, aa));
*/

      /* Only add in high pam scores and convert back to
	            standard pam scores (by subtracting HighPass): */
      if (p_index >= SMatrix->highpass)
         motif->score += (p_index - SMatrix->highpass);
   }  /*  end of col */

   motif->score = motif->score / sqrt((float) motif->domain);

}  /*  end of score_mot */
/*======================================================================
	Put position-based sequence weights in the weights array
	col = block column, pos = sequence position
	Version for motifj:
======================================================================*/
void pb_weights(block, weights)
struct motif_struct *block;
double weights[MAXFREQ];
{
   struct pb_counts *pb;
   double factor, dtemp;
   int iseq, seq, pos, col, aa, width, offset, nident;
   char c1, c[2];

   offset = (int) (block->domain - block->distance1 - block->distance2 - 1) / 2;
   width = block->domain;
   nident = 0;

   pb = (struct pb_counts *) malloc(width*sizeof(struct pb_counts));

   for (col = 0; col < width; col++)
   { 
      pb[col].diffaas = 0.0;
      for (aa = 0; aa < 20; aa++)
         pb[col].naas[aa] = (double) 0.0;
   }
  
   for (col = 0; col < width; col++)
   {
      for (iseq = 0; iseq < block->freq; iseq++)
      {
         seq = block->seq_no[iseq];
	 pos = block->pos[iseq] + col - offset;
	 if (pos >= 0 && pos < Len[seq] )
	 {
	     aa = Seq[seq][pos];
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
   }

   for (col = 0; col < width; col++)
   { 
      for (aa = 0; aa < 20; aa++)
      {
         if (pb[col].naas[aa] > 0.0)
         {
            pb[col].diffaas += 1;	/* # of different types of aas in col */
         }
	 if (pb[col].naas[aa] >= (double) Signif)
	 {
	     nident += 1; 	/* # of conserved cols in block */
	 }
      }
   }

   for (iseq = 0; iseq < block->freq; iseq++)
   {
      seq = block->seq_no[iseq];
      weights[iseq] = 0.0;
      for (col = 0; col < width; col++)
      {
	 pos = block->pos[iseq] + col - offset;
	 if (pos >= 0 && pos < Len[seq] )
	 {
	     aa = Seq[seq][pos];
             dtemp = pb[col].diffaas * pb[col].naas[aa];
             if (dtemp > 0.0) weights[iseq] += 1.0 / dtemp;
         }
      } 
   }
 
   /*  Normalize weights to add to #seqs == block->freq  */
   dtemp = 0.0;
   for (iseq = 0; iseq < block->freq; iseq++) dtemp += weights[iseq];
   factor = 1.0;
   if (dtemp > 0.0) factor = (double) block->freq/dtemp;
   for (iseq = 0; iseq < block->freq; iseq++) weights[iseq] *= factor;

   free(pb);
}  /* end of pb_weights */
/*===================================================================
    Realign sequences in the motif & return some information.
    Assigns positions to sequences not having the motif,
    may assign new positions to sequences with motif.
    Probably should re-compute column and block scores, but doesn't
    (motomat does).
=====================================================================*/
void realign(mot, dat)
struct motif_struct *mot;
struct realign_struct *dat;
{
   int x, i, j, k, count, offset, pos, position, datseq, motseq;
   int iscore, aa1, aa2, aas[20][MAX_DOMAIN_WIDTH];
   long int score;

   offset = (int) (mot->domain - mot->distance1 - mot->distance2 - 1) / 2;
   count = mot->freq - mot->dups;

   /*  Initialize  */
   for (x=0; x < NumSeqs; x++)
   {  
     dat->sequence[x] = 0;
     dat->max_score[x] = dat->max2_score[x] = (long int) 0;
     dat->max_pos[x] = dat->max2_pos[x] = 0;
     dat->motif_pos[x] = 0;
   }

  /*-----------------------------------------------------------------------*/
  /*    Order sequences with motif before those without */
  /*   seq_no[] may appear more than once in mot->, but only once in dat->
     EG: suppose NumSeqs = 5, mot->dups = 1, mot->freq = 5, count = 4
          x	mot->seq_no[x]		dat->sequence[x]
	  0	0			0
	  1	1			1
	  2	1			2
	  3	2			4  count = 4 seqs with motif
          4	4 mot->freq = 5		3  NumSeqs = 5
	-------------
	  5     3	(sequences without the motif)
  -------------------------------------------------------------------------*/
  j = 0;         				/* j is sequence counter */
  k = 0;					/* k is duplication counter */
  for (x = 0; x < NumSeqs + mot->dups; x++)
  {
      if (mot->seq_no[j + k] == x - k)
	dat->sequence[j++] = x - k;
      else
      {
	 if (j > 0 && mot->seq_no[j + k] == x - k - 1) k++;
	 else dat->sequence[x-k-j + count] = x-k;
      }
  }
 
  /*-----------------------------------------------------------------------*/
  /*  Compute Smith's "sliding sequence column scores",  aas[aa][i] is
      the average pam score of aa in col i of all occurences of motif */
   for (i=0; i < mot->domain; i++)
   {
      for (aa1=0; aa1 < 20; aa1++)
      {
         score = 0;
         for (x=0; x < mot->freq; x++)
         {
            pos = mot->pos[x] + i - offset;
            if (pos >= 0 && pos < Len[ mot->seq_no[x] ]  )
            {
               aa2 = Seq[ mot->seq_no[x] ][pos];
               score += (long int) (SMatrix->scores[aa1][aa2] + 1) * 100 - SMatrix->highpass;
            }
         }
         if (score > 0) aas[aa1][i] = (int) (score / mot->freq);
         else               aas[aa1][i] = 0;
      }
   }


  /*-----------------------------------------------------------------------*/
  /*   Now realign each sequence in dat->sequence[] order  
       The scores used for the realignment are never modified  */
  for (x = 0; x < NumSeqs; x++) 
  {
      datseq = dat->sequence[x];
      dat->max_score[x] = 0;
      dat->max_pos[x] = -1;

      /* Move along entire length of protein, looking for maximum fit
	 (also records next best fit in max2_pos): */
      for (position = 0; position < Len[datseq]; position++)
      {
	score = 0;
	for (i = 0; i < mot->domain; i++)
        {
          iscore = mot->scores[i] - SMatrix->highpass;
          if (iscore <= 0) iscore = 0;
          pos = position + i - offset;
	  if (pos >= 0 && pos < Len[datseq])
	  {
	      aa1 = Seq[datseq][pos];
	      if ( (aa1 >= 0 && aa1 < 20) && iscore)
		score += (long int) aas[aa1][i] * iscore;
	  }
        }
	if (score > dat->max_score[x]) 
        {
	  dat->max2_score[x] = dat->max_score[x];
	  dat->max2_pos[x] = dat->max_pos[x];
	  dat->max_score[x] = score;
	  dat->max_pos[x] = position;
	}
      } /* end of for position */

      position = dat->max_pos[x];
      /*---  Add info about sequences without the motif            JGH*/
      /*  This assumes the first "count" sequence numbers in
          dat->sequence[] are the unique numbers among the first
          mot->freq entries of mot->seq_no[] (see example above) */
      if (x >= count) 
      { 
	 mot->seq_no[x + mot->dups] = datseq; 
	 mot->pos[x + mot->dups] = position;
      } 
      /*---- Update info about sequences with motif that were realigned */
      else 
      {
          /*  Find first corresponding mot->seq_no[] sequences 
              (can occur more than once if mot->dups > 0  */
          motseq = -1;
          i = 0;
          while (motseq < 0 && i < mot->freq)
          {
             if (mot->seq_no[i] == datseq)
             {
                motseq = datseq;
                dat->motif_pos[x] = mot->pos[i];
                mot->pos[i] = position;
             }
             i++;
         }
      }
 }   /*  end of sequence x */

}  /* end of realign */
/*========================================================================
       Write an output file for motomat
==========================================================================*/
void write_motifs(total_motifs, motif)
int total_motifs;
struct motif_struct *motif;
{
  FILE *mot;
  char temp;
  int i, j;

  /* Open an output file for the motifs:  Batch_Filename.mot               */
  if (strlen(Batch_Filename) != 0)
  {
     Mot_Filename[0] = '\0';
     strncat(Mot_Filename, Batch_Filename+Batsplit->dir_len,
             Batsplit->name_len);
     Mot_Filename[Batsplit->name_len] = '\0';
     strcat(Mot_Filename, ".mot");
  }
  else
     strcpy(Mot_Filename, "motifj.mot");
  if ( (mot=fopen(Mot_Filename, "w+b")) == NULL)
  {
     printf("\nUnable to open file %s", Mot_Filename);               /*JGH*/
     exit(-1);
  }
  i = VERSION;
  fwrite(&i, sizeof(int), 1, mot);
  fwrite(&RunType, sizeof(int), 1, mot);
  fwrite(&Signif, sizeof(int), 1, mot);
  fwrite(&Dups, sizeof(int), 1, mot);
  fwrite(&Distance, sizeof(int), 1, mot);
  fwrite(&NumSeqs, sizeof(int), 1, mot);
  fwrite(&total_motifs, sizeof(int), 1, mot);
  fwrite(Len, sizeof(int), NumSeqs, mot);
  fwrite(Seqname, SNAMELEN*sizeof(char), NumSeqs, mot);
  for (i=0; i<NumSeqs; i++)
     for (j=0; j<Len[i]; j++)
     {
	temp = Seq[i][j];
	fwrite(num_to_aachar(temp), sizeof(char), 1, mot);
     }
  fwrite(motif, sizeof(struct motif_struct), total_motifs, mot);
  fwrite(Title, MAXLINE*sizeof(char), 1, mot);
  fclose(mot);	
}   /* end of write_motifs  */
/*======================================================================
     Write the sequences out to a .pros file
========================================================================*/
void write_pros()
{
  FILE *out;
  char *ctemp, filename[FNAMELEN];
  int i, j, min, max, imin, imax;

  /* Find min/max sequence lengths   */
  imin = imax = max = -1; min = 30000;
  for (i=0; i<NumSeqs; i++)
  {
      if (Len[i] < min) {min = Len[i]; imin = i;}
      if (Len[i] > max) {max = Len[i]; imax = i;}
  }
  /* Open an output file for the pros:  Batch_Filename.pros               */
  if (strlen(Batch_Filename) != 0)
  {
     filename[0] = '\0';
     strncat(filename, Batch_Filename+Batsplit->dir_len,
             Batsplit->name_len);
     filename[Batsplit->name_len] = '\0';
     strcat(filename, ".motifj.pros");
  }
  else
     strcpy(filename, "motifj.pros");
  if ( (out=fopen(filename, "w+b")) == NULL)
  {
     printf("\nUnable to open file %s", filename);               /*JGH*/
  }
  else
  {
     for (i=0; i<NumSeqs; i++)
     {
        fprintf(out, ">%s %d", Seqname[i], Len[i]);
        if (i == imin)      fprintf(out, " MINIMUM\n");
        else if (i == imax) fprintf(out, " MAXIMUM\n");
             else           fprintf(out, "\n");
        for (j=0; j<Len[i]; j++)
        {
	   ctemp = num_to_aachar(Seq[i][j]);
	   fprintf(out, "%1s", ctemp);
	   if ( ((j+1) % 80) == 0) fprintf(out, "\n");
        }
        fprintf(out, "\n");
      }
      fclose(out);
   }
}   /* end of write_pros */
