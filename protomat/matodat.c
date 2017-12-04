/*=====================================================================*/
/*(C) Copyright 1991, Fred Hutchinson Cancer Research Center           */
/*        matodat.c  reads BLIMPS MATRIX search output file
	    eg:  matodat PS00011.hom PS00011.lis PS00011A.blk blocks.dat
			  or
		 matodat PS00011.hom PS00011.lis none
			  or
		 matodat PS00011.hom none none

      If .lis file is provided, looks for those sequences in the .hom file.
      If .blk and .dat files are provided, adds score of 99.5% result
	not in .lis file to .blk and appends block to .dat file.
      Always appends statistics to file named matodat.dat.

       [nl][nl]Consensus...
       [nl]Target...
       [nl][nl]
       [nl][entry# for 13][sp][ID & DE for 60][score for 5]...etc.
       Results start on 7th line.
       Swiss-Prot ID is in columns 16-25.
       Score is in columns 75-79.
*/
/*---------------------------------------------------------------------*/
/*   5/6/91   J. Henikoff
  >>>>>>>>>>>>>>>>>>>>   Blocks 9.0 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  2/ 9/98 Fixed problem re-writing BL line (when input blocks come from
	  addseqs)
  2/20/98 Increased NSCORE from 1000 to 1500
  2/26/98 Treat PS=P in .lis file as fragments; change actually made
	   to get_ids() in motmisc.c
  3/ 1/98 Fix 99.5% on BL line (% wasn't printed)
=====================================================================*/
#include "motifj.h"
#include <math.h>

#define NSCORE 1500		/* Maximum # of scores in .hom file */
#define MAXBUCKET 50		/* Maximum # of frequency buckets */
#define SWISS 38303		/* Number of sequences in SwissProt 29 */

int read_hom();
struct db_id *makedbid();
int get_ids();
void append_block();

char AC[9];				/* Global variables */
int LisSeq, FragSeq, InBlock, Repeats;
double Swiss, Scores, Aligns;

/*======================================================================*/
void main(argc, argv)
int argc;
char *argv[];
{
   FILE *fhom, *flis, *fdat, *fblk;
   char homfile[FNAMELEN], lisfile[FNAMELEN], *ptr, *ptr1;
   char datfile[FNAMELEN], blkfile[FNAMELEN], ctemp[FNAMELEN+20];
   int nids, nscores;
   struct db_id *ids;

   printf("\nMATODAT: (C) Copyright 1991, Fred Hutchinson Cancer ");
   printf("Research Center");
   if (argc < 3)
   {
      printf("\nUSAGE to calibrate blocks:\n");
      printf("     matodat PS00011.bli PS00011.lis PS00011A.blk blocks.dat\n");
      printf("  or matodat PS00011.bli PS00011.lis none\n");
      printf("  or matodat PS00011.bli none none\n");
   }
/*------------- arg 1:  .hom file -------------------------------*/
   if (argc > 1)
      strcpy(homfile, argv[1]);
   else
   {
      printf("\nEnter name of BLIMPS results file: ");
      gets(homfile);
   }
   if ( (fhom=fopen(homfile, "r")) == NULL)
   {
      printf("\nCannot open file %s\n", homfile);
      exit(-1);
   }
   else printf("\nReading %s...", homfile);
   ptr = strrchr(homfile,'/');			/* look for file name */
   if (ptr != NULL)  ptr1 = strtok(ptr+1, ".");
   else  ptr1 = strtok(homfile, ".");
   strcpy(AC, ptr1);
   if (strlen(AC) < 8) strcat(AC, "*");

/*------------- arg 2:  .lis file -------------------------------*/
   if (argc > 2)
      strcpy(lisfile, argv[2]);
   else
   {
      printf("\nEnter name of file containing list of sequences: ");
      gets(lisfile);
   }
   if (!strlen(lisfile)) flis = NULL;
   else if ( (flis=fopen(lisfile, "r")) == NULL)
   {
      printf("\nCannot open file %s", lisfile);
      printf("\nNo list of sequences will be checked.");
   }
/*------------- arg 3:  .blk file -------------------------------*/
   if (argc > 3)
      strcpy(blkfile, argv[3]);
   else
   {
      printf("\nEnter name of block file: ");
      gets(blkfile);
   }
   if (!strlen(blkfile)) fblk = NULL;
   else if ((fblk=fopen(blkfile, "r")) == NULL)
   {
      printf("\nCannot open file %s\n", blkfile);
      printf("\nNo blocks will be updated.");
   }
/*------------- arg 4:  blocks.dat file -------------------------------*/
   if (fblk != NULL)
   {
      if (argc > 4) strcpy(datfile, argv[4]);
      else
      {
	 printf("\nEnter name of blocks database output file: ");
	 gets(datfile);
      }
      if (!strlen(datfile)) fdat = NULL;
      else     /*  blocks database has been specified */
      {
	 if ((fdat=fopen(datfile, "rt")) == NULL)   /* does it exist? */
	 {
	    printf("\nInitializing blocks database %s", datfile);
            sprintf(ctemp, "cp matodat.stp %s", datfile);
            system(ctemp);
	 }
	if ((fdat=fopen(datfile, "a")) == NULL)   /* already exists */
	{
	    printf("\nCannot open file %s\n", datfile);
	    printf("\nNo blocks database will be appended.");
	 }
      }
   }
   else fdat = NULL;
/*---------------  Get the sequences in the .lis file -----------------*/
   nids = 0;
   ids = makedbid();
   if (flis != NULL)
      printf("\n%d IDs in %s", (nids=get_ids(flis, ids)), lisfile);

/*---------------------------------------------------------------------*/
   nscores = read_hom(fhom, fblk, fdat, ids);
   printf("\n%d scores in %s\n", nscores, homfile);
   fclose(fhom);
   if (fblk != NULL) fclose(fblk);
   if (fdat != NULL) fclose(fdat);
   printf("\n");
   exit(0);
}  /*  end of main */
/*======================================================================*/
int read_hom(fhom, fblk, fdat, ids)
FILE *fhom, *fblk, *fdat;
struct db_id *ids;
{
   FILE *fout;			/* statistics */
   int tns, tn[NSCORE], minscore, maxscore, i, rank, tscore;
   int tp[NSCORE], lfreq[MAXBUCKET], lminrank, lmaxrank;
   int iminrank, imaxrank, tpmedian, tp995, pearson;
   int found, freq[MAXBUCKET], minseq, maxseq, tps, x, lower, upper;
   int fragment, frags, inblock, inblocks;
   double tenth, roc;
   char line[MAXLINE], query[FNAMELEN], db[FNAMELEN], id[12], score[8];
   char ctemp[MAXLINE];
   struct db_id *did;

   /*------------------Intialize -------------------------------------*/
   for (i=0; i<MAXBUCKET; i++)
   { freq[i] = lfreq[i] = 0; }
   tns = 0; minscore=9999; maxscore=0; rank=0;
   minseq = 9999; maxseq = 0; tps = frags = inblocks = 0;
   tp[0] = tn[0] = 0;
   lminrank = iminrank = 9999; lmaxrank = imaxrank = 0;
   query[0] = db[0] = ctemp[0] = '\0';

   /*----------- Count the true positive sequences -------------------*/
   LisSeq = FragSeq = InBlock = 0;
   did = ids->next;
   while (did != NULL) 
   {
      if (!did->frag && !did->block) LisSeq++;
      else if (did->frag)            FragSeq++;
      else if (did->block)           InBlock++;
      did = did->next;
   }

   /*---------Get the 99.5% score ------------------------------------*/
   Swiss = (double) SWISS;
   lower = (int) (SWISS-LisSeq-FragSeq-InBlock)*.005;
   printf("\n%d true positives, %d fragments, %d in block,",
          LisSeq, FragSeq, InBlock);
   printf(" default 99.5th rank = %d", lower);
   roc = 0.0;
   Scores = Aligns = 0.0;
   Repeats = NO;

   /*------------- Read the search results ---------------------------*/
   while (!feof(fhom) && fgets(line,MAXLINE,fhom) != NULL)
   {
      if (strlen(line) > 11 && strncmp(line, "Block", 5) == 0)
      {   /*  "Block File: query" */
	 strcpy(query, &line[11]);
      }
      else if (strlen(line) > 12 && strncmp(line, "Matrix", 6) == 0)
      {   /*  "Matrix File: query" */
	 strcpy(query, &line[12]);
      } 
      else if (strlen(line) > 18 && strncmp(line, "Target", 6) == 0)
      {   /*  "Target File (s) : db" */
	 strcpy(db, &line[18]);
      }
      else if (strlen(line) > 17 && strncmp(line, "Records", 7) == 0)
      {  /* "Records Searched: n" */
         strcpy(ctemp, &line[17]);
         Swiss = atof(ctemp);
	 if (Swiss > 0)
         {
            lower = (int) (Swiss-LisSeq-FragSeq-InBlock)*.005;
          printf("\n99.5th rank adjusted to %d based on %.0f records searched",
                		lower, Swiss);
         }
      }
      else if (strlen(line) > 12 && strncmp(line, "Scores", 6) == 0)
      {  /* "Scores Done: n" */
         strcpy(ctemp, &line[12]);
         Scores = atof(ctemp);
         if (Scores > Swiss)
         {
             Repeats = YES;
             printf("\nWARNING: Blimps search was done with repeats on,");
             printf("\n         true positives may be counted more than once.");
         }
      }
      else if (strlen(line) > 16 && strncmp(line, "Alignments", 10) == 0)
      {  /* "Alignments Done: n" */
         strcpy(ctemp, &line[16]);
         Aligns = atof(ctemp);
      }
      else if (strlen(line) > 96)
      {
	 rank += 1;
/*   Old way:  got id from description field
         if (line[14] == '>')
	    strncpy(id, &line[15], IDLEN);  
         else
	    strncpy(id, &line[14], IDLEN);
*/
         /*   Get id */
         if (line[0] == '>')
	    strncpy(id, &line[1], IDLEN);  
         else
	    strncpy(id, &line[0], IDLEN);
         id[IDLEN] = '\0';
         for (i=0; i < IDLEN; i++)
            if (id[i] == ' ') id[i] = '\0';
	 strncpy(score, &line[75], 4); score[4] = '\0';		/* score */
	 tscore = atoi(score);
	 /*----- Change % to $ for comparison & correct length ---*/
	 for (i=0; i<strlen(id); i++)
	 {
	    if (id[i] == '%') id[i] = '$';
	    else if (id[i] == ' ' || id[i] == '\n') id[i] = '\0';
	 }
         /*  now check to see if id is in the .lis file...  */
	 found = fragment = inblock = NO;
	 did = ids->next;
	 while (did != NULL && 
		strcmp(did->entry, id) <= 0)
	 {
            if (!did->found && strcmp(did->entry, id) == 0)
	    {
	       did->rank = rank;
	       did->score = tscore;
               did->found = found = YES;
               if (!did->frag && !did->block)
               {
                  i = (int) tscore/100;
	          if (i >= 0 && i < MAXBUCKET) ++lfreq[i];
	          tp[tps] = tscore;
	          tps++;
	          if (rank < lminrank) lminrank = rank;
	          if (rank > lmaxrank) lmaxrank = rank;
	          if (tscore < minseq) minseq = tscore;
	          if (tscore > maxseq) maxseq = tscore;
	          printf("\nTP: rank %.3d=%.10s, score=%.4d",
                          rank, id, tscore);
	       }
               else if (did->frag)
               {
                  fragment = YES;
                  frags++;
	          printf("\n   FRAGMENT: rank %.3d=%.10s, score=%.4d",
                                rank, id, tscore);
               }
               else if (did->block)
               {
                  inblock = YES;
                  inblocks++;
	          printf("\n   IN BLOCK: rank %.3d=%.10s, score=%.4d",
                                rank, id, tscore);
               }
            }
	    did = did->next;
	 }

	 if (!found && tns < NSCORE)      /*  Seqs not in list */
	 {
            i = (int) tscore/100;
	    if (i >= 0 && i < MAXBUCKET) ++freq[i];
	    tn[tns] = tscore;
	    if (tn[tns] < minscore) minscore = tn[tns];
	    if (tn[tns] > maxscore) maxscore = tn[tns];
	    if (rank < iminrank) iminrank = rank;
	    if (rank > imaxrank) imaxrank = rank;
	    if (tns==0)   printf("\n   1st TN rank=%d score=%d", rank, tn[0]);
	    else if ( ((tns+1)%100)==0 )
	       printf("\n   %dth TN rank=%d score=%d", tns+1, rank, tn[tns]);
	    else if ((tns+1) == lower)
	       printf("\n   %dth (99.5) rank=%d TN score=%d",
                       tns+1, rank, tn[tns]);
	    /*  Print names of other top seqs */
	    if (rank <= LisSeq)
	    {
	       printf("\n   TN: rank %.3d=%.10s, score=%.4d",
			 rank, id, tscore);
	    }
	    tns++;
            roc += (double) tps / (double) LisSeq;
	 }
      } /* end of search result */
   } /* end of file */
   if (tns > 0) roc /= (double) tns;
   else roc = 1.0;

   /*--------------------------Print the results--------------------*/
   if (lower > tns-1) lower = tns-1;	/* adjust 99.55 TN if nec. */
   if (lower < 0) lower = 0;
   /*-- For median tp, use LisSeq/2, not tps/2, because some tp
	sequences may not have been found --*/
   i = (int) LisSeq/2;
   if (tps == 0) tpmedian=0;	/* no true positives found */
   else if (i >= tps) 		/* less than half true positives found */
       tpmedian = tn[tns-1];  		/* assume the worst */
   else tpmedian = tp[i];
   printf("\n%d of %d true positive scores found,", tps, LisSeq);
   printf(" %d of %d in block,", inblocks, InBlock);
   printf(" %d of %d fragments,\n       %d other scores", frags, FragSeq, tns);
   if (tps > 0)
      printf("\nSummary of true positive scores:  Maximum=%d, Median=%d, Minimum=%d",
	     tp[0], tpmedian, tp[tps-1]);
   else printf ("\nNo true positive sequences found.");
   printf("\nSummary of other scores: Maximum=%d, Median=%d, Minimum=%d, 99th=%d",
      tn[0], tn[(int) tns/2], tn[tns-1],
      tn[(int) tns/100]);
   printf("\nFrequency distribution of all scores:");
   printf("(*+=TP, x/=TN)");
   if (minseq < minscore) minscore = minseq;
   if (maxseq > maxscore) maxscore = maxseq;
   for (i=(int) minscore/100; i<= (int) maxscore/100; i++)
   {
      printf("\n %.4d-%.4d  %.4d  ", i*100, i*100+99, lfreq[i]+freq[i]);
      for (x=0; x< (int) lfreq[i]/5; x++)  printf("*");
      if ((freq[i]-5*x) > 0) printf("+");
      for (x=0; x< (int) freq[i]/5; x++)  printf("x");
      if ((freq[i]-5*x) > 0) printf("/");
   }
/*----------Append statistics to matodat.dat ----------------------------*/
   if ( (fout=fopen("matodat.dat", "a")) == NULL)
   {
      printf("\nCannot open file matodat.dat\n");
      printf("\nNo statistics will be saved.");
   }
   if (fout != NULL)
   {
      fprintf(fout, "%.7s %.1s %d", AC, AC+7, LisSeq);
      /*  Median CI for z=2: order stat = (n+1)/2 - z*sqrt(n)/2.  */
      if (tps > 5)
      {
	 upper = (float) tps; upper=sqrt(upper);
	 upper += (float) (tps+1)/2;
      }
      else upper = (float) tps-1;
      if (upper < 0) upper= 0; if (upper > tps-1) upper=tps-1;
      /*  Info for seqs in the lis file (true positives)*/
      if (tps > 0)
         fprintf(fout, " %d %d %d %d %d %d %d",
	   tps, lminrank, lmaxrank, tp[(int) upper],
	   tp[tps-1], tpmedian, tp[0]);
      else
         fprintf(fout, " 0 0 0 0 0 0 0");
      /*  Now add decile scores for seqs in the lis file */
      tenth = (float) tps/10;
      for (i=9; i>=0; i--)
      {
	 x = (int) tenth/2 + i*tenth;
	 if (x<0) x = 0; if (x>tps-1) x=tps-1;
/*	 fprintf(fout, " %d", tp[x]); */
      }
      /*  Info for seqs not in the lis file (true negatives) */
      if (tns > 0)
         fprintf(fout, " %d %d %d %d %d %d %d",
	   tns, iminrank, imaxrank, tn[lower],
	   tn[tns-1], tn[(int) tns/2], tn[0]);
      else
         fprintf(fout, " 0 0 0 0 0 0 0");
      /* Number of tps above 99.5% of tns */
      tp995=0;
      for (i=0; i<tps; i++)
	 if (tp[i] >= tn[lower]) tp995++;
      printf("\n%d of %d true positive scores at or above", tp995, tps);
      printf("\nthe 99.5 percent true negative score of %d", tn[lower]);
      pearson = LisSeq - tps;
      i = tps - 1;
      while (pearson >= 0 && pearson < tns &&
             i >= 0 && i < tps &&
             tn[pearson] > tp[i])
      { pearson++; i--; }
      printf("\n%d = Pearson equivalence number (%d/%d)",
              pearson, LisSeq - pearson, LisSeq);
      fprintf(fout, " %d %d", tp995, pearson);
      printf("\nROC area = %f", roc);
      fprintf(fout, " %f", roc);
   }

/*---------------Put .5% score in BL line of .blk file and append it
      to blocks.dat, tack on block data to matodat.dat --------------*/
/*------- Put median TP score in block -------------------*/
   if (fblk != NULL && fdat != NULL)
      append_block(fblk, fdat, fout, tn[lower], tpmedian);
   if (fout != NULL)
   { fprintf(fout, "\n"); fclose(fout); }
   return(tns+tps+frags+inblocks);
}  /* end of read_hom */
/*======================================================================
       Append the block to the database, changing some of the text
       and adding the score.

BL could look like:
From motomat: "BL   FWN motif=[18,7,17] motomat=[1,80,-10] width=11 seqs=19"
From addseqs: 
"BL   FWN motif=[18,7,17] motomat=[1,80,-10] width=11; seqs=20; 99.5%=559; strength=846"
Want: "BL   FWN motif; width=11; seqs=20; 99.5%=575; strength=1445
=======================================================================*/
void append_block(fblk, fdat, fout, lower, upper)
FILE *fblk, *fdat, *fout;
int lower, upper;
{
   char line[MAXLINE], saveline[MAXLINE], *ptr, *ptr1, ntemp[10];
   char savemot[20], stemp[MAXLINE];
   int lastword, len, width, seqs, minprev, maxprev, dups, strength;

   width=seqs=minprev=maxprev=dups = strength = -1;

   while ( fgets(line, MAXLINE, fblk) != NULL)
   {
	 if (strncmp(line, "AC   ", 5)==0)
	 {       /* Change PS to BL */
	    if (strncmp(line+5, "PS", 2)==0)
	    { line[5]='B'; line[6]='L';  }
	    strcpy(saveline,line);	/* pick up minprev & maxprev */
	    ptr=strtok(saveline+5, "(");
	    if (ptr != NULL) ptr=strtok(NULL, ",");
	    if (ptr != NULL)
	    {  minprev = atoi(ptr);
	       ptr=strtok(NULL, ")");
	       if (ptr != NULL) maxprev = atoi(ptr);
	    }
	 }
	 else if (strncmp(line, "DE   ", 5)==0)
	 {
	    strcpy(saveline, line); line[5] = '\0';
	    lastword = 0;     /* 1="active" */
			      /* 2="active site" or "domain" or "protein" */
			      /* 3="signature" */
	    ptr = strtok(saveline+5, " \t\n");
	    while (ptr != NULL)
	    {
	       if (lastword == 1 && strncmp(ptr, "site", 4) != 0)
		  strcat(line, "active ");
	       if (strcmp(ptr, "active")==0)
		  lastword = 1;
	       else if (strncmp(ptr, "site", 4)==0)
	       {
		  if (lastword == 1)
		  {  strcat(line, "proteins "); lastword = 2; }
		  else
		  {  strcat(line, "site proteins "); lastword = 0; }
	       }
	       else if (strncmp(ptr, "signature", 9) == 0)
	       {
		  if (lastword == 2)  strcat(line, " ");
		  else                strcat(line, "proteins ");
		  lastword = 3;
	       }
	       else if (strncmp(ptr, "1", 1) == 0 && lastword == 3)
	       { strcat(line, " "); lastword = 0; }
	       else if (strncmp(ptr, "domain", 6)==0)
	       { strcat(line, "domain proteins "); lastword = 2; }
	       else if (strncmp(ptr, "protein", 7)==0)
	       { strcat(line, ptr); strcat(line, " "); lastword = 2;}
	       else if (strncmp(ptr, "repeat", 6)==0)
	       { strcat(line, "repeat proteins "); lastword = 0; }
	       else if (strncmp(ptr, "sequence", 8)==0)
	       { strcat(line, "sequence proteins "); lastword = 0; }
	       else if (strncmp(ptr, "putative", 8)==0)
	       { strcat(line, " "); lastword = 0; }
	       else
	       {  strcat(line, ptr); strcat(line, " "); lastword = 0;}
	       ptr = strtok(NULL,  " \t\n");
	    }
	    /*---- Remove extra spaces at end of line & add period ---*/
	    len = strlen(line) - 1;
	    while (len >= 0 && (line[len] == '.' || line[len] == ' '))
	       line[len--] = '\0';
	    if (line[len] != '.') strcat(line, ".");
	    strcat(line, "\n");
	 }
	 if (strncmp(line, "BL   ", 5)==0 ||
             strncmp(line, "MA   ", 5)==0  )
	 {
            /*   Get width & #seqs */
	    strcpy(saveline, line);
	    ptr=strstr(saveline, "width=");
	    ptr1=strtok(ptr,"=");
	    if (ptr1 != NULL) ptr1 = strtok(NULL, ";");
	    width=atoi(ptr1);

	    strcpy(saveline, line);
	    ptr=strstr(saveline, "seqs=");
	    ptr1=strtok(ptr,"=");
	    if (ptr1 != NULL) ptr1 = strtok(NULL, ";");
	    seqs=atoi(ptr1);

	    /*  Remove motif & motomat parameters before writing it out */
            /*  123456789
                BL   XXX motif     */
	    strcpy(saveline, line); line[8] = '\0';
	    ptr = strtok(saveline+8, " \t\n");
	    while (ptr != NULL)
	    {
	       if (strncmp(ptr, "motif", 5)==0)
	       {
		  strcat(line, " motif; ");
		  strcpy(savemot, ptr);
	       }
               else
               {
                  if (strncmp(ptr, "gibbs", 5)==0)
                  {
                     strcat(line, " gibbs; ");
		     strcpy(savemot, ptr);
                  }
               }
	       ptr = strtok(NULL,  " \t\n");
	    }
            if (lower > 0) strength = (int) ((float) 1000*upper/lower);
            else           strength = 0;
            sprintf(stemp, "width=%d; seqs=%d; 99.5%%=%d; strength=%d\n",
                    width, seqs, lower, strength);
	    strcat(line, stemp); 

	    /*---------------- pick up dups--------------*/
	    ptr=strtok(savemot, ",");
	    if (ptr !=NULL)
	    {
	       ptr=strtok(NULL, ",");
	       if (ptr != NULL) dups=atoi(ptr);
	    }
	 }  /*  end of BL */
	 fputs(line, fdat);
   }  /* end of fblk */

   /*   Put stuff in matodat.dat */
   fprintf(fout, " %d %d %d %d", minprev, maxprev, width, dups);
}  /* end of append_block */
