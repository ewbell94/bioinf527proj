/*=====================================================================*/
/*(C) Copyright 1991 by Fred Hutchinson Cancer Research Center         */
/*        uextract.c   Un-indexed extract from GENBANK, EMBL, .UNI     */
/*   Reads through ENTIRE list of requested entries for each db entry. */
/*     USE: uextract lisname dbname [-f -ofilename]                    */
/*        dbname = name of database file, program will determine its type
	  lisname = name of a file containing a list of entries to
	    be extracted.  For the various file types this should be:
	      GENBANK:  LOCUS name (FASTA type 1)
	      EMBL:     ID name (FASTA type 3)
	      PIR:	ENTRY name (NBRF/CODATA format, FASTA type 2))
	      VMS:	name immediately following the ; on the first line.
			(NBRF/VMS format, FASTA type 5)
	      UNI:      name immediately following the > on the first line.
			(FASTA type 0)
           -f => extract fragments (not extracted by default)
           -ofilename => put all sequences in file filename (extracted
                         to separate files by default)
           -n => don't execute motifj (executed by default if MOTIF
                 is on title line of .lis file

    * The "lisname" file should have a title line starting with ">".
      The second line may contain the name of a directory where the
      extracted sequences will be stored; it will be created if it
      does not exist. Subsequent lines should contain the sequence
      entry IDs, one per line. EG:
		>PS00094 ;Cytosine methylases
		c:\proteins\PS00094\
		MTB1$BACSH
		MTB1$BACSU
		MTB1$BREEP

     * One output file is created for each entry in "universal" format
	(>title $, then sequence, then *). The name of the file
	is the 1st 8 characters of the entry name followed by .dna for
	GENBANK and by .pro for EMBL and UNI. EG: "MTB1$BAC.pro".
	NOTE FOR DOS: IF THE ENTRY NAME IS > 8 CHARACTERS THE FILE NAME
	MAY NOT BE UNIQUE SINCE IT IS TRUNCATED!!
     * Another file with a ".lst" extension is created containing the
       names of sequences in lisname which were found and do not appear
       to be fragments (the word FRAGMENT does not appear in the
       description for the sequence). In addition, if the IDs. in the
       lisname  file are of SWISS-PROT format, aaaa$aaaaa, then only
       the longest of all the sequences with the same characters before
       the $ and with the first three characters following the $ the
       same is written to the .lst file. For example, if lisname contains
       MTB1$BACSU and MTB1$BACSH and MTB1$BREEP, and MTB1$BACSU is longer
       than MTB1$BACSH, then MTB1$BACSU and MTB1$BREEP will be written
       to .lst, but not MTB1$BACSH. See lst_list().
     * A statistics file named "uextract.dat" is created.
  KNOWN PROBLEMS:  *If input database is not sorted by ID, may miss some
		    requested entries.
		   *If input list file has extension .lst, overwrites it.
		   *VMS format does not have fragment information.
--------------------------------------------------------------------------*/
/*
>>>>>>>>>>>>>>>>>>>> Blocks 7.x, 8.0 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 2/19/95  Don't create directory if filename is specified.
          Don't activate -n if filename is specified, works okay now.
	  Don't write .lst or rewrite .lis files if -n (don't run motif)
11/14/95  Added usage;  update title signif, etc. before writing .lis/.lst
=========================================================================*/
#include <sys/types.h>
/*#include <sys/dir.h>
*/
#include <dirent.h>
#include "motifj.h"

int get_ids();
int lst_list();
int lis_list();
/*--------------  Routines from motmisc.obj --------------------------*/
void init_dbs();
int type_dbs();
int extract_seqs();
struct db_id *makedbid();
char *dir_unix();
struct split_name *split_names();

char Pros[FNAMELEN];

/*======================================================================*/
void main(argc, argv)
int argc;
char *argv[];
{
   char infile[FNAMELEN], lisfile[FNAMELEN], title[MAXLINE], stemp[MAXLINE];
   char foutname[FNAMELEN], stitle[MAXLINE];
   char lstname[FNAMELEN], *temp, *ptr;
   char runtype[3], signif[5], dups[5], distance[5];
   struct db_info *dbs[MAXDB];
   struct db_id *ids;
   int arg, totseqs, nids, test, totlst, reseqs, frag, motif, i, j;
   FILE *fin, *flis, *flst, *fdat, *fout;
   struct split_name *lissplit;

   printf("\nUEXTRACT: (C) Copyright 1991 by Fred Hutchinson");
   printf(" Cancer Research Center\n");
   if (argc < 3)
   {
      printf("USAGE: uextract <lisfile> <dbfile> [-f -o<outfile> -n]\n");
      printf("       <lisfile> = file listing sequence names to extract\n");
      printf("                   in PROTOMAT format\n");
      printf("       <dbfile>  = sequence database\n");
      printf("       -f          to extract FRAGMENTs\n");
      printf("       -o<outfile> to have all sequences put in <outfile>\n");
      printf("       -n          to NOT execute motifj automatically\n");
   }
/*-------------  arg 1.  List of entries to extract ---------------------*/
   if (argc > 1)
      strcpy(lisfile, argv[1]);
   else
   {
      printf("\nEnter name of file containing list of entries to extract: ");
      gets(lisfile);
   }
   if ( (flis=fopen(lisfile, "r")) == NULL)
   {
      printf("\nCannot open file %s\n", lisfile);
      exit(-1);
   }
   lissplit=split_names(lisfile);
/*----------------- arg 2 database name --------------------------------*/
   if (argc > 2)
      strcpy(infile, argv[2]);
   else
   {
      printf("\nEnter name of database file to extract entries from: ");
      gets(infile);
   }
   if ( (fin=fopen(infile, "r")) == NULL)
   {
      printf("\nCannot open file %s\n", infile);
      exit(-1);
   }
/*----------------- args 3+ extract fragments---------------------------*/
   frag = NO;
   foutname[0] = '\0';
   motif = YES;
   if (argc > 3)
   {
      for (arg=3; arg < argc; arg++)
      {
         if (argv[arg][0] == '-')
            switch(argv[arg][1])
            {  
               case 'f': frag = YES; break;
               case 'o': strcpy(foutname, argv[arg]+2); break;
               case 'n': motif = NO; break;
               default: break;
            }
      }
   }
   if (strlen(foutname))
   {
/*    motif = NO;		 motifj won't work with this option */
      if ( (fout=fopen(foutname, "w+t")) == NULL)
      {
         printf("\nCannot open file %s\n", foutname);
         exit(-1);
      }
      else
         printf("\nExtracting all sequences to %s", foutname);
   }
   else fout = NULL;

/*-------------  First line of extract file may have a title ----------*/
   fgets(title, MAXLINE, flis);
   if (strlen(title) && title[0] != '>')
   {
      title[0] = '\0';
      rewind(flis);
   }
   else       /*  Don't print the MOTIFJ= part, confuses people ! */
   {
      strcpy(stitle, title);  ptr = strstr(stitle, "MOTIFJ=[");
      if (ptr != NULL) stitle[strlen(stitle)-strlen(ptr)] = '\0';
      printf("\n%s", stitle);
/*      if (fout != NULL) fprintf(fout, "%s\n", stitle);  gibbsj loops! */
   }
/*------------Second line may have a directory name --------------*/
   Pros[0] = '\0';
   if (getcwd(Pros, FNAMELEN) != NULL) strcat(Pros, "/");           /*DOS*/
   fgets(stemp, MAXLINE, flis);
   if (fout == NULL && strlen(stemp) && stemp[0] != '>' &&
       strstr(stemp, "/") != NULL)                                /*DOS*/
         strcpy(Pros,dir_unix(stemp));      /* create the directory if nec. */

/*--------------- Get extract list --------------------------------------*/
   ids = makedbid();
   nids = get_ids(flis, ids);
   fclose(flis);
   printf("\n%d sequences requested for extract from %s", nids, infile);
   if (fout == NULL)
      printf("\n  and deposit into directory %s", Pros);
   else
      printf("\n  and deposit into file %s", foutname);

/*------------------- Extract the sequences ---------------------------*/
   if (nids > 0)
   {
      init_dbs(dbs);		       /* load database infor. */
      totseqs = extract_seqs(nids, dbs, fin, ids, Pros, fout, frag);
   }
   else  totseqs = 0;
   printf("\n%d sequences extracted\n",  totseqs);

   fclose(fin);


   /*-----------If they want to run motif, do a lot of stuff now ------*/
   if (motif && totseqs > 0)
   {
      /*--------  Make the .lst file in the current directory----------*/
      lstname[0] = '\0';
      strncat(lstname, lisfile+lissplit->dir_len, lissplit->name_len);
      lstname[lissplit->name_len] = '\0';
      strcat(lstname, ".lst");
      if ( (flst=fopen(lstname, "w+t")) == NULL)
      {
	 printf("\nCannot open file %s\n", lstname);
	 exit(-1);
      }
      printf("%s\n",title);

      if (strlen(title) > 2) fprintf(flst, "%s", title);
      if (strlen(Pros) > 2)  fprintf(flst, "%s\n", Pros);
      totlst = lst_list(flst, ids);
      printf("\n%d entries written to %s\n", totlst, lstname);
      fclose(flst);

      /*-----  Update the title line for the .lis file ------------*/
      /*  NOTE: should do this in the .lst file, too, but requires totlst */
      /*   See if any info. from PROTOMOT; just use runtype & dups    */
      temp = strstr(title, "MOTIFJ=[");
      if (temp != NULL)
      {
         ptr = strtok(temp, "[");
         ptr = strtok(NULL, ",");
         strcpy(runtype, ptr);
         ptr = strtok(NULL, ",");
         strcpy(signif, ptr);
         ptr = strtok(NULL, ",");
         strcpy(dups, ptr);
         test = atoi(dups);
         if (test < 0) test = 0;		/* Default is zero dups */
         kr_itoa(test, dups, 10);
         ptr = strtok(NULL, "]");
         strcpy(distance, ptr);
      }
      else				/* use defaults */
      {
         strcpy(runtype, "4");
         strcpy(dups, "0");
      }
      /*----  redo distance & signif;  don't use PROTOMOT stuff!---*/
      strcpy(distance, "17");		/* Force distance=17 */
      /*--     totlst=#seqs in .lst; totseqs=#seqs in .lis     --*/
      if (strcmp(runtype, "4") == 0)
	       test = (int) ((totlst+1)/2) + 1; /*  start at NUMSEQS/2 +1 */
      else test = totlst;
      if (test <= MINSEQS)  test = MINSEQS + 1;
      kr_itoa(test, signif, 10);
      /* Re-write MOTIFJ=[] stuff in title now!   */
      sprintf(title, "%s MOTIFJ=[%s,%s,%s,%s];$\n", 
               stitle, runtype, signif, dups, distance);

      /*---- Rewrite the .lis file ------------------------------------*/
      if ( (flis=fopen(lisfile, "w+t")) == NULL)
      {
	 printf("\nCannot open file %s for update\n", lisfile);
      }
      else
      {
	 if (strlen(title) > 2) fprintf(flis, "%s", title);
	 if (strlen(Pros) > 2)  fprintf(flis, "%s\n", Pros);
	 reseqs = lis_list(flst, ids);
	 printf("\n%d entries re-written to %s\n", reseqs, lisfile);
	 fclose(flis);
      }

      /*-----------------Write stats to the .dat file --------------------*/
      if ( (fdat=fopen("uextract.dat", "a")) != NULL)
      {
            fprintf(fdat, "%s %d %d %s %s %s\n",
	        lisfile, totseqs, totlst, signif, dups, distance);
            fclose(fdat);
      }
      /*------Execute motifj----------------------------------------*/
      /*------motifj will read either a list of sequences or a file
         of sequences if the file name is preceded with a "-" --------*/
      if (fout == NULL) strcpy(stemp, lstname);
      else
      {  strcpy(stemp, "-"); strcat(stemp, foutname);
            fclose(fout);
      }
      execlp("motifj", "motifj", runtype, stemp,
		signif, dups, distance, NULL);
   }  /*  end of if motif  */

   exit(0);
}  /*  end of main */

/*======================================================================*/
/*   Build the .lst file for motifj:
     The sequences in the .lis file are all included EXCEPT:
	1. If one was not extracted (id->found=NO)
	2. If one was found to be a fragment (id->frag=1)
	3. If two or more sequences appear to be from SWISS-PROT
	   (id contains a $ or _) and their ids are identical before the
	   $ and for 3 characters following the $, only the longest of
	   the group is included.                                       */
/*======================================================================*/
int lst_list(flst, ids)
FILE *flst;
struct db_id *ids;
{
   struct db_id *id, *save;
   char sid[20], tid[20], sid1[10];
   int maxlen, nlst, doll;
   
   sid1[0] = sid[0] = '\0';
   nlst = 0;
   id = ids->next;				/* initialize 1st set */
   /*----  Skip to first non-fragment, non-P sequence ----*/
   while (id != NULL && 
          (id->frag || !id->found || strcmp(id->ps, "P")==0)) 
                 id = id->next;
   if (id == NULL) return(0);   /*  No .lst sequences */

   strcpy(sid, id->entry);  maxlen = id->len;  save=id;
   /*-----  Copy up to $ or _, plus 3 more characters ---------*/
   strcpy(sid1, sid);
   doll = strcspn(sid1, "$_");
      if ( (doll+4) < strlen(sid1) ) sid1[doll+4] = '\0';
   /*---- Now do the rest of the entries -------------------*/
   while (id != NULL &&
             (id->frag || !id->found || strcmp(id->ps, "P")==0)) 
                   id=id->next;
   while (id != NULL)
   {
      strcpy(tid, id->entry);
      doll = strcspn(tid, "$_");
      if ( (doll+4) < strlen(tid) ) tid[doll+4] = '\0';
      if (strcmp(tid, sid1) != 0)      /*  New set */
      {
	 nlst += 1;
	 /*  Print the previous set */
	 fprintf(flst, "%-12s", sid);
	 if (strlen(save->ps)) fprintf(flst, "  PS=%s", save->ps);
	 if (save->len > 0) fprintf(flst, "  LENGTH=%-6d", save->len);
	 fprintf(flst, "\n");
	 save->lst = YES;
	 strcpy(sid, id->entry);        /*  Initialize the next set */
	 save = id;
	 strcpy(sid1, tid);
	 maxlen = id->len;
      }
      else if (strcmp(tid, sid1) == 0 && id->len > maxlen)
      {
	 strcpy(sid, id->entry); save = id;
	 maxlen = id->len;
      }
      /*----  Skip to next non-fragment, non-P sequence ------*/
      id = id->next;
      while (id != NULL && 
             (id->frag || !id->found || strcmp(id->ps, "P")==0)) 
                id = id->next;
   }
   nlst += 1;
   /*  Print the last set */
   fprintf(flst, "%-12s", sid);
   if (strlen(save->ps)) fprintf(flst,"  PS=%s", save->ps);
   if (save->len > 0) fprintf(flst,"  LENGTH=%-6d", save->len);
   fprintf(flst,"\n");
   save->lst = YES;

   return(nlst);

}   /*  end of lst_list */
/*======================================================================*/
int lis_list(flis, ids)
FILE *flis;
struct db_id *ids;
{
  
   struct db_id *id;
   int nlis, i;

   nlis = 0;
   id = ids->next;
   while (id != NULL)
   {
       nlis += 1;
       /*   Change % back to $ */
       for (i=0; i<strlen(id->entry); i++)
	  if (id->entry[i] == '%') id->entry[i] = '$';
       fprintf(flis, "%-12s", id->entry);
       if (strlen(id->ps)) fprintf(flis, "  PS=%s", id->ps);
       if (id->len > 0 ) fprintf(flis, "  LENGTH=%-6d", id->len);
       if (id->frag)     fprintf(flis, "  FRAGMENT");
       if (id->lst)      fprintf(flis, "  LST");
       if (strlen(id->pir)) fprintf(flis, " PIR=%s", id->pir);
       fprintf(flis, "\n");
       id = id->next;
   }
   return(nlis);
}   /*  end of lis_list */
