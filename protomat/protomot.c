/*=====================================================================*/
/*(C) Copyright 1991 by Fred Hutchinson Cancer Research Center         */
/*        protomot.c Extracts entries from PROSITE.DAT                 */
/*      USE: protomot prosite AC name [swiss pros]
	 prosite: name of prosite.dat file
	 AC:      accession # in prosite.dat file, or "all" for all
		  accession numbers in the file
	 name:    short (7 chars) name Prefix for output files
	 swiss:   name of corresponding SWISS-PROT database
	 pros:    directory for deposition of proteins                */
/*    Reads prosite.dat as an EMBL database.  Looks for patterns, and
      associated protein sequences:
	->start = ID   name; PATTERN.
		  AC   Prosite accession #.
	->desc  = DE   description.
		  PA   pattern.  (can be > 1 line)
		  MA   matrix. (can be > 1 line)
	->seq   = DR   P07507, MGP$BOVIN , T; (etc.)
	->end   = //                                                   */
/*   Outputs file "name.lis" containing a list of Swiss-Prot IDs
     for true positive and false negative sequences in the PROSITE
     entry. The title line of the ".lis" file contains parameter
     information for the MOTIFJ program. If the swiss name is provided,
     executes UEXTRACT.                                                */
/*---------------------------------------------------------------------*/
/*    8/9/90   J. Henikoff
>>>>>>>>>>>>>>>>>>>>>>>> Blocks 8.0 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    3/10/94  Added check for MAX-REPEAT. 
    7/31/94  Look for MA entries in Prosite (new with release 12)
   12/ 5/95 Modified format of /FALSE_NEG on NR line (changed w/ release 13)
>>>>>>>>>>>>>>>>>>>>>>>> Blocks 9.0 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   12/12/95 Changed dups from MAX-REPEAT - 1 to nhit - npos; print nhit/npos.
>>>>>>>>>>>>>>>>>>>>>>>> Blocks 10.0 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    2/27/98 Create repeats.dat file: (MAX-REPEAT - 1)
=========================================================================*/
#include <sys/types.h>
/*#include <sys/dir.h>*/
#include <dirent.h> 
#include "motifj.h"

void init_dbs();
int type_dbs();
int strnjcmp();

int get_pat();
int find_swiss();
int screen_pat();
int classify();

char PatName[15], Prefix[15], Pros[FNAMELEN];

/*======================================================================*/
void main(argc, argv)
int argc;
char *argv[];
{
   FILE *fdat, *frep;
   char infile[FNAMELEN], defname[FNAMELEN], swiss[FNAMELEN];
   struct db_info *dbs[MAXDB];
   int totpats;

   printf("\nPROTOMOT: (C) Copyright 1991 by Fred Hutchinson Cancer Research Center");
   printf("\nFor PROSITE, please cite:  A. Bairoch, PROSITE:");
   printf(" A Dictionary of Protein\n Sites and Patterns,");
   printf(" U. Geneva, 5th Release (1990).\n");

/*------------- arg 1:  prosite.dat file -------------------------------*/
   if (argc > 1)
      strcpy(infile, argv[1]);
   else
   {
      strcpy(defname, "/cdrom/prosite/prosite.dat");
      printf("\nEnter name of ProSite data file [%s]:\n", defname);
      gets(infile);
      if ((int) strlen(infile) < 2)
	 strcpy(infile, defname);
   }
   if ( (fdat=fopen(infile, "r")) == NULL)
   {
      printf("\nCannot open file %s\n", infile);
      exit(-1);
   }
/*------------------- arg 2:  Accession # in prosite.dat ----------------*/
   if (argc > 2)
      strcpy(PatName, argv[2]);
   else
   {
      printf("\nEnter accession number of ProSite pattern [all]: ");
      gets(PatName);
   }
   if ((int) strlen(PatName) < 2)   strcpy(PatName, "all");
   else
   {
      if (PatName[0] == 'p') PatName[0] = 'P';
      if (PatName[1] == 's') PatName[1] = 'S';
   }
/*------------------- arg 3:  Name Prefix for output files --------------*/
   strcpy(defname, PatName);
   if (argc > 3)
      strcpy(Prefix, argv[3]);
   else if (strcmp(PatName, "all") != 0)
   {
      printf("\nEnter name (1-7 characters) for this group of proteins [%s]: ",
	      defname);
      gets(Prefix);
   }
   if ((int) strlen(Prefix) < 2)
      strcpy(Prefix, defname);
/*------------------- arg 4:  Location of Swiss-Prot database -----------*/
   if (argc > 4)
      strcpy(swiss, argv[4]);
   else
   {
      printf("\nEnter name SWISS-PROT database file [none]:\n");
      gets(swiss);
   }
   if ((int) strlen(swiss) < 2)   strcpy(swiss, "none");
/*------------------- arg 5:  Directory for extracted proteins ----------*/
   getcwd(defname, 40);  strcat(defname, "/");
   if (strcmp(PatName, "all") != 0)
   {
      strcat(defname, Prefix); strcat(defname, "/");
   }
   if (argc > 5)
      strcpy(Pros, argv[5]);
   else if (strcmp(swiss, "none") != 0)
   {
      printf("\nEnter directory into which proteins extracted from %s",
	      swiss);
      printf("\n should be deposited [%s]: ", defname);
      gets(Pros);
   }
   if ((int) strlen(Pros) < 2)   strcpy(Pros, defname);
/*           Make sure directory ends with a backslash; expected by
		 uextract and motifj */
   if (Pros[(int) strlen(Pros)-1] != '/') strcat(Pros, "/");

/*---------- Append to repeats.dat file---------------------------------*/
   if ( (frep=fopen("repeats.dat", "a+t")) == NULL)
      printf("Cannot open repeats.dat\n");

/*--------------------- Load DB format information ---------------------*/
   init_dbs(dbs);

/*---------------------------------------------------------------------*/
   totpats = get_pat(dbs, fdat, frep);
   printf("\n%d patterns processed\n", totpats);
   fclose(fdat); fclose(frep);
   if ((strcmp(swiss, "none") != 0) && (strcmp(PatName, "all") != 0))
   {
      strcpy(defname, Prefix); strcat(defname,".lis");
      printf("\nStarting extract of sequences from %s...\n", swiss);
      execlp("uextract", "uextract", defname, swiss, NULL);
   }
   printf("\n");
   exit(0);
}  /*  end of main */
/*======================================================================*/
int get_pat(dbs, fdat, frep)
struct db_info *dbs[];
FILE *fdat, *frep;
{
   int nseq, npat, nswiss, i, db, goodpat, done;
   char line[MAXLINE], title[MAXLINE], id[25], *ptr;
   char pattern[3*MAXLINE], counts[3*MAXLINE], ac[10], filename[12];
   char repeats[3*MAXLINE];
   FILE *flis;
/* FILE *fpat;           Removed pattern file */

   nseq = npat = 0;
   pattern[0] = title[0] = id[0] = counts[0] = repeats[0] = '\0';
   done = NO;
   db = type_dbs(fdat, dbs);
   if (db >= 0 && db < MAXDB)
   {
      printf("\nLooking for pattern %s...", PatName);
      do
      {
	 if (strncmp(line, dbs[db]->start, (int) strlen(dbs[db]->start)) == 0)
	 {
	    nseq += 1;
	    id[0] = ac[0] = pattern[0] = counts[0] = repeats[0] = '\0';
	    ptr = strtok(&line[dbs[db]->title_offset], " ");
	    strcat(id, ptr);
	    ptr = strtok(NULL, ".");
	    if (strnjcmp(ptr, "PATTERN", 7) == 0 ||
                strnjcmp(ptr, "MATRIX", 6)  == 0)
	    {
	       while(fgets(line, MAXLINE, fdat) != NULL &&
		  strncmp(line, "DR", 2) != 0 &&
		  strncmp(line, dbs[db]->end, (int) strlen(dbs[db]->end)) != 0 )
	       {
		  if (strncmp(line, "AC", 2) == 0)
		  {
		     strcpy(ac, &line[dbs[db]->title_offset]);
		     ac[7] = '\0';
		     strcpy(title, ">");
		     strcat(title, ac);
		     strcat(title, " ;");	/* leave space for A,B...*/
		     strcat(title, id);
		  }
		  else if (strncmp(line, "DE", 2) == 0)
		  {
		     strcat(title, &line[dbs[db]->title_offset]);
		     ptr = strtok(title, "\r\n");
		     strcat(title, ";");
		  }
		  else if (strncmp(line, "DO", 2) == 0)
		  {
		     strcat(title, &line[dbs[db]->title_offset]);
		     ptr = strtok(title, "\r\n");
		     strcat(title, ";");
		  }
		  else if (strncmp(line, "PA", 2) == 0)
		  {
		     strcat(pattern, &line[dbs[db]->title_offset]);
		     ptr = strtok(pattern, "\r\n");
		  }
		  else if (strncmp(line, "NR", 2) == 0)
		  {
		     strcat(counts, &line[dbs[db]->title_offset]);
		     ptr = strtok(counts, "\r\n");
		  }
		  else if (strncmp(line, "CC", 2) == 0 &&
                           strstr(line, "/MAX-REPEAT") != NULL)
		  {
		     strcat(repeats, &line[dbs[db]->title_offset]);
		     ptr = strtok(repeats, "\r\n");
		  }
	       }   /*  now at DR or end */
/*    Fix up title;  change $ to %  */
	       for (i=0; i< (int) strlen(title); i++)
		  if (title[i] == '$') title[i] = '%';

/*    Check to see if we want this entry, & write it out if we do (we
      either look for one entry, or for all entries)  */
	       goodpat = NO;
	       if ( (strcmp(PatName, "all") == 0  ||
		     strcmp(PatName, ac) == 0) &&
		     (screen_pat(pattern, title, counts, repeats, frep) == YES) )
	       {
		  goodpat = YES;
		  npat++;
		  strcat(title, "$");    /*  finish off title */

	       }
	       if ( strcmp(PatName, "all") != 0 && goodpat) done = YES;
	       if (goodpat)
	       {
/*      Pattern file no longer produced - PATMAT can't use it
		  if ((strcmp(PatName, "all")) == 0)
		     strcpy(filename, ac);
		  else  strcpy(filename, Prefix);
		  strcat(filename, ".pat");
		  printf("\nCreating pattern file %s for pattern %s",
			 filename, ac);
		  if ( (fpat=fopen(filename, "w+t")) == NULL)
		  {
		     printf("\nCannot open file %s\n", filename);
		     return(-1);
		  }
		  fprintf(fpat, "%s\n", title);
		  fprintf(fpat, "%s\n", pattern);
		  fclose(fpat);
*/
/*------ Process the sequences, if there are any-------------------------*/
		  nswiss = 0;
		  if (strncmp(line, "DR", 2) == 0)
		  {
		     if ((strcmp(PatName, "all")) == 0)
			strcpy(filename, ac);
		     else  strcpy(filename, Prefix);
		     strcat(filename,".lis");
		     printf("\nCreating batch file %s for entry %s",
			     filename, ac);
		     if ( (flis=fopen(filename, "w+t")) == NULL)
		     {
			printf("\nCannot open file %s\n", filename);
			exit(-1);
		     }
		     fprintf(flis, "%s\n", title);    /* Title */
		     fprintf(flis, "%s", Pros);       /* Directory */
		     if (strcmp(PatName, "all") == 0)
			fprintf(flis,"%s/",ac);
		     fprintf(flis, "\n");
		     printf("\nCollecting sequences for %s...", ac);
		     nswiss += find_swiss(&line[dbs[db]->seq_offset], flis);
		     while(!feof(fdat) && fgets(line, MAXLINE, fdat) != NULL &&
		       strncmp(line, dbs[db]->end, strlen(dbs[db]->end)) != 0 &&
		       strncmp(line, "DR", 2) == 0)
		     {
			nswiss += find_swiss(&line[dbs[db]->seq_offset], flis);
		     }  /*  end of sequence */
		  }
		  fclose(flis);
		  printf("Found %d sequences", nswiss);
/*-----Not possible to execute uextract here & return;...-------*/
	       }  /*  end of if goodpat  */
	    }  /*  end of PATTERN */
	 }  /*  end of start of entry */
      }  while (!feof(fdat) && fgets(line, MAXLINE, fdat) != NULL && !done);
   }   /*  end of if db valid */
   printf("\n");

   return(npat);
}  /* end of get_pat */
/*===================================================================*/
/*  find_swiss                                                       */
/*     Finds XXXX..X$YYYY...Y in a string which can have up to 3     */
/*     sets of this form:  "P07507, MGP$BOVIN , T;"                  */
/*===================================================================*/
int find_swiss(line, flis)
char *line;
FILE *flis;
{
   int nswiss, len;
   char *ptr, swiss_key[20], temp[20];

   nswiss = 0;
   ptr = strtok(line, ",");			/* skip past "PO7507" */
   while (ptr != NULL)
   {
	 ptr = strtok(NULL, ",");	/*  now ptr = " MGP$BOVIN " */
	 if (ptr != NULL)
	 {
	       strcpy(temp, ptr);
	    /*  get rid of leading spaces   */
	       len = strspn(temp, " ");
	       strcpy(swiss_key, &temp[len]);
	    /*  get rid of trailing spaces  */
	       len = strcspn(swiss_key, " ");
	       swiss_key[len] = '\0';
	       ptr = strtok(NULL, ";");   /* now ptr = " T" */
	    /* get true positives (T), false negatives (N) and
		potential sequences (P)  */
	       if (ptr != NULL &&
		   (strcmp(ptr, " T") == 0 || strcmp(ptr, " N") == 0 ||
		    strcmp(ptr, " P") == 0) )
	       {
		  nswiss += 1;
		  fprintf(flis, "%-12s  PS=%s\n", swiss_key, ptr+1);
	       }
	       ptr = strtok(NULL, ",");     /*  now ptr = "P08493" */
	 }
   }
   return(nswiss);
}  /*  end of find_swiss */
/*=====================================================================*/
/*    Pattern must have at least one substitution, which means at
      least one [] or {}. It cannot have the words "signature 2" or
      "signature 3" anywhere in its description. It must have at least
      the value MINSEQS as true positives on the NR line; ie
      ... /POSITIVE=nn(mm) ...  mm >= MINSEQS                          */
/*=====================================================================*/
int screen_pat(pattern, title, counts, repeats, frep)
char *pattern;
char *title;
char *counts, *repeats;
FILE *frep;
{
   char *temp, *ptr, mem[40], tempc[MAXLINE];
   int npos, nhit, dups, runtype, distance, type, n, nl, d1, d2;
   int a1, a2, a3, nrep;

/*------ First see if we want this entry at all ---------------------*/
/*       Check for second or third signatures if all entries selected */
if (strcmp(PatName, "all") == 0)
{
   if ( strstr(title, "_2;") != NULL)         return(NO);
   if ( strstr(title, "_3;") != NULL)         return(NO);
   if ( strstr(title, "signature 2") != NULL) return(NO);
   if ( strstr(title, "signature 3") != NULL) return(NO);
}

/*       How many true positive sequences are there?                  */
   strcpy(tempc, counts);
   temp = strstr(tempc, "/POSITIVE=");
   if (temp == NULL)                          return(NO);
   nhit = npos = 0;
   ptr = strtok(temp, "=");
   ptr = strtok(NULL, "(");
   if (ptr != NULL)
   {
      nhit = atoi(ptr);
      ptr = strtok(NULL, ")");
      npos = atoi(ptr);
   }
/*       How many false negative sequences are there?                  */
   strcpy(tempc, counts);
   temp = strstr(tempc, "/FALSE_NEG=");
   if (temp != NULL)
   {
      ptr = strtok(temp, "=");
      if (ptr != NULL)
      {
	 nhit = nhit + atoi(ptr);
	 npos = npos + atoi(ptr);
      }
   }
/*       How many repeats are there?                  */
   strcpy(tempc, repeats);
   temp = strstr(tempc, "/MAX-REPEAT=");
   if (temp != NULL)
   {
      ptr = strtok(temp, "=");
      ptr = strtok(NULL, "(");
      if (ptr != NULL) nrep = atoi(ptr);
      if (frep != NULL) fprintf(frep, "%s %6d\n", Prefix, (nrep - 1) );
   }
   if (npos < MINSEQS)                        return(NO);

/*---------- Determine MOTIFJ parameters & add them to the title line ----*/
   runtype = 4;
   dups = nhit - npos; 
   printf("\nnhit = %d, npos = %d, dups = %d, nhit/npos = %f\n",
          nhit, npos, dups, (float) nhit / npos);
/*   dups = nrep - 1;         Why did we do this???? */
   distance = MAX_DISTANCE;
   d1 = d2 = nl = n = 0;
   a1 = a2 = a3 = NO;
/*   NOTE:  strtok doesn't work with the intervening call to classify! */
   while ( (n=strcspn(&pattern[nl], "-.")) > 0)
   {
      strncpy(&tempc[0], &pattern[nl], n);
      tempc[n] = '\0';
      nl += (n+1);
      type = classify(tempc);
      if (type < 0)			/* ambiguous distance; restart */
      {
	 a1 = a2 = a3 = NO; d1 = d2 = 0;
      }
      else if (type == 0)		/* amino acid */
      {
	 if (a1 == NO)      a1 = YES;
	 else if (a2 == NO) a2 = YES;
	 else if (a3 == NO) a3 = YES;
      }
      else				/* distance */
      {
	 if (a2 == YES)      d2 += type;
	 else if (a1 == YES) d1 += type;
      }
      if (a1 == YES && a2 == YES && a3 == YES)
      {
/*	 distance = UMIN(distance, UMAX(d1, d2)+1); */
	 if (d2 > d1) d2 = d1; d1++;
	 if (d1 < distance) distance = d1;
	 a1 = a2 = a3 = NO; d1 = d2 = 0;
      }
   }
   sprintf(mem, " MOTIFJ=[%d,%d,%d,%d];",
	   runtype, npos, dups, distance);
   strcat(title, mem);

   return(YES);
}  /*  end of screen_pat */
/*=====================================================================*/
/*   Classify pattern token into one of three categories:
       -1   ambiguous distance if x(m,n), [..](m,n), {..}(m,n)
	1+  distance if x, x(n), [..], [..](n), {..}, {..}(n)
	0   amino acid
========================================================================*/
int classify(token)
char *token;
{
   char *ptr, tempc[30];

   if (strstr(token, "x") != NULL ||
       strstr(token, "[") != NULL ||
       strstr(token, "{") != NULL )
   {
      if (strstr(token, ",") != NULL)         return(-1);   /* x(m,n) */
      if (strstr(token, "(") == NULL)         return(1);    /* just x */
      strcpy(tempc, token);
      ptr = strtok(tempc, "(");
      ptr = strtok(NULL, ")");
      return(atoi(ptr));
   }
   else return (0);					    /*  AA    */
}  /*  end of classify */
