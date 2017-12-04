/*    blclump.c  Compute Steve's sequence weights for blocks
           blclump <input blocks file> <output blocks file> <type> <scale>
		type = P  (position-based weights)
		     = V  (Voronoi weights)
		     = A  (Vingron & Argos weights)
                     = Cn (n% clustering)
 		scale = max (max == 100)
			n (sum == number of sequences)
			# (sum == #)
--------------------------------------------------------------------
 2/10/94 J. Henikoff (from blcons.c)
 2/21/94 Added Voronoi option
 3/ 4/94 Added cluster % option. Position-based columns are normalized.
 5/23/94 Added Vingron & Argos option.
 8/ 6/94 Added option to scale weights to sum to number of sequences
 4/ 5/95 Added call to free_block()
 1/ 6/96 Increased MAXWIDTH from 60 to 100 for Shmuel; still too inflexible.
12/ 5/96 Ignore X instead of treating it like another character in position().
 2/ 2/97 Increased MAXROWS from 400 to 800 to accomodate current blocks.
         Increased MAXWIDTH from 100 to 400 for Shmuel.
 2/ 2/97 Set MAXPWIDTH to 1500 for postion-based only; other methods 
         use MAXWIDTH = 100.
====================================================================*/

#define BLWEIGHT_C_
#define EXTERN
#define YES 1
#define NO 0
#define MAXNAME 80	/* Maximum file name length */
#define MAXAA 26	/* Dimension of aa_ arrays */
#define MAXWIDTH 100	/* Maximum block width */
#define MAXPWIDTH 1500	/* Maximum block width for postion-based method */
#define MAXROWS 800	/* Maximum number of rows in a block */
#define MAXINT 32000
/*   Already defined for Solaris */
#define RAND_MAX 2147483647		/* 2^31 - 1 */
#define NUMOFMUTANTS 10000		/* Voronoi randomizations */
#define ROUND(x) ((x >= 0.0) ? (int) (x+0.5) : (int) (x-0.5))
#define ORD(x) ((int) x - 65)
#define INDEX(n, col, row) (col*n - (col*(col+3))/2 - 1 + row)

#include "blocksprogs.h"

/*
 * Local variables and data structures
 */

void position();
void voronoi();
void noteaa();				/* used by voronoi() */
void calculate();			/* used by voronoi() */
void adjustweights();			/* used by voronoi() */
void cluster();
void argos();
void scale_weights();
void put_block();

struct vor {				/* used by voronoi() */
  int aapresent[MAXWIDTH][MAXAA];
  int aanumber[MAXWIDTH];
};

struct pair {				/* used by cluster() */
   int score, cluster;
};

/*=======================================================================*/
/*
 * main
 *   controls flow of program
 *   Parameters: argc, argv
 *   Error codes:
 */

int main(argc, argv)
     int argc;
     char *argv[];
{
  FILE *bfp, *ofp;
  Block *block;
  int i, wtype, clus, stype;
  char bdbname[MAXNAME], conname[MAXNAME], ctemp[10];

/* ------------1st arg = blocks database -------------------------------*/
   if (argc > 1)
      strcpy(bdbname, argv[1]);
   else
   {
      printf("\nEnter name of blocks database: ");
      gets(bdbname);
   }
   if ( (bfp=fopen(bdbname, "r")) == NULL)
   {
      printf("\nCannot open file %s\n", bdbname);
      exit(-1);
   }
/* ------------2nd arg = weighted blocks database ---------------------*/
   if (argc > 2)
      strcpy(conname, argv[2]);
   else
   {
      printf("\nEnter name of new weighted blocks database: ");
      gets(conname);
   }
   if ( (ofp=fopen(conname, "w")) == NULL)
   {
      printf("\nCannot open file %s\n", conname);
      exit(-1);
   }
/* ------------3rd arg = weighting type -------------------------------*/
   wtype = 0;			/* Position-based is default */
   clus = -1;
   ctemp[0] = '\0';
   if (argc > 3)
      strcpy(ctemp, argv[3]);
   else
   {
      printf("\nEnter weighting scheme (P=position-based, V=Voronoi, ");
      printf("\n A=Vingron & Argos, Cn = n% cluster [P]: ");
      gets(ctemp);
   }
   if (strlen(ctemp))
   {
      if (ctemp[0] == 'V' || ctemp[0] == 'v') wtype = 1;
      else if (ctemp[0] == 'A' || ctemp[0] == 'a') wtype = 3;
           else if (ctemp[0] == 'C' || ctemp[0] == 'c')
                {
                   wtype = 2;  clus = atoi(ctemp+1);
                   if (clus < 0) clus = 0;
                   if (clus > 100) clus = 100;
                }
   }
   if (wtype < 0 || wtype > 3) wtype = 0;
   switch (wtype)
   {
        case 0:
  	   printf("blweight: Calculating position-based weights\n");
	   break;
	case 1:
  	   printf("blweight: Calculating Voronoi weights\n");
	   break;
	case 2:
  	   printf("blweight: Calculating %d percent cluster weights\n", clus);
	   break;
	case 3:
  	   printf("blweight: Calculating Vingron & Argos weights\n");
	   break;
	default:
	   break;
   }
   /*---------------------Scale ----------------------------------*/
   stype = 0;
   if (argc > 4)
      strcpy(ctemp, argv[4]);
   else
   {
      printf("\nEnter scale (M=>max=100, N=>sum=#seq, ###=>sum=###) [M]: ");
      gets(ctemp);
   }
   if (strlen(ctemp))
   {
      if (ctemp[0] == 'm' || ctemp[0] == 'M') stype = 0;
      else if (ctemp[0] == 'n' || ctemp[0] == 'N') stype = 1;
           else 
           {
              stype = atoi(ctemp);
              if (stype < 100) stype = 100;
              if (stype > 999) stype = 999;
           }
   }
   switch (stype)
   {
        case 0:
  	   printf("blweight: Integer weights, maximum weight = 100\n");
	   break;
	case 1:
  	   printf("blweight: Decimal weights, sum = number of sequences\n");
	   break;
	default:
  	   printf("blweight: Integer weights, sum = %d\n", stype);
	   break;
   }
   /*-----------------------------------------------------------------*/

  while ((block = read_a_block(bfp)) != NULL)
  {
     switch (wtype)
     {
        case 0:
  	   position(block);
	   break;
	case 1:
	   voronoi(block);
	   break;
	case 2:
	   cluster(clus, block);
	   break;
	case 3:
	   argos(block);
	   break;
	default:
	   break;
     }
     scale_weights(block, stype);
     put_block(block, ofp, stype);
     free_block(block);
  }
   
  fclose(bfp); fclose(ofp);
  exit(0);

}  /* end of main */
/*=======================================================================
      Compute Steve's position-based sequence weights      
	Doesn't count - (0), X (23) or * (24), other chars are 1-22
========================================================================*/
void position(block)
Block *block;
{
   double diffaas[MAXPWIDTH], naas[MAXPWIDTH][MAXAA+1], dtemp;
   int seq, pos, aa, width;

   for (pos = 0; pos < MAXPWIDTH; pos++)
   { 
      diffaas[pos] = 0.0;
      for (aa = 0; aa < MAXAA+1; aa++)
         naas[pos][aa] = (double) 0.0;
   }
  
   width = block->sequences[0].length;
   if (width > MAXPWIDTH)
   {
        width = MAXPWIDTH;
	printf("ERROR: Block %s is too wide (%d), truncated to %d\n",
		block->number, block->sequences[0].length, width);
   }
   for (pos = 0; pos < width; pos++)
      for (seq = 0; seq < block->num_sequences; seq++)
      {
         if (block->residues[seq][pos] >= 1 &&
             block->residues[seq][pos] <= 22)
         {
            naas[pos][block->residues[seq][pos]] += 1;
         }
         else printf("Residue not counted: %d\n", block->residues[seq][pos]);
      }

   for (pos = 0; pos < width; pos++)
   { 
      for (aa = 1; aa <= 22; aa++)
      {
         if (naas[pos][aa] > 0.0)
         {
            diffaas[pos] += 1;	/* # of different types of aas in pos */
         }
      }
   }

   for (seq = 0; seq < block->num_sequences; seq++)
   {
      block->sequences[seq].weight = 0.0;
      for (pos = 0; pos < width; pos++)
      {
         aa = block->residues[seq][pos];
         dtemp = diffaas[pos] * naas[pos][aa];
         if (dtemp > 0.0)
           block->sequences[seq].weight += 1.0 / dtemp;
      } 
   }
}  /* end of position */
/*=========================================================================
    For each pair of segments in the block, count the number of positions
    that differ.
===========================================================================*/
void argos(block)
Block *block;
{
   int width, pos, seq, seq1, seqs, diff[MAXROWS][MAXROWS];
   double total;

   width = block->sequences[0].length;
   if (width > MAXWIDTH)
   {
        width = MAXWIDTH;
	printf("ERROR: Block %s is too wide (%d), truncated to %d\n",
		block->number, block->sequences[0].length, width);
   }
   seqs = block->num_sequences;
   if (seqs > MAXROWS)
   {
      seqs = MAXROWS;
	printf("ERROR: Block %s is too deep (%d), truncated to %d sequences\n",
		block->number, block->num_sequences, seqs);
   }
   for (seq = 0; seq < seqs; seq++)
      for (seq1 = 0; seq1 < seqs; seq1++)
         diff[seq][seq1] = 0;
   for (seq = 0; seq < seqs; seq++)
      for (seq1 = 0; seq1 < seqs; seq1++)
         for (pos = 0; pos < width; pos++)
            if (block->residues[seq][pos] !=
                block->residues[seq1][pos])
                    diff[seq][seq1] = diff[seq][seq1] + 1;
 
   total = 0.0;
   for (seq = 0; seq < seqs; seq++)
      block->sequences[seq].weight = 0.0;
   for (seq = 0; seq < seqs; seq++)
      for (seq1 = 0; seq1 < seqs; seq1++)
      {
         block->sequences[seq].weight += (double) diff[seq][seq1];
         total += (double) diff[seq][seq1];
      }
   for (seq = 0; seq < seqs; seq++)
   {
      block->sequences[seq].weight /= total;
   }
   
}   /* end of argos */
/*===================================================================
======================================================================*/
void voronoi(block)
Block *block;
{
  struct vor *v;
  int irow, seqs;
  double owned[MAXROWS];

  v = (struct vor *) malloc(sizeof(struct vor));
  noteaa(block, v);
  seqs = block->num_sequences;
  if (seqs > MAXROWS)
  {
     seqs = MAXROWS;
     printf("ERROR: Block %s is too deep (%d), truncated to %d sequences\n",
	block->number, block->num_sequences, seqs);
  }
  for (irow = 0; irow < seqs; irow++)
      owned[irow]  =  0.0; /* set # owned to zero */
  calculate(block, v, owned);
  for (irow = 0; irow < seqs; irow++)
     block->sequences[irow].weight = owned[irow];
}  /* end of voronoi */
/*===================================================================
 examine the alignment and store the amino acids that occur at
 each col in the array aapresent. Note how many different aa  
 occur at each col and store in the array aanumber.            
======================================================================*/
void noteaa(block, v)
Block *block;
struct vor *v;
{ /* noteaa */
 int acol, arow;
 int alreadypresent;      /* used to ensure that only one occurrence */
                          /* of each aa goes into aapresent          */
 int aap;                /* used to count in the array aanumber */
 int cols, rows;

 cols = block->sequences[0].length;
 if (cols > MAXWIDTH)
 {
     cols = MAXWIDTH;
     printf("ERROR: Block %s is too wide (%d), truncated to %d\n",
		block->number, block->sequences[0].length, cols);
 }
 rows = block->num_sequences;
 if (rows > MAXROWS)
 {
     rows = MAXROWS;
     printf("ERROR: Block %s is too deep (%d), truncated to %d sequences\n",
	block->number, block->num_sequences, rows);
 }

 for (acol = 0; acol < cols; acol++) {
  v->aapresent[acol][0] = block->residues[0][acol] ;
  v->aanumber[acol]  =  1;
 }
 for (acol = 0; acol < cols; acol++) {
  for (arow = 1; arow < rows; arow++) {
   alreadypresent  =  NO;
   for (aap = 0; aap < v->aanumber[acol]; aap++)
     if (block->residues[arow][acol] == v->aapresent[acol][aap])
       alreadypresent  =  YES;
   if (!alreadypresent) {
    v->aapresent[acol][v->aanumber[acol]] = block->residues[arow][acol];
    v->aanumber[acol] += 1; /* another aa occurs */
   } /* if */
  } /* for arow */
  v->aapresent[acol][v->aanumber[acol]] = '\0';
 } /* for acol */
}   /* end of noteaa */
/*===================================================================
 for every mutant needed,  generate the mutant, find its distance 
 to the other sequences in the alignment,  keeping short distance
 and the identities of those sequences that lie closest. Then   
 adjust weights.                                               
======================================================================*/
void calculate(block, v, owned)
Block *block;
struct vor *v;
double owned[MAXROWS];
{
 double mutant[MAXWIDTH][MAXAA];
 int icol, irow, mut, z;
 double total, value, rrand;
 int cols, rows;

 cols = block->sequences[0].length;
 if (cols > MAXWIDTH)
 {
     cols = MAXWIDTH;
     printf("ERROR: Block %s is too wide (%d), truncated to %d\n",
		block->number, block->sequences[0].length, cols);
 }
 rows = block->num_sequences;
 if (rows > MAXROWS)
 {
   rows = MAXROWS;
   printf("ERROR: Block %s is too deep (%d), truncated to %d rows\n",
		block->number, block->num_sequences, rows);
 }

 total = 0.0;
 for (mut = 0; mut < NUMOFMUTANTS; mut++) {
  for (icol = 0; icol < cols; icol++) { /* make the mutant */
   value = 0;
   for (z = 0; z < v->aanumber[icol]; z++) {
    rrand = (double) ran0() / (RAND_MAX + 1.0);  /* between 0.0 and 1.0 */
    if (rrand != 0.0) value = -log(rrand);
    else value = 0.0;
    mutant[icol][v->aapresent[icol][z]] = value;
    total = total + value;
   }
   for (z = 0; z < v->aanumber[icol]; z++)
    mutant[icol][v->aapresent[icol][z]] =
      mutant[icol][v->aapresent[icol][z]] / total;
  } /* for icol */ /* that finishes making the mutant */
    adjustweights(block, mutant, owned);
 } /* for mut,  i.e. for every mutant needed */
}  /*  end of calculate */
/*===================================================================
======================================================================*/
void adjustweights(block, mutant, owned)
Block *block;
double mutant[MAXWIDTH][MAXAA], owned[MAXROWS];
{
 double distance, shortdistance;
 int jrow, jcol, rowset[MAXROWS], r;
 int numinset; 			/* number of elements in the set */
 int cols, rows;

 cols = block->sequences[0].length;
 if (cols > MAXWIDTH)
 {
     cols = MAXWIDTH;
     printf("ERROR: Block %s is too wide (%d), truncated to %d\n",
		block->number, block->sequences[0].length, cols);
 }
 rows = block->num_sequences;
 if (rows > MAXROWS)
 {
   rows = MAXROWS;
   printf("ERROR: Block %s is too deep (%d), truncated to %d rows\n",
		block->number, block->num_sequences, rows);
 }

 shortdistance  =  MAXINT;
 for (r = 0; r < rows; r++) rowset[r] = NO;
 for (jrow = 0; jrow < rows; jrow++) {
  distance  =  0;
  /* identity distances are used */
  for (jcol = 0; jcol < cols; jcol++)
   distance = distance + mutant[jcol][block->residues[jrow][jcol]];
  if (distance <= shortdistance)
  {
   if (distance < shortdistance)
   {
     for (r = 0; r < rows; r++) rowset[r] = NO;
     rowset[jrow]  =  YES;		/* initialize rowset here */
     shortdistance  =  distance;
     numinset  =  1;
   }
   else
   {
     rowset[jrow] = YES;		/* add row to rowset here */
     numinset  =  numinset + 1;
   }
  } /* if distance <= */
 } /* for jrow */ /* found the shortest distance and those rows in the */
                    /* alignment that have that distance. now weights.  */
 /* calculate a single weight */
 for (jrow = 0; jrow < rows; jrow++)
  if (rowset[jrow])
     owned[jrow] += ((double) 1.0 / numinset);
}  /* end of adjustweights */
/*======================================================================*/
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

/*=========================================================================
===========================================================================*/
/*======================================================================*/
/*    Cluster sequences in a block based on the number of               */
/*    identities within the block. Sets Block.cluster & Block.ncluster  */
/*      1. Compute number of identities for each possible pair of seqs. */
/*         Results stored in lower half of matrix (pairs).              */
/*      2. Use clustering threshold % of # of AAs in trimmed block.     */
/*      3. Cluster recursively by traversing cols, rows of matrix.      */
/*UNIX NOTE:  Program aborts when running under UNIX at free(pairs),
   so use the fixed size declaration pairs & remove the malloc() &
   free() calls when compiling for UNIX                                 */
/*======================================================================*/
void cluster(clus, block)
int clus;
Block *block;
{
   int iclus, npair, threshold, s1, s2, l1, l2, px, i, i1, i2;
   int nclus[MAXROWS], minclus, oldclus, width, nseq;
   struct pair pairs[MAXROWS*(MAXROWS-1)/2]; 
   int icluster[MAXROWS];

   width = block->sequences[0].length;
/*
   if (width > MAXWIDTH)
   {
        width = MAXWIDTH;
	printf("ERROR: Block %s is too wide (%d), truncated to %d\n",
		block->number, block->sequences[0].length, width);
   }
*/
   nseq = block->num_sequences;
   if (nseq > MAXROWS)
   {
     nseq = MAXROWS;
     printf("ERROR: Block %s is too deep (%d), truncated to %d nseq\n",
		block->number, block->num_sequences, nseq);
   }
   npair = nseq*(nseq-1)/2;
   threshold = (int) (clus*(width))/100;

/*    Compute scores for all possible pairs of sequences            */
   for (s1=0; s1<nseq-1; s1++)   		/* col = 0, n-2     */
   {
      l1 = 0;
      for (s2=s1+1; s2<nseq; s2++)	/* row = col+1, n-1 */
      {
	 l2 = 0;
	 px = INDEX(nseq, s1, s2);
	 pairs[px].score = 0;
	 pairs[px].cluster = -1;
	 for (i=0; i<=width; i++)
	 {
	    i1 = l1+i;  i2 = l2+i;
	    if (i1 >= 0 && i1 < width &&
		i2 >= 0 && i2 < width &&
		block->residues[s1][i1] == block->residues[s2][i2])
		   pairs[px].score += 1;
	 }
      }  /* end of s2 */
   }  /* end of s1 */

/*  Print scores */
/*   printf("\nThreshold=%d", threshold);
   for (s2=1; s2<nseq; s2++)
   {
      printf ("\n");
      for (s1=0; s1<s2; s1++)
      {
	 px = INDEX(nseq, s1, s2);
	 printf(" %.3d", pairs[px].score);
      }
    }
*/

/*-------Cluster if score exceeds threshold by scanning cols (s1) */
   for (s1=0; s1<nseq; s1++)
   {
      icluster[s1] = -1;			/* clear out old values */
      nclus[s1] = 0;
   }
   iclus = 0;        				/* cluster number */
   for (s1=0; s1<nseq-1; s1++)   		/* col = 0, n-2     */
      for (s2=s1+1; s2<nseq; s2++)	/* row = col+1, n-1 */
      {
	 px = INDEX(nseq, s1, s2);
	 if (pairs[px].score >= threshold)	/*  cluster this pair */
	 {
	    if (icluster[s1] < 0)          /* s1 not yet clustered */
	    {
	       if (icluster[s2] < 0)       /* new cluster */
	       {
		  icluster[s1] = iclus++;
		  icluster[s2] = icluster[s1];
	       }
	       else  				/* use s2's cluster  */
		  icluster[s1] =  icluster[s2];
	    }
	    /*  use s1's cluster if it has one and s2 doesn't */
	    else if (icluster[s1] >= 0 && icluster[s2] < 0)
	       icluster[s2] = icluster[s1];
	    /* merge the two clusters into the lower number */
	    else if (icluster[s1] >= 0 && icluster[s2] >= 0)
	    {
	       minclus = icluster[s1]; oldclus = icluster[s2];
	       if (icluster[s2] < icluster[s1])
	       {
		  minclus = icluster[s2]; oldclus = icluster[s1];
	       }
	       for (i1=0; i1<nseq; i1++)
		 if (icluster[i1] == oldclus)
		     icluster[i1] = minclus;
	    }
	 }  /* end of if pairs */
      }  /* end of s2 */

   /*---  Set ncluster, get rid of negative cluster numbers --*/
   for (s1=0; s1<nseq; s1++)
   {
      if (icluster[s1] < 0) 
	  icluster[s1] = iclus++;
   }
   for (s1=0; s1<nseq; s1++)
	  nclus[icluster[s1]] += 1;

   /*----------  Now compute the weight for each sequence ------------*/
   for (s1 = 0; s1 < nseq; s1++)
   {
      block->sequences[s1].weight = (double) 1.0 / nclus[icluster[s1]];
   }
}  /* end of cluster */
/*======================================================================*/
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
 * put_block
     Same as blimps output_block(), but outputs weight in floating
     point if stype == 1 => weights sum to # of sequences
 *   Outputs a block data structure to the given file. 
 *   Parameters: 
 *     Block *block:  the block to print
 *     FILE  *obfp:   the ouput block file pointer
 ======================================================================*/

void put_block(block, obfp, stype)
     Block *block;
     FILE *obfp;
     int stype;
{
  int i,j,k;

  fprintf(obfp, "ID   %s\n", block->id);
  fprintf(obfp, "AC   %s\n", block->ac);
  fprintf(obfp, "DE   %s\n", block->de);
  fprintf(obfp, "BL   %s", block->bl);

  /* print clusters */
  for (i=0; i<block->num_clusters; i++) {
    fprintf(obfp, "\n");
    for (j=0; j<block->clusters[i].num_sequences; j++) {
      fprintf(obfp, "%10.10s (%4d) ", 
	      block->clusters[i].sequences[j].name,
	      block->clusters[i].sequences[j].position);
      for(k=0; k<block->clusters[i].sequences[j].length; k++) {
	fprintf(obfp, "%c", 
		aa_btoa[block->clusters[i].sequences[j].sequence[k]]);
      }
      if (stype == 1)
      fprintf(obfp, " %7.3f\n", block->clusters[i].sequences[j].weight);
      else
      fprintf(obfp, " %3d\n", round(block->clusters[i].sequences[j].weight));
    } /* end of cluster */
  }  
  /* end of the block */
  fprintf(obfp, "//\n");
  fflush(obfp);
} /* end of put_block */
