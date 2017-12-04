/*    blalign.c    Print blocks as multiple alignment      
        blalign <blocks database> -[p|s|f]
         <blocks database> = name of blocks db file, aligns each group
         p|s|f = p for Posfai style, s for short style, f for fasta style
--------------------------------------------------------------------
 4/21/94 J. Henikoff (from blcovar.c)
 4/26/94 Added option for "Posfai" or abbreviated style of output
 6/23/94 Print number of sequences and block widths
 8/12/96 Added fasta style output
 2/19/97 Added output width parameter = OutWidth
====================================================================*/

#define BLALIGN_C_
#define EXTERN
#define MAXNAME 80	/* Maximum file name length */
#define MAXAA 25	/* Dimension of aa_ arrays */
#define MAXWIDTH 60	/* Maximum block width */
#define OUTWIDTH 66	/* Default output line */
#define OUTMAX 240	/* Maximum output line width */
#define MAXLEN 4000	/* Maximum sequence length */
#define FLZERO 0.0000001
#define INDEX(n, col, row) (col*n + row - (col*(col+1))/2)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "blocksprogs.h"


/*
 * Local variables and data structures
 */

/*  first entry has no block, just nblock & totwidth, other entries
    in list have just pointers to the blocks */
struct block_list {		/* list of blocks for a family */
   int nblock;				/* number of blocks in family */
   int nseq;				/* number of sequences in blocks */
   int totwidth;			/* total width of the blocks */
   Block *block;
   struct block_list *next;
};
struct seqseq {			/* reorder sequences */
   int seq;
   int pos;
};
struct align {			/* alignment line */
   char line[MAXLEN];
};

struct block_list *make_blist();
void insert_blist();
void free_blist();
void order_seq();
int print_align();
void format_seqs();
void posfai_seqs();
void short_seqs();
void fasta_seqs();

int OutWidth;			/* width of output */

/*=======================================================================
 * main
 =======================================================================*/

int main(argc, argv)
     int argc;
     char *argv[];
{
  FILE *bfp;
  Block *block;
  struct block_list *blist;
  int i, np, style;
  char bdbname[MAXNAME], conname[MAXNAME], save_number[10], ctemp[10];

   if (argc < 2)
   {
      printf("BLALIGN: Copyright 1997 Fred Hutchinson Cancer Research ");
      printf("Center\nPrints blocks as a multiple alignment.\n");
      printf("USAGE:  blalign <file of blocks> [-p|s|f] [width]\n");
      printf("          p=Posfai style, s=short style, f=fasta style\n");
      printf("          defaults are short style, width = 66\n");
   }

/* ------------1st arg = blocks database -------------------------------*/
   if (argc > 1)
      strcpy(bdbname, argv[1]);
   else
   {
      printf("\nEnter name of file of blocks: ");
      gets(bdbname);
   }
   if ( (bfp=fopen(bdbname, "r")) == NULL)
   {
      printf("\nCannot open file %s\n", bdbname);
      exit(-1);
   }
/*----------2nd arg = p or s -------------------------------------------*/
   style = 0; 			/* short*/
   if (argc > 2 && argv[2][0] == '-')
   {
	if      (argv[2][1] == 'p' || argv[2][1] == 'P') style = 1;
	else if (argv[2][1] == 'f' || argv[2][1] == 'F') style = 2;
   }
   
/*-----------3rd arg = output width -----------------------------------*/
   OutWidth = OUTWIDTH;
   if (argc > 3) OutWidth = atoi(argv[3]);
   if (OutWidth < 1 || OutWidth > OUTMAX) OutWidth = OUTWIDTH;

/*----------------------------------------------------------------------*/
  save_number[0] = '\0';
  blist = NULL;
  while ((block = read_a_block(bfp)) != NULL)
  {
     if (strncmp(save_number, block->number, 7) != 0)
     {       /*  new family  */
        /*  process last family */
        if (blist != NULL && blist->nblock > 0)
        {
           np = print_align(blist, style);
           free_blist(blist);
        }
        strncpy(save_number, block->number, 7); save_number[7] = '\0';
        if (style < 2) printf(">%s %s\n", save_number, block->de);
        blist = make_blist();
        insert_blist(blist, block);
     }
     else    
     {     /*  same family */
        insert_blist(blist, block);
     }
  }
  /*  process final family */
  if (blist != NULL && blist->nblock > 0)
  {
     np = print_align(blist, style);
     free_blist(blist);
  }

  fclose(bfp);
  exit(0);

}  /* end of main */
/*=======================================================================
     routines for a list of blocks
========================================================================*/
struct block_list *make_blist()
{
   struct block_list *new;
   
   new = (struct block_list *) malloc (sizeof(struct block_list));
   new->nblock = new->nseq = new->totwidth = 0;
   new->block = NULL;
   new->next = NULL;

   return(new);
}  /* end of make_blist */

void insert_blist(blist, block)
struct block_list *blist;
Block *block;
{
   struct block_list *cur;

   cur = blist;
   while (cur->next != NULL)
      cur = cur->next;
   cur->next = make_blist();
   cur->next->block = block;
   blist->nblock += 1;
   blist->nseq = block->num_sequences;
   blist->totwidth += block->sequences[0].length;
   
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
      free(cur);
      last->next = NULL;
      free_blist(last);
   }
   else free(blist);

}  /* end of free_blist */
/*=======================================================================
	blist = list of blocks for a group
========================================================================*/
int print_align(blist, style)
struct block_list *blist;
int style;
{
   int i, np, seq;
   int nseq;
   struct block_list *bcur, *bfirst;
   struct seqseq *sseq;

   /*     Initialize */
   nseq = blist->nseq;
   sseq = (struct seqseq *) malloc (nseq * sizeof(struct seqseq));
  
   /*   Sequences aren't in the same order in different blocks */
   /*   They are printed in the order of the first block       */

   /*   Do the first block                                     */
   bfirst = blist->next;
   if (bfirst->block == NULL) return(0);

   for (seq = 0; seq < nseq; seq++)
   {
      sseq[seq].seq = seq;
      sseq[seq].pos = 0;
   }
   format_seqs(sseq, bfirst->block);
   np = 1;

   bcur = bfirst->next;
   while (bcur != NULL && bcur->block != NULL)
   {
      /*  sequences may not be in the same order */
      order_seq(sseq, bfirst->block, bcur->block);
      format_seqs(sseq, bcur->block);
      np++;
      bcur = bcur->next;
   } /* end of block */

   if (style < 2)
      printf("%d sequences are included in %d blocks\n", 
           blist->nseq, blist->nblock);
   if (style == 1)      posfai_seqs(sseq, blist);
   else if (style == 2) fasta_seqs(sseq, blist);
   else                 short_seqs(sseq, blist);
   return(np);

}  /* end of print_align */

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
/*=======================================================================
        Formats one block: determines how many gaps and unaligned
        positions are required to line up the block in each sequence.

	These values are stored as sequences[].type (gaps) and
        sequences[].max_length (unaligned).
========================================================================*/
void format_seqs(sseq, b)
struct seqseq *sseq;
Block *b;
{
   int nseq, seq, seq1, pos, posgap, posun;
   int between;

   nseq = b->num_sequences;
   /*    find maximum distance since last block in all sequences */
   between = 0;
   for (seq = 0; seq < nseq ; seq++)
   {
      seq1 = sseq[seq].seq;
      posun = b->sequences[seq1].position - sseq[seq].pos;
      if (posun > between) between = posun;
   }
   for (seq=0; seq < nseq; seq++)
   {
      seq1 = sseq[seq].seq;
      /*   unaligned region between blocks */
      posun = b->sequences[seq1].position - sseq[seq].pos;
      /*   gap to accomodate longest unaligned region */
      posgap = between - posun;
      sseq[seq].pos = b->sequences[seq1].position + b->sequence_length;
      /*    Reusing these fields !!! */
      b->sequences[seq1].type = posgap;
      b->sequences[seq1].max_length = posun;
   }
}   /* end of format_seqs */
/*==================================================================
      Prints the alignment in the style of Posfai (every position)
====================================================================*/
void posfai_seqs(sseq, blist)
struct seqseq *sseq;
struct block_list *blist;
{
   int alen, nseq, seq, seq1, pos, posgap, posun, curpos;
   char ctemp[OUTMAX];
   struct block_list *bcur;
   Block *b;
   struct align *out;

   nseq = blist->nseq;
   out = (struct align *) malloc(nseq * sizeof(struct align));
   if (out == NULL)
   {   printf("\nposfai_seqs: OUT OF MEMORY"); exit(-1); }
   for (seq=0; seq < nseq; seq++) out[seq].line[0] = '\0';

   bcur = blist->next;
   while (bcur != NULL && bcur->block != NULL)
   {
      b = bcur->block;
      for (seq=0; seq < nseq; seq++)
      {
         order_seq(sseq, blist->next->block, bcur->block);
         seq1 = sseq[seq].seq;
         /*   unaligned region between blocks */
         posun = b->sequences[seq1].max_length;
         /*   gap to accomodate longest unaligned region */
         posgap = b->sequences[seq1].type;
         for (pos=0; pos < posgap; pos++)
            strcat(out[seq].line, " ");
         for (pos=0; pos < posun; pos++)
            strcat(out[seq].line, ".");
         sprintf(ctemp, "%5d ", b->sequences[seq1].position);
         strcat(out[seq].line, ctemp);
         for (pos=0; pos < b->sequence_length; pos++)
            ctemp[pos] = aa_btoa[b->residues[seq1][pos]];
         ctemp[b->sequence_length] = '\0';
         strcat(out[seq].line, ctemp);
      }
      bcur = bcur->next;
   }

   alen = strlen(out[0].line);
   b = blist->next->block;
   for (curpos = 0; curpos < alen; curpos += OutWidth)
   {
      for (seq=0; seq < nseq; seq++)
      {
         printf("\n%10s ", b->sequences[seq].name);
         for (pos=curpos; pos < curpos+OutWidth; pos++)
            if (pos < alen)   printf("%c", out[seq].line[pos]);
      }
      printf("\n");
   }
   printf("\n");
}   /* end of posfai_seqs */
/*==================================================================
      Prints the alignment in short style (doesn't print spaces
      between blocks)
====================================================================*/
void short_seqs(sseq, blist)
struct seqseq *sseq;
struct block_list *blist;
{
   int alen, nseq, seq, seq1, pos, posgap, posun, curpos;
   char ctemp[OUTMAX], header[OUTMAX];
   struct block_list *bcur;
   Block *b;
   struct align *out;

   nseq = blist->nseq;
   out = (struct align *) malloc(nseq * sizeof(struct align));
   if (out == NULL)
   {   printf("\nshort_seqs: OUT OF MEMORY"); exit(-1); }
   for (seq=0; seq < nseq; seq++) out[seq].line[0] = '\0';
   header[0] = '\0';

   curpos = 11;   /* seq name width */
   bcur = blist->next;
   while (bcur != NULL && bcur->block != NULL)
   {
      b = bcur->block;

      /*  Page width required to print each block is
              name=11 + distance from last block=8 + offset=6 +
              max block width for blockmaker=55 => 80 max  */
      /*   Print full line */
      if ((curpos + 14 + b->sequence_length) > OutWidth)
      {
         printf("\n            %s", header);
         header[0] = '\0';
         for (seq=0; seq < nseq; seq++)
         {
            printf("\n%10s ", blist->next->block->sequences[seq].name);
            printf("%s", out[seq].line);
            out[seq].line[0] = '\0';
         }
         printf("\n");
         curpos = 11;
      }
      curpos += (14 + b->sequence_length);

      /*    Format the block  */
      strcat(header, b->number);
      sprintf(ctemp, ", width = %d", b->sequence_length);
      strcat(header, ctemp);
      for (seq=0; seq < nseq; seq++)
      {
         order_seq(sseq, blist->next->block, bcur->block);
         seq1 = sseq[seq].seq;
         /*   unaligned region between blocks */
         posun = b->sequences[seq1].max_length;
         /*  print distance bewteen blocks for this sequence after 1st block*/
         if (bcur != blist->next)
         {
            sprintf(ctemp, " (%4d) ", posun);
            strcat(out[seq].line, ctemp);
         }
         sprintf(ctemp, "%5d ", b->sequences[seq1].position);
         strcat(out[seq].line, ctemp);
         for (pos=0; pos < b->sequence_length; pos++)
            ctemp[pos] = aa_btoa[b->residues[seq1][pos]];
         ctemp[b->sequence_length] = '\0';
         strcat(out[seq].line, ctemp);
      }
      /*   Pad out the header line */
      alen = strlen(header);
      for (pos = alen; pos < strlen(out[0].line); pos++)
         strcat(header, " ");
      bcur = bcur->next;
   }
   /*   Print last block */
   printf("\n            %s", header);
   for (seq=0; seq < nseq; seq++)
   {
      printf("\n%10s ", blist->next->block->sequences[seq].name);
      printf("%s", out[seq].line);
   }

   printf("\n");
}   /* end of short_seqs */
/*==================================================================
      Prints the portion of the sequences in blocks in fasta format
====================================================================*/
void fasta_seqs(sseq, blist)
struct seqseq *sseq;
struct block_list *blist;
{
   int alen, nseq, seq, seq1, pos, posgap, posun, curpos;
   char ctemp[OUTMAX], header[OUTMAX];
   struct block_list *bcur;
   Block *b;
   struct align *out;

   nseq = blist->nseq;
   out = (struct align *) malloc(nseq * sizeof(struct align));
   if (out == NULL)
   {   printf("\nfasta_seqs: OUT OF MEMORY"); exit(-1); }
   for (seq=0; seq < nseq; seq++) out[seq].line[0] = '\0';
   header[0] = '\0';

   for (seq=0; seq < nseq; seq++)
   {
      printf(">%10s from blocks\n", blist->next->block->sequences[seq].name);
      curpos = 0;
      bcur = blist->next;
      while (bcur != NULL && bcur->block != NULL)
      {
         b = bcur->block;
         order_seq(sseq, blist->next->block, b);
         seq1 = sseq[seq].seq;
         printf("%5d ", b->sequences[seq1].position);
         for (pos=0; pos < b->sequence_length; pos++)
            printf("%c", aa_btoa[b->residues[seq1][pos]]);
         printf("\n");
         bcur = bcur->next;
      } /* end of blocks */
   }   /*  end of a sequence  */
}   /* end of fasta_seqs */
