# Recalculation of BLOSUM using modern sequence data

Files in this Repository:
+ **1234.txt** results from homologsearch.py for seed 1234
+ **5678.txt** results from homologsearch.py for seed 5678
+ **9001.txt** results from homologsearch.py for seed 9001
+ **astral-scopedom-seqres-gd-sel-gs-bib-49-2.06.fa** sequences from the Astral database for testing
+ **blosum62.csv** results for blosum62
+ **blosum62.mat** the BLOSUM62 matrix used in NWalign
+ **homologsearch.py** runs NW comparisons and generates confusion matrices for each substitution matrix
+ **iij2mat.py** converts files output by blosum.c into format readable by NWalign
+ **makefile** makefile for NWalign.cpp
+ **new99.csv** results for new99
+ **new99.mat** our new generated matrix used in NWalign
+ **NWalign.cpp** Needleman-Wunsch alignment implemented in C++ for efficiency
+ **NWalign.h** header file for NWalign.cpp
+ **NWalign.py** original NW implementation from which NWalign.cpp was made
+ **randblos.csv** results for randblos
+ **randblos.mat** randomly generated substitution matrix based off of BLOSUM62
+ **randblos.py** randomly generates a subsitution matrix given value frequencies present in BLOSUM62
+ **README.md** this file
+ **seeds.txt** the three seeds used to run triplicate
+ **tests.R** R script that generates mean accuracy, precision, and recall and performs t-tests of all these parameters between the three matrices