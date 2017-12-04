//
//  NWalign.h
//  NWalign
//
//  Created by 包雨薇 on 17/11/19.
//  Copyright © 2017年 包雨薇. All rights reserved.
//

#ifndef NWalign_h
#define NWalign_h


#include <string>
#include <utility>
#include <deque>
#include <unordered_map>
using namespace std;


struct grid {
    char source = 'd';
    int myScore = 0;
    int wScore = 0;
    int nScore = 0;
    pair<int, int> wOrigin = make_pair(0, 0);
    pair<int, int> nOrigin = make_pair(0, 0);
};


int GO = -11;
int GX = -1;
int aligLn = 0;
int idenLn = 0;
bool printS = false;
bool printM = false;
bool quiet = false;
string S1;
string S2;
string MXFILE = "blosum62.mat";
unordered_map<string, int> SCORE;
vector<vector<grid>> myMap;


void getMode(int argc, char * argv[]);
void parseScoreMX();
void calculateMap();
void compareGrid(int i, int j);
void printMap();
void backTrace(deque<char> &seq1, deque<char> &seq2);
void printPath(const deque<char> &seq1, deque<char> &seq2);
#endif /* NWalign_h */
