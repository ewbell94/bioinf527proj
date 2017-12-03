//
//  main.cpp
//  NWalign
//
//  Created by 包雨薇 on 17/11/19.
//  Copyright © 2017年 包雨薇. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <deque>
#include <vector>
#include <utility>
#include <limits>
#include <getopt.h>
#include <assert.h>
#include "NWalign.h"
using namespace std;

int main(int argc, char * argv[]) {

    // parse command line
    getMode(argc, argv);

    // parse scoring matrix
    parseScoreMX();

    // calculate the matrix
    calculateMap();
    if (printM && !quiet) {
        printMap();
    }
    

    // back trace
    deque<char> seq1;
    deque<char> seq2;
    backTrace(seq1, seq2);
    
    // output
    if (printS && !quiet) {
        printPath(seq1, seq2);
    }
    
    
    
    return 0;
}


void getMode(int argc, char * argv[]) {
    opterr = true;
    int choice;
    int option_index = 0;
    option long_options[] = {
        { "seq1", required_argument, nullptr, 'p'},
        { "seq2", required_argument, nullptr, 'q'},
        { "mx", optional_argument, nullptr, 'm'},
        { "go", optional_argument, nullptr, 'O'},
        { "gx", optional_argument, nullptr, 'X'},
        { "printS", no_argument, nullptr, 'S'},
        { "pringM", no_argument, nullptr, 'M'},
        { "quiet", no_argument, nullptr, 'Q'},
        { "help", no_argument, nullptr, 'h'},
        { nullptr, 0, nullptr, '\0'}
    };

    while ((choice = getopt_long(argc, argv, "p:q:m:O:X:SMQh",
                                 long_options, &option_index)) != -1) {
        switch(choice) {
            case 'p':{
//                string s1file = optarg;
//                ifstream infile1;
//                infile1.open(s1file);
//                if (!infile1.is_open()) {
//                    cout << "error opening infile" << endl;
//                    exit(1);
//                }
//                string inLine1 = "";
//                while (getline(infile1, inLine1)) {
//                    S1 += inLine1;
//                }
//                infile1.close();
                S1 = optarg;
                break;
            }
            case 'q':{
//                string s2file = optarg;
//                ifstream infile2;
//                infile2.open(s2file);
//                if (!infile2.is_open()) {
//                    cout << "error opening infile" << endl;
//                    exit(1);
//                }
//                string inLine2 = "";
//                while (getline(infile2, inLine2)) {
//                    S2 += inLine2;
//                }
//                infile2.close();
                S2 = optarg;
                break;
            }
            case 'm':{
                string a = optarg;
                MXFILE = a;
                break;
            }
            case 'O':{
                GO = stoi(optarg);
                break;
            }
            case 'X':{
                GX = stoi(optarg);
                break;
            }
            case 'S':{
                printS = true;
                break;
            }
            case 'M':{
                printM = true;
                break;
            }
            case 'Q':{
                quiet = true;
                break;
            }
            case 'h':{
                cout << "seq1/p: Sequence 1(just a sequence!NO FILE) " << endl;
                cout << "seq2/q: Sequence 2(just a sequence!NO FILE) " << endl;
                cout << "mx/m: Scoring matrix, default blosum62.mat " << endl;
                cout << "go/O: Gap opening penalty, default -11 " << endl;
                cout << "gx/X: Gap extension penalty, default -1 " << endl;
                cout << "printS/S: Print the aligned sequences " << endl;
                cout << "printM/M: Print the calculated matrix " << endl;
                cout << "quiet/Q: Print the sequence id only " << endl;
                exit(0);
            }
            default:{
                cerr << "Error: invalid option " << '\n';
                exit(1);
            }
        }
    }
    
    return;
}




void parseScoreMX() {
    ifstream infile;
    infile.open(MXFILE);
    if (!infile.is_open()) {
        cout << "error opening infile" << endl;
        exit(1);
    }

    // the first row
    vector<char> header;
    string inHeader;
    getline(infile, inHeader);
    stringstream sh(inHeader);
    char chh;
    while (sh >> chh) {
        header.push_back(chh);
    }
    
    // the rest rows
    string inLine = "";
    while (getline(infile, inLine)) {
        stringstream ss(inLine);
        char ch;
        ss >> ch;
        string num = "";
        int i = 0;
        while (ss >> num) {
            string myPair = "";
            myPair.push_back(ch);
            myPair.push_back(header[i]);
            SCORE[myPair] = stoi(num);
            ++i;
        }
    }
    infile.close();
    
    
    //test
//    cout << "pring out the scoring matrix..." << endl;
//    for (auto it = SCORE.begin(); it != SCORE.end(); ++it) {
//        cout << it->first << "\t" << it->second << endl;
//    }
    
    return;
}





void calculateMap() {
    myMap.resize(S2.size()+1);
    
    // the first row
    myMap[0].resize(S1.size()+1);
    myMap[0][0].myScore = 0;
    myMap[0][0].wScore = 0;
    myMap[0][0].nScore = 0;
    myMap[0][0].source = 'x';
    myMap[0][0].wOrigin = make_pair(0, 0);
    myMap[0][0].nOrigin = make_pair(0, 0);
    
    myMap[0][1].myScore = GO;
    myMap[0][1].wScore = GO;
    myMap[0][1].nScore = numeric_limits<int>::min()-GO;
    myMap[0][1].source = 'w';
    myMap[0][1].wOrigin = make_pair(0, 1);
    myMap[0][1].nOrigin = make_pair(0, 1);
    for (int i = 2; i < int(S1.size()+1); ++i) {
        myMap[0][i].wScore = myMap[0][i-1].myScore + GX;
        myMap[0][i].nScore = numeric_limits<int>::min()-GO;
        myMap[0][i].source = 'w';
        myMap[0][i].wOrigin = make_pair(0, 1);
        myMap[0][i].nOrigin = make_pair(0, i);
        myMap[0][i].myScore = myMap[0][i].wScore;
    }

    // the rest rows
    for (int i = 1; i < int(S2.size()+1); ++i) {
        myMap[i].resize(S1.size()+1);
        // the first col
        if (i == 1) {
            myMap[i][0].nScore = myMap[i-1][0].myScore + GO;
        } else {
            myMap[i][0].nScore = myMap[i-1][0].myScore + GX;
        }
        myMap[i][0].source = 'n';
        myMap[i][0].wOrigin = make_pair(i, 0);
        myMap[i][0].nOrigin = make_pair(1, 0);
        myMap[i][0].wScore = numeric_limits<int>::min()-GO;
        myMap[i][0].myScore = myMap[i][0].nScore;
        
        
        // everything else
        for (int j = 1; j < int(S1.size()+1); ++j) {
            compareGrid(i, j);
        }
    }
    return;
}








void compareGrid(int i, int j) {
    assert(i > 0);
    assert(j > 0);

    // calculate wScore
    if (myMap[i][j-1].wScore + GX > myMap[i][j-1].myScore + GO) {
        myMap[i][j].wScore = myMap[i][j-1].wScore + GX;
        myMap[i][j].wOrigin = myMap[i][j-1].wOrigin;
    } else {
        myMap[i][j].wScore = myMap[i][j-1].myScore + GO;
        myMap[i][j].wOrigin = make_pair(i, j);
    }
    
    // calculate nScore
    if (myMap[i-1][j].nScore + GX > myMap[i-1][j].myScore + GO) {
        myMap[i][j].nScore = myMap[i-1][j].nScore + GX;
        myMap[i][j].nOrigin = myMap[i-1][j].nOrigin;
    } else {
        myMap[i][j].nScore = myMap[i-1][j].myScore + GO;
        myMap[i][j].nOrigin = make_pair(i, j);
    }
    
    // west
    myMap[i][j].myScore = myMap[i][j].wScore;
    myMap[i][j].source = 'w';
    
    // north
    if (myMap[i][j].nScore > myMap[i][j].myScore) {
        myMap[i][j].myScore = myMap[i][j].nScore;
        myMap[i][j].source = 'n';
    }
    
    // diagonal
    string myPair = "";
    myPair.push_back(S2[i-1]);
    myPair.push_back(S1[j-1]);
    
    if (SCORE.find(myPair) == SCORE.end()) {
        cerr << "SEQUENCE CHAR NOT FOUND IN SCORING MATRIX" << endl;
        exit(1);
    }
    int diagonal = myMap[i-1][j-1].myScore + SCORE[myPair];
    if (diagonal > myMap[i][j].myScore) {
        myMap[i][j].myScore = diagonal;
        myMap[i][j].source = 'd';
    }
    return;
}






void printMap() {
    cout << "########## PRINT THE CALCULATED MATRIX ##########" << endl;
    // pring the first row
    cout << "\t\t" ;
    for (int i = 0; i < int(S1.size()+1); ++i) {
        cout << S1[i] << "\t";
    }
    cout << endl;
    
    // pring the rest
    for (int i = 0; i < int(myMap.size()); ++i) {
        // the first col
        if (i == 0) {
            cout << "\t";
        } else {
            cout << S2[i-1] << "\t";
        }
        
        // the rest
        for (int j = 0; j < int(myMap[i].size()); ++j) {
            cout << myMap[i][j].myScore << "\t";
        }
        cout << endl;
    }
    cout << "############## END OF PRINTING ##############" << endl << endl;
    return;
}





void backTrace(deque<char> &seq1, deque<char> &seq2) {
    int i = int(S2.size());
    int j = int(S1.size());
    while (myMap[i][j].source != 'x') {
        if (myMap[i][j].source == 'd') {
            aligLn++;
            seq1.push_front(S1[j-1]);
            seq2.push_front(S2[i-1]);
            if (S1[j-1] == S2[i-1]) {
                idenLn++;
            }
            --i;
            --j;
        } else if (myMap[i][j].source == 'w') {
            pair<int, int> temp = myMap[i][j].wOrigin;
            while (i != temp.first || j != temp.second) {
                seq1.push_front(S1[j-1]);
                seq2.push_front('-');
                --j;
            }
            seq1.push_front(S1[j-1]);
            seq2.push_front('-');
            --j;
            
        } else if (myMap[i][j].source == 'n') {
            pair<int, int> temp = myMap[i][j].nOrigin;
            while (i != temp.first || j != temp.second) {
                seq1.push_front('-');
                seq2.push_front(S2[i-1]);
                --i;
            }
            seq1.push_front('-');
            seq2.push_front(S2[i-1]);
            --i;
        }
    }

    
    if (!quiet) {
        cout << "############## ALIGNMENT SUMMARY ##############" << endl;
        cout << "Seq1 length: " << S1.size() << endl;
        cout << "Seq2 length: " << S2.size() << endl;
        cout << "Aligned length: " << aligLn << endl;
        cout << "Identical length: " << idenLn << endl;
        cout << "Sequence Identity: " << double(idenLn)/double(S2.size()) << endl;
        cout << "################ END OF SUMMARY ################" << endl << endl;
    } else {
        cout << double(idenLn)/double(S2.size()) << endl;
    }
    
    return;
}





void printPath(const deque<char> &seq1, deque<char> &seq2) {
    assert(seq1.size() == seq2.size());
    cout << "########## PRINT THE ALIGNED SEQUENCES ##########" << endl;
    int i = 0;
    while (i < int(seq1.size())) {
        for (int j = i; j < i+80; ++j) {
            if (j >= int(seq1.size())) {
                break;
            }
            cout << seq1[j];
        }
        cout << endl;
        for (int j = i; j < i+80; ++j) {
            if (j >= int(seq1.size())) {
                break;
            }
            if (seq1[j] == seq2[j]) {
                cout << ":";
            } else {
                cout << " ";
            }
        }
        cout << endl;
        for (int j = i; j < i+80; ++j) {
            if (j >= int(seq1.size())) {
                break;
            }
            cout << seq2[j];
        }
        cout << endl << endl;
        
        i+= 80;
    }
    cout << "############## END OF PRINTING ##############" << endl << endl;
    return;
}




