/**********************
 *  Purpose: Represent and manipulate an amino acid sequence.
 *  Author: Chris Horan
 *  Date: April 10th, 2020
 *  
 * 
 */

#include "Protein.h"
#include <iostream>

using namespace std;

int main() {
    string aaSeq = "";
    cout << "Provide an amino acid sequence: ";
    cin >> aaSeq;
    Protein p(aaSeq);
    cout<<p.seq()<<endl;
    cout<<"This protein contains "<<p.size()<< " residues." << endl;
    cout<<"This protein weighs "<<p.protWeight()<<" kiloDaltons."<<endl;
    return 0;
}