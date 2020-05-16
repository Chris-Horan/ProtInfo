/**********************
 *  Purpose: Represent a protein.
 *  Author: Chris Horan
 *  Date: Feb 12th, 2020
 * 
 * TODO:
 *  -toUpper(string)
 *  -toUpper(char)
 */

#ifndef PROTEIN_H
#define PROTEIN_H

#include <string>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

/**********************
 * Represents a protein.
 * 
 * Functionality: 
 * Add a residue: append(char), void
 * Returns the sequence of residues: seq(), string
 */
class Protein {
    private:
        /**************
         * Library of data on amino acids.
         * Index:
         *  [0]: Full name of AA
         *  [1]: Molecular weight of residue in Daltons.
         *  [2]: Volume in Angstroms of Amino Acid.
         */
        map<char, vector<string>> aaLib = {
            {'A', {"Alanine",                   "89.09",    "88.6"}},
            {'R', {"Arginine",                  "174.20",   "173.4"}},
            {'N', {"Asparagine",                "132.12",   "114.1"}},
            {'D', {"Aspartic Acid",             "133.10",   "111.1"}},
            {'C', {"Cysteine",                  "121.15",   "108.5"}},
            {'Q', {"Glutamine",                 "146.15",   "143.8"}},
            {'E', {"Glutamic Acid",             "147.13",   "138.4"}},
            {'G', {"Glycine",                   "75.07",    "60.1"}},
            {'H', {"Histidine",                 "155.16",   "153.2"}},
            {'I', {"Isoleucine",                "131.17",   "166.7"}},
            {'L', {"Leucine",                   "131.17",   "166.7"}},
            {'K', {"Lysine",                    "146.19",   "168.6"}},
            {'M', {"Methionine",                "149.21",   "162.9"}},
            {'F', {"Phenylalanine",             "165.19",   "189.9"}},
            {'P', {"Proline",                   "115.13",   "112.7"}},
            {'S', {"Serine",                    "105.09",   "89.0"}},
            {'T', {"Threonine",                 "119.12",   "116.1"}},
            {'W', {"Tryptophan",                "204.23",   "227.8"}},
            {'Y', {"Tyrosine",                  "181.19",   "193.6"}},
            {'V', {"Valine",                    "117.15",   "140.0"}},
            {'B', {"Asparagine/Aspartic Acid",  "132.67",   "112.6"}},
            {'Z', {"Glutamine/Glutamic Acid",   "146.76",   "146.6"}}
        };

        struct AminoAcid {
            AminoAcid* prev;
            char ID;
            AminoAcid* next;
        };

        AminoAcid* cTerm;
        AminoAcid* nTerm;
        int len;

        void init() {
            cTerm = nullptr;
            nTerm = nullptr;
            len = 0;
        }

        void clear(){
            AminoAcid* curr = cTerm;
            AminoAcid* tmp;
            while(curr != nullptr) {
                tmp = curr;
                curr = curr->next;
                delete tmp;
            }
        }

        //Returns the AminoAcid object reference at the specified index
        AminoAcid* aaAt(int n) {
            if(n > len) return nullptr;
            AminoAcid* curr = cTerm;
            for(int i = 0; i<n; i++) {
                curr = curr->next;
            }
            return curr;
        }

        string toUpper(string old) {
            string newS = "";
            char curr = '\0';
            for(int i = 0; i<old.size(); i++) {
                curr = old[i];
                curr = toUpper(curr);
                newS += curr;
            }
            return newS;
        }

        char toUpper(char old) {
            char newC = old;
            if(!((newC > 64)&&(newC < 91))) {
                newC -= 32;
            }
            return newC;
        }

        bool checkChars(string seq) {
            string allowed = "ARNDBCQEZGHILKMFPSTWYV";
            int size = seq.size();
            for(int i = 0; i<size; i++) {
                if(allowed.find(seq[i])==string::npos)
                    return false;
            }
            return true;
        }

        bool checkChars(const char seq) {
            string allowed = "ARNDBCQEZGHILKMFPSTWYV";
            if(allowed.find(seq)==string::npos)
                return false;
            return true;
        }

    public:
        //Standard constructor
        Protein() {
            init();
        }

        //Constructor with sequence
        Protein(string sequence) {
            init();
            sequence = toUpper(sequence);
            if(!checkChars(sequence))
                return;
            int i;
            for(i = 0; i < sequence.size(); i++) {
                append(sequence[i]);
            }
            len = i;
        }

        //Standard destructor
        ~Protein() {
            clear();
            delete cTerm;
            delete nTerm;
        }

        //Append a new amino acid to your protein.
        void append(const char abbr) {
            if(!checkChars(abbr))
                return;
            AminoAcid* newAA = new AminoAcid{nTerm,abbr,nullptr};
            if(cTerm == nullptr) {
                cTerm = newAA;
                nTerm = newAA;
            } else {
                nTerm->next = newAA;
                nTerm = newAA;
            }
            len++;
        }

        // Insert a residue at the specified index.
        void insert(const char abbr, int n) {
            if(!checkChars(abbr))
                return;
            if(n > len || n < 0) return; //ERROR
            if(n == len) {
                append(abbr);
                return;
            }
            if(n == 0) {
                AminoAcid* newAA = new AminoAcid{nullptr,abbr,cTerm};
                cTerm->prev = newAA;
                cTerm = newAA;
                len++;
                return;
            }
            AminoAcid* curr = cTerm;
            for(int i = 0; i < n - 1; i++) {
                curr = curr->next;
            }
            AminoAcid* newAA = new AminoAcid{curr,abbr,curr->next};
            curr->next = newAA;
            newAA->next->prev = newAA;
            len++;
        }

        //Returns the length of the protein as an integer.
        int size() {
            return len;
        }

        //Returns the molecular weight of the protein in Daltons
        //  as a double.
        double protWeight() {
            if(cTerm==nullptr) return 0;
            AminoAcid* curr = cTerm;
            double weight = 0;
            int cntr = 0;
            double water = 18.0153;
            while(curr != nullptr) {
                weight += stod(aaLib[(curr->ID)][1]);
                cntr++;
                curr = curr->next;
            }
            water = water * (cntr - 1);
            return (weight - water) / 1000;
        }

        //Returns the sequence of residues C-N as a string.
        string seq() {
            AminoAcid* curr = cTerm;
            string prot = "";
            while(curr != nullptr) {
                prot += curr->ID;
                curr = curr->next;
            }
            return prot;
        }

        //Returns the full name of the residue 
        string residueAt(int n) {
            AminoAcid* aa = aaAt(n);
            return aaLib[(aa->ID)][0];
        }
};
#endif