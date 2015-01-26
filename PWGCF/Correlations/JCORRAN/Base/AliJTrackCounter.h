/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

//===========================================================
// AliJTrackCounter.h
//   Created  Fri Mar 11 15:24:11 EET 2011  by classmaker
//   Jan Rak
//===========================================================

#ifndef ALIJTRACKCOUNTER_H
#define ALIJTRACKCOUNTER_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

class AliJTrackCounter : public TObject {

    public:
        AliJTrackCounter();    //constructor

        AliJTrackCounter(int inind0, int inind1, int inlist0, int inlist1, double indphi);

        virtual ~AliJTrackCounter(){;}    //destructor


        void Reset();
        bool Exists() const {return fpt>=0 ? true : false;}
        void Store(int inind, double inpt, int inBin){findex=inind; fpt=inpt; fptBin=inBin;}
        int GetIndex() const {return findex;}
        void SetIndex(int inIdx) {findex=inIdx;}
        double Pt() const {return fpt;} // BS to compatible with AliJBaseTrack
        double GetPt() const {return fpt;}
        int GetPtBin() const {return fptBin;}

        int GetPairTrackID(int ip ) const {return (ip==0||ip==1) ? fpairTrackID[ip] : -1; }
        int GetPairListID(int ip ) const {return (ip==0||ip==1) ? fpairListID[ip] : -1; }
        double GetPairDPhi() const {return fdphi;}

        void Print(Option_t* option = "") const{
           // We must have option here because Print overrides the same function from TObject
           // TODO: make some sensible use of option
            cout<<"LPindex="<<findex << option <<" LPpt="<<fpt <<" bin= "<<fptBin 
                <<" fpairTrackID[0]="<<fpairTrackID[0] <<" fpairTrackID[1]="<<fpairTrackID[1] <<" pairDPHI="<< fdphi<< endl;    
        }

       AliJTrackCounter& operator=(const AliJTrackCounter& counter);

    protected:

        int findex;  // comment me
        
        int fpairTrackID[2];  // comment me
        int fpairListID[2];  // comment me
        int fptBin;  // comment me
        double fpt;  // comment me
        double fdphi;  // comment me

        ClassDef(AliJTrackCounter,1)
};

#endif






















