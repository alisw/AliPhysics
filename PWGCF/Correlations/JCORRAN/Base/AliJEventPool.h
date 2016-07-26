//===========================================================
// AliJEventPool.h
//
//===========================================================

#ifndef ALIJEVENTPOOL_H
#define ALIJEVENTPOOL_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

#include <AliJConst.h>

class TClonesArray;
class AliJBaseTrack;
class AliJPhoton;
class AliJTrack;
class AliJCard;
class AliJCorrelations;
class AliJHistos;
class TH1D;

#define   MAXNOEVENT 1000    // Maximum no of events in pools (400 used for QM anal.) 

class AliJEventPool {

    public:
      AliJEventPool(AliJCard *cardin, AliJHistos *histosin, AliJCorrelations *coin, particleType particle );
      virtual ~AliJEventPool( );
      AliJEventPool(const AliJEventPool& obj);
      AliJEventPool& operator=(const AliJEventPool& obj);
  

    public:
        void Mix( TClonesArray *triggList, 
                corrFillType cFTyp, 
                float cent, float Z, float thisMult, int iev, bool leadingParticle = false);

       //void MixRNDM( AliJEventPool *cross, void (AliJCorrelations::*fillHisto)(fillType, int, AliJBaseTrack*, AliJBaseTrack*) );

        void AcceptList(TClonesArray *inList, float cent, float Z, float inMult, int iev);

        void Mysample(TH1D *fromh, TH1D *toh );
        void PrintOut(){for(int i=0;i<kMaxNoCentrBin;i++)
            cout<<"c: "<<i<<" mixed "<<fnoMix[i]<<" accepted "<<fnoMixCut[i]<<" "<<(fnoMix[i]>0?fnoMixCut[i]*1.0/fnoMix[i]:0)<< endl;}

    protected:

        int   fevent[kMaxNoCentrBin][MAXNOEVENT];  // comment me
        float fZVertex[kMaxNoCentrBin][MAXNOEVENT];  // comment me
        float fcentrality[kMaxNoCentrBin][MAXNOEVENT];  // comment me
        float fmult[kMaxNoCentrBin][MAXNOEVENT];  // comment me
        long  flastAccepted[kMaxNoCentrBin];  // comment me
        long fwhereToStore[kMaxNoCentrBin];   // comment me
        long fnoMix[kMaxNoCentrBin];  // comment me
        long fnoMixCut[kMaxNoCentrBin];   // comment me

        TClonesArray   *fLists[kMaxNoCentrBin][MAXNOEVENT]; // mix lists
        AliJCard  *fcard;  // card
        AliJCorrelations *fcorrelations; // correlation object
        AliJHistos *fhistos;  // histos
        //AliJBaseTrack *ftk; // track
        //AliJBaseTrack *ftk1; // track
        //AliJBaseTrack *ftk2; // track
        particleType fthisPoolType; // pool type

        TClonesArray  *fpoolList;  // pool list

        //int   trials[MAXNOEVENT];

};

#endif






















