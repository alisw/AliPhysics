/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Implementation for correlation analysis

//===========================================================
// AliJIaaCorrelations.h
//===========================================================

#ifndef ALIJIAACORRELATIONS_H
#define ALIJIAACORRELATIONS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <TRandom3.h>

#include "AliJIaaHistos.h"
#include "../AliJAcceptanceCorrection.h"
#include "../AliJCorrelationInterface.h"

using namespace std;

class AliJIaaCorrelations;
class AliJIaaHistos;
class AliJBaseTrack;
class AliJCard;

class AliJIaaCorrelations : public AliJCorrelationInterface{
  
public:

    AliJIaaCorrelations( AliJCard *cardIn, AliJIaaHistos *histosIn);

    virtual ~AliJIaaCorrelations(){;}    //destructor
    AliJIaaCorrelations();
    AliJIaaCorrelations(const AliJIaaCorrelations& in);
    AliJIaaCorrelations& operator=(const AliJIaaCorrelations& obj);

    void PrintOut(){cout<<"Real correl = "<<fnReal<<"  mixed = "<<fnMix<<endl;}

    void FillHisto(corrFillType cFTyp, fillType fTyp,    int cBin, int zBin, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2);
    void FillCorrelationHistos (fillType fTyp,    int cBin, int zBin, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2);

    void SetSamplingInclusive(){fsamplingMethod = 1;}
    void SetAcceptanceCorrection(AliJAcceptanceCorrection *accCorr){fAcceptanceCorrection = accCorr;}

    double DeltaPhi(double phi1, double phi2);

protected:
  
    AliJCard*   fcard; // card
    AliJIaaHistos* fhistos;  // histos
    AliJAcceptanceCorrection *fAcceptanceCorrection;  // acceptance correction container
    int fnReal; // comment me
    int fnMix; // comment me
    int fsamplingMethod; // comment me
    double fmaxEtaRange; // comment me

    TRandom3 *frandom; // comment me

    double fptt;  // pT of the trigger particle in the correlation loop
    double fpta;  // pT ot the asociated particle in the correlation loop
    double fTrackPairEfficiency;  // pair efficiency for the tracks in the correlation loop
    int fpttBin;  // Bin index for the trigger pT bin
    int fptaBin;  // Bin index for the associated pT bin
    double fPhiTrigger;  // Azimuthal angle of the trigger particle
    double fPhiAssoc;  // Asimuthal angle of the associated particle
    double fEtaTrigger;  // Pseodurapidity of the trigger particle
    double fEtaAssoc;  // Pseudorapidity of the associated particle
    double fDeltaPhi;  // Difference of the azimuthal angles of trigger and associated particles
    double fDeltaPhiPiPi;  // The same as above but measured from -pi to pi
    double fDeltaEta;  // Difference of the pseudorapidities of the trigger and associated particles
    double fXlong;  // The xlong value of the trigger and associated particles

    bool fNearSide;  // true if near side correlation, false if away side
    bool fNearSide3D; // near side defined by the half ball around the trigger
    int fEtaGapBin;  // Bin index for the current eta gap
    int fPhiGapBinNear;  // Bin index for the phi gap in the near side
    int fPhiGapBinAway;  // Bin index for the phi gap in the away side
    int fRGapBinNear;  // Bin index for the R gap in the near side
    int fRGapBinAway;  // Bin index for the R gap in the away side
    int fCentralityBin;  // Bin index for the centrality bin
    int fXlongBin;  // Bin index for xlong bin

    bool fIsLikeSign; // True = like sign correlation, false = unlike sign correlation

    double fGeometricAcceptanceCorrection;   // Acceptance correction due to the detector geometry
    double fGeometricAcceptanceCorrection3D; // Acceptance correction due to the detector geometry for 3D near side

private:
  
    void FillDeltaEtaHistograms(fillType fTyp, int ZBin); // deltaEta histogram filler
    void FillDeltaPhiHistograms(fillType fTyp); // deltaPhi histogram filler
    void FillDeltaEtaDeltaPhiHistograms(fillType fTyp); // deltaEta deltaPhi histogram filler
    void FillPtaHistograms(fillType fTyp); // pTa histogram filler
    int GetBinIndex(int assocType, TLorentzVector *vTrigger, TLorentzVector *vAssoc);

};

#endif // ALIJIAACORRELATIONS_H
