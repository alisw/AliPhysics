/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Implementation for correlation analysis

//===========================================================
// AliJJtCorrelations.h
//
// Author: Jussi Viinikainen
//===========================================================

#ifndef ALIJJTCORRELATIONS_H
#define ALIJJTCORRELATIONS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <TRandom3.h>

#include "AliJJtHistograms.h"
#include "../AliJAcceptanceCorrection.h"
#include "../AliJCorrelationInterface.h"

using namespace std;

class AliJJtCorrelations;
class AliJJtHistograms;
class AliJBaseTrack;
class AliJCard;

class AliJJtCorrelations : public AliJCorrelationInterface{
  
public:
  
  AliJJtCorrelations(); // default constructor
  AliJJtCorrelations( AliJCard *cardIn, AliJJtHistograms *histosIn); // constructor
  AliJJtCorrelations(const AliJJtCorrelations& in); // copy constructor
  virtual ~AliJJtCorrelations(); //destructor
  AliJJtCorrelations& operator=(const AliJJtCorrelations& obj); // equal sign operator
  
  void PrintOut(){cout<<"Number of events = "<<fnReal<<"  Mixed events = "<<fnMix<<endl;} // Event count print
  
  void FillHisto(corrFillType cFTyp, fillType fTyp, int cBin, int zBin, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2); // correlation histogram filler based on correlation type
  void FillCorrelationHistograms (fillType fTyp, int cBin, int zBin, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2); // correlation histogram filler
  
  void SetSamplingInclusive(){fsamplingMethod = 1;} // Setter for inclusive sampling
  void UseZVertexAcceptance(bool useZ){fUseZVertexBinsAcceptance = useZ;} // Setter for fUseZVertexBinsAcceptance
  void SetAcceptanceCorrection(AliJAcceptanceCorrection *accCorr){fAcceptanceCorrection = accCorr;} // Setter for acceptance correction
  
  
protected:
  
  AliJCard* fcard; // Card with binning information etc.
  AliJJtHistograms* fhistos;  // Histograms needed in the analysis
  AliJAcceptanceCorrection *fAcceptanceCorrection;  // Acceptance correction container
  int fnReal; // Number of events
  int fnMix; // Number of mixed events
  int fsamplingMethod; // Sampling method flag (flat or inclusive)
  double fmaxEtaRange; // Eta range used in the analysis
  
  TRandom3 *frandom; // Random number generator
  
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
  int fZBin; // z-vertex bin
  
  bool fNearSide;  // true if near side correlation, false if away side
  bool fNearSide3D; // near side defined by the half ball around the trigger
  int fEtaGapBin;  // Bin index for the current eta gap
  int fPhiGapBinNear;  // Bin index for the phi gap in the near side
  int fRGapBinNear;  // Bin index for the R gap in the near side
  int fCentralityBin;  // Bin index for the centrality bin
  int fXlongBin;  // Bin index for xlong bin

  bool fIsLikeSign; // True = like sign correlation, false = unlike sign correlation
  
  bool fUseZVertexBinsAcceptance; // false = integrate over z-vertex bins, true = Do acceptance correction in z-vertex bins
  
private:
  
  void FillJtHistograms(fillType fTyp, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2, bool fill2DBackground); // jT histogram filler
  void FillDeltaEtaHistograms(fillType fTyp); // deltaEta histogram filler
  void FillDeltaEtaDeltaPhiHistograms(fillType fTyp); // deltaEta deltaPhi histogram filler
  void FillPtaHistograms(fillType fTyp); // pTa histogram filler
  void FillJtDistributionHistograms(fillType fTyp, int assocType, TLorentzVector *vTrigger, TLorentzVector *vAssoc, AliJTH1D &hDistribution, AliJTH1D &hDistributionLikeSign, AliJTH1D &hDistributionUnlikeSign, AliJTH1D &hInvariantMass, AliJTH1D &hInvariantMassLikeSign, AliJTH1D &hInvariantMassUnlikeSign); // jT distribution filler
  void FillJtBackgroundHistograms(int assocType, int gapType, TLorentzVector *vTrigger, TLorentzVector *vAssoc, AliJTH1D &hBackground, AliJTH1D &hBackgroundLikeSign, AliJTH1D &hBackgroundUnlikeSign, AliJTH1D &hPtAssoc, AliJTH2D &hBackground2D, bool fill2DBackground); // jT background filler
  int GetBinIndex(int assocType, TLorentzVector *vTrigger, TLorentzVector *vAssoc); // Bin index getter
};

#endif






















