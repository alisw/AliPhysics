/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

//===========================================================
// AliJCorrelations.h
//   Created  Thu Apr 17 12:40:29 EEST 2008  by classmaker
//   Jan Rak
//===========================================================

#ifndef ALIJCORRELATIONS_H
#define ALIJCORRELATIONS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <TRandom3.h>  //FK//

#include  <AliJConst.h>

using namespace std;

class AliJCorrelations;
class AliJHistos;
class AliJBaseTrack;
class AliJCard;

class AliJCorrelations {
  
public:
  
  AliJCorrelations( AliJCard *cardIn, AliJHistos *histosIn);
  
  virtual ~AliJCorrelations(){;}    //destructor
  AliJCorrelations();
  AliJCorrelations(const AliJCorrelations& in);
  AliJCorrelations& operator=(const AliJCorrelations& obj);
  
  void PrintOut(){cout<<"Real correl = "<<fnReal<<"  mixed = "<<fnMix<<endl;}
  
  void FillHisto(corrFillType cFTyp, fillType fTyp,    int cBin, int zBin, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2);
  void FillAzimuthHistos (fillType fTyp,    int cBin, int zBin, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2);
  
  double GetGeoAccCorrFlat(double deltaEta);
  double GetGeoAccCorrIncl(double deltaEta);
  
  void SetSampligInclusive(){fsamplingMethod = 1;}
  
  
  double DeltaPhi(double phi1, double phi2);
  
  //double DeltaPhi(double phi1, double phi2) {
  //    double res =  atan2(sin(phi1-phi2), cos(phi1-phi2));
  //    //return res>-kJPi/3.0 ? res : kJTwoPi+res ;
  //    return res > dPhiRange ? res : kJTwoPi+res ;
  //}
  
protected:
  
  AliJCard*   fcard; // card
  AliJHistos* fhistos;  // histos
  int fnReal; // comment me
  int fnMix; // comment me
  int fsumTriggerAndAssoc; // comment me
  int fsamplingMethod; // comment me
  int fIsHeavyIon; // comment me
  double fawayPhiGap; // comment me
  double fDPhiUERegion[2]; // comment me
  double fRGap[30]; // comment me
  double fmaxEtaRange; // comment me
  
  int fRSignalBin; // comment me
  
  TRandom3 *frandom; // comment me
  
  double fptt;  // pT of the trigger particle in the correlation loop
  double fpta;  // pT ot the asociated particle in the correlation loop
  double fTrackPairEfficiency;  // pair efficiency for the tracks in the correlation loop
  bool fIsIsolatedTrigger;  // Tells whether the trigger is isolated or not
  int fpttBin;  // Bin index for the trigger pT bin
  int fptaBin;  // Bin index for the associated pT bin
  double fPhiTrigger;  // Azimuthal angle of the trigger particle
  double fPhiAssoc;  // Asimuthal angle of the associated particle
  double fDeltaPhi;  // Difference of the azimuthal angles of trigger and associated particles
  double fDeltaPhiPiPi;  // The same as above but measured from -pi to pi
  double fDeltaEta;  // Difference of the pseudorapidities of the trigger and associated particles
  double fXlong;  // The xlong value of the trigger and associated particles
  
  bool fNearSide;  // true if near side correlation, false if away side
  int fEtaGapBin;  // Bin index for the current eta gap
  int fPhiGapBinNear;  // Bin index for the phi gap in the near side
  int fPhiGapBinAway;  // Bin index for the phi gap in the away side
  int fRGapBinNear;  // Bin index for the R gap in the near side
  int fRGapBinAway;  // Bin index for the R gap in the away side
  int fCentralityBin;  // Bin index for the centrality bin
  int fXlongBin;  // Bin index for xlong bin
  
  double fGeometricAcceptanceCorrection;  // Acceptance correction due to the detector geometry
  
private:
  
  void FillPairPtAndCosThetaStarHistograms(fillType fTyp, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2);
  void FillXeHistograms(fillType fTyp);
  void FillJtHistograms(fillType fTyp, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2, bool fill2DBackground);
  void FillDeltaEtaHistograms(fillType fTyp, int zBin);
  void FillDeltaPhiHistograms(fillType fTyp);
  void FillPtaHistograms(fillType fTyp);
  void FillIAAAndMoonHistograms(fillType fTyp, int zBin);
};

#endif






















