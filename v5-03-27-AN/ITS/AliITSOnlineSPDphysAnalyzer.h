#ifndef ALIITSONLINESPDPHYSANALYZER_H
#define ALIITSONLINESPDPHYSANALYZER_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// This class is used in the detector algorithm framework //
// to process the data stored in special container files  //
// (see AliITSOnlineSPDphys).                             //
////////////////////////////////////////////////////////////

#include <TString.h>

class AliITSOnlineSPDphys;
class AliITSOnlineCalibrationSPDhandler;
class TGraph;
class TH2F;

class AliITSOnlineSPDphysAnalyzer {

 public:
  AliITSOnlineSPDphysAnalyzer(const Char_t *fileName, AliITSOnlineCalibrationSPDhandler* handler, Bool_t readFromGridFile=kFALSE);
  AliITSOnlineSPDphysAnalyzer(AliITSOnlineSPDphys* physObj, AliITSOnlineCalibrationSPDhandler* handler);
  AliITSOnlineSPDphysAnalyzer(const AliITSOnlineSPDphysAnalyzer& handle);
  ~AliITSOnlineSPDphysAnalyzer();

  AliITSOnlineSPDphysAnalyzer& operator=(const AliITSOnlineSPDphysAnalyzer& handle);

  void       SetCalibHandler(AliITSOnlineCalibrationSPDhandler *handler) {fHandler=handler;}
  void       SetParam(const Char_t *pname, const Char_t *pval);
  void       ReadParamsFromLocation(const Char_t *dirName);

  UInt_t     ProcessDeadPixels();
  UInt_t     ProcessNoisyPixels();
  UInt_t     ProcessNoisyPixels(UInt_t eq, UInt_t nrEvts);

  UInt_t     GetNrEnoughStatChips();
  UInt_t     GetNrDeadChips();
  UInt_t     GetNrInefficientChips();
  UInt_t     GetNrNeedsMoreStatChips();

  AliITSOnlineSPDphys* GetOnlinePhys() {return fPhysObj;}
  UInt_t     GetEqNr() const;
  UInt_t     GetNrEvents() const;

  TH2F*      GetHitMapTot();
  TH2F*      GetPhysicalHitMapTot();
  TH2F*      GetHitMapChip(UInt_t hs, UInt_t chip);

 private:
  TString  fFileName; // container file name
  enum     calibvals{kMINTH,kMEANTH,kDAC,kUNIMA,kNOISE,kDELAY};  // calib types
  AliITSOnlineSPDphys               *fPhysObj; // container obj
  AliITSOnlineCalibrationSPDhandler *fHandler; // calib helper obj
  void     Init(Bool_t readFromGridFile=kFALSE);    // initialization
  void     Exponent(Double_t &val, Int_t &valExp) const;

  UInt_t   fNrEnoughStatChips;    // nr of enough stat chips
  UInt_t   fNrDeadChips;          // nr of dead chips
  UInt_t   fNrInefficientChips;   // nr of inefficient chips

  Double_t fNrEqHits;         // total nr of hits for associated eq
  Bool_t   fbDeadProcessed;   // flag to tell if ProcessDeadPixels has been called

  // dead noisy parameters:
  Double_t fThreshNoisy;       // at what confidence level do we search for noisy
  Double_t fThreshDead;        // at what confidence level do we search for dead
  UInt_t   fMinEventsForNoisy; // min nr of events required to try noisy algorithm
  UInt_t   fMinEventsForDead;  // min nr of events required to try dead algorithm
  Float_t  fDefinitelyNoisyRatio; // if a pixel fires more than this ratio of the events, it must be noisy
  Double_t fMinNrEqHitsForDeadChips; // minimum nr of hits for eq to assign dead chip
  Double_t fRatioToMeanForInefficientChip; // ratio to mean nr of hits per chip to assign ineff. chip
  

};

#endif
