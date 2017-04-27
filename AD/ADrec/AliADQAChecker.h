// -*- C++ -*-
#ifndef ALIADQACHECKER_H
#define ALIADQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/*
  Checks the quality of the data
  by comparing with reference data
  which should be loaded from QA ref DB
*/

// --- ROOT system ---
class TFile;
class TH1;
class TH2;
class TObjArray;
class TVirtualPad;
class TLegend;

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQACheckerBase.h"

class AliADLoader;
class AliCDBManager;
class AliCDBStorage;
class AliADQAParam;

class AliADQAChecker: public AliQACheckerBase {
public:
  AliADQAChecker();
  virtual ~AliADQAChecker(); // destructor

  virtual void   Init(const AliQAv1::DETECTORINDEX_t det);

  AliADQAParam *GetQAParam() const;
  void SetLowEventCut(Int_t nEvents) {fLowEventCut = nEvents;}
  void SetORvsANDCut(Double_t cut) {fORvsANDCut = cut;}
  void SetBGvsBBCut(Double_t cut) {fBGvsBBCut = cut;}
  void SetSatMedCut(Double_t cut) {fSatMed = cut;}
  void SetSatHighCut(Double_t cut) {fSatHigh = cut;}
  void SetSatHugeCut(Double_t cut) {fSatHuge = cut;}
  void SetMaxPedDiffCut(Double_t cut) {fMaxPedDiff = cut;}

protected:
  virtual void Check( Double_t * test, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam * recoParam);
  Double_t CheckRaws(TObjArray * list) const;
  Double_t CheckPedestals(TObjArray * list) const;
  Double_t CheckEsds(TObjArray * list) const;

  virtual void   MakeImage( TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode) ;
  virtual void SetQA(AliQAv1::ALITASK_t index, Double_t * value) const ;

private:
  AliADQAChecker(const AliADQAChecker& qac); // cpy ctor
  AliADQAChecker &operator=(const AliADQAChecker& qac); // assignment operator

  AliADQAParam *fQAParam;
  Int_t    fLowEventCut; // Minimum number of events required by the QA checker
  Float_t fORvsANDCut; // AD OR vs AD AND counters cut
  Float_t fBGvsBBCut; // AD beam-gas vs beam-beam counters cut
  Float_t fSatMed; //Medium saturation cut
  Float_t fSatHigh; //High saturation cut
  Float_t fSatHuge; //Very high saturation cut
  Int_t fMaxPedDiff; //Pedestal difference cut
  Float_t fMaxPedWidth; //Pedestal width cut
  Int_t fChargeChannelZoomMin; //Min for Zoom on charge
  Int_t fChargeChannelZoomMax; //Max for Zoom on charge
  Float_t fTimeRatioBBZoomMin; //Min for Zoom time/BB ratio
  Float_t fTimeRatioBBZoomMax; //Max for Zoom time/BB ratio
  Float_t fTimeRatioBGZoomMin; //Min for Zoom time/BG ratio
  Float_t fTimeRatioBGZoomMax; //Max for Zoom time/BG ratio
  Float_t fChargeTrendMin;
  Float_t fChargeTrendMax;
  Float_t fMaxNoTimeRate;
  Float_t fMaxNoFlagRate;
  Float_t fMaxBBVariation;
  Float_t fMaxBGVariation;
  Float_t fAsynchronBB;
  Float_t fAsynchronBG;
  
  TVirtualPad* pCharge;
  TVirtualPad* pChargeZoom;
  TVirtualPad* pTime;
  TVirtualPad* pTimeRatio;
  TVirtualPad* pMeanTime;
  TVirtualPad* pClockFg;
  TVirtualPad* pPed;
  TVirtualPad* pMaxCh;
  TVirtualPad* pCoinc;
  TVirtualPad* pTrigger;
  TVirtualPad* pDecision;
  TVirtualPad* pChargeTrend;
  
  TH1* histo;
  TH1* histoBlue;
  TH1* histoRed;
  TH1* histoDenominator;
  TH1* histoRatio;
  TH2* histo2D;
  TH2* histo2DCopy;
  
  TLegend *myLegend1;
  TLegend *myLegend2;
  TLegend *myLegend3;

  ClassDef(AliADQAChecker,7);  // description
};

#endif // AliADQAChecker_H
