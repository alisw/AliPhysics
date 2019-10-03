#ifndef AliTwoPlusOneContainer_H
#define AliTwoPlusOneContainer_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// data container for 2+1 particle analysis

#include "TNamed.h"
#include "AliUEHist.h"
#include "THn.h"
// #include "THn.h" // in cxx file causes .../THn.h:257: error: conflicting declaration ‘typedef class THnT<float> THnF’

class AliVParticle;

class TList;
class TSeqCollection;
class TObjArray;
class TH1F;
class TH2F;
class TH3F;

class AliTwoPlusOneContainer : public TNamed
{
 public:
  AliTwoPlusOneContainer(const char* name = "AliTwoPlusOneContainer", const char* uEHist_name = "TwoPlusOne", const char* binning = 0, Double_t alpha = 0.2);

  virtual ~AliTwoPlusOneContainer();
  
  enum PlotKind {kSameNS = 0, kSameAS, kMixedNS, kMixedAS, kMixedCombNS, kMixedCombAS, k1plus1, kBackgroundSameNS, kBackgroundSameAS, kMixed1plus1, kParticleDist, kMixedMixedCombNS, kMixedMixedCombAS, kMixedBackgroundSameNS, kMixedBackgroundSameAS};

  Int_t FillCorrelations(Double_t centrality, Float_t zVtx, AliTwoPlusOneContainer::PlotKind step, TObjArray* triggerNear, TObjArray* triggerAway, TObjArray* assocNear, TObjArray* assocAway, Double_t weight, Bool_t is1plus1, Bool_t isBackgroundSame, Bool_t applyEfficiency);
  void FillParticleDist(Double_t centrality, Float_t zVtx, TObjArray* particleDist, Double_t weight, Bool_t applyEfficiency);
  
  AliUEHist* GetData() {return fTwoPlusOne;}
  TH1F* GetAsymmetry() {return fAsymmetry;}
  TH1F* GetAsymmetryMixed() {return fAsymmetryMixed;}
  TH2F* GetTriggerPt() {return fTriggerPt;}
  Double_t getTriggerPt1Min() {return fTriggerPt1Min;}
  Double_t getTriggerPt1Max() {return fTriggerPt1Max;}
  Double_t getTriggerPt2Min() {return fTriggerPt2Min;}
  Double_t getTriggerPt2Max() {return fTriggerPt2Max;}
  Double_t getPtAssocMin() {return fPtAssocMin;}
  Double_t getPtAssocMax() {return fPtAssocMax;}
  
  void SetData(AliUEHist* obj) {fTwoPlusOne = obj; }
  void SetAlpha(Double_t value) {fAlpha = value;}
  void SetUseLeadingPt(Bool_t flag) { fUseLeadingPt = flag; }
  void SetUseAllT1(Bool_t flag) { fUseAllT1 = flag; }
  void SetUseBackgroundSameOneSide(Bool_t flag) { fUseBackgroundSameOneSide = flag; }
  void SetUseSmallerPtAssoc(Bool_t flag) { fUseSmallerPtAssoc = flag; }
  void SetEfficiencyCorrection(THnF* hist)   { fEfficiencyCorrection = hist;   }

  AliTwoPlusOneContainer(const AliTwoPlusOneContainer &c);
  AliTwoPlusOneContainer& operator=(const AliTwoPlusOneContainer& c);
  virtual void Copy(TObject& c) const;
  virtual Long64_t Merge(TCollection* list);
 
protected:
  void DeleteContainers();
  Double_t getEfficiency(Double_t pt, Double_t eta, Double_t centrality, Double_t zVtx);
  
  AliUEHist* fTwoPlusOne;	     //a 7 dim histogram which actually contains all the data

  TH1F* fAsymmetry;                  //asymmetry of the same event
  TH1F* fAsymmetryMixed;             //asymmetry of the mixed event

  TH2F* fTriggerPt;                  //2 dim histogramm with fine binning to describe the pT distribution of the trigger particles

  Double_t fTriggerPt1Min;           //minimum energy for the first trigger particle
  Double_t fTriggerPt1Max;           //maximum energy for the first trigger particle
  Double_t fTriggerPt2Min;           //minimum energy for the second trigger particle
  Double_t fTriggerPt2Max;           //maximum energy for the second trigger particle
  Double_t fPtAssocMin;              //minimum energy for the associated particles
  Double_t fPtAssocMax;              //maximum energy for the associated particle
  Double_t fAlpha;                   //minimum energy for the first trigger particle
  Bool_t fUseLeadingPt;        //decides if all particles of a cone are used as trigger particles or only the leading particles within alpha (apply this on near and away side)
  Bool_t fUseAllT1;                  //use all possible T1 combinations
  Bool_t fUseBackgroundSameOneSide;  //uses only one side for background same
  Bool_t fUseSmallerPtAssoc;         //use only associated particles with less pT than the trigger particles
  THnF* fEfficiencyCorrection;   // if non-0 this efficiency correction is applied on the fly to the filling for trigger particles. The factor is multiplicative, i.e. should contain 1/efficiency
  Int_t fMergeCount;	             // counts how many objects have been merged together
  
  ClassDef(AliTwoPlusOneContainer, 10)  // underlying event histogram container
};


#endif
