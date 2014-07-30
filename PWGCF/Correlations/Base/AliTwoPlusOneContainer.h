#ifndef AliTwoPlusOneContainer_H
#define AliTwoPlusOneContainer_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// data container for 2+1 particle analysis

#include "TNamed.h"
#include "AliUEHist.h"
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
  AliTwoPlusOneContainer(const char* name = "AliTwoPlusOneContainer", const char* binning = 0, Double_t alpha = 0.2);

  virtual ~AliTwoPlusOneContainer();
  
  enum PlotKind {kSameNS = 0, kSameAS, kMixedNS, kMixedAS, kMixedCombNS, kMixedCombAS, k1plus1, kBackgroundSameNS, kBackgroundSameAS};

  void FillCorrelations(Double_t centrality, Float_t zVtx, AliTwoPlusOneContainer::PlotKind step, TObjArray* triggerNear, TObjArray* triggerAway, TObjArray* assocNear, TObjArray* assocAway, Double_t weight, Bool_t is1plus1, Bool_t isBackgroundSame);
  
  AliUEHist* GetData() {return fTwoPlusOne;}
  Double_t getTriggerPt1Min() {return fTriggerPt1Min;}
  Double_t getTriggerPt1Max() {return fTriggerPt1Max;}
  Double_t getTriggerPt2Min() {return fTriggerPt2Min;}
  Double_t getTriggerPt2Max() {return fTriggerPt2Max;}
  Double_t getPtAssocMin() {return fPtAssocMin;}
  Double_t getPtAssocMax() {return fPtAssocMax;}
  
  void SetData(AliUEHist* obj) {fTwoPlusOne = obj; }
  void SetAlpha(Double_t value) {fAlpha = value;}

  AliTwoPlusOneContainer(const AliTwoPlusOneContainer &c);
  AliTwoPlusOneContainer& operator=(const AliTwoPlusOneContainer& c);
  virtual void Copy(TObject& c) const;
  virtual Long64_t Merge(TCollection* list);
 
protected:
  void DeleteContainers();
  
  AliUEHist* fTwoPlusOne;	     //a 6 dim histogram which actually contains all the data
  
  Double_t fTriggerPt1Min;           //minimum energy for the first trigger particle
  Double_t fTriggerPt1Max;           //maximum energy for the first trigger particle
  Double_t fTriggerPt2Min;           //minimum energy for the second trigger particle
  Double_t fTriggerPt2Max;           //maximum energy for the second trigger particle
  Double_t fPtAssocMin;              //minimum energy for the associated particles
  Double_t fPtAssocMax;              //maximum energy for the associated particle
  Double_t fAlpha;                   //minimum energy for the first trigger particle
  Int_t fMergeCount;	             // counts how many objects have been merged together
  
  ClassDef(AliTwoPlusOneContainer, 2)  // underlying event histogram container
};


#endif
