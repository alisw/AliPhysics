/* This file is property of and copyright                                 *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliAnalysisTaskGammaJet.h
/// @author Svein Lindal
/// @brief  Class used to run isolation studies of conversion gamma/pions
 

#ifndef ALIANACONVISOLATION_CXX
#define ALIANACONVISOLATION_CXX

#include "TObject.h"
#include "Rtypes.h"
#include "TF1.h"
class TH2F;
class TH1F;
class AliAODConversionPhoton;
class TClonesArray;
class TList;

class AliAnaConvIsolation : public TObject {

public:
  
  AliAnaConvIsolation(); 
  AliAnaConvIsolation(Float_t coneSize, Float_t maxPtThreshold, Float_t sumPtThreshold, Float_t maxPtFraction, Float_t sumPtFraction);
  virtual ~AliAnaConvIsolation();
  

  void CreateHistograms();


  ///Set And get cone size
  void SetConeSize(Float_t cs) {fConeSize = cs;}
  Float_t GetConeSize() const {return fConeSize;}


  //Set and get max pt threshold
  void SetMaxPtThreshold(Float_t cs) {fMaxPtThreshold = cs;}
  Float_t GetPtThreshold() const {return fMaxPtThreshold;}


  //Set and get sum pt threshold
  void SetSumPtThreshold(Float_t cs) {fSumPtThreshold = cs;}
  Float_t GetPtSumThreshold() const {return fSumPtThreshold;}

  //Set and get max Pt fraction threshold
  void SetMaxPtFraction(Float_t cs) {fMaxPtFraction = cs;}
  Float_t GetPtFraction() const {return fMaxPtFraction;}

  //Set and get sum pt fraction threshold
  void SetSumPtFraction(Float_t cs) {fSumPtFraction = cs;}
  Float_t GetPtSumFraction() const {return fSumPtFraction;}

  //Set min pt for a particle to be counted in bg
  void SetMinPt( Float_t minpt) {fMinPt = minpt; }

  //Get isolation curve
  TF1 * GetIsolationCurve() const { return fIsoCurve; }

  //Set function of isolation curve
  void SetIsolationCurve(TString  curve) { 
    fCurveFunction = curve;
  }

  //Fill histograms
  void FillHistograms(Float_t pt, Float_t ptMax, Float_t ptSum, Bool_t isolated, Int_t nTracks);

  //Get list of histograms
  TList * GetHistograms() const { return fHistograms;}

  //Is particle isolated
  Bool_t IsIsolated( const AliAODConversionPhoton * const particle, const TClonesArray * const tracks, Bool_t &leading);
  Bool_t IsIsolated( AliAODConversionPhoton * const particle, const TClonesArray * const tracks, const Int_t nSpawn, const Int_t * const spawn, Bool_t &leading );

 private:

  //Is eta - phi distance smaller than conesize ?
  Bool_t IsInCone(Float_t dEta, Float_t dPhi, Float_t coneSize) const {
    return ( (dEta*dEta + dPhi*dPhi) < coneSize*coneSize);
  }

  ///Evaluate whether particle is isolated according to criterie
  Bool_t EvaluateIsolationCriteria(Float_t ptSum, Float_t pt) const;

  TF1 * fIsoCurve; ///Curve defining if particle is isolated or not 
  TString fCurveFunction; ///Funtion defining curve

  Float_t fConeSize; //Size of isolation cone
  Float_t fMinPt; //min pt for bg particles
  Float_t fMaxPtThreshold; //max pt threshold
  Float_t fSumPtThreshold; //sum pt threhold
  Float_t fMaxPtFraction;  //max pt fraction threshold
  Float_t fSumPtFraction; //sum pt fraction threshold

  TList * fHistograms; //list of histograms

  TH2F * fhMaxPtInCone[2]; //histogram of max pt in cone
  TH2F * fhSumPtInCone[2]; //histogram of sum pt in cone
  TH2F * fhSumPtVsMaxPt[2];//sum pt vs max pt
  
  // TH1F * fHistSumPt[2][2]; 
  // TH1F * fHistMaxPt[2][2];
  
  TH1F * fhPtCandidates[2]; //pt distribution of isolation candidates
  TH1F * fhTrackMult[2];    //Track multiplicity of events with / wo isolated particles
  

  Float_t fHistogramMaxPt; //Upper pt limit in histograms


  AliAnaConvIsolation(const AliAnaConvIsolation&); // not implemented
  AliAnaConvIsolation& operator=(const AliAnaConvIsolation&); // not implemented
  ClassDef(AliAnaConvIsolation, 2); // example of analysis
};

#endif
