#ifndef ALIMUONPAIRLIGHT_H
#define ALIMUONPAIRLIGHT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 06/09/2006

/// \ingroup evaluation
/// \class AliMUONPairLight
/// \brief Compact information for the generated muon pairs 
/// 
/// Compact information for the generated muon pairs in the MUON arm 
/// useful at the last stage of the analysis chain
/// Pairs are built with two AliMUONTrackLight objects 
/// Using the class AliMUONTrackLight this class combines the decay
/// information ("history") of the reconstructed tracks and fills
/// a series of flags for the formed reconstructed dimuon:
/// fIsCorrelated, fCreationProcess, fIsFeedDown, ...
/// for information about the dimuon, use PrintInfo with the appropriate
/// printflag
/// To be used together with AliMUONTrackLight
///
/// \author This class was prepared by INFN Cagliari, July 2006
/// (authors: H.Woehri, A.de Falco)

// MUON classes
#include "AliMUONTrackLight.h"

// ROOT classes
//#include "TLorentzVector.h"
class TLorentzVector;

class AliMUONPairLight : public TObject { 
public:
  AliMUONPairLight();
  AliMUONPairLight(AliMUONPairLight &dimuCopy);
  virtual ~AliMUONPairLight();
  virtual void SetMuons(AliMUONTrackLight mu0, AliMUONTrackLight mu1);
  AliMUONTrackLight* GetMuon(Int_t index) ; 
  Int_t GetMuonMotherPDG(Int_t imuon, Int_t mother=0) ; 
  Bool_t IsOpenCharm() {return (TMath::Abs(fMu0.GetParentFlavour(0))==4 && TMath::Abs(fMu1.GetParentFlavour(0))==4 && fIsCorrelated && !IsAResonance());}
  Bool_t IsOpenBeauty() {return (TMath::Abs(fMu0.GetParentFlavour(0))==5 && TMath::Abs(fMu1.GetParentFlavour(0))==5 && fIsCorrelated  && !IsAResonance());}
  Bool_t IsAResonance(); 
  Bool_t IsCorrelated() const {return fIsCorrelated;}   
  Int_t GetCauseOfCorrelation() const {return fCauseOfCorrelation;}
  Bool_t IsFeedDown() const {return fIsFeedDown;}   
  Bool_t IsOneTrackNotAMuon() {return (!( fMu0.IsAMuon() && fMu1.IsAMuon() )) ;}
  Int_t GetCharge() {return fMu0.GetCharge() + fMu1.GetCharge();}
  Int_t GetCreationProcess() const  {return fCreationProcess;}
  void SetCorrelated(Bool_t answer) {fIsCorrelated = answer; }
  void SetCauseOfCorrelation(Int_t pdg) {fCauseOfCorrelation = pdg; }
  void SetFeedDown(Int_t answer) {fIsFeedDown = answer;}
  TLorentzVector GetPRec(){return fMu0.GetPRec()+fMu1.GetPRec();}
  TLorentzVector GetPGen(){return fMu0.GetPGen()+fMu1.GetPGen();}
  Double_t GetOpeningAngle(); 
  virtual void PrintInfo(Option_t* opt);//"H" single muons' decay histories
                                           //"K" dimuon kinematics
                                           //"F" dimuon flags
                                           //"A" all variables
protected: 
  AliMUONTrackLight fMu0;   /// first muon  
  AliMUONTrackLight fMu1;   /// second muon
  Int_t fCreationProcess; ///0: pair creation, 1: gluon splitting, 2: flavour excitation, 3: same fragmented mother, -1: resonance
  Bool_t fIsCorrelated; ///tells if the two muons are of correlated origin
  Int_t fCauseOfCorrelation; ///pdg of common mother
  Int_t fIsFeedDown;     /// tells if the process is from feeddown 

  void SetProcess();///checks if muons are correlated and assigns

  ClassDef(AliMUONPairLight,1)  /// Dimuon carrying info at generation and reconstruction level
}; 
#endif
