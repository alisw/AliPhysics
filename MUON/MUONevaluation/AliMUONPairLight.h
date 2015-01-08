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
  AliMUONPairLight& operator=(const AliMUONPairLight&);
  virtual void SetMuons(const AliMUONTrackLight& mu0, const AliMUONTrackLight& mu1);
  AliMUONTrackLight* GetMuon(Int_t index) ;   
  Int_t GetMuonMotherPDG(Int_t imuon, Int_t mother=0) ; 

  /// returns kTRUE if the creation process of the pair was "open charm" (kFALSE... otherwise)
  Bool_t IsOpenCharm() {return (TMath::Abs(fMu0.GetParentFlavour(0))==4 && TMath::Abs(fMu1.GetParentFlavour(0))==4 && fIsCorrelated && !IsAResonance());}
  /// returns kTRUE if the creation process of the pair was "open beauty" (kFALSE... otherwise)
  Bool_t IsOpenBeauty() {return (TMath::Abs(fMu0.GetParentFlavour(0))==5 && TMath::Abs(fMu1.GetParentFlavour(0))==5 && fIsCorrelated  && !IsAResonance());}
  Bool_t IsAResonance(); 
  /// returns kTRUE if at least one of the first hadronised parent is a pi or a K (kFALSE... otherwise)
  Bool_t IsOneMuonFromPionOrKaon(){return (fMu0.IsParentPionOrKaon(0) || fMu1.IsParentPionOrKaon(0));}
  /// Return the info if the two muons are of correlated origin
  Bool_t IsCorrelated() const {return fIsCorrelated;}   
  /// Return the pdg of common mother
  Int_t GetCauseOfCorrelation() const {return fCauseOfCorrelation;}
  /// Return the info if the process is from feeddown 
  Bool_t IsFeedDown() const {return fIsFeedDown;}   
  /// returns kTRUE if at least one of the reconstructed tracks is not a muon (kFALSE... otherwise)
  Bool_t IsOneTrackNotAMuon() {return (!( fMu0.IsAMuon() && fMu1.IsAMuon() )) ;}
  /// returns the charge of the created pair
  Int_t GetCharge() {return fMu0.GetCharge() + fMu1.GetCharge();}
  /// \brief Return the info ablout creation process
  ///0: pair creation, 1: gluon splitting, 2: flavour excitation, 3: same fragmented mother, -1: resonance
  Int_t GetCreationProcess() const  {return fCreationProcess;}
  /// Set the info ablout creation process
  void SetCorrelated(Bool_t answer) {fIsCorrelated = answer; }
  /// Set the pdg of common mother
  void SetCauseOfCorrelation(Int_t pdg) {fCauseOfCorrelation = pdg; }
  /// Set the info if the process is from feeddown
  void SetFeedDown(Int_t answer) {fIsFeedDown = answer;}
  /// returns a TLorentzVector containing the reconstructed kinematics of the pair
  TLorentzVector GetPRec(){return fMu0.GetPRec()+fMu1.GetPRec();}
  /// returns a TLorentzVector containing the generated kinematics of the pair
  TLorentzVector GetPGen(){return fMu0.GetPGen()+fMu1.GetPGen();}
  Double_t GetOpeningAngle(); 
  Bool_t IsDimuonFromCorrPiK();
  virtual void PrintInfo(const Option_t* opt);//"H" single muons' decay histories
                                           //"K" dimuon kinematics
                                           //"F" dimuon flags
                                           //"A" all variables

protected: 
  /// Checks if muons are correlated and assigns
  void SetProcess();

  AliMUONTrackLight fMu0;   ///< first muon  
  AliMUONTrackLight fMu1;   ///< second muon
  Int_t fCreationProcess; ///<0: pair creation, 1: gluon splitting, 2: flavour excitation, 3: same fragmented mother, -1: resonance
  Bool_t fIsCorrelated; ///<tells if the two muons are of correlated origin
  Int_t fCauseOfCorrelation; ///<pdg of common mother
  Int_t fIsFeedDown;     ///< tells if the process is from feeddown 

  ClassDef(AliMUONPairLight,1)  /// Dimuon carrying info at generation and reconstruction level
}; 
#endif
