#ifndef AliAODEventInfo_H
#define AliAODEventInfo_H

/* $Id$ */ 

#include "TClonesArray.h"

#include "AliAODEvent.h"

class AliAODHeader;

// AliAODEventInfo: a class for AODs for the MUON Arm of the ALICE Experiment
// Author: P. Cortese, Universita' del Piemonte Orientale in Alessandria and
// INFN of Torino - Italy
//
// This class provides additional information about the AliAODEvent, some
// of this is specific to MUON. The information stored in this class will 
// follow the evolution of the framework (i.e. some data members  may be 
// moved into the header in the future).
//
//

// 2007/11/06 v1.00 Initial version
// 2007/12/18 v1.01 More compact information for what regards trigger info
// 2008/02/01 v1.02 Apply coding conventions

class AliAODEventInfo : public TNamed {
public:
  AliAODEventInfo();
  ~AliAODEventInfo();

  void SetBeamEnergy(Double_t BeamEnergy){ fBeamEnergy=BeamEnergy; };
  Double_t EBeam() const { return fBeamEnergy; };
  Double_t SqrtS() const { return 2*fBeamEnergy; };

  // -- Trigger information (to be updated in future)
  void SelectTriggerBits(UChar_t muonSingleLPtL0, UChar_t muonSingleHPtL0,
			 UChar_t muonLikeLPtL0, UChar_t muonLikeHPtL0,
			 UChar_t muonUnlikeLPtL0, UChar_t muonUnlikeHPtL0);

  UChar_t GetBitSingleLPtL0() const { return fMuonSingleLPtL0; }
  UChar_t GetBitSingleHPtL0() const { return fMuonSingleHPtL0; }
  UChar_t GetBitLikeLPtL0() const { return fMuonLikeLPtL0; }
  UChar_t GetBitLikeHPtL0() const { return fMuonLikeHPtL0; }
  UChar_t GetBitUnlikeLPtL0() const { return fMuonUnlikeLPtL0; }
  UChar_t GetBitUnlikeHPtL0() const { return fMuonUnlikeHPtL0; }

  Bool_t MuonSingleLPtL0() const; // Test trigger pattern for MUON_Single_LPt_L0
  Bool_t MuonSingleHPtL0() const; // Test trigger pattern for MUON_Single_HPt_L0
  Bool_t MuonLikeLPtL0() const; // Test trigger pattern for MUON_Like_LPt_L0
  Bool_t MuonLikeHPtL0() const; // Test trigger pattern for MUON_Like_HPt_L0
  Bool_t MuonUnlikeLPtL0() const; // Test trigger pattern for MUON_Unlike_LPt_L0
  Bool_t MuonUnlikeHPtL0() const; // Test trigger pattern for MUON_Unlike_HPt_L0

  Int_t GetNDimuons() const { return (fDi!=0) ? ((TClonesArray*)fDi.GetObject())->GetSize() : 0;}
  Int_t NDimu() const { return GetNDimuons();}

  // Pointers
  void SetEv(AliAODEvent *ev){ fEv=ev; }
  void SetEi(AliAODEventInfo *ei){ fEi=ei; }
  void SetHe(AliAODHeader *he){ fHe=he; }
  void SetTr(TClonesArray *tr){ fTr=tr; }
  void SetDi(TClonesArray *di){ fDi=di; }

  AliAODEvent *GetEv() { return (fEv!=0) ? (AliAODEvent*)fEv.GetObject() : 0; }
  AliAODEventInfo* GetEi() { return (fEi!=0) ? (AliAODEventInfo*)fEi.GetObject() : 0; }
  TClonesArray *GetDi() { return (fDi!=0) ? (TClonesArray*)fDi.GetObject() : 0; } // Get dimuon array

  AliAODEvent *Ev() { return (fEv!=0) ? (AliAODEvent*)fEv.GetObject() : 0; }
  TClonesArray *Di() { return (fDi!=0) ? (TClonesArray*)fDi.GetObject() : 0; } // Get dimuon array
  AliAODEventInfo* Ei() { return (fEi!=0) ? (AliAODEventInfo*)fEi.GetObject() : 0; }

//  AliAODHeader *GetHe() { return (fHe!=0) ? (AliAODHeader*)fHe.GetObject() : 0; }
//  AliAODTrack *GetTr();

  Bool_t IsHeaderAccessible(const Char_t *msg=0) const;

  protected:
  // Missing in AliAODHeader and added here
  Double_t fBeamEnergy; // Add beam energy not present in AliAODHeader
  UChar_t fMuonSingleLPtL0; // Decode trigger info
  UChar_t fMuonSingleHPtL0; // Decode trigger info
  UChar_t fMuonLikeLPtL0;   // Decode trigger info
  UChar_t fMuonLikeHPtL0;   // Decode trigger info
  UChar_t fMuonUnlikeLPtL0; // Decode trigger info
  UChar_t fMuonUnlikeHPtL0; // Decode trigger info

  // Data members to provide automatic access
  TRef fEv; // Event
  TRef fEi; // EventInfo
  TRef fHe; // Header
  TRef fTr; // Tracks
  TRef fDi; // Dimuons

  ClassDef(AliAODEventInfo,1)  // Additional header for MUON arm
};

#endif
