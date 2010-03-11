#ifndef ALIMUONINFOSTOREMC_H
#define ALIMUONINFOSTOREMC_H

/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliMuonInfoStoreRD
// class used to extract and store info of MC particle
// Author: X-M. Zhang, zhang@clermont.in2p3.fr
//                     zhangxm@iopp.ccnu.edu.cn
//***********************************************************

#include <TString.h>
#include <TParticle.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliAODTrack.h"
#include "AliMuonInfoStoreRD.h"

class AliMuonInfoStoreMC : public AliMuonInfoStoreRD {
 public:

  AliMuonInfoStoreMC();
  AliMuonInfoStoreMC(AliAODTrack *trkAOD, TClonesArray *mcClArr, Bool_t full=kFALSE);
  AliMuonInfoStoreMC(AliESDMuonTrack *trkESD, AliESDEvent *esd, AliMCEventHandler *mcH, Bool_t full=kFALSE);
  AliMuonInfoStoreMC(AliESDMuonTrack *trkESD, AliMCEventHandler *mcH, Bool_t full=kFALSE);
  AliMuonInfoStoreMC(AliESDMuonTrack *trkESD, AliStack *stack);
  AliMuonInfoStoreMC(const AliMuonInfoStoreMC &src);
  AliMuonInfoStoreMC& operator=(const AliMuonInfoStoreMC &src);
  virtual ~AliMuonInfoStoreMC();

  Int_t ParentFlavour(Int_t i=0) const;
  Bool_t IsMotherAResonance(Int_t i) const;
  Bool_t IsAMuon() const { return (fSource>=0 && fSource!=4); }

  TLorentzVector LorentzP()      const { return fLorentzP; }
  Int_t MuonSource()             const { return fSource; }
  Int_t TrackIndex()             const { return fTrackIndex; }
  Int_t TrackPDGCode()           const { return fTrackPDGCode; }
  Int_t NParents()               const { return fNParents; }
  Int_t ParentIndex(Int_t i=0)   const { return (i<fNParents ? fParentIndex[i] : -1); }
  Int_t ParentPDGCode(Int_t i=0) const { return (i<fNParents ? fParentPDGCode[i] : 0); }
  Int_t QuarkIndex(Int_t i=0)    const { return (i<4 ? fQuarkIndex[i] : -1); }
  Int_t QuarkPDGCode(Int_t i=0)  const { return (i<4 ? fQuarkPDGCode[i] : 0); }
  Bool_t IsOscillation()         const { return fOscillation; }
  Double_t Weight()              const { return fWeight; }

  static const char* StdBranchName() { return fgkStdBranchName.Data(); }
  const static Int_t NSources()      { return fgkNSources;             }

 private:

  AliAODMCParticle* FindTrackRef(AliAODTrack*     const trkAOD, TClonesArray* const mcClArr);
  TParticle*        FindTrackRef(AliESDMuonTrack* const trkESD, AliMCEventHandler* const mcH);
  TParticle*        FindTrackRef(AliESDMuonTrack* const trkESD, AliESDEvent* const esd, AliMCEventHandler* const mcH);
  void SetMCInfo(AliAODMCParticle *pMC, TClonesArray* const mcClArr);
  void SetMCInfo(TParticle *pMC, AliMCEventHandler *mcH);
  void FillHistoryParents(Int_t lineM, TClonesArray* const mcClArr);
  void FillHistoryParents(Int_t lineM, AliStack* const stack);
  void FillHistoryQuarks(Int_t lineM, TClonesArray *mcClArr);
  void FillHistoryQuarks(Int_t lineM, AliStack *stack);
  Int_t SelectHFMuon();

  Bool_t IsDiquark(Int_t pdg);
  void ResetQuarkInfo();

  static const TString fgkStdBranchName;  // Standard branch name
  static const Int_t   fgkNSources;       // num. of muon sources

  Bool_t fIsFull;            // whether to use full mode (Pb-Pb)
  TLorentzVector fLorentzP;  // lorentz momentum of particle
  Int_t fTrackIndex;         // index of the MC particle
  Int_t fTrackPDGCode;       // PDG code of the MC particle
  Int_t fSource;  // = 0, mu<-b 
                  // = 1, mu<-c 
                  // = 2, primary mu
                  // = 3, secondary mu
                  // = 4, not mu
                  // = 5, unidentified track

  Int_t fParentIndex[5];    // index of parents
  Int_t fParentPDGCode[5];  // PDG code of parents
  Int_t fNParents;          // num. of parents
  Bool_t fOscillation;      // flag of oscillation

  Int_t fQuarkIndex[4];    // index of quarks
  Int_t fQuarkPDGCode[4];  // PDG code of quarks

  Double_t fWeight;  // for PbPb collisoions

  ClassDef(AliMuonInfoStoreMC, 3);
};

#endif
