#ifndef ALIMUONINFOSTOREMC_H
#define ALIMUONINFOSTOREMC_H

/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

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
#include "AliMCEvent.h"

class AliMuonInfoStoreMC : public AliMuonInfoStoreRD {
 public:

  AliMuonInfoStoreMC();
  AliMuonInfoStoreMC(AliAODTrack     *trkAOD, AliMCEvent *mcEvent, Bool_t full=kFALSE);
  AliMuonInfoStoreMC(AliESDMuonTrack *trkESD, AliMCEvent *mcEvent, Bool_t full=kFALSE);
  AliMuonInfoStoreMC(const AliMuonInfoStoreMC &src);
  AliMuonInfoStoreMC& operator=(const AliMuonInfoStoreMC &src);
  virtual ~AliMuonInfoStoreMC();

  Int_t  ParentFlavour(Int_t i=0)    const;
  Bool_t IsMotherAResonance(Int_t i) const;

  TLorentzVector LorentzP()         const { return fLorentzP; }
  Int_t    Source()                 const { return fSource; }
  Int_t    TrackIndex()             const { return fTrackIndex; }
  Int_t    TrackPDGCode()           const { return fTrackPDGCode; }
  Int_t    ParentsN()               const { return fNParents; }
  Int_t    ParentIndex(Int_t i=0)   const { return (i<fNParents ? fParentIndex[i] : -1); }
  Int_t    ParentPDGCode(Int_t i=0) const { return (i<fNParents ? fParentPDGCode[i] : 0); }
  Int_t    QuarkIndex(Int_t i=0)    const { return (i<4 ? fQuarkIndex[i] : -1); }
  Int_t    QuarkPDGCode(Int_t i=0)  const { return (i<4 ? fQuarkPDGCode[i] : 0); }
  Bool_t   IsOscillation()          const { return fOscillation; }
  Double_t Weight()                 const { return fWeight; }

  static const char* StdBranchName() { return fgkStdBranchName.Data(); }
  static Int_t SourcesN()            { return fgkSourcesN;             }

 private:

  void SetMCInfoAOD(AliMCEvent *mcEvent, Int_t label);
  void SetMCInfoESD(AliMCEvent *mcEvent, Int_t label);
  void FillHistoryQuarksAOD(AliMCEvent *mcEvent, Int_t lineM);
  void FillHistoryQuarksESD(AliMCEvent *mcEvent, Int_t lineM);
  Int_t SelectHFMuon();

  Bool_t IsDiquark(Int_t pdg);
  void ResetQuarkInfo();

  static const TString fgkStdBranchName;  // Standard branch name
  static const Int_t   fgkSourcesN;       // num. of muon sources

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

  ClassDef(AliMuonInfoStoreMC, 5);
};

#endif
