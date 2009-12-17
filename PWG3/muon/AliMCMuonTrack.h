#ifndef ALIMCMUONTRACK_H
#define ALIMCMUONTRACK_H

#include <TParticle.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliAODTrack.h"
#include "AliAODMuonTrack.h"

class AliMCMuonTrack : public AliAODMuonTrack {
 public:

  AliMCMuonTrack();
  AliMCMuonTrack(AliAODTrack *trkAOD, TClonesArray *mcClArr, Bool_t full=kFALSE);
  AliMCMuonTrack(AliESDMuonTrack *trkESD, AliESDEvent *esd, AliMCEventHandler *mcH, Bool_t full=kFALSE);
  AliMCMuonTrack(AliESDMuonTrack *trkESD, AliStack *stack);
  virtual ~AliMCMuonTrack();

  Int_t GetParentFlavour(Int_t i=0) const;
  Bool_t IsMotherAResonance(Int_t i) const;
  Bool_t IsAMuon() const { return (fSource>=0 && fSource!=4); }

  TLorentzVector GetPGen() const { return fPGen; }
  Int_t GetSource() const { return fSource; }
  Int_t GetTrackIndex() const { return fTrackIndex; }
  Int_t GetTrackPDGCode() const { return fTrackPDGCode; }
  Int_t GetNParents() const { return fNParents; }
  Int_t GetParentIndex(Int_t i=0) const { return (i<fNParents ? fParentIndex[i] : -1); }
  Int_t GetParentPDGCode(Int_t i=0) const { return (i<fNParents ? fParentPDGCode[i] : 0); }
  Int_t GetQuarkIndex(Int_t i=0) const { return (i<4 ? fQuarkIndex[i] : -1); }
  Int_t GetQuarkPDGCode(Int_t i=0) const { return (i<4 ? fQuarkPDGCode[i] : 0); }
  Bool_t IsOscillation() const { return fOscillation; }
  Double_t GetWeight() const { return fWeight; }

 private:

  AliAODMCParticle* FindTrackRef(AliAODTrack *trkAOD, TClonesArray *mcClArr);
  TParticle* FindTrackRef(AliESDMuonTrack *trkESD, AliESDEvent *esd, AliMCEventHandler *mcH);  // do not implement on official train
  void SetMCInfo(AliAODMCParticle *pMC, TClonesArray *mcClArr);
  void SetMCInfo(TParticle *pMC, AliMCEventHandler *mcH);
  void FillHistoryParents(Int_t lineM, TClonesArray *mcClArr);
  void FillHistoryParents(Int_t lineM, AliStack *stack);
  void FillHistoryQuarks(Int_t lineM, TClonesArray *mcClArr);
  void FillHistoryQuarks(Int_t lineM, AliStack *stack);
  //AliMUONTrack* CovESDtoMuonTrack(AliESDMuonTrack &trkESD);  // do not implement on official train
  Int_t SelectHFMuon();

  Bool_t IsDiquark(Int_t pdg);
  void ResetQuarkInfo();

  static const Double_t fgkSigmaCut;

  Bool_t fIsFull;
  static const Int_t fgkNParentsMax = 5;
  TLorentzVector fPGen;
  Int_t fTrackIndex;
  Int_t fTrackPDGCode;
  Int_t fSource;  // = 0, mu<-b 
                  // = 1, mu<-c 
                  // = 2, primary mu
                  // = 3, secondary mu
                  // = 4, not mu

  Int_t fParentIndex[fgkNParentsMax];
  Int_t fParentPDGCode[fgkNParentsMax];
  Int_t fNParents;
  Bool_t fOscillation;

  Int_t fQuarkIndex[4];
  Int_t fQuarkPDGCode[4];

  Double_t fWeight;  // for PbPb collisoions

  ClassDef(AliMCMuonTrack, 2);
};

#endif
