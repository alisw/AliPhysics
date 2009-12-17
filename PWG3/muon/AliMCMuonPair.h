#ifndef ALIMCMUONPAIR_H
#define ALIMCMUONPAIR_H

#include <TLorentzVector.h>

#include "AliMCMuonTrack.h"
#include "AliAODMuonPair.h"

class AliMCMuonPair : public AliAODMuonPair {
 public:

  AliMCMuonPair();
  AliMCMuonPair(AliMCMuonTrack *trk0, AliMCMuonTrack *trk1, Bool_t full=kFALSE);
  virtual ~AliMCMuonPair();

  AliMCMuonTrack* GetTrack(Int_t i) const { return (i<2 ? (AliMCMuonTrack*)(fTrk[i].GetObject()) : 0x0); }

  TLorentzVector GetPGen() const { return fPGen; }
  Int_t GetSource() const { return fSource; }

 private:

  void FindDimuonSourceFast();
  void FindDimuonSourceFull();

  Bool_t fIsFull;
  TLorentzVector fPGen;
  Int_t fSource;  // = 0, BBdiff
                  // = 1, Bchain
                  // = 2, DDdiff
                  // = 3, Dchain
                  // = 4, Resonance
                  // = 5, UnCorr bkg

  ClassDef(AliMCMuonPair, 1);
};

#endif
