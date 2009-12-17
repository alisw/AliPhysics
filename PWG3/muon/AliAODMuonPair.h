#ifndef ALIAODMUONPAIR_H
#define ALIAODMUONPAIR_H

#include <TObject.h>
#include <TRef.h>
#include <TLorentzVector.h>

#include "AliAODMuonTrack.h"

class AliAODMuonPair : public TObject {
 public:

  AliAODMuonPair();
  AliAODMuonPair(AliAODMuonTrack *trk0, AliAODMuonTrack *trk1);
  virtual ~AliAODMuonPair();

  AliAODMuonTrack* GetTrack(Int_t i) const { return (i<2 ? (AliAODMuonTrack*)(fTrk[i].GetObject()) : 0x0); }

  TLorentzVector GetP()                     const { return fP;                           }
  Int_t          GetCharge()                const { return fCharge;                      }
  Int_t          GetMuMatchTrigger(Int_t i) const { return (i<2 ? fTrigger[i] : -9999 ); }
  Double_t       GetMuDCA(Int_t i)          const { return (i<2 ? fDca[i]     : -9999.); }
  Double_t       GetMuChi2(Int_t i)         const { return (i<2 ? fChi2[i]    : -9999.); }

 protected:

  TRef fTrk[2];
  void FillPairInfo();

 private:


  TLorentzVector fP;
  Int_t fCharge;
  Int_t fTrigger[2];
  Double_t fDca[2];
  Double_t fChi2[2];

  ClassDef(AliAODMuonPair, 3);
};

#endif
