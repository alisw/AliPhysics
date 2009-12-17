#ifndef ALIAODMUONTRACK_H
#define ALIAODMUONTRACK_H

#include <TObject.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

#include "AliAODTrack.h"
#include "AliESDMuonTrack.h"

class AliAODMuonTrack : public TObject {
 public:

  AliAODMuonTrack();
  AliAODMuonTrack(AliAODTrack *trk);
  AliAODMuonTrack(AliESDMuonTrack *trk);
  virtual ~AliAODMuonTrack();

  Bool_t SelectSingleMuon(Double_t cuts[10]);

  TLorentzVector GetP()       const { return fP;       }
  Int_t          GetCharge()  const { return fCharge;  }
  Int_t          GetTrigger() const { return fTrigger; }
  Double_t       GetDCA()     const { return fDca;     }
  Double_t       GetChi2()    const { return fChi2;    }
  Double_t       GetCentr()   const { return fCentr;   }

 private:

  void FillTrackInfo(AliAODTrack *trk);
  void FillTrackInfo(AliESDMuonTrack *trk);

  TLorentzVector fP;
  Short_t fCharge;
  Int_t fTrigger;
  Double_t fDca;
  Double_t fChi2;
  Double_t fCentr;  // used for PbPb conllisions

  ClassDef(AliAODMuonTrack, 3);
};

#endif
