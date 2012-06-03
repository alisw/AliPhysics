#ifndef ALIEMCALPARTICLE_H
#define ALIEMCALPARTICLE_H

// $Id$

#include <TLorentzVector.h>
#include <TObject.h>

#include "AliVTrack.h"

class AliVCluster;

class AliEmcalParticle: public TObject {
 public:
  AliEmcalParticle();
  AliEmcalParticle(TObject *particle, Int_t id = -1);
  ~AliEmcalParticle();
  AliEmcalParticle(const AliEmcalParticle &p); 
  AliEmcalParticle &operator=(const AliEmcalParticle &p);

  Short_t           Charge()               const { if (fTrack) return fTrack->Charge(); else return 0               ; }
  Double_t          GetTrackPhiOnEMCal()   const { if (fTrack) return fTrack->GetTrackPhiOnEMCal(); else return -999; }
  Double_t          GetTrackEtaOnEMCal()   const { if (fTrack) return fTrack->GetTrackEtaOnEMCal(); else return -999; }
  Bool_t            IsEMCAL()              const { if (fCluster) return kTRUE; 
                                                   if (fTrack)   return fTrack->IsEMCAL(); 
                                                   return kFALSE; }
  AliVCluster*      GetCluster()           const { return fCluster                                                  ; }
  AliVTrack*        GetTrack()             const { return fTrack                                                    ; }
  Bool_t            IsCluster()            const { return (Bool_t) fCluster != 0                                    ; }
  Bool_t            IsTrack()              const { return (Bool_t) fTrack   != 0                                    ; }
  Int_t             IdInCollection()       const { return fId                                                       ; }
  Double_t          M()                    const { return 0.13957                                                   ; }

  Int_t             GetMatchedObjId(UShort_t i = 0)         const { return fNMatched > i ? fMatchedIds[i]  : -1 ; }
  Double_t          GetMatchedObjDistance(UShort_t i = 0)   const { return fNMatched > i ? fMatchedDist[i] : -1 ; }
  UShort_t          GetNumberOfMatchedObj()                 const { return fNMatched                            ; }

  TLorentzVector   &GetLorentzVector(const Double_t *vertex = 0)  const;

  void              AddMatchedObj(Int_t id, Double_t d);
  void              ResetMatchedObjects();
  void              SetIdInCollection(Int_t id)          { fId = id                                                                      ; }
  void              SetMatchedObj(Int_t id, Double_t d)  { ResetMatchedObjects(); fMatchedIds[0] = id; fMatchedDist[0] = d; fNMatched = 1; }

 protected:
  static const UShort_t fSizeMatched = 10;        //!size of matched clusters array

  AliVTrack        *fTrack;                       //!track
  AliVCluster      *fCluster;                     //!cluster
  UShort_t          fMatchedIds[fSizeMatched];    //!ids of matched clusters, ordered from the closest to the farthest
  Double_t          fMatchedDist[fSizeMatched];   //!distances of matched clusters
  UShort_t          fNMatched;                    //!number of matched objects 
  Int_t             fId;                          //!id in collection

  ClassDef(AliEmcalParticle, 1) // Emcal particle class
};
#endif
