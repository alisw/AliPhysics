#ifndef ALIEMCALPARTICLE_H
#define ALIEMCALPARTICLE_H

// $Id$

#include <TLorentzVector.h>
#include <TMath.h>
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"

class AliEmcalParticle: public AliVParticle {
 public:
  AliEmcalParticle();
  AliEmcalParticle(TObject *particle, Int_t id = -1, Double_t vx=0, Double_t vy=0, Double_t vz=0);
  AliEmcalParticle(const AliEmcalParticle &p); 
  AliEmcalParticle &operator=(const AliEmcalParticle &p);
  virtual ~AliEmcalParticle();

  // AliVParticle interface
  Double_t          Px()        const { return fPt*TMath::Cos(fPhi);  }
  Double_t          Py()        const { return fPt*TMath::Sin(fPhi);  };
  Double_t          Pz()        const { return fPt*TMath::SinH(fEta); }
  Double_t          Pt()        const { return fPt ; }
  Double_t          P()         const { return TMath::Sqrt(Px()*Px()+Py()*Py()+Pz()*Pz()); }
  Bool_t            PxPyPz(Double_t p[3]) const { p[0] = Px(); p[1] = Py(); p[2] = Pz(); return 1; }
  Double_t          Xv()        const { return 0; }
  Double_t          Yv()        const { return 0; }
  Double_t          Zv()        const { return 0; }
  Bool_t            XvYvZv(Double_t x[3]) const { x[0] = Xv(); x[1] = Yv(); x[2] = Zv(); return 1; }
  Double_t          OneOverPt() const { return 1./fPt; }
  Double_t          Phi()       const { return fPhi; }
  Double_t          Theta()     const { return 0.; }
  Double_t          E()         const { if (fTrack) return fTrack->E(); return fCluster->E(); }
  Double_t          M()         const { if (fTrack) return fTrack->M(); return 0; }
  Double_t          Eta()       const { return fEta; }
  Double_t          Y()         const { if (fTrack) return fTrack->Y(); return fEta; }
  Short_t           Charge()    const { if (fTrack) return fTrack->Charge(); else return 0; }
  Int_t   GetLabel()    const { if (fTrack) return fTrack->GetLabel(); return fCluster->GetLabel(); }
  Int_t   PdgCode()     const { return 0; }
  const Double_t *PID() const { return 0; }

  AliVCluster*      GetCluster()           const { return fCluster                                                  ; }
  Int_t             GetMatchedObjId(UShort_t i = 0)         const { return fNMatched > i ? fMatchedIds[i]  : -1     ; }
  Double_t          GetMatchedObjDistance(UShort_t i = 0)   const { return fNMatched > i ? fMatchedDist[i] : -1     ; }
  UShort_t          GetNumberOfMatchedObj()                 const { return fNMatched                                ; }
  AliVTrack*        GetTrack()             const { return fTrack                                                    ; }
  Double_t          GetTrackPhiOnEMCal()   const { if (fTrack) return fTrack->GetTrackPhiOnEMCal(); else return -999; }
  Double_t          GetTrackEtaOnEMCal()   const { if (fTrack) return fTrack->GetTrackEtaOnEMCal(); else return -999; }
  Int_t             IdInCollection()       const { return fId                                                       ; }
  Bool_t            IsCluster()            const { return (Bool_t) fCluster != 0                                    ; }
  Bool_t            IsEMCAL()              const { if (fCluster) return kTRUE; 
                                                   if (fTrack)   return fTrack->IsEMCAL(); 
                                                   return kFALSE; }
  Bool_t            IsTrack()              const { return (Bool_t) fTrack   != 0                                    ; }
  Bool_t            IsMC(Int_t minLabel=0) const { if (fTrack) return (TMath::Abs(fTrack->GetLabel()) > minLabel); 
                                                   return (fCluster->GetLabel() > minLabel); }

  void              AddMatchedObj(Int_t id, Double_t d);
  void              ResetMatchedObjects();
  void              SetIdInCollection(Int_t id)          { fId = id                                                                      ; }
  void              SetMatchedObj(Int_t id, Double_t d)  { ResetMatchedObjects(); fMatchedIds[0] = id; fMatchedDist[0] = d; fNMatched = 1; }

 protected:
  TLorentzVector   &GetLorentzVector(const Double_t *vertex = 0)  const;

  static const UShort_t fSizeMatched = 99;        //!size of matched clusters array

  AliVTrack        *fTrack;                       //!track
  AliVCluster      *fCluster;                     //!cluster
  UShort_t          fMatchedIds[fSizeMatched];    //!ids of matched clusters, ordered from the closest to the farthest
  Double_t          fMatchedDist[fSizeMatched];   //!distances of matched clusters
  UShort_t          fNMatched;                    //!number of matched objects 
  Int_t             fId;                          //!id in original collection
  Double_t          fPhi;                         //!phi
  Double_t          fEta;                         //!eta
  Double_t          fPt;                          //!pt

  ClassDef(AliEmcalParticle, 1) // Emcal particle class
};
#endif
