#ifndef ALIUNICOREVENTALICEESD_H
#define ALIUNICOREVENTALICEESD_H

/* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2008

#include <cmath>
#include "TVector2.h"
#include "AliESDEvent.h"
#include "AliUnicorEvent.h"

//=============================================================================
class AliUnicorEventAliceESD : public AliUnicorEvent {

 public:
              AliUnicorEventAliceESD(AliESDEvent *esd=0);
              AliUnicorEventAliceESD(const AliUnicorEventAliceESD &ev): AliUnicorEvent(ev), fESD(ev.fESD) {}
  virtual     ~AliUnicorEventAliceESD();
  AliUnicorEventAliceESD &operator=(const AliUnicorEventAliceESD &source) {fESD = source.fESD; return *this;}
  Double_t    Etamin() const {return -0.75;}
  Double_t    Etamax() const {return  0.75;}
  void        AttachTree(TTree *tr) {fESD->ReadFromTree(tr);}
  Bool_t      Good() const;
  Double_t    Centrality() const {return 0.9999*exp(-NParticles()/20.0);} // OK for pp
  void        RP(Double_t &qx, Double_t &qy) const {AliUnicorEvent::RP(qx,qy,2);}
  Double_t    RPphi() const {Double_t qx,qy; RP(qx,qy); return atan2(qy,qx);}
  Double_t    Zver() const {return fESD->GetPrimaryVertexTPC()->GetZv()/10.0;}
  Int_t       NParticles() const {return fESD->GetNumberOfTracks();}

  Bool_t      ParticleGood(Int_t i, Int_t pidi=0) const;
  Double_t    ParticleP(Int_t i)     const {return fESD->GetTrack(i)->GetTPCInnerParam()->P();}
  Double_t    ParticleTheta(Int_t i) const {return fESD->GetTrack(i)->GetTPCInnerParam()->Theta();}
  Double_t    ParticlePhi(Int_t i)   const {return TVector2::Phi_mpi_pi(fESD->GetTrack(i)->GetTPCInnerParam()->Phi());}
  Double_t    ParticleDedx(Int_t i)  const {return fESD->GetTrack(i)->GetTPCsignal()/50.0;}
  Bool_t      PairGood(Double_t p0, Double_t the0, Double_t phi0, 
		       Double_t p1, Double_t the1, Double_t phi1) const;
  // alternative: GetConstrainedParam, GetInnerParam, GetTPCInnerParam 
  void        SetESD(AliESDEvent * const esd) {fESD = esd;}
  AliESDEvent *GetESD() const {return fESD;}

 protected:
  AliESDEvent *fESD;   //! pointer to the actual source of data

  ClassDef(AliUnicorEventAliceESD,0)
};
#endif 
//=============================================================================
