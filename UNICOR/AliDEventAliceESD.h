// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2008

#ifndef ALIDEVENTALICEESD_H
#define ALIDEVENTALICEESD_H

#include <cmath>
#include "TVector2.h"
#include "AliDEvent.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"

//=============================================================================
class AliDEventAliceESD : public AliDEvent, public AliESDEvent {

 public:
              AliDEventAliceESD();
  virtual    ~AliDEventAliceESD();
  Double_t    Etamin() const {return -0.75;}
  Double_t    Etamax() const {return  0.75;}
  void        AttachTree(TTree *tr) {ReadFromTree(tr);}
  Bool_t      Good() const;
  Double_t    Centrality() {return 0.9999*exp(-NParticles()/20.0);} // OK for pp
  void        RP(Double_t &qx, Double_t &qy) const {AliDEvent::RP(qx,qy,2);}
  Double_t    RPphi() const {Double_t qx,qy; RP(qx,qy); return atan2(qy,qx);}
  Double_t    Zver() const {return AliESDEvent::GetPrimaryVertex()->GetZv()/10.0;}
  Int_t       NParticles() const {return AliESDEvent::GetNumberOfTracks();}

  Bool_t      ParticleGood(Int_t i, Int_t pidi=0) const;
  Double_t    ParticleP(Int_t i)     const {return AliESDEvent::GetTrack(i)->GetTPCInnerParam()->P();}
  Double_t    ParticleTheta(Int_t i) const {return AliESDEvent::GetTrack(i)->GetTPCInnerParam()->Theta();}
  Double_t    ParticlePhi(Int_t i)   const {return TVector2::Phi_mpi_pi(AliESDEvent::GetTrack(i)->GetTPCInnerParam()->Phi());}
  Double_t    ParticleDedx(Int_t i)  const {return AliESDEvent::GetTrack(i)->GetTPCsignal()/47.0;}
  Bool_t      PairGood(Double_t p0, Double_t the0, Double_t phi0, 
		       Double_t p1, Double_t the1, Double_t phi1) const;
  // alternative: GetTPCInnerParam, GetConstrainedParam
  ClassDef(AliDEventAliceESD,0)
};
#endif 
//=============================================================================
