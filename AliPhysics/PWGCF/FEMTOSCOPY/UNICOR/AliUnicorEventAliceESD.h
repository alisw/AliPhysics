#ifndef ALIUNICOREVENTALICEESD_H
#define ALIUNICOREVENTALICEESD_H

/* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2008

#include <cmath>
#include "TVector2.h"
#include "AliESDEvent.h"
//#include "AliPhysicsSelection.h"
#include "AliExternalTrackParam.h"
#include "AliUnicorEvent.h"

//class AliPhysicsSelection;
class TGraph;

//=============================================================================
class AliUnicorEventAliceESD : public AliUnicorEvent {

 public:
              AliUnicorEventAliceESD(AliESDEvent *esd=0);
	      //AliUnicorEventAliceESD(const AliUnicorEventAliceESD &ev): AliUnicorEvent(ev), fViper(ev.fViper), fESD(ev.fESD), fPhysicsSelection(ev.fPhysicsSelection){}
              AliUnicorEventAliceESD(const AliUnicorEventAliceESD &ev): AliUnicorEvent(ev), fViper(ev.fViper), fESD(ev.fESD) {}
  virtual     ~AliUnicorEventAliceESD();
  AliUnicorEventAliceESD &operator=(const AliUnicorEventAliceESD &source) {
    if(&source == this) return *this;
    AliUnicorEvent::operator=(source);
    fViper=source.fViper; fESD=source.fESD; return *this;}
  Double_t    Etamin() const {return -0.75;}
  Double_t    Etamax() const {return  0.75;}
  void        AttachTree(TTree *tr) {fESD->ReadFromTree(tr);}
  Bool_t      Good() const;
  //  Double_t    Centrality() const {return 0.9999*exp(-NGoodParticles()/1000.0);} // 7 for pp 900 GeV, 12 for pp 7 TeV, 1000  PbPb
  Double_t    Centrality() const; 
  void        RP(Double_t &qx, Double_t &qy) const {AliUnicorEvent::RP(qx,qy,2);}
  Double_t    RPphi() const {Double_t qx,qy; RP(qx,qy); return atan2(qy,qx);}
  Double_t    Zver() const {return fESD->GetPrimaryVertex()->GetZ()/12.0;}
  Int_t       NParticles() const {return fESD->GetNumberOfTracks();}

  Bool_t      ParticleGood(Int_t i, Int_t pidi=0) const;
  Double_t    ParticleP(Int_t i)     const {return GetTrackParam(i)->P();}
  Double_t    ParticleTheta(Int_t i) const {return GetTrackParam(i)->Theta();}
  Double_t    ParticlePhi(Int_t i)   const {return TVector2::Phi_mpi_pi(GetTrackParam(i)->Phi());}
  Double_t    ParticleDedx(Int_t i)  const {return fESD->GetTrack(i)->GetTPCsignal()/50.0;}
  Bool_t      PairGood(double p0, double the0, double phi0, double z0, 
		       double p1, double the1, double phi1, double z1) const;
  void        SetESD(AliESDEvent * const esd) {fESD = esd;}
  AliESDEvent *GetESD() const {return fESD;}
  //const AliExternalTrackParam *GetTrackParam(Int_t i) const {return fESD->GetTrack(i);}
  //const AliExternalTrackParam *GetTrackParam(Int_t i) const {return fESD->GetTrack(i)->GetConstrainedParam();}
  //const AliExternalTrackParam *GetTrackParam(Int_t i) const {return fESD->GetTrack(i)->GetInnerParam();} // not at vtx!
  const AliExternalTrackParam *GetTrackParam(Int_t i) const {return fESD->GetTrack(i)->GetTPCInnerParam();}

 protected:
  TGraph      *fViper;                    // V0 centrality percentile
  AliESDEvent *fESD;                      //! pointer to the actual source of data
  //  AliPhysicsSelection *fPhysicsSelection; //! interaction event filter

  ClassDef(AliUnicorEventAliceESD,0)
};
#endif 
//=============================================================================
