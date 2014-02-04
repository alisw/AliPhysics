/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//====================================================================================================================================================
//
//      ESD description of an ALICE muon forward track, combining the information of the Muon Spectrometer and the Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliESDMuonGlobalTrack.h"
#include "AliESDEvent.h"

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TDatabasePDG.h"

ClassImp(AliESDMuonGlobalTrack)

//====================================================================================================================================================

AliESDMuonGlobalTrack::AliESDMuonGlobalTrack():
  AliVParticle(),
  fCharge(0),
  fMatchTrigger(0),
  fPx(0), 
  fPy(0), 
  fPz(0), 
  fPt(0), 
  fP(0), 
  fEta(0), 
  fRapidity(0),
  fChi2(0),
  fChi2MatchTrigger(0),
  fLabel(-1),
  fESDEvent(0)
{

  //  Default constructor

}

//====================================================================================================================================================

AliESDMuonGlobalTrack::AliESDMuonGlobalTrack(Double_t px, Double_t py, Double_t pz):
  AliVParticle(),
  fCharge(0),
  fMatchTrigger(0),
  fPx(0), 
  fPy(0), 
  fPz(0), 
  fPt(0), 
  fP(0), 
  fEta(0), 
  fRapidity(0),
  fChi2(0),
  fChi2MatchTrigger(0),
  fLabel(-1),
  fESDEvent(0)
{

  //  Constructor with kinematics

  SetPxPyPz(px, py, pz);

}

//====================================================================================================================================================

AliESDMuonGlobalTrack::AliESDMuonGlobalTrack(const AliESDMuonGlobalTrack& muonTrack):
  AliVParticle(muonTrack),
  fCharge(muonTrack.fCharge),
  fMatchTrigger(muonTrack.fMatchTrigger),
  fPx(muonTrack.fPx), 
  fPy(muonTrack.fPy), 
  fPz(muonTrack.fPz), 
  fPt(muonTrack.fPt), 
  fP(muonTrack.fP), 
  fEta(muonTrack.fEta), 
  fRapidity(muonTrack.fRapidity),
  fChi2(muonTrack.fChi2),
  fChi2MatchTrigger(muonTrack.fChi2MatchTrigger),
  fLabel(muonTrack.fLabel),
  fESDEvent(muonTrack.fESDEvent)
{

  // Copy constructor
  
}

//====================================================================================================================================================

AliESDMuonGlobalTrack& AliESDMuonGlobalTrack::operator=(const AliESDMuonGlobalTrack& muonTrack) {

  // Assignment operator

  if (this == &muonTrack) return *this;

  // Base class assignement
  AliVParticle::operator=(muonTrack);

  fCharge           = muonTrack.fCharge;
  fMatchTrigger     = muonTrack.fMatchTrigger;
  fPx               = muonTrack.fPx; 
  fPy               = muonTrack.fPy; 
  fPz               = muonTrack.fPz; 
  fPt               = muonTrack.fPt; 
  fP                = muonTrack.fP;
  fEta              = muonTrack.fEta;
  fRapidity         = muonTrack.fRapidity;
  fChi2             = muonTrack.fChi2;
  fChi2MatchTrigger = muonTrack.fChi2MatchTrigger;
  fLabel            = muonTrack.fLabel;
  fESDEvent         = muonTrack.fESDEvent;

  return *this;

}

//====================================================================================================================================================

void AliESDMuonGlobalTrack::Copy(TObject &obj) const {
  
  // This overwrites the virtual TObject::Copy()
  // to allow run time copying without casting
  // in AliESDEvent

  if (this==&obj) return;
  AliESDMuonGlobalTrack *robj = dynamic_cast<AliESDMuonGlobalTrack*>(&obj);
  if (!robj) return; // not an AliESDMuonGlobalTrack
  *robj = *this;

}

//====================================================================================================================================================

void AliESDMuonGlobalTrack::SetPxPyPz(Double_t px, Double_t py, Double_t pz) {

  Double_t mMu = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
  Double_t eMu = TMath::Sqrt(mMu*mMu + px*px + py*py + pz*pz);

  TLorentzVector kinem(px, py, pz, eMu);

  fPx       =  kinem.Px();
  fPy       =  kinem.Py();
  fPz       =  kinem.Pz();
  fP        =  kinem.P();
  fPt       =  kinem.Pt();
  fEta      =  kinem.Eta();
  fRapidity =  kinem.Rapidity(); 

}

//====================================================================================================================================================
