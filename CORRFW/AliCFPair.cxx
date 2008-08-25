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


////////////////////////////////////////////////
// Class to handle pairs of tracks
// Useful for resonance analysis
// Derives from AliVParticle => 
// usable in Correction Framework
////////////////////////////////////////////////
// author : renaud.vernet@cern.ch
////////////////////////////////////////////////

#include "AliCFPair.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliESDEvent.h"
#include "TMath.h"
#include "AliAODv0.h"

ClassImp(AliCFPair)

AliCFPair::AliCFPair(AliVParticle* t1, AliVParticle* t2) :
  AliVParticle(),
  fIsV0(0),
  fTrackNeg(t1),
  fTrackPos(t2),
  fESDV0(0x0),
  fAODV0(0x0),
  fLabel(-1),
  fV0PDG(0)
{
  //  
  // 2-track ctor
  //
}
AliCFPair::AliCFPair(AliESDv0* v0, AliESDEvent* esd) :
  AliVParticle(),
  fIsV0(1),
  fTrackNeg(esd->GetTrack(v0->GetNindex())),
  fTrackPos(esd->GetTrack(v0->GetPindex())),
  fESDV0(v0),
  fAODV0(0x0),
  fLabel(-1),
  fV0PDG(0)
{
  //  
  // V0 ctor
  //
}
AliCFPair::AliCFPair(AliAODv0* v0) :
  AliVParticle(),
  fIsV0(1),
  fTrackNeg((AliVParticle*)v0->GetSecondaryVtx()->GetDaughter(1)),
  fTrackPos((AliVParticle*)v0->GetSecondaryVtx()->GetDaughter(0)),
  fESDV0(0x0),
  fAODV0(v0),
  fLabel(-1),
  fV0PDG(0)
{
  //  
  // V0 ctor
  //
}
AliCFPair::AliCFPair(const AliCFPair& c) :
  AliVParticle(c),
  fIsV0(c.fIsV0),
  fTrackNeg(c.fTrackNeg),
  fTrackPos(c.fTrackPos),
  fESDV0(c.fESDV0),
  fAODV0(c.fAODV0),
  fLabel(c.fLabel),
  fV0PDG(c.fV0PDG)
{
  // 
  // Copy constructor.
  // 
}
AliCFPair& AliCFPair::operator=(const AliCFPair& c) {
  // 
  // assignment operator.
  // 
  
  if (this!=&c) {
    AliVParticle::operator=(c);
    fIsV0 = c.fIsV0;
    fTrackNeg = c.fTrackNeg ;
    fTrackPos = c.fTrackPos ;
    fESDV0 = c.fESDV0 ;
    fAODV0 = c.fAODV0 ;
    fLabel = c.fLabel ;
    fV0PDG = c.fV0PDG ;
  }
  return *this;
}
Bool_t AliCFPair::PxPyPz(Double_t p[3]) const {
  //
  // sets pair total momentum in vector p
  //
  if (fIsV0) {
    if      (fESDV0) fESDV0->GetPxPyPz(p[0],p[1],p[2]);
    else if (fAODV0) fAODV0->PxPyPz(p);
    else Error("PxPyPz","Pointer to V0 not found");
  }
  else if (fTrackNeg && fTrackPos) {
    p[0]=fTrackNeg->Px()+fTrackPos->Px();
    p[1]=fTrackNeg->Py()+fTrackPos->Py();
    p[2]=fTrackNeg->Pz()+fTrackPos->Pz();
  }
  else Error("PxPyPz","Could not find V0 nor track pointers");
  return kTRUE;
}

Double32_t AliCFPair::P() const {
  //
  // returns pair total momentum norm
  //
  Double32_t mom[3];
  PxPyPz(mom);
  return TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
}

Double32_t AliCFPair::Px() const {
  //
  // returns pair X-projected momentum
  //
  Double32_t mom[3];
  PxPyPz(mom);
  return mom[0];
}
Double32_t AliCFPair::Py() const {
  //
  // returns pair Y-projected momentum
  //
  Double32_t mom[3];
  PxPyPz(mom);
  return mom[1];
}
Double32_t AliCFPair::Pz() const {
  //
  // returns pair Z-projected momentum
  //
  Double32_t mom[3];
  PxPyPz(mom);
  return mom[2];
}
Double32_t AliCFPair::Pt() const {
  //
  // returns pair transverse (XY) momentum
  //
  Double32_t mom[3];
  PxPyPz(mom);
  return sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
}
Double32_t AliCFPair::E() const {
  //
  // returns pair total energy according to ESD-calculated mass
  //
  Double32_t mom[3];
  PxPyPz(mom);
  Double32_t mass=M() ;
  return TMath::Sqrt(mass*mass + mom[0]*mom[0]+ mom[1]*mom[1]+ mom[2]*mom[2]);
}
Bool_t AliCFPair::XvYvZv(Double_t x[3]) const {
  //
  // sets pair position to x
  // since this class is designed for resonances, the assumed pair position
  // should be the same for both tracks. neg track position is kept here
  //
  
  if (fIsV0) {
    if      (fESDV0) fESDV0->GetXYZ(x[0],x[1],x[2]);
    else if (fAODV0) fAODV0->XvYvZv(x);
    else Error("PxPyPz","Pointer to V0 not found");
  }
  else if (fTrackNeg) fTrackNeg->PxPyPz(x);
  else Error("PxPyPz","Could not find V0 nor track pointers");

  return kTRUE;
}
Double32_t AliCFPair::Xv() const {
  //
  // returns pair X-projected position
  //
  Double32_t pos[3];
  XvYvZv(pos);
  return pos[0];
}
Double32_t AliCFPair::Yv() const {
  //
  // returns pair Y-projected position
  //
  Double32_t pos[3];
  XvYvZv(pos);
  return pos[1];
}
Double32_t AliCFPair::Zv() const {
  //
  // returns pair Z-projected position
  //
  Double32_t pos[3];
  XvYvZv(pos);
  return pos[2];
}
Double32_t AliCFPair::Phi() const {
  //
  // returns pair phi angle (in transverse plane)
  //
  return TMath::Pi()+TMath::ATan2(-Py(),-Px()); 
}
Double32_t AliCFPair::Theta() const { 
  //
  // returns pair theta angle (in YZ plane)
  //
  return (Pz()==0)?TMath::PiOver2():TMath::ACos(Pz()/P());
}
Double32_t AliCFPair::Eta() const {
  //
  // returns pair pseudo-rapidity
  //
  Double32_t pmom = P();
  Double32_t pz = Pz();
  if (pmom != TMath::Abs(pz)) return 0.5*TMath::Log((pmom+pz)/(pmom-pz));
  else                        return 999;
}
Double32_t AliCFPair::Y() const {
  //
  // returns pair rapidity
  //
  Double32_t e  = E();
  Double32_t pz = Pz();

  if (e == pz || e == -pz) {
    printf("GetRapidity : ERROR : rapidity for 4-vector with E = Pz -- infinite result");
    return 999;
  }
  if (e < pz) {
    printf("GetRapidity : ERROR : rapidity for 4-vector with E = Pz -- infinite result");
    return 999;
  }
  Double32_t y = 0.5 * log((e + pz) / (e - pz));
  return y;
} 
Double32_t AliCFPair::M() const {
  //
  // returns pair invariant mass
  // in case of a V0, returns the current mass hypothesis
  // otherwise returns ESD-calculated mass
  //

  Double32_t minv = 0. ;

  if (fIsV0) {
    if      (fESDV0) {
      fESDV0->ChangeMassHypothesis(fV0PDG);
      minv = (Double32_t)fESDV0->GetEffMass();
    }
    else if (fAODV0) {
      switch (fV0PDG) {
      case 310  : minv = fAODV0->MassK0Short()    ; break ;
      case 3122 : minv = fAODV0->MassLambda()     ; break ;
      case -3122: minv = fAODV0->MassAntiLambda() ; break ;
      default:    minv = -1.              ; break ;
      }
    }
    else Error("PxPyPz","Pointer to V0 not found");
  }
  else if (fTrackNeg && fTrackPos) {
    Double32_t p  = P() ;
    Double32_t e = fTrackNeg->E() + fTrackPos->E() ;
    minv = TMath::Sqrt(e*e-p*p);
  }
  else Error("M","Could not find V0 nor track pointers");
  
  return minv ;
}
