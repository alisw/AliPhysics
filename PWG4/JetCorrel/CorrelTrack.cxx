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
/* $Id: $ */

//_____________________________________________________
// Container class for global tracks
// Use CorrelKFTrack_t to reconstruct parent with AliKF
//-- Author: Paul Constantin

#include "CorrelTrack.h"

using namespace std;

CorrelTrack_t::CorrelTrack_t() : CorrelParticle_t(), fTPCx(-999.), fTPCy(-999.), fTPCz(-999.) {
  // default constructor:
}

CorrelTrack_t::CorrelTrack_t(Float_t pt, Float_t p, Float_t e, Float_t m, cPartType_t i, Float_t x, Float_t y, Float_t z) : 
  CorrelParticle_t(pt,p,e,m,i), fTPCx(x), fTPCy(y), fTPCz(z){
  // constructor:
}

CorrelTrack_t::CorrelTrack_t(const CorrelTrack_t &p) : 
  CorrelParticle_t(p.fPt,p.fPhi,p.fEta,p.fMass,p.fID), fTPCx(p.fTPCx), fTPCy(p.fTPCy), fTPCz(p.fTPCz) {
  // copy constructor:
}

CorrelTrack_t* CorrelTrack_t::Copy(){
  // creates and returns a deep object copy
  CorrelTrack_t *copy = new CorrelTrack_t;
  copy->fPt    = this->Pt()*this->Q();
  copy->fPhi   = this->Phi();
  copy->fEta   = this->Eta();
  copy->fMass  = this->M();
  copy->fID    = this->ID();
  copy->fTPCx  = this->X();
  copy->fTPCy  = this->Y();
  copy->fTPCz  = this->Z();
  return copy;
}

Float_t CorrelTrack_t::Dist(CorrelTrack_t * const that) const {
  return TMath::Sqrt((this->X()-that->X())*(this->X()-that->X()) +
		     (this->Y()-that->Y())*(this->Y()-that->Y()) +
		     (this->Z()-that->Z())*(this->Z()-that->Z()));
}

void CorrelTrack_t::Show() const {
  std::cout<<" Track pT="<<Pt()<<" phi="<<Phi()<<" eta="<<Eta()<<" m="<<M()<<" id="<<ID()
	   <<" x="<<X()<<" y="<<Y()<<" z="<<Z()<<std::endl;
}
