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

//__________________________________________________
// Container class for a track with AliKF parameters
//-- Author: Paul Constantin

#include "CorrelKFTrack.h"

using namespace std;

CorrelKFTrack_t::CorrelKFTrack_t() : CorrelParticle_t(), fParam(NULL), fCovar(NULL){
  // default constructor:
}

CorrelKFTrack_t::CorrelKFTrack_t(Float_t pt, Float_t p, Float_t e, Float_t m, cPartType_t i, 
				   Double_t* par, Double_t* cov) : 
  CorrelParticle_t(pt,p,e,m,i), fParam(par), fCovar(cov){
  // constructor:
}

CorrelKFTrack_t::CorrelKFTrack_t(const CorrelKFTrack_t &p) : 
  CorrelParticle_t(p.fPt,p.fPhi,p.fEta,p.fMass,p.fID), fParam(p.fParam), fCovar(p.fCovar){
  // copy constructor:
}

CorrelKFTrack_t& CorrelKFTrack_t::operator=(const CorrelKFTrack_t& rhs){
  fPt      = rhs.Pt()*rhs.Q();
  fPhi     = rhs.Phi();
  fEta     = rhs.Eta();
  fMass    = rhs.M();
  fID      = rhs.ID();
  fParam   = rhs.Param();
  fCovar   = rhs.Covar();
  return *this;
}

CorrelKFTrack_t* CorrelKFTrack_t::Copy(){
  // creates and returns a copy object
  CorrelKFTrack_t *copy = new CorrelKFTrack_t;
  copy->fPt    = this->Pt()*this->Q();
  copy->fPhi   = this->Phi();
  copy->fEta   = this->Eta();
  copy->fMass  = this->M();
  copy->fID    = this->ID();
  copy->fParam = this->Param();
  copy->fCovar = this->Covar();
  return copy;
}

void CorrelKFTrack_t::Show() const {
  // printout method
  std::cout<<" Electron pT="<<Pt()<<" phi="<<Phi()<<" eta="<<Eta()<<" m="<<M()<<" id="<<ID()
	   <<" param="<<Param()<<" covar="<<Covar()<<std::endl;
}
