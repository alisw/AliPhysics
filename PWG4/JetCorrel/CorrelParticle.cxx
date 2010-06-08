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

//________________________________________________
// Main container class - stores generic particle.
//-- Author: Paul Constantin

#include "CorrelParticle.h"

using namespace std;

CorrelParticle_t::CorrelParticle_t() : fPt(-999.), fPhi(-999.), fEta(-999.), fMass(-999.), fID(t_unknown){
  // default constructor
}

CorrelParticle_t::CorrelParticle_t(Float_t pt, Float_t p, Float_t t, Float_t m, cPartType_t i) : fPt(pt), fPhi(p), fEta(t), fMass(m), fID(i){
  // constructor
}

CorrelParticle_t::CorrelParticle_t(const CorrelParticle_t& p) : fPt(p.fPt), fPhi(p.fPhi), fEta(p.fEta), fMass(p.fMass), fID(p.fID){
  // copy constructor
}

CorrelParticle_t* CorrelParticle_t::Copy(){
  // creates and returns a deep object copy
  CorrelParticle_t *copy = new CorrelParticle_t;
  copy->fPt    = this->Pt()*this->Q();
  copy->fPhi   = this->Phi();
  copy->fEta   = this->Eta();
  copy->fMass  = this->M();
  copy->fID    = this->ID();
  return copy;
}

void CorrelParticle_t::Show() const {
  // printout method
  std::cout<<" Particle pT="<<Pt()<<" phi="<<Phi()<<" eta="<<Eta()<<" m="<<M()<<" id="<<ID()<<std::endl;
}
