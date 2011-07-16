/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD class for photons
//     Author: Yves Schutz, CERN
//-------------------------------------------------------------------------

#include <TLorentzVector.h>
#include "AliAODPhoton.h"

ClassImp(AliAODPhoton)


//______________________________________________________________________________
AliAODPhoton::AliAODPhoton() :
    AliVParticle(),
    fMomentum(0)
{
  // constructor
}

AliAODPhoton::AliAODPhoton(Double_t px, Double_t py, Double_t pz, Double_t e):
    AliVParticle(),
    fMomentum(0)
{
  // constructor
    fMomentum = new TLorentzVector(px, py, pz, e);
}

AliAODPhoton::AliAODPhoton(TLorentzVector & p):
    AliVParticle(),
    fMomentum(0)
{
  // constructor
    fMomentum = new TLorentzVector(p);
}


//______________________________________________________________________________
AliAODPhoton::~AliAODPhoton() 
{
  // destructor
    delete fMomentum;
}

//______________________________________________________________________________
AliAODPhoton::AliAODPhoton(const AliAODPhoton& photon) :
    AliVParticle(photon),
    fMomentum(0)
{
  // Copy constructor
    fMomentum = new TLorentzVector(*photon.fMomentum);
    
}

//______________________________________________________________________________
AliAODPhoton& AliAODPhoton::operator=(const AliAODPhoton& photon)
{
  // Assignment operator
  if(this!=&photon) {
  }

  return *this;
}

void AliAODPhoton::Print(Option_t* /*option*/) const 
{
  // Print information of all data members
  printf("Photon 4-vector:\n");
  printf("     E  = %13.3f\n", E() );
  printf("     Px = %13.3f\n", Px());
  printf("     Py = %13.3f\n", Py());
  printf("     Pz = %13.3f\n", Pz());
}
