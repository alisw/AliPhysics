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
//     AOD class for di-jets
//     The present version is for test purposes only
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TLorentzVector.h>
#include "AliAODDiJet.h"

ClassImp(AliAODDiJet)


//______________________________________________________________________________
AliAODDiJet::AliAODDiJet() :
    AliAODJet(),
    fJetR(0),
    fJet1(0),
    fJet2(0)
{
  // constructor
}

AliAODDiJet::AliAODDiJet(Double_t px, Double_t py, Double_t pz, Double_t e):
    AliAODJet(px, py, pz, e), 
    fJetR(new TRefArray(2)),
    fJet1(0),
    fJet2(0)
{
  // another constructor
}

AliAODDiJet::AliAODDiJet(TLorentzVector & p):
    AliAODJet(p),
    fJetR(new TRefArray(2)),
    fJet1(0),
    fJet2(0)
{
  // constructor
}


//______________________________________________________________________________
AliAODDiJet::~AliAODDiJet() 
{
  // destructor
}
