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
//     AOD track base class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include "AliAODJet.h"

ClassImp(AliAODJet)


//______________________________________________________________________________
AliAODJet::AliAODJet() 
{
  // constructor
}

//______________________________________________________________________________
AliAODJet::~AliAODJet() 
{
  // destructor
}

//______________________________________________________________________________
AliAODJet::AliAODJet(const AliAODJet& jet) :
  AliVirtualParticle(jet)
{
  // Copy constructor

}

//______________________________________________________________________________
AliAODJet& AliAODJet::operator=(const AliAODJet& jet)
{
  // Assignment operator
  if(this!=&jet) {
  }

  return *this;
}

