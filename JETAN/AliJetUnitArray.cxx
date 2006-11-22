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


/* $Id$ */

//_________________________________________________________________________
//  Unit used by UA1 algorithm
// --
//*-- Author: Sarah Blyth (LBL/UCT)
// --
// Revised Version for JETAN - 30/03/2006
// -- Magali Estienne (IReS)

#include "AliJetUnitArray.h"


ClassImp(AliJetUnitArray)

AliJetUnitArray::AliJetUnitArray()
{
  // Default constructor
  fUnitEnergy         = 0.0;
  fUnitEta            = 0.0;
  fUnitPhi            = 0.0;
  fUnitID             = 0;
  fUnitFlag           = kOutJet;
}  

AliJetUnitArray::~AliJetUnitArray()
{
  // Destructor 
}
	
Bool_t AliJetUnitArray::operator>(AliJetUnitArray unit) const
{
  // Greater than operator used by sort
  if( fUnitEnergy > unit.GetUnitEnergy())
    return kTRUE;
  else 
    return kFALSE;
}

Bool_t AliJetUnitArray::operator<( AliJetUnitArray unit) const
{
  // Less than operator used by sort
  if( fUnitEnergy < unit.GetUnitEnergy())
    return kTRUE;
  else
    return kFALSE;
}

Bool_t AliJetUnitArray::operator==( AliJetUnitArray unit) const
{
  // equality operator used by sort
  if( fUnitEnergy == unit.GetUnitEnergy())
    return kTRUE;
  else
    return kFALSE;
}
