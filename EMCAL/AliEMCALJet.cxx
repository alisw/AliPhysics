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

//*-- Author: Andreas Morsch (CERN)

#include "AliEMCALJet.h"
#include "Ecommon.h"

ClassImp(AliEMCALJet)

//____________________________________________________________________________
AliEMCALJet::AliEMCALJet()
{
// Default constructor
}

AliEMCALJet::AliEMCALJet(Float_t energy, Float_t phi, Float_t eta)
{
// Constructor
    fEnergy = energy;
    fPhi    = phi;
    fEta    = eta;
}

//____________________________________________________________________________
AliEMCALJet::~AliEMCALJet()
{
// Destructor
}











