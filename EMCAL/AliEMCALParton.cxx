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

/* $Id:  */
/*
  $Log$
  Revision 1.1.2.1  2002/08/28 15:06:50  alibrary
  Updating to v3-09-01

  Revision 1.1  2002/08/21 10:29:29  schutz
  New classes (by Renan)

*/

//*-- Author: Renan Cabrera (Creighton U.)

#include "AliEMCALParton.h"
#include "Ecommon.h"
  
ClassImp(AliEMCALParton)   
    
//____________________________________________________________________________
AliEMCALParton::AliEMCALParton()
{
  // Default constructor
}

AliEMCALParton::AliEMCALParton(Float_t energy, Float_t phi, Float_t eta)
{
  // Constructor
  fEnergy = energy;
  fPhi    = phi;
  fEta    = eta;
}

//____________________________________________________________________________

AliEMCALParton::~AliEMCALParton()
{
  // Destructor
}
