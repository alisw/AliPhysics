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


/*
 
$Log$
Revision 1.1.1.1  2003/05/29 18:56:35  horner
Initial import - Mark

 
*/

//_________________________________________________________________________
//  Unit used by UA1 algorithm
//
//*-- Author: Sarah Blyth (LBL/UCT)
//


#include "AliEMCALJetFinderAlgoUA1Unit.h"


ClassImp(AliEMCALJetFinderAlgoUA1Unit)

AliEMCALJetFinderAlgoUA1Unit::AliEMCALJetFinderAlgoUA1Unit()
    {
      fUnitEnergy         = 0.0;
      fUnitEta            = 0.0;
      fUnitPhi            = 0.0;
      fUnitID             = 0;
      fUnitFlag           = kOutJet;
    }  

AliEMCALJetFinderAlgoUA1Unit::~AliEMCALJetFinderAlgoUA1Unit()
{

}
	
Bool_t AliEMCALJetFinderAlgoUA1Unit::operator>(AliEMCALJetFinderAlgoUA1Unit unit)
{
  if( fUnitEnergy > unit.GetUnitEnergy())
    return kTRUE;
  else 
    return kFALSE;
}

Bool_t AliEMCALJetFinderAlgoUA1Unit::operator<( AliEMCALJetFinderAlgoUA1Unit unit)
{
  if( fUnitEnergy < unit.GetUnitEnergy())
     return kTRUE;
  else
     return kFALSE;
}

Bool_t AliEMCALJetFinderAlgoUA1Unit::operator==( AliEMCALJetFinderAlgoUA1Unit unit)
    {
      if( fUnitEnergy == unit.GetUnitEnergy())
	 return kTRUE;
      else
	 return kFALSE;
    }
