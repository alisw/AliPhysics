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

//_________________________________________________________________________
// PHOSRecPoint base class deriving from AliRecPoint
//*-- Author : Gines MARTINEZ  SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h"
#include "AliPHOSRecPoint.h"

ClassImp(AliPHOSRecPoint)


//____________________________________________________________________________
AliPHOSRecPoint::AliPHOSRecPoint()
  : AliRecPoint()
{
  fGeom =   AliPHOSGeometry::GetInstance() ;
  fPHOSMod = 0;
}

//____________________________________________________________________________
AliPHOSRecPoint::~AliPHOSRecPoint()
{
  // dtor
}

//____________________________________________________________________________
Int_t AliPHOSRecPoint::GetPHOSMod()
{ 
  if(fPHOSMod > 0) return fPHOSMod ;
  Int_t relid[4] ;
  
  AliPHOSDigit * digit   ;
  digit = (AliPHOSDigit *) fDigitsList[0] ;
  AliPHOSGeometry * PHOSGeom =  (AliPHOSGeometry *) fGeom ;

  PHOSGeom->AbsToRelNumbering(digit->GetId(), relid) ;
  fPHOSMod = relid[0];
  return fPHOSMod ;
}
